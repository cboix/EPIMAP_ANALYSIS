import java.io.BufferedReader;
import java.io.*;
import java.text.*;
import java.util.*;

// import org.apache.commons.math3.stat.inference.*;
// import org.apache.commons.math3.distribution.*;
// import org.apache.commons.math3.linear.*;

public class StateMWGwasPeakHyper {

  public static Hashtable htbinom = new Hashtable();

  /** Returns the log of the binomial coefficient N choose ni */
  public static double logbinomcoeff(int ni, int N) {
    String sz = ni + ";" + N;
    Double dobj = (Double) htbinom.get(sz);
    double dsum;

    if (dobj != null) {
      dsum = ((Double) dobj).doubleValue();
    } else {
      dsum = 0;
      int dmax = Math.max(ni, N - ni);
      int dmin = Math.min(ni, N - ni);

      // the log of the part of the numerator not cancelled by
      // the larger factorial in the denominator
      for (int nj = dmax + 1; nj <= N; nj++) {
        dsum += Math.log(nj);
      }

      // subtract off the log of the denominator of the smaller term
      for (int nj = 2; nj <= dmin; nj++) {
        dsum -= Math.log(nj);
      }
      // store it
      htbinom.put(sz, new Double(dsum));
    }

    return dsum;
  }

  /**
   * Returns the probability of seeing more than x objects of type A, when there are nA objects of
   * type A nB objects of type B, and nm objects total drawn This can be used to compute a more
   * accurate p-values than a 1-cumulative probability calculation
   */
  public static double hypergeometrictail(int nx, int nA, int nB, int nm) {
    if (nx < 0) {
      return 1;
    }

    double dprob = 0;

    int nminx = Math.max(nx + 1, 0); // the first element to the right of x or 0
    int nmaxx = Math.min(nA, nm); // the max of type A there can be

    if (nminx > nmaxx) {
      // if nx approaches infinity tail should be 0
      return 0;
    }

    double dsum =
        -logbinomcoeff(nm, nA + nB) + logbinomcoeff(nminx, nA) + logbinomcoeff(nm - nminx, nB);

    double dlogprob = dsum;
    for (int ni = nminx + 1; ni <= nmaxx; ni++) {
        dsum += Math.log(nA - ni + 1) - Math.log(nB - nm + ni) + Math.log(nm - ni + 1) - Math.log(ni);

        if (dsum >= dlogprob) {
            dlogprob = dsum + Math.log(1 + Math.pow(Math.E, dlogprob - dsum));
      } else {
        dlogprob = dlogprob + Math.log(1 + Math.pow(Math.E, dsum - dlogprob));
      }
    }

    dprob = Math.pow(Math.E, dlogprob);

    if (dprob <= 0) {
      return 0;
    } else if (dprob >= 1) {
      return 1;
    } else {
      return dprob;
    }
  }

  static class RecBackFore {

    double dval;
    boolean bfore;

    RecBackFore(double dval, boolean bfore) {

      this.dval = dval;
      this.bfore = bfore;
    }
  }

  static class RecBackForeCompare implements Comparator {

    public int compare(Object o1, Object o2) {
      RecBackFore r1 = (RecBackFore) o1;
      RecBackFore r2 = (RecBackFore) o2;

      if (r1.dval > r2.dval) {
        return -1;
      } else if (r1.dval < r2.dval) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  static class Rec {
    int nbegin;
    double dpval;

    Rec(int nbegin, double dpval) {

      this.nbegin = nbegin;
      this.dpval = dpval;
    }
  }

  static class RecCompare implements Comparator {

    public int compare(Object o1, Object o2) {
      Rec r1 = (Rec) o1;
      Rec r2 = (Rec) o2;

      if (r1.dpval < r2.dpval) {
        return -1;
      } else if (r1.dpval > r2.dpval) {
        return 1;
      } else if (r1.nbegin < r2.nbegin) {
        return -1;
      } else if (r1.nbegin > r2.nbegin) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  // Main load datasets + format:
  public static void main(String[] args) throws IOException {
    String[] chroms = {
      "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
      "18", "19", "20", "21", "22", "X"
    };

    // /broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/n15/enh/clust/n15_enh_k300_cls_300/n15_enh_k300_
    String prefix = args[0];
    String outprefix = args[1];
    String gwascatalog = args[2];
    // Extend to X bp?
    int extend = Integer.parseInt(args[3]);
    int NUMCLS = Integer.parseInt(args[4]);
    int halfext = extend / 2;
    // "gwascatalog_dec21_2018.txt"; 

    // OLD ARGS:
    //String szmark = args[0];
    // int ntype = Integer.parseInt(args[1]);

    // Distance between two SNPs, keep top.
    int DIST = 1000000; 

    // Set up output file:
    PrintWriter pw = null;
    String extstr = String.valueOf(extend);
    pw = new PrintWriter(new FileWriter(outprefix + "_" + extstr + "_enrich.tsv"));

    // MannWhitneyUTest theMannWhitneyUTest = new MannWhitneyUTest();

    NumberFormat nf = NumberFormat.getInstance();
    nf.setMaximumFractionDigits(2);
    nf.setGroupingUsed(false);
    Random theRandom = new Random(2433);

    HashSet hscategory = new HashSet();
    HashMap hmcategory = new HashMap();
    ArrayList alcategory = new ArrayList();
    // NOTE: GWAS should be formatted as 2014 version
    BufferedReader brgwas = new BufferedReader(new FileReader(gwascatalog));
    brgwas.readLine();
    String szLine;

    int ncategory = 0;
    while ((szLine = brgwas.readLine()) != null) {
      StringTokenizer st = new StringTokenizer(szLine, "\t");
      st.nextToken();
      st.nextToken();
      st.nextToken();
      st.nextToken();
      st.nextToken();
      String szpubmed = st.nextToken();
      st.nextToken();
      st.nextToken();
      st.nextToken();
      st.nextToken();
      String sztrait = st.nextToken();

      if (!hscategory.contains(szpubmed + "\t" + sztrait)) {
        alcategory.add(szpubmed + "\t" + sztrait);
        hscategory.add(szpubmed + "\t" + sztrait);
        hmcategory.put(szpubmed + "\t" + sztrait, new Integer(ncategory));
        ncategory++;
      }
    }
    brgwas.close();

    for (int ncls = 0; ncls < NUMCLS; ncls++){
      String szcls = String.valueOf(ncls);
      System.out.println(szcls);
      String szfilename = prefix + szcls + ".bed"; 

      int ngwasbackgroundcountsHIT = 0;
      int[] ngwasforegroundcountsHIT = new int[alcategory.size()];
      int ngwasbackgroundcountsALL = 0;
      int[] ngwasforegroundcountsALL = new int[alcategory.size()];

      ArrayList gwasbackgroundcounts = new ArrayList();
      ArrayList[] gwasforegroundcounts = new ArrayList[alcategory.size()];
      for (int na = 0; na < gwasforegroundcounts.length; na++)
        gwasforegroundcounts[na] = new ArrayList();

      // Check file exists:
      File f = null;
      f = new File(szfilename);
      if (!f.exists()) continue;

      for (int nchrom = 0; nchrom < chroms.length; nchrom++) {
        String szcurrchrom = "chr" + chroms[nchrom];
        boolean[] stateA = new boolean[250000000];  // BOOLEAN OF MAX CHR LENGTH
        brgwas = new BufferedReader(new FileReader(gwascatalog));

        BufferedReader brchrom = null;
        brchrom = new BufferedReader(new FileReader(szfilename)); 

        String szLineA;
        while ((szLineA = brchrom.readLine()) != null) {
          StringTokenizer ststate = new StringTokenizer(szLineA, "\t");
          String szchrom = ststate.nextToken();
          if (szchrom.equals(szcurrchrom)) {
            int nbegin = Integer.parseInt(ststate.nextToken());
            int nend = Integer.parseInt(ststate.nextToken());
            if (extend != 0){
                int nmid = (nend + nbegin) / 2;
                nbegin = nmid - halfext;
                nend = nmid + halfext;
                if (nbegin < 0){ nbegin = 0; }
            }
            for (int nk = nbegin; nk < nend; nk++) {
              stateA[nk] = true;
            }
          }
        }
        brchrom.close();

        HashSet hsvisited = new HashSet();
        HashMap hmvisited = new HashMap();

        brgwas.readLine();
        while ((szLine = brgwas.readLine()) != null) {
          StringTokenizer st = new StringTokenizer(szLine, "\t", true);
          st.nextToken();
          st.nextToken(); // tab

          String szchrom = "chr" + st.nextToken();
          st.nextToken(); // tab

          if (szchrom.equals(szcurrchrom)) {
            // Note: Tokens read tabs as well, so read 2x # of lines.
            String szbegin = st.nextToken();
            st.nextToken(); // tab
            String szend = st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            String szpubmedid = st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            String sztrait = st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            String szval = st.nextToken();
            if (!szval.equals("\t")) st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            st.nextToken();
            st.nextToken(); // tab
            String szpval = st.nextToken();
            if (szpval.equals("NS")) continue;
            if (szpval.equals("E")) continue;
            if (szpval.equals("NR")) continue;

            ArrayList alvisit = (ArrayList) hmvisited.get(szpubmedid + "\t" + sztrait);
            if (alvisit == null) {
              alvisit = new ArrayList();
              hmvisited.put(szpubmedid + "\t" + sztrait, alvisit);
            }

            int nbegin = Integer.parseInt(szbegin);
            // System.out.println(szend + " " + szpubmedid + " " + sztrait + " " + szval + " " + szpval);
            double dpval = Double.parseDouble(szpval);
            alvisit.add(new Rec(nbegin, dpval));
          }
        }
        brgwas.close();

        Iterator itrcategory = (Iterator) hmvisited.keySet().iterator();

        while (itrcategory.hasNext()) {
          String szkey = (String) itrcategory.next();
          ArrayList alpval = (ArrayList) hmvisited.get(szkey);

          ncategory = ((Integer) hmcategory.get(szkey)).intValue();

          Rec[] therecA = new Rec[alpval.size()];

          for (int nm = 0; nm < therecA.length; nm++) {
            therecA[nm] = (Rec) alpval.get(nm);
          }

          Arrays.sort(therecA, new RecCompare());

          ArrayList alkeep = new ArrayList();

          for (int na = 0; na < therecA.length; na++) {
            boolean bok = true;
            for (int nb = 0; nb < alkeep.size(); nb++) {
              Rec currRec = (Rec) alkeep.get(nb);
              if (Math.abs(currRec.nbegin - therecA[na].nbegin) < DIST) bok = false;
            }

            if (bok) {
              alkeep.add(therecA[na]);
              int nwindow = therecA[na].nbegin;

              if (stateA[nwindow]) {
                ngwasbackgroundcountsHIT++;
                ngwasforegroundcountsHIT[ncategory]++;
              }
              ngwasbackgroundcountsALL++;
              ngwasforegroundcountsALL[ncategory]++;
            }
          }
        }
      }

      for (int ni = 0; ni < ngwasforegroundcountsHIT.length; ni++) {
        if (ngwasforegroundcountsHIT[ni] >= 1) {
          double dpval =
              StateMWGwasPeakHyper.hypergeometrictail(
                  ngwasforegroundcountsHIT[ni] - 1,
                  ngwasbackgroundcountsHIT,
                  ngwasbackgroundcountsALL - ngwasbackgroundcountsHIT,
                  ngwasforegroundcountsALL[ni]);
          double dfold =
              (ngwasforegroundcountsHIT[ni] / (double) ngwasforegroundcountsALL[ni])
                  / (ngwasbackgroundcountsHIT / (double) ngwasbackgroundcountsALL);
          pw.println(dpval + "\t" + szcls + "\t" + alcategory.get(ni) + "\t" + ngwasforegroundcountsHIT[ni] + "\t" + ngwasforegroundcountsALL[ni] + "\t" + dfold);
          //System.out.println(dpval + "\t" + szcls + "\t" + alcategory.get(ni) + "\t" + ngwasforegroundcountsHIT[ni] + "\t" + ngwasforegroundcountsALL[ni] + "\t" + dfold);
        }
      }
    }
    pw.close();
  }
}
