import java.io.*;
import java.text.*;
import java.util.*;

// import org.apache.commons.math3.stat.inference.*;
// import org.apache.commons.math3.distribution.*;
// import org.apache.commons.math3.linear.*;

public class StateMWGwasPeakHyperWithRandom {

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

  static class RandomRecCompare implements Comparator {
    public int compare(Object o1, Object o2) {
      Rec r1 = (Rec) o1;
      Rec r2 = (Rec) o2;
      if (r1.drandom < r2.drandom) {
        return -1;
      } else if (r1.drandom > r2.drandom) {
        return 1;
      } else {
        return 0;
      }
    }
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
      // computing the increase in probability mass
      // numerator has nA!/(ni!(nA-ni)!) * nB!/((nm-ni)!(nB-nm+ni)!)
      // denominator has (nA+nB)!/(nm!(nA+nB-nm)!)

      // numerator has nA!/((ni-1)!(nA-ni+1)!) * nB!/((nm-ni+1)!(nB-nm+ni-1)!)
      // denominator has (nA+nB)!/(nm!(nA+nB-nm)!)
      // cancelling gives
      // 1/(ni!(nA-ni)!) * 1/((nm-ni)!(nB-nm+ni)!) over
      // 1/((ni-1)!(nA-ni+1)!) * 1/((nm-ni+1)!(nB-nm+ni-1)!)
      dsum += Math.log(nA - ni + 1) - Math.log(nB - nm + ni) + Math.log(nm - ni + 1) - Math.log(ni);

      // log(a+b+c+d+e)
      // log(e) + log(a+b+c+d+e) - log(e)
      // log(e) + log((a+b+c+d+e)/e)
      // log(e) + log(1+(a+b+c+d)/e)
      // log(e) + log(1+Math.exp(log(a+b+c+d)-log(e)))

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

  /*
  static class Rec
  {
  int nbegin;
  double dpval;

  Rec(int nbegin, double dpval)
  {

  this.nbegin = nbegin;
  this.dpval = dpval;
  }
  }
  */

  static class Rec {
    int nbegin;
    double dpval;
    String szid;
    String szchrom;
    double drandom;

    Rec(String szid, String szchrom, int nbegin, double dpval, double drandom) {
      this.szid = szid;
      this.szchrom = szchrom;
      this.nbegin = nbegin;
      this.dpval = dpval;
      this.drandom = drandom;
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

  public static void main(String[] args) throws IOException {

    // Chromosomes
    String[] chroms = {
      "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
      "18", "19", "20", "21", "22", "X"
    };

    String szmark = args[0];
    int ntype = Integer.parseInt(args[1]);

    int NUMSTATES = 25;
    int WINDOW = 25;
    int DIST = 1000000;

    int SEED = Integer.parseInt(args[2]);
    ;
    Random theRandom = new Random(SEED);

    PrintWriter pw = null;
    if (ntype == 0) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_observed_" + szmark + "_" + SEED + "_narrow.txt"));
    } else if (ntype == 1) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_imputed_" + szmark + "_" + SEED + "_narrow.txt"));
    } else if (ntype == 2) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_observed_" + szmark + "_" + SEED + "_gapped.txt"));
    } else if (ntype == 3) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_imputed_" + szmark + "_" + SEED + "_gapped.txt"));
    } else if (ntype == 4) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_observed_" + szmark + "_" + SEED + "_broad.txt"));
    } else if (ntype == 5) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_observed_"
                      + szmark
                      + "_"
                      + SEED
                      + "_hotspot.all.peaks.txt"));
    } else if (ntype == 6) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_observed_"
                      + szmark
                      + "_"
                      + SEED
                      + "_hotspot.fdr0.01.peaks.txt"));
    } else if (ntype == 7) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_observed_"
                      + szmark
                      + "_"
                      + SEED
                      + "_macs2.narrowPeak.txt"));
    } else if (ntype == 8) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_observed_"
                      + szmark
                      + "_"
                      + SEED
                      + "_hotspot.broad.bed.txt"));
    } else if (ntype == 9) {
      pw =
          new PrintWriter(
              new FileWriter(
                  "FULL_GWAS_PEAKHYPER/random_observed_"
                      + szmark
                      + "_"
                      + SEED
                      + "_hotspot.fdr0.01.broad.bed.txt"));
    }

    // Random theRandom = new Random(SEED);

    /*
    {
    pw = new PrintWriter(new FileWriter("FULL_GWAS_MW/average_"+szmark+"_2sided.txt"));
    }
    */

    // MannWhitneyUTest theMannWhitneyUTest = new MannWhitneyUTest();

    NumberFormat nf = NumberFormat.getInstance();
    nf.setMaximumFractionDigits(2);
    nf.setGroupingUsed(false);
    // Random theRandom = new Random(2433);

    HashSet hscategory = new HashSet();
    HashMap hmcategory = new HashMap();
    ArrayList alcategory = new ArrayList();
    BufferedReader brgwas = new BufferedReader(new FileReader("gwascatalog_sep12_2014.txt"));
    brgwas.readLine();
    String szLine;

    int ncategory = 0;
    ArrayList algwasrecs = new ArrayList();

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

    int NUMCELLS = 129;
    for (int ncell = 1; ncell <= NUMCELLS; ncell++)
    // for (int ncell = 124; ncell <= 124; ncell++)
    {
      System.out.println(ncell);

      if (ncell == 60) continue;
      if (ncell == 64) continue;

      int ngwasbackgroundcountsHIT = 0;
      int[] ngwasforegroundcountsHIT = new int[alcategory.size()];
      int ngwasbackgroundcountsALL = 0;
      int[] ngwasforegroundcountsALL = new int[alcategory.size()];

      ArrayList gwasbackgroundcounts = new ArrayList(); // int[NUMSTATES+1];
      ArrayList[] gwasforegroundcounts = new ArrayList[alcategory.size()]; // [NUMSTATES+1];
      for (int na = 0; na < gwasforegroundcounts.length; na++)
        gwasforegroundcounts[na] = new ArrayList();

      String szcell;

      if (ncell <= 9) szcell = "E00" + ncell;
      else if (ncell <= 99) szcell = "E0" + ncell;
      else szcell = "E" + ncell;

      File f = null;

      // IMPUTEDPEAKS/E003-H3K18ac.imputed.gappedPeak.bed.gz
      // IMPUTEDPEAKS/E003-H3K79me1.imputed.gappedPeak.bed.gz
      // IMPUTEDPEAKS/E003-H3K18ac.imputed.narrowPeak.bed.gz
      // IMPUTEDPEAKS/E003-H3K79me1.imputed.narrowPeak.bed.g

      if (ntype == 0) {
        f = new File( "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/" + szcell + "-" + szmark + ".narrowPeak.gz");
        // "CONVERTED_ROADMAPRAWPVAL/chr"+chroms[0]+"_"+szcell+"-"+szmark+".pval.signal.bigwig.bedgraph.gz.wig.gz");
      } else if (ntype == 1) {
        f = new File("IMPUTEDPEAKS/" + szcell + "-" + szmark + ".imputed.narrowPeak.bed.gz");
        // f = new
        // File("IMPUTE_ALL/chr"+chroms[0]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
        // f = new
        // File("TIER2_IMPUTE_ALL/chr"+chroms[0]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
      } else if (ntype == 2) {
        f = new File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/" + szcell + "-" + szmark + ".gappedPeak.gz");
        ///
        // broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/E001-H3K9me3.gappedPeak.gz*
        // f = new
        // File("AVG_CONVERTED_ROADMAPRAWPVAL/chr"+chroms[0]+"_"+szmark+"_"+szcell+".wig.gz");
        // AVG_CONVERTED_ROADMAPRAWPVAL/chr10_H3K27me3_E001.wig.gz
      } else if (ntype == 3) {
        f = new File("IMPUTEDPEAKS/" + szcell + "-" + szmark + ".imputed.gappedPeak.bed.gz");
      } else if (ntype == 4) {
        f = new File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/" + szcell + "-" + szmark + ".broadPeak.gz");
      } else if (ntype == 5) {
        f = new File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/" + szcell + "-" + szmark + ".hotspot.all.peaks.bed.gz");
      } else if (ntype == 6) {
        f = new File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/" + szcell + "-" + szmark + ".hotspot.fdr0.01.peaks.bed.gz");
      } else if (ntype == 7) {
        f = new File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/" + szcell + "-" + szmark + ".macs2.narrowPeak.gz");
      } else if (ntype == 8) {
        f = new File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/" + szcell + "-" + szmark + ".hotspot.broad.bed.gz");
      } else if (ntype == 9) {
        f = new File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/" + szcell + "-" + szmark + ".hotspot.fdr0.01.broad.bed.gz");
      }

      // "CONVERTED_ROADMAPRAWPVAL/chr"+chroms[0]+"_"+szcell+"-"+szmark+".pval.signal.bigwig.bedgraph.gz.wig.gz");

      // IMPUTE_ALL/chr10_impute_E001_DNase_1_100000_20.wig.gz

      if (!f.exists()) continue;

      // System.out.println(szcell);

      /*
             HashMap hmscore = new HashMap();
          //for (int nchrom = 0; nchrom < chroms.length; nchrom++)
          {
          BufferedReader brchrom;//  = Util.getBufferedReader("IMPUTE_ALL/chr"+chroms[0]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
          if (ntype == 0)
          {
          //		f = new File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/"+szcell+"-"+szmark+".narrowPeak.gz");
          brchrom = Util.getBufferedReader("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/"+szcell+"-"+szmark+".narrowPeak.gz");

          //brchrom  = Util.getBufferedReader("CONVERTED_ROADMAPRAWPVAL/chr"+chroms[nchrom]+"_"+szcell+"-"+szmark+".pval.signal.bigwig.bedgraph.gz.wig.gz");
          }
          else if (ntype == 1)
          {
          brchrom = Util.getBufferedReader("IMPUTEDPEAKS/"+szcell+"-"+szmark+".imputed.narrowPeak.bed.gz");
          //		    brchrom  = Util.getBufferedReader("IMPUTE_ALL/chr"+chroms[nchrom]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
          //brchrom  = Util.getBufferedReader("TIER2_IMPUTE_ALL/chr"+chroms[nchrom]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
          }
          else if (ntype == 2)
          {
          brchrom = Util.getBufferedReader("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/"+szcell+"-"+szmark+".gappedPeak.gz");
          //brchrom  = Util.getBufferedReader("AVG_CONVERTED_ROADMAPRAWPVAL/chr"+chroms[nchrom]+"_"+szmark+"_"+szcell+".wig.gz");
          }
          else if (ntype == 3)
          {
          brchrom = Util.getBufferedReader("IMPUTEDPEAKS/"+szcell+"-"+szmark+".imputed.gappedPeak.bed.gz");

          }
          else if (ntype == 4)
          {
          brchrom = Util.getBufferedReader("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/"+szcell+"-"+szmark+".broadPeak.gz");
          }

          //brchrom.readLine();
          //brchrom.readLine();

          while ((szLine = brchrom.readLine())!=null)
          {
          Integer obj = (Integer) hmscore.get(new Double(szLine));

          if (obj == null)
          {
          hmscore.put(new Double(szLine), new Integer(1));
          }
          else
          {
          hmscore.put(new Double(szLine), new Integer(1+obj.intValue()));
          }

          }
          brchrom.close();

          }

          double[] dvalsA = new double[hmscore.size()];
          Iterator itrkey = (Iterator) hmscore.keySet().iterator();
          for (int nk = 0; nk < dvalsA.length; nk++)
          {
          dvalsA[nk] = ((Double) itrkey.next()).doubleValue();
          }
          Arrays.sort(dvalsA);


          HashMap hmbeginrank = new HashMap();
          HashMap hmendrank = new HashMap();

          int nbeginrank = 1;
          for (int ni = 0; ni < dvalsA.length; ni++)
          {
          int ntotalrank = ((Integer) hmscore.get(new Double(dvalsA[ni]))).intValue();

          hmbeginrank.put(new Double(dvalsA[ni]),new Integer(nbeginrank));
          hmendrank.put(new Double(dvalsA[ni]),new Integer(nbeginrank+ntotalrank));
          nbeginrank += ntotalrank;
          }
      //int nbeginrank = ((Integer) hmbeginrank.get(stateA[nwindow] )).intValue();
      //int nendrank = ((Integer) hmendrank.get(stateA[nwindow])).intValue();
      */

      for (int nchrom = 0; nchrom < chroms.length; nchrom++) {
        // System.out.println("\t"+nchrom);
        boolean[] stateA = new boolean[250000000];
        brgwas = new BufferedReader(new FileReader("gwascatalog_sep12_2014.txt"));

        String szcurrchrom = "chr" + chroms[nchrom];

        BufferedReader brchrom = null;
        if (ntype == 0) {
          //		f = new
          // File("/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/"+szcell+"-"+szmark+".narrowPeak.gz");
          brchrom =
              Util.getBufferedReader(
                  "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/"
                      + szcell
                      + "-"
                      + szmark
                      + ".narrowPeak.gz");

          // brchrom  =
          // Util.getBufferedReader("CONVERTED_ROADMAPRAWPVAL/chr"+chroms[nchrom]+"_"+szcell+"-"+szmark+".pval.signal.bigwig.bedgraph.gz.wig.gz");
        } else if (ntype == 1) {
          brchrom =
              Util.getBufferedReader(
                  "IMPUTEDPEAKS/" + szcell + "-" + szmark + ".imputed.narrowPeak.bed.gz");
          //		    brchrom  =
          // Util.getBufferedReader("IMPUTE_ALL/chr"+chroms[nchrom]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
          // brchrom  =
          // Util.getBufferedReader("TIER2_IMPUTE_ALL/chr"+chroms[nchrom]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
        } else if (ntype == 2) {
          brchrom =
              Util.getBufferedReader(
                  "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/"
                      + szcell
                      + "-"
                      + szmark
                      + ".gappedPeak.gz");
          // brchrom  =
          // Util.getBufferedReader("AVG_CONVERTED_ROADMAPRAWPVAL/chr"+chroms[nchrom]+"_"+szmark+"_"+szcell+".wig.gz");
        } else if (ntype == 3) {
          brchrom =
              Util.getBufferedReader(
                  "IMPUTEDPEAKS/" + szcell + "-" + szmark + ".imputed.gappedPeak.bed.gz");

        } else if (ntype == 4) {
          brchrom =
              Util.getBufferedReader(
                  "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/"
                      + szcell
                      + "-"
                      + szmark
                      + ".broadPeak.gz");
        } else if (ntype == 5) {
          brchrom =
              Util.getBufferedReader(
                  "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/"
                      + szcell
                      + "-"
                      + szmark
                      + ".hotspot.all.peaks.bed.gz");
        } else if (ntype == 6) {
          brchrom =
              Util.getBufferedReader(
                  "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/"
                      + szcell
                      + "-"
                      + szmark
                      + ".hotspot.fdr0.01.peaks.bed.gz");
        } else if (ntype == 7) {
          brchrom =
              Util.getBufferedReader(
                  "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/combrep/"
                      + szcell
                      + "-"
                      + szmark
                      + ".macs2.narrowPeak.gz");
        } else if (ntype == 8) {
          brchrom =
              Util.getBufferedReader(
                  "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/"
                      + szcell
                      + "-"
                      + szmark
                      + ".hotspot.broad.bed.gz");
        } else if (ntype == 9) {
          brchrom =
              Util.getBufferedReader(
                  "/broad/compbio/anshul/projects/roadmap/peaks/stdnames30M/broadPeak/"
                      + szcell
                      + "-"
                      + szmark
                      + ".hotspot.fdr0.01.broad.bed.gz");
        }

        /*
           if (ntype == 0)
           {
           brchrom  = Util.getBufferedReader("CONVERTED_ROADMAPRAWPVAL/"+szcurrchrom+"_"+szcell+"-"+szmark+".pval.signal.bigwig.bedgraph.gz.wig.gz");
           }
           else if (ntype ==1)
           {
           brchrom  = Util.getBufferedReader("IMPUTE_ALL/chr"+chroms[nchrom]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
        //     brchrom  = Util.getBufferedReader("TIER2_IMPUTE_ALL/chr"+chroms[nchrom]+"_impute_"+szcell+"_"+szmark+"_1_100000_20.wig.gz");
           }
           else
           {
           brchrom  = Util.getBufferedReader("AVG_CONVERTED_ROADMAPRAWPVAL/chr"+chroms[nchrom]+"_"+szmark+"_"+szcell+".wig.gz");
           }

           brchrom.readLine();
           brchrom.readLine();
           */
        //	    BufferedReader brstate  =
        // Util.getBufferedReader("TIER2_REORDER_MODEL/states25/"+szcell+"_25_imputed12marks_mnemonics.bed.gz");
        //   BufferedReader brstate  =
        // Util.getBufferedReader("/broad/compbio/anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/n25/"+szcell+"_25_coreMarks_segments.bed");

        ///
        // broad/compbio/anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/n15/reordered/"+szcell+"_15_coreMarks_segments.bed
        //       BufferedReader brstate  =
        // Util.getBufferedReader("/broad/compbio/anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/n15/reordered/"+szcell+"_15_coreMarks_segments.bed");
        String szLineA;

        // int nk = 0;
        while ((szLineA = brchrom.readLine()) != null) {
          StringTokenizer ststate = new StringTokenizer(szLineA, "\t");

          String szchrom = ststate.nextToken();
          if (szchrom.equals(szcurrchrom)) {
            int nbegin = Integer.parseInt(ststate.nextToken());
            int nend = Integer.parseInt(ststate.nextToken());
            for (int nk = nbegin; nk < nend; nk++) {
              stateA[nk] = true;
            }
          }

          // stateA[nk] = Double.parseDouble(szLineA);
          // nk++;

          /*

          //String szchrom = ststate.nextToken();
          if (szchrom.equals(szcurrchrom))
          {
          int nbegin = Integer.parseInt(ststate.nextToken())/WINDOW;
          int nend = Integer.parseInt(ststate.nextToken())/WINDOW;

          StringTokenizer stu = new StringTokenizer(ststate.nextToken(),"_");

          //int nstate = Integer.parseInt(stu.nextToken());
          int nstate = Integer.parseInt(stu.nextToken().substring(1));

          //if ((nstate == 14)||(nstate == 15))
          //nstate = 13;
          for (int nk = nbegin; nk < nend; nk++)
          {
          stateA[nk] = nstate;
          }
          }
          */
        }
        brchrom.close();

        // HashSet hsvisited = new HashSet();
        // HashMap hmvisited = new HashMap();

        brgwas.readLine();
        while ((szLine = brgwas.readLine()) != null) {
          StringTokenizer st = new StringTokenizer(szLine, "\t", true);
          st.nextToken();
          st.nextToken(); // tab

          String szchrom = st.nextToken();
          st.nextToken(); // tab

          if (szchrom.equals(szcurrchrom)) {
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
            // System.out.println(szpval);

            // ArrayList alvisit = (ArrayList) hmvisited.get(szpubmedid+"\t"+sztrait);
            // if (alvisit == null)
            // {
            //	alvisit = new ArrayList();
            //	hmvisited.put(szpubmedid+"\t"+sztrait,alvisit);
            // }

            int nbegin = Integer.parseInt(szbegin);
            double dpval = Double.parseDouble(szpval);
            algwasrecs.add(
                new Rec(
                    szpubmedid + "\t" + sztrait, szchrom, nbegin, dpval, theRandom.nextDouble()));

            if (!hscategory.contains(szpubmedid + "\t" + sztrait)) {
              alcategory.add(szpubmedid + "\t" + sztrait);

              hscategory.add(szpubmedid + "\t" + sztrait);
              hmcategory.put(szpubmedid + "\t" + sztrait, new Integer(ncategory));
              ncategory++;
            }

            // alvisit.add(new Rec(nbegin, dpval));
          }

          /*

          if (!hsvisited.contains(szbegin+"\t"+szpubmedid+"\t"+sztrait))
          {
          ncategory = ((Integer) hmcategory.get(szpubmedid+"\t"+sztrait)).intValue();
          int nwindow = Integer.parseInt(szbegin)/WINDOW;
          gwasbackgroundcounts[stateA[nwindow]]++;
          gwasforegroundcounts[ncategory][stateA[nwindow]]++;

          hsvisited.add(szbegin+"\t"+szpubmedid+"\t"+sztrait, new Rec(szbegin, dpval));
          }
          */
        }
        // }
        brgwas.close();

        Rec[] theGWASRecA = new Rec[algwasrecs.size()];

        for (int na = 0; na < theGWASRecA.length; na++) {
          theGWASRecA[na] = (Rec) algwasrecs.get(na);
        }
        Arrays.sort(theGWASRecA, new RecCompare());

        HashMap hmbyID = new HashMap();

        ArrayList algwasrecsCLUMPED = new ArrayList();
        for (int nc = 0; nc < theGWASRecA.length; nc++) {
          ArrayList al = (ArrayList) hmbyID.get(theGWASRecA[nc].szid);

          if (al == null) {
            al = new ArrayList();
            hmbyID.put(theGWASRecA[nc].szid, al);
            al.add(theGWASRecA[nc]);
            algwasrecsCLUMPED.add(theGWASRecA[nc]);
          } else {
            boolean bok = true;
            for (int nk = 0; nk < al.size(); nk++) {
              Rec currrec = (Rec) al.get(nk);
              if ((theGWASRecA[nc].szchrom.equals(currrec.szchrom))
                  && (Math.abs(currrec.nbegin - theGWASRecA[nc].nbegin) < DIST)) {
                bok = false;
              }
            }

            if (bok) {
              al.add(theGWASRecA[nc]);
              algwasrecsCLUMPED.add(theGWASRecA[nc]);
            }
          }
        }

        Rec[] algwasrecsCLUMPEDA = new Rec[algwasrecsCLUMPED.size()];
        for (int na = 0; na < algwasrecsCLUMPEDA.length; na++) {
          Rec currRec = (Rec) algwasrecsCLUMPED.get(na);
          algwasrecsCLUMPEDA[na] =
              new Rec(
                  currRec.szid, currRec.szchrom, currRec.nbegin, currRec.dpval, currRec.drandom);
        }

        Arrays.sort(algwasrecsCLUMPEDA, new RandomRecCompare());
        for (int na = 0; na < algwasrecsCLUMPEDA.length; na++) {
          algwasrecsCLUMPEDA[na].szid = ((Rec) algwasrecsCLUMPED.get(na)).szid;
        }

        HashSet hsvisited = new HashSet();
        HashMap hmvisited = new HashMap();

        // algwasrecsCLUMPEDA[na].
        for (int ne = 0; ne < algwasrecsCLUMPEDA.length; ne++) {
          if (algwasrecsCLUMPEDA[ne].szchrom.equals(szcurrchrom)) {

            ArrayList alvisit = (ArrayList) hmvisited.get(algwasrecsCLUMPEDA[ne].szid);
            if (alvisit == null) {
              alvisit = new ArrayList();
              hmvisited.put(algwasrecsCLUMPEDA[ne].szid, alvisit);
            }

            // int nbegin = Integer.parseInt(algwasrecsCLUMPEDA[na].szbegin);
            // double dpval = Double.parseDouble(szpval);
            alvisit.add(algwasrecsCLUMPEDA[ne]);
            // new Rec(algwasrecsCLUMPEDA[ne].nbegin, algwasrecsCLUMPEDA[ne].dpval));// dpval));
          }
        }

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
              int nwindow = therecA[na].nbegin; // /WINDOW;
              //	   if ((ncell==116)&&(szkey.contains("19838193")))
              // System.out.println(szcurrchrom+"\t"+therecA[na].nbegin+"\t"+stateA[nwindow]);
              // int nbeginrankval = ((Integer) hmbeginrank.get(stateA[nwindow])).intValue();
              // int nendrank = ((Integer) hmendrank.get(stateA[nwindow])).intValue();
              // int ndiff = nendrank-nbeginrankval+1;
              // int nrank = nbeginrankval +((int) (theRandom.nextDouble()*ndiff));

              if (stateA[nwindow]) {
                ngwasbackgroundcountsHIT++;
                ngwasforegroundcountsHIT[ncategory]++;
              }
              ngwasbackgroundcountsALL++;
              ngwasforegroundcountsALL[ncategory]++;

              /*
              if (stateA[nwindow])
              {
              gwasbackgroundcounts.add(new Integer(1));//new Integer(nrank));//[getrank(stateA[nwindow],hmscore); //++;
              gwasforegroundcounts[ncategory].add(new Integer(1));// new Integer(nrank));
              }
              else
              {
              gwasbackgroundcounts.add(new Integer(0));//new Integer(nrank));//[getrank(stateA[nwindow],hmscore); //++;
              gwasforegroundcounts[ncategory].add(new Integer(0));// new Integer(nrank));
              }
              */

              // System.out.println(ncategory+"\t"+gwasforegroundcounts[ncategory].size());

            }
          }
        }

        //      for (int ni = 0; ni <gwasforegroundcounts.length; ni++)
        //  System.out.println(ni+"\t"+gwasforegroundcounts[ni].size());
      }

      for (int ni = 0; ni < ngwasforegroundcountsHIT.length; ni++) {
        if (ngwasforegroundcountsHIT[ni] >= 1) {
          //   double dhyper =
          // hypergeometrictail(gwasforegroundcounts[ni][nj]-1,gwasbackgroundcounts[nj] ,
          // nsumbg-gwasbackgroundcounts[nj] ,nsumfg);
          // * Returns the probability of seeing more than x  objects of type A, when there are nA
          // objects of type A
          //    * nB objects of type B, and nm objects total drawn
          // * This can be used to compute a more accurate p-values than a 1-cumulative probability
          // calculation

          double dpval =
              StateMWGwasPeakHyper.hypergeometrictail(
                  ngwasforegroundcountsHIT[ni] - 1,
                  ngwasbackgroundcountsHIT,
                  ngwasbackgroundcountsALL - ngwasbackgroundcountsHIT,
                  ngwasforegroundcountsALL[ni]);
          double dfold =
              (ngwasforegroundcountsHIT[ni] / (double) ngwasforegroundcountsALL[ni])
                  / (ngwasbackgroundcountsHIT / (double) ngwasbackgroundcountsALL);
          // System.out.println(dpval+"\t"+szcell+"\t"+alcategory.get(ni)+"\t"+ngwasforegroundcountsHIT[ni]+"\t"+ngwasforegroundcountsALL[ni]+"\t"+dfold);
          pw.println(
              dpval
                  + "\t"
                  + szcell
                  + "\t"
                  + alcategory.get(ni)
                  + "\t"
                  + ngwasforegroundcountsHIT[ni]
                  + "\t"
                  + ngwasforegroundcountsALL[ni]
                  + "\t"
                  + dfold);
        }
      }

      /*
             double dsumbg = 0;
             int ncountbg = 0;
             for (int nj = 0; nj < gwasbackgroundcounts.size(); nj++)
             {
             dsumbg += ((Integer) gwasbackgroundcounts.get(nj)).intValue();
             ncountbg++;
             }

             for (int ni = 0; ni < gwasforegroundcounts.length; ni++)
             {
          //		System.out.println("==="+ni+"\t"+gwasforegroundcounts[ni].size());
          double dsumfg = 0;
          int ncountfg = 0;
          for (int nj = 0; nj < gwasforegroundcounts[ni].size(); nj++)
          {
          dsumfg += ((Integer) gwasforegroundcounts[ni].get(nj)).intValue();//[nj];
          ncountfg++;
          }


          if (gwasforegroundcounts[ni].size() >= 1)
          {
          //   double davgbg= dsumbg/(double) ncountbg;
          //double davgfg = dsumfg/(double) ncountfg;

          for (int nj = 0; nj < gwasforegroundcounts[ni].size(); nj++)
          {
          gwasbackgroundcounts.remove(gwasforegroundcounts[ni].get(nj));

          }
          //double[] d1 = new

          double[] backA = new double[gwasbackgroundcounts.size()];
          double dsumbackA = 0;


          ArrayList alrecbackfore = new ArrayList();
          for (int nk = 0; nk < backA.length; nk++)
          {
          backA[nk] = ((Double) gwasbackgroundcounts.get(nk)).doubleValue();
          dsumbackA += backA[nk];
          alrecbackfore.add(new RecBackFore(backA[nk],false));
          }

          double dsumforeA = 0;
          double[] foreA = new double[gwasforegroundcounts[ni].size()];
          for (int nk = 0; nk < foreA.length; nk++)
          {
          foreA[nk] = ((Double) gwasforegroundcounts[ni].get(nk)).doubleValue();
          dsumforeA += foreA[nk];
          alrecbackfore.add(new RecBackFore(foreA[nk], true));
          }


          RecBackFore[] theRecBackFore = new RecBackFore[alrecbackfore.size()];
          for (int na = 0; na < theRecBackFore.length; na++)
          {
          theRecBackFore[na] = (RecBackFore) alrecbackfore.get(na);
          }

          Arrays.sort(theRecBackFore, new RecBackForeCompare());

          int na = 0;
          double dsumforerank =0;
          double dsumbackrank = 0;
          int nb;
          while (na < theRecBackFore.length)
          {
          nb = na +1;

          while ((nb < theRecBackFore.length) && (theRecBackFore[na].dval == theRecBackFore[nb].dval))
          {
              nb++;
          }

          double dcurrrank = (na+nb-1)/2.0;
          for (int nk = na; nk < nb; nk++)
          {
              if (theRecBackFore[nk].bfore)
              {
                  dsumforerank += dcurrrank;
              }
              else
              {
                  dsumbackrank += dcurrrank;
              }
          }

          na = nb;
          }

      // double davgbg = dsumbackA/backA.length;
      //double davgfg = dsumforeA/foreA.length;

      double davgbg = dsumbackrank/backA.length;
      double davgfg = dsumforerank/foreA.length;

      // double dpval = theMannWhitneyUTest.mannWhitneyUTest(backA,foreA)/2;
      double dpval = theMannWhitneyUTest.mannWhitneyUTest(backA,foreA);
      */

      /*

      if ((szcell.equals("E124"))&&(((String) alcategory.get(ni)).contains("23740937")))
      {
      for (int nd = 0; nd < backA.length;nd++)
      {
      System.out.println("back\t"+nd+"\t"+backA[nd]);
      }

      for (int nd = 0; nd < foreA.length;nd++)
      {
      System.out.println("fore\t"+nd+"\t"+foreA[nd]);
      }

      System.out.println(szcell+"\t"+alcategory.get(ni)+"\t"+davgbg+"\t"+davgfg+"\t"+backA.length+"\t"+foreA.length+"\t"+dpval);
      }
      */

      /*
         String szbetter;
         if (davgfg < davgbg)
         szbetter = "1";
         else
         szbetter = "0";
      // if ((dpval < 0.01)&&(davgfg> davgbg))
      //if (dpval < 0.01)
      {

      pw.println(dpval+"\t"+szbetter+"\t"+szcell+"\t"+alcategory.get(ni)+"\t"+nf.format(davgbg)+"\t"+nf.format(davgfg)+"\t"+backA.length+"\t"+foreA.length+"\t");
      //			           dpval+"\t"+szbetter);

      //pw.println(szcell+"\t"+alcategory.get(ni)+"\t"+davgbg+"\t"+davgfg+"\t"+backA.length+"\t"+foreA.length+"\t"+dpval);
      //System.out.println(szcell+"\t"+alcategory.get(ni)+"\t"+davgbg+"\t"+davgfg+"\t"+backA.length+"\t"+foreA.length+"\t"+dpval);

      }

      for (int nj = 0; nj < gwasforegroundcounts[ni].size(); nj++)
      {
      gwasbackgroundcounts.add(gwasforegroundcounts[ni].get(nj));

      }
      */

    }

    //	    double dmututalinfo = 0;
    // for (int nj = 1; nj < gwasforegroundcounts[ni].length; nj++)
    // {

    /**
     * Returns the probability of seeing more than x objects of type A, when there are nA objects of
     * type A nB objects of type B, and nm objects total drawn This can be used to compute a more
     * accurate p-values than a 1-cumulative probability calculation
     */
    //    public static double hypergeometrictail(int nx, int nA, int nB, int nm)

    // double dhyper = hypergeometrictail(gwasforegroundcounts[ni][nj]-1,gwasbackgroundcounts[nj] ,
    // nsumbg-gwasbackgroundcounts[nj] ,nsumfg);
    // if (dhyper < 0.001)
    // System.out.println(szcell+"\t"+alcategory.get(ni)+"\t"+nj+"\t"+dhyper+"\t"+gwasforegroundcounts[ni][nj]+"\t"+gwasbackgroundcounts[nj]
    //                   +"\t"+(nsumbg-gwasbackgroundcounts[nj])+"\t"+nsumfg);
    /*
    double dxy0 = gwasforegroundcounts[ni][nj]/dsumbg;
    double dxy1 = (gwasbackgroundcounts[nj]-gwasforegroundcounts[ni][nj])/dsumbg;
    double dx = gwasbackgroundcounts[nj]/dsumbg;
    double dy0 = dsumfg/dsumbg;
    double dy1 =  1-dy0;//gwasbackgroundcounts[nj]/dsumbg;

    if (dxy0 > 0)
    dmututalinfo += (dxy0 *Math.log(dxy0/(dx*dy0)));//+

    if (dxy1 > 0)
    dmututalinfo +=(dxy1 *Math.log(dxy1/(dx*dy1)));
    */
    //   gwasforegroundcounts[ni][nj] /= dsumfg;

    pw.close();
  }
}
