import java.io.*;
import java.util.*;

public class ConvertPVALtoFDR {

  public static void main(String[] args) throws IOException {
    //	1:29am jernst@node1390 /broad/compbio/jernst/compbio-hp/SIGNALPREDICT6 $ sort -k1,1 -g
    // FULL_GWAS_PEAKHYPER/random_observed_H3K4me1_100_narrow.txt | more
    // 1.5758141804633378E-5   E030    19343178        Height  9       20      5.3013093289689035

    String szLine;

    ArrayList alrandomvals = new ArrayList();
    String szmark = args[0];
    String sztype = args[1];

    int NUMRANDOM = 100;

    for (int nseed = 100; nseed <= 199; nseed++) {
      BufferedReader br =
          new BufferedReader(
              new FileReader(
                  "FULL_GWAS_PEAKHYPER/random_observed_"
                      + szmark
                      + "_"
                      + nseed
                      + "_"
                      + sztype
                      + ".txt"));

      while ((szLine = br.readLine()) != null) {
        StringTokenizer st = new StringTokenizer(szLine, "\t");
        alrandomvals.add(new Double(Double.parseDouble(st.nextToken())));
      }
      br.close();
    }

    double[] alrandomA = new double[alrandomvals.size()];
    for (int nk = 0; nk < alrandomA.length; nk++) {
      alrandomA[nk] = ((Double) alrandomvals.get(nk)).doubleValue();
    }
    Arrays.sort(alrandomA);

    BufferedReader brobserved =
        new BufferedReader(
            new FileReader("FULL_GWAS_PEAKHYPER/observed_" + szmark + "_" + sztype + ".txt"));
    ArrayList alobs = new ArrayList();
    while ((szLine = brobserved.readLine()) != null) {
      StringTokenizer st = new StringTokenizer(szLine, "\t");
      alobs.add(new Double(Double.parseDouble(st.nextToken())));
    }
    brobserved.close();

    double[] alobsA = new double[alobs.size()];
    for (int nk = 0; nk < alobsA.length; nk++) {
      alobsA[nk] = ((Double) alobs.get(nk)).doubleValue();
    }
    Arrays.sort(alobsA);

    int[] countrandombetter = new int[alobsA.length];
    int nrandomindex = 0;
    for (int nk = 0; nk < alobsA.length; nk++) {
      while ((nrandomindex < alrandomA.length) && (alobsA[nk] >= alrandomA[nrandomindex])) {
        nrandomindex++;
      }
      countrandombetter[nk] = nrandomindex;
    }

    double[] fdr = new double[countrandombetter.length];
    for (int nk = 0; nk < fdr.length; nk++) {
      fdr[nk] = (countrandombetter[nk] / (double) NUMRANDOM) / (double) (nk + 1);
    }

    for (int nk = fdr.length - 2; nk >= 0; nk--) {

      fdr[nk] = Math.min(fdr[nk], fdr[nk + 1]);
    }

    for (int na = 0; na < fdr.length; na++) {
      System.out.println(na + "\t" + (-Math.log(alobsA[na]) / Math.log(10)) + "\t" + fdr[na]);
    }
  }
}
