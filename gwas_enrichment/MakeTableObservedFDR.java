import java.io.*;
import java.util.*;

public class MakeTableObservedFDR {

  public static void main(String[] args) throws IOException {
    HashSet hssigstudy = new HashSet();
    HashMap hmimpute = new HashMap();
    HashMap hmobserve = new HashMap();

    ArrayList alcells = new ArrayList();

    String szmark = args[0];
    String sztype = args[1];
    String sztype2 = sztype;
    if (sztype.equals("broad")) sztype2 = "gapped";
    //	BufferedReader br = new BufferedReader(new FileReader("imputed_H3K27ac_2sided.txt"));
    // BufferedReader brimpute = new BufferedReader(new
    // FileReader("imputed_"+szmark+"_"+sztype2+".txt"));//"observed_H3K27ac_2sided.txt"));

    BufferedReader brobserve =
        new BufferedReader(new FileReader("observed_" + szmark + "_" + sztype + ".txt"));
    String szLine;

    ArrayList alpvals = new ArrayList();
    ArrayList alfdrs = new ArrayList();
    //	BufferedReader brfdrlookup  = new BufferedReader(new FileReader("pval_fdr_"+szmark+".txt"));
    BufferedReader brfdrlookup =
        new BufferedReader(new FileReader("pval_fdr_" + szmark + "_" + sztype + ".txt"));

    while ((szLine = brfdrlookup.readLine()) != null) {
      StringTokenizer st = new StringTokenizer(szLine, "\t");
      st.nextToken();
      alpvals.add(new Double(st.nextToken()));
      alfdrs.add(new Double(st.nextToken()));
    }
    brfdrlookup.close();

    BufferedReader brcells = new BufferedReader(new FileReader("cellorder.txt"));
    while ((szLine = brcells.readLine()) != null) {
      alcells.add(szLine);
    }
    brcells.close();

    // PrintWriter pwimpute = new PrintWriter(new
    // FileWriter("table_imputed_"+szmark+"_"+sztype2+".txt"));
    PrintWriter pwobserve =
        new PrintWriter(new FileWriter("fdr_table_observe_" + szmark + "_" + sztype + ".txt"));

    HashSet hscells = new HashSet();

    //	PrintWriter pwimputenh = new PrintWriter(new
    // FileWriter("noheader_table_imputed_"+szmark+"_2sided.txt"));
    // PrintWriter pwobservenh = new PrintWriter(new
    // FileWriter("noheader_table_observe_"+szmark+"_2sided.txt"));

    /*
    while ((szLine =brimpute.readLine())!=null)
    {
        //0.06432427754895576     1       E001    23563607        Height  7388.11 6468.98 14694   74
        StringTokenizer st  =new StringTokenizer(szLine,"\t");

        double dpval = Double.parseDouble(st.nextToken());
        //String szbetter = st.nextToken();
        String szcell = st.nextToken();
        String szid = st.nextToken()+"\t"+st.nextToken();

        if (dpval <0.001)//&&(szbetter.equals("1")))
    	hssigstudy.add(szid);



        //if (szbetter.equals("1"))
        	hmimpute.put(szcell+"\t"+szid, new Double(-Math.log(dpval)/Math.log(10)));
    	//else
    	//hmimpute.put(szcell+"\t"+szid, new Double(Math.log(dpval)/Math.log(10)));
    }
    brimpute.close();
    */

    while ((szLine = brobserve.readLine()) != null) {
      // 0.06432427754895576     1       E001    23563607        Height  7388.11 6468.98 14694   74
      StringTokenizer st = new StringTokenizer(szLine, "\t");

      double dpval = Double.parseDouble(st.nextToken());
      // String szbetter = st.nextToken();
      String szcell = st.nextToken();
      String szid = st.nextToken() + "\t" + st.nextToken();

      hscells.add(szcell);
      if (dpval < 0.001) // &&(szbetter.equals("1")))
      hssigstudy.add(szid);

      // if (szbetter.equals("1"))
      hmobserve.put(szcell + "\t" + szid, new Double(-Math.log(dpval) / Math.log(10)));
      // else
      // hmobserve.put(szcell+"\t"+szid, new Double(Math.log(dpval)/Math.log(10)));
    }
    brobserve.close();

    int NUMCELL = 129;

    Iterator itr = hssigstudy.iterator();
    pwobserve.print("pmid\ttrait");
    // pwimpute.print("pmid\ttrait");

    for (int ncellindex = 0; ncellindex < alcells.size(); ncellindex++) {
      String szcell = (String) alcells.get(ncellindex);

      // pwimpute.print("\t"+szcell);
      pwobserve.print("\t" + szcell);
    }
    // pwimpute.println();

    pwobserve.print("\tpval\tfdr");
    pwobserve.println();

    while (itr.hasNext()) {

      String szid = (String) itr.next();
      pwobserve.print(szid);
      //  pwimpute.print(szid);

      double dmaxlogpval = -1;

      for (int ncellindex = 0; ncellindex < alcells.size(); ncellindex++) {
        String szcell = (String) alcells.get(ncellindex);

        if (hmobserve.get(szcell + "\t" + szid) == null) {
          // if (ncellindex > 0)
          // pwobservenh.print("\tn");
          if (hscells.contains(szcell)) pwobserve.print("\t0");
          else pwobserve.print("\tn");
        } else {
          // pwobservenh.print("\t"+((Double)hmobserve.get(szcell+"\t"+szid)).doubleValue());
          double dval = ((Double) hmobserve.get(szcell + "\t" + szid)).doubleValue();
          pwobserve.print("\t" + dval);

          if (dval >= dmaxlogpval) {
            dmaxlogpval = dval;
          }
        }

        /*
                      if (hmimpute.get(szcell+"\t"+szid)==null)
            {
                          pwimpute.print("\tn");
            //if (ncellindex >0)
                          //pwimputenh.print("\tn");

            }
                      else
            {
                          pwimpute.print("\t"+((Double)hmimpute.get(szcell+"\t"+szid)).doubleValue());
            //  pwimputenh.print("\t"+((Double)hmimpute.get(szcell+"\t"+szid)).doubleValue());

            }
        */

      }

      int nh = 0;
      while (dmaxlogpval < ((Double) alpvals.get(nh)).doubleValue()) {
        //	       alpvals.add(new Double(st.nextToken()));
        // alfdrs.add(new Double(st.nextToken()));
        nh++;
      }
      double dfdr = ((Double) alfdrs.get(nh)).doubleValue();
      pwobserve.println("\t" + dmaxlogpval + "\t" + dfdr);
      // pwimpute.println();

      // pwobservenh.println();
      //  pwimputenh.println();

    }

    // pwimputenh.close();
    //        pwobservenh.close();

    // pwimpute.close();
    pwobserve.close();
  }
}
