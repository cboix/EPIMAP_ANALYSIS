import java.io.*;
import java.util.*;
import weka.core.*;
import weka.classifiers.functions.*;
import weka.classifiers.*;
import java.text.*;
// From https://github.com/dnaase/Bis-tools/blob/master/recombination_valley_paper/LinkingRM.java //
// Also see http://www.biolchem.ucla.edu/labs/ernst/roadmaplinking/ //

public class LinkingRM
{
    /**
     * Maximum number of base pairs to consider from TSS
     */
    static int OUTERRANGE =1000000;// 125000;

    /**
     * Minimum distance required from TSS
     */
    static int INNERRANGE = 0;//5000;
    static int SMOOTHRANGE = 5000;

    /**
     * The window size being used
     */
    static int WINDOW = 200;

    static int MISSINGVALUE = -1;

    static int NUMGROUP = 8;

    /**
     *
     */
    static String STATEASSIGNDIR = "/broad/compbio/anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/n15/reordered/STATEBYLINE";
	// "POSTERIORS15_100";

    /**
     *
     */
    static String SIGNALDIR = "/broad/compbio/jernst/compbio-hp/SIGNALPREDICT6/CONVERTED_ROADMAPRAWPVAL"; 
	//"BINNEDDATA";

    /**
     *
     */
    static String EXPRESSIONFILE = "57epigenomes.RPKM.pc";
	// "cellheader_refseqexpressnoprobe.txt";

    //  static String EXPRESSIONFILEHEADER = "EG.name.txt";

    /**
     *
     */
    static String OUTPUTDIR = "LOGISTICLINKS";

    /**
     *
     */
    static String GENEANNOTATIONFILE = "Ensembl_v65.Gencode_v10.ENSG.gene_info";
	// "/seq/compbio-hp/jernst/ENCODEPRODUCTION2/refgenes_12_15_09.txt";

    static String[] chroms= {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"};

    //CTCF, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3,H3K9ac,H4K20me1, WCE
    static String[] marks = {"H3K27ac","H3K4me1","H3K4me2","H3K9ac","DNase"};
    //static int[] marksuse = {1,4,5};

    static double RIDGE = 1.0;

    ///////////////////////////////////////////////////////////////////////////////////

    static class Rec
    {
        String szchrom;
        int ntssbin;
        boolean bposstrand;

        Rec(String szchrom, int ntssbin, boolean bposstrand)
        {
            this.szchrom =szchrom;
            this.ntssbin = ntssbin;
            this.bposstrand = bposstrand;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////

    static class RandomRec
    {
        int nindex;
        float drval;
        RandomRec(int nindex, float drval)
        {
            this.nindex = nindex;
            this.drval = drval;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////

    /**
     * 
     */
    static public class RandomRecCompare implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
            RandomRec r1 = (RandomRec) o1;
            RandomRec r2 = (RandomRec) o2;
            if (r1.drval < r2.drval)
	    {
	       return -1;
	    }
            else if (r1.drval > r2.drval)
	    {
	       return 1;
	    }
            else
	    {
	       return 0;
	    }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////

    static class LinkRec
    {
	String szchrom;
	int nbegin;
	int nend;
	String szgene;
	double dprob;
	int ndist;

	public LinkRec(int ndist, double dprob,String szchrom, int nbegin, int nend, String szgene)
	{
	    this.ndist = ndist;
	    this.dprob = dprob;
	    this.szchrom = szchrom;
	    this.nbegin = nbegin;
	    this.nend = nend;
	    this.szgene = szgene;
	}

    }

    //////////////////////////////////////////////////////////////////////////////////
    /**
     *
     */
    static public class LinkRecCompare implements Comparator
    {
        public int compare(Object o1, Object o2)
        {
	    LinkRec r1 = (LinkRec) o1;
	    LinkRec r2 = (LinkRec) o2;

	    int ncompareval = r1.szchrom.compareTo(r2.szchrom);

	    if (ncompareval != 0)
	    {
		return ncompareval;
	    }
	    else if (r1.nbegin < r2.nbegin)
	    {
		return -1;
	    }
	    else if (r1.nbegin > r2.nbegin)
      	    {
		return 1;
	    }
	    else
	    {
		return r1.szgene.compareTo(r2.szgene);
	    }
	}
    }



    /////////////////////////////////////////////////////////////////////////////////

    /**
     * Computes the correlation between two vectors
     */
    public static float correlation(float[] xvalues, float[] yvalues,boolean[] present)
    {
       float  dsumx = 0,
      	      dsumy = 0,
	      dsumxsq = 0,
	      dsumysq = 0,
	      dsumxy = 0,
	      dvarx,
	      dvary,
	      dcoeff;

       int numvalues = 0;

       for (int nindex = 0; nindex < xvalues.length; nindex++)
       {
	   if (present[nindex])
	  {
          dsumx += xvalues[nindex];
	  dsumy += yvalues[nindex];
	  dsumxsq += xvalues[nindex]*xvalues[nindex];
	  dsumysq += yvalues[nindex]*yvalues[nindex];
	  dsumxy  += xvalues[nindex]*yvalues[nindex];
	  numvalues++;
	  }
       }
       // numvalues = xvalues.length;

       if (numvalues==0)
       {
          dcoeff = 0;
       }
       else
       {
          dvarx = dsumxsq - dsumx*dsumx/numvalues;
          dvary = dsumysq - dsumy*dsumy/numvalues;

          if (dvarx*dvary <= 0)
          {
             dcoeff = 0;
	  }
          else
          {
	     dcoeff = (dsumxy - dsumx*dsumy/numvalues)/(float) Math.sqrt(dvarx*dvary);
	  }
       }
       return dcoeff;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    public static void main(String[] args) throws Exception
    {
       String SZCELL = args[0];
       String STATEWANT = args[1];

       NumberFormat nf = NumberFormat.getInstance();
       nf.setMaximumFractionDigits(3);
       double dthreshold = Double.parseDouble(args[2]);

       int numsmooth = SMOOTHRANGE/WINDOW;

       //stores a mapping of all genes to chromosome, TSS bin, and strand
       BufferedReader brgene = new BufferedReader(new FileReader(GENEANNOTATIONFILE));
       PrintWriter pw = new PrintWriter(new FileWriter(OUTPUTDIR+"/links_"+SZCELL+"_"+STATEWANT+"_"+dthreshold+".txt"));

       HashMap hmgene = new HashMap();

       String szLine;
       //int[] marksuse = new int[marksuseNames.length];

       //builds a set of valid chromosomes
       HashSet hsvalidchrom = new HashSet();
       for (int na = 0; na < chroms.length; na++)
       {
          hsvalidchrom.add("chr"+chroms[na]);
       }


       //reads in the set of genes
       //genes must be on a valid chromosome
       //format of file is
       //filler col, gene, chrom, strand, begin, end
       while ((szLine = brgene.readLine())!= null)
       {
          int ntss;
          String[] szA = szLine.split("\\s");
          String szstrand = szA[4];
          String szgene = szA[0];
          String szchrom = "chr"+szA[1];

	  if (hsvalidchrom.contains(szchrom))
	  {
	      //gene is on a valid chromosome
	     String szbegin = szA[2];
	     String szend = szA[3];
	     if (szstrand.equals("1"))
             {
	        ntss  = Integer.parseInt(szbegin);
             }
             else if (szstrand.equals("-1"))
             {
                ntss = Integer.parseInt(szend);
             }
	     else
	     {
		 throw new IllegalArgumentException("Invalid strand of "+szstrand+" given!");
	     }
             int ntssbin = ntss/WINDOW;

	     //takes the most outer most bin on valid chrom
	     //if a gene has multiple transcription start sites uses the one with
	     //the lowest coordinate value if on the positive strand, and the greatest
	     //if on the negative strand
	     Rec theRec = (Rec) hmgene.get(szgene);
	     if (theRec == null)
             {
      	        hmgene.put(szgene,new Rec(szchrom,ntssbin,szstrand.equals("1")));
	     }
	     else if ((szstrand.equals("1"))&&(ntssbin < theRec.ntssbin))
	     {
	        hmgene.put(szgene,new Rec(szchrom,ntssbin,szstrand.equals("1")));
	     }
	     else if ((szstrand.equals("-1"))&&(ntssbin > theRec.ntssbin))
             {
      	        hmgene.put(szgene,new Rec(szchrom,ntssbin,szstrand.equals("1")));
	     }
	  }
       }
       brgene.close();

       BufferedReader brexpress = new BufferedReader(new FileReader(EXPRESSIONFILE));
       // BufferedReader brexpressheader = new BufferedReader(new FileReader(EXPRESSIONFILEHEADER));
   
       ArrayList alcells = new ArrayList();
       //brexpressheader.readLine(); //flush E000
       StringTokenizer stu = new StringTokenizer(brexpress.readLine(),"\t");
       stu.nextToken();
       stu.nextToken();
       while (stu.hasMoreTokens())
       {
	   alcells.add(stu.nextToken());
       }
       // StringTokenizer stheader = new StringTokenizer(brexpress.readLine(),"\t"); //expression header line
       //stheader.nextToken();//get rid of header first column
       //int ncellindex = 0;
       String[] cells = new String[alcells.size()];
       for (int ncellindex = 0; ncellindex < cells.length; ncellindex++)
       {
	   cells[ncellindex] = (String) alcells.get(ncellindex);
       }


       HashSet  hsunique = new HashSet(); //only stores unique TSS gene-expression location
       ArrayList aluniqueExpress  = new ArrayList(); //stores the expression for unique locations
       ArrayList aluniqueID = new ArrayList(); //contains unique genes

       //reads in the expression file       
       //only keeps genes with unique start bin or expression
       //as sometimes multiple genes are associated with the same probe with the same gene start
       while ((szLine = brexpress.readLine())!=null)
       {
          StringTokenizer st = new StringTokenizer(szLine,"\t");
          String szID = st.nextToken();
          String szexpressvals = "";
          float[] fexpressvals = new float[cells.length];
	  int ncell = 0;
	  st.nextToken();//flush E000
          while (st.hasMoreTokens())
	  {
	     String sztoken = st.nextToken();
	     szexpressvals += sztoken+"\t";
	     fexpressvals[ncell] = Float.parseFloat(sztoken);
	     ncell++;
	  }

	  Rec theRec = (Rec) hmgene.get(szID);
	   
	  if ((theRec != null)&&(!hsunique.contains(theRec.szchrom+"\t"+theRec.ntssbin+"\t"+theRec.bposstrand+"\t"+szexpressvals)))
          {
	     aluniqueID.add(szID); //stores ID 
             aluniqueExpress.add(fexpressvals); //stores expression vals
	     
             hsunique.add(theRec.szchrom+"\t"+theRec.ntssbin+"\t"+theRec.bposstrand+"\t"+szexpressvals);
          }
       }
       brexpress.close();



       int numgenes = aluniqueExpress.size();

       /*  
       for (int ncell = 0; ncell < cells.length; ncell++)
       {
          for (int nchrom = 0; nchrom < chroms.length; nchrom++)
          {
	     BufferedReader br = new BufferedReader(new FileReader(SIGNALDIR+"/"+cells[ncell]+"_chr"+chroms[nchrom]+".txt"));
	     String szmarkheader = br.readLine();
	     if ((ncell ==0)&&(nchrom ==0))
	     {
		 String[] headerA = szmarkheader.split("\\s");
		 for (int nj = 0; nj < marksuseNames.length; nj++)
		 {
		     int nk = 0;
		     boolean bfound = false;
		     while ((nk < headerA.length)&&(!bfound))
		     {
		        if (headerA[nk].equals(marksuseNames[nj]))
		        {
			    marksuse[nj] = nk;
			    bfound = true;
			}
			nk++;
		     }
		     if (!bfound)
		     {
			 throw new IllegalArgumentException("Did not find mark "+marksuseNames[nj]);
		     }
		 }
	     }	    
	  }
       }
       */

       HashMap[] hmloctosignal = new HashMap[marks.length];
       for (int nk = 0; nk < hmloctosignal.length; nk++)
       {
          hmloctosignal[nk] = new HashMap();
       }



       HashMap[] hmloctopresent = new HashMap[marks.length];
       for (int nk = 0; nk < hmloctopresent.length; nk++)
       {
	   hmloctopresent[nk] = new HashMap();
       }


       for (int nchrom = 0; nchrom < chroms.length; nchrom++)
       {	
	   System.out.println(chroms[nchrom]);
	   //going to iterate over each chromosome

          String szchrom = "chr"+chroms[nchrom];
          BufferedReader brmax = Util.getBufferedReader(STATEASSIGNDIR+"/"+SZCELL+"_15_coreMarks_chr"+chroms[nchrom]+"_statebyline.txt.gz");
	  brmax.readLine();
	  brmax.readLine();
								   //"/max_"+SZCELL+"_chr"+chroms[nchrom]+".txt"));
	  BufferedReader[][] brsignal = new BufferedReader[marks.length][cells.length];
	  for (int nmark = 0; nmark < marks.length; nmark++)
	  {
             for (int ncell = 0; ncell < cells.length; ncell++)
	     {
		 String szinfile = SIGNALDIR+"/chr"+chroms[nchrom]+"_"+cells[ncell]+"-"+marks[nmark]+".pval.signal.bigwig.bedgraph.gz.wig.gz";
		 File f = new File(szinfile);
		 if (f.exists())
		 {
	         ///broad/compbio/jernst/compbio-hp/SIGNALPREDICT6/CONVERTED_ROADMAPRAWPVAL/chr10_E001-H3K4me3.pval.signal.bigwig.bedgraph.gz.wig.gz
		     brsignal[nmark][ncell] = Util.getBufferedReader(szinfile);
						     //SIGNALDIR+"/chr"+chroms[nchrom]+"_"+cells[ncell]+"-"+marks[nmark]+".pval.signal.bigwig.bedgraph.gz.wig.gz"));
                    brsignal[nmark][ncell].readLine();
                    brsignal[nmark][ncell].readLine();
		 }
		 else
		 {
		     brsignal[nmark][ncell] = null;
		 }
	     }
          }

	  String szLineState;
	  int nline = 0;
	  while ((szLineState = brmax.readLine())!=null)
          {
	     if (szLineState.equals(STATEWANT))
	     {
	        //we found the state we are predicting for
	        float[][] signalvals = new float[marks.length][cells.length];
		boolean[][] bpresent = new boolean[marks.length][cells.length];

	        for (int ncell = 0; ncell < cells.length; ncell++)
	        {
		    for (int nmark = 0; nmark < marks.length; nmark++)
		    {
			//String szLineData = brsignal[ncell].readLine();
		   
			//  	   String[] szvalsA = szLineData.split("\\s");
			//for (int nmark = 0; nmark < marks.length; nmark++)
			//{	
		       //stores in the normalized signal values              
			if (brsignal[nmark][ncell] != null)
			{
			    bpresent[nmark][ncell] = true;
			    //signalvals[nmark][ncell] = MISSINGVALUE;

			   float dsum = 0;
			   int nread = 0;
			   for (int ni = 0; ni < NUMGROUP; ni++)
			   {
			       String szLineGroup = brsignal[nmark][ncell].readLine();
			       if (szLineGroup != null)
			       {
			          dsum += Double.parseDouble(szLineGroup);
				  nread++;
			       }
			       else
			       {
				   break;
			       }
			   }  
			   signalvals[nmark][ncell] = dsum/nread;// Float.parseFloat(szvalsA[marks[nmark]]);
		       //}
			}

		    }
		}

		for (int nmark = 0; nmark < hmloctosignal.length; nmark++)
	        {
	           hmloctosignal[nmark].put(szchrom+"\t"+nline, signalvals[nmark]);
                   hmloctopresent[nmark].put(szchrom+"\t"+nline, bpresent[nmark]);

		}		    		 
	          //will store all signals associated with this
	     }
	     else
	     {
	        for (int nmark = 0; nmark < marks.length; nmark++)
		{
		 //skip over the signal values without storing
	           for (int ncell = 0; ncell < cells.length; ncell++)
	           {
		       for (int ni = 0; ni < NUMGROUP; ni++)
		       {
			   if (brsignal[nmark][ncell] != null)
			   {
			       brsignal[nmark][ncell].readLine();
			   }
		       }
		   }
		}
	     }
	     nline++;
	  }

	  brmax.close();

	  for (int nmark = 0; nmark < marks.length; nmark++)
	  {
	     for (int ncell = 0; ncell < brsignal.length; ncell++)
	     {
	        brsignal[nmark][ncell].close();
	     }
	  }
       } //finished reading through all the chromosomes

       ///////////////////////////////////////////////////////////////////////////////

       ArrayList allinkrec = new ArrayList();

	int nbegin = -OUTERRANGE/WINDOW;
	int nend = OUTERRANGE/WINDOW;
	for (int npos = nbegin; npos <= nend; npos++)
	{
	    //position determines position relative to TSS
	    if (npos % 100 == 0)
	    System.out.println(npos);

	   if (Math.abs(npos) >=INNERRANGE/WINDOW)//&&(Math.abs(npos)<=OUTERRANGE/WINDOW))
	   {
	      //only considering position if suff
	      //going to build a classifier with all positive those within NUMSMOOTH intervals which are positive
	      //negative for those in randomization

              FastVector attributes = new FastVector();
           
              for (int ni = 0; ni < marks.length; ni++)
      	      {
	         Attribute exp = new Attribute("mark"+ni);
                 attributes.addElement(exp);
	      }
	      FastVector classregset = new FastVector();
	      classregset.addElement("0");
	      classregset.addElement("1");
	      Attribute classreg = new Attribute("class",classregset);
	      attributes.addElement(classreg);

	      //int[] randomindex = new int[numgenes];
	      RandomRec[] theRandomRecs = new RandomRec[numgenes];
	      Random theRandom = new Random(4543);


	      //randomizes the ordering of genes
   	      for (int ngene = 0; ngene < numgenes; ngene++)
	      {
                 theRandomRecs[ngene] = new RandomRec(ngene,theRandom.nextFloat());
	      }
	      Arrays.sort(theRandomRecs, new RandomRecCompare());

	      
	      //for (int nj = 0; nj < theRandomRecs.length; nj++)
	      //{
              //   randomindex[nj] = theRandomRecs[nj].nindex;
	      //}	   

	      Classifier theLinkClassifier = new Logistic();
	      ((Logistic) theLinkClassifier).setRidge(RIDGE);
	      Instances theInstances = new Instances("express",attributes,0);
	      theInstances.setClassIndex(theInstances.numAttributes()-1);
	
	      ArrayList alapply = new ArrayList();
	      ArrayList alapplyID = new ArrayList();
	      ArrayList alapplyRecs = new ArrayList();

	      for (int nsmooth = npos - numsmooth; nsmooth <= npos + numsmooth; nsmooth++)
	      {	      
		  //going through all positions within numsmooth to 
		  //smooth the range

                 for (int ngene = 0; ngene < numgenes; ngene++)
                 {
                    Rec theRec = (Rec) hmgene.get((String) aluniqueID.get(ngene));
	            int nshiftpos;

	            if (theRec.bposstrand)
	            {
			//position after shifting for smoothing done in a strand aware way
                       nshiftpos = theRec.ntssbin - nsmooth;
		    }
	            else
	            {
                       nshiftpos = theRec.ntssbin + nsmooth;
		    }

	            float[] dsignalvals  = (float[]) hmloctosignal[0].get(theRec.szchrom+"\t"+nshiftpos);
	            if (dsignalvals != null)
	            {	            		      
			//we are in the desired state at this position


		       float[] dexpressvals = (float[]) aluniqueExpress.get(ngene);
		       float[] drandomexpressvals = (float[]) aluniqueExpress.get(theRandomRecs[ngene].nindex);

		       Instance theInstance = new DenseInstance(marks.length+1);
		       //creates a feature vector for the correlation between the mark and expression
		       theInstance.setDataset(theInstances);
		       //gives the label one to the vector
		       theInstance.setValue(marks.length,"1");
		       theInstance.setWeight(1);

		       if (npos == nsmooth)
		       {
			   //if we are at the center position, then saving this entry to apply later
		          alapply.add(theInstance);
		  	  alapplyID.add(aluniqueID.get(ngene));
			  alapplyRecs.add(theRec);
		       }

		       //creating the random instance vector
		       Instance theRandomInstance;
		       theRandomInstance = new DenseInstance(marks.length+1);
		       theRandomInstance.setDataset(theInstances);
		       theRandomInstance.setValue(marks.length,"0");
		       theRandomInstance.setWeight(1);		       
		     
	               for (int nmark = 0; nmark < marks.length; nmark++)
	               {
			   //going through each mark

			   //computing the correlation with gene expression based on the real signal
		          dsignalvals  = (float[]) hmloctosignal[nmark].get(theRec.szchrom+"\t"+nshiftpos);
			  boolean[] bpresent = (boolean[]) hmloctopresent[nmark].get(theRec.szchrom+"\t"+nshiftpos);

			  float dcorr = correlation(dsignalvals,dexpressvals,bpresent);			
			  theInstance.setValue(nmark, new Float(dcorr));
			  //System.out.println("observed\t"+dcorr);

			  //computing the correlation with gene expression based on the random signal
			  dcorr = correlation(dsignalvals,drandomexpressvals,bpresent);
			  theRandomInstance.setValue(nmark, new Float(dcorr));			
			  //System.out.println("random\t"+dcorr);
		       }

		       //adds in the training instances
		       theInstances.add(theInstance);    	     
		       theInstances.add(theRandomInstance);		     
		    }
		 }
	      }

	      //training the classifier
	      theLinkClassifier.buildClassifier(theInstances);

	      //System.out.println("***"+theInstances.size()+"****"+theLinkClassifier);
	      //applying the classifier
	      for (int ni = 0; ni < alapply.size(); ni++)
	      {
	         String szchrom;
	         int ntssbin;
	         boolean bposstrand;

	         //applies the classifier to the training data
	         double[] vals =theLinkClassifier.distributionForInstance((Instance) alapply.get(ni));
	         Rec theRec = (Rec) alapplyRecs.get(ni);
		 //System.out.println(vals[1]);
	         if (vals[1] >=dthreshold/(1.0+dthreshold))
	         {   //x/(1-x)=k --> (k(1-x))=x --> k=x+kx --> k/(1+k)
	            if (theRec.bposstrand)
	            {
			allinkrec.add(new LinkRec(-npos*WINDOW,vals[1],theRec.szchrom, (theRec.ntssbin-npos)*WINDOW, (theRec.ntssbin-npos+1)*WINDOW,(String) alapplyID.get(ni)));
	            }
	            else
	            {
			allinkrec.add(new LinkRec(-npos*WINDOW,vals[1],theRec.szchrom, (theRec.ntssbin+npos)*WINDOW, (theRec.ntssbin+npos+1)*WINDOW,(String) alapplyID.get(ni)));
		    }
		 }
	      }
	   }
	}


	LinkRec[] allinkrecA = new LinkRec[allinkrec.size()];
	for (int nk = 0; nk < allinkrecA.length; nk++)
	{
	    allinkrecA[nk] = (LinkRec) allinkrec.get(nk);
	}
	Arrays.sort(allinkrecA, new LinkRecCompare());

	String szprev = "";
	for (int nk = 0; nk < allinkrecA.length; nk++)
	{
	    LinkRec currLinkRec = allinkrecA[nk];
	    double dratio = currLinkRec.dprob/(1-currLinkRec.dprob);
	    String szcurr = currLinkRec.szchrom+"\t"+currLinkRec.nbegin+"\t"+currLinkRec.nend+"\t"+currLinkRec.szgene+"\t"+nf.format(dratio)+"\t"+currLinkRec.ndist;

	    if (!szprev.equals(szcurr))
	    {
	       pw.println(szcurr);
	    }
	    szprev= szcurr;
	}
       
	pw.close();
    }
}


