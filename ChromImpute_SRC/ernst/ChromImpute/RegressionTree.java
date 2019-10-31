package ernst.ChromImpute;

import java.io.*;
import java.util.*;
import java.text.*;
    
public class RegressionTree
{

    /**
     * array list holding an array of pointers for xvaluescount at various depths
     */
    ArrayList xvaluescountAL;

    /**
     * array list holding an array of pointers for y sums at various depths
     */
    ArrayList yvaluessumAL;

    /**
     * array list holding an array of pointers to arrays of the indicies of locations assigned to the right
     */
    ArrayList rightlocationsAL;

    /**
     * stores the number of positions in the array actually being used at each depth
     */
    ArrayList rightlocationsSizeAL;

    /**
     * Pointer to the root of the regression tree
     */
    TreeNode theTree;

    /**
     * For each feature maps x-values to an index
     */
    HashMap[] hmdatatoindex;

    /**
     * Stores for each feature the set of x-values observed for it sorted
     */
    float[][] xvals;

    /**
     * Random number generator used to break ties on selecting features 
     */
    Random theRandom;

    /**
     * Stores the data to use in training, the data is organized features then instance values
     */
    float[][] data;

    /**
     * Stores the target output values for each data instance
     */
    float[] output;

    /**
     * A numberformat object for outputting the contents of the tree
     */
    NumberFormat nf;

    /**
     * Minimum number of locations
     */
    int nminnumlocations;


    ////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Record for a node in regression tree
     * Has pointers to the left and right node, an index of a feature to split on,
     * a split value, and an average value associated with the elements of the node
     */
    static class TreeNode
    {
	TreeNode left;
	TreeNode right;
	int nsplitfeatureindex;
	double dsplitval;
	double dmean;
    }

    ////////////////////////////////////////////////////////////////////////////////////////

    static class Rec
    {
	float dfeatureval;
	float doutval;

	Rec(float dfeatureval, float doutval)
	{
	    this.dfeatureval = dfeatureval;
	    this.doutval = doutval;
	}
    }


    ////////////////////////////////////////////////////////////////////////////////////////


    /**
     *
     */
    public RegressionTree(ArrayList dataAL, ArrayList outAL, int nminnumlocations)
    {
	//initializes the minimum num locations
       this.nminnumlocations = nminnumlocations;

       //stores the total number of instances in the training data
       int numinstances = dataAL.size();

       //stores the total number of features in the training data
       int numfeatures = ((float[]) dataAL.get(0)).length;

       //stores the contents of dataL and outAL into data and output
       output = new float[numinstances];
       data = new float[numfeatures][numinstances]; // the first dimension of data is the feature then the instance
       for (int ni = 0; ni < numinstances; ni++)
       {
          float[] fvals = (float[]) dataAL.get(ni);
          for (int nj = 0; nj < fvals.length; nj++)
          {
             data[nj][ni] = fvals[nj];
          }
	  output[ni] = ((Float) outAL.get(ni)).floatValue();
	}

        //inititalizes the output numberformat
	nf = NumberFormat.getInstance(Locale.ENGLISH);
	nf.setMaximumFractionDigits(2);
	nf.setGroupingUsed(false);

	//inititalizes the Randomization for breaking ties on feature to split on
	this.theRandom = new Random(456);

	//a vector of the indicies of the data associated with a node of the tree
	int[] datalocations = new int[data[0].length];

	//stores for each feature how often it occured
	int[][] xvaluescount  = new int[data.length][];

	//stores for each feature the sum of the y-values corresponding to it
	double[][] yvaluessum = new double[data.length][];

	//initializes the tree 
	initTree(datalocations, xvaluescount, yvaluessum);

	//initialize arraylist containing array pointers of right locations
        rightlocationsAL = new ArrayList();

	//inititalize arraylist with size of corresponding right location
        rightlocationsSizeAL = new ArrayList();

	//inititalize arraylist with x-value count storage
	xvaluescountAL = new ArrayList();

	//initialize arraylist with y-value sum  storage
	yvaluessumAL = new ArrayList();

	buildTree(theTree,datalocations,datalocations.length,xvaluescount,yvaluessum,0);

	//helping the garbage collection now to clear out the memory allocated in storing
	//data to build the regression tree previously stored in class variables
	xvaluescountAL = null;
	yvaluessumAL = null;
	this.data =null;
	this.output = null;
        hmdatatoindex = null;
        xvals = null;
	rightlocationsAL = null;
        rightlocationsSizeAL = null;
    }


    ////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Constructor that builds a RegressionTree by reading contents of brtreefile
     */ 
    public RegressionTree(BufferedReader brtreefile) throws IOException
    {
	nf = NumberFormat.getInstance(Locale.ENGLISH);
	nf.setMaximumFractionDigits(2);
	nf.setGroupingUsed(false);
	theTree = new TreeNode();
	loadTree(brtreefile, theTree);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * This procedure loads a regression tree stored in brtreefile
     */
    public void loadTree(BufferedReader brtreefile, TreeNode theNode) throws IOException
    {
	StringTokenizer st = new StringTokenizer(brtreefile.readLine(),"\t");
	String sztoken1 = st.nextToken();
	if (sztoken1.equals("n"))
	{
	    //this is a leaf node reads the predicted value associated with it
	    theNode.dmean = Double.parseDouble(st.nextToken());
	}
	else
	{
	    //loads the index of the feature to split on
	    theNode.nsplitfeatureindex = Integer.parseInt(sztoken1);

	    //loads the value of the split feature
	    theNode.dsplitval = Double.parseDouble(st.nextToken());

	    //creates a node to the left
	    theNode.left = new TreeNode();
	    loadTree(brtreefile, theNode.left);

	    //creates a node to the right
	    theNode.right = new TreeNode();
	    loadTree(brtreefile, theNode.right);
	}
    }
 
    /////////////////////////////////////////////////////////////////////////

    /**
     * Compares records on the feature values
     */    
    static class RecCompare implements Comparator, Serializable
    {
	public int compare(Object o1, Object o2)
	{
	    Rec r1 = (Rec) o1;
	    Rec r2 = (Rec) o2;

	    if (r1.dfeatureval < r2.dfeatureval)
	    {
		return -1;
	    }
	    else if (r1.dfeatureval > r2.dfeatureval)
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	}
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * Compares records based on the output feature value
     */
    static class OutCompare implements Comparator, Serializable
    {
	public int compare(Object o1, Object o2)
	{
	    Rec r1 = (Rec) o1;
	    Rec r2 = (Rec) o2;

	    if (r1.doutval < r2.doutval)
	    {
		return -1;
	    }
	    else if (r1.doutval > r2.doutval)
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	}
    }


    ////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Record for y-tally
     */
    static class RecYtally
    {
       int ntally;
       double dysum;

       RecYtally(int ntally, double dysum)
       {
	  this.ntally = ntally;
	  this.dysum = dysum;
       }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////


    /**
     * Stores information on the best split found so far
     */
    static class BestSplit
    {
	int nbestfeatureindex; //index of best feature to split on
	double dbestdataval; //value of the best split feature
	int nbestdataindex;  //last index of the left side if vals sorted
	double dbestleftmean; //best split feature average to the left
	double dbestrightmean; //best split feature average to the right
	//double dbestSS; //bestSS for split

        BestSplit(int nbestfeatureindex,int nbestdataindex, double dbestdataval,double dbestleftmean,double dbestrightmean)//, double dbestSS)
	{
	    this.nbestfeatureindex = nbestfeatureindex;
	    this.nbestdataindex = nbestdataindex;
	    this.dbestdataval = dbestdataval;
	    this.dbestleftmean = dbestleftmean;
	    this.dbestrightmean = dbestrightmean;
	    //this.dbestSS = dbestSS;
	}
    }


    ////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Inititalizes the regression tree that will then be built
     */
    public void initTree(int[] datalocations, int[][] xvaluescount, double[][] yvaluessum)
    {

	//Initializes the root of the tree
	theTree = new TreeNode();

	//will store for each feature unique x-values in sorted order
	xvals = new float[data.length][];

	for (int nindex = 0; nindex < datalocations.length; nindex++)
	{
	    //initializes the root of the tree to be associated with all data points
	    datalocations[nindex] = nindex;
	} 

	//array of HashMap from data values to index
	hmdatatoindex = new HashMap[data.length];

	for (int nfeature = 0; nfeature < data.length; nfeature++)
	{
	    //going through each feature
	   hmdatatoindex[nfeature] = new HashMap();

	   //stores a map of an x-feature value to an index
	   HashMap hmdatatoindex_nfeature = hmdatatoindex[nfeature];

	   //stores for each x-feature value a record with the tally of the number of times that
	   //feature occured and the total sum of the y-values when that x-feature occured
	   HashMap hmtallysum = new HashMap();

	   //gets a pointer to all array values for the feature
	   float[] data_nfeature = data[nfeature];

	   //going through all data locations
	   for (int ndataindex = 0; ndataindex < datalocations.length; ndataindex++)
	   {
	       //converts the x-feature value to an object
	       //Float objd = new Float(data_nfeature[ndataindex]);
	      Float objd = Float.valueOf(data_nfeature[ndataindex]);

	      //gets the tally associated with this x-feature
	      RecYtally theRec = (RecYtally) hmtallysum.get(objd);
	    
	      if (theRec == null)
	      {
		  //first time we've encountered this data value
		  hmtallysum.put(objd, new RecYtally(1, output[ndataindex]));
	      }
	      else
	      {
		  //we've encountered this data value before incrementing count
	         hmtallysum.put(objd, new RecYtally(theRec.ntally+1,output[ndataindex]+theRec.dysum));
	      }
	   }

	   //stores all the unique x-values in hmtally sum into xvals_nfeature and sorts them
	   xvals[nfeature] = new float[hmtallysum.size()];
           float[] xvals_nfeature = xvals[nfeature];
	   Iterator hmtallyxitr = hmtallysum.keySet().iterator();

	   for (int ni = 0; ni < xvals_nfeature.length; ni++)
	   {
	      Float objd = (Float) hmtallyxitr.next();
	      xvals_nfeature[ni] = objd.floatValue();
	   }
	   Arrays.sort(xvals_nfeature);

	   //stores into xvaluescount and yvaluessum the corresponding occurence tally and y-sum for each xfeature position
	   xvaluescount[nfeature] = new int[xvals_nfeature.length];
	   yvaluessum[nfeature] = new double[xvals_nfeature.length];

	   int[] xvaluescount_nfeature = xvaluescount[nfeature];
	   double[] yvaluessum_nfeature = yvaluessum[nfeature];

	   for (int ni = 0; ni < xvals_nfeature.length; ni++)
	   {
	       //Float objd = new Float(xvals_nfeature[ni]);
	      Float objd = Float.valueOf(xvals_nfeature[ni]);
	      RecYtally theRec = (RecYtally) hmtallysum.get(objd);

	      xvaluescount_nfeature[ni] = theRec.ntally;
	      yvaluessum_nfeature[ni] = theRec.dysum;
	      //hmdatatoindex_nfeature.put(new Float(xvals_nfeature[ni]), Integer.valueOf(ni));
	      hmdatatoindex_nfeature.put(Float.valueOf(xvals_nfeature[ni]), Integer.valueOf(ni));
	   }	   
	}
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Goes through each feature and each x-value for the feature and finds the best one split on returns the split to Bestsplit
     */
    public BestSplit evaluateSplits(TreeNode theTreeNode, int ndatalocationssize, int[][] xvaluescount, double[][] yvaluessum)
    {
	//variables associated with the best split
	int nbestdataindex = -1;
	int nbestfeatureindex = 0;
	double dbestdataval = -1;
	double dbestleftmean = 0;
	double dbestrightmean = 0;
	double dbestrandomval = 2;

	double dbestSS = Float.MAX_VALUE;

	int numfeatures = data.length;

	for (int nfeature = 0; nfeature < numfeatures; nfeature++)
	{
	    //going through each feature 

	    //HashMap hmtallysum = new HashMap();

	   //initally assume all values are to the right of the split value
	   double dsumlefty = 0;
	   double dsumrighty = 0;
	   int ncountleft = 0;
	   int ncountright = ndatalocationssize;

	   int[] xvaluescount_nfeature = xvaluescount[nfeature];
	   double[] yvaluessum_nfeature = yvaluessum[nfeature];
	   float[] xvals_nfeature = xvals[nfeature];

	   for (int ni = 0; ni < yvaluessum_nfeature.length; ni++)
	   {
	       //computing the sum of all y-values ot the right
	      dsumrighty += yvaluessum_nfeature[ni];
	   }

	   //computes the average value to the right
	   double davgrighty = dsumrighty/(double) ncountright;

	       	     
           double dSS = -davgrighty*dsumrighty;

	   //gets the last index to xvals_nfeature
	   int nxvalslengthm1 = xvals[nfeature].length -1;

	   if (dSS <= dbestSS)
	   {
	      double drandomval = theRandom.nextDouble();
              if ((dSS < dbestSS) || ((dSS==dbestSS)&&(drandomval < dbestrandomval)))
              {
		  //best split is assigning everything to the left node, the right node will not be used
	         dbestSS = dSS;
	         dbestrandomval = drandomval;
	         nbestfeatureindex = nfeature;
	         nbestdataindex = ncountright-1;
	         dbestleftmean = dsumrighty/(double) ncountright;
	         dbestrightmean =Double.MAX_VALUE;
                 dbestdataval = xvals_nfeature[nxvalslengthm1];
	      }
	   }

	   for (int ni = 0; ni < nxvalslengthm1; ni++)
	   {	   
	       //gets the current count of the number of instances associated with this feature
	      int ntally = xvaluescount_nfeature[ni];
	      if (ntally > 0)
	      {
		  //only going to pursue this x-value as a split if some observed value has been associated with it
	         double dysum = yvaluessum_nfeature[ni];

		 //increases the counts associated with splitting to the left
  	         ncountleft += ntally;
		 dsumlefty += dysum;

		 //decreases the counts 
		 ncountright -= ntally;
		 dsumrighty -= dysum;

		 //computes the updated averages for the left and right
	         double davglefty = dsumlefty/(double) ncountleft;
	         davgrighty = dsumrighty/(double) ncountright;
	       	     
                 dSS = -davglefty*dsumlefty-davgrighty*dsumrighty;
		 //derivation
		 //sum i = 1 to n (x_i - avg(x_L))^2 - sum i = n+1 to R (x_i - avg(X_R))^2
		 //sum i = 1 to n (-2*x_i*avg(x_L) + avg(x_L)^2) + sum i = n+1 to R (-2*x_i*avg(x_R) + avg(x_R)^2)
		 //avg(x_L)*(sum i = 1 to n (-2*x_i + avg(x_L))) + avg(x_R)*(sum i = n+1 to R (-2*x_i + avg(x_R)))
		 //avg(x_L)*(- sum i = 1 to n (x_i)) + avg(x_R)*(- sum i = n+1 to R (x_i))

		 if (dSS <= dbestSS)
		 {
		     //computes a random value  to use in breaking ties
		    double drandomval = theRandom.nextDouble();

                    if (((dSS < dbestSS)||((dSS==dbestSS)&&(drandomval<dbestrandomval)))&&(ncountleft>= nminnumlocations)&& (ncountright>= nminnumlocations))
	            {
			//only accepts splits which improves the bestSS value or is tie and has a lower random value
			//also require the counts to the left and right are greater than the nminnumlocations parameter

			//updates the best sum of square value
	               dbestSS = dSS;
		       dbestrandomval = drandomval;

		       //updates best tie value
	               nbestfeatureindex = nfeature;

		       //index of last position in sorted order that would be less than split val
	               nbestdataindex = ncountleft-1;

	               dbestleftmean = davglefty;//dsumlefty/(double) ncountleft;
                       dbestrightmean = davgrighty;// dsumrighty/(double) ncountright;

		       //value corresponding to the best split feature
                       dbestdataval = xvals_nfeature[ni];
		    }
		 }
	      }
	   }
	}

	//returns the best split found so far
        return new BestSplit(nbestfeatureindex, nbestdataindex, dbestdataval, dbestleftmean, dbestrightmean);//,dbestSS);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * Procedure for constructing the regression tree
     */
    public void buildTree(TreeNode currNode,int[] datalocations, int ndatalocationssize, int[][] xvaluescount,double[][] yvaluessum,int ndepth)
    {
       if (ndatalocationssize >= 2*nminnumlocations)
       { 
	   //it is possible to find a split so each leaf node has at least minnum locations
	  BestSplit theBestSplit = evaluateSplits(currNode,ndatalocationssize, xvaluescount, yvaluessum);	   

	  //sets the split feature index and value 
  	  currNode.nsplitfeatureindex = theBestSplit.nbestfeatureindex;
    	  currNode.dsplitval = theBestSplit.dbestdataval;

	  //stores the number of features which is number of rows of data
	  int numfeatures = data.length;

	  //stores the number of data elements going to the left
	  int nsizeleft = theBestSplit.nbestdataindex+1;

	  //stores the number of data elements going to the right
	  int nsizeright = ndatalocationssize -theBestSplit.nbestdataindex-1;

	  int[] leftlocations = datalocations;
	  int[] rightlocations;
	  
	  if (ndepth < rightlocationsAL.size())
	  {
	      //we've been to this depth before
	      //get current array for depth
	     rightlocations = (int[]) rightlocationsAL.get(ndepth);
	     if (rightlocations.length < nsizeright)
	     {
		 //check if need to reallocated to make larger
	        rightlocations = new int[nsizeright];

		//stores the new array back into the array list
	        rightlocationsAL.set(ndepth,rightlocations);
	     }
	     //updates the number of elements at the right side of this depth
	     rightlocationsSizeAL.set(ndepth, Integer.valueOf(nsizeright));
	  }
	  else
	  {
	      //first time at this depth allocating an array for it to store indicies
	     rightlocations = new int[nsizeright];
	     rightlocationsAL.add(rightlocations);

	     //stores the size of the array
	     rightlocationsSizeAL.add(Integer.valueOf(nsizeright));
	  }

	  int[][] xvaluescountright;
	  double[][] yvaluessumright;

	  int[][] xvaluescountleft;
	  double[][] yvaluessumleft;
	  
          if ((nsizeleft >0)&&(nsizeright >0))
	  {	      
	      //we are going to need to perform a split since both sides are >0

	     int nleftindex = 0;
	     int nrightindex = 0;
	     float[] data_nsplitfeatureindex = data[currNode.nsplitfeatureindex];

	     if (nsizeleft <= nsizeright)
	     {
		 //left size is smaller than right side
	        if (ndepth < xvaluescountAL.size())
	        {
		    //we've been at this depth before
		    //getting the stored contents
		   xvaluescountleft = (int[][]) xvaluescountAL.get(ndepth);
		   yvaluessumleft = (double[][]) yvaluessumAL.get(ndepth);

	           for (int nfeature = 0; nfeature < xvaluescountleft.length; nfeature++)
		   {
		      int[] xvaluescountleft_nfeature = xvaluescountleft[nfeature];
		      double[] yvaluessumleft_nfeature = yvaluessumleft[nfeature];

		      for (int nindex = 0; nindex < xvaluescountleft_nfeature.length; nindex++)
		      {
			  //resets the count and y-sum values to 0 for all x-features values
		         xvaluescountleft_nfeature[nindex] = 0;
		         yvaluessumleft_nfeature[nindex]= 0;
		      }
		   }
		}
		else
	        {
		    //first time at this depth allocating feature array
	           xvaluescountleft = new int[xvaluescount.length][];
	           yvaluessumleft = new double[yvaluessum.length][];

	           for (int ni = 0; ni < xvaluescountleft.length; ni++)
	           {
		       //creating an x-value count for each x-split value
	              xvaluescountleft[ni] = new int[xvaluescount[ni].length];

		      //creating a y-sum value for each x-split value
	              yvaluessumleft[ni] = new double[yvaluessum[ni].length];
	           }

		   //stores the allocated arrays
		    xvaluescountAL.add(xvaluescountleft);
		    yvaluessumAL.add(yvaluessumleft);
		 }

		 for (int ni = 0; ni < ndatalocationssize; ni++)
	         {
	            int nmappedindex = datalocations[ni];
		    float fcurrval = data_nsplitfeatureindex[nmappedindex];
	            if (fcurrval <= currNode.dsplitval)
	            {
			//value is less than split value which occurs less than half the time

		       float output_nmappedindex = output[nmappedindex];
		       for (int nfeature = 0; nfeature < numfeatures; nfeature++)
		       {
			   //updating the left side xvalues and y-sum for including this value on the left side
			  float fval = data[nfeature][nmappedindex]; //gets the actual value for the feature
		          //int ndataxvalindex = ((Integer) hmdatatoindex[nfeature].get(new Float(fval))).intValue();
		          int ndataxvalindex = ((Integer) hmdatatoindex[nfeature].get(Float.valueOf(fval))).intValue();
			  //updating count of how often to the left side we have this feature value
		          xvaluescountleft[nfeature][ndataxvalindex]++;
			  //updating corresponding sum of y for this feature value
		          yvaluessumleft[nfeature][ndataxvalindex] += output_nmappedindex;			  
		       }

	   	       leftlocations[nleftindex] = nmappedindex;
		       nleftindex++;
		    }
	            else
	            {
		       rightlocations[nrightindex] = nmappedindex;
		       nrightindex++;
		    }
		 }

		 //the right side inherits the parents and will then subtract out the left side
		 xvaluescountright = xvaluescount;
		 yvaluessumright = yvaluessum;
		 for (int nfeature = 0; nfeature < numfeatures; nfeature++)
		 {
		    int[] xvaluescountright_nfeature = xvaluescountright[nfeature];
		    double[] yvaluessumright_nfeature = yvaluessumright[nfeature];
		    double[] yvaluessumleft_nfeature = yvaluessumleft[nfeature];
		    int[] xvaluescountleft_nfeature = xvaluescountleft[nfeature];

		    for (int nxvalindex = 0; nxvalindex < xvaluescountleft_nfeature.length; nxvalindex++)
		    { 
			//subtract out the left size for each feature
		       xvaluescountright_nfeature[nxvalindex] -= xvaluescountleft_nfeature[nxvalindex]; 
		       yvaluessumright_nfeature[nxvalindex] -= yvaluessumleft_nfeature[nxvalindex];
		    }
		 }
	     }
	     else
	     {
		 //left side is larger going to iterate through right
	        if (ndepth < xvaluescountAL.size())
	        {
		    //we've been to this depth before
		   xvaluescountright = (int[][]) xvaluescountAL.get(ndepth);
		   yvaluessumright = (double[][]) yvaluessumAL.get(ndepth);

		   for (int nfeature = 0; nfeature < xvaluescountright.length; nfeature++)
		   {
		      int[] xvaluescountright_nfeature = xvaluescountright[nfeature];
		      double[] yvaluessumright_nfeature = yvaluessumright[nfeature];

		      for (int nindex = 0; nindex < xvaluescountright_nfeature.length; nindex++)
		      {
			  //clearing out the counts for the feature values
		         xvaluescountright_nfeature[nindex] = 0;
		         yvaluessumright_nfeature[nindex]= 0;
		      }
		   }
		}
		else
	        {
		    //first time at this depth allocating space for the x-count values and y-sum values
		   xvaluescountright = new int[xvaluescount.length][];
		   yvaluessumright = new double[yvaluessum.length][];

	           for (int ni = 0; ni < xvaluescountright.length; ni++)
		   {
	              xvaluescountright[ni] = new int[xvaluescount[ni].length];
		      yvaluessumright[ni] = new double[yvaluessum[ni].length];
	           }

	           xvaluescountAL.add(xvaluescountright);
	           yvaluessumAL.add(yvaluessumright);
		}

		for (int ni = 0; ni < ndatalocationssize; ni++)
	        {
		    //going through all the data location indicies
	           int nmappedindex = datalocations[ni];
		   float fcurrval = data_nsplitfeatureindex[nmappedindex];

	           if (fcurrval <= currNode.dsplitval)
	           {
		       //the value is less than or equal to the split value
		       //storing this is a left index value
	              leftlocations[nleftindex] = nmappedindex;		    
	              nleftindex++;
	           }
	           else
	           {
		       
		      float output_nmappedindex = output[nmappedindex];
		   
		      //computing the right tally and output sum for each feature value
		      for (int nfeature = 0; nfeature < numfeatures; nfeature++)
		      {
		         float fval = data[nfeature][nmappedindex];
		         //int ndataxvalindex = ((Integer) hmdatatoindex[nfeature].get(new Float(fval))).intValue();
		         int ndataxvalindex = ((Integer) hmdatatoindex[nfeature].get(Float.valueOf(fval))).intValue();
		         xvaluescountright[nfeature][ndataxvalindex]++;
		         yvaluessumright[nfeature][ndataxvalindex] += output_nmappedindex;			 
		      }

		      //storing this in the right index
		      rightlocations[nrightindex] = nmappedindex;
		      nrightindex++;
		   }
		}

		//left side inherits parent and now decrementing right side values
		xvaluescountleft = xvaluescount;
		yvaluessumleft = yvaluessum;

		for (int nfeature = 0; nfeature < numfeatures; nfeature++)
		{
		    int[] xvaluescountleft_nfeature = xvaluescountleft[nfeature];
		    double[] yvaluessumleft_nfeature = yvaluessumleft[nfeature];
		    int[] xvaluescountright_nfeature = xvaluescountright[nfeature];
		    double[] yvaluessumright_nfeature = yvaluessumright[nfeature];

		    for (int nxvalindex = 0; nxvalindex < xvaluescountleft_nfeature.length; nxvalindex++)
		    { 
			xvaluescountleft_nfeature[nxvalindex] -= xvaluescountright_nfeature[nxvalindex]; 
			yvaluessumleft_nfeature[nxvalindex] -= yvaluessumright_nfeature[nxvalindex];
		    }
		}
	     }

	     
	     //creates a new left node
	     TreeNode leftnode = new TreeNode();
	     leftnode.left = null;
	     leftnode.right = null;
	     leftnode.dmean = theBestSplit.dbestleftmean;
             currNode.left = leftnode;

	     //creates a new right node
             TreeNode rightnode = new TreeNode();
	     rightnode.left = null;
             rightnode.right = null;
	     rightnode.dmean = theBestSplit.dbestrightmean;
	     currNode.right = rightnode;

	     //recurisvely build tree for left and right sides	     
	     buildTree(leftnode, leftlocations,nsizeleft,  xvaluescountleft, yvaluessumleft, ndepth + 1);
	     buildTree(rightnode, rightlocations,nsizeright, xvaluescountright, yvaluessumright, ndepth +1);
	  }
       }
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Converts the contents of theTree to a string
     */
    public String toString()    
    {
	StringBuffer sbtreestring = new StringBuffer();
	traverse(theTree,sbtreestring);
	return sbtreestring.toString();

    }


    /////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Converts dval to a string with NumberFormat nf and then parses the leading 0
     */
    public String numformat(double dval)
    {

        String szformat = nf.format(dval);

        if (szformat.startsWith("0."))
        {
           return szformat.substring(1);
	}
        else
        {
	   return szformat;
        }
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Traverses the tree appending the conents to sbtreestring
     */
    public void traverse(TreeNode theTreeNode, StringBuffer sbtreestring)
    {
	if (theTreeNode != null)
	{
	    if (theTreeNode.left == null)
	    {
		//at a leaf node outputting the leaf value
		sbtreestring.append("n\t"+numformat(theTreeNode.dmean)+"\n");
	    }
	    else
	    {
		//at a non-leaf node outputting the split index and value, then subtraversing left and then right trees
		sbtreestring.append(theTreeNode.nsplitfeatureindex+"\t"+numformat(theTreeNode.dsplitval)+"\n");
	        traverse(theTreeNode.left,sbtreestring);
	        traverse(theTreeNode.right,sbtreestring);
	    }
	}
    }


    ////////////////////////////////////////////////////////////////////////////

    /**
     * Given a vector of values in the Instance determines the corresponding 
     * root value
     * if theTree is null return 0
     */
    public double classifyInstance(float[] theInstance)
    {
	TreeNode ptr = theTree;
	double dmean = 0;

	while (ptr != null)
	{
	    //walks through the tree to find the value to classify it to
	    dmean = ptr.dmean;
	    if (theInstance[ptr.nsplitfeatureindex] <= ptr.dsplitval)
	    {
		//values less than or equal to the split feature go to the left
		ptr = ptr.left;
	    }
	    else
	    {
		ptr = ptr.right;
	    }
	}

	return dmean;
    }

    ////////////////////////////////////////////////////////////////////////////////

    
}
