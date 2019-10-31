package ernst.ChromImpute;

import java.io.*;
import java.text.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

// import weka.core.matrix.*; //UNCOMMENT for linear regression uses external Weka package

public class ChromImpute {

  // static int NUMLINES = 100000; // number of lines currently stored in memory
  // static int NUMLINES = 50000; // number of lines currently stored in memory
  static int NUMLINES = 10000; // number of lines currently stored in memory

  static int DEFAULTCHROMHMMBIN = 200;

  /** Default value of dna-methylation */
  static boolean DEFAULTDNAMETHYL = false;

  static int DEFAULTMINNUMLOCATIONS = 20;

  static double DEFAULT_EVALPERCENT1 = 1;

  static double DEFAULT_EVALPERCENT2 = 5;

  static String DEFAULTEXTENSION = ".wig.gz";

  /** */
  static int DEFAULTMAXKNN = 10;

  /** */
  int nmaxknn = DEFAULTMAXKNN;

  /** */
  static int DEFAULTNUMBAGS = 1;

  /** */
  int numbags;

  /** */
  String szclassifierdir = null;

  /** */
  String szdistancedir = null;

  /** */
  int nroundval = 10;

  /** */
  float froundval = (float) nroundval;

  // boolean bprintfile = false;

  /** */
  boolean busesamecellfeatures = true;

  /** */
  boolean buseorderfeatures = true;

  // boolean bloadtrainfile = false;

  /** */
  int nholdoutcellrequest = -1;

  /** */
  int nbagrequest = -1;

  // boolean bnormalize = false;

  /**
   * Contains the indicies of positions to sample for each classifier, first dimension is
   * classifiers flattened based on bag, second index is the samples
   */
  int[][] samples;

  /** */
  int nincrementnarrow = ChromImpute.DEFAULTINCREMENTNARROW;
  // 1; //ChromImpute.DEFAULTINCREMENT; //increment for the main cell type

  /** */
  // static int DEFAULTINCREMENT = 1;

  /** */
  // static int DEFAULTMAXOFFSET = 20;

  /** */
  int nincrementwide = ChromImpute.DEFAULTINCREMENTWIDE; // 20;

  /** */
  int nmaxoffsetnarrow = ChromImpute.DEFAULTMAXOFFSETNARROW;

  static int DEFAULTKNNWINDOW = 20;
  static int DEFAULTINCREMENTNARROW = 1;
  static int DEFAULTINCREMENTWIDE = 20;
  static int DEFAULTMAXOFFSETWIDE = 400;
  static int DEFAULTMAXOFFSETNARROW = 20;
  static int DEFAULTMINTOTALENSEMBLE = 0;

  /** */
  int nmaxoffsetwide = ChromImpute.DEFAULTMAXOFFSETWIDE; // 400;

  /** */
  int nknnoffset = ChromImpute.DEFAULTKNNWINDOW;;

  /** */
  static int DEFAULTRES = 25;

  /** */
  static int DEFAULTNUMSAMPLES = 100000;

  /** The directory containing the input signal files to impute */
  String szinputdir;

  /** A table of the cell and mark of the data file to impute */
  String szimputeinfoINfile;

  /** The directory to which the signal files should be written */
  String szoutdir;

  /** The output cell type to which we should impute values for */
  String szoutcell;

  /** The output mark to which we should impute values */
  String szoutmark;

  // String szvalidatefile;
  //
  // The file we are predict for currently exists
  //

  //
  // A table of the cell and mark of the data file to impute
  //
  // String szimputeinfoOUTfile;

  /** Maps each chromosome to size */
  HashMap hmchromsize;

  /** Maps each chromosome to index */
  HashMap hmchromindex;

  /**
   * A file containing information on the chromosomes to include in the imputation and there file
   * lengths
   */
  String szchrominfo;

  /** The resolution at which the imputation should be done */
  int nresolution = ChromImpute.DEFAULTRES;

  /** output file for single imputation */
  String szoutfile;

  /** the maximum chromosome size */
  long nmaxchromsize;

  /** List of chromosomes */
  ArrayList alchrom;

  /** Seed for the random number generator */
  static int DEFAULTSEED = 100;

  /** min number of elements need to in a leaf node */
  int nminnumlocations;

  /** */
  int nseed = ChromImpute.DEFAULTSEED;

  /** The set of all marks that will be used in the imputation */
  HashSet hsmarks;

  /** The set of all cell types that will be used in the imputation */
  HashSet hscells;

  /** A mapping from mark to the list cell types for which there is data for */
  HashMap hmmarkcell;

  /** A mapping from mark,cell combo to a file for which there is data for */
  HashMap hmmarkcellfile;

  // ArrayList altargetrecs;

  // total bins of the chromosomes
  /** total number of locations excluding edge positions */
  long ntotalbins;

  /** */
  int[] noverallpos;

  /** */
  int[] nsampleindexA;

  /** */
  // double dsamplingfraction;

  /** */
  RegressionLinear[] theClassifierLinearA;

  /** */
  RegressionTree[] theClassifierA;

  /** */
  GZIPOutputStream pwout;

  /** */
  GZIPOutputStream pwattributes;

  /** */
  int numclassifiers;

  /** */
  String szpioneermark = null;

  /** */
  HashSet hspioneermarks;

  /** */
  String szchromwant;

  /** */
  String szchromwantgenerate;

  /** */
  NumberFormat nf1;

  // add params to for generating features and applying
  /** */
  boolean bdnamethyl = false;

  // PARAMS used for DNA METHYLATION
  /** */
  String szmethylheader = null;

  /** */
  String szmethylinfo = null;

  /** */
  String szmethylDIR = null;

  /** edge value for methylation */
  static double DEDGEVAL = 0.5;

  /** */
  static int NMETHYLSCALE = 100;

  /** */
  int ntotalsizeDNAmethylvalid;

  /** */
  String[] dnamethylheader;

  /** */
  int[] regularcelltodnamethylindex;

  /** */
  HashMap hmchromdnamethylfile;

  /** */
  double[] methylavg;

  /** */
  ArrayList[] attributes;

  /** */
  ArrayList theTrainInstances;

  /** */
  ArrayList theTrainInstancesOutput;

  /** */
  // ArrayList[] theTestInstances;

  /** */
  GZIPOutputStream[] trainPW;

  /** */
  int[][] samplesdnamethylbinindex;

  /** */
  int[][] samplesdnamethylchromindex;

  /** */
  float[][][] samplesdnamethylvals;

  /** */
  int numdnamethylcells;

  /** */
  String[] marksA;

  /** */
  String[] cellsA;

  /** */
  boolean[][] bmarkcell;

  /** */
  boolean[][] bcellmark;

  /** */
  boolean[] bdnamethylcell;

  /** */
  // ArrayList alinstances;

  /** */
  int numcells;

  /** */
  int nummarks;

  /** */
  int ntargetcell = -2;

  /** */
  int ntargetmark = -2;

  /** */
  ArrayList[][] distmarkcellAL;

  int numsamples;

  String szextension;

  boolean bprintonefile;

  boolean bprintbrowserheader;

  /** if non-null specifies the mark to convert data for */
  String szconvertmark;

  /** if non-null specifies the cell types to convert data for */
  String szconvertcell;

  /** file suffix for imputed data to eval */
  String szevalobserveddir;

  /** file suffix for observed data to eval */
  String szevalobservedfile;

  /** file suffix for imputed data to eval */
  String szevalimputedir;

  /** file suffix for imputed data to eval */
  String szevalimputefile;

  double devalpercent1;

  double devalpercent2;

  private boolean BLINEAR = false;
  // private boolean BLINEAR = true; // linear regression not officially supported, weka calls to
  // support it is commented out
  // flag should stay false unless changing the commented out code

  boolean bmethylavggenome = false;

  // boolean bmethylavgchrom = false;

  boolean btieglobal = false;

  String szevaloutfile = null;

  String szpeakevalfile = null;

  /** signal treshold to exceed to make 1 call when outputing binarization */
  double dsignalthresh;

  /** the binsize of ChromHMM */
  int nbinsize;

  /** whether to include partial lines */
  boolean bpartial;

  /** whether to use names in file when making chromhmm output */
  boolean busenames;

  ////////////////////////////////////////////////////////////////////////////////////////////////

  static class RegressionLinear {
    double[] coeffs;

    RegressionLinear(BufferedReader br) throws IOException {
      String szLine;

      ArrayList al = new ArrayList();
      while ((szLine = br.readLine()) != null) {
        al.add(szLine);
      }
      coeffs = new double[al.size()];

      for (int nindex = 0; nindex < al.size(); nindex++) {
        coeffs[nindex] = Double.parseDouble((String) al.get(nindex));
      }

      br.close();
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  static class BaseDistRec {
    double ddist;
    int ncell;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////

  public static class BaseDistRecCompare implements Comparator, Serializable {
    public int compare(Object o1, Object o2) {
      BaseDistRec r1 = (BaseDistRec) o1;
      BaseDistRec r2 = (BaseDistRec) o2;

      if (r1.ddist < r2.ddist) {
        return -1;
      } else if (r1.ddist > r2.ddist) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /** Constructor for exportToChromHMM */
  public ChromImpute(
      String szchrominfo,
      String szinputdir,
      String szimputeinfoINfile,
      String szoutdir,
      int nresolution,
      double dsignalthresh,
      int nbinsize,
      boolean bpartial,
      boolean busenames)
      throws IOException, IllegalArgumentException {

    this.szchrominfo = szchrominfo;
    this.szinputdir = szinputdir;
    this.szimputeinfoINfile = szimputeinfoINfile;
    this.szoutdir = szoutdir;
    this.nresolution = nresolution;

    this.dsignalthresh = dsignalthresh;
    this.nbinsize = nbinsize;
    this.bpartial = bpartial;
    this.busenames = busenames;

    exportToChromHMM();
  }

  /**
   * Loads information on the data imputation File is of the form input: load cell type, load mark,
   * data file
   */
  public void loadExportInfo() throws IOException {
    // stores the list of all marks
    hsmarks = new HashSet();

    // stores the list of all cell types involved in the imputation
    hscells = new HashSet();

    // stores for each mark an array list with all cell types associated with the mark
    hmmarkcell = new HashMap();

    // stores for a mark-cell combination the associated input data file
    hmmarkcellfile = new HashMap();

    //
    BufferedReader brimputeinfo = Util.getBufferedReader(szimputeinfoINfile);
    String szLine;

    while ((szLine = brimputeinfo.readLine()) != null) {
      StringTokenizer st = new StringTokenizer(szLine, "\t");
      if (st.countTokens() == 0) continue;
      // reading cell and mark in the file
      String szcell = st.nextToken();
      String szmark = st.nextToken();

      // adding mark to the set of marks available
      hsmarks.add(szmark);

      // adding mark to the set of cell types
      hscells.add(szcell);

      ArrayList alcells = (ArrayList) hmmarkcell.get(szmark);
      if (alcells == null) {
        // first time with this mark, creates arraylist associated with it
        alcells = new ArrayList();
        hmmarkcell.put(szmark, alcells);
      }
      // adding cell to the list of cell
      alcells.add(szcell);

      // gets the input file
      if (busenames) {
        String szinfile = st.nextToken();
        // maps a mark cell combination to a file
        hmmarkcellfile.put(szmark + "\t" + szcell, szinfile);
      }

      // reads in all the mark-cell-data file information
    }
    brimputeinfo.close();

    // stores the count of the number of marks and cells
    nummarks = hsmarks.size();
    numcells = hscells.size();

    // stores the marks and cells into the marksA and cellsA array
    // and then sorts them
    marksA = new String[nummarks];
    cellsA = new String[numcells];
    Iterator itrcells = (Iterator) hscells.iterator();
    Iterator itrmark = (Iterator) hsmarks.iterator();

    for (int nmark = 0; nmark < marksA.length; nmark++) {
      marksA[nmark] = (String) itrmark.next();
    }

    for (int ncell = 0; ncell < cellsA.length; ncell++) {
      cellsA[ncell] = (String) itrcells.next();
    }

    Arrays.sort(marksA);
    Arrays.sort(cellsA);
  }

  public void exportToChromHMM() throws IllegalArgumentException, IOException {

    int numgroup = 0;
    if (nbinsize % nresolution != 0) {
      throw new IllegalArgumentException(
          "ChromHMM binsize should be evenly divisible by ChromImpute resolution, but current values of "
              + nbinsize + " " + nresolution + " do not satisfy that");
    } else {
      numgroup = nbinsize / nresolution;
    }

    loadChromInfo();
    loadExportInfo();

    NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
    nf.setMaximumFractionDigits(2);
    nf.setGroupingUsed(false);

    boolean bbinary = (dsignalthresh < Double.MAX_VALUE);

    for (int ncell = 0; ncell < cellsA.length; ncell++) {
      String szcell = cellsA[ncell];

      Iterator itrchrom = (Iterator) hmchromsize.keySet().iterator();

      while (itrchrom.hasNext()) {
        String szchrom = (String) itrchrom.next();

        int nchromsize = ((Integer) hmchromsize.get(szchrom)).intValue();

        int numfulllines = ((int) (nchromsize / (numgroup * nresolution))) * numgroup;

        GZIPOutputStream pw;

        if (bbinary) {
          pw =
              new GZIPOutputStream(
                  new FileOutputStream(szoutdir + "/" + szcell + "_" + szchrom + "_binary.txt.gz"));
        } else {
          pw =
              new GZIPOutputStream(
                  new FileOutputStream(szoutdir + "/" + szcell + "_" + szchrom + "_signal.txt.gz"));
        }

        String szout = szcell + "\t" + szchrom + "\n";

        byte[] btformat = szout.getBytes();
        pw.write(btformat, 0, btformat.length);

        szout = marksA[0];
        for (int nmark = 1; nmark < marksA.length; nmark++) {
          szout += "\t" + marksA[nmark];
        }
        szout += "\n";

        btformat = szout.getBytes();
        pw.write(btformat, 0, btformat.length);

        BufferedReader[] brA = new BufferedReader[marksA.length];
        for (int nmark = 0; nmark < marksA.length; nmark++) {
          String szinfile;
          if (busenames) {
            szinfile =
                szinputdir
                    + "/"
                    + szchrom
                    + "_"
                    + ((String) hmmarkcellfile.get(marksA[nmark] + "\t" + szcell));
          } else {
            szinfile =
                szinputdir + "/" + szchrom + "_impute_" + szcell + "_" + marksA[nmark] + ".wig.gz";
          }

          File f = new File(szinfile);

          if (f.exists()) {
            brA[nmark] = Util.getBufferedReader(szinfile);
          } else {
            throw new IllegalArgumentException("Could not find file " + szinfile);
          }

          brA[nmark].readLine();
          brA[nmark].readLine();
        }

        String szLine;
        boolean bcontinue = true;

        int[] ntotallines = new int[marksA.length];
        while (bcontinue) {
          boolean bfound = false;
          szout = "";

          for (int nmark = 0; nmark < brA.length; nmark++) {
            double dsum = 0;
            int ntally = 0;

            for (int nsmall = 0; nsmall < numgroup; nsmall++) {
              szLine = brA[nmark].readLine();
              if (szLine == null) {
                bcontinue = false;
                break;
              } else {
                ntotallines[nmark]++;
                bfound = true;
                dsum += Double.parseDouble(szLine);
                ntally++;
              }
            }

            if (bfound) {
              double davg = dsum / (double) ntally;
              String szval;
              if (bbinary) {
                if (davg >= dsignalthresh) {
                  szval = "1";
                } else {
                  szval = "0";
                }
              } else {
                szval = nf.format(davg);
              }

              if (nmark >= 1) {
                szout += "\t" + szval;
              } else {
                szout += szval;
              }
            }
          }

          if ((bfound) && (bpartial || (ntotallines[0] <= numfulllines))) {
            szout += "\n";
            btformat = szout.getBytes();
            pw.write(btformat, 0, btformat.length);
          }
        }
        pw.finish();
        pw.close();

        for (int ni = 0; ni < brA.length; ni++) {
          brA[ni].close();
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  /** Constructor for Eval */
  public ChromImpute(
      String szevalobserveddir,
      String szevalobservedfile,
      String szevalimputedir,
      String szevalimputefile,
      String szchrominfo,
      double devalpercent1,
      double devalpercent2,
      boolean bprintbrowserheader,
      boolean bprintonefile,
      String szevaloutfile,
      String szpeakevalfile)
      throws IOException {
    this.szevalobserveddir = szevalobserveddir;
    this.szevalobservedfile = szevalobservedfile;
    this.szevalimputedir = szevalimputedir;
    this.szevalimputefile = szevalimputefile;
    this.szchrominfo = szchrominfo;
    this.devalpercent1 = devalpercent1;
    this.devalpercent2 = devalpercent2;
    this.bprintbrowserheader = bprintbrowserheader;
    this.bprintonefile = bprintonefile;
    this.szevaloutfile = szevaloutfile;
    this.szpeakevalfile = szpeakevalfile;

    loadChromInfo();

    if (szpeakevalfile != null) {
      evaluatePeaks();
    } else {
      evaluate();
    }
    // evaluateROC();
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  /** Constructor for global correlation */
  public ChromImpute(
      String szchrominfo,
      String szinputdir,
      String szimputeinfoINfile,
      String szoutmark,
      String szoutcell,
      String szoutdir,
      int nresolution,
      String szextension)
      throws IOException {
    this.szchrominfo = szchrominfo;
    this.szinputdir = szinputdir;
    this.szimputeinfoINfile = szimputeinfoINfile;
    this.szoutmark = szoutmark;
    this.szoutcell = szoutcell;
    this.szoutdir = szoutdir;
    this.nresolution = nresolution;
    this.szextension = szextension;

    loadChromInfo(); // reads file with the length info for each chromosome
    loadImputeInfo(); // reads information on data to impute file from

    if ((szoutcell != null) && (!hscells.contains(szoutcell))) {
      throw new IllegalArgumentException(
          "Sample " + szoutcell + " was not defined in the sample mark table.");
    }

    if ((szoutmark != null) && (!hsmarks.contains(szoutmark))) {
      throw new IllegalArgumentException(
          "Mark " + szoutmark + " was not defined in the sample mark table.");
    }

    computeGlobalCorrelations();
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  /** Constructor for Apply */
  // TODO verbose imputation apply
  public ChromImpute(
      String szchrominfo,
      String szinputdir,
      String szdistancedir,
      String szimputeinfoINfile,
      String szoutdir,
      String szoutcell,
      String szoutmark, // String szimputeinfoOUTfile,
      int nresolution,
      String szoutfile,
      // String szvalidatefile,
      String szclassifierdir,
      int nmaxknn,
      String szpioneermark,
      boolean busesamecellfeatures,
      boolean buseorderfeatures,
      String szchromwant,
      boolean bdnamethyl,
      int nmintotalensemble,
      int numbags,
      boolean bprintbrowserheader,
      boolean bprintonefile,
      String szmethylheader,
      String szmethylinfo,
      String szmethylDIR,
      int nmaxoffsetnarrow,
      int nmaxoffsetwide,
      int nincrementnarrow,
      int nincrementwide,
      int nknnoffset,
      boolean bmethylavggenome,
      boolean bmethylavgchrom,
      boolean btieglobal
      // String szmethylheader, String szmethylinfo, String szmethylDIR,
      // double dedgeval, double dmethylscale
      ) throws Exception {
    this.nmaxknn = nmaxknn;
    this.bdnamethyl = bdnamethyl;
    this.bmethylavggenome = bmethylavggenome;
    // this.bmethylavgchrom = bmethylavgchrom;
    this.btieglobal = btieglobal;

    if ((!bmethylavggenome) && (!bmethylavgchrom)) {
      if (szchromwant == null) {
        this.bmethylavggenome = true;
      }

      /*
      if (szchromwant != null)
      {
      this.bmethylavgchrom = true;
      }
      else
      {
      this.bmethylavggenome = true;
      }
      */
    }

    // this.numbags = numbags;

    this.nmaxoffsetnarrow = nmaxoffsetnarrow;
    this.nmaxoffsetwide = nmaxoffsetwide;
    this.nincrementnarrow = nincrementnarrow;
    this.nincrementwide = nincrementwide;
    this.nknnoffset = nknnoffset;

    this.bprintbrowserheader = bprintbrowserheader;
    this.bprintonefile = bprintonefile;

    this.busesamecellfeatures = busesamecellfeatures;
    this.buseorderfeatures = buseorderfeatures;
    this.szchrominfo = szchrominfo;
    this.szinputdir = szinputdir;
    this.szdistancedir = szdistancedir;
    this.szpioneermark = szpioneermark;
    hspioneermarks = new HashSet();
    if (szpioneermark != null) {
      StringTokenizer stpioneer = new StringTokenizer(szpioneermark, ",");
      while (stpioneer.hasMoreTokens()) {
        hspioneermarks.add(stpioneer.nextToken());
      }
    }
    this.szimputeinfoINfile = szimputeinfoINfile;
    this.szoutdir = szoutdir;

    this.szoutcell = szoutcell;
    this.szoutmark = szoutmark;
    // this.szimputeinfoOUTfile = szimputeinfoOUTfile;
    this.nresolution = nresolution;
    // this.szoutfile = szoutfile;
    // this.szvalidatefile = szvalidatefile;

    this.szmethylheader = szmethylheader;
    this.szmethylinfo = szmethylinfo;
    this.szmethylDIR = szmethylDIR;

    this.szclassifierdir = szclassifierdir;
    this.szchromwant = szchromwant;

    nf1 = NumberFormat.getInstance(Locale.ENGLISH);
    nf1.setMaximumFractionDigits(1);
    nf1.setGroupingUsed(false);
    if (szoutfile != null) {
      this.szoutfile = szoutfile; // +"_"+numbags+"_"+numsamples+"_"+nminnumlocations+".wig";
    } else {
      this.szoutfile =
          "impute_" + szoutcell + "_" + szoutmark
              + ".wig"; // +"_"+numbags+"_"+numsamples+"_"+nminnumlocations+".wig";
    }

    if (bdnamethyl) {
      loadDNAMethylHeader();
    }
    loadChromInfo();
    loadImputeInfo();

    if (!hscells.contains(szoutcell)) {
      throw new IllegalArgumentException(
          "Sample " + szoutcell + " was not defined in the sample mark table.");
    }

    if ((!hsmarks.contains(szoutmark)) && (!bdnamethyl)) {
      throw new IllegalArgumentException(
          "Mark " + szoutmark + " was not defined in the sample mark table.");
    }

    if (bdnamethyl) {
      this.numbags =
          (int) Math.max(numbags, Math.ceil(nmintotalensemble / (double) numdnamethylcells));
    } else {
      ArrayList alcells = (ArrayList) hmmarkcell.get(szoutmark);
      this.numbags =
          (int) Math.max(numbags, Math.ceil(nmintotalensemble / (double) alcells.size()));
      // fixed to update numbags count to subtract
      // System.out.println("numbags is "+this.numbags+" "+nmintotalensemble+" "+alcells.size());
    }

    // loadTargets();
    loadDistInfo();
    executeApply();
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  /** Constructor for train features */
  // TODO FIX HERE
  public ChromImpute(
      String szchrominfo,
      String szinputdir, // String szdistancedir,
      String szimputeinfoINfile,
      String szoutdir,
      String szoutcell,
      String szoutmark, // String szimputeinfoOUTfile,
      // int nresolution,String szoutfile,//String szvalidatefile,
      int nholdoutcellrequest,
      int nbagrequest,
      // boolean bprintfile, boolean bloadtrainfile,
      // int nseed,
      int nmaxknn,
      String szpioneermark,
      boolean busesamecellfeatures,
      boolean buseorderfeatures,
      boolean bdnamethyl,
      int nminnumlocations,
      int nmintotalensemble,
      int numbags,
      String szmethylheader)
      throws Exception {
    this.szmethylheader = szmethylheader;
    // this.nmaxknn = nmaxknn;
    this.bdnamethyl = bdnamethyl;
    this.busesamecellfeatures = busesamecellfeatures;
    this.buseorderfeatures = buseorderfeatures;
    this.nholdoutcellrequest = nholdoutcellrequest;
    this.nbagrequest = nbagrequest;
    // this.szchrominfo = szchrominfo;
    this.szinputdir = szinputdir;
    // this.szdistancedir = szdistancedir;
    this.szimputeinfoINfile = szimputeinfoINfile;
    this.szpioneermark = szpioneermark;
    this.nminnumlocations = nminnumlocations;
    hspioneermarks = new HashSet();
    if (szpioneermark != null) {
      StringTokenizer stpioneer = new StringTokenizer(szpioneermark, ",");
      while (stpioneer.hasMoreTokens()) {
        hspioneermarks.add(stpioneer.nextToken());
      }
    }
    this.szoutdir = szoutdir;
    this.szoutcell = szoutcell;
    this.szoutmark = szoutmark;
    // this.szimputeinfoOUTfile = szimputeinfoOUTfile;
    // this.nresolution = nresolution;
    // this.szoutfile = szoutfile;
    // this.szvalidatefile = szvalidatefile;
    // this.bprintfile = bprintfile;
    // this.bloadtrainfile = bloadtrainfile;
    // this.nseed = nseed;

    nf1 = NumberFormat.getInstance(Locale.ENGLISH);
    nf1.setMaximumFractionDigits(1);
    nf1.setGroupingUsed(false);

    // Specific function for DNase-seq and Methylation
    if (bdnamethyl) {
      loadDNAMethylHeader();
    }

    // loadChromInfo();   // Reads in chromsizes
    loadImputeInfo(); // reads information on data to impute file from
    // Looks like this:
    // hscells cells list
    // hsmarks marks list 
    // hmmarkcell for each mark all cells assoc with it
    // hmmarkcellfile file for a mark-cell combo

    // Check that the call is to a real cell and a real mark (separately):
    if (!hscells.contains(szoutcell)) {
      throw new IllegalArgumentException(
          "Sample " + szoutcell + " was not defined in the sample mark table.");
    }

    if ((!hsmarks.contains(szoutmark)) && (!bdnamethyl)) {
      throw new IllegalArgumentException(
          "Mark " + szoutmark + " was not defined in the sample mark table.");
    }

    if (bdnamethyl) {
      this.numbags =
          (int) Math.max(numbags, Math.ceil(nmintotalensemble / (double) numdnamethylcells));
    } else {
      ArrayList alcells = (ArrayList) hmmarkcell.get(szoutmark);
      //removed in v0.9.6 //
      // if ((alcells == null) || (alcells.size() == 1)) {
      //   this.numbags = numbags;
      // } else {
      this.numbags =
          (int) Math.max(numbags, Math.ceil(nmintotalensemble / (double) alcells.size()));
      // }
      // System.out.println("numbags is "+this.numbags+" "+nmintotalensemble+" "+alcells.size());
      // TODO ADDED THIS FOR DEBUGGING
      // System.out.println("numbags is "+this.numbags+" "+nmintotalensemble+" "+alcells.size());
    }
    executeTrain();
  }

  /////////////////////////////////////////////////////////////////////////////////////

  /** Constructor for generate features */
  public ChromImpute(
      String szchromwantgenerate,
      String szchrominfo,
      String szinputdir,
      String szdistancedir,
      String szimputeinfoINfile,
      String szoutdir,
      String szoutcell,
      String szoutmark, // String szimputeinfoOUTfile,
      int nresolution,
      String szoutfile,
      // int nholdoutcellrequest,int nbagrequest,
      // boolean bprintfile, boolean bloadtrainfile,
      int nseed,
      int nmaxknn,
      boolean bdnamethyl,
      int numsamples,
      int nmintotalensemble,
      int numbags,
      String szmethylheader,
      String szmethylinfo,
      String szmethylDIR,
      int nmaxoffsetnarrow,
      int nmaxoffsetwide,
      int nincrementnarrow,
      int nincrementwide,
      int nknnoffset,
      boolean bmethylavggenome,
      boolean bmethylavgchrom,
      boolean btieglobal)
      throws Exception {

    this.szchromwantgenerate = szchromwantgenerate;
    this.bmethylavggenome = bmethylavggenome;
    // this.bmethylavgchrom = bmethylavgchrom;
    this.btieglobal = btieglobal;

    if ((!bmethylavggenome) && (!bmethylavgchrom)) {
      if (szchromwantgenerate == null) {
        this.bmethylavggenome = true;
      }

      /*
      if (szchromwantgenerate != null)
      {
      this.bmethylavgchrom = true;
      }
      else
      {
      this.bmethylavggenome = true;
      }
      */
    }

    this.szchrominfo = szchrominfo;
    this.szinputdir = szinputdir;
    this.szdistancedir = szdistancedir;
    this.szimputeinfoINfile = szimputeinfoINfile;
    this.szoutdir = szoutdir;
    this.szoutcell = szoutcell;
    this.szoutmark = szoutmark;
    // this.szimputeinfoOUTfile = szimputeinfoOUTfile;
    this.nresolution = nresolution;
    this.szoutfile = szoutfile;
    // this.nholdoutcellrequest = nholdoutcellrequest;
    // this.nbagrequest = nbagrequest;
    // this.bprintfile = true;
    // this.bloadtrainfile = false;

    // this.bprintfile = bprintfile;
    // this.bloadtrainfile = bloadtrainfile;
    this.nseed = nseed;
    this.nmaxknn = nmaxknn;
    this.bdnamethyl = bdnamethyl;
    this.numsamples = numsamples;
    // this.numbags = numbags;

    // this.busesamecellfeatures =busesamecellfeatures;
    // this.buseorderfeatures = buseorderfeatures;

    this.szmethylheader = szmethylheader;
    this.szmethylinfo = szmethylinfo;
    this.szmethylDIR = szmethylDIR;
    this.nmaxoffsetnarrow = nmaxoffsetnarrow;
    this.nmaxoffsetwide = nmaxoffsetwide;
    this.nincrementnarrow = nincrementnarrow;
    this.nincrementwide = nincrementwide;
    this.nknnoffset = nknnoffset;

    // this.szpioneermark = szpioneermark;
    // hspioneermarks  =new HashSet();
    // if (szpioneermark != null)
    // {
    // StringTokenizer stpioneer = new StringTokenizer(szpioneermark,",");
    // while (stpioneer.hasMoreTokens())
    //  {
    //    hspioneermarks.add(stpioneer.nextToken());
    // }
    // }

    // szoutfile =
    // "impute_"+szoutcell+"_"+szoutmark+".wig";//+"_"+numbags+"_"+numsamples+"_"+nminnumlocations+".wig";

    nf1 = NumberFormat.getInstance(Locale.ENGLISH);
    nf1.setMaximumFractionDigits(1);
    nf1.setGroupingUsed(false);

    if (bdnamethyl) {
      loadDNAMethylHeader();
    }
    loadChromInfo();

    if ((szchromwantgenerate != null) && (hmchromsize.get(szchromwantgenerate) == null)) {
      throw new IllegalArgumentException(
          szchromwantgenerate + " is not a valid chromosome as listed in the file " + szchrominfo + "!");
    }

    loadImputeInfo();

    if ((!hsmarks.contains(szoutmark)) && (!bdnamethyl)) {
      throw new IllegalArgumentException(
          "Mark " + szoutmark + " was not defined in the sample mark table.");
    }

    if (bdnamethyl) {
      this.numbags =
          (int) Math.max(numbags, Math.ceil(nmintotalensemble / (double) numdnamethylcells));
    } else {
      ArrayList alcells = (ArrayList) hmmarkcell.get(szoutmark);
      this.numbags =
          (int) Math.max(numbags, Math.ceil(nmintotalensemble / (double) alcells.size()));
      // System.out.println("numbags is "+this.numbags+" "+nmintotalensemble+" "+alcells.size());
    }

    // loadTargets();
    loadDistInfo();

    executeGenerateTraining();
  }

  /////////////////////////////////////////////////////////////////////////////////////
  /** Constructor for converting into fixed resolution form */
  public ChromImpute(
      String szchrominfo,
      String szinputdir,
      String szimputeinfoINfile,
      String szoutdir,
      int nresolution,
      String szchromwant,
      String szconvertmark,
      String szconvertcell)
      throws IOException {
    this.szchrominfo = szchrominfo;
    this.szinputdir = szinputdir;
    this.szimputeinfoINfile = szimputeinfoINfile;
    this.szoutdir = szoutdir;
    this.nresolution = nresolution;
    this.szchromwant = szchromwant;
    this.szconvertmark = szconvertmark;
    this.szconvertcell = szconvertcell;

    nf1 = NumberFormat.getInstance(Locale.ENGLISH);
    nf1.setMaximumFractionDigits(1);
    nf1.setGroupingUsed(false);

    if (bdnamethyl) {
      loadDNAMethylHeader();
    }
    loadChromInfo(); // reads file with the length info for each chromosome
    loadImputeInfo(); // reads information on data to impute file from

    if ((szconvertcell != null) && (!hscells.contains(szconvertcell))) {
      throw new IllegalArgumentException(
          "Sample " + szconvertcell + " was not defined in the sample mark table.");
    }

    if ((szconvertmark != null) && (!hsmarks.contains(szconvertmark))) {
      throw new IllegalArgumentException(
          "Mark " + szconvertmark + " was not defined in the sample mark table.");
    }

    convertData(); // converts data into desired format
  }

  //////////////////////////////////////////////////////////////////////////////////
  /** */
  public void loadDNAMethylHeader() throws IOException {
    BufferedReader brheaderfile = Util.getBufferedReader(szmethylheader);
    String szLine = brheaderfile.readLine();
    StringTokenizer st = new StringTokenizer(szLine, "\t");
    numdnamethylcells = st.countTokens() - 1; // first column is position

    dnamethylheader = new String[numdnamethylcells];
    st.nextToken(); // flush position element

    // reads in the methylation header
    for (int nkindex = 0; nkindex < numdnamethylcells; nkindex++) {
      dnamethylheader[nkindex] = st.nextToken();
    }

    brheaderfile.close();
  }

  /////////////////////////////////////////////////////////////////////////////////////////////

  /** Reads in a tab delimited file which gives each chromosome and its length */
  public void loadChromInfo() throws IOException {
    // System.out.println("Reading "+szchrominfo);
    // System.out.println(bmethylavggenome+"\t"+bmethylavgchrom+"\t"+szchromwant+"\t"+szchromwantgenerate);
    BufferedReader brchrom = Util.getBufferedReader(szchrominfo);

    hmchromsize = new HashMap(); // maps each chrom to a size
    hmchromindex = new HashMap(); // maps each chromosome to an index
    alchrom = new ArrayList(); // stores the chromosomes in an array list
    ArrayList alchromALL = new ArrayList();

    String szLine;
    long ntotalsize = 0;
    int nindex = 0;

    while ((szLine = brchrom.readLine()) != null) {
      StringTokenizer st = new StringTokenizer(szLine, "\t");
      if (st.countTokens() == 0) continue;
      String szchrom = st.nextToken();

      alchromALL.add(szchrom);

      if ((szchromwant == null) || (szchromwant.equals(szchrom))) {
        int nsize = Integer.parseInt(st.nextToken());

        alchrom.add(szchrom);
        hmchromsize.put(szchrom, Integer.valueOf(nsize));
        if (nsize > nmaxchromsize) {
          // chromosome size larger than the maximum found
          nmaxchromsize = nsize;
        }

        // increment total size of the chromosome
        ntotalsize += nsize;
        hmchromindex.put(szchrom, Integer.valueOf(nindex));
        nindex++;
      }
    }
    brchrom.close();

    // determines the total number of bins excluding thoses on the ends in which do not have
    // sufficient number of bins
    int nmaxparam = Math.max(nmaxoffsetwide, nknnoffset);
    ntotalbins = ntotalsize / nresolution - 2 * nmaxparam * nindex;

    if (bdnamethyl) {

      // int numchrom = nindex;
      // System.out.println("numchrom is "+numchrom);
      // BufferedReader brheaderfile = Util.getBufferedReader(szmethylheader);
      // szLine = brheaderfile.readLine();
      // StringTokenizer st = new StringTokenizer(szLine,"\t ");

      ntotalsizeDNAmethylvalid = 0;
      //
      // numdnamethylcells = st.countTokens()-1; //first column is position
      methylavg = new double[numdnamethylcells];
      double[] methylcount = new double[numdnamethylcells];

      // has for each chrom the name of the DNA methylation file
      BufferedReader brdnamethylinfo = Util.getBufferedReader(szmethylinfo);

      hmchromdnamethylfile = new HashMap();
      StringTokenizer st;
      ArrayList alchromdnamethylinfo = new ArrayList();
      while ((szLine = brdnamethylinfo.readLine()) != null) {
        st = new StringTokenizer(szLine, "\t");
        if (st.countTokens() == 0) continue;
        String szchrom = st.nextToken();
        alchromdnamethylinfo.add(szchrom);
        String szfile = st.nextToken();
        hmchromdnamethylfile.put(szchrom, szfile);
      }
      brdnamethylinfo.close();

      // dnamethylheader = new String[numdnamethylcells];
      // st.nextToken(); //flush position element

      // reads in the methylation header
      // for (int nkindex = 0; nkindex < numdnamethylcells; nkindex++)
      // {
      //   dnamethylheader[nkindex] = st.nextToken();
      // }

      // has for each chrom the name of the DNA methylation file
      // for (int nchromindex = 0; nchromindex < numchrom; nchromindex++)
      for (int nchromindex = 0; nchromindex < alchromALL.size(); nchromindex++) {
        // String szchrom = (String) alchromALL.get(nchromindex);
        String szchrom = (String) alchromdnamethylinfo.get(nchromindex);
        // String szchrom = (String) alchromdnamethyl.get(nchromindex);

        // don't need to read this file if not averaging genomewide, and only one apply chromosome
        // selected
        if ((!bmethylavggenome) && (szchromwant != null) && (!szchromwant.equals(szchrom)))
          continue;

        String szfile = (String) hmchromdnamethylfile.get(szchrom);

        if (szfile == null) {
          throw new IllegalArgumentException(
              "Chromosome " + szchrom + " not found in DNA methylation info file " + szmethylinfo);
        }
        // System.out.println("opening "+szmethylDIR+"/"+szfile);
        // System.out.println("expecting "+methylavg.length+" columns");
        BufferedReader brdnamethyldata = Util.getBufferedReader(szmethylDIR + "/" + szfile);

        int nchromsize = 0;
        if (szchromwant == null) {
          // always go here when training
          // getting this as long as not applying to specific chromosome
          Integer intobjchromsize = (Integer) hmchromsize.get(szchrom);
          // if (intobjchromsize == null) //updated in v0.9.6 to remove this as wan't doing anything
          // && (szchromwant != null))
          // {
          // don't need chromosome size
          //   throw new IllegalArgumentException("Chromosome "+szchrom+" found in DNA methylation
          // file but not present in "+szchrominfo);
          // }
          nchromsize = ((Integer) intobjchromsize).intValue();
        }

        while ((szLine = brdnamethyldata.readLine()) != null) {
          // first value is the coordinate
          StringTokenizer stdnamethyl = new StringTokenizer(szLine, "\t ");
          // assuming 1-based here
          int ncoord = Integer.parseInt(stdnamethyl.nextToken()) - 1;
          int nbin = ncoord / nresolution;
          for (int nk = 0; nk < methylavg.length; nk++) {
            if (!stdnamethyl.hasMoreTokens()) {
              throw new IllegalArgumentException(
                  "Did not find the expected "
                      + methylavg.length
                      + " number of data columns in line "
                      + szLine);
            }

            double dval = Double.parseDouble(stdnamethyl.nextToken());
            if (dval >= 0) {
              // this is a valid methylation; missing are negative
              // bmethylavggenome, bmethylavgchrom
              if ((bmethylavggenome)
                  || (((szchromwant != null) && szchrom.equals(szchromwant))
                      || ((szchromwantgenerate != null) && szchrom.equals(szchromwantgenerate)))) {
                // contribues to average if genome-wide
                // or specified chromosome and matching it
                methylavg[nk] += dval;
                methylcount[nk]++;
              }
            }
          }

          // chrom size 50
          // array index 0 to 49
          // if 20 is distance can't start from 0 to 19 or 30 to 49

          if ((nbin >= nmaxparam) && (nbin < nchromsize - nmaxparam)) {
            ntotalsizeDNAmethylvalid++;
          }
        }
      }

      // computes the average methylation value
      for (int nk = 0; nk < methylavg.length; nk++) {
        methylavg[nk] /= methylcount[nk];
        // System.out.println("methylavg "+nk+" "+methylavg[nk]);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////

  /** */
  public void generateSamples() {
    // System.out.println("numcells is "+numcells+" "+ntargetcell+" "+nbagrequest+"
    // "+nholdoutcellrequest+" "+numbags);

    samples = new int[numcells * numbags][numsamples]; // index of samples
    nsampleindexA = new int[numcells * numbags]; //
    noverallpos = new int[numcells * numbags];

    int nclassifierindex = 0;
    // for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++)
    {
      for (int nbag = 0; nbag < numbags; nbag++) {
        // System.out.println("Gensamples\t"+nbag);
        // if (((nbagrequest == -1)||(nbagrequest ==nbag))&&
        //   ((nholdoutcellrequest == -1)||(nholdoutcellrequest == nholdoutcell)))
        // {
        // currently using the same set of positions for each classifier of the same bag
        Random theRandom = new Random(nbag + nseed); // (1000*nholdoutcell+nbag+nseed);

        int[] samples_nclassifierindex = samples[nclassifierindex];

        for (int ni = 0; ni < numsamples; ni++) {
          samples_nclassifierindex[ni] = ni;
        }

        // sampling without replacement data locations wanted
        for (int ni = numsamples; ni < ntotalbins; ni++) {
          if (theRandom.nextDouble() < numsamples / ((double) ni + 1)) {
            samples_nclassifierindex[theRandom.nextInt(numsamples)] = ni;
          }
        }

        Arrays.sort(samples_nclassifierindex);
        nsampleindexA[nclassifierindex] = 0;
        noverallpos[nclassifierindex] = 0;
        // }

        // System.out.println("generate samples "+nbag+"
        // "+samples_nclassifierindex[samples_nclassifierindex.length-1]+" "+ntotalbins);
        nclassifierindex++;
      }
    }

    // using samne set of loci for each bag
    for (int nholdoutcell = 1; nholdoutcell < numcells; nholdoutcell++) {
      for (int nbag = 0; nbag < numbags; nbag++) {
        int[] samples_nclassifierindex = samples[nclassifierindex];

        int[] samples_nbag = samples[nbag];
        for (int ni = 0; ni < numsamples; ni++) {
          samples_nclassifierindex[ni] = samples_nbag[ni];
        }

        nsampleindexA[nclassifierindex] = 0;
        noverallpos[nclassifierindex] = 0;
        nclassifierindex++;
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  /** */
  /*
     public void generateSamples()
     {
     System.out.println("numcells is "+numcells+" "+ntargetcell+" "+nbagrequest+" "+nholdoutcellrequest+" "+numbags);

     samples = new int[numcells*numbags][numsamples]; //index of samples
     nsampleindexA = new int[numcells*numbags]; //
     noverallpos = new int[numcells*numbags];

     int nclassifierindex = 0;
     for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++)
     {
     System.out.println("Gensamples\t"+nholdoutcell);
     for (int nbag = 0; nbag < numbags; nbag++)
     {
  //if (((nbagrequest == -1)||(nbagrequest ==nbag))&&
  //   ((nholdoutcellrequest == -1)||(nholdoutcellrequest == nholdoutcell)))
  //{
  //currently using the same set of positions for each classifier of the same bag
  Random theRandom = new Random(nbag+nseed);//(1000*nholdoutcell+nbag+nseed);

  int[] samples_nclassifierindex = samples[nclassifierindex];

  for (int ni = 0; ni < numsamples; ni++)
  {
  samples_nclassifierindex[ni] = ni;
  }

  //sampling without replacement data locations wanted
  for (int ni = numsamples; ni < ntotalbins; ni++)
  {
  if (theRandom.nextDouble() < numsamples/((double) ni +1))
  {
  samples_nclassifierindex[theRandom.nextInt(numsamples)] = ni;
  }
  }

  Arrays.sort(samples_nclassifierindex);
  nsampleindexA[nclassifierindex] = 0;
  noverallpos[nclassifierindex] = 0;
  //}
  nclassifierindex++;
     }
     }
     }
     */

  ////////////////////////////////////////////////////////////////////////////////////////////

  /** */
  public void generateSamplesDNAmethyl(String[] chroms) throws IOException {
    // System.out.println("numcells is "+numcells+" "+ntargetcell+" "+nbagrequest+"
    // "+nholdoutcellrequest+" "+numbags);
    int[] samples = new int[numsamples];
    nsampleindexA = new int[numcells * numbags];
    noverallpos = new int[numcells * numbags];

    samplesdnamethylbinindex = new int[numbags][numsamples];
    samplesdnamethylchromindex = new int[numbags][numsamples];
    samplesdnamethylvals = new float[numbags][numsamples][]; // [numdnamethylcells];

    for (int nbag = 0; nbag < numbags; nbag++) {

      int[] samplesdnamethylbinindex_nbag = samplesdnamethylbinindex[nbag];
      int[] samplesdnamethylchromindex_nbag = samplesdnamethylchromindex[nbag];
      float[][] samplesdnamethylvals_nbag = samplesdnamethylvals[nbag];

      Random theRandom = new Random(nbag + nseed); // (1000*nholdoutcell+nbag+nseed);

      for (int ni = 0; ni < numsamples; ni++) {
        samples[ni] = ni;
      }

      // sampling without replacement data locations wanted
      for (int ni = numsamples; ni < ntotalsizeDNAmethylvalid; ni++) {
        if (theRandom.nextDouble() < numsamples / ((double) ni + 1)) {
          samples[theRandom.nextInt(numsamples)] = ni;
        }
      }

      Arrays.sort(samples);

      int noverallindex = 0;
      int nsampleindex = 0;
      int nmaxparam = Math.max(nmaxoffsetwide, nknnoffset);

      // copies the chromosomes into an array
      /*
      Iterator chromitr = hmchromsize.keySet().iterator();
      String[] chroms = new String[hmchromsize.keySet().size()];
      int nchromindex = 0;
      while (chromitr.hasNext())
      {
      chroms[nchromindex++] = (String) chromitr.next();
      }
      Arrays.sort(chroms);
      */

      for (int nchromindex = 0; nchromindex < chroms.length; nchromindex++) {
        String szchrom = (String) chroms[nchromindex];

        String szfile = (String) hmchromdnamethylfile.get(szchrom);
        BufferedReader brdnamethyldata = Util.getBufferedReader(szmethylDIR + "/" + szfile);

        int intobjchromsize = (Integer) hmchromsize.get(szchrom);
        int nchromsize = ((Integer) intobjchromsize).intValue();
        String szLine;

        while ((szLine = brdnamethyldata.readLine()) != null) {
          // first value is the coordinate
          StringTokenizer stdnamethyl = new StringTokenizer(szLine, "\t ");
          int ncoord = Integer.parseInt(stdnamethyl.nextToken()) - 1;
          int nbin = ncoord / nresolution;

          // chrom size 50
          // array index 0 to 49
          // if 20 is distance can't start from 0 to 19 or 30 to 49

          if ((nbin >= nmaxparam) && (nbin < nchromsize - nmaxparam)) {
            // valid position
            if ((nsampleindex < samples.length) && (noverallindex == samples[nsampleindex])) {
              // we are on a position we wanted to sample
              // storing the bin and chromosome corresponding to it

              if ((szchromwantgenerate != null) && (!szchromwantgenerate.equals(szchrom))) {
                samplesdnamethylvals_nbag[nsampleindex] = null;
                samplesdnamethylbinindex_nbag[nsampleindex] = -1;
                // if in parallel then this will become new chromosome
                samplesdnamethylchromindex_nbag[nsampleindex] = -1; // nchromindex;
              } else {
                samplesdnamethylbinindex_nbag[nsampleindex] = nbin;
                // if in parallel then this will become new chromosome
                if (szchromwantgenerate != null) {
                  samplesdnamethylchromindex_nbag[nsampleindex] = 0; // nchromindex;
                } else {
                  samplesdnamethylchromindex_nbag[nsampleindex] = nchromindex;
                }

                samplesdnamethylvals_nbag[nsampleindex] = new float[numdnamethylcells];
                // copying the methylation values into the sample replacing negative sample values
                // with the average
                float[] samplesdnamethylvals_nbag_nsampleindex =
                    samplesdnamethylvals_nbag[nsampleindex];

                for (int nk = 0; nk < samplesdnamethylvals_nbag_nsampleindex.length; nk++) {
                  float fval = Float.parseFloat(stdnamethyl.nextToken());
                  if (fval < 0) {
                    fval = (float) methylavg[nk];
                  }
                  samplesdnamethylvals_nbag_nsampleindex[nk] =
                      ((int) (nroundval * NMETHYLSCALE * fval + 0.5)) / (froundval);
                }
              }
              nsampleindex++; // we had a hit with this position increment sample count
            }
            noverallindex++; // this a position that counts to overall
          }
        }
        brdnamethyldata.close();

        if ((szchromwantgenerate != null) && (szchromwantgenerate.equals(szchrom))) {
          while (nsampleindex < samplesdnamethylbinindex_nbag.length) {
            samplesdnamethylbinindex_nbag[nsampleindex] = -1;
            samplesdnamethylchromindex_nbag[nsampleindex] = -1;
            nsampleindex++;
          }
          break;
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  static class TrainFileRec {
    String szprefix;
    BufferedReader br;

    public TrainFileRec(String szprefix, BufferedReader br) {
      this.szprefix = szprefix;
      this.br = br;
    }
  }

  static class TrainFileRecCompare implements Comparator, Serializable {

    public int compare(Object o1, Object o2) {
      TrainFileRec r1 = (TrainFileRec) o1;
      TrainFileRec r2 = (TrainFileRec) o2;
      return (r1.szprefix.compareTo(r2.szprefix));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** */
  public void executeTrain() throws Exception {

    NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
    nf.setMaximumFractionDigits(2);
    nf.setGroupingUsed(false);

    // Check that the target mark exists and put it as ntargetmark
    ntargetmark = -1;
    for (int nmark = 0; nmark < marksA.length; nmark++) {
      if (marksA[nmark].equals(szoutmark)) {
        ntargetmark = nmark;
        break;
      }
    }

    if ((!bdnamethyl) && (ntargetmark == -1)) {
      throw new IllegalArgumentException(szoutmark + " was not found as a mark!");
    }

    // Check that the target cell exists and put it as ntargetcell
    if (szoutcell != null) {
      ntargetcell = -1;
      for (int ncell = 0; ncell < cellsA.length; ncell++) {
        if (cellsA[ncell].equals(szoutcell)) {
          ntargetcell = ncell;
          break;
        }
      }

      if (ntargetcell == -1) {
        throw new IllegalArgumentException(szoutcell + " was not found as a cell!");
      }
    }

    // nholdoutcellrequest defaults to -1: TODO figure out what exactly it does for us
    if ((nholdoutcellrequest != -1)
        && (!bdnamethyl && (!bmarkcell[ntargetmark][nholdoutcellrequest]))) {
      throw new IllegalArgumentException(
          "No data available for\tmark " + ntargetmark + "\t cell " + nholdoutcellrequest);
    }

    // getting list of marks indices target cell

    // attributes = new ArrayList[numcells*numbags];

    /*
    int nclassifierindex = 0;
    for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++)
    {
    for (int nbag = 0; nbag < numbags; nbag++)
    {
    if (((nbagrequest == -1)||(nbagrequest ==nbag))&&
    ((nholdoutcellrequest == -1)||((nholdoutcellrequest == nholdoutcell)))&&
    (((bdnamethyl&&bdnamethylcell[nholdoutcell])||(!bdnamethyl&&bmarkcell[ntargetmark][nholdoutcell]))))
    {
    theTrainInstances[nclassifierindex] = new ArrayList();
    theTrainInstancesOutput[nclassifierindex] = new ArrayList();
    }
    nclassifierindex++;
    }
    }
    */

    // loads from training data
    // nclassifierindex = 0;

    // TODO printing to debug as well
    // System.out.println("IN LOADFILE\t" + nbagrequest + "\t" + nholdoutcellrequest
            // + "\t" + numcells + "\t" + ntargetcell + "\t");
            //+ btrain);
    // updates numbags based on nbagrequest amount since we know it needs to be larger
    if (nbagrequest >= numbags) { numbags = nbagrequest + 1; }

    for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
      if (nholdoutcell != ntargetcell) {
        for (int nbag = 0; nbag < numbags; nbag++) {
          if (((nbagrequest == -1) || (nbagrequest == nbag))
              && ((nholdoutcellrequest == -1) || (nholdoutcellrequest == nholdoutcell))
              && (((bdnamethyl && bdnamethylcell[nholdoutcell])
                  || (!bdnamethyl && bmarkcell[ntargetmark][nholdoutcell]))))
          // (bdnamethyl||bmarkcell[ntargetmark][nholdoutcell])&&(theTrainInstances[nclassifierindex]!=null))
          {
            theTrainInstances = new ArrayList(); // [numcells*numbags];
            theTrainInstancesOutput = new ArrayList(); // [numcells*numbags];

            // BufferedReader brtraindata;
            BufferedReader brattributes;

            BufferedReader[] brtraindataA;
            // TODO what happens when we need to open a chromosome instead?
            String infile = szinputdir + "/traindata_" + szoutmark // +theTargetRec.szmark
                + "_" + nholdoutcell + "_" + nbag + ".txt.gz";
            File f = new File(infile);
            // System.out.println("Checking for file " + infile);

            if (f.exists()) {
              brtraindataA = new BufferedReader[1];
              brtraindataA[0] = Util.getBufferedReader(infile);
              // System.out.println("File exists and was read.");
            } else {
              // Check if directory is there:
              File dir = new File(szinputdir);
              String[] files = dir.list();

              // TODO issue if some files are by chromosome and some are not!
              if (files == null) {
                throw new IllegalArgumentException("Directory " + szinputdir + " was not found!");
              }
              ArrayList albrtraindataRecs = new ArrayList();
              String suffix = "_traindata_" + szoutmark + "_" + nholdoutcell + "_" + nbag + ".txt.gz";
              // System.out.println("Looking for files with suffix " + suffix);
              // Look at file listing for suffix matches
              for (int nk = 0; nk < files.length; nk++) {
                if (files[nk].endsWith(suffix)) {
                  int nindexpos = files[nk].indexOf(suffix);
                  String szprefix = files[nk].substring(0, nindexpos);
                  // TODO Print when we have found a file that works:
                  //System.out.println("prefix is " + szprefix + " for file " + files[nk]);

                  // going to sort so
                  albrtraindataRecs.add(
                      new TrainFileRec(
                          szprefix, Util.getBufferedReader(szinputdir + "/" + files[nk])));
                }
              }

              TrainFileRec[] brtraindataRecsA = new TrainFileRec[albrtraindataRecs.size()];
              for (int nk = 0; nk < brtraindataRecsA.length; nk++) {
                brtraindataRecsA[nk] = (TrainFileRec) albrtraindataRecs.get(nk);
              }
              Arrays.sort(brtraindataRecsA, new TrainFileRecCompare());
              brtraindataA = new BufferedReader[brtraindataRecsA.length];
              for (int nk = 0; nk < brtraindataRecsA.length; nk++) {
                brtraindataA[nk] = brtraindataRecsA[nk].br;
              }
            }

            String attrfile = szinputdir + "/attributes_" + szoutmark // +theTargetRec.szmark
                        + "_" + nholdoutcell + "_" + nbag + ".txt.gz";
            // System.out.println("attributes file is " + attrfile);
            brattributes = Util.getBufferedReader(attrfile);

            String useattrfile = szoutdir + "/useattributes_" + szoutcell + "_" + szoutmark // theTargetRec.szcell+"_"+theTargetRec.szmark
                            + "_" + nholdoutcell + "_" + nbag + ".txt.gz";
            // System.out.println("useattributes file is " + useattrfile);
            GZIPOutputStream pwuseattributes;
            pwuseattributes =
                new GZIPOutputStream(
                    new FileOutputStream(useattrfile));

            // System.out.println("OPENING\t" + szinputdir + "/traindata_" + szoutcell + "_" + szoutmark // theTargetRec.szcell+"_"+theTargetRec.szmark
            //         + "_" + nholdoutcell + "_" + nbag + ".txt.gz");
            // System.out.println("VALS\t" + ntargetcell + "\t" + nholdoutcell + "\t" + nholdoutcellrequest + "\t" + nbag + "\t" + nbagrequest);
            // loadInstancesFromFile(brtraindata, brattributes, pwuseattributes, ntargetcell);
            loadInstancesFromFile(brtraindataA, brattributes, pwuseattributes, ntargetcell);
            for (int nk = 0; nk < brtraindataA.length; nk++) {
              brtraindataA[nk].close();
            }
            pwuseattributes.finish();
            pwuseattributes.close();

            // if (BLINEAR) {
            //          theClassifierLinearA = new RegressionLinear[theTrainInstances.length];
            // } else {
            //          theClassifierA = new RegressionTree[theTrainInstances.length];
            // }

            if (((float[]) theTrainInstances.get(0)).length >= 1) {
              // only classify if attributes available
              if (BLINEAR) {
                //// System.out.println("getting theTrainInstances "+nclassifier+" "+nholdoutcell+"
                // "+nbag);
                //// if ((theTrainInstances[nclassifier] != null)&&(((float[])
                // theTrainInstances[nclassifier].get(0)).length >=1))
                //// {

                RegressionLinear theClassifier;
                double[][] data =
                    new double[theTrainInstances.size()]
                        [1 + ((float[]) theTrainInstances.get(0)).length];
                double[][] y = new double[data.length][1];
                for (int na = 0; na < data.length; na++) {
                  float[] al_na = (float[]) theTrainInstances.get(na);
                  double[] data_na = data[na];
                  for (int nb = 0; nb < al_na.length; nb++) {
                    data_na[nb + 1] = (double) al_na[nb];
                  }
                  data_na[0] = 1;
                }

                // ArrayList theTrainInstancesOutput_nclassifier =
                // theTrainInstancesOutput[nclassifier];
                for (int na = 0; na < y.length; na++) {
                  // y[na][0] = ((Float) theTrainInstancesOutput_classifier.get(na)).floatValue();
                  y[na][0] = ((Float) theTrainInstancesOutput.get(na)).floatValue();
                }

                double[] beta = null;
                /*
                //UNCOMMENT FOR LINEAR REGRESSION
                Matrix m_data = new Matrix(data);
                Matrix m_y = new Matrix(y);
                double dridge = 1;
                //double dridge = 10000;

                LinearRegression lr = new LinearRegression(m_data,m_y,dridge);
                beta = lr.getCoefficients();
                */

                // if ((nholdoutcellrequest !=-1)|| (nbagrequest != -1))
                {
                  GZIPOutputStream pw =
                      new GZIPOutputStream(
                          new FileOutputStream(
                              szoutdir
                                  + "/linearclassifier_"
                                  + szoutcell
                                  + "_"
                                  + szoutmark
                                  + "_"
                                  + nholdoutcell
                                  + "_"
                                  + nbag
                                  + ".txt.gz"));

                  StringBuffer sb = new StringBuffer();
                  for (int nk = 0; nk < beta.length; nk++) sb.append(beta[nk] + "\n");
                  byte[] btformat = sb.toString().getBytes();
                  pw.write(btformat, 0, btformat.length);
                  pw.finish();
                  pw.close();
                }

              } else {
                // TODO CONTINUED HERE (not LINEAR REGRESSION)
                // System.out.println(
                //         "getting theTrainInstances " + nclassifier + " " + nholdoutcell + " " + nbag);
                // if ((theTrainInstances[nclassifier] != null)&&(((float[])
                // theTrainInstances[nclassifier].get(0)).length >=1))
                // {
                // System.out.println(nclassifier+" before train
                // "+theTrainInstances[nclassifier].size());

                // System.out.println("here before regression tree classifier");
                // not doing equivalent for linear
                RegressionTree theClassifier =
                    new RegressionTree(
                        theTrainInstances, theTrainInstancesOutput, nminnumlocations);

                String classifierfile = szoutdir + "/classifier_" + szoutcell + "_" + szoutmark
                        + "_" + nholdoutcell + "_" + nbag + ".txt.gz";
                // System.out.println("Printing theClassifier to file " + classifierfile);
                // not doing equivalent for linear
                GZIPOutputStream pw =
                    new GZIPOutputStream(
                        new FileOutputStream(classifierfile));
                byte[] btformat = theClassifier.toString().getBytes();
                pw.write(btformat, 0, btformat.length);
                pw.finish();
                pw.close();
              }
            }
          }
          // nclassifierindex++;
        }
      }
      // else
      // {
      //  nclassifierindex += numbags;
      // }

    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** */
  public void executeApply() throws Exception {

    NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
    nf.setMaximumFractionDigits(2);
    nf.setGroupingUsed(false);

    ntargetmark = -1;
    for (int nmark = 0; nmark < marksA.length; nmark++) {
      if (marksA[nmark].equals(szoutmark)) {
        ntargetmark = nmark;
        break;
      }
    }

    if ((!bdnamethyl) && (ntargetmark == -1)) {
      throw new IllegalArgumentException(szoutmark + " was not found as a mark!");
    }

    if (szoutcell != null) {
      ntargetcell = -1;
      for (int ncell = 0; ncell < cellsA.length; ncell++) {
        if (cellsA[ncell].equals(szoutcell)) // theTargetRec.szcell))
        {
          ntargetcell = ncell;
          break;
        }
      }

      if (ntargetcell == -1) {
        throw new IllegalArgumentException(
            szoutcell + " was not found as a cell!"); // theTargetRec.szcell+" was not found!");
      }
    }

    // if ((nholdoutcellrequest != -1) &&
    // (!bdnamethyl&&(!bmarkcell[ntargetmark][nholdoutcellrequest])))
    // {
    //   throw new IllegalArgumentException("No data available for\tmark "+ntargetmark+"\t cell
    // "+nholdoutcellrequest);
    // }

    // getting list of marks indicies target cell
    ArrayList almarkstargetcell = new ArrayList();
    for (int nmark = 0; nmark < bmarkcell.length; nmark++) {
      if (bmarkcell[nmark][ntargetcell]) {
        almarkstargetcell.add(Integer.valueOf(nmark));
      }
    }

    // theTrainInstances = new ArrayList[numcells*numbags];
    // theTrainInstancesOutput = new ArrayList[numcells*numbags];
    // theTestInstances  = new ArrayList[numcells*numbags];
    attributes = new ArrayList[numcells * numbags];
    // trainPW = new GZIPOutputStream[numcells*numbags];

    loadClassifiers();

    if (numclassifiers == 0) {
      throw new IllegalArgumentException(
          "No previously trained classifiers for mark "
              + szoutmark + " were found available to load!");
    }

    // System.out.println("Target marks");
    int[] markstargetcellA = new int[almarkstargetcell.size()];
    for (int nmark = 0; nmark < markstargetcellA.length; nmark++) {
      markstargetcellA[nmark] = ((Integer) almarkstargetcell.get(nmark)).intValue();
    }

    int nclassifierindex = 0;
    for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
      for (int nbag = 0; nbag < numbags; nbag++) {
        // theTestInstances[nclassifierindex] = new ArrayList();
        nclassifierindex++;
      }
    }

    System.out.println("Allocating space: " + ChromImpute.NUMLINES + " x " + nummarks + " x " +  numcells );
    float[][][] databinfirst = new float[ChromImpute.NUMLINES][nummarks][numcells]; // bin, mark, then cell

    System.out.println("Allocating space: " + nummarks + " x " + nummarks + " x " +ChromImpute.NUMLINES );
    float[][][] databinlast = new float[nummarks][numcells][ChromImpute.NUMLINES]; // mark, cell, then bin

    System.out.println("Allocating space: " + numcells + " x " + nummarks);
    float[][][] databinlastCELLMARK = new float[numcells][nummarks][];

    // stores for each mark and cell total signal and number of entries
    BufferedReader[][] brA = new BufferedReader[marksA.length][cellsA.length];

    // System.out.println("bnormalize is "+bnormalize);

    // getting chroms into array array
    Iterator chromitr = hmchromsize.keySet().iterator();
    String[] chroms = new String[hmchromsize.keySet().size()];
    int nchromindex = 0;
    while (chromitr.hasNext()) {
      chroms[nchromindex++] = (String) chromitr.next();
    }
    Arrays.sort(chroms);

    if (bprintonefile) {
      System.out.println("Will print all chromosomes to " + szoutdir + "/" + szoutfile + ".gz");
      pwout = new GZIPOutputStream(new FileOutputStream(szoutdir + "/" + szoutfile + ".gz"));
      if (bprintbrowserheader) {
        String szHeader = "track type=wiggle_0 name=" + szoutcell + "_" + szoutmark + "_imputed\n";

        byte[] btformat = szHeader.getBytes();
        pwout.write(btformat, 0, btformat.length);
      }
    }

    for (int nchrom = 0; nchrom < chroms.length; nchrom++) {
      long ltime = System.currentTimeMillis();
      // iterating through each chromosome generating training instances for that chromosome
      String szchrom = chroms[nchrom];
      System.out.println("Adding chromosome " + szchrom);

      // TODO This checks if the chromosome is already completed.
      // TODO do gzip -t test on chromfiles before we run.
      String chromfile = szoutdir + "/" + szchrom + "_" + szoutfile + ".gz";
      File chrf = new File(chromfile);
      if (chrf.exists()) {
        System.out.println("SKIPPING: " + chromfile + " - already exists.");
      } else {
        ArrayList almethylvals = null;
        ArrayList almethylcoord = null;
        int[] methylindexA = null;

        // will apply the classifiers
        if (!bprintonefile) {
          System.out.println("COMPUTING: Will print to " + chromfile);
          pwout = 
              new GZIPOutputStream(
                  new FileOutputStream(chromfile));
        }

        // prints header lines for data
        String szHeader = "";
        if ((bprintbrowserheader) && (!bprintonefile)) {
          szHeader = "track type=wiggle_0 name=" + szoutcell + "_" + szoutmark + "_imputed\n";
        }
        if (bdnamethyl) {
          szHeader += "variableStep chrom=" + szchrom + "\n";
        } else {
          szHeader += "fixedStep chrom=" + szchrom + " start=1 step=" + nresolution
              + " span=" + nresolution + "\n";
        }
        byte[] btformat = szHeader.getBytes();
        pwout.write(btformat, 0, btformat.length);

        if (bdnamethyl) {
          String szfile = (String) hmchromdnamethylfile.get(szchrom);
          if (szfile == null) {
            throw new IllegalArgumentException(
                "Chromosome " + szchrom + " was not found in the DNA methylation info file");
          }
          BufferedReader brdnamethyldata =
              Util.getBufferedReader(
                  szmethylDIR + "/" + szfile); // assumes methylation data in separate files
          // int intobjchromsize = (Integer) hmchromsize.get(szchrom);
          String szLine;

          almethylvals = new ArrayList();
          almethylcoord = new ArrayList();
          methylindexA = new int[1];
          methylindexA[0] = 0; // used to pass by reference for methylation data

          while ((szLine = brdnamethyldata.readLine()) != null) {
            // first value is the coordinate
            StringTokenizer stdnamethyl = new StringTokenizer(szLine, "\t ");
            int ncoord =
                Integer.parseInt(stdnamethyl.nextToken()) - 1; // converts coordinate to 0 base
            almethylcoord.add(Integer.valueOf(ncoord)); //

            float[] methylvals =
                new float[numdnamethylcells]; // allocates an array with the dna methylati val

            almethylvals.add(methylvals);

            // chrom size 50
            // array index 0 to 49
            // if 20 is distance can't start from 0 to 19 or 30 to 49

            // going through the methylation values
            for (int nk = 0; nk < methylvals.length; nk++) {
              float fval = Float.parseFloat(stdnamethyl.nextToken());
              if (fval < 0) {
                // methylation values less than 0 get replaced with average
                fval = (float) methylavg[nk];
              }

              // rounding methylation values to nearest tength after multiplying by 100
              methylvals[nk] = ((int) (nroundval * NMETHYLSCALE * fval + 0.5)) / (froundval);
            }
          }
          brdnamethyldata.close();
        }

        // int nchromsize = ((Integer) hmchromsize.get(szchrom)).intValue()/nresolution+1;
        int nchromsize = (((Integer) hmchromsize.get(szchrom)).intValue() - 1) / nresolution + 1;

        float[] predictvals = null;
        boolean[] presentvals = null;

        if (bdnamethyl) {
          predictvals = new float[almethylvals.size()]; // number of positions of dna methylation
          presentvals =
              new boolean
                  [predictvals
                      .length]; // stores whether a prediction was actually made for the position or
          // if it should be treated as missing
        } else {
          predictvals = new float[nchromsize];
        }

        // allocating memory
        System.out.println("Allocating memory");
        System.out.println("-- Loop is " + marksA.length + " x " + cellsA.length);
        for (int nmark = 0; nmark < marksA.length; nmark++) {
          for (int ncell = 0; ncell < cellsA.length; ncell++) {
            String szinfile = (String) hmmarkcellfile.get(marksA[nmark] + "\t" + cellsA[ncell]);

            if (szinfile != null) {
              // System.out.println("READING IN " + szinputdir + "/" + szchrom + "_" + szinfile + ".wig.gz");
              brA[nmark][ncell] =
                  Util.getBufferedReader(szinputdir + "/" + szchrom + "_" + szinfile + ".wig.gz");

              brA[nmark][ncell].readLine();
              brA[nmark][ncell].readLine();
              // if (databinfirst == null)
              // {
              //   databinfirst = new float[ChromImpute.NUMLINES][nummarks][numcells];
              //   databinlast = new float[nummarks][numcells][ChromImpute.NUMLINES];
              //   databinlastCELLMARK = new float[numcells][nummarks][];
              // }
              databinlastCELLMARK[ncell][nmark] = databinlast[nmark][ncell];
            }
          }
        }

        int nmaxPARAM = Math.max(nknnoffset, nmaxoffsetwide);
        int n2nmaxoffset = 2 * nmaxPARAM;

        ltime = System.currentTimeMillis();
        // reading in the initial portion of the data
        for (int nbeginread = 0; nbeginread < n2nmaxoffset; nbeginread++) {
          float[][] databinfirst_nbeginread = databinfirst[nbeginread];
          for (int nmark = 0; nmark < marksA.length; nmark++) {
            float[] databinfirst_nbeginread_nmark = databinfirst_nbeginread[nmark];
            float[][] databinlast_nmark = databinlast[nmark];
            boolean[] bmarkcell_nmark = bmarkcell[nmark];
            BufferedReader[] brA_nmark = brA[nmark];

            for (int ncell = 0; ncell < cellsA.length; ncell++) {
              if (bmarkcell_nmark[ncell]) {
                String szLine = brA_nmark[ncell].readLine();
                if (szLine == null) break;

                // rounding value here so will be consistent with format of data from loaded instances
                float fval = ((int) (nroundval * Float.parseFloat(szLine) + 0.5)) / froundval;
                databinfirst_nbeginread_nmark[ncell] = fval;
                databinlast_nmark[ncell][nbeginread] = fval;
              }
            }
          }
        }

        // System.out.println("B "+(System.currentTimeMillis()-ltime));
        // ltime = System.currentTimeMillis();
        int nbin = 0;

        // System.out.println("Reading files:");
        for (int nbeginread = n2nmaxoffset;
            nbeginread < nchromsize;
            nbeginread += (ChromImpute.NUMLINES - n2nmaxoffset)) {
            // going to read out the rest of the chromosome in chunks
            // System.out.println("Reading next chunk:");

            boolean bok = true;
            ltime = System.currentTimeMillis();
            for (nbin = n2nmaxoffset; nbin < ChromImpute.NUMLINES && bok; nbin++) {
                float[][] databinfirst_nbin = databinfirst[nbin];
                for (int nmark = 0; nmark < marksA.length; nmark++) {
                    // System.out.println("here 3");
                    float[] databinfirst_nbin_nmark = databinfirst_nbin[nmark];
                    float[][] databinlast_nmark = databinlast[nmark];
                    boolean[] bmarkcell_nmark = bmarkcell[nmark];
                    BufferedReader[] brA_nmark = brA[nmark];

                    String szLine;

                    for (int ncell = 0; ncell < cellsA.length; ncell++) {
                        if (bmarkcell_nmark[ncell]) {
                            if ((szLine = brA_nmark[ncell].readLine()) != null) {
                                // rounding value here so will be consistent with format of data from loaded
                                // instances
                                float fval = ((int) (nroundval * Float.parseFloat(szLine) + 0.5)) / froundval;
                                databinfirst_nbin_nmark[ncell] = fval;
                                databinlast_nmark[ncell][nbin] = fval;
                            } else {
                                bok = false;
                            }
                        }
                    }
                }
            }

          // added in to update way handling last prediction position on chromosome
          if ((!bdnamethyl) && (!bok)) {
            // comment out to revert to old way
            nbin--;
          }
          // System.out.println("CALLING TEST
          // with\t"+nbin+"\t"+nbeginread+"\t"+n2nmaxoffset+"\t"+(nbeginread-n2nmaxoffset));

          // System.out.println("Processing the chunks read");
          // process the chunk read
          if (bdnamethyl) {
            generateInstanceDataTestDNAMethyl(
                databinlast,
                databinlastCELLMARK,
                databinfirst,
                nbin,
                nbeginread - n2nmaxoffset,
                // bmarkcell, bcellmark,
                markstargetcellA,
                ntargetcell,
                // false,
                predictvals,
                presentvals,
                almethylvals,
                almethylcoord,
                methylindexA);
          } else {
            generateInstanceDataTest( // ERROR IS HERE
                databinlast,
                databinlastCELLMARK,
                databinfirst,
                nbin,
                nbeginread - n2nmaxoffset,
                markstargetcellA,
                ntargetcell,
                ntargetmark,
                predictvals);
          }

          // System.out.println("Copy last elements to the front");
          // copies the last n2nmaxoffset elements to the front
          for (int ni = 0; ni < n2nmaxoffset; ni++) {
            float[][] databinfirst_ni = databinfirst[ni];
            float[][] databinfirst_nend = databinfirst[databinfirst.length - n2nmaxoffset + ni];
            for (int nmark = 0; nmark < databinfirst_ni.length; nmark++) {
              float[][] databinlast_nmark = databinlast[nmark];
              float[] databinfirst_ni_nmark = databinfirst_ni[nmark];
              float[] databinfirst_nend_nmark = databinfirst_nend[nmark];
              for (int ncell = 0; ncell < databinfirst_ni_nmark.length; ncell++) {
                // copys the contents
                float fval = databinfirst_nend_nmark[ncell];
                databinfirst_ni_nmark[ncell] = fval;
                databinlast_nmark[ncell][ni] = fval;
              }
            }
          }
          // System.out.println("E "+(System.currentTimeMillis()-ltime));
        }
        // }

        for (int na = 0; na < predictvals.length; na++) {
          if (bdnamethyl) {
            String szval;
            int nmethylcoord = ((Integer) almethylcoord.get(na)).intValue() + 1;
            if (presentvals[na]) {
              // we are making a prediction for the dna methylation
              // also have to rescale based on NMETHYLSCALE value
              szval =
                  nmethylcoord
                      + "\t"
                      + nf.format(predictvals[na] / (NMETHYLSCALE * numclassifiers))
                      + "\n";
            } else {
              // going to use default edge val
              szval = nmethylcoord + "\t" + nf.format(DEDGEVAL) + "\n";
            }
            btformat = szval.getBytes();
            pwout.write(btformat, 0, btformat.length);
          } else {
            // prints the average
            String szval = nf.format(predictvals[na] / numclassifiers) + "\n";
            btformat = szval.getBytes();
            pwout.write(btformat, 0, btformat.length);
          }
        }

        // generates data for the current chromosome

        for (int na = 0; na < brA.length; na++) {
          for (int nb = 0; nb < brA[na].length; nb++) {
            if (brA[na][nb] != null) {
              brA[na][nb].close();
            }
          }
        }

        if (!bprintonefile) {
          pwout.finish();
          pwout.close();
        }
      }
    }

    if (bprintonefile) {
      pwout.finish();
      pwout.close();
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * Loads information on the data imputation File is of the form input: load cell type, load mark,
   * data file
   */
  public void loadImputeInfo() throws IOException {
    // stores the list of all marks
    hsmarks = new HashSet();

    // stores the list of all cell types involved in the imputation
    hscells = new HashSet();

    // stores for each mark an array list with all cell types associated with the mark
    hmmarkcell = new HashMap();

    // stores for a mark-cell combination the associated input data file
    hmmarkcellfile = new HashMap();

    //
    BufferedReader brimputeinfo = Util.getBufferedReader(szimputeinfoINfile);
    String szLine;

    while ((szLine = brimputeinfo.readLine()) != null) {
      StringTokenizer st = new StringTokenizer(szLine, "\t");
      if (st.countTokens() == 0) continue;
      // reading cell and mark in the file
      String szcell = st.nextToken();
      String szmark = st.nextToken();

      if (szmark.startsWith("OrderNarrow_"))
        throw new IllegalArgumentException(
            "Mark "
                + szmark
                + " name starts with 'OrderNarrow_' which is a reserved prefix, please change the mark name.");
      if (szmark.startsWith("OrderGlobal_"))
        throw new IllegalArgumentException(
            "Mark "
                + szmark
                + " name starts with 'OrderGlobal_' which is a reserved prefix, please change the mark name.");
      if (szmark.contains("_by_"))
        throw new IllegalArgumentException(
            "Mark "
                + szmark
                + " contains '_by_' which is a reserved substring, please change the mark name.");

      // skips an entry if szpioneermark is not null
      // hspioneer does not contains the mark
      // and this is for the target out cell
      if ((szpioneermark != null)
          && (!hspioneermarks.contains(szmark))
          && (szcell.equals(szoutcell))) continue;

      // gets the input file
      String szinfile = st.nextToken();

      // adding mark to the set of marks available
      hsmarks.add(szmark);

      // adding mark to the set of cell types
      hscells.add(szcell);

      ArrayList alcells = (ArrayList) hmmarkcell.get(szmark);
      if (alcells == null) {
        // first time with this mark, creates arraylist associated with it
        alcells = new ArrayList();
        hmmarkcell.put(szmark, alcells);
      }
      // adding cell to the list of cell
      alcells.add(szcell);

      // maps a mark cell combination to a file
      hmmarkcellfile.put(szmark + "\t" + szcell, szinfile);

      // reads in all the mark-cell-data file information
    }
    brimputeinfo.close();

    // stores the count of the number of marks and cells
    nummarks = hsmarks.size();
    numcells = hscells.size();

    // updates the max-nearest neighbor based on the number of cells
    nmaxknn = Math.min(nmaxknn, numcells - 1);

    // stores the marks and cells into the marksA and cellsA array
    // and then sorts them
    marksA = new String[nummarks];
    cellsA = new String[numcells];
    Iterator itrcells = (Iterator) hscells.iterator();
    Iterator itrmark = (Iterator) hsmarks.iterator();

    for (int nmark = 0; nmark < marksA.length; nmark++) {
      marksA[nmark] = (String) itrmark.next();
    }

    for (int ncell = 0; ncell < cellsA.length; ncell++) {
      cellsA[ncell] = (String) itrcells.next();
    }

    Arrays.sort(marksA);
    Arrays.sort(cellsA);

    if (bdnamethyl) {
      regularcelltodnamethylindex =
          new int[numcells]; // contains a mapping to the methyl cell index if being used
      bdnamethylcell = new boolean[numcells]; // true if the cell is being used
      for (int nregularcellindex = 0; nregularcellindex < cellsA.length; nregularcellindex++) {
        // bdnamethylcell[na] = false;
        for (int nmethylcellindex = 0;
            nmethylcellindex < dnamethylheader.length;
            nmethylcellindex++) {
          if (cellsA[nregularcellindex].equals(dnamethylheader[nmethylcellindex])) {
            regularcelltodnamethylindex[nregularcellindex] = nmethylcellindex;
            bdnamethylcell[nregularcellindex] = true;
            break;
          }
        }
      }
    }

    // stores for each mark and cell type combination if it is present
    bmarkcell = new boolean[nummarks][numcells];
    bcellmark = new boolean[numcells][nummarks];

    for (int nmark = 0; nmark < marksA.length; nmark++) {
      // iterates through each mark and marks which cell types the mark has been mapped in
      String szmark = marksA[nmark];
      // marksA[nmark] = szmark;
      ArrayList alcells = (ArrayList) hmmarkcell.get(szmark);

      int nalcellsize = alcells.size();

      boolean[] bmarkcell_nmark = bmarkcell[nmark];
      for (int nindex = 0; nindex < nalcellsize; nindex++) {
        // gets the corresponding cell index
        int njindex = Arrays.binarySearch(cellsA, alcells.get(nindex));
        bmarkcell_nmark[njindex] = true;
        bcellmark[njindex][nmark] = true;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** */
  public void loadDistInfo() throws IOException {
    //
    distmarkcellAL = new ArrayList[nummarks][numcells];
    HashMap hmCellToID = new HashMap();

    // inititalizes for each mark, and each cell an array list to have in sorted order the closest
    // cells
    for (int nmark = 0; nmark < distmarkcellAL.length; nmark++) {
      for (int ncell = 0; ncell < numcells; ncell++) {
        distmarkcellAL[nmark][ncell] = new ArrayList();
      }
    }

    // maps each cell identifer to index
    for (int ncell = 0; ncell < cellsA.length; ncell++) {
      hmCellToID.put(cellsA[ncell], Integer.valueOf(ncell));
    }

    for (int nmark = 0; nmark < bmarkcell.length; nmark++) {
      // added for RNA-seq don't need the target mark one since not using
      if ((szoutmark != null) && (szoutmark.equals(marksA[nmark]))) continue;

      for (int ncell = 0; ncell < bmarkcell[nmark].length; ncell++) {
        if (bmarkcell[nmark][ncell]) {
          // we have the mark in the cell type
          BufferedReader brdist =
              new BufferedReader(
                  new FileReader(
                      szdistancedir + "/" + cellsA[ncell] + "_" + marksA[nmark] + ".txt"));
          String szLine;

          // reading in the nearest other cell types in order
          while ((szLine = brdist.readLine()) != null) {
            StringTokenizer st = new StringTokenizer(szLine, "\t");
            String szcell = st.nextToken();
            Integer objID = (Integer) hmCellToID.get(szcell);
            if (objID != null) {
              // adding the cell identifier into the list
              distmarkcellAL[nmark][ncell].add(objID);
            }
          }
          brdist.close();
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////

  public void executeGenerateTraining() throws Exception {

    NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
    nf.setMaximumFractionDigits(2);
    nf.setGroupingUsed(false);

    ntargetmark = -1;
    for (int nmark = 0; nmark < marksA.length; nmark++) {
      if (marksA[nmark].equals(szoutmark)) // theTargetRec.szmark))
      {
        ntargetmark = nmark;
        break;
      }
    }

    if ((!bdnamethyl) && (ntargetmark == -1)) {
      throw new IllegalArgumentException(szoutmark + " was not found as mark!");
    }

    /*
    if (szoutcell != null)//theTargetRec.szcell != null) {
        ntargetcell = -1;
        for (int ncell = 0; ncell < cellsA.length; ncell++) {
            if (cellsA[ncell].equals(szoutcell)) //theTargetRec.szcell)) {
                ntargetcell = ncell;
                break;
            }
        }

        if (ntargetcell == -1) {
            throw new IllegalArgumentException(
                szoutcell+" was not found!"); //theTargetRec.szcell+" was not found!");
        }
    }


    if ((nholdoutcellrequest != -1) && (!bdnamethyl&&(!bmarkcell[ntargetmark][nholdoutcellrequest])))
    {
    throw new IllegalArgumentException("No data available for\tmark "+ntargetmark+"\t cell "+nholdoutcellrequest);
    }
    */

    String[] chroms;

    Iterator chromitr = hmchromsize.keySet().iterator();
    chroms = new String[hmchromsize.keySet().size()];
    int nchromindex = 0;
    while (chromitr.hasNext()) {
      chroms[nchromindex++] = (String) chromitr.next();
    }
    Arrays.sort(chroms);

    if (bdnamethyl) {
      generateSamplesDNAmethyl(chroms);

      if (szchromwantgenerate != null) {
        // System.out.println("before while");
        // paralllelize this
        int nchrom = 0;
        while ((nchrom < chroms.length) && (!szchromwantgenerate.equals(chroms[nchrom]))) {
          int nclassifierindex = 0;
          int nmaxPARAM = Math.max(nknnoffset, nmaxoffsetwide);
          int n2nmaxoffset = 2 * nmaxPARAM;

          // int nchromsize = (((Integer)
          // hmchromsize.get(chroms[nchrom])).intValue()-1)/nresolution+1;

          for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
            // don't need to do this for each held out cell
            for (int nbag = 0; nbag < numbags; nbag++) {
              int[] samplesdnamethylchromindex_nbag = samplesdnamethylchromindex[nbag];
              int nsampleindexA_nclassifierindex = nsampleindexA[nclassifierindex];

              while ((nsampleindexA_nclassifierindex < samplesdnamethylchromindex_nbag.length)
                  && (samplesdnamethylchromindex_nbag[nsampleindexA_nclassifierindex] < 0)) {
                nsampleindexA[nclassifierindex]++;
                // also updating our local variable for it
                nsampleindexA_nclassifierindex++;
              }

              // System.out.println("==> "+nsampleindexA[nclassifierindex]+"
              // "+nsampleindexA[nclassifierindex]+"
              // "+samplesdnamethylchromindex[nbag][nsampleindexA[nclassifierindex]]);
              nclassifierindex++;
            }
          }
          nchrom++;
        }
      }

    } else {
      generateSamples();

      if (szchromwantgenerate != null) {
        int nchrom = 0;
        while ((nchrom < chroms.length) && (!szchromwantgenerate.equals(chroms[nchrom]))) {
          int nclassifierindex = 0;
          int nmaxPARAM = Math.max(nknnoffset, nmaxoffsetwide);
          int n2nmaxoffset = 2 * nmaxPARAM;

          int nchromsize =
              (((Integer) hmchromsize.get(chroms[nchrom])).intValue() - 1) / nresolution + 1;

          for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
            // don't need to do this for each held out cell
            for (int nbag = 0; nbag < numbags; nbag++) {
              int[] samples_nclassifierindex = samples[nclassifierindex];
              int nsampleindexA_nclassifierindex = nsampleindexA[nclassifierindex];

              for (int npos = nmaxPARAM; npos <= nchromsize - 1 - nmaxPARAM; npos++)
              // for (int npos = nmaxPARAM; npos <= nchromsize  - nmaxPARAM; npos++)
              {
                // updated to nchromsize - 1 for backwards compatibility

                while ((nsampleindexA_nclassifierindex < samples_nclassifierindex.length)
                    &&
                    //        (samples_nclassifierindex[nsampleindexA_nclassifierindex] <=
                    // noverallpos[nclassifierindex]))
                    (samples_nclassifierindex[nsampleindexA_nclassifierindex]
                        == noverallpos[nclassifierindex])) {
                  nsampleindexA[nclassifierindex]++;
                  // also updating our local variable for it
                  nsampleindexA_nclassifierindex++;
                }
                // moving which base we have sampled for this classifier
                noverallpos[nclassifierindex]++;
              }
              // System.out.println("===> "+chroms[nchrom]+" "+nclassifierindex+"
              // "+nsampleindexA[nclassifierindex]+" "+noverallpos[nclassifierindex]);
              nclassifierindex++;
            }
          }

          for (int na = 0; na < noverallpos.length; na++) {
            noverallpos[na]++;
            // if (na <= 2)
            // System.out.println("position
            // "+na+"\t"+nchrom+"\t"+noverallpos[na]+"\t"+nsampleindexA[na]);
          }
          nchrom++;
        }
      }
    }

    if (szchromwantgenerate != null) {
      chroms = new String[1];
      chroms[0] = szchromwantgenerate;
    }

    // theTrainInstances = new ArrayList[numcells*numbags];
    // theTrainInstancesOutput = new ArrayList[numcells*numbags];
    // theTestInstances  = new ArrayList[numcells*numbags];

    attributes = new ArrayList[numcells * numbags];
    trainPW = new GZIPOutputStream[numcells * numbags];

    // System.out.println("Target marks");
    /*
    int[] markstargetcellA = new int[almarkstargetcell.size()];
    for (int nmark = 0; nmark < markstargetcellA.length; nmark++)
    {
    markstargetcellA[nmark] = ((Integer) almarkstargetcell.get(nmark)).intValue();
    }
    */

    boolean bprintattributes =
        ((szchromwantgenerate == null) || (szchromwantgenerate.equals((String) alchrom.get(0))));
    int nclassifierindex = 0;

    for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
      for (int nbag = 0; nbag < numbags; nbag++) {
        if ((bdnamethyl && bdnamethylcell[nholdoutcell])
            || (!bdnamethyl && bmarkcell[ntargetmark][nholdoutcell])) {
          // we have the target mark in this cell type will create features for it

          if (bprintattributes) {
            pwattributes =
                new GZIPOutputStream(
                    new FileOutputStream(
                        szoutdir
                            + "/attributes_"
                            + szoutmark // +theTargetRec.szmark
                            + "_"
                            + nholdoutcell
                            + "_"
                            + nbag
                            + ".txt.gz"));
          }

          if (bdnamethyl) {
            generateFeaturesDNAMethyl( // bmarkcell, markstargetcellA,
                nholdoutcell, nclassifierindex, bprintattributes);
          } else {
            generateFeatures( // bmarkcell, markstargetcellA,
                ntargetmark, nholdoutcell, nclassifierindex, bprintattributes);
          }

          if (bprintattributes) {
            pwattributes.finish();
            pwattributes.close();
          }

          // theTrainInstances[nclassifierindex] = new ArrayList();
          // theTrainInstancesOutput[nclassifierindex] = new ArrayList();
          // theTestInstances[nclassifierindex] = new ArrayList();
          // Different file for chromosome segmented or not
          String outfile;
          if (szchromwantgenerate == null) {
            outfile = szoutdir + "/traindata_" + szoutmark // +theTargetRec.szmark
                + "_" + nholdoutcell + "_" + nbag + ".txt.gz";
          } else {
            outfile = szoutdir + "/" + szchromwantgenerate + "_traindata_" 
                + szoutmark + "_" + nholdoutcell + "_" + nbag + ".txt.gz";
          }
          System.out.println("Printing to " + outfile);
          trainPW[nclassifierindex] =
              new GZIPOutputStream(
                      new FileOutputStream(outfile));
        }
        nclassifierindex++;
      }
    }

    // float[][][] databinfirst = null;  //bin, mark, then cell
    // float[][][] databinlast = null;   //mark, cell, then bin
    // float[][][] databinlastCELLMARK = null;

    float[][][] databinfirst = new float[ChromImpute.NUMLINES][nummarks][numcells];
    float[][][] databinlast = new float[nummarks][numcells][ChromImpute.NUMLINES];
    float[][][] databinlastCELLMARK = new float[numcells][nummarks][];

    // stores for each mark and cell total signal and number of entries
    BufferedReader[][] brA = new BufferedReader[marksA.length][cellsA.length];

    for (int nchrom = 0; nchrom < chroms.length; nchrom++) {

      // System.out.println("===> "+chroms.length+" "+nchrom+" "+chroms[nchrom]);//+" "+
      // for (int nk = 0; nk < nsampleindexA.length; nk++)
      // {
      //    System.out.println(nsampleindexA[nk]+" "+noverallpos[nk]);
      // nsampleindexA[0]+" "+noverallpos[0]+" "+nsampleindexA[1]+" "+noverallpos[1]);
      // }
      // long ltime = System.currentTimeMillis();
      // iterating through each chromosome generating training instances for that chromosome
      String szchrom = chroms[nchrom];

      // ArrayList almethylvals = null;
      // ArrayList almethylcoord = null;
      // int[] methylindexA = null;

      // int nchromsize = ((Integer) hmchromsize.get(szchrom)).intValue()/nresolution+1;
      int nchromsize = (((Integer) hmchromsize.get(szchrom)).intValue() - 1) / nresolution + 1;

      // float[] predictvals = null;
      // boolean[] presentvals =null;

      // if (!bloadtrainfile)
      // {

      // generating the features from the raw data
      for (int nmark = 0; nmark < marksA.length; nmark++) {
        for (int ncell = 0; ncell < cellsA.length; ncell++) {
          String szinfile = (String) hmmarkcellfile.get(marksA[nmark] + "\t" + cellsA[ncell]);

          if (szinfile != null) {

            // also means bmarkcell will be true
            brA[nmark][ncell] =
                Util.getBufferedReader(szinputdir + "/" + szchrom + "_" + szinfile + ".wig.gz");

            brA[nmark][ncell].readLine();
            brA[nmark][ncell].readLine();
            // if (databinfirst == null)
            // {
            //   databinfirst = new float[ChromImpute.NUMLINES][nummarks][numcells];
            //   databinlast = new float[nummarks][numcells][ChromImpute.NUMLINES];
            //   databinlastCELLMARK = new float[numcells][nummarks][];
            // }
            databinlastCELLMARK[ncell][nmark] = databinlast[nmark][ncell];
          }
        }
      }

      int nmaxPARAM = Math.max(nknnoffset, nmaxoffsetwide);
      int n2nmaxoffset = 2 * nmaxPARAM;

      // System.out.println("A "+(System.currentTimeMillis()-ltime));
      // ltime = System.currentTimeMillis();

      // reads in the initial portion up to twice the maximum offset
      for (int nbeginread = 0; nbeginread < n2nmaxoffset; nbeginread++) {
        float[][] databinfirst_nbeginread = databinfirst[nbeginread];
        for (int nmark = 0; nmark < marksA.length; nmark++) {
          float[] databinfirst_nbeginread_nmark = databinfirst_nbeginread[nmark];
          float[][] databinlast_nmark = databinlast[nmark];
          boolean[] bmarkcell_nmark = bmarkcell[nmark];
          BufferedReader[] brA_nmark = brA[nmark];

          for (int ncell = 0; ncell < cellsA.length; ncell++) {
            if (bmarkcell_nmark[ncell]) {
              // only read in data if we have mark in cell type
              String szLine = brA_nmark[ncell].readLine();
              if (szLine == null) break;

              float fval = ((int) (nroundval * Float.parseFloat(szLine) + 0.5)) / froundval;
              databinfirst_nbeginread_nmark[ncell] = fval;
              databinlast_nmark[ncell][nbeginread] = fval;
            }
          }
        }
      }

      // System.out.println("B "+(System.currentTimeMillis()-ltime));
      // ltime = System.currentTimeMillis();
      int nbin = 0;

      for (int nbeginread = n2nmaxoffset;
          nbeginread < nchromsize;
          nbeginread += (ChromImpute.NUMLINES - n2nmaxoffset)) {

        boolean bok = true;
        // ltime = System.currentTimeMillis();
        // reads in the additional data elements needed
        for (nbin = n2nmaxoffset; nbin < ChromImpute.NUMLINES && bok; nbin++) {
          float[][] databinfirst_nbin = databinfirst[nbin];
          for (int nmark = 0; nmark < marksA.length; nmark++) {
            float[] databinfirst_nbin_nmark = databinfirst_nbin[nmark];
            float[][] databinlast_nmark = databinlast[nmark];
            boolean[] bmarkcell_nmark = bmarkcell[nmark];
            BufferedReader[] brA_nmark = brA[nmark];

            String szLine;

            for (int ncell = 0; ncell < cellsA.length; ncell++) {
              if (bmarkcell_nmark[ncell]) {
                if ((szLine = brA_nmark[ncell].readLine()) != null) {
                  // stores rounded value
                  float fval = ((int) (nroundval * Float.parseFloat(szLine) + 0.5)) / froundval;
                  databinfirst_nbin_nmark[ncell] = fval;
                  databinlast_nmark[ncell][nbin] = fval;
                } else {
                  bok = false;
                  // nbin--;
                }
              }
            }
          }
        }

        // updated nbin handling
        if ((!bdnamethyl) && (!bok)) {
          nbin--;
        }
        // System.out.println("C "+(System.currentTimeMillis()-ltime)+"\t"+nbin);
        // ltime = System.currentTimeMillis();

        nclassifierindex = 0;
        for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
          // if (//((nholdoutcellrequest == -1)||(nholdoutcellrequest==nholdoutcell))&&
          // (bdnamethyl&&bdnamethylcell[nholdoutcell]) ||
          //  (!bdnamethyl&&bmarkcell[ntargetmark][nholdoutcell]))
          // {
          // System.out.println("Holding out\t"+nholdoutcell+"\t"+nbeginread);
          // these calls will do all bags for this classifier index
          if ((bdnamethyl) && bdnamethylcell[nholdoutcell]) {
            // nbin has the number of lines read for currnet holding
            // nbeginread-n2nmaxoffset represents the overall starting position
            generateInstanceDataTrainDNAMethyl(
                databinlast,
                databinlastCELLMARK,
                databinfirst,
                nbin,
                nbeginread - n2nmaxoffset,
                // bmarkcell,bcellmark, markstargetcellA,ntargetcell,
                nholdoutcell,
                nclassifierindex,
                nchrom); // true,null,nchrom);
          } else if (!bdnamethyl && bmarkcell[ntargetmark][nholdoutcell]) {
            // nbin has the number of lines read from current holding
            // nbeginread-n2nmaxoffset represents the overall starting position
            generateInstanceDataTrain(
                databinlast,
                databinlastCELLMARK,
                databinfirst,
                nbin,
                nbeginread - n2nmaxoffset,
                // bmarkcell, bcellmark, markstargetcellA, ntargetcell,
                ntargetmark,
                nholdoutcell,
                nclassifierindex); // ,true,null);
          }
          // }
          nclassifierindex += numbags;
        }

        for (int ni = 0; ni < n2nmaxoffset; ni++) {
          // copys the contents of the last n2nmaxoffset to the front of the array
          float[][] databinfirst_ni = databinfirst[ni];
          float[][] databinfirst_nend = databinfirst[databinfirst.length - n2nmaxoffset + ni];
          for (int nmark = 0; nmark < databinfirst_ni.length; nmark++) {
            float[][] databinlast_nmark = databinlast[nmark];
            float[] databinfirst_ni_nmark = databinfirst_ni[nmark];
            float[] databinfirst_nend_nmark = databinfirst_nend[nmark];
            for (int ncell = 0; ncell < databinfirst_ni_nmark.length; ncell++) {
              float fval = databinfirst_nend_nmark[ncell];
              databinfirst_ni_nmark[ncell] = fval;
              databinlast_nmark[ncell][ni] = fval;
            }
          }
        }

        // System.out.println("E "+(System.currentTimeMillis()-ltime));
      }

      if (!bdnamethyl) {
        for (int na = 0; na < noverallpos.length; na++) {
          // for backwards compatibility treating chromosome as one extra position longer
          noverallpos[na]++;
        }
      }

      // generates data for the current chromosome
      // updated in v0.9.7 to close after each chromosome
      // closes out the bufferedreaders
      for (int na = 0; na < brA.length; na++) {
        for (int nb = 0; nb < brA[na].length; nb++) {
          if (brA[na][nb] != null) {
            brA[na][nb].close();
          }
        }
      }
    }

    // System.out.println("IN FINISH");
    // closes out the print writers
    for (int ni = 0; ni < trainPW.length; ni++) {
      if (trainPW[ni] != null) {
        trainPW[ni].finish();
        trainPW[ni].close();
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////

  /** Stores into attributes[nclassiferindex] the contents of szattributefile */
  public void loadFeatures(int nclassifierindex, String szattributefile) throws IOException {
    BufferedReader br = Util.getBufferedReader(szattributefile);
    attributes[nclassifierindex] = new ArrayList();
    String szLine;
    while ((szLine = br.readLine()) != null) {
      attributes[nclassifierindex].add(szLine);
    }
    br.close();
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /** Procedure generates the features that will be printed to the training files */
  public void generateFeatures( // boolean[][] bmarkcell,
      // int[] markstargetcellA,
      int ntargetmark, int nholdoutcell, int nclassifierindex, boolean bprintattributes)
      throws IOException {
    // System.out.println("Generating Features nclassifierindex = "+nclassifierindex);
    attributes[nclassifierindex] = new ArrayList(); // could get away with just using feature count

    int nfeaturecount = 0;

    for (int nmarktargetcell = 0; nmarktargetcell < marksA.length; nmarktargetcell++) {
      if ((nmarktargetcell != ntargetmark)
          && (bmarkcell[ntargetmark][nholdoutcell])
          && (bmarkcell[nmarktargetcell][nholdoutcell])) {
        // both the order by mark and the target mark are available in the replacement hold out cell
        // which is a different mark from the target
        int nhit = 1;
        for (int ncell = 0; ncell < bmarkcell[ntargetmark].length; ncell++) {
          // going to figure out the maximum number of knn features that can be made
          if ((ncell != nholdoutcell)) {
            // a nearest cell can't be the same as the replacement held out cell
            if ((bmarkcell[ntargetmark][ncell]) && (bmarkcell[nmarktargetcell][ncell])) {
              // both the target mark and mark being ordered by need to be available to use this
              // cell type
              if (nhit <= nmaxknn + 1) {
                // allowing an additional cell type above maxknn since might need to skip over one
                // cell type when that is the feature being predicted

                // creates a feature for each mark and ranked position
                // excluding cell type held out and cell type predicting here

                // adds in a feature based on the narrow distance
                if (bprintattributes) {
                  String szfeature =
                      nfeaturecount
                          + "\t"
                          + "OrderNarrow_"
                          + marksA[ntargetmark]
                          + "_"
                          + nhit
                          + "_by_"
                          + marksA[nmarktargetcell]
                          + "\n";
                  byte[] btformat = szfeature.getBytes();
                  pwattributes.write(btformat, 0, btformat.length);
                  nfeaturecount++;
                  szfeature =
                      nfeaturecount
                          + "\t"
                          + "OrderGlobal_"
                          + marksA[ntargetmark]
                          + "_"
                          + nhit
                          + "_by_"
                          + marksA[nmarktargetcell]
                          + "\n";
                  btformat = szfeature.getBytes();
                  pwattributes.write(btformat, 0, btformat.length);
                  nfeaturecount++;
                }

                int[] samples_nclassifierindex = samples[nclassifierindex];

                // System.out.println("Classifier "+nclassifierindex+"\t"+nfeaturecount+"\t"+
                //
                // "OrderNarrow_"+marksA[ntargetmark]+"_"+nhit+"_by_"+marksA[nmarktargetcell]);

                attributes[nclassifierindex].add(
                    "OrderNarrow_"
                        + marksA[ntargetmark]
                        + "_"
                        + nhit
                        + "_by_"
                        + marksA[nmarktargetcell]);

                // adds in a feature based on the global distance
                // System.out.println("Classifier "+nclassifierindex+"\t"+nfeaturecount+"\t"+
                //
                // "OrderGlobal_"+marksA[ntargetmark]+"_"+nhit+"_by_"+marksA[nmarktargetcell]);

                attributes[nclassifierindex].add(
                    "OrderGlobal_"
                        + marksA[ntargetmark]
                        + "_"
                        + nhit
                        + "_by_"
                        + marksA[nmarktargetcell]);

                nhit++;
              }
            }
          }
        }
      }
    }

    // System.out.println("Feature count order by\t"+nfeaturecount+"\t"+marksA.length);

    for (int nmark = 0; nmark < marksA.length; nmark++) {
      // int nmark = markstargetcellA[nmarkindex];
      // System.out.println("IMPUTEHERE "+nmark+" "+ ntargetmark+" "+nholdoutcell+"
      // "+bmarkcell[nmark][nholdoutcell]);
      if ((bmarkcell[nmark][nholdoutcell]) && (nmark != ntargetmark)) {
        // creating feature values for the replacement cell type
        // mark can't be the target mark

        if (bprintattributes) {
          String szfeature = nfeaturecount + "\t" + marksA[nmark] + "_center\n";
          byte[] btformat = szfeature.getBytes();
          pwattributes.write(btformat, 0, btformat.length);
          nfeaturecount++;
        }

        // creates a feature for the exact position

        attributes[nclassifierindex].add(marksA[nmark] + "_center");

        // creates features for the positions at the more spacing up to nmaxoffsetnarrow
        // puts them in position adding left and right
        for (int ni = nincrementnarrow; ni <= nmaxoffsetnarrow; ni += nincrementnarrow) {
          // DNase signal to the left, right, and cumulative

          // System.out.println("Classifier "+nclassifierindex+"\t"+marksA[nmark]+" left
          // "+(ni*nresolution));
          // System.out.println("Classifier "+nclassifierindex+"\t"+marksA[nmark]+" right
          // "+(ni*nresolution));

          if (bprintattributes) {
            String szfeature =
                nfeaturecount + "\t" + marksA[nmark] + "_left_" + (ni * nresolution) + "\n";
            byte[] btformat = szfeature.getBytes();
            pwattributes.write(btformat, 0, btformat.length);
            nfeaturecount++;
            szfeature =
                nfeaturecount + "\t" + marksA[nmark] + "_right_" + (ni * nresolution) + "\n";
            btformat = szfeature.getBytes();
            pwattributes.write(btformat, 0, btformat.length);
            nfeaturecount++;
          }

          attributes[nclassifierindex].add(marksA[nmark] + "_left_" + (ni * nresolution));
          attributes[nclassifierindex].add(marksA[nmark] + "_right_" + (ni * nresolution));
        }

        for (int ni = nmaxoffsetnarrow + nincrementwide;
            ni <= nmaxoffsetwide;
            ni += nincrementwide) {
          // System.out.println("Classifier "+nclassifierindex+"\t"+marksA[nmark]+" left
          // "+(ni*nresolution));
          // System.out.println("Classifier "+nclassifierindex+"\t"+marksA[nmark]+" right
          // "+(ni*nresolution));

          if (bprintattributes) {
            String szfeature =
                nfeaturecount + "\t" + marksA[nmark] + "_left_" + (ni * nresolution) + "\n";
            byte[] btformat = szfeature.getBytes();
            pwattributes.write(btformat, 0, btformat.length);
            nfeaturecount++;

            szfeature =
                nfeaturecount + "\t" + marksA[nmark] + "_right_" + (ni * nresolution) + "\n";
            btformat = szfeature.getBytes();
            pwattributes.write(btformat, 0, btformat.length);
            nfeaturecount++;
          }

          attributes[nclassifierindex].add(marksA[nmark] + "_left_" + (ni * nresolution));
          attributes[nclassifierindex].add(marksA[nmark] + "_right_" + (ni * nresolution));
        }
      }
    }

    // System.out.println("Classifier "+nclassifierindex+"\t"+"target "+marksA[ntargetmark]);
    // System.out.println("DEFOUT\t"+ntargetmark+"\t"+nholdoutcell+"\t"+nfeaturecount+"\t"+attributes[nclassifierindex].size());
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  public void generateFeaturesDNAMethyl( // boolean[][] bmarkcell,
      // int[] markstargetcellA,
      int nholdoutcell, int nclassifierindex, boolean bprintattributes) throws IOException {
    // System.out.println("Generating Features nclassifierindex = "+nclassifierindex);
    attributes[nclassifierindex] = new ArrayList();

    int nfeaturecount = 0;

    for (int nmarktargetcell = 0; nmarktargetcell < marksA.length; nmarktargetcell++) {
      if ((bdnamethylcell[nholdoutcell]) && (bmarkcell[nmarktargetcell][nholdoutcell])) {
        int nhit = 1;
        for (int ncell = 0; ncell < bdnamethylcell.length; ncell++) {
          if ((ncell != nholdoutcell)) {
            if ((bdnamethylcell[ncell]) && (bmarkcell[nmarktargetcell][ncell])) {

              if (nhit <= nmaxknn + 1) {
                // creates a feature for each mark and ranked position
                // excluding cell type held out and cell type predicting here
                // Attribute feature = new
                // Attribute("Order_"+marksA[nmark]+"_"+nhit+"_by_"+marksA[markstargetcellA[nmarktargetcell]]);

                if (bprintattributes) {
                  String szfeature =
                      nfeaturecount
                          + "\t"
                          + "OrderNarrow_DNAMethyl_"
                          + nhit
                          + "_by_"
                          + marksA[nmarktargetcell]
                          + "\n";
                  byte[] btformat = szfeature.getBytes();
                  pwattributes.write(btformat, 0, btformat.length);
                  nfeaturecount++;

                  szfeature =
                      nfeaturecount
                          + "\t"
                          + "OrderGlobal_DNAMethyl_"
                          + nhit
                          + "_by_"
                          + marksA[nmarktargetcell]
                          + "\n";
                  btformat = szfeature.getBytes();
                  pwattributes.write(btformat, 0, btformat.length);
                  nfeaturecount++;
                }

                // System.out.println("Classifier "+nclassifierindex+"\t"+nfeaturecount+"\t"+
                //                 "OrderNarrow_DNAMethyl_"+nhit+"_by_"+marksA[nmarktargetcell]);

                attributes[nclassifierindex].add(
                    "OrderNarrow_DNAMethyl_" + nhit + "_by_" + marksA[nmarktargetcell]);

                // System.out.println("Classifier "+nclassifierindex+"\t"+nfeaturecount+"\t"+
                //                 "OrderGlobal_DNAMethyl_"+nhit+"_by_"+marksA[nmarktargetcell]);

                attributes[nclassifierindex].add(
                    "OrderGlobal_DNAMethyl_" + nhit + "_by_" + marksA[nmarktargetcell]);

                nhit++;
              }
            }
          }
        }
      }
    }

    // System.out.println("Feature count order by\t"+nfeaturecount);

    for (int nmark = 0; nmark < marksA.length; nmark++) {
      // int nmark = markstargetcellA[nmarkindex];
      // System.out.println("IMPUTEHERE "+nmark+" "+ ntargetmark+" "+nholdoutcell+"
      // "+bmarkcell[nmark][nholdoutcell]);
      if (bmarkcell[nmark][nholdoutcell]) {
        if (bprintattributes) {
          String szfeature = nfeaturecount + "\t" + marksA[nmark] + "_center\n";
          byte[] btformat = szfeature.getBytes();
          pwattributes.write(btformat, 0, btformat.length);
          nfeaturecount++;
        }

        attributes[nclassifierindex].add(marksA[nmark] + "_center");

        for (int ni = nincrementnarrow; ni <= nmaxoffsetnarrow; ni += nincrementnarrow) // ++)
        {

          // System.out.println("Classifier "+nclassifierindex+"\t"+marksA[nmark]+" left
          // "+(ni*nresolution));
          // System.out.println("Classifier "+nclassifierindex+"\t"+marksA[nmark]+" right
          // "+(ni*nresolution));

          if (bprintattributes) {
            String szfeature =
                nfeaturecount + "\t" + marksA[nmark] + "_left_" + (ni * nresolution) + "\n";
            byte[] btformat = szfeature.getBytes();
            pwattributes.write(btformat, 0, btformat.length);
            nfeaturecount++;
            szfeature =
                nfeaturecount + "\t" + marksA[nmark] + "_right_" + (ni * nresolution) + "\n";
            btformat = szfeature.getBytes();
            pwattributes.write(btformat, 0, btformat.length);
            nfeaturecount++;
          }

          attributes[nclassifierindex].add(marksA[nmark] + "_left_" + (ni * nresolution));
          attributes[nclassifierindex].add(marksA[nmark] + "_right_" + (ni * nresolution));
        }

        for (int ni = nmaxoffsetnarrow + nincrementwide;
            ni <= nmaxoffsetwide;
            ni += nincrementwide) // ++)
        {
          // System.out.println("Classifier "+nclassifierindex+"\t"+marksA[nmark]+" left
          // "+(ni*nresolution));
          // System.out.println("Classifier "+nclassifierindex+"\t"+marksA[nmark]+" right
          // "+(ni*nresolution));

          if (bprintattributes) {
            String szfeature =
                nfeaturecount + "\t" + marksA[nmark] + "_left_" + (ni * nresolution) + "\n";
            byte[] btformat = szfeature.getBytes();
            pwattributes.write(btformat, 0, btformat.length);

            nfeaturecount++;

            szfeature =
                nfeaturecount + "\t" + marksA[nmark] + "_right_" + (ni * nresolution) + "\n";
            btformat = szfeature.getBytes();
            pwattributes.write(btformat, 0, btformat.length);

            nfeaturecount++;
          }

          attributes[nclassifierindex].add(marksA[nmark] + "_left_" + (ni * nresolution));
          attributes[nclassifierindex].add(marksA[nmark] + "_right_" + (ni * nresolution));
        }
      }
    }

    // System.out.println("Classifier "+nclassifierindex);
    // System.out.println("DEFOUT\t"+nholdoutcell+"\t"+nfeaturecount+"\t"+attributes[nclassifierindex].size());
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  static class RecMax {
    float dval;
    int nline;

    public RecMax(float dval, int nline) {
      this.dval = dval;
      this.nline = nline;
    }
  }

  static class RecMaxCompare implements Comparator, Serializable {

    public int compare(Object o1, Object o2) {
      RecMax r1 = (RecMax) o1;
      RecMax r2 = (RecMax) o2;

      if (r1.dval > r2.dval) {
        return -1;
      } else if (r1.dval < r2.dval) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  static class ROCRec {
    int ntotal;
    int nhit;
  }

  /** For comparing observed to imputed data */
  public void evaluateROC() throws IOException {

    double dkeepreal = 0;
    double dkeepimpute = 0;

    NumberFormat nf4 = NumberFormat.getInstance(Locale.ENGLISH);
    nf4.setGroupingUsed(false);
    nf4.setMaximumFractionDigits(4);

    HashMap hmtally = new HashMap();
    int nlinetotal = 0;

    // figuring out the needed cut-off for storing devalpercent2 of real data
    for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
      BufferedReader brreal;
      brreal =
          Util.getBufferedReader(
              szevalobserveddir + "/" + alchrom.get(nchrom) + "_" + szevalobservedfile);

      String szLine;
      brreal.readLine();
      brreal.readLine();

      while ((szLine = brreal.readLine()) != null) {
        double dval = Double.parseDouble(szLine);

        // Integer intcount = (Integer) hmtally.get(new Double(dval));
        Integer intcount = (Integer) hmtally.get(Double.valueOf(dval));

        if (intcount == null) {
          hmtally.put(Double.valueOf(dval), Integer.valueOf(1));
          // hmtally.put(new Double(dval), new Integer(1));
        } else {
          hmtally.put(Double.valueOf(dval), Integer.valueOf(1 + intcount.intValue()));
          // hmtally.put(new Double(dval), new Integer(1+intcount.intValue()));
        }
        nlinetotal++;
      }
      brreal.close();
    }

    Iterator itrkeys = hmtally.keySet().iterator();

    String szkey;

    int nindex2 = 0;
    double[] dvals = new double[hmtally.size()];
    while (itrkeys.hasNext()) {
      dvals[nindex2] = ((Double) itrkeys.next()).doubleValue();
      nindex2++;
    }
    Arrays.sort(dvals);

    int nhits = 0;
    int nmaxindex = dvals.length - 1;
    while ((nmaxindex >= 0) && (nhits / (double) nlinetotal < devalpercent2 / (double) 100)) {
      // nhits += ((Integer) hmtally.get(new Double(dvals[nmaxindex]))).intValue();
      nhits += ((Integer) hmtally.get(Double.valueOf(dvals[nmaxindex]))).intValue();
      dkeepreal = dvals[nmaxindex];
      nmaxindex--;
    }

    // hmtally = new HashMap();
    nlinetotal = 0;

    HashSet hsMax = new HashSet();
    int nlinereal = 0;
    ArrayList al = new ArrayList();

    // storing all observed data above the observed data keep threshold
    for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
      BufferedReader brreal;
      brreal =
          Util.getBufferedReader(
              szevalobserveddir + "/" + alchrom.get(nchrom) + "_" + szevalobservedfile);
      String szLine;
      brreal.readLine();
      brreal.readLine();

      while ((szLine = brreal.readLine()) != null) {
        double dval = Double.parseDouble(szLine);
        if (dval >= dkeepreal) {
          al.add(new RecMax((float) dval, nlinereal));
        }
        nlinereal++;
      }
      brreal.close();
    }

    int numtop1 = (int) (nlinereal * (devalpercent1 / (double) 100));

    RecMax[] alMax = new RecMax[al.size()];
    for (int nk = 0; nk < alMax.length; nk++) {
      alMax[nk] = (RecMax) al.get(nk);
    }
    Arrays.sort(alMax, new RecMaxCompare());

    for (int ni = 0; ni < numtop1; ni++) {
      hsMax.add(Integer.valueOf(alMax[ni].nline));
      // hsMax.add(new Integer(alMax[ni].nline));
    }

    HashMap hmRec = new HashMap();

    nlinereal = 0;

    if (bprintonefile) {
      BufferedReader brimpute = Util.getBufferedReader(szevalimputedir + "/" + szevalimputefile);
      // brimpute.readLine();
      if (bprintbrowserheader) {
        brimpute.readLine();
      }

      // BufferedReader brreal = null;

      String szLine;

      // nlines = 0;

      while ((szLine = brimpute.readLine()) != null) {
        if ((szLine.startsWith("fixed"))) // ||(szLine.startsWith("variable")))
        {

          // if (brreal != null)
          // {
          //   brreal.close();
          // }

          // StringTokenizer stchromline = new StringTokenizer(szLine, " \t=");
          // stchromline.nextToken();//fixed/variable
          // stchromline.nextToken();//chrom
          // String szcurrchrom = stchromline.nextToken();

          // brreal =
          // Util.getBufferedReader(szevalobserveddir+"/"+szcurrchrom+"_"+szevalobservedfile);

          // brreal.readLine();
          // brreal.readLine();
          szLine = brimpute.readLine();
        }

        // double dval1 = Double.parseDouble(brreal.readLine());
        double dval2 = Double.parseDouble(szLine);

        ROCRec theROCRec = (ROCRec) hmRec.get("" + dval2);

        if (theROCRec == null) {
          theROCRec = new ROCRec();
          hmRec.put("" + dval2, theROCRec);
        }

        theROCRec.ntotal++;
        // if (hsMax.contains(new Integer(nlinereal)))
        if (hsMax.contains(Integer.valueOf(nlinereal))) {
          theROCRec.nhit++;
        }

        // nlines++;
        nlinereal++;
      }
      brimpute.close();
    } else {
      for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
        // BufferedReader brreal;

        // brreal =
        // Util.getBufferedReader(szevalobserveddir+"/"+alchrom.get(nchrom)+"_"+szevalobservedfile);

        // brreal.readLine();
        // brreal.readLine();

        // int nindex = 0;
        String szLine;

        BufferedReader brimpute =
            Util.getBufferedReader(
                szevalimputedir + "/" + alchrom.get(nchrom) + "_" + szevalimputefile);
        brimpute.readLine();
        if (bprintbrowserheader) {
          brimpute.readLine();
        }
        // nlines = 0;

        while ((szLine = brimpute.readLine()) != null) {
          double dval2 = Double.parseDouble(szLine);

          ROCRec theROCRec = (ROCRec) hmRec.get("" + dval2);

          if (theROCRec == null) {
            theROCRec = new ROCRec();
            hmRec.put("" + dval2, theROCRec);
          }

          theROCRec.ntotal++;
          // if (hsMax.contains(new Integer(nlinereal)))
          if (hsMax.contains(Integer.valueOf(nlinereal))) {
            theROCRec.nhit++;
            // (dval2 >= dkeepimpute)
            // stores imputed
            // al.add(new RecMax((float) dval2,nlinereal));
          }

          // nlines++;
          nlinereal++;
        }
        brimpute.close();
        // brreal.close();
      }
    }

    Iterator itr = hmRec.keySet().iterator();
    double[] keys = new double[hmRec.size()];
    int nkey = 0;
    while (itr.hasNext()) {
      keys[nkey] = Double.parseDouble((String) itr.next());
      nkey++;
    }
    Arrays.sort(keys);

    int ntp = 0;
    int nfp = 0;
    int ntotalfalse = nlinereal - numtop1;
    for (int nindex = keys.length - 1; nindex >= 0; nindex--) {
      ROCRec theROCRec = (ROCRec) hmRec.get("" + keys[nindex]);

      ntp += theROCRec.nhit;
      nfp += theROCRec.ntotal - theROCRec.nhit;
      double dfprate = nfp / (double) ntotalfalse;
      double dtprate = ntp / (double) numtop1;

      System.out.println(
          keys[nindex]
              + "\t"
              + nf4.format(dfprate)
              + "\t"
              + nf4.format(dtprate)
              + "\t"
              + ((ntp + nfp) / (double) (numtop1 + ntotalfalse)));
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** For comparing observed to imputed data */
  public void evaluatePeaks() throws IOException {

    String szLine;

    NumberFormat nf4 = NumberFormat.getInstance(Locale.ENGLISH);
    nf4.setGroupingUsed(false);
    nf4.setMaximumFractionDigits(4);

    HashMap hmROCRec = new HashMap();
    int ntotaltrue = 0;
    int nlinetotalall = 0;

    for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
      boolean[] bpeakA = new boolean[(int) nmaxchromsize / nresolution + 1];
      BufferedReader brpeaks = Util.getBufferedReader(szpeakevalfile);
      String szcurrchrom = (String) alchrom.get(nchrom);
      while ((szLine = brpeaks.readLine()) != null) {
        StringTokenizer st = new StringTokenizer(szLine, "\t");
        String szchrom = st.nextToken();
        if (szchrom.equals(szcurrchrom)) {
          int nbegin = Integer.parseInt(st.nextToken()) / nresolution;
          int nend = (Integer.parseInt(st.nextToken()) - 1) / nresolution;

          for (int nk = nbegin; nk <= nend; nk++) {
            bpeakA[nk] = true;
            ntotaltrue++;
          }
        }
      }
      brpeaks.close();

      if (bprintonefile) {

        // "fixedStep chrom="+szchrom+" start=1 step="+nresolution+" span="+nresolution+"\n";
        BufferedReader brdata = Util.getBufferedReader(szevalimputedir + "/" + szevalimputefile);
        // String szLine;
        brdata.readLine();

        if (bprintbrowserheader) {
          brdata.readLine();
        }

        boolean bfound = false;
        // String szprefix = "fixedStep chrom="+alchrom.get(nchrom);
        while ((szLine = brdata.readLine()) != null) {
          if (szLine.startsWith("fixed")) {
            StringTokenizer stchromline = new StringTokenizer(szLine, " \t=");
            stchromline.nextToken(); // fixed/variable
            stchromline.nextToken(); // chrom
            String szchrom = stchromline.nextToken();
            if (szcurrchrom.equals(szchrom)) {
              bfound = true;
            }
          }
        }

        if (bfound) {
          int nline = 0;
          while (((szLine = brdata.readLine()) != null) && (!szLine.startsWith("fixed"))) {
            double dval1 = Double.parseDouble(szLine);
            ROCRec theROCRec = (ROCRec) hmROCRec.get("" + dval1);
            if (theROCRec == null) {
              theROCRec = new ROCRec();
              hmROCRec.put("" + dval1, theROCRec);
            }

            theROCRec.ntotal++;
            if (bpeakA[nline]) {
              theROCRec.nhit++;
            }
            nline++;
          }
          nlinetotalall += nline;
          brdata.close();
        }
      } else {
        BufferedReader brdata =
            Util.getBufferedReader(
                szevalimputedir + "/" + alchrom.get(nchrom) + "_" + szevalimputefile);
        // String szLine;
        brdata.readLine();

        if (bprintbrowserheader) {
          brdata.readLine();
        }

        int nline = 0;
        while ((szLine = brdata.readLine()) != null) {

          double dval1 = Double.parseDouble(szLine);
          ROCRec theROCRec = (ROCRec) hmROCRec.get("" + dval1);
          if (theROCRec == null) {
            theROCRec = new ROCRec();
            hmROCRec.put("" + dval1, theROCRec);
          }

          theROCRec.ntotal++;
          if (bpeakA[nline]) {
            theROCRec.nhit++;
          }
          nline++;
        }
        nlinetotalall += nline;
        brdata.close();
      }
    }

    Iterator itr = hmROCRec.keySet().iterator();
    double[] keys = new double[hmROCRec.size()];
    int nkey = 0;
    while (itr.hasNext()) {
      keys[nkey] = Double.parseDouble((String) itr.next());
      nkey++;
    }
    Arrays.sort(keys);

    int ntp = 0;
    int nfp = 0;
    int ntotalfalse = nlinetotalall - ntotaltrue;
    double dauctop = 0;
    double doldfprate = 0;
    double doldtprate = 0;
    for (int nindex = keys.length - 1; nindex >= 0; nindex--) {
      ROCRec theROCRec = (ROCRec) hmROCRec.get("" + keys[nindex]);

      ntp += theROCRec.nhit;
      nfp += theROCRec.ntotal - theROCRec.nhit;
      double dfprate = nfp / (double) ntotalfalse;
      double dtprate = ntp / (double) ntotaltrue;

      dauctop += (dfprate - doldfprate) * (dtprate + doldtprate) / 2.0;
      doldfprate = dfprate;
      doldtprate = dtprate;
      // System.out.println(keys[nindex]+"\t"+nf4.format(dfprate)+"\t"+nf4.format(dtprate)+"\t"+((ntp+nfp)/(double)(numtop1+ntotalfalse)));
    }

    String szheaderline = "AUC_Predict_Peaks_" + szpeakevalfile + "_with_IMPUTE";

    String szoutputline = nf4.format(dauctop);

    if (szevaloutfile == null) {
      System.out.println(szheaderline);
      System.out.println(szoutputline);
    } else {
      PrintWriter pweval = new PrintWriter(szevaloutfile);
      pweval.println(szheaderline);
      pweval.println(szoutputline);
      pweval.close();
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /** For comparing observed to imputed data */
  public void evaluate() throws IOException {

    double dkeepreal = 0;
    double dkeepimpute = 0;

    NumberFormat nf4 = NumberFormat.getInstance(Locale.ENGLISH);
    nf4.setGroupingUsed(false);
    nf4.setMaximumFractionDigits(4);

    HashMap hmtally = new HashMap();
    int nlinetotal = 0;

    // figuring out the needed cut-off for storing devalpercent2 of real data
    for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
      BufferedReader brreal;
      brreal =
          Util.getBufferedReader(
              szevalobserveddir + "/" + alchrom.get(nchrom) + "_" + szevalobservedfile);

      String szLine;
      brreal.readLine();
      brreal.readLine();

      while ((szLine = brreal.readLine()) != null) {
        double dval = Double.parseDouble(szLine);

        // Integer intcount = (Integer) hmtally.get(new Double(dval));
        Integer intcount = (Integer) hmtally.get(Double.valueOf(dval));

        if (intcount == null) {
          // hmtally.put(new Double(dval), new Integer(1));
          hmtally.put(Double.valueOf(dval), Integer.valueOf(1));
        } else {
          // hmtally.put(new Double(dval), new Integer(1+intcount.intValue()));
          hmtally.put(Double.valueOf(dval), Integer.valueOf(1 + intcount.intValue()));
        }
        nlinetotal++;
      }
      brreal.close();
    }

    Iterator itrkeys = hmtally.keySet().iterator();

    String szkey;

    int nindex2 = 0;
    double[] dvals = new double[hmtally.size()];
    while (itrkeys.hasNext()) {
      dvals[nindex2] = ((Double) itrkeys.next()).doubleValue();
      nindex2++;
    }
    Arrays.sort(dvals);

    int nhits = 0;
    int nmaxindex = dvals.length - 1;
    while ((nmaxindex >= 0) && (nhits / (double) nlinetotal < devalpercent2 / (double) 100)) {
      // nhits += ((Integer) hmtally.get(new Double(dvals[nmaxindex]))).intValue();
      nhits += ((Integer) hmtally.get(Double.valueOf(dvals[nmaxindex]))).intValue();
      dkeepreal = dvals[nmaxindex];
      nmaxindex--;
    }

    hmtally = new HashMap();
    nlinetotal = 0;

    if (bprintonefile) {
      BufferedReader brimpute;

      brimpute = Util.getBufferedReader(szevalimputedir + "/" + szevalimputefile);

      brimpute.readLine();
      if (bprintbrowserheader) {
        brimpute.readLine();
      }

      String szLine;
      while ((szLine = brimpute.readLine()) != null) {
        if (szLine.startsWith("fixed")) // ||(szLine.startsWith("variable")))
        {
          continue;
        }

        double dval = Double.parseDouble(szLine);

        // Integer intcount = (Integer) hmtally.get(new Double(dval));
        Integer intcount = (Integer) hmtally.get(Double.valueOf(dval));

        if (intcount == null) {
          // hmtally.put(new Double(dval), new Integer(1));
          hmtally.put(Double.valueOf(dval), Integer.valueOf(1));
        } else {
          // hmtally.put(new Double(dval), new Integer(1+intcount.intValue()));
          hmtally.put(Double.valueOf(dval), Integer.valueOf(1 + intcount.intValue()));
        }
        nlinetotal++;
      }
      brimpute.close();
    } else {
      // figuring out the needed cut-off for storing devalpercent1 of imputed data
      for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
        BufferedReader brimpute;

        brimpute =
            Util.getBufferedReader(
                szevalimputedir + "/" + alchrom.get(nchrom) + "_" + szevalimputefile);

        brimpute.readLine();
        if (bprintbrowserheader) {
          brimpute.readLine();
        }

        String szLine;
        while ((szLine = brimpute.readLine()) != null) {
          double dval = Double.parseDouble(szLine);

          // Integer intcount = (Integer) hmtally.get(new Double(dval));
          Integer intcount = (Integer) hmtally.get(Double.valueOf(dval));

          if (intcount == null) {
            hmtally.put(Double.valueOf(dval), Integer.valueOf(1));
            // hmtally.put(new Double(dval), new Integer(1));
          } else {
            // hmtally.put(new Double(dval), new Integer(1+intcount.intValue()));
            hmtally.put(Double.valueOf(dval), Integer.valueOf(1 + intcount.intValue()));
          }
          nlinetotal++;
        }
        brimpute.close();
      }
    }

    itrkeys = hmtally.keySet().iterator();
    nindex2 = 0;
    dvals = new double[hmtally.size()];
    while (itrkeys.hasNext()) {
      dvals[nindex2] = ((Double) itrkeys.next()).doubleValue();
      nindex2++;
    }
    Arrays.sort(dvals);

    nhits = 0;
    nmaxindex = dvals.length - 1;
    while ((nmaxindex >= 0) && (nhits / (double) nlinetotal < devalpercent2 / (double) 100)) {
      // nhits += ((Integer) hmtally.get(new Double(dvals[nmaxindex]))).intValue();
      nhits += ((Integer) hmtally.get(Double.valueOf(dvals[nmaxindex]))).intValue();
      dkeepimpute = dvals[nmaxindex];
      nmaxindex--;
    }

    // System.out.println("BOTH_"+devalpercent1+"\tIMPUTE_"+devalpercent1+"_OBSERVED_"+devalpercent2+"\tOBSERVED_"+devalpercent1+"_IMPUTE_"+devalpercent2+"\tCorrelation\t"+
    //
    // "OBSERVED_"+devalpercent1+"_AUC_PREDICT_IMPUTE\tIMPUTE_"+devalpercent1+"_AUC_PREDICT_OBSERVE\t");
    // int nlines = 0;

    HashSet hsMax = new HashSet();
    HashSet hsMax2 = new HashSet();

    int nlinereal = 0;
    ArrayList al = new ArrayList();

    // storing all observed data above the observed data keep threshold
    for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
      BufferedReader brreal;
      brreal =
          Util.getBufferedReader(
              szevalobserveddir + "/" + alchrom.get(nchrom) + "_" + szevalobservedfile);
      String szLine;
      brreal.readLine();
      brreal.readLine();

      while ((szLine = brreal.readLine()) != null) {
        double dval = Double.parseDouble(szLine);
        if (dval >= dkeepreal) {
          al.add(new RecMax((float) dval, nlinereal));
        }
        nlinereal++;
      }
      brreal.close();
    }

    int numtop1 = (int) (nlinereal * (devalpercent1 / (double) 100));
    int numtop2 = (int) (nlinereal * (devalpercent2 / (double) 100));

    RecMax[] alMax = new RecMax[al.size()];
    for (int nk = 0; nk < alMax.length; nk++) {
      alMax[nk] = (RecMax) al.get(nk);
    }
    Arrays.sort(alMax, new RecMaxCompare());

    for (int ni = 0; ni < numtop1; ni++) {
      hsMax.add(Integer.valueOf(alMax[ni].nline));
      // hsMax.add(new Integer(alMax[ni].nline));
    }

    for (int ni = 0; ni < numtop2; ni++) {
      // hsMax2.add(new Integer(alMax[ni].nline));
      hsMax2.add(Integer.valueOf(alMax[ni].nline));
    }

    HashSet hsMaxImpute = new HashSet();

    int nlineimpute = 0;
    ArrayList alimpute = new ArrayList();

    if (bprintonefile) {
      BufferedReader brimpute;

      brimpute = Util.getBufferedReader(szevalimputedir + "/" + szevalimputefile);

      brimpute.readLine();
      if (bprintbrowserheader) {
        brimpute.readLine();
      }

      String szLine;

      while ((szLine = brimpute.readLine()) != null) {
        if (szLine.startsWith("fixed")) // ||(szLine.startsWith("variable")))
        {
          continue;
        }

        double dval = Double.parseDouble(szLine);

        if (dval >= dkeepimpute) {
          alimpute.add(new RecMax((float) dval, nlineimpute));
        }

        nlineimpute++;
      }
      brimpute.close();

    } else {
      for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
        BufferedReader brimpute;
        brimpute =
            Util.getBufferedReader(
                szevalimputedir + "/" + alchrom.get(nchrom) + "_" + szevalimputefile);

        String szLine;
        brimpute.readLine();
        brimpute.readLine();

        while ((szLine = brimpute.readLine()) != null) {
          double dval = Double.parseDouble(szLine);
          if (dval >= dkeepimpute) {
            alimpute.add(new RecMax((float) dval, nlineimpute));
          }
          nlineimpute++;
        }
        brimpute.close();
      }
    }

    RecMax[] alMaxImpute = new RecMax[alimpute.size()];
    for (int nk = 0; nk < alMaxImpute.length; nk++) {
      alMaxImpute[nk] = (RecMax) alimpute.get(nk);
    }
    Arrays.sort(alMaxImpute, new RecMaxCompare());

    for (int ni = 0; ni < numtop1; ni++) {
      hsMaxImpute.add(Integer.valueOf(alMaxImpute[ni].nline));
      // hsMaxImpute.add(new Integer(alMaxImpute[ni].nline));
    }

    // al = new ArrayList();
    nlinereal = 0;

    double dsumx = 0;
    double dsumy = 0;
    double dsumxsq = 0;
    double dsumysq = 0;
    double dsumxy = 0;

    HashMap hmROCRecImpute = new HashMap();
    HashMap hmROCRecReal = new HashMap();

    if (bprintonefile) {
      BufferedReader brimpute = Util.getBufferedReader(szevalimputedir + "/" + szevalimputefile);
      // brimpute.readLine();
      if (bprintbrowserheader) {
        brimpute.readLine();
      }

      BufferedReader brreal = null;

      String szLine;

      // nlines = 0;

      while ((szLine = brimpute.readLine()) != null) {
        if ((szLine.startsWith("fixed"))) // ||(szLine.startsWith("variable")))
        {

          if (brreal != null) {
            brreal.close();
          }

          StringTokenizer stchromline = new StringTokenizer(szLine, " \t=");
          stchromline.nextToken(); // fixed/variable
          stchromline.nextToken(); // chrom
          String szcurrchrom = stchromline.nextToken();

          brreal =
              Util.getBufferedReader(
                  szevalobserveddir + "/" + szcurrchrom + "_" + szevalobservedfile);

          brreal.readLine();
          brreal.readLine();
          szLine = brimpute.readLine();
        }

        double dval1 = Double.parseDouble(brreal.readLine());
        double dval2 = Double.parseDouble(szLine);

        ROCRec theROCRecImpute = (ROCRec) hmROCRecImpute.get("" + dval2);
        if (theROCRecImpute == null) {
          theROCRecImpute = new ROCRec();
          hmROCRecImpute.put("" + dval2, theROCRecImpute);
        }

        theROCRecImpute.ntotal++;
        // if (hsMax.contains(new Integer(nlinereal)))
        if (hsMax.contains(Integer.valueOf(nlinereal))) {
          theROCRecImpute.nhit++;
        }

        ROCRec theROCRecReal = (ROCRec) hmROCRecReal.get("" + dval1);
        if (theROCRecReal == null) {
          theROCRecReal = new ROCRec();
          hmROCRecReal.put("" + dval1, theROCRecReal);
        }

        theROCRecReal.ntotal++;

        // if (hsMaxImpute.contains(new Integer(nlinereal)))
        if (hsMaxImpute.contains(Integer.valueOf(nlinereal))) {
          theROCRecReal.nhit++;
        }

        dsumx += dval1;
        dsumxsq += dval1 * dval1;
        dsumy += dval2;
        dsumysq += dval2 * dval2;
        dsumxy += dval1 * dval2;

        // if (dval2 >= dkeepimpute)
        // {
        // stores imputed
        // al.add(new RecMax((float) dval2,nlinereal));
        // }
        // nlines++;
        nlinereal++;
      }
      brimpute.close();
      brreal.close();
    } else {
      for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
        BufferedReader brreal;

        brreal =
            Util.getBufferedReader(
                szevalobserveddir + "/" + alchrom.get(nchrom) + "_" + szevalobservedfile);

        brreal.readLine();
        brreal.readLine();

        // int nindex = 0;
        String szLine;

        BufferedReader brimpute =
            Util.getBufferedReader(
                szevalimputedir + "/" + alchrom.get(nchrom) + "_" + szevalimputefile);
        brimpute.readLine();
        if (bprintbrowserheader) {
          brimpute.readLine();
        }
        // nlines = 0;

        while ((szLine = brimpute.readLine()) != null) {
          double dval1 = Double.parseDouble(brreal.readLine());
          double dval2 = Double.parseDouble(szLine);

          ROCRec theROCRecImpute = (ROCRec) hmROCRecImpute.get("" + dval2);
          if (theROCRecImpute == null) {
            theROCRecImpute = new ROCRec();
            hmROCRecImpute.put("" + dval2, theROCRecImpute);
          }

          theROCRecImpute.ntotal++;
          // if (hsMax.contains(new Integer(nlinereal)))
          if (hsMax.contains(Integer.valueOf(nlinereal))) {
            theROCRecImpute.nhit++;
          }

          ROCRec theROCRecReal = (ROCRec) hmROCRecReal.get("" + dval1);
          if (theROCRecReal == null) {
            theROCRecReal = new ROCRec();
            hmROCRecReal.put("" + dval1, theROCRecReal);
          }

          theROCRecReal.ntotal++;
          // if (hsMaxImpute.contains(new Integer(nlinereal)))
          if (hsMaxImpute.contains(Integer.valueOf(nlinereal))) {
            theROCRecReal.nhit++;
          }

          dsumx += dval1;
          dsumxsq += dval1 * dval1;
          dsumy += dval2;
          dsumysq += dval2 * dval2;
          dsumxy += dval1 * dval2;

          // if (dval2 >= dkeepimpute)
          // {
          // stores imputed
          //  al.add(new RecMax((float) dval2,nlinereal));
          // }

          // nlines++;
          nlinereal++;
        }
        brimpute.close();
        brreal.close();
      }
    }

    // alMax = new RecMax[al.size()];
    // for (int nk = 0; nk < alMax.length; nk++)
    // {
    //   alMax[nk] = (RecMax) al.get(nk);
    // }
    // Arrays.sort(alMax, new RecMaxCompare());
    // sorts the imputed data matches

    int ntopmatchboth = 0;
    int ntopmatchreal = 0;
    int ntopmatchimpute = 0;

    // counts matches with the observed data
    // checks if smaller imputed overlaps with smaller observed
    for (int ni = 0; ni < numtop1; ni++) {
      // if (hsMax.contains(new Integer(alMaxImpute[ni].nline)))
      if (hsMax.contains(Integer.valueOf(alMaxImpute[ni].nline))) {
        ntopmatchboth++;
      }
    }

    // checks if larger imputed overlaps with smaller observed
    for (int ni = 0; ni < numtop2; ni++) {
      // if (hsMax.contains(new Integer(alMaxImpute[ni].nline)))
      if (hsMax.contains(Integer.valueOf(alMaxImpute[ni].nline))) {
        ntopmatchreal++;
      }
    }

    // checks if smaller imputed overlaps greater observed
    for (int ni = 0; ni < numtop1; ni++) {
      // if (hsMax2.contains(new Integer(alMaxImpute[ni].nline)))
      if (hsMax2.contains(Integer.valueOf(alMaxImpute[ni].nline))) {
        ntopmatchimpute++;
      }
    }

    // computes the correlation coefficient
    double dcorr;
    double dvarx = dsumxsq - dsumx * dsumx / nlinereal;
    double dvary = dsumysq - dsumy * dsumy / nlinereal;
    double dvarxdvary = dvarx * dvary;
    if (dvarxdvary <= 0) {
      dcorr = 0;
    } else {
      dcorr = (dsumxy - dsumx * dsumy / nlinereal) / Math.sqrt(dvarxdvary);
    }

    Iterator itr = hmROCRecImpute.keySet().iterator();
    double[] keys = new double[hmROCRecImpute.size()];
    int nkey = 0;
    while (itr.hasNext()) {
      keys[nkey] = Double.parseDouble((String) itr.next());
      nkey++;
    }
    Arrays.sort(keys);

    int ntp = 0;
    int nfp = 0;
    int ntotalfalse = nlinereal - numtop1;
    double dauctopreal = 0;
    double doldfprate = 0;
    double doldtprate = 0;
    for (int nindex = keys.length - 1; nindex >= 0; nindex--) {
      ROCRec theROCRec = (ROCRec) hmROCRecImpute.get("" + keys[nindex]);

      ntp += theROCRec.nhit;
      nfp += theROCRec.ntotal - theROCRec.nhit;
      double dfprate = nfp / (double) ntotalfalse;
      double dtprate = ntp / (double) numtop1;

      dauctopreal += (dfprate - doldfprate) * (dtprate + doldtprate) / 2.0;
      doldfprate = dfprate;
      doldtprate = dtprate;
      // System.out.println(keys[nindex]+"\t"+nf4.format(dfprate)+"\t"+nf4.format(dtprate)+"\t"+((ntp+nfp)/(double)(numtop1+ntotalfalse)));
    }

    itr = hmROCRecReal.keySet().iterator();
    keys = new double[hmROCRecReal.size()];
    nkey = 0;
    while (itr.hasNext()) {
      keys[nkey] = Double.parseDouble((String) itr.next());
      nkey++;
    }
    Arrays.sort(keys);

    ntp = 0;
    nfp = 0;
    ntotalfalse = nlinereal - numtop1;
    double dauctopimpute = 0;
    doldfprate = 0;
    doldtprate = 0;
    for (int nindex = keys.length - 1; nindex >= 0; nindex--) {
      ROCRec theROCRec = (ROCRec) hmROCRecReal.get("" + keys[nindex]);

      ntp += theROCRec.nhit;
      nfp += theROCRec.ntotal - theROCRec.nhit;
      double dfprate = nfp / (double) ntotalfalse;
      double dtprate = ntp / (double) numtop1;

      dauctopimpute += (dfprate - doldfprate) * (dtprate + doldtprate) / 2.0;
      doldfprate = dfprate;
      doldtprate = dtprate;
      // System.out.println(keys[nindex]+"\t"+nf4.format(dfprate)+"\t"+nf4.format(dtprate)+"\t"+((ntp+nfp)/(double)(numtop1+ntotalfalse)));
    }

    String szheaderline =
        "BOTH_"
            + devalpercent1
            + "\tIMPUTE_"
            + devalpercent1
            + "_OBSERVED_"
            + devalpercent2
            + "\tOBSERVED_"
            + devalpercent1
            + "_IMPUTE_"
            + devalpercent2
            + "\tCorrelation\t"
            + "IMPUTE_"
            + devalpercent1
            + "_AUC_PREDICT_OBSERVE\tOBSERVED_"
            + devalpercent1
            + "_AUC_PREDICT_IMPUTE";

    String szoutputline =
        nf4.format(100 * ntopmatchboth / (double) numtop1)
            + "\t"
            + nf4.format(100 * ntopmatchimpute / (double) numtop1)
            + "\t"
            + nf4.format(100 * ntopmatchreal / (double) numtop1)
            + "\t"
            + nf4.format(dcorr)
            + "\t"
            + nf4.format(dauctopimpute)
            + "\t"
            + nf4.format(dauctopreal);

    if (szevaloutfile == null) {
      System.out.println(szheaderline);
      System.out.println(szoutputline);
    } else {
      PrintWriter pweval = new PrintWriter(szevaloutfile);
      pweval.println(szheaderline);
      pweval.println(szoutputline);
      pweval.close();
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** Reads in classifiers saved in files */
  public void loadClassifiers() throws IOException {

    numclassifiers = 0;

    if (BLINEAR) {
        theClassifierLinearA = new RegressionLinear[numcells * numbags];
    } else {
        theClassifierA = new RegressionTree[numcells * numbags];
    }
    System.out.println("Made classifier of size: " + numcells * numbags + " - " + numcells + " x " + numbags);

    int nclassifierindex = 0;
    for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
      if (nholdoutcell != ntargetcell) {
        for (int nbag = 0; nbag < numbags; nbag++) {
          String szFile;

          if (BLINEAR) {
            szFile =
                szclassifierdir
                    + "/"
                    + "linearclassifier_"
                    + szoutcell
                    + "_"
                    + szoutmark
                    + "_"
                    + nholdoutcell
                    + "_"
                    + nbag
                    + ".txt.gz";
          } else {
            szFile =
                szclassifierdir
                    + "/"
                    + "classifier_"
                    + szoutcell
                    + "_"
                    + szoutmark
                    + "_"
                    + nholdoutcell
                    + "_"
                    + nbag
                    + ".txt.gz";
          }

          File f = new File(szFile);

          if (f.exists()) {
            String szattributefile =
                szclassifierdir
                    + "/"
                    + "useattributes_"
                    + szoutcell
                    + "_"
                    + szoutmark
                    + "_"
                    + nholdoutcell
                    + "_"
                    + nbag
                    + ".txt.gz";
            loadFeatures(nclassifierindex, szattributefile);

            if (BLINEAR) {
              theClassifierLinearA[nclassifierindex] =
                  new RegressionLinear(Util.getBufferedReader(szFile));
            } else {
              theClassifierA[nclassifierindex] = new RegressionTree(Util.getBufferedReader(szFile));
            }
            // increment total classifier count of available classifiers
            numclassifiers++;
          }
          nclassifierindex++; // increments the classifier index
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  public void insertionsortFirstK(BaseDistRec[] theArray, int nk) {
    if (nk > theArray.length) nk = theArray.length;
    for (int nindex = 1; nindex < nk; nindex++) {
      BaseDistRec insertRec = theArray[nindex];
      double dval = insertRec.ddist;
      int ncell = insertRec.ncell;
      int nholepos = nindex;

      BaseDistRec theArray_nholepos = theArray[nholepos];

      while (nholepos > 0) {
        // while there is position to the left
        BaseDistRec theArray_nholeposm1 = theArray[nholepos - 1];
        if (dval < theArray_nholeposm1.ddist) {
          // element to the right is greater than element to the left we are swapping
          theArray_nholepos.ddist = theArray_nholeposm1.ddist;
          theArray_nholepos.ncell = theArray_nholeposm1.ncell;
          theArray_nholepos = theArray_nholeposm1;
          nholepos--;
        } else {
          break;
        }
      }

      if (nholepos < nindex) {
        // we've made a swap, storing back the predicted values
        theArray_nholepos.ddist = dval;
        theArray_nholepos.ncell = ncell;
      }
    }
    // first k now sorted

    // now seeing if any more past first k can crack that
    for (int nindex = nk; nindex < theArray.length; nindex++) {
      BaseDistRec insertRec = theArray[nindex];
      double dval = insertRec.ddist;
      int ncell = insertRec.ncell;
      // only cares about getting sorted order from position k correct
      int nholepos = nk;

      //
      BaseDistRec theArray_nholepos = theArray[nholepos];

      while (nholepos > 0) {
        BaseDistRec theArray_nholeposm1 = theArray[nholepos - 1];

        if (dval < theArray_nholeposm1.ddist) {
          theArray_nholepos.ncell = theArray_nholeposm1.ncell;
          theArray_nholepos.ddist = theArray_nholeposm1.ddist;
          nholepos--;
          theArray_nholepos = theArray_nholeposm1;
        } else {
          break;
        }
      }

      if (nholepos < nk) {
        // we've made a change storing the value to insert into the array
        theArray_nholepos.ncell = ncell;
        theArray_nholepos.ddist = dval;
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////

  public void insertionsortFirstK(BaseDistRec[] theArray, int nk, int[] globalRank) {
    if (nk > theArray.length) nk = theArray.length;

    for (int nindex = 1; nindex < nk; nindex++) {
      BaseDistRec insertRec = theArray[nindex];
      double dval = insertRec.ddist;
      int ncell = insertRec.ncell;
      int nholepos = nindex;

      BaseDistRec theArray_nholepos = theArray[nholepos];

      while (nholepos > 0) {
        // while there is position to the left
        BaseDistRec theArray_nholeposm1 = theArray[nholepos - 1];
        if ((dval < theArray_nholeposm1.ddist)
            || ((dval == theArray_nholeposm1.ddist)
                && (globalRank[ncell] < globalRank[theArray_nholeposm1.ncell]))) {
          // element to the right is greater than element to the left we are swapping
          theArray_nholepos.ddist = theArray_nholeposm1.ddist;
          theArray_nholepos.ncell = theArray_nholeposm1.ncell;
          theArray_nholepos = theArray_nholeposm1;
          nholepos--;
        } else {
          break;
        }
      }

      if (nholepos < nindex) {
        // we've made a swap, storing back the predicted values
        theArray_nholepos.ddist = dval;
        theArray_nholepos.ncell = ncell;
      }
    }
    // first k now sorted

    // now seeing if any more past first k can crack that
    for (int nindex = nk; nindex < theArray.length; nindex++) {
      BaseDistRec insertRec = theArray[nindex];
      double dval = insertRec.ddist;
      int ncell = insertRec.ncell;
      // only cares about getting sorted order from position k correct
      int nholepos = nk;

      //
      BaseDistRec theArray_nholepos = theArray[nholepos];

      while (nholepos > 0) {
        BaseDistRec theArray_nholeposm1 = theArray[nholepos - 1];

        if ((dval < theArray_nholeposm1.ddist)
            || ((dval == theArray_nholeposm1.ddist)
                && (globalRank[ncell] < globalRank[theArray_nholeposm1.ncell]))) {
          theArray_nholepos.ncell = theArray_nholeposm1.ncell;
          theArray_nholepos.ddist = theArray_nholeposm1.ddist;
          nholepos--;
          theArray_nholepos = theArray_nholeposm1;
        } else {
          break;
        }
      }

      if (nholepos < nk) {
        // we've made a change storing the value to insert into the array
        theArray_nholepos.ncell = ncell;
        theArray_nholepos.ddist = dval;
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** Stores into array list */
  public void generateInstanceDataTrain(
      float[][][] databinlast,
      float[][][] databinlastCELLMARK,
      float[][][] databinfirst,
      int numbin,
      int nbegin,
      // boolean[][] bmarkcell, boolean[][] bcellmark,int[] markstargetcellA,
      int ntargetmark,
      int nholdoutcell,
      int nclassifierindex)
      throws Exception {

    // features are target cell type

    // long lstartoverall = System.currentTimeMillis();
    // long ltimeprint = 0;
    // long ltimeyes = 0;

    // NumberFormat nf2 = NumberFormat.getInstance();
    // nf2.setMaximumFractionDigits(2);
    // nf2.setGroupingUsed(false);

    int[][] savedorder = new int[marksA.length][nmaxknn + 1];

    boolean[] bmarkcell_ntargetmark = bmarkcell[ntargetmark];
    // long ltimeinner = 0;
    // long ltimesort = 0;
    // long ltimemark = 0;

    // int nusecelltarget = nholdoutcell;

    // System.out.println("HOLDOUT is "+nusecelltarget+" "+nholdoutcell);

    // initializing the global ordering
    int[][] global = new int[distmarkcellAL.length][numcells];
    int[][] globalRank = new int[distmarkcellAL.length][numcells];
    for (int na = 0; na < global.length; na++) {
      for (int nb = 0; nb < global[na].length; nb++) {
        // no element at that position if neg
        global[na][nb] = -1;
        globalRank[na][nb] = -1;
      }
    }

    // for (int ntargetmarkindex = 0; ntargetmarkindex < markstargetcellA.length;
    // ntargetmarkindex++)
    for (int nmarktarget = 0; nmarktarget < marksA.length; nmarktarget++) {
      // going through all marks in the target cell
      // int nmarktarget = markstargetcellA[ntargetmarkindex];

      // gets actual index of target mark
      if (bmarkcell[nmarktarget][nholdoutcell]) {
        // probably don't need to do this if nmarktarget is same ntargetmark
        // we have this mark in the holdout cell

        int npos = 0;

        for (int nindex = 0; nindex < distmarkcellAL[nmarktarget][nholdoutcell].size(); nindex++) {
          int ncell = ((Integer) distmarkcellAL[nmarktarget][nholdoutcell].get(nindex)).intValue();

          if ((ncell != nholdoutcell)
              && (bmarkcell[nmarktarget][ncell])
              && (bmarkcell_ntargetmark[ncell])) {
            // only going to use the cell if both the target and mark target are there
            // and is not the holdout cell
            global[nmarktarget][npos] = ncell;
            globalRank[nmarktarget][ncell] = npos;

            npos++;
          }
        }
      }
    }

    // long ltimeinit = System.currentTimeMillis()-lstartoverall;

    for (int nbag = 0; nbag < numbags; nbag++) {
      // if ((nbagrequest == -1)||(nbagrequest == nbag))
      // {
      // no bag request or we have the bag request we want

      // ArrayList theCurrInstances;
      // ArrayList theCurrInstancesOutput = null;

      // theCurrInstances       = theTrainInstances[nclassifierindex];
      // theCurrInstancesOutput = theTrainInstancesOutput[nclassifierindex];

      // now moving to a new classifier index
      GZIPOutputStream trainPW_nclassifierindex = trainPW[nclassifierindex];
      int[] samples_nclassifierindex = samples[nclassifierindex];

      int numfeatures = attributes[nclassifierindex].size();
      // System.out.println("getting
      // numfeatures\t"+nclassifierindex+"\t"+attributes.length+"\t"+numfeatures);

      String[] theInstance = new String[numfeatures];

      float[] databinlast_ntargetmark_nholdoutcell = databinlast[ntargetmark][nholdoutcell];

      // used for locally sorting the nearest cells
      BaseDistRec[] theBaseDistRecA = new BaseDistRec[bmarkcell[0].length];

      for (int nindex = 0; nindex < theBaseDistRecA.length; nindex++) {
        theBaseDistRecA[nindex] = new BaseDistRec();
      }

      int nmaxPARAM = Math.max(nknnoffset, nmaxoffsetwide);
      for (int nbinindex = nmaxPARAM; nbinindex < numbin - nmaxPARAM; nbinindex++) {
        // considering all positions that don't go off the edge
        int nsampleindexA_nclassifierindex = nsampleindexA[nclassifierindex];
        // corresponds to the sample index of the classifier

        // moved databifirstinto if then
        // while instead of if since can have multiple samples

        while ((nsampleindexA_nclassifierindex < samples_nclassifierindex.length)
            && (samples_nclassifierindex[nsampleindexA_nclassifierindex]
                <= noverallpos[nclassifierindex]))
        // changed from == to <= for back compatibility
        // (samples_nclassifierindex[nsampleindexA_nclassifierindex] ==
        // noverallpos[nclassifierindex]))
        {
          // sample position matches overall position for classifier
          // if (nclassifierindex == 0)
          //  System.out.println("PQR\t"+nbinindex+"\t"+nsampleindexA_nclassifierindex+"\t"+
          // samples_nclassifierindex[nsampleindexA_nclassifierindex]+"\t"+noverallpos[nclassifierindex]+"\t"+numbin+"\t"+nmaxPARAM+"\t"+(numbin-nmaxPARAM));

          float[][] databinfirst_nbinindex = databinfirst[nbinindex];
          // this is a position we've sampled

          int ninstanceindex = 0;
          // if (bmarkcell[ntargetmark][nholdoutcell])
          // {
          // for (int ntargetmarkindex = 0; ntargetmarkindex < markstargetcellA.length;
          // ntargetmarkindex++)
          for (int nmarktarget = 0; nmarktarget < marksA.length; nmarktarget++) {
            // going through all marks in the target cell
            // int nhit = 0;
            // int nmarktarget = markstargetcellA[ntargetmarkindex];

            // gets actual index of target mark
            if ((bmarkcell[nmarktarget][nholdoutcell]) && (nmarktarget != ntargetmark)) {
              boolean[] bmarkcell_nmarktarget =
                  bmarkcell[nmarktarget]; // array stating if mark is present in cell
              float[][] databinlast_nmarktarget = databinlast[nmarktarget];
              float[] databinlast_nmarktarget_nholdoutcell = databinlast_nmarktarget[nholdoutcell];
              // float[] databinfirst_nbinindex_nmarktarget = databinfirst_nbinindex[nmarktarget];

              int[] global_nmarktarget = global[nmarktarget];
              int[] globalRank_nmarktarget = globalRank[nmarktarget];
              // array with data for target mark and cell
              // long ltimestart = System.currentTimeMillis();

              for (int ncell = 0; ncell < bmarkcell_nmarktarget.length; ncell++) {
                if ((ncell != nholdoutcell)
                    && (bmarkcell_nmarktarget[ncell])
                    && (bmarkcell_ntargetmark[ncell])) {
                  float ddist = 0;
                  float[] databinlast_nmarktarget_ncell = databinlast_nmarktarget[ncell];
                  int nchange = nbinindex - nknnoffset;
                  // computing euclidean distance for all points within nknnoffset positions
                  for (int npos = -nknnoffset; npos <= nknnoffset; npos++) {
                    float ddiff =
                        (databinlast_nmarktarget_ncell[nchange]
                            - databinlast_nmarktarget_nholdoutcell[nchange]);
                    nchange++;
                    ddist += ddiff * ddiff; // takes the squared error of the signal difference
                  }

                  BaseDistRec theBaseDistRecA_ncell = theBaseDistRecA[ncell];
                  theBaseDistRecA_ncell.ddist = ddist;
                  theBaseDistRecA_ncell.ncell = ncell;
                } else {
                  // cell type is not eligible gives maximum value
                  BaseDistRec theBaseDistRecA_ncell = theBaseDistRecA[ncell];
                  theBaseDistRecA_ncell.ddist = Double.MAX_VALUE;
                  theBaseDistRecA_ncell.ncell = ncell;
                }
              }

              // ltimeinner += System.currentTimeMillis()-ltimestart;
              // int[] savedorder_ntargetmarkindex = savedorder[ntargetmarkindex];
              // getting out the saved order of top ranking cell types for this mark
              int[] savedorder_nmarktarget = savedorder[nmarktarget];

              // long ltimesortstart = System.currentTimeMillis();
              // initializes the sorted order to approximately what was used at the last position
              // not exactly since doing in place swaps
              for (int nrankpos = 0; nrankpos <= nmaxknn; nrankpos++) {
                int nbestindex = savedorder_nmarktarget[nrankpos];
                BaseDistRec temprec = theBaseDistRecA[nbestindex];
                theBaseDistRecA[nbestindex] = theBaseDistRecA[nrankpos];
                theBaseDistRecA[nrankpos] = temprec;
              }

              if (btieglobal) {
                insertionsortFirstK(theBaseDistRecA, nmaxknn + 1, globalRank_nmarktarget);
              } else {
                insertionsortFirstK(theBaseDistRecA, nmaxknn + 1);
              }

              for (int nrankpos = 0; nrankpos <= nmaxknn; nrankpos++) {
                // savedorder_ntargetmarkindex[nrankpos] = theBaseDistRecA[nrankpos].ncell;
                savedorder_nmarktarget[nrankpos] = theBaseDistRecA[nrankpos].ncell;
              }

              float[] databinfirst_nbinindex_ntargetmark = databinfirst_nbinindex[ntargetmark];
              int numknn = 0;

              for (int ncellrank = 0; ncellrank <= nmaxknn; ncellrank++) {
                if (theBaseDistRecA[ncellrank].ddist < Double.MAX_VALUE) {
                  // only taking positions which are valid
                  // gets the index of the closest instance
                  int theBaseDistRecA_ncellrank_ncell = theBaseDistRecA[ncellrank].ncell;
                  int global_nmarktarget_ncellrank = global_nmarktarget[ncellrank];

                  // prints the local then global nearest feature
                  // gives distance and then rank order
                  theInstance[ninstanceindex++] =
                      numformat(databinfirst_nbinindex_ntargetmark[theBaseDistRecA_ncellrank_ncell])
                          + "|"
                          + theBaseDistRecA_ncellrank_ncell;

                  // if (global_nmarktarget[ncellrank] ==-1)
                  // {
                  //  System.out.println("global_nmarktarget -1\t"+nmarktarget+"\t"+ncellrank);
                  // }
                  theInstance[ninstanceindex++] =
                      numformat(databinfirst_nbinindex_ntargetmark[global_nmarktarget_ncellrank])
                          + "|"
                          + global_nmarktarget_ncellrank;
                } else {
                  break;
                }
              }
              // ltimesort += System.currentTimeMillis()-ltimesortstart;
            }
          }
          // }

          // long ltimemarkstart = System.currentTimeMillis();

          boolean[] bcellmark_nholdoutcell = bcellmark[nholdoutcell];
          float[][] databinlastCELLMARK_nholdoutcell = databinlastCELLMARK[nholdoutcell];

          // for (int nmarkindex = 0; nmarkindex < markstargetcellA.length; nmarkindex++)
          for (int nmark = 0; nmark < marksA.length; nmark++) {
            // int nmark = markstargetcellA[nmarkindex];
            // going through all marks available in target cell

            if ((bcellmark_nholdoutcell[nmark]) && (nmark != ntargetmark)) {
              // we've had the mark and is not the target mark
              float[] databinlast_nmark_nholdoutcell = databinlastCELLMARK_nholdoutcell[nmark];

              // getting the center value
              theInstance[ninstanceindex++] = numformat(databinlast_nmark_nholdoutcell[nbinindex]);

              int nleftindex = nbinindex - nincrementnarrow;
              int nrightindex = nbinindex + nincrementnarrow;

              for (int ni = nincrementnarrow; ni <= nmaxoffsetnarrow; ni += nincrementnarrow) {
                // adding values to the left and right
                if (ninstanceindex >= theInstance.length) {
                  // gives some debug info
                  System.out.println(
                      "OUT theInstance[ninstanceindex++]\t"
                          + theInstance.length
                          + " \t"
                          + ninstanceindex
                          + "\t"
                          + nbinindex
                          + "\t"
                          + nincrementnarrow
                          + "\t"
                          + nmaxoffsetnarrow);
                }

                if (nleftindex >= databinlast_nmark_nholdoutcell.length) {
                  // gives some debug info
                  System.out.println(
                      "OUT databinlast_nmark_ntargetcell nleftindex\t"
                          + nleftindex
                          + "\t"
                          + nbinindex
                          + "\t"
                          + nincrementnarrow);
                }

                // adds the left narrow
                theInstance[ninstanceindex++] =
                    numformat(databinlast_nmark_nholdoutcell[nleftindex]);
                nleftindex -= nincrementnarrow;

                // adds the right narrow
                theInstance[ninstanceindex++] =
                    numformat(databinlast_nmark_nholdoutcell[nrightindex]);
                nrightindex += nincrementnarrow;
              }

              // starts at adding increment wide to narrow position
              nleftindex = nbinindex - nmaxoffsetnarrow - nincrementwide;
              nrightindex = nbinindex + nmaxoffsetnarrow + nincrementwide;

              for (int ni = nmaxoffsetnarrow + nincrementwide;
                  ni <= nmaxoffsetwide;
                  ni += nincrementwide) {
                // DNase signal to the left, right, and cumulative
                theInstance[ninstanceindex++] =
                    numformat(databinlast_nmark_nholdoutcell[nleftindex]);
                nleftindex -= nincrementwide;

                theInstance[ninstanceindex++] =
                    numformat(databinlast_nmark_nholdoutcell[nrightindex]);
                nrightindex += nincrementwide;
              }
            }
          }

          // ltimemark += System.currentTimeMillis()-ltimemarkstart;

          // long ltimeprintstart = System.currentTimeMillis();
          // writing the data to train a classifier out to disk

          StringBuffer sb = new StringBuffer();
          for (int nel = 0; nel < theInstance.length; nel++) {
            sb.append(theInstance[nel] + "\t");
          }

          // adds the target value to the string buffer
          sb.append(databinlast_ntargetmark_nholdoutcell[nbinindex] + "\n");

          byte[] btformat = sb.toString().getBytes();
          trainPW_nclassifierindex.write(btformat, 0, btformat.length);

          // ltimeprint += System.currentTimeMillis()-ltimeprintstart;
          // ltimeyes += System.currentTimeMillis()-ltimeyesstart;

          // moving to the next position for this classifier
          nsampleindexA[nclassifierindex]++;
          // also updating our local variable for it
          nsampleindexA_nclassifierindex++;
        }
        // moving which base we have sampled for this classifier
        noverallpos[nclassifierindex]++;
      }
      // }
      nclassifierindex++;
    }

    // System.out.println("Leaving
    // generateInstanceDataTrain\t"+(System.currentTimeMillis()-lstartoverall)+"\t"+ltimeinner+"\t"+ltimesort+"\t"
    //       +ltimeprint+"\t"+ltimeinit+"\t"+ltimeyes+"\t"+ltimemark);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  /** Stores into array list */
  public void generateInstanceDataTrainDNAMethyl(
      float[][][] databinlast,
      float[][][] databinlastCELLMARK,
      float[][][] databinfirst,
      int numbin,
      int nbegin,
      // boolean[][] bmarkcell, boolean[][] bcellmark, //int[] basemarksA, int[] basecellsA,
      // int[] markstargetcellA,
      int nholdoutcell,
      int nclassifierindex, // float[][] avgsignal,
      int nchromindex)
      throws Exception {
    // generateData(data,nbin,basemarksA,basecellsA, markstargetcellA,ntargetcell,ntargetmark);
    // features are target cell type

    // long lstartoverall = System.currentTimeMillis();
    // long ltimeprint = 0;
    // long ltimeyes = 0;

    // NumberFormat nf2 = NumberFormat.getInstance();
    // nf2.setMaximumFractionDigits(2);
    // nf2.setGroupingUsed(false);

    // int[][] savedorder = new int[markstargetcellA.length][nmaxknn+1];
    int[][] savedorder = new int[marksA.length][nmaxknn + 1];

    // long ltimeinner = 0;
    // long ltimesort = 0;
    // long ltimemark = 0;

    // int nusecelltarget = nholdoutcell;

    // System.out.println("HOLDOUT is "+nusecelltarget+" "+nholdoutcell);

    int[][] global = new int[distmarkcellAL.length][numcells];
    int[][] globalRank = new int[distmarkcellAL.length][numcells];
    for (int na = 0; na < global.length; na++) {
      for (int nb = 0; nb < global[na].length; nb++) {
        global[na][nb] = -1;
        globalRank[na][nb] = -1;
      }
    }

    // for (int ntargetmarkindex = 0; ntargetmarkindex < markstargetcellA.length;
    // ntargetmarkindex++)
    for (int nmarktarget = 0; nmarktarget < marksA.length; nmarktarget++) {
      // going through all marks in the target cell  //nmarktarget == nmarktargetcell
      // int nmarktarget = markstargetcellA[ntargetmarkindex];

      // gets actual index of target mark
      if (bmarkcell[nmarktarget][nholdoutcell]) {
        int npos = 0;
        for (int nindex = 0; nindex < distmarkcellAL[nmarktarget][nholdoutcell].size(); nindex++) {
          int ncell = ((Integer) distmarkcellAL[nmarktarget][nholdoutcell].get(nindex)).intValue();

          if ((ncell != nholdoutcell)
              && (bmarkcell[nmarktarget][ncell])
              && (bdnamethylcell[ncell])) {
            global[nmarktarget][npos] = ncell;
            globalRank[nmarktarget][ncell] = npos;
            npos++;
          }
        }
      }
    }

    // long ltimeinit = System.currentTimeMillis()-lstartoverall;

    for (int nbag = 0; nbag < numbags; nbag++) {
      // if ((nbagrequest == -1)||(nbagrequest == nbag))
      // {
      // no bag request or we have the bag request we want

      // ArrayList theCurrInstances;
      // ArrayList theCurrInstancesOutput = null;

      // theCurrInstances       = theTrainInstances[nclassifierindex];
      // theCurrInstancesOutput = theTrainInstancesOutput[nclassifierindex];
      GZIPOutputStream trainPW_nclassifierindex = trainPW[nclassifierindex];
      // int[] samples_nclassifierindex = nsampleindexA[nclassifierindex];

      int numfeatures = attributes[nclassifierindex].size();
      // System.out.println("getting
      // numfeatures\t"+nclassifierindex+"\t"+attributes.length+"\t"+numfeatures);

      String[] theInstance = new String[numfeatures];

      BaseDistRec[] theBaseDistRecA = new BaseDistRec[bmarkcell[0].length];

      for (int nindex = 0; nindex < theBaseDistRecA.length; nindex++) {
        theBaseDistRecA[nindex] = new BaseDistRec();
      }

      int[] samplesdnamethylbinindex_nbag = samplesdnamethylbinindex[nbag];
      int[] samplesdnamethylchromindex_nbag = samplesdnamethylchromindex[nbag];

      int nmaxPARAM = Math.max(nknnoffset, nmaxoffsetwide);
      // System.out.println("Doing "+nchromindex+" "+(nmaxPARAM+nbegin)+" to
      // "+(nbegin+numbin-nmaxPARAM-1));
      for (int nbinindex = nmaxPARAM; nbinindex < numbin - nmaxPARAM; nbinindex++) {
        int nsampleindexA_nclassifierindex = nsampleindexA[nclassifierindex];

        // moved databifirstinto if then

        while ((nsampleindexA_nclassifierindex < samplesdnamethylbinindex_nbag.length)
            && (samplesdnamethylbinindex_nbag[nsampleindexA_nclassifierindex] == nbegin + nbinindex)
            && (samplesdnamethylchromindex_nbag[nsampleindexA_nclassifierindex] == nchromindex)) {
          // we are on the chromosome position we have sampled

          // float[][] databinfirst_nbinindex = databinfirst[nbinindex];
          // this is a position we've sampled
          // long ltimeyesstart = System.currentTimeMillis();

          int ninstanceindex = 0;

          if (bdnamethylcell[nholdoutcell]) {
            // we have dna methylation as a target in this held out cell
            // for (int ntargetmarkindex = 0; ntargetmarkindex < markstargetcellA.length;
            // ntargetmarkindex++)
            for (int nmarktarget = 0; nmarktarget < marksA.length; nmarktarget++) {
              // going through all marks in the target cell
              int nhit = 0;
              // int nmarktarget = markstargetcellA[ntargetmarkindex];
              // gets actual index of target mark
              if (bmarkcell[nmarktarget][nholdoutcell]) {
                // order by mark is also available in held out cell
                boolean[] bmarkcell_nmarktarget =
                    bmarkcell[nmarktarget]; // array stating if mark is present in cell
                float[][] databinlast_nmarktarget = databinlast[nmarktarget];
                float[] databinlast_nmarktarget_nholdoutcell =
                    databinlast_nmarktarget[nholdoutcell];
                // float[] databinfirst_nbinindex_nmarktarget = databinfirst_nbinindex[nmarktarget];

                int[] global_nmarktarget = global[nmarktarget];
                int[] globalRank_nmarktarget = globalRank[nmarktarget];
                // array with data for target mark and cell

                // long ltimestart = System.currentTimeMillis();

                for (int ncell = 0; ncell < bmarkcell_nmarktarget.length; ncell++) {
                  if ((ncell != nholdoutcell)
                      && (bmarkcell_nmarktarget[ncell])
                      && (bdnamethylcell[ncell])) // ncell != ntargetcell
                  {
                    float ddist = 0;
                    float[] databinlast_nmarktarget_ncell = databinlast_nmarktarget[ncell];
                    int nchange = nbinindex - nknnoffset;
                    for (int npos = -nknnoffset; npos <= nknnoffset; npos++) {
                      float ddiff =
                          (databinlast_nmarktarget_ncell[nchange]
                              - databinlast_nmarktarget_nholdoutcell[nchange]);
                      nchange++;
                      ddist += ddiff * ddiff; // takes the squared error of the signal difference
                    }

                    BaseDistRec theBaseDistRecA_ncell = theBaseDistRecA[ncell];
                    theBaseDistRecA_ncell.ddist = ddist;
                    theBaseDistRecA_ncell.ncell = ncell;
                  } else {
                    // cell type is not eligible gives maximum value
                    BaseDistRec theBaseDistRecA_ncell = theBaseDistRecA[ncell];
                    theBaseDistRecA_ncell.ddist = Double.MAX_VALUE;
                    theBaseDistRecA_ncell.ncell = ncell;
                  }
                }

                // ltimeinner += System.currentTimeMillis()-ltimestart;

                // int[] savedorder_ntargetmarkindex = savedorder[ntargetmarkindex];
                int[] savedorder_nmarktarget = savedorder[nmarktarget];

                // long ltimesortstart = System.currentTimeMillis();
                // initializes the sorted order to approximately what was used at the last position
                // using just in place swaps
                for (int nrankpos = 0; nrankpos <= nmaxknn; nrankpos++) {
                  int nbestindex = savedorder_nmarktarget[nrankpos];
                  BaseDistRec temprec = theBaseDistRecA[nbestindex];
                  theBaseDistRecA[nbestindex] = theBaseDistRecA[nrankpos];
                  theBaseDistRecA[nrankpos] = temprec;
                }

                if (btieglobal) {
                  insertionsortFirstK(theBaseDistRecA, nmaxknn + 1, globalRank_nmarktarget);
                } else {
                  insertionsortFirstK(theBaseDistRecA, nmaxknn + 1); // , nbinindex);
                }

                for (int nrankpos = 0; nrankpos <= nmaxknn; nrankpos++) {
                  savedorder_nmarktarget[nrankpos] = theBaseDistRecA[nrankpos].ncell;
                }

                float[] samplesdnamethylvals_nbag_nsampleindexA_nclassifierindex =
                    samplesdnamethylvals[nbag][nsampleindexA[nclassifierindex]];

                int numknn = 0;

                for (int ncellrank = 0; ncellrank <= nmaxknn; ncellrank++) {
                  if (theBaseDistRecA[ncellrank].ddist < Double.MAX_VALUE) {
                    // we have a valid a cell type for ordering
                    int theBaseDistRecA_ncellrank_ncell = theBaseDistRecA[ncellrank].ncell;
                    int global_nmarktarget_ncellrank = global_nmarktarget[ncellrank];

                    // maps the regular cell index to the dna methylation index for the local and
                    // then global ordering
                    theInstance[ninstanceindex++] =
                        numformat(
                                samplesdnamethylvals_nbag_nsampleindexA_nclassifierindex[
                                    regularcelltodnamethylindex[theBaseDistRecA_ncellrank_ncell]])
                            + "|"
                            + theBaseDistRecA_ncellrank_ncell;

                    // if (global_nmarktarget[ncellrank] ==-1)
                    //   System.out.println("global_nmarktarget -1\t"+nmarktarget+"\t"+ncellrank);

                    theInstance[ninstanceindex++] =
                        numformat(
                                samplesdnamethylvals_nbag_nsampleindexA_nclassifierindex[
                                    regularcelltodnamethylindex[global_nmarktarget_ncellrank]])
                            + "|"
                            + global_nmarktarget_ncellrank;
                  } else {
                    break;
                  }
                }
                // ltimesort += System.currentTimeMillis()-ltimesortstart;
              }
            }
          }

          // long ltimemarkstart = System.currentTimeMillis();

          boolean[] bcellmark_nholdoutcell = bcellmark[nholdoutcell];
          float[][] databinlastCELLMARK_nholdoutcell = databinlastCELLMARK[nholdoutcell];

          // for (int nmarkindex = 0; nmarkindex < markstargetcellA.length; nmarkindex++)
          for (int nmark = 0; nmark < marksA.length; nmark++) {
            // int nmark = markstargetcellA[nmarkindex];
            // going through all marks

            if ((bcellmark_nholdoutcell[nmark]) && (nmark != ntargetmark)) {
              // checks if mark is available in the hold out cell and not target mark
              float[] databinlast_nmark_nholdoutcell = databinlastCELLMARK_nholdoutcell[nmark];
              theInstance[ninstanceindex++] = numformat(databinlast_nmark_nholdoutcell[nbinindex]);

              int nleftindex = nbinindex - nincrementnarrow;
              int nrightindex = nbinindex + nincrementnarrow;

              // reads in same cell type narrow features
              for (int ni = nincrementnarrow; ni <= nmaxoffsetnarrow; ni += nincrementnarrow) {
                // DNase signal to the left, right, and cumulative
                if (ninstanceindex >= theInstance.length)
                  System.out.println(
                      "OUT theInstance[ninstanceindex++]\t"
                          + theInstance.length
                          + " \t"
                          + ninstanceindex
                          + "\t"
                          + nbinindex
                          + "\t"
                          + nincrementnarrow
                          + "\t"
                          + nmaxoffsetnarrow);
                if (nleftindex >= databinlast_nmark_nholdoutcell.length)
                  System.out.println(
                      "OUT databinlast_nmark_nholdoutcell nleftindex\t"
                          + nleftindex
                          + "\t"
                          + nbinindex
                          + "\t"
                          + nincrementnarrow);

                theInstance[ninstanceindex++] =
                    numformat(databinlast_nmark_nholdoutcell[nleftindex]);
                nleftindex -= nincrementnarrow;
                theInstance[ninstanceindex++] =
                    numformat(databinlast_nmark_nholdoutcell[nrightindex]);
                nrightindex += nincrementnarrow;
              }

              nleftindex = nbinindex - nmaxoffsetnarrow - nincrementwide;
              nrightindex = nbinindex + nmaxoffsetnarrow + nincrementwide;
              // reads in same cell type wide features
              for (int ni = nmaxoffsetnarrow + nincrementwide;
                  ni <= nmaxoffsetwide;
                  ni += nincrementwide) {
                // DNase signal to the left, right, and cumulative
                theInstance[ninstanceindex++] =
                    numformat(databinlast_nmark_nholdoutcell[nleftindex]);
                nleftindex -= nincrementwide;
                theInstance[ninstanceindex++] =
                    numformat(databinlast_nmark_nholdoutcell[nrightindex]);
                nrightindex += nincrementwide;
              }
            }
          }
          // ltimemark += System.currentTimeMillis()-ltimemarkstart;

          // long ltimeprintstart = System.currentTimeMillis();
          // writing the data to train a classifier out to disk

          // outputs the contents of the feature
          StringBuffer sb = new StringBuffer();

          for (int nel = 0; nel < theInstance.length; nel++) {
            sb.append(theInstance[nel] + "\t");
          }
          sb.append(
              samplesdnamethylvals[nbag][nsampleindexA[nclassifierindex]][
                      regularcelltodnamethylindex[nholdoutcell]]
                  + "\n");
          // adds in the target feature value

          byte[] btformat = sb.toString().getBytes();
          trainPW_nclassifierindex.write(btformat, 0, btformat.length);
          // ltimeprint += System.currentTimeMillis()-ltimeprintstart;

          // increment the overall sampling position
          nsampleindexA[nclassifierindex]++;
          nsampleindexA_nclassifierindex++;
          // ltimeyes += System.currentTimeMillis()-ltimeyesstart;
        }
      }
      // }
      nclassifierindex++;
    }

    // System.out.println("Leaving
    // generateInstanceDataTrain\t"+(System.currentTimeMillis()-lstartoverall)+"\t"+ltimeinner+"\t"+ltimesort+"\t"
    //       +ltimeprint+"\t"+ltimeinit+"\t"+ltimeyes+"\t"+ltimemark);
  }

  /////////////////////////////////////////////////////////////////////////////////////////

  public String numformat(float fval) {
    String szformat = nf1.format(fval);

    if (szformat.startsWith("0.")) // assuming no negatives
    {
      return szformat.substring(1);
    } else {
      return szformat;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////

  public void loadInstancesFromFile( // int nclassifierindex
      // BufferedReader brtraindata,
      BufferedReader[] brtraindataA,
      BufferedReader brattributes,
      GZIPOutputStream pwuseattributes,
      int ntargetcell)
      throws IOException {

    // System.out.println("IN LOADINSTANCESFROMFILE
    // \t"+nclassifierindex+"\t"+theTrainInstances[nclassifierindex]);
    // ArrayList theCurrInstances = theTrainInstances[nclassifierindex];
    // ArrayList theCurrInstancesOutput = theTrainInstancesOutput[nclassifierindex];

    String szLine;
    boolean bfirst = true;
    ArrayList almaxranks = new ArrayList();
    String szprevmark = "";
    boolean bopen = true;
    int nlastrank = -1;

    int numfeatures = 0;

    ArrayList alfeatureSets = new ArrayList();
    ArrayList albpresentfeatureSets = new ArrayList();
    ArrayList alcurrset = new ArrayList();
    ArrayList alotherfeatures = new ArrayList();
    // ArrayList alotherfeaturesMark = new ArrayList();
    ArrayList albpresentotherfeatures = new ArrayList();

    HashSet hstargetmarks = new HashSet();

    boolean[] bcellmark_ntargetcell = bcellmark[ntargetcell];

    // getting into hstargetmarks all marks that are marks in the target cell and pioneer if we are
    // using pioneer marks
    for (int nmark = 0; nmark < bcellmark_ntargetcell.length; nmark++) {
      if ((bcellmark_ntargetcell[nmark])
          && ((szpioneermark == null) || (hspioneermarks.contains(marksA[nmark])))) {
        hstargetmarks.add(marksA[nmark]);
      }
    }

    // System.out.println("====>"+hstargetmarks.size());

    // OrderNarrow_H3K4me1_1_by_H3K27me3
    while ((szLine = brattributes.readLine()) != null) {
      // reading attribute file
      // temporary change space
      StringTokenizer st = new StringTokenizer(szLine, "\t");
      // StringTokenizer st = new StringTokenizer(szLine,"\t_ "); //-- comment out in v0.9.6
      st.nextToken(); // gets rid of feature index
      String szorder = st.nextToken();
      // if (szorder.startsWith("Order")) -- comment out in v0.9.6
      if ((szorder.startsWith("OrderNarrow_")) || (szorder.startsWith("OrderGlobal_"))) {
        szorder = szorder.substring(12); // added in v0.9.6
        // st.nextToken(); //skips output mark //rem
        int nbyindex = szorder.indexOf("_by_");
        StringTokenizer stu = new StringTokenizer(szorder.substring(0, nbyindex), "_");
        // need to prevent _by_ in mark names

        String szrank = null;
        while (stu.hasMoreTokens()) {
          szrank = stu.nextToken();
        }

        // int nrank = Integer.parseInt(st.nextToken()); //gets rank position
        int nrank = Integer.parseInt(szrank);
        // st.nextToken(); //skips by
        // String szmark = st.nextToken(); //gets order by mark
        // String szmark = st.nextToken(); //gets order by mark
        String szmark = szorder.substring(nbyindex + 4, szorder.length());

        if ((!szmark.equals(szprevmark)) && (!bfirst)) {
          // were on a different order by mark saving the last rank
          almaxranks.add(Integer.valueOf(nlastrank));

          if ((hstargetmarks.contains(szprevmark)) && (buseorderfeatures)) {
            // going to accept this feature mark is in target and we are using order by features
            albpresentfeatureSets.add(Boolean.valueOf(true));
          } else {
            // not using this set of features
            albpresentfeatureSets.add(Boolean.valueOf(false));
          }

          // storing set of features all associated with this order by amount
          alfeatureSets.add(alcurrset);

          // resets the current set of features
          alcurrset = new ArrayList();
        }

        alcurrset.add(szLine);

        // updates the previous mark and last rank
        szprevmark = szmark;
        bfirst = false;
        nlastrank = nrank;
      } else {

        // System.out.println(szLine);

        // will assume all other features come before regular features
        alotherfeatures.add(szLine);

        // temporary change space
        // st = new StringTokenizer(szLine,"\t_ "); //removed in v0.9.5
        st = new StringTokenizer(szLine, "\t"); // added in v0.9.6
        st.nextToken(); // skip index
        String szLineAtt = st.nextToken();
        String szmarkdirect;
        if (szLineAtt.endsWith("_center")) {
          szmarkdirect = szLineAtt.substring(0, szLineAtt.length() - 7);
        } else {
          int nlastindex =
              Math.max(szLineAtt.lastIndexOf("_left_"), szLineAtt.lastIndexOf("_right_"));
          szmarkdirect = szLineAtt.substring(0, nlastindex);
        }

        // String szmarkdirect = st.nextToken();//gets
        // alotherfeaturesMark.add(szmarkdirect); //adds mark to the set of marks being considered
        if ((hstargetmarks.contains(szmarkdirect)) && (busesamecellfeatures)) {
          // if mark is in target cell and using same cell type features increment number of
          // features selected
          numfeatures++;
          albpresentotherfeatures.add(Boolean.valueOf(true));
        } else {
          albpresentotherfeatures.add(Boolean.valueOf(false));
        }
      }
    }

    // adds in from the last feature from the order by set
    if ((hstargetmarks.contains(szprevmark)) && buseorderfeatures) {
      albpresentfeatureSets.add(Boolean.valueOf(true));
    } else {
      albpresentfeatureSets.add(Boolean.valueOf(false));
    }

    alfeatureSets.add(alcurrset);

    almaxranks.add(Integer.valueOf(nlastrank));

    brattributes.close();

    int nbrindex = 0;
    BufferedReader brtraindata = brtraindataA[nbrindex];
    nbrindex++;
    szLine = brtraindata.readLine();

    StringTokenizer st = new StringTokenizer(szLine, "\t");
    // st = new StringTokenizer(szLine, "\t");

    int numdistclass = 2; // parameter that we have two classes narrow and wide

    int[] maxcount = new int[almaxranks.size()];

    // going through all the max rank sets for order by
    for (int nmaxrankindex = 0; nmaxrankindex < maxcount.length; nmaxrankindex++) {
      // gets the numranks available for the current mark
      int numranks = ((Integer) almaxranks.get(nmaxrankindex)).intValue();
      if (((Boolean) albpresentfeatureSets.get(nmaxrankindex)).booleanValue()) {
        // using this order by feature set
        int nnarrowcount = 0;

        for (int nrank = 0; nrank < numranks; nrank++) {
          // going through each element
          String sztoken = st.nextToken();
          StringTokenizer stpipe = new StringTokenizer(sztoken, "|");
          stpipe.nextToken();
          int ncell = Integer.parseInt(stpipe.nextToken());
          if (ncell != ntargetcell) {
            // in order to accept ncell can't be target count
            if (nnarrowcount < nmaxknn) {
              // only accept the feature if a knn-feature
              nnarrowcount++;
              numfeatures +=
                  numdistclass; // increment feature count for num associated with this rank
              maxcount[nmaxrankindex] =
                  nnarrowcount; // updates the highest knn number for this rank position
            }
          }

          // flush out the other order by feature not needed to get count
          for (int nj = 2; nj <= numdistclass; nj++) {
            st.nextToken();
          }
        }
      } else {
        // not available order set flushing out the corresponding elements
        for (int nrank = 0; nrank < numranks; nrank++) {
          for (int nj = 1; nj <= numdistclass; nj++) {
            st.nextToken();
          }
        }
      }
    }

    int nuseindex = 0;
    for (int nmaxrankindex = 0; nmaxrankindex < maxcount.length; nmaxrankindex++) {
      // going through all rank sets
      ArrayList alset = (ArrayList) alfeatureSets.get(nmaxrankindex);

      for (int nk = 0; nk < maxcount[nmaxrankindex] * numdistclass; nk++) {
        // outputting knn features actually using
        StringTokenizer sta = new StringTokenizer((String) alset.get(nk), "\t");
        sta.nextToken(); // gets rid of index

        String szfeature = nuseindex + "\t" + sta.nextToken() + "\n";
        byte[] btformat = szfeature.getBytes();
        pwuseattributes.write(btformat, 0, btformat.length);

        nuseindex++;
      }
    }

    // printing out same cell type features
    for (int notherfeature = 0; notherfeature < alotherfeatures.size(); notherfeature++) {
      if (((Boolean) albpresentotherfeatures.get(notherfeature)).booleanValue()) {
        // feature is valid

        StringTokenizer sta =
            new StringTokenizer((String) alotherfeatures.get(notherfeature), "\t");
        sta.nextToken(); // get rid of index

        // update use index and then print feature
        String szfeature = nuseindex + "\t" + sta.nextToken() + "\n";
        byte[] btformat = szfeature.getBytes();
        pwuseattributes.write(btformat, 0, btformat.length);
        nuseindex++;
      }
    }

    do {
      // allocates memory for the instance
      float[] theInstance = new float[numfeatures];
      st = new StringTokenizer(szLine, "\t");
      int nnarrowindex = 0;
      int nglobalindex = 1;

      for (int nmaxrankindex = 0; nmaxrankindex < almaxranks.size(); nmaxrankindex++) {

        int numranks = ((Integer) almaxranks.get(nmaxrankindex)).intValue();
        if (((Boolean) albpresentfeatureSets.get(nmaxrankindex)).booleanValue()) {
          int nnarrowcount = 0;
          int nglobalordercount = 0;
          float dKNNsum = 0;
          float dKNNsumglobal = 0;
          for (int nrank = 0; nrank < numranks; nrank++) {
            String sztoken = st.nextToken();
            StringTokenizer stpipe = new StringTokenizer(sztoken, "|");
            float fval = Float.parseFloat(stpipe.nextToken());
            int ncell = Integer.parseInt(stpipe.nextToken());
            if (ncell != ntargetcell) {
              if (nnarrowcount < nmaxknn) {
                dKNNsum += fval;
                nnarrowcount++;
                double davg = dKNNsum / nnarrowcount;
                theInstance[nnarrowindex] = ((int) (nroundval * davg + 0.5)) / froundval;
                nnarrowindex += numdistclass; // increments to next place for narrowindex vals
              }
            }

            sztoken = st.nextToken();
            stpipe = new StringTokenizer(sztoken, "|");
            fval = Float.parseFloat(stpipe.nextToken());
            ncell = Integer.parseInt(stpipe.nextToken());

            if (ncell != ntargetcell) {
              if (nglobalordercount < nmaxknn) {
                dKNNsumglobal += fval;
                nglobalordercount++;
                double davg = dKNNsumglobal / nglobalordercount;
                theInstance[nglobalindex] = ((int) (nroundval * davg + 0.5)) / froundval;
                nglobalindex += numdistclass; // increments to next place for globalindex vals
              }
            }
          }
        } else {
          // not using this feature set flushing out values
          for (int nrank = 0; nrank < numranks; nrank++) {
            for (int nj = 1; nj <= numdistclass; nj++) {
              st.nextToken();
            }
          }
        }
      }

      int nfeature = nnarrowindex;
      int numpresentotherfeatures = albpresentotherfeatures.size();
      for (int notherfeature = 0; notherfeature < numpresentotherfeatures; notherfeature++) {
        // adding in the same cell type features
        String sztoken = st.nextToken();
        if (((Boolean) albpresentotherfeatures.get(notherfeature)).booleanValue()) {
          theInstance[nfeature++] = Float.parseFloat(sztoken);
        }
      }

      theTrainInstances.add(theInstance);

      // storing the output value
      theTrainInstancesOutput.add(Float.valueOf(st.nextToken()));
      // theTrainInstancesOutput.add(new Float(st.nextToken()));
      szLine = brtraindata.readLine();
      if ((szLine == null) && (nbrindex < brtraindataA.length)) {
        brtraindata = brtraindataA[nbrindex];
        szLine = brtraindata.readLine();
        nbrindex++;
      }
    } while (szLine != null);

    // System.out.println("LEAVING IN LOADINSTANCESFROMFILE
    // \t"+nclassifierindex+"\t"+theCurrInstances.size()+"\t"+theCurrInstancesOutput.size());
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  /** Stores into array list */
  public void generateInstanceDataTest(
      float[][][] databinlast,
      float[][][] databinlastCELLMARK,
      float[][][] databinfirst,
      int numbin,
      int nbegin,
      // boolean[][] bmarkcell, boolean[][] bcellmark,
      int[] markstargetcellA,
      int ntargetcell,
      int ntargetmark,
      // boolean btrain,
      float[] predictvals)
      throws Exception {
    // features are target cell type

    // ArrayList theCurrInstances;
    // ArrayList theCurrInstancesOutput=null;

    int[][] savedorder = new int[markstargetcellA.length][nmaxknn + 1];
    boolean[] bmarkcell_ntargetmark = bmarkcell[ntargetmark];

    // nusecelltarget = ntargetcell;

    BaseDistRec[][] theBaseDistRecA = new BaseDistRec[markstargetcellA.length][bmarkcell[0].length];

    float[][] databinlastCELLMARK_ntargetcell = databinlastCELLMARK[ntargetcell];

    for (int nindex = 0; nindex < theBaseDistRecA.length; nindex++) {
      for (int njindex = 0; njindex < theBaseDistRecA[nindex].length; njindex++) {
        theBaseDistRecA[nindex][njindex] = new BaseDistRec();
      }
    }

    // boolean bweka = false;

    // does classification without going through weka

    // long ltime1=0,ltime1_13=0,ltime1_135=0,ltime1_14=0,ltime1_142=0,
    //  ltime1_143=0,ltime1_145=0,ltime1_147=0, ltime1_15=0,ltime2=0,
    // ltime3=0,ltime1_1=0,ltime1_2=0,ltime1_3=0,ltime2_5 = 0;

    double[][] storeddist = new double[databinlast.length][databinlast[0].length]; // mark, cell

    float[][] theInstanceA = new float[attributes.length][];
    for (int na = 0; na < attributes.length; na++) {
      if (attributes[na] != null) {
        theInstanceA[na] = new float[attributes[na].size()];
      }
    }

    // loading global nearest cell type for each mark
    int[][] global = new int[distmarkcellAL.length][numcells];
    int[][] globalRank = new int[distmarkcellAL.length][numcells];
    for (int na = 0; na < global.length; na++) {
      for (int nb = 0; nb < global[na].length; nb++) {
        global[na][nb] = -1;
        globalRank[na][nb] = -1;
      }
    }

    for (int ntargetmarkindex = 0; ntargetmarkindex < markstargetcellA.length; ntargetmarkindex++) {
      // going through all marks in the target cell  //nmarktarget == nmarktargetcell
      int nmarktarget = markstargetcellA[ntargetmarkindex];
      // gets actual index of target mark

      if (bmarkcell[nmarktarget][ntargetcell]) {
        // we have the mark target in the desired cell type

        int npos = 0;
        for (int nindex = 0; nindex < distmarkcellAL[nmarktarget][ntargetcell].size(); nindex++) {
          int ncell = ((Integer) distmarkcellAL[nmarktarget][ntargetcell].get(nindex)).intValue();

          // for cell to be eligible can't be target cell and needs both target and order by mark
          if ((ncell != ntargetcell)
              && (bmarkcell[nmarktarget][ncell])
              && (bmarkcell[ntargetmark][ncell])) {
            global[nmarktarget][npos] = ncell;
            globalRank[nmarktarget][ncell] = npos;
            npos++;
          }
        }
      }
    }

    int nmaxPARAM = Math.max(nmaxoffsetwide, nknnoffset);

    // initialize the stored distance for each order mark compared to each other cell type starting
    // from nmaxPARAM-nknnoffset
    for (int ntargetmarkindex = 0; ntargetmarkindex < markstargetcellA.length; ntargetmarkindex++) {
      int nmarktarget = markstargetcellA[ntargetmarkindex]; // gets the index of a target mark

      float[][] databinlast_nmarktarget = databinlast[nmarktarget];
      double[] storeddist_nmarktarget = storeddist[nmarktarget];

      for (int ncell = 0; ncell < databinlast_nmarktarget.length; ncell++) {
        double ddist = 0;

        float[] databinlast_nmarktarget_ncell = databinlast_nmarktarget[ncell];
        float[] databinlast_nmarktarget_ntargetcell = databinlast_nmarktarget[ntargetcell];

        for (int npos = -nknnoffset; npos <= nknnoffset; npos++) {
          int nchange = nmaxPARAM + npos;

          double ddiff =
              (databinlast_nmarktarget_ncell[nchange]
                  - databinlast_nmarktarget_ntargetcell[nchange]);
          ddist += ddiff * ddiff;
        }

        storeddist_nmarktarget[ncell] = ddist;
      }
    }

    int numbin_minus_nmaxPARAM = numbin - nmaxPARAM;
    for (int nbinindex = nmaxPARAM; nbinindex < numbin_minus_nmaxPARAM; nbinindex++) {
      // long lbegin1 = System.currentTimeMillis();

      int nposleft = nbinindex - nknnoffset - 1;
      int nposright = nbinindex + nknnoffset;

      float[][] databinfirst_nbinindex = databinfirst[nbinindex];
      float[][] databinfirst_nposleft = null;
      float[][] databinfirst_nposright = null;

      if (nbinindex > nmaxPARAM) {
        databinfirst_nposleft = databinfirst[nposleft];
        databinfirst_nposright = databinfirst[nposright];
      }

      for (int ntargetmarkindex = 0;
          ntargetmarkindex < markstargetcellA.length;
          ntargetmarkindex++) {
        // not the first time through
        // long lbegin1_1 = System.currentTimeMillis();
        BaseDistRec[] theBaseDistRecA_ntargetmarkindex = theBaseDistRecA[ntargetmarkindex];

        int nmarktarget = markstargetcellA[ntargetmarkindex]; // gets the index of a target mark

        boolean[] bmarkcell_nmarktarget =
            bmarkcell[nmarktarget]; // gets what cells it is present in
        // float[] databinfirst_nbinindex_nmarktarget = databinfirst_nbinindex[nmarktarget]; //gets
        // the data associated with it

        double[] storeddist_nmarktarget = storeddist[nmarktarget];

        int[] globalRank_nmarktarget = globalRank[nmarktarget];

        if (nbinindex > nmaxPARAM) {
          // we are not on the first position
          float[] databinfirst_nposleft_nmarktarget = databinfirst_nposleft[nmarktarget];
          float[] databinfirst_nposright_nmarktarget = databinfirst_nposright[nmarktarget];
          float databinfirst_nposleft_nmarktarget_ntargetcell =
              databinfirst_nposleft_nmarktarget[ntargetcell];
          float databinfirst_nposright_nmarktarget_ntargetcell =
              databinfirst_nposright_nmarktarget[ntargetcell];

          int bmarkcell_nmarktarget_length = bmarkcell_nmarktarget.length;

          for (int ncell = 0; ncell < bmarkcell_nmarktarget_length; ncell++) {
            BaseDistRec theBaseDistRecA_ntargetmarkindex_ncell =
                theBaseDistRecA_ntargetmarkindex[ncell];
            theBaseDistRecA_ntargetmarkindex_ncell.ncell = ncell; // initialize cell index

            if ((ncell != ntargetcell)
                && (bmarkcell_nmarktarget[ncell])
                && (bmarkcell_ntargetmark[ncell])) {
              // this is a valid cell
              double ddiffleft =
                  (databinfirst_nposleft_nmarktarget[ncell]
                      - databinfirst_nposleft_nmarktarget_ntargetcell);
              double ddiffright =
                  (databinfirst_nposright_nmarktarget[ncell]
                      - databinfirst_nposright_nmarktarget_ntargetcell);

              // removes distribution from position now outside left boundary
              // and adds distribution position now inside right boundary
              theBaseDistRecA_ntargetmarkindex_ncell.ddist =
                  storeddist_nmarktarget[ncell] - ddiffleft * ddiffleft + ddiffright * ddiffright;

              // updated in 0.9.7 to -100
              if (theBaseDistRecA_ntargetmarkindex_ncell.ddist < -100) {
                System.out.println(
                    "NEG DIST\t"
                        + ddiffleft
                        + "\t"
                        + ddiffright
                        + "\t"
                        + numbin
                        + "\t"
                        + nbegin
                        + "\t"
                        + nbinindex
                        + "\t"
                        + ncell
                        + "\t"
                        + nmarktarget
                        + "\t"
                        + ntargetcell
                        + "\t"
                        + theBaseDistRecA_ntargetmarkindex_ncell.ddist
                        + "\t"
                        + databinfirst_nposleft_nmarktarget[ncell]
                        + "\t"
                        + databinfirst_nposright_nmarktarget[ncell]
                        + "\t"
                        + storeddist_nmarktarget[ncell]);

                throw new Exception();
              }

              // updates the stored distance for this cell type
              storeddist_nmarktarget[ncell] = theBaseDistRecA_ntargetmarkindex_ncell.ddist;
            } else {
              theBaseDistRecA_ntargetmarkindex_ncell.ddist = Double.MAX_VALUE;
            }
          }
        } else {
          // first time through going through all other cells
          for (int ncell = 0; ncell < bmarkcell_nmarktarget.length; ncell++) {

            BaseDistRec theBaseDistRecA_ntargetmarkindex_ncell =
                theBaseDistRecA_ntargetmarkindex[ncell];
            theBaseDistRecA_ntargetmarkindex_ncell.ncell = ncell; // initializing cell value

            if ((ncell != ntargetcell)
                && (bmarkcell_nmarktarget[ncell])
                && (bmarkcell_ntargetmark[ncell])) {
              // ncell not equal to target cell and both order by and target mark available in cell
              theBaseDistRecA_ntargetmarkindex_ncell.ddist = storeddist_nmarktarget[ncell];
            } else {
              // not an acceptable neighbor setting to max value
              theBaseDistRecA_ntargetmarkindex_ncell.ddist = Double.MAX_VALUE;
            }
          }
        }

        int[] savedorder_ntargetmarkindex = savedorder[ntargetmarkindex];

        // initializes the sorted order to approximately what was used last time using in place
        // swaps only
        for (int nrankpos = 0; nrankpos <= nmaxknn; nrankpos++) {
          int nbestindex = savedorder_ntargetmarkindex[nrankpos];
          BaseDistRec temprec = theBaseDistRecA_ntargetmarkindex[nbestindex];
          theBaseDistRecA_ntargetmarkindex[nbestindex] = theBaseDistRecA_ntargetmarkindex[nrankpos];
          theBaseDistRecA_ntargetmarkindex[nrankpos] = temprec;
        }

        // sorts the first maxknn+1 elements
        if (btieglobal) {
          insertionsortFirstK(
              theBaseDistRecA_ntargetmarkindex, nmaxknn + 1, globalRank_nmarktarget);
        } else {
          insertionsortFirstK(theBaseDistRecA_ntargetmarkindex, nmaxknn + 1); // , nbinindex);
        }

        for (int nrankpos = 0; nrankpos <= nmaxknn; nrankpos++) {
          savedorder_ntargetmarkindex[nrankpos] = theBaseDistRecA_ntargetmarkindex[nrankpos].ncell;
        }

        // long lbegin2 = System.currentTimeMillis();
        // ltime1 += (lbegin2-lbegin1);
      }

      int nclassifierindex = 0;

      // long lbegin2_5 = System.currentTimeMillis();
      // ltime2_5 += System.currentTimeMillis()-lbegin2_5;

      for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
        if (nholdoutcell != ntargetcell) {
          // NOTE: uncomment to debug which files it breaks on:
          // System.out.println("Computing cell: " + nholdoutcell 
          //         + " while targetcell is: " + ntargetcell);
          // System.out.println("current nclassifierindex is: " + nclassifierindex + " where length of instance was: " + theInstanceA.length);
          
          // theCurrInstances = theTestInstances[nclassifierindex];
          // if (((BLINEAR)&&(theClassifierLinearA[nclassifierindex]==null))||((!BLINEAR)&&
          // (theClassifierA[nclassifierindex] == null)))
          // {
          // System.out.println("Warning no classifier in "+nclassifierindex+" "+BLINEAR+"
          // "+theClassifierLinearA[nclassifierindex]);
          // }
          // else

          if (((!BLINEAR) && (theClassifierA[nclassifierindex] != null))
              || ((BLINEAR) && (theClassifierLinearA[nclassifierindex] != null))) {
            // may want to move this outside for loop
            float[] theInstance = theInstanceA[nclassifierindex]; 
            int ninstanceindexnarrow = 0;
            int ninstanceindexglobal = 1;
            int numdistclasses = 2;
            // long lbegin2_9 = System.currentTimeMillis();

            if (buseorderfeatures) {
              // if we are using order features
              for (int ntargetmarkindex = 0;
                  ntargetmarkindex < markstargetcellA.length;
                  ntargetmarkindex++) {
                int nmarktarget =
                    markstargetcellA[ntargetmarkindex]; // gets the index of a target mark

                boolean[] bmarkcell_nmarktarget =
                    bmarkcell[nmarktarget]; // gets what cells it is present in
                if ((bmarkcell_nmarktarget[nholdoutcell]) && (nmarktarget != ntargetmark)) {
                  // order by mark is different than targetmark and is available in the held out
                  // cell type

                  int[] global_nmarktarget = global[nmarktarget]; // gets the global ordering

                  // float[] databinfirst_nbinindex_nmarktarget =
                  // databinfirst_nbinindex[nmarktarget]; //gets the data associated with it
                  BaseDistRec[] theBaseDistRecA_ntargetmarkindex =
                      theBaseDistRecA[ntargetmarkindex];

                  float[] databinfirst_nbinindex_ntargetmark = databinfirst_nbinindex[ntargetmark];
                  int numknn = 0;

                  float dKNNsum = 0;
                  // for (int ncellrank = 0; numknn< nmaxknn; ncellrank++)
                  int ncellrank = 0;
                  while (numknn < nmaxknn) {
                    BaseDistRec theBaseDistRecA_ntargetmarkindex_ncellrank =
                        theBaseDistRecA_ntargetmarkindex[ncellrank];

                    if (theBaseDistRecA_ntargetmarkindex_ncellrank.ddist < Double.MAX_VALUE) {
                      if (theBaseDistRecA_ntargetmarkindex_ncellrank.ncell != nholdoutcell) {
                        // we don't use target cell used in training to be consistent in the number
                        // of possible knn features
                        dKNNsum +=
                            databinfirst_nbinindex_ntargetmark[
                                theBaseDistRecA_ntargetmarkindex_ncellrank.ncell];
                        numknn++; // increment knn info
                        double davg = dKNNsum / numknn;
                        theInstance[ninstanceindexnarrow] =
                            ((int) (nroundval * davg + 0.5)) / froundval;

                        ninstanceindexnarrow +=
                            numdistclasses; // increment narrow index of classifier
                      }
                    } else {
                      break;
                    }
                    ncellrank++; // increment rank position
                  }

                  numknn = 0;
                  dKNNsum = 0;
                  ncellrank = 0;
                  // for (int ncellrank = 0; numknn< nmaxknn; ncellrank++)
                  while (numknn < nmaxknn) {
                    int global_nmarktarget_ncellrank = global_nmarktarget[ncellrank];
                    if (global_nmarktarget_ncellrank != -1) {
                      // global at a valid distance
                      if (global_nmarktarget_ncellrank != nholdoutcell) {
                        // not tthe held out cell
                        dKNNsum += databinfirst_nbinindex_ntargetmark[global_nmarktarget_ncellrank];
                        numknn++;
                        double davg = dKNNsum / numknn;
                        theInstance[ninstanceindexglobal] =
                            ((int) (nroundval * davg + 0.5)) / froundval;
                        ninstanceindexglobal += numdistclasses;
                      }
                    } else {
                      break;
                    }
                    ncellrank++;
                  }
                }
              }
            }

            int ninstanceindex = ninstanceindexnarrow;

            if (busesamecellfeatures) {
              // we are using same cell features
              boolean[] bcellmark_nholdoutcell = bcellmark[nholdoutcell];

              for (int nmarkindex = 0; nmarkindex < markstargetcellA.length; nmarkindex++) {
                int nmark = markstargetcellA[nmarkindex];

                if ((bcellmark_nholdoutcell[nmark]) && (nmark != ntargetmark)) {
                  // mark was available in the holdout cell and not target mark

                  float[] databinlast_nmark_ntargetcell = databinlastCELLMARK_ntargetcell[nmark];

                  theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nbinindex];

                  int nleftindex = nbinindex - nincrementnarrow;
                  int nrightindex = nbinindex + nincrementnarrow;

                  for (int ni = nincrementnarrow; ni <= nmaxoffsetnarrow; ni += nincrementnarrow) {
                    // gets all values within nmaxoffsetnarrow window, incrementing by increment
                    // narrow
                    theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nleftindex];
                    nleftindex -= nincrementnarrow;
                    theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nrightindex];
                    nrightindex += nincrementnarrow;
                  }

                  nleftindex = nbinindex - nmaxoffsetnarrow - nincrementwide;
                  nrightindex = nbinindex + nmaxoffsetnarrow + nincrementwide;

                  for (int ni = nmaxoffsetnarrow + nincrementwide;
                      ni <= nmaxoffsetwide;
                      ni += nincrementwide) {
                    // getting the values in the wider index window
                    // TODO: Error here for some runs:
                    theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nleftindex];
                    nleftindex -= nincrementwide;
                    theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nrightindex];
                    nrightindex += nincrementwide;
                  }
                }
              }
            }

            // long lbegin3 = System.currentTimeMillis();
            // ltime2 += (lbegin3-lbegin2_9);

            int nbinindex_plus_nbegin = nbinindex + nbegin;
            for (int nbag = 0; nbag < numbags; nbag++) {

              double dpredictval = 0;

              if (BLINEAR) {
                double[] coeffs = theClassifierLinearA[nclassifierindex + nbag].coeffs;
                dpredictval = coeffs[0];
                for (int ni = 1; ni < coeffs.length; ni++) {
                  dpredictval += coeffs[ni] * theInstance[ni - 1];
                }
              } else {
                RegressionTree.TreeNode ptr = theClassifierA[nclassifierindex + nbag].theTree;
                // walks through the tree to find the value to classify it to
                // TODO Introduced try catch to see what the error is:
                try {
                    // go through tree, finding best end node:
                    while (ptr != null) {
                        dpredictval = ptr.dmean;
                        if (theInstance[ptr.nsplitfeatureindex] <= ptr.dsplitval) {  // TODO ERROR IS HERE
                            ptr = ptr.left;
                        } else {
                            ptr = ptr.right;
                        }
                    }
                } catch(Exception e){
                    System.out.println("Exception occurred");
                    // ptr is a treeNode, with: 
                    // nsplitfeatureindex, left, right, dsplitval, dmean
                    System.out.println("nbag: " + nbag + " out of numbags: " + numbags);
                    System.out.println("nclassind: " + nclassifierindex);
                    // Predicted value:
                    System.out.println("dpredval: " + dpredictval);
                    // theInstance.length MUST be larger than pr.nsplit.
                    System.out.println("theInstance.length: " + theInstance.length);
                    System.out.println("ptr.nsplit: " + ptr.nsplitfeatureindex);
                    System.out.println("ptr.dsplitval: " + ptr.dsplitval); 
                    System.out.println(e);
                    // Print values:
                    System.out.println(
                            "OUT of bounds " 
                            + nbinindex_plus_nbegin + "\t" + predictvals.length
                            + "\t" + nbinindex + "\t" + nbegin
                            + "\t" + numbin + "\t" + nmaxoffsetwide);
                    System.exit(1);
                }
              }
              if (nbinindex_plus_nbegin >= predictvals.length) {
                System.out.println(
                    "OUT of bounds " + nbinindex_plus_nbegin + "\t" + predictvals.length
                        + "\t" + nbinindex + "\t" + nbegin
                        + "\t" + numbin + "\t" + nmaxoffsetwide);
                System.exit(1);
              }

              predictvals[nbinindex_plus_nbegin] += dpredictval;
            }
            // long lbegin4 = System.currentTimeMillis();
            // ltime3 += (lbegin4-lbegin3);
          }
          nclassifierindex += numbags;
        }
      }
    } // goes through each bin
    // System.out.println("TIMEHEREGLOBAL\t"+ltime1+"\t"+ltime2+"\t"+ltime2_5+"\t"+ltime3);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** Stores into array list */
  public void generateInstanceDataTestDNAMethyl(
      float[][][] databinlast,
      float[][][] databinlastCELLMARK,
      float[][][] databinfirst,
      int numbin,
      int nbegin,
      // boolean[][] bmarkcell, boolean[][] bcellmark,
      int[] markstargetcellA,
      int ntargetcell,
      // boolean btrain,
      float[] predictvals,
      boolean[] presentvals,
      ArrayList almethylvals,
      ArrayList almethylcoord,
      int[] methylindexA)
      throws Exception {

    // ArrayList theCurrInstances;
    // ?ArrayList theCurrInstancesOutput=null;

    int[][] savedorder = new int[markstargetcellA.length][nmaxknn + 1];

    // int nusecelltarget;

    // nusecelltarget = ntargetcell;

    BaseDistRec[][] theBaseDistRecA =
        new BaseDistRec[markstargetcellA.length][bmarkcell[0].length]; // numbasedist];

    float[][] databinlastCELLMARK_ntargetcell = databinlastCELLMARK[ntargetcell];

    for (int nindex = 0; nindex < theBaseDistRecA.length; nindex++) {
      for (int njindex = 0; njindex < theBaseDistRecA[nindex].length; njindex++) {
        theBaseDistRecA[nindex][njindex] = new BaseDistRec();
      }
    }

    // boolean bweka = false;

    // does classification without going through weka

    // long ltime1=0,ltime1_13=0,ltime1_135=0,ltime1_14=0,ltime1_142=0,
    // ltime1_143=0,ltime1_145=0,ltime1_147=0, ltime1_15=0,ltime2=0,
    // ltime3=0,ltime1_1=0,ltime1_2=0,ltime1_3=0,ltime2_5 = 0;

    double[][] storeddist = new double[databinlast.length][databinlast[0].length];

    float[][] theInstanceA = new float[attributes.length][];
    for (int na = 0; na < attributes.length; na++) {
      if (attributes[na] != null) {
        theInstanceA[na] = new float[attributes[na].size()];
      }
    }

    int[][] global = new int[distmarkcellAL.length][numcells];
    int[][] globalRank = new int[distmarkcellAL.length][numcells];
    for (int na = 0; na < global.length; na++) {
      for (int nb = 0; nb < global[na].length; nb++) {
        global[na][nb] = -1;
        globalRank[na][nb] = -1;
      }
    }

    for (int ntargetmarkindex = 0; ntargetmarkindex < markstargetcellA.length; ntargetmarkindex++) {
      // going through all marks in the target cell  //nmarktarget == nmarktargetcell
      int nmarktarget = markstargetcellA[ntargetmarkindex];

      // gets actual index of target mark
      if (bmarkcell[nmarktarget][ntargetcell]) {
        int npos = 0;

        for (int nindex = 0; nindex < distmarkcellAL[nmarktarget][ntargetcell].size(); nindex++) {
          int ncell = ((Integer) distmarkcellAL[nmarktarget][ntargetcell].get(nindex)).intValue();

          if ((ncell != ntargetcell) && (bmarkcell[nmarktarget][ncell]) && bdnamethylcell[ncell]) {
            global[nmarktarget][npos] = ncell;
            globalRank[nmarktarget][ncell] = npos;
            npos++;
          }
        }
      }
    }

    int nmaxPARAM = Math.max(nmaxoffsetwide, nknnoffset);

    for (int ntargetmarkindex = 0; ntargetmarkindex < markstargetcellA.length; ntargetmarkindex++) {
      int nmarktarget = markstargetcellA[ntargetmarkindex]; // gets the index of a target mark

      float[][] databinlast_nmarktarget = databinlast[nmarktarget];
      double[] storeddist_nmarktarget = storeddist[nmarktarget];

      for (int ncell = 0; ncell < databinlast_nmarktarget.length; ncell++) {
        double ddist = 0;
        float[] databinlast_nmarktarget_ncell = databinlast_nmarktarget[ncell];
        float[] databinlast_nmarktarget_ntargetcell = databinlast_nmarktarget[ntargetcell];

        for (int npos = -nknnoffset; npos <= nknnoffset; npos++) {
          int nchange = nmaxPARAM + npos;

          double ddiff =
              (databinlast_nmarktarget_ncell[nchange]
                  - databinlast_nmarktarget_ntargetcell[nchange]);
          ddist += ddiff * ddiff;
        }

        storeddist_nmarktarget[ncell] = ddist;
      }
    }

    int numbin_minus_nmaxPARAM = numbin - nmaxPARAM;
    for (int nbinindex = nmaxPARAM; nbinindex < numbin_minus_nmaxPARAM; nbinindex++) {
      // long lbegin1 = System.currentTimeMillis();

      int nposleft = nbinindex - nknnoffset - 1;
      int nposright = nbinindex + nknnoffset;

      // float[][] databinfirst_nbinindex = databinfirst[nbinindex];
      float[][] databinfirst_nposleft = null;
      float[][] databinfirst_nposright = null;

      if (nbinindex > nmaxPARAM) {
        databinfirst_nposleft = databinfirst[nposleft];
        databinfirst_nposright = databinfirst[nposright];
      }

      for (int ntargetmarkindex = 0;
          ntargetmarkindex < markstargetcellA.length;
          ntargetmarkindex++) {
        // long lbegin1_1 = System.currentTimeMillis();
        BaseDistRec[] theBaseDistRecA_ntargetmarkindex = theBaseDistRecA[ntargetmarkindex];

        int nmarktarget = markstargetcellA[ntargetmarkindex]; // gets the index of a target mark

        boolean[] bmarkcell_nmarktarget =
            bmarkcell[nmarktarget]; // gets what cells it is present in
        // float[] databinfirst_nbinindex_nmarktarget = databinfirst_nbinindex[nmarktarget]; //gets
        // the data associated with it

        double[] storeddist_nmarktarget = storeddist[nmarktarget];
        int[] globalRank_nmarktarget = globalRank[nmarktarget];

        if (nbinindex > nmaxPARAM) {
          float[] databinfirst_nposleft_nmarktarget = databinfirst_nposleft[nmarktarget];
          float[] databinfirst_nposright_nmarktarget = databinfirst_nposright[nmarktarget];
          float databinfirst_nposleft_nmarktarget_ntargetcell =
              databinfirst_nposleft_nmarktarget[ntargetcell];
          float databinfirst_nposright_nmarktarget_ntargetcell =
              databinfirst_nposright_nmarktarget[ntargetcell];

          int bmarkcell_nmarktarget_length = bmarkcell_nmarktarget.length;

          for (int ncell = 0; ncell < bmarkcell_nmarktarget_length; ncell++) {
            BaseDistRec theBaseDistRecA_ntargetmarkindex_ncell =
                theBaseDistRecA_ntargetmarkindex[ncell];
            theBaseDistRecA_ntargetmarkindex_ncell.ncell = ncell;

            if ((ncell != ntargetcell)
                && (bmarkcell_nmarktarget[ncell])
                && (bdnamethylcell[ncell])) // (bmarkcell_ntargetmark[ncell]))
            {
              double ddiffleft =
                  (databinfirst_nposleft_nmarktarget[ncell]
                      - databinfirst_nposleft_nmarktarget_ntargetcell);
              double ddiffright =
                  (databinfirst_nposright_nmarktarget[ncell]
                      - databinfirst_nposright_nmarktarget_ntargetcell);

              theBaseDistRecA_ntargetmarkindex_ncell.ddist =
                  storeddist_nmarktarget[ncell] - ddiffleft * ddiffleft + ddiffright * ddiffright;

              // removing 0.9.7 to -100
              if (theBaseDistRecA_ntargetmarkindex_ncell.ddist < -100) {
                System.out.println(
                    "NEG DIST\t"
                        + ddiffleft + "\t" + ddiffright + "\t"
                        + numbin + "\t" + nbegin + "\t"
                        + nbinindex + "\t" + ncell + "\t"
                        + nmarktarget + "\t" + ntargetcell + "\t"
                        + theBaseDistRecA_ntargetmarkindex_ncell.ddist + "\t"
                        + databinfirst_nposleft_nmarktarget[ncell] + "\t"
                        + databinfirst_nposright_nmarktarget[ncell] + "\t"
                        + storeddist_nmarktarget[ncell]);
                throw new Exception();
              }

              storeddist_nmarktarget[ncell] = theBaseDistRecA_ntargetmarkindex_ncell.ddist;

            } else {
              theBaseDistRecA_ntargetmarkindex_ncell.ddist = Double.MAX_VALUE;
            }
          }
        } else {
          for (int ncell = 0; ncell < bmarkcell_nmarktarget.length; ncell++) {
            BaseDistRec theBaseDistRecA_ntargetmarkindex_ncell =
                theBaseDistRecA_ntargetmarkindex[ncell];
            theBaseDistRecA_ntargetmarkindex_ncell.ncell = ncell;

            if ((ncell != ntargetcell)
                && (bmarkcell_nmarktarget[ncell])
                && (bdnamethylcell[ncell])) {
              theBaseDistRecA_ntargetmarkindex_ncell.ddist = storeddist_nmarktarget[ncell];
            } else {
              theBaseDistRecA_ntargetmarkindex_ncell.ddist = Double.MAX_VALUE;
            }
          }
        }

        int[] savedorder_ntargetmarkindex = savedorder[ntargetmarkindex];

        // initializes the sorted order to approximately what was used at the last position using in
        // place swaps
        for (int nrankpos = 0; nrankpos < nmaxknn + 1; nrankpos++) {
          int nbestindex = savedorder_ntargetmarkindex[nrankpos];
          BaseDistRec temprec = theBaseDistRecA_ntargetmarkindex[nbestindex];
          theBaseDistRecA_ntargetmarkindex[nbestindex] = theBaseDistRecA_ntargetmarkindex[nrankpos];
          theBaseDistRecA_ntargetmarkindex[nrankpos] = temprec;
        }

        if (btieglobal) {
          insertionsortFirstK(
              theBaseDistRecA_ntargetmarkindex, nmaxknn + 1, globalRank_nmarktarget);
        } else {
          insertionsortFirstK(theBaseDistRecA_ntargetmarkindex, nmaxknn + 1); // , nbinindex);
        }

        for (int nrankpos = 0; nrankpos < nmaxknn + 1; nrankpos++) {
          savedorder_ntargetmarkindex[nrankpos] = theBaseDistRecA_ntargetmarkindex[nrankpos].ncell;
        }

        // long lbegin2 = System.currentTimeMillis();
        // ltime1 += (lbegin2-lbegin1);
      }

      // long lbegin2_5 = System.currentTimeMillis();

      // ltime2_5 += System.currentTimeMillis()-lbegin2_5;

      int nsizealmethylcoord = almethylcoord.size();
      int ncurrbin;

      System.out.println("methylindexA[0]: " + methylindexA[0]);
      System.out.println("nsizealmethylcoord: " + nsizealmethylcoord);
      // Run methyl predict?
      while ((methylindexA[0] < nsizealmethylcoord)
          && ((ncurrbin = (((Integer) almethylcoord.get(methylindexA[0])).intValue() / nresolution))
              <= nbegin + nbinindex)) {
        // iterates while we have more methylation index positions, and  there positions are less
        // than or equal to the current bin
        if (ncurrbin == (nbegin + nbinindex)) {
          int nclassifierindex = 0;
          for (int nholdoutcell = 0; nholdoutcell < numcells; nholdoutcell++) {
            // this is a cell type we can use for order by information
            if (nholdoutcell != ntargetcell) {
              // this is the bin of the
              float[] fmethylvals =
                  (float[])
                      almethylvals.get(
                          methylindexA[0]); // gets dna methylation values at that position
              // theCurrInstances = theTestInstances[nclassifierindex];

              // if (theClassifierA[nclassifierindex] == null)
              if (((!BLINEAR) && (theClassifierA[nclassifierindex] != null))
                  || ((BLINEAR) && (theClassifierLinearA[nclassifierindex] != null))) {

                float[] theInstance =
                    theInstanceA[
                        nclassifierindex]; // new float[numfeatures]; //need to move this outside
                // for loop
                int ninstanceindexnarrow = 0;

                int ninstanceindexglobal = 1;
                int numdistclasses = 2;
                // long lbegin2_9 = System.currentTimeMillis();

                if (buseorderfeatures) {
                  for (int ntargetmarkindex = 0;
                      ntargetmarkindex < markstargetcellA.length;
                      ntargetmarkindex++) {
                    int nmarktarget =
                        markstargetcellA[ntargetmarkindex]; // gets the index of a target mark

                    boolean[] bmarkcell_nmarktarget =
                        bmarkcell[nmarktarget]; // gets what cells it is present in
                    if (bmarkcell_nmarktarget[nholdoutcell]) {
                      int[] global_nmarktarget = global[nmarktarget];
                      // float[] databinfirst_nbinindex_nmarktarget =
                      // databinfirst_nbinindex[nmarktarget]; //gets the data associated with it
                      BaseDistRec[] theBaseDistRecA_ntargetmarkindex =
                          theBaseDistRecA[ntargetmarkindex];

                      // Get  the value here
                      int numknn = 0;

                      float dKNNsum = 0;
                      int ncellrank = 0;
                      // for (int ncellrank = 0; numknn< nmaxknn; ncellrank++)
                      while (numknn < nmaxknn) {
                        BaseDistRec theBaseDistRecA_ntargetmarkindex_ncellrank =
                            theBaseDistRecA_ntargetmarkindex[ncellrank];

                        if (theBaseDistRecA_ntargetmarkindex_ncellrank.ddist < Double.MAX_VALUE) {
                          if (theBaseDistRecA_ntargetmarkindex_ncellrank.ncell != nholdoutcell) {
                            dKNNsum +=
                                fmethylvals[
                                    regularcelltodnamethylindex[
                                        theBaseDistRecA_ntargetmarkindex_ncellrank.ncell]];

                            numknn++;
                            double davg = dKNNsum / numknn;
                            theInstance[ninstanceindexnarrow] =
                                ((int) (nroundval * davg + 0.5)) / froundval;
                            ninstanceindexnarrow += numdistclasses;
                          }
                        } else {
                          break;
                        }
                        ncellrank++;
                      }

                      numknn = 0;
                      dKNNsum = 0;
                      // for (int ncellrank = 0; numknn< nmaxknn; ncellrank++)
                      ncellrank = 0;
                      while (numknn < nmaxknn) {
                        int global_nmarktarget_ncellrank = global_nmarktarget[ncellrank];
                        if (global_nmarktarget_ncellrank != -1) {
                          if (global_nmarktarget_ncellrank != nholdoutcell) {
                            dKNNsum +=
                                fmethylvals[
                                    regularcelltodnamethylindex[global_nmarktarget_ncellrank]];
                            numknn++;
                            double davg = dKNNsum / numknn;
                            theInstance[ninstanceindexglobal] =
                                ((int) (nroundval * davg + 0.5)) / froundval;
                            ninstanceindexglobal += numdistclasses;
                          }
                        } else {
                          break;
                        }
                        ncellrank++;
                      }
                    }
                  }
                }

                int ninstanceindex = ninstanceindexnarrow;

                if (busesamecellfeatures) {
                  boolean[] bcellmark_nholdoutcell = bcellmark[nholdoutcell];
                  for (int nmarkindex = 0; nmarkindex < markstargetcellA.length; nmarkindex++) {
                    int nmark = markstargetcellA[nmarkindex];

                    if (bcellmark_nholdoutcell[nmark]) {
                      float[] databinlast_nmark_ntargetcell =
                          databinlastCELLMARK_ntargetcell[nmark];

                      theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nbinindex];

                      int nleftindex = nbinindex - nincrementnarrow;
                      int nrightindex = nbinindex + nincrementnarrow;

                      for (int ni = nincrementnarrow;
                          ni <= nmaxoffsetnarrow;
                          ni += nincrementnarrow) {
                        theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nleftindex];
                        nleftindex -= nincrementnarrow;

                        theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nrightindex];
                        nrightindex += nincrementnarrow;
                      }

                      nleftindex = nbinindex - nmaxoffsetnarrow - nincrementwide;
                      nrightindex = nbinindex + nmaxoffsetnarrow + nincrementwide;

                      for (int ni = nmaxoffsetnarrow + nincrementwide;
                          ni <= nmaxoffsetwide;
                          ni += nincrementwide) {
                        theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nleftindex];
                        nleftindex -= nincrementwide;
                        theInstance[ninstanceindex++] = databinlast_nmark_ntargetcell[nrightindex];
                        nrightindex += nincrementwide;
                      }
                    }
                  }
                }

                // long lbegin3 = System.currentTimeMillis();
                // ltime2 += (lbegin3-lbegin2_9);

                int nbinindex_plus_nbegin = nbinindex + nbegin;

                for (int nbag = 0; nbag < numbags; nbag++) {
                  double dpredictval = 0;

                  if (BLINEAR) {

                    double[] coeffs = theClassifierLinearA[nclassifierindex + nbag].coeffs;
                    dpredictval = coeffs[0];
                    for (int ni = 1; ni < coeffs.length; ni++) {
                      dpredictval += coeffs[ni] * theInstance[ni - 1];
                    }

                    predictvals[methylindexA[0]] += dpredictval;
                    presentvals[methylindexA[0]] = true;

                  } else {
                    RegressionTree.TreeNode ptr = theClassifierA[nclassifierindex + nbag].theTree;
                    while (ptr != null) {
                      // walks through the tree to find the value to classify it to
                      dpredictval = ptr.dmean;

                      if (theInstance[ptr.nsplitfeatureindex] <= ptr.dsplitval) {
                        ptr = ptr.left;
                      } else {
                        ptr = ptr.right;
                      }
                    }

                    predictvals[methylindexA[0]] += dpredictval;
                    presentvals[methylindexA[0]] = true;
                  }
                }
                // long lbegin4 = System.currentTimeMillis();
                // ltime3 += (lbegin4-lbegin3);
              } // else
              nclassifierindex += numbags;
            }
          }
        }
        methylindexA[0]++;
      }
    } // goes through each bin
    // System.out.println("TIMEHEREGLOBAL\t"+ltime1+"\t"+ltime2+"\t"+ltime2_5+"\t"+ltime3);//+"\t"+nglobalcount);

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  /** converts data into the desired binned resolution. and outputs a single chromosome per file */
  public void convertData() throws IOException {
    Iterator itrmark = (Iterator) hsmarks.iterator();

    String[][] markcellA =
        new String[hsmarks.size()][]; // stores for each mark all cell types mark is available for
    String[] markA = new String[hsmarks.size()]; // note not marksA

    int nmark = 0;

    while (itrmark.hasNext()) {
      // iterating over all marks
      String szmark = (String) itrmark.next();
      markA[nmark] = szmark;
      ArrayList alcell = (ArrayList) hmmarkcell.get(szmark);
      markcellA[nmark] = new String[alcell.size()];
      for (int nindex = 0; nindex < alcell.size(); nindex++) {
        // storing cell associated with mark
        markcellA[nmark][nindex] = (String) alcell.get(nindex);
      }
      nmark++;
    }

    NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
    nf.setMaximumFractionDigits(2);
    nf.setGroupingUsed(false);

    float[][] data = new float[alchrom.size()][];
    int nchrom;

    for (nchrom = 0; nchrom < data.length; nchrom++) {
      // allocates data array for each chromosome
      String szcurrchrom = (String) alchrom.get(nchrom);
      nchrom = ((Integer) hmchromindex.get(szcurrchrom)).intValue();
      // gets the size of the chromosome
      data[nchrom] =
          new float[(((Integer) hmchromsize.get(szcurrchrom)).intValue() - 1) / nresolution + 1];
    }
    // iterating over each chromosome

    // iterating over marks
    // System.out.println("\t"+markcellA.length+"\t"+markcellA[0].length);
    for (nmark = 0; nmark < markcellA.length; nmark++) {
      if ((szconvertmark != null) && (!markA[nmark].equals(szconvertmark))) {
        continue;
      }

      for (int ncell = 0; ncell < markcellA[nmark].length; ncell++) {

        if ((szconvertcell != null) && (!markcellA[nmark][ncell].equals(szconvertcell))) {
          continue;
        }

        // iterating over cell types of the marks
        // System.out.println("\t"+nmark+"\t"+ncell);
        for (int ni = 0; ni < data.length; ni++) {
          float[] data_ni = data[ni];
          for (int na = 0; na < data_ni.length; na++) {
            // resetting the contents of data to 0
            data_ni[na] = 0;
          }
        }

        // System.out.println("\t"+markcellA[nmark][ncell]);

        // iterating over cell types
        float fval = 0;
        int nstep = 1;
        int nspan = 1;
        int nposition = 0;
        String szcurrchrom = null;
        String szinfile =
            (String) hmmarkcellfile.get(markA[nmark] + "\t" + markcellA[nmark][ncell]);

        if ((szinfile.toLowerCase(Locale.ENGLISH).endsWith(".wig"))
            || (szinfile.toLowerCase(Locale.ENGLISH).endsWith(".wig.gz"))) {

          BufferedReader br;
          if (szchromwant == null) {
            File f = new File(szinputdir + "/" + szinfile);
            if (!f.exists()) {
              System.out.println(
                  "WARNING " + szinputdir + "/" + szinfile + " not found and not converted!");
              continue;
            }
            br = Util.getBufferedReader(szinputdir + "/" + szinfile);
          } else {
            File f = new File(szinputdir + "/" + szchromwant + "_" + szinfile);
            if (!f.exists()) {
              System.out.println(
                  "WARNING "
                      + szinputdir
                      + "/"
                      + szchromwant
                      + "_"
                      + szinfile
                      + " not found and not converted!");
              continue;
            }

            br = Util.getBufferedReader(szinputdir + "/" + szchromwant + "_" + szinfile);
          }

          boolean bdeclared = false;
          boolean bvariable = true;
          String szLine;

          while ((szLine = br.readLine()) != null) {
            if (szLine.startsWith("variableStep")) {
              bdeclared = true;
              StringTokenizer st = new StringTokenizer(szLine, " \t=");
              st.nextToken(); // variable
              st.nextToken(); // chrom
              szcurrchrom = st.nextToken();
              bvariable = true;
              if (st.hasMoreTokens()) {
                st.nextToken();
                nspan = Integer.parseInt(st.nextToken());
              }
            } else if (szLine.startsWith("fixedStep")) {
              bdeclared = true;
              bvariable = false;
              StringTokenizer st = new StringTokenizer(szLine, " \t=");
              st.nextToken(); // fixed
              st.nextToken(); // chrom
              szcurrchrom = st.nextToken();
              st.nextToken(); // start
              nposition = Integer.parseInt(st.nextToken()) - 1;
              st.nextToken();
              nstep = Integer.parseInt(st.nextToken());
              if (st.hasMoreTokens()) {
                st.nextToken();
                nspan = Integer.parseInt(st.nextToken());
              }
            } else if ((!szLine.startsWith("#"))
                && (!szLine.toLowerCase(Locale.ENGLISH).startsWith("browser"))
                && (!szLine.toLowerCase(Locale.ENGLISH).startsWith("track"))) {
              if (!bdeclared) {
                throw new IllegalArgumentException(
                    szinfile + " is a wig file but variable or fixed not declared");
              }

              Integer objInt = (Integer) hmchromindex.get(szcurrchrom);
              StringTokenizer st = new StringTokenizer(szLine, "\t ");
              if (bvariable) {
                if (st.countTokens() != 2) {
                  throw new IllegalArgumentException(
                      szinfile
                          + " is a variable step wig expecting 2 tokens but found "
                          + st.countTokens()
                          + " token, in line "
                          + szLine);
                }
              } else {
                if (st.countTokens() != 1) {
                  throw new IllegalArgumentException(
                      szinfile
                          + " is a fixed step wig expecting 1 token but found "
                          + st.countTokens()
                          + " token, in line "
                          + szLine);
                }
              }
              if (objInt != null) {
                nchrom = (objInt).intValue();

                if (bvariable) {

                  nposition = Integer.parseInt(st.nextToken()) - 1;
                  fval = Float.parseFloat(st.nextToken());
                } else {
                  fval = Float.parseFloat(szLine);
                }

                nchrom = ((Integer) hmchromindex.get(szcurrchrom)).intValue();
                float[] data_nchrom = data[nchrom];
                int nbegin = nposition / nresolution;

                // nposition
                int nend = (nposition + nspan) / nresolution;
                int nshort =
                    (nbegin + 1) * nresolution
                        - nposition
                        - nspan; // how short we are from filling the first bin
                if (nshort > 0) {
                  // span does not even fill the first position
                  data_nchrom[nbegin] +=
                      fval
                          * (nresolution - (nposition - nbegin * nresolution) - nshort)
                          / ((float) nresolution);
                } else {
                  data_nchrom[nbegin] +=
                      fval
                          * (nresolution - (nposition - nbegin * nresolution))
                          / ((float) nresolution);
                }

                if (nend > nbegin) {
                  // these are full bins
                  for (int nbin = nbegin + 1; nbin < nend; nbin++) {
                    // store full value
                    data_nchrom[nbin] += fval;
                  }

                  // amount to add for last bin
                  // if (nend >= data_nchrom.length)
                  // {
                  // updated
                  //   System.out.println("Invalid nend value "+nend+"\t"+szLine+"\t"+szcurrchrom);
                  // }

                  if (nend
                      < data_nchrom
                          .length) // added in 1.0.0 to handle chromosome length multiple of
                  // resolution
                  {
                    data_nchrom[nend] +=
                        fval * (nposition + nspan - nresolution * nend) / ((float) nresolution);
                  }
                }

                if (!bvariable) {
                  nposition += nstep;
                }
              }
            }
          }
          br.close();
        } else if ((szinfile.toLowerCase(Locale.ENGLISH).endsWith(".bedgraph"))
            || (szinfile.toLowerCase(Locale.ENGLISH).endsWith(".bedgraph.gz"))) {
          BufferedReader br;
          if (szchromwant == null) {
            File f = new File(szinputdir + "/" + szinfile);
            if (!f.exists()) {
              System.out.println(
                  "WARNING " + szinputdir + "/" + szinfile + " not found and not converted!");
              continue;
            }
            br = Util.getBufferedReader(szinputdir + "/" + szinfile);
          } else {
            File f = new File(szinputdir + "/" + szchromwant + "_" + szinfile);
            if (!f.exists()) {
              System.out.println(
                  "WARNING "
                      + szinputdir
                      + "/"
                      + szchromwant
                      + "_"
                      + szinfile
                      + " not found and not converted!");
              continue;
            }
            br = Util.getBufferedReader(szinputdir + "/" + szchromwant + "_" + szinfile);
          }

          String szLine;
          while ((szLine = br.readLine()) != null) {
            StringTokenizer st =
                new StringTokenizer(
                    szLine, "\t "); // added space here for delimiter in bedgraph files
            if ((!szLine.startsWith("#"))
                && (!szLine.toLowerCase(Locale.ENGLISH).startsWith("browser"))
                && (!szLine.toLowerCase(Locale.ENGLISH).startsWith("track"))) {
              // adding error checking the input is in expected format
              if (st.countTokens() != 4) {
                throw new IllegalArgumentException(
                    "Found a line "
                        + szLine
                        + " in a bedgraph file with "
                        + st.countTokens()
                        + " tokens while expecting 4");
              }
              szcurrchrom = st.nextToken();

              Integer objInt = (Integer) hmchromindex.get(szcurrchrom);
              if (objInt != null) {

                nchrom = (objInt).intValue();

                int nactualbegin = Integer.parseInt(st.nextToken());
                int nactualend = Integer.parseInt(st.nextToken()); // last base not include

                fval = Float.parseFloat(st.nextToken());
                int nbegin = nactualbegin / nresolution;
                int nend = nactualend / nresolution;
                int nshort = (nbegin + 1) * nresolution - nactualend;

                float[] data_nchrom = data[nchrom];

                if (nshort > 0) {
                  data_nchrom[nbegin] +=
                      fval
                          * (nresolution - (nactualbegin - nbegin * nresolution) - nshort)
                          / ((float) nresolution);
                } else {
                  data_nchrom[nbegin] +=
                      fval
                          * (nresolution - (nactualbegin - nresolution * nbegin))
                          / ((float) nresolution);
                }

                if (nend > nbegin) {
                  for (int nbin = nbegin + 1; nbin < nend; nbin++) {
                    data_nchrom[nbin] += fval;
                  }

                  if (nend
                      < data_nchrom
                          .length) // added in 1.0.0 to handle chromosome length multiple of
                  // resolution
                  {
                    data_nchrom[nend] +=
                        fval * (nactualend - nresolution * nend) / ((float) nresolution);
                  }
                }
              }
            }
          }
          br.close();
        }

        for (int ni = 0; ni < data.length; ni++) {
          GZIPOutputStream pw =
              new GZIPOutputStream(
                  new FileOutputStream(
                      szoutdir + "/" + alchrom.get(ni) + "_" + szinfile + ".wig.gz"));

          String szheader =
              "track type=wiggle_0 name="
                  + markcellA[nmark][ncell]
                  + "_"
                  + markA[nmark]
                  + "_observed\n"
                  + "fixedStep  chrom="
                  + alchrom.get(ni)
                  + " start=1 step="
                  + nresolution
                  + " span="
                  + nresolution
                  + "\n";
          byte[] btformat = szheader.getBytes();
          pw.write(btformat, 0, btformat.length);

          float[] data_ni = data[ni];
          int nlen =
              (((Integer) hmchromsize.get(alchrom.get(ni))).intValue() - 1) / nresolution
                  + 1; // changed by JE to subtract 1
          for (int nj = 0; nj < nlen; nj++) {
            String szvals = nf.format(data_ni[nj]) + "\n";
            btformat = szvals.getBytes();
            pw.write(btformat, 0, btformat.length);
          }
          pw.finish();
          pw.close();
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  public void computeGlobalCorrelations() throws IOException {

    if (szoutcell != null) {
      if (hmmarkcellfile.get(szoutmark + "\t" + szoutcell) != null) {
        computePairCorrelation(szoutmark, szoutcell);
      } else {
        System.out.println("Not found file for " + szoutmark + " " + szoutcell);
      }
      //          computeAllPairsCorrelation(szoutmark, szoutcell);
    } else if (szoutmark != null) {
      computeAllPairCorrelation(szoutmark);
    } else {
      for (int nmark = 0; nmark < marksA.length; nmark++) {
        computeAllPairCorrelation(marksA[nmark]);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////

  public void computeAllPairCorrelation(String sztargetmark) throws IOException {

    ArrayList alothercells = new ArrayList();

    for (int ncell = 0; ncell < cellsA.length; ncell++) {
      if (hmmarkcellfile.get(sztargetmark + "\t" + cellsA[ncell]) != null) {
        alothercells.add(cellsA[ncell]);
      }
    }

    long numvalues = 0;

    double[] dsumx = new double[alothercells.size()];
    double[] dsumxsq = new double[alothercells.size()];
    double[][] dsumxy = new double[alothercells.size()][alothercells.size()];
    double[] dvalsA = new double[alothercells.size()];
    BufferedReader[] brothercells = new BufferedReader[alothercells.size()];

    for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
      String szchrom = (String) alchrom.get(nchrom);

      // System.out.println(szchrom+"\t"+alothercells.size()+"\t"+sztargetmark);
      for (int ni = 0; ni < alothercells.size(); ni++) {
        brothercells[ni] =
            Util.getBufferedReader(
                szinputdir
                    + "/"
                    + szchrom
                    + "_"
                    + hmmarkcellfile.get(sztargetmark + "\t" + alothercells.get(ni))
                    + szextension);

        brothercells[ni].readLine();
        brothercells[ni].readLine();
      }

      String szLine;

      while ((szLine = brothercells[0].readLine()) != null) {
        dvalsA[0] = Double.parseDouble(szLine);
        for (int ni = 1; ni < dvalsA.length; ni++) {
          dvalsA[ni] = Double.parseDouble(brothercells[ni].readLine());
        }

        for (int ni = 0; ni < dvalsA.length; ni++) {
          double dval = dvalsA[ni];
          dsumx[ni] += dval;
          dsumxsq[ni] += dval * dval;
          double[] dsumxy_ni = dsumxy[ni];
          for (int nj = ni + 1; nj < dvalsA.length; nj++) {
            dsumxy_ni[nj] += dval * dvalsA[nj];
          }
        }
        numvalues++;

        // if (numvalues % 100000 == 0)
        //  System.out.println(numvalues);
      }

      for (int nj = 0; nj < brothercells.length; nj++) {
        brothercells[nj].close();
      }
    }

    //       double[] dist = new double[alothercells.size()-];
    for (int ncell1 = 0; ncell1 < alothercells.size(); ncell1++) {
      PrintWriter pw =
          new PrintWriter(szoutdir + "/" + alothercells.get(ncell1) + "_" + sztargetmark + ".txt");
      GlobalCorrRec[] recA = new GlobalCorrRec[alothercells.size() - 1];

      int nrecindex = 0;
      for (int ncell2 = 0; ncell2 < alothercells.size(); ncell2++) {
        if (ncell1 == ncell2) continue;

        double dcoeff;
        double dvarx = dsumxsq[ncell1] - dsumx[ncell1] * dsumx[ncell1] / numvalues;
        double dvary = dsumxsq[ncell2] - dsumx[ncell2] * dsumx[ncell2] / numvalues;
        double dvarxdvary = dvarx * dvary;

        if (dvarxdvary <= 0) {
          dcoeff = 0;
        } else {
          if (ncell1 < ncell2)
            dcoeff =
                (dsumxy[ncell1][ncell2] - dsumx[ncell1] * dsumx[ncell2] / numvalues)
                    / Math.sqrt(dvarxdvary);
          else
            dcoeff =
                (dsumxy[ncell2][ncell1] - dsumx[ncell1] * dsumx[ncell2] / numvalues)
                    / Math.sqrt(dvarxdvary);
        }

        double ddist = 1 - dcoeff;
        recA[nrecindex++] = new GlobalCorrRec((String) alothercells.get(ncell2), ddist);
      }
      Arrays.sort(recA, new GlobalCorrRecCompare());

      for (int ni = 0; ni < recA.length; ni++) {
        pw.println(recA[ni].szcell + "\t" + recA[ni].ddist);
      }
      pw.close();
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////

  // updated in v.0.9.6 to make static
  static class RecOtherCell {
    double dsumy;
    double dsumysq;
    double dsumxy;
    BufferedReader brothercells;
  }

  public void computePairCorrelation(String sztargetmark, String sztargetcell) throws IOException {
    PrintWriter pw = new PrintWriter(szoutdir + "/" + sztargetcell + "_" + sztargetmark + ".txt");

    ArrayList alothercells = new ArrayList();

    for (int ncell = 0; ncell < cellsA.length; ncell++) {
      if (!cellsA[ncell].equals(sztargetcell)) {
        if (hmmarkcellfile.get(sztargetmark + "\t" + cellsA[ncell]) != null) {
          alothercells.add(cellsA[ncell]);
        }
      }
    }

    long numvalues = 0;
    double dsumx = 0;
    // double[] dsumy = new double[alothercells.size()];
    double dsumxsq = 0;
    // double[] dsumysq = new double[alothercells.size()];
    // double[] dsumxy = new double[alothercells.size()];

    // BufferedReader[] brothercells = new BufferedReader[alothercells.size()];

    RecOtherCell[] theRecOtherCellA = new RecOtherCell[alothercells.size()];
    for (int nk = 0; nk < theRecOtherCellA.length; nk++) {
      theRecOtherCellA[nk] = new RecOtherCell();
    }

    for (int nchrom = 0; nchrom < alchrom.size(); nchrom++) {
      String szchrom = (String) alchrom.get(nchrom);
      // System.out.println(szchrom+"\t"+alothercells.size());
      for (int ni = 0; ni < alothercells.size(); ni++) {
        // theRecOtherCellA[ni] = new RecOtherCell();
        BufferedReader brothercells =
            Util.getBufferedReader(
                szinputdir
                    + "/"
                    + szchrom
                    + "_"
                    + hmmarkcellfile.get(sztargetmark + "\t" + alothercells.get(ni))
                    + szextension);

        brothercells.readLine();
        brothercells.readLine();
        theRecOtherCellA[ni].brothercells = brothercells;
      }

      BufferedReader brtargetcell =
          Util.getBufferedReader(
              szinputdir
                  + "/"
                  + szchrom
                  + "_"
                  + hmmarkcellfile.get(sztargetmark + "\t" + sztargetcell)
                  + szextension);
      brtargetcell.readLine();
      brtargetcell.readLine();

      String szLine;

      while ((szLine = brtargetcell.readLine()) != null) {
        double dtargetval = Double.parseDouble(szLine);

        dsumx += dtargetval;
        dsumxsq += dtargetval * dtargetval;

        for (int nj = 0; nj < theRecOtherCellA.length; nj++) {
          RecOtherCell theRecOtherCell = theRecOtherCellA[nj];
          double dval = Double.parseDouble(theRecOtherCell.brothercells.readLine());

          // double dval = Double.parseDouble(szLine);
          theRecOtherCell.dsumy += dval;
          theRecOtherCell.dsumysq += dval * dval;
          theRecOtherCell.dsumxy += dtargetval * dval;
        }
        numvalues++;
      }

      brtargetcell.close();
      for (int nj = 0; nj < theRecOtherCellA.length; nj++) {
        theRecOtherCellA[nj].brothercells.close();
      }
    }

    GlobalCorrRec[] recA = new GlobalCorrRec[alothercells.size()];
    for (int ni = 0; ni < recA.length; ni++) {

      double dcoeff;
      double dvarx = dsumxsq - dsumx * dsumx / numvalues;
      double dvary =
          theRecOtherCellA[ni].dsumysq
              - theRecOtherCellA[ni].dsumy * theRecOtherCellA[ni].dsumy / numvalues;
      double dvarxdvary = dvarx * dvary;
      if (dvarxdvary <= 0) {
        dcoeff = 0;
      } else {
        dcoeff =
            (theRecOtherCellA[ni].dsumxy - dsumx * theRecOtherCellA[ni].dsumy / numvalues)
                / Math.sqrt(dvarxdvary);
      }

      double ddist = 1 - dcoeff;
      recA[ni] = new GlobalCorrRec((String) alothercells.get(ni), ddist);
    }

    Arrays.sort(recA, new GlobalCorrRecCompare());

    for (int ni = 0; ni < recA.length; ni++) {
      pw.println(recA[ni].szcell + "\t" + recA[ni].ddist);
    }

    pw.close();
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  static class GlobalCorrRec {

    String szcell;
    double ddist;

    GlobalCorrRec(String szcell, double ddist) {
      this.szcell = szcell;
      this.ddist = ddist;
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  static class GlobalCorrRecCompare implements Comparator, Serializable {
    public int compare(Object o1, Object o2) {
      GlobalCorrRec r1 = (GlobalCorrRec) o1;
      GlobalCorrRec r2 = (GlobalCorrRec) o2;

      if (r1.ddist < r2.ddist) return -1;
      else if (r1.ddist > r2.ddist) return 1;
      else return 0;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  public static void main(String[] args) {

    /** This specifies the input directory */
    String szinputdir = null;

    /** This specifies the distance directory */
    String szdistancedir = null;

    /** This specifies the directory with the saved classifiers */
    String szclassifierdir = null;

    /** This specifies the directory where the imputed files should be written */
    String szoutdir = null;

    /** the file with information on the input for imputation */
    String szimputeinfoINfile = null;

    /** the file with information on the output for the imputation */
    // String szimputeinfoOUTfile = null;

    /** This specifies a single cell type to impute */
    String szoutcell = null;

    /** This specifies a single mark to impute */
    String szoutmark = null;

    /** contains info on the chromosomes */
    String szchrominfo = null;

    /** output file */
    String szoutfile = null;

    String szmethylheader = null;

    String szmethylinfo = null;

    String szmethylDIR = null;

    /** This parameter determines the output resolution */
    int nresolution = ChromImpute.DEFAULTRES;
    // -o outputfile - specifies the output imputation files
    // by default a full matrix is written when no file is specified
    // with the file name MARK_CELL

    // outputfile is cell, mark, outfile, target depth

    boolean bok = true;

    int nargindex = 0;

    String szcommand = "";

    if (args.length >= 1) {
      szcommand = args[nargindex++];
    }

    String szpioneermark = null;
    String szchromwant = null;
    String szchromwantgenerate = null;
    boolean busesamecellfeatures = true;
    boolean buseorderfeatures = true;

    int nmaxoffsetnarrow = ChromImpute.DEFAULTMAXOFFSETNARROW;
    int nmaxoffsetwide = ChromImpute.DEFAULTMAXOFFSETWIDE;
    int nincrementnarrow = ChromImpute.DEFAULTINCREMENTNARROW;
    int nincrementwide = ChromImpute.DEFAULTINCREMENTWIDE;
    int nknnoffset = ChromImpute.DEFAULTKNNWINDOW;

    if (szcommand.equalsIgnoreCase("Version")) {
      System.out.println("This is version 1.0.3 of ChromImpute");
    } else if (szcommand.equalsIgnoreCase("ExportToChromHMM")) {
      boolean bpartial = false;
      boolean busenames = false;
      int nbinsize = ChromImpute.DEFAULTCHROMHMMBIN;
      double dsignalthresh = Double.MAX_VALUE;
      // the convert command takes a target set of files
      // and converts them consistent with the desired resolution
      while (nargindex < args.length - 4) {
        // imputation file is CELL type, then mark, then data file
        // the -r option specifies the output resolution to be used

        String szoption;

        szoption = args[nargindex++];

        if (szoption.equals("-r")) {
          nresolution = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-g")) {
          dsignalthresh = Double.parseDouble(args[nargindex++]);
        } else if (szoption.equals("-b")) {
          nbinsize = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-partial")) {
          bpartial = true;
        } else if (szoption.equals("-usenames")) {
          busenames = true;
        } else {
          bok = false;
        }
      }

      if (nargindex == args.length - 4) {
        // input data
        szinputdir = args[nargindex++];
        szimputeinfoINfile = args[nargindex++];
        szchrominfo = args[nargindex++];
        szoutdir = args[nargindex];

        File f = new File(szoutdir);
        if (!f.exists()) {
          if (!f.mkdirs()) {
            // throw new IllegalArgumentException(szoutdir+" does not exist and could not be
            // created!");
            System.out.println("ERROR: " + szoutdir + " does not exist and could not be created!");
            System.exit(1);
          }
        }
      } else {
        bok = false;
      }

      if (!bok) {
        System.out.println(
            "USAGE: java ChromImpute ExportToChromHMM [-b chromhmmbinsize][-g signalthresh][-partial][-r resolution][-usenames] CHROMIMPUTEDIR inputinfofile chrominfofile CHROMHMMDIR");
        System.exit(1);
      } else {
        try {
          new ChromImpute(
              szchrominfo,
              szinputdir,
              szimputeinfoINfile,
              szoutdir,
              nresolution,
              dsignalthresh,
              nbinsize,
              bpartial,
              busenames);
        } catch (Exception ex) {
          ex.printStackTrace(System.out);
          System.exit(1);
        }
      }

    } else if (szcommand.equalsIgnoreCase("Convert")) {

      String szconvertmark = null;
      String szconvertcell = null;

      // the convert command takes a target set of files
      // and converts them consistent with the desired resolution
      while (nargindex < args.length - 4) {
        // imputation file is CELL type, then mark, then data file
        // default output is the full imputation grid
        // the -r option specifies the output resolution to be used

        // CELL type
        // files should be organized CELL type then directory

        // input is a cell type/mark text file indicating
        // if that cell type and mark should also
        String szoption;

        szoption = args[nargindex++];

        if (szoption.equals("-r")) {
          nresolution = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-c")) {
          szchromwant = args[nargindex++];
        } else if (szoption.equals("-m")) {
          szconvertmark = args[nargindex++];
        } else if (szoption.equals("-l")) {
          szconvertcell = args[nargindex++];
        } else {
          bok = false;
        }
      }

      if (nargindex == args.length - 4) {
        // input data
        szinputdir = args[nargindex++];
        szimputeinfoINfile = args[nargindex++];
        szchrominfo = args[nargindex++];
        szoutdir = args[nargindex];

        File f = new File(szoutdir);
        if (!f.exists()) {
          if (!f.mkdirs()) {
            // throw new IllegalArgumentException(szoutdir+" does not exist and could not be
            // created!");
            System.out.println("ERROR: " + szoutdir + " does not exist and could not be created!");
            System.exit(1);
          }
        }
      } else {
        bok = false;
      }

      if (!bok) {
        System.out.println(
            "USAGE: java ChromImpute Convert [-c chrom][-l convertsample][-m convertmark][-r resolution] INPUTDIR inputinfofile chrominfofile CONVERTEDDIR");
        System.exit(1);
      } else {
        try {
          new ChromImpute(
              szchrominfo,
              szinputdir,
              szimputeinfoINfile,
              szoutdir,
              nresolution,
              szchromwant,
              szconvertmark,
              szconvertcell);
        } catch (Exception ex) {
          ex.printStackTrace(System.out);
          System.exit(1);
        }
      }
    } else if (szcommand.equalsIgnoreCase("Apply")) {
      // String szvalidatefile = null;
      int nholdoutcellrequest = -1;
      int nbagrequest = -1;
      int nmaxknn = ChromImpute.DEFAULTMAXKNN;

      boolean bdnamethyl = ChromImpute.DEFAULTDNAMETHYL;
      boolean bprintonefile = false;
      boolean bprintbrowserheader = true; // ChromImpute.DEFAULTPRINTBROWSERHEADER;
      int numbags = ChromImpute.DEFAULTNUMBAGS;
      int nmintotalensemble = ChromImpute.DEFAULTMINTOTALENSEMBLE;

      boolean bmethylavggenome = false;
      boolean bmethylavgchrom = false;
      boolean btieglobal = false;
      // String szchromwant = null;

      while (nargindex < args.length - 8) {
        // imputation file is CELL type, then mark, then data file
        // default output is the full imputation grid
        // the -i option allows imputing a subset of files
        // the -s options allow specifying the specific cell type and mark to impute
        // the -t option allows specifying a target imputation depth, by default average is used
        // the -r option specifies the output resolution to be used

        // optionally specify a specific subset of marks

        // input are signal files in Wig or BigWig format
        // output file
        // The program takes a design matrix of
        // cell types and mark wig signals then imputes the
        // missing chromatin mark files

        // CELL type
        // files should be organized CELL type then directory

        // input is a cell type/mark text file indicating
        // if that cell type and mark should also
        String szoption;

        szoption = args[nargindex++];

        if (szoption.equals("-dnamethyl")) {
          bdnamethyl = true;
          szmethylinfo = args[nargindex++];
          szmethylDIR = args[nargindex++];
          szmethylheader = args[nargindex++];
        } else if (szoption.equals("-methylavggenome")) {
          bmethylavggenome = true;
        } else if (szoption.equals("-methylavgchrom")) {
          bmethylavgchrom = true;
        } else if (szoption.equals("-printonefile")) {
          bprintonefile = true;
        } else if (szoption.equals("-noprintbrowserheader")) {
          bprintbrowserheader = false;
        } else if (szoption.equals("-o")) {
          szoutfile = args[nargindex++];
        } else if (szoption.equals("-p")) {
          szpioneermark = args[nargindex++];
        } else if (szoption.equals("-a")) {
          nmintotalensemble = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-b")) {
          numbags = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-r")) {
          nresolution = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-s")) // here just for back compatiblity
        {
          // this specifies the target cell and mark
          szoutcell = args[nargindex++];
          szoutmark = args[nargindex++];
        } else if ((szoption.equals("-nosame")) || (szoption.equals("-markonly"))) {
          busesamecellfeatures = false;
        } else if ((szoption.equals("-noorder")) || (szoption.equals("-sampleonly"))) {
          buseorderfeatures = false;
        } else if (szoption.equals("-w")) {
          nmaxoffsetnarrow = Integer.parseInt(args[nargindex++]);
          nmaxoffsetwide = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-i")) {
          nincrementnarrow = Integer.parseInt(args[nargindex++]);
          nincrementwide = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-n")) {
          nknnoffset = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-c")) {
          szchromwant = args[nargindex++];
        } else if (szoption.equals("-k")) {
          nmaxknn = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-tieglobal")) {
          btieglobal = true;
        } else {
          bok = false;
        }
      }

      if ((bmethylavggenome) && (bmethylavgchrom)) {
        bok = false;
      }

      if (bmethylavgchrom && (szchromwant == null)) {
        System.out.println("-methylavgchrom flag present but no chromosome specified");
        bok = false;
      }

      // System.out.println(nargindex+"\t"+(args.length-4));
      if (nargindex == args.length - 8) {
        szinputdir = args[nargindex++];
        szdistancedir = args[nargindex++];
        szclassifierdir = args[nargindex++];
        szimputeinfoINfile = args[nargindex++];
        szchrominfo = args[nargindex++];
        szoutdir = args[nargindex++];

        File f = new File(szoutdir);
        if (!f.exists()) {
          if (!f.mkdirs()) {
            // throw new IllegalArgumentException(szoutdir+" does not exist and could not be
            // created!");
            System.out.println("ERROR: " + szoutdir + " does not exist and could not be created!");
            System.exit(1);
          }
        }
        szoutcell = args[nargindex++];
        szoutmark = args[nargindex++];

      } else {
        bok = false;
      }

      if (!bok) {
        System.out.println(
            "USAGE: java ChromImpute Apply [-a mintotalensemble][-b numbags][-c chrom][-sampleonly][-dnamethyl infofile directory header]"
                + "[-i incrementnarrow incrementwide][-k maxknn][-markonly][-methylavggenome|-methylavgchrom][-n knnwindow]"
                + "[-noprintbrowserheader][-o outputfile][-p selectedmarks][-printonefile][-r resolution][-tieglobal][-w windownarrow windowwide] "
                + "CONVERTEDDIR DISTANCEDIR PREDICTORDIR inputinfofile chrominfo OUTPUTIMPUTEDIR sample mark");
        System.exit(1);

      } else {
        try {
          new ChromImpute(
              szchrominfo,
              szinputdir,
              szdistancedir,
              szimputeinfoINfile,
              szoutdir,
              szoutcell,
              szoutmark, // szimputeinfoOUTfile,
              nresolution,
              szoutfile, // szvalidatefile,
              szclassifierdir,
              nmaxknn,
              szpioneermark,
              busesamecellfeatures,
              buseorderfeatures,
              szchromwant,
              bdnamethyl,
              nmintotalensemble,
              numbags,
              bprintbrowserheader,
              bprintonefile,
              szmethylheader,
              szmethylinfo,
              szmethylDIR,
              nmaxoffsetnarrow,
              nmaxoffsetwide,
              nincrementnarrow,
              nincrementwide,
              nknnoffset,
              bmethylavggenome,
              bmethylavgchrom,
              btieglobal);

        } catch (Exception ex) {
          ex.printStackTrace(System.out);
          System.exit(1);
        }
      }
    } else if ((szcommand.equalsIgnoreCase("ComputeGlobalDist"))
        || (szcommand.equalsIgnoreCase("ComputeGlobalCorr"))) {

      String szoption;
      String szextension = ChromImpute.DEFAULTEXTENSION;

      while (nargindex < args.length - 4) {
        szoption = args[nargindex++];

        if (szoption.equals("-s")) {
          // this specifies the target cell and mark
          szoutcell = args[nargindex++];
          szoutmark = args[nargindex++];
        } else if (szoption.equals("-m")) {
          szoutmark = args[nargindex++];
        } else if (szoption.equals("-x")) {
          szextension = args[nargindex++];
        } else if (szoption.equals("-r")) {
          nresolution = Integer.parseInt(args[nargindex++]);
        } else {
          bok = false;
        }
      }

      if (nargindex == args.length - 4) {
        szinputdir = args[nargindex++];
        szimputeinfoINfile = args[nargindex++];
        szchrominfo = args[nargindex++];
        szoutdir = args[nargindex++];

        File f = new File(szoutdir);
        if (!f.exists()) {
          if (!f.mkdirs()) {
            // throw new IllegalArgumentException(szoutdir+" does not exist and could not be
            // created!");
            System.out.println("ERROR: " + szoutdir + " does not exist and could not be created!");
            System.exit(1);
          }
        }
      } else {
        bok = false;
      }

      if (!bok) {
        System.out.println(
            "USAGE: java ChromImpute ComputeGlobalDist [-m mark][-r resolution][-s sample mark][-x extension] "
                + "CONVERTEDDIR inputinfofile chrominfo DISTANCEDIR");
        System.exit(1);
      } else {
        try {
          new ChromImpute(
              szchrominfo,
              szinputdir,
              szimputeinfoINfile,
              szoutmark,
              szoutcell,
              szoutdir,
              nresolution,
              szextension);
        } catch (Exception ex) {
          ex.printStackTrace(System.out);
          System.exit(1);
        }
      }
    } else if (szcommand.equalsIgnoreCase("Eval")) {

      String szoption;

      String szevalobserveddir = null;
      String szevalobservedfile = null;
      String szevalimputedir = null;
      String szevalimputefile = null;
      String szevaloutfile = null;
      String szpeakevalfile = null;

      double devalpercent1 = ChromImpute.DEFAULT_EVALPERCENT1;
      double devalpercent2 = ChromImpute.DEFAULT_EVALPERCENT2;
      boolean bprintbrowserheader = true;
      boolean bprintonefile = false;

      while (nargindex < args.length - 5) {
        szoption = args[nargindex++];

        if (szoption.equals("-o")) {
          szevaloutfile = args[nargindex++];
        } else if (szoption.equals("-p")) {
          devalpercent1 = Double.parseDouble(args[nargindex++]);
          devalpercent2 = Double.parseDouble(args[nargindex++]);
        } else if (szoption.equals("-noprintbrowserheader")) {
          bprintbrowserheader = false;
        } else if (szoption.equals("-f")) {
          szpeakevalfile = args[nargindex++];
        } else if (szoption.equals("-printonefile")) {
          bprintonefile = true;
        }
      }

      if (nargindex == args.length - 5) {
        szevalobserveddir = args[nargindex++];
        szevalobservedfile = args[nargindex++];
        szevalimputedir = args[nargindex++];
        szevalimputefile = args[nargindex++];
        szchrominfo = args[nargindex++];
      } else {
        bok = false;
      }

      if (!bok) {
        System.out.println(
            "USAGE: java ChromImpute Eval [-f peakevalfile][-noprintbrowserheader][-o outfile][-p percent1 percent2][-printonefile] CONVERTEDDIR ConvertedFile IMPUTEDIR ImputeFile chrominfo");
        System.exit(1);
      } else {
        try {
          new ChromImpute(
              szevalobserveddir,
              szevalobservedfile,
              szevalimputedir,
              szevalimputefile,
              szchrominfo,
              devalpercent1,
              devalpercent2,
              bprintbrowserheader,
              bprintonefile,
              szevaloutfile,
              szpeakevalfile);
        } catch (Exception ex) {
          ex.printStackTrace(System.out);
          System.exit(1);
        }
      }
    } else if ((szcommand.equalsIgnoreCase("GenerateTrainData"))
        || (szcommand.equalsIgnoreCase("GenerateFeatures"))) // for backwards compability
    {
      // generate features
      // train the model

      // String szvalidatefile = null;
      int nholdoutcellrequest = -1;
      int nbagrequest = -1;

      int nseed = ChromImpute.DEFAULTSEED;

      int nmaxknn = ChromImpute.DEFAULTMAXKNN;
      int numsamples = ChromImpute.DEFAULTNUMSAMPLES;

      boolean bdnamethyl = ChromImpute.DEFAULTDNAMETHYL;
      // int nminnumlocations = ChromImpute.DEFAULTMINNUMLOCATIONS;

      int numbags = ChromImpute.DEFAULTNUMBAGS;
      int nmintotalensemble = ChromImpute.DEFAULTMINTOTALENSEMBLE;

      boolean bmethylavggenome = false;
      boolean bmethylavgchrom = false;
      boolean btieglobal = false;
      // String szchromwant = null;

      while (nargindex < args.length - 6) {
        // imputation file is CELL type, then mark, then data file
        // default output is the full imputation grid
        // the -i option allows imputing a subset of files
        // the -s options allow specifying the specific cell type and mark to impute
        // the -t option allows specifying a target imputation depth, by default average is used
        // the -r option specifies the output resolution to be used

        // optionally specify a specific subset of marks

        // input are signal files in Wig or BigWig format
        // output file
        // The program takes a design matrix of
        // cell types and mark wig signals then imputes the
        // missing chromatin mark files

        // CELL type
        // files should be organized CELL type then directory

        // input is a cell type/mark text file indicating
        // if that cell type and mark should also
        String szoption;

        szoption = args[nargindex++];

        // if (szoption.equals("-i"))
        // {
        // this specifies a file with a list of cell types and marks to generate outputs
        //   szimputeinfoOUTfile = args[nargindex++];
        // }
        // else
        if (szoption.equals("-dnamethyl")) {
          bdnamethyl = true;
          szmethylinfo = args[nargindex++];
          szmethylDIR = args[nargindex++];
          szmethylheader = args[nargindex++];
        } else if (szoption.equals("-methylavggenome")) {
          bmethylavggenome = true;
        } else if (szoption.equals("-methylavgchrom")) {
          bmethylavgchrom = true;
        }
        // else if (szoption.equals("-o"))
        // {
        //  szoutfile = args[nargindex++];
        // }
        else if (szoption.equals("-f")) {
          numsamples = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-a")) {
          nmintotalensemble = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-b")) {
          numbags = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-d")) {
          nseed = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-r")) {
          nresolution = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-c")) {
          szchromwantgenerate = args[nargindex++];
        }
        /*
           else if ((szoption.equals("-s"))||(szoption.equals("-m")))//-s for backward compatiblity
           {
        //this specifies the target cell and mark
        //szoutcell = args[nargindex++];
        szoutmark = args[nargindex++];
           }
           */
        else if (szoption.equals("-k")) {
          nmaxknn = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-w")) {
          nmaxoffsetnarrow = Integer.parseInt(args[nargindex++]);
          nmaxoffsetwide = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-i")) {
          nincrementnarrow = Integer.parseInt(args[nargindex++]);
          nincrementwide = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-n")) {
          nknnoffset = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-tieglobal")) {
          btieglobal = true;
        }
        // -w windownarrow windowwide
        // else if (szoption.equals("-q"))
        // {
        //      nholdoutcellrequest = Integer.parseInt(args[nargindex++]);
        //      nbagrequest =  Integer.parseInt(args[nargindex++]);
        // }
        else {
          bok = false;
        }
      }

      // System.out.println(nargindex+"\t"+(args.length-4));
      if (nargindex == args.length - 6) {
        szinputdir = args[nargindex++];
        szdistancedir = args[nargindex++];
        szimputeinfoINfile = args[nargindex++];
        szchrominfo = args[nargindex++];
        szoutdir = args[nargindex++];

        File f = new File(szoutdir);
        if (!f.exists()) {
          if (!f.mkdirs()) {
            // throw new IllegalArgumentException(szoutdir+" does not exist and could not be
            // created!");
            System.out.println("ERROR: " + szoutdir + " does not exist and could not be created!");
            System.exit(1);
          }
        }

        szoutmark = args[nargindex++];
      } else {
        bok = false;
      }

      if ((bmethylavggenome) && (bmethylavgchrom)) {
        bok = false;
      }

      if (bmethylavgchrom && (szchromwantgenerate == null)) {
        System.out.println("-methylavgchrom flag present but no chromosome specified");
        bok = false;
      }

      if (!bok) {
        System.out.println(
            "USAGE: java ChromImpute GenerateTrainData [-a mintotalensemble][-b numbags][-c chrom][-d seed]"
                + "[-dnamethyl infofile directory header][-f numsamples][-i incrementnarrow incrementwide]"
                + "[-k maxknn][-methylavggenome|-methylavgchrom][-n knnwindow][-r resolution][-tieglobal][-w windownarrow windowwide]"
                + " CONVERTEDDIR DISTANCEDIR inputinfofile chrominfo TRAINDATADIR mark");
        System.exit(1);
      } else {
        try {
          new ChromImpute(
              szchromwantgenerate,
              szchrominfo,
              szinputdir,
              szdistancedir,
              szimputeinfoINfile,
              szoutdir,
              szoutcell,
              szoutmark, // szimputeinfoOUTfile,
              nresolution,
              szoutfile,
              // nholdoutcellrequest, nbagrequest,//true,false,
              nseed,
              nmaxknn,
              bdnamethyl,
              numsamples,
              nmintotalensemble,
              numbags,
              szmethylheader,
              szmethylinfo,
              szmethylDIR,
              nmaxoffsetnarrow,
              nmaxoffsetwide,
              nincrementnarrow,
              nincrementwide,
              nknnoffset,
              bmethylavggenome,
              bmethylavgchrom,
              btieglobal);
        } catch (Exception ex) {
          ex.printStackTrace(System.out);
          System.exit(1);
        }
      }
    } else if (szcommand.equalsIgnoreCase("Train")) // impute
    {
      // train the model

      // String szvalidatefile = null;
      int nholdoutcellrequest = -1;
      int nbagrequest = -1;
      // boolean bloadtrainfile = true; //false;
      int nmaxknn = ChromImpute.DEFAULTMAXKNN;
      // int nseed = ChromImpute.DEFAULTSEED;
      boolean bdnamethyl = ChromImpute.DEFAULTDNAMETHYL;
      int nminnumlocations = ChromImpute.DEFAULTMINNUMLOCATIONS;
      int numbags = ChromImpute.DEFAULTNUMBAGS;
      int nmintotalensemble = ChromImpute.DEFAULTMINTOTALENSEMBLE;

      while (nargindex < args.length - 5) {
        // imputation file is CELL type, then mark, then data file
        // default output is the full imputation grid
        // the -i option allows imputing a subset of files
        // the -s options allow specifying the specific cell type and mark to impute
        // the -t option allows specifying a target imputation depth, by default average is used
        // the -r option specifies the output resolution to be used

        // optionally specify a specific subset of marks

        // input are signal files in Wig or BigWig format
        // output file
        // The program takes a design matrix of
        // cell types and mark wig signals then imputes the
        // missing chromatin mark files

        // CELL type
        // files should be organized CELL type then directory

        // input is a cell type/mark text file indicating
        // if that cell type and mark should also
        String szoption;

        szoption = args[nargindex++];

        /*
           if (szoption.equals("-i"))
           {
        //this specifies a file with a list of cell types and marks to generate outputs
        szimputeinfoOUTfile = args[nargindex++];
           }
           else if (szoption.equals("-o"))
           {
           szoutfile = args[nargindex++];
           }
           else if (szoption.equals("-f"))
           {
           NUMSAMPLES = Integer.parseInt(args[nargindex++]);
           }
           */
        if (szoption.equals("-b")) {
          numbags = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-a")) {
          nmintotalensemble = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-m")) {
          nminnumlocations = Integer.parseInt(args[nargindex++]);
        }
        // else if (szoption.equals("-d"))
        // {
        //  nseed = Integer.parseInt(args[nargindex++]);
        // }
        // else if (szoption.equals("-r"))
        // {
        // nresolution = Integer.parseInt(args[nargindex++]);
        // }
        else if (szoption.equals("-p")) {
          szpioneermark = args[nargindex++];
        }
        /*
           else if (szoption.equals("-s"))
           {
        //this specifies the target cell and mark
        szoutcell = args[nargindex++];
        szoutmark = args[nargindex++];
           }
           */
        else if (szoption.equals("-q")) {
          nholdoutcellrequest = Integer.parseInt(args[nargindex++]);
        } else if (szoption.equals("-g")) {
          nbagrequest = Integer.parseInt(args[nargindex++]);
        }
        /*
        else if (szoption.equals("-v"))
        {
        szvalidatefile = args[nargindex++];
        }
        */
        else if (szoption.equals("-k")) {
          nmaxknn = Integer.parseInt(args[nargindex++]);
        }
        // else if (szoption.equals("-generate"))
        // {
        //  bloadtrainfile = false;
        // }
        else if (szoption.equals("-dnamethyl")) {
          bdnamethyl = true;
          szmethylheader = args[nargindex++];
        } else if ((szoption.equals("-nosame")) || (szoption.equals("-markonly"))) {
          busesamecellfeatures = false;
        } else if ((szoption.equals("-noorder")) || (szoption.equals("-sampleonly"))) {
          buseorderfeatures = false;
        } else {
          bok = false;
        }
      }

      // System.out.println(nargindex+"\t"+(args.length-4));
      if (nargindex == args.length - 5) {
        szinputdir = args[nargindex++];
        // szdistancedir = args[nargindex++];
        szimputeinfoINfile = args[nargindex++];
        // szchrominfo = args[nargindex++];
        szoutdir = args[nargindex++];

        File f = new File(szoutdir);
        if (!f.exists()) {
          if (!f.mkdirs()) {
            // throw new IllegalArgumentException(szoutdir+" does not exist and could not be
            // created!");
            System.out.println("ERROR: " + szoutdir + " does not exist and could not be created!");
            System.exit(1);
          }
        }

        szoutcell = args[nargindex++];
        szoutmark = args[nargindex++];
      } else {
        bok = false;
      }

      if (!bok) {
        System.out.println(
            "USAGE: java ChromImpute Train [-a mintotalensemble][-b numbags][-sampleonly][-dnamethyl header][-g bagrequest][-k maxknn][-m minnumpoints][-markonly][-p selectedmarks][-q samplerequest] TRAINDATADIR inputinfofile PREDICTORDIR sample mark");
        System.exit(1);
      } else {
        try {
          // for (int nholdoutcellindex = 0; nholdoutcellindex <= nholdoutcellrequest;
          // nholdoutcellindex++)
          {
            try {
              System.out.println("Started reading file:");

              new ChromImpute(
                  szchrominfo,
                  szinputdir,
                  szimputeinfoINfile,
                  szoutdir,
                  szoutcell,
                  szoutmark,
                  // szimputeinfoOUTfile, nresolution, szoutfile,//szvalidatefile,
                  nholdoutcellrequest,
                  nbagrequest,
                  nmaxknn,
                  szpioneermark,
                  busesamecellfeatures,
                  buseorderfeatures,
                  bdnamethyl,
                  nminnumlocations,
                  nmintotalensemble,
                  numbags,
                  szmethylheader);

            } catch (IllegalArgumentException ex) {
              System.out.println("ERROR\t" + ex.getMessage());
              System.exit(1);
            }
          }
        } catch (Exception ex) {
          ex.printStackTrace(System.out);
          System.exit(1);
        }
      }
    } else {
      System.out.println(
          "Need to specify the mode Convert|ComputeGlobalDist|ExportToChromHMM|GenerateTrainData|Train|Apply|Eval|Version ");
      System.exit(1);
    }

    System.exit(0);
  }
}
