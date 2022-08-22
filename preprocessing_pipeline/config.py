import os
# ===========================================================
# Main Directories:
domain = os.popen('hostname -d').read().split()[0]
HOME = os.environ.get("HOME")
TMP = os.environ.get("TMP")
if domain == 'broadinstitute.org':
    DBDIR = "/broad/compbio/cboix/EPIMAP_ANALYSIS/db"
    BINDIR = "/broad/compbio/cboix/EPIMAP_ANALYSIS/bin"
    SFTDIR = "/broad/compbio/cboix/software/bin"
else:
    DBDIR = HOME+"/EPIMAP_ANALYSIS/db"
    SFTDIR = HOME+"/bin"
    BINDIR = HOME + "/EPIMAP_ANALYSIS/bin"
ANNDIR = DBDIR + "/Annotation"
ENCODE_PASS = "R7GOMSZL:roctnxvjxljlfrwg"

# Software:
SAMTL = "/broad/software/free/Linux/redhat_6_x86_64/pkgs/" + \
    "samtools/samtools_1.3.1/bin/samtools"
SPPDIR = "/broad/compbio/cboix/software/phantompeakqualtools-master"
HOTSPOTDIR = "/broad/compbio/cboix/software/hotspot/hotspot-distr"
MARKDUP = "/seq/software/picard/1.802/bin/MarkDuplicates.jar"
CHROMIMPUTE = "/broad/compbio/cboix/software/ChromImpute/ChromImpute.jar"
CHMM = "/broad/compbio/cboix/software/ChromHMM/ChromHMM.jar"
MAPDIR = ANNDIR+"/umap"
CIDIR = DBDIR+"/ChromImpute"
CONVERTED_DATADIR = CIDIR+"/converted"
CHMMDIR = DBDIR+"/ChromHMM"
BINARIZED_DIR = CHMMDIR+"/binarized"
# MACS2.1 should be in path (@ $SFTDIR/bin)
# if [[ -z $(which macs2) ]]; then echo 'ERROR: macs2 not in $PATH'; exit 1; fi

# Data directories:
VLDIR = DBDIR+"/Validation"
SAMPDIR = DBDIR+"/Validation/Sampling"
DHSDIR = DBDIR+"/DNase-seq"
DIRS = [DBDIR, BINDIR, ANNDIR, MAPDIR, CIDIR, CONVERTED_DATADIR,
        VLDIR, SAMPDIR, DHSDIR, CHMMDIR, BINARIZED_DIR]

# WCE mapping:
WCE_MAPPING = ANNDIR+"/control_mapping_20181023.tsv"

# Genome Annotation Files:
CHROMSIZES_noY = ANNDIR+"/hg19.chrom.sizes_noY"
CHROMSIZES = ANNDIR+"/hg19.chrom.sizes"
CHROMSIZES_hg38 = ANNDIR+"/hg38.chrom.sizes"
CHROMARM = ANNDIR+"/chromArm.bed"
LOCHAIN = ANNDIR+"/hg38ToHg19.over.chain.gz"
LOCHAIN_REVERSE = ANNDIR+"/hg19ToHg38.over.chain.gz"
GTFDIR = ANNDIR+"/GENCODE"
# GTF_SUFFIX = "v25.primary_assembly"  # For hg38 only
GTF_SUFFIX = "v27lift37.primary_assembly"  # For hg19
GENCODE = GTFDIR+"/gencode."+GTF_SUFFIX+".annotation.gtf.gz"
GENE_GTF = GTFDIR+"/Gene."+GTF_SUFFIX+".bed"
GENE_COLS = GTFDIR+"/Gene."+GTF_SUFFIX+".cols.tsv"

# Mappability regions
MAP36 = MAPDIR+"/k36.Umap.MultiTrackMappability.filtered.bed"
MAP50 = MAPDIR+"/k50.Umap.MultiTrackMappability.filtered.bed"
MAP100 = MAPDIR+"/k100.Umap.MultiTrackMappability.filtered.bed"
BLACKLIST = ANNDIR+"/hg19_blacklist_ENCFF000KJP.bed"

# Variables:
GENOMESIZE = 'hs'  # for MACS2.1
CHRS_noY = ['chr' + str(i+1) for i in range(22)] + ['chrX']
CHRS = ['chr' + str(i+1) for i in range(22)] + ['chrX', 'chrY']
WINDOWS = ["CONTROL2_5k", "2_5k", "5k", "50k", "100k", "150k", "200k"]

# Conda - activate and set temporary pythonpath with just the miniconda libs:
MACS_ENV = 'set +u; source activate macs2_env; \
    PYTHONPATH=/broad/compbio/cboix/software/' + \
    'miniconda2/lib/python2.7/site-packages; set -u'
# ===========================================================


# Helper functions to make settings for processing:
def get_extension(link):
    splt = link.split('.')
    linkext = splt[-2] + '.' + splt[-1]
    fullext = splt[-1]
    if linkext == 'tagAlign.gz':
        return('tagAlign')
    elif fullext == 'bam':
        return('bam')
    else:
        raise AssertionError("Unexpected value of extension 'LINKEXT'", linkext)


# TODO: Sometimes MAPQ threshold fails: look Phred vals
def get_mapq_thresh(proc_status):
    if proc_status == "alignments":
        # print("""Setting no MAPQ threshold for alignments because they have
        #       already been filtered. Some filtered files are Phred64, which
        #       causes problems with using the consistent MAPQ=30.""")
        return(0)
    else:
        return(30)


# Read fraglenfile for fragment length:
def get_fraglen(fraglenfile, replace=73):
    try:
        with open(fraglenfile, 'r') as f:
            for line in f:
                line = line.split("\t")
                fraglen = round(int(line[2]) / 2)
                if fraglen < 10:
                    fraglen = replace
    except:
        fraglen = replace
    return(fraglen)


def get_sval(svalfile):
    with open(svalfile, 'r') as f:
        for line in f:
            sval = line.split()[0]
    return(sval)


def check_mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def read_isPaired(input):
    with open(input, 'r') as f:
        for line in f:
            line = int(line.split("\n")[0])
    return(line)
