#!/usr/bin/R
# ---------------------------------------
# Load the distance matrices from fdist:
# Imputed, Observed, and Mixed
# Additionally, load in annotation
# ---------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
library(circlize)
library(ComplexHeatmap)

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    nregions = 20000
    dataset = 'cor'
    print(paste0("No arguments supplied. Defaulting to loading nregions = ", 
                 nregions, ' and dataset ', dataset))
} else {        
    nregions = args[1]
    dataset = args[2]
}

# Set the metric printable name:
if (dataset == 'spearman'){
    metric = 'Spearman Rho'
} else if (dataset == 'normdist') {
    metric = 'Normalized Euclidean Distance'
} else if (dataset == 'jaccard') {
    metric = 'Jaccard Distance'
} else if (dataset == 'cor'){
    metric = 'Correlation'
} else if (dataset == 'rtd'){
    metric = 'Rogers-Tanimoto Distance'
}

# -----------------------
# Load datasets/matrices:
# -----------------------
# Locations/files:
fnames <- list.files(path='ChromImpute/region_distance/',pattern=paste0('all_fixeddist_',nregions,'_*'))
marks <- c(sub('.tsv.gz','',sub(paste0('all_fixeddist_',nregions,'_'),'',fnames)),'Full')
NMARKS <- length(marks) - 1

# For getting names from file prefix
rownames(tracktab) = tracktab$prefix

# For split obsimp:
splitmat = function(mat, set, mdf){
    mat = mat[set, set]
    colnames(mat) = mdf[set, 1]
    rownames(mat) = mdf[set, 1]
    return(mat)
}

# Imputed + observed distance matrices:
# Read in long format:
ll <- list()
obsll <- list()
full.ll <- list()
ct.list <- c()
obsct.list = c()
fullct.list = c()
for (i in 1:NMARKS){
    mark <- marks[i]
    df <- read.delim(gzfile(paste0('ChromImpute/region_distance/',fnames[i])),sep="\t",header=F)
    if(nrow(df) > 1){
        print(paste0('Evaluating ',mark))
        names(df) = c('prefix','against','nregions','fixed','cor','normdist', 'spearman','jaccard','rtd')

        cwide = spread(df[,c('prefix','against',dataset)], against, dataset)
        cmat = as.matrix(cwide[,-1])
        colnames(cmat) = tracktab[colnames(cmat), "uqsample"]
        rownames(cmat) = tracktab[as.character(cwide$prefix), "uqsample"]
        rn = unique(sort(rownames(cmat)))
        cmat = cmat[rn,rn]
        if (dataset %in% c('cor','spearman')){
            cmat = 1 - cmat 
        }

        # Split up matrices:
        nam = colnames(cmat)
        mdf = t(sapply(nam, function(x){strsplit(x, "_")[[1]]}))
        iset = which(mdf[,2] == 'imp')
        oset = which(mdf[,2] == 'obs')
        cimat = splitmat(cmat, iset, mdf)
        comat = splitmat(cmat, oset, mdf)

        # Add each cell list
        ct.list <- sort(unique(c(ct.list, colnames(cimat))))
        obsct.list <- sort(unique(c(obsct.list, colnames(comat))))
        fullct.list<- sort(unique(c(fullct.list, colnames(cmat))))
        
        # Add matrices - correlation as 1 - cor
        ll[[mark]] <- cimat
        obsll[[mark]] <- comat
        full.ll[[mark]] <- cmat
    }
}


# Remove cells with no ATAC-seq imputed (bad imputation):
atacimp = as.character(rownames(ll[['ATAC-seq']]))
ct.list = ct.list[ct.list %in% atacimp]
ct.list = ct.list[ct.list %in% cellorder]

# Fused distance matrix by equal weight of all:
N <- length(ct.list)
full  <- matrix(0, nrow=N, ncol=N, dimnames=list(ct.list,ct.list))
idm = full
tmp = full * NA
mainmarks = c('H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9me3')
NMAIN = length(mainmarks)
for (i in 1:NMAIN){
    mark <- mainmarks[i]
    mat <- ll[[mark]]
    if (!is.null(mat)){
        rn = as.character(colnames(mat))
        rn = rn[rn %in% ct.list]
        tmp[rn,rn] = mat[rn,rn]
        idm = idm + 1.0 * (!is.na(tmp[ct.list, ct.list]))
        # Weighted sum:
        # full = full + tmp[ct.list, ct.list] / NMARKS
        tmp[is.na(tmp)] = 0
        full = full + tmp[ct.list, ct.list]
        tmp = tmp * NA
    }
}
full = full / idm
