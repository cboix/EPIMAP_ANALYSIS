#!/usr/bin/R
# ---------------------------------------
# Load the distance matrices:
# Imputed, Observed, and Mixed
# Additionally, load in annotation
# ---------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
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

# -----------------------
# Load datasets/matrices:
# -----------------------
# Locations/files:
fnames <- list.files(path='ChromImpute',pattern='imp_distance_all_*')
marks <- c(sub('.tsv','',sub('imp_distance_all_','',fnames)),'Full')
ofnames <- list.files(path='ChromImpute',pattern='^distance_all_*') # Observed
omarks <- c(sub('.tsv','',sub('distance_all_','',ofnames)),'Full')
keep.id = which(omarks %in% marks)
omarks = omarks[keep.id]
ofnames = ofnames[keep.id]
NMARKS <- length(marks) - 1

# Imputed + observed distance matrices:
ll <- list()
obsll <- list()
full.ll <- list()
ct.list <- c()
obsct.list = c()
fullct.list = c()
for (i in 1:NMARKS){
    mark <- marks[i]
    mat <- read.delim(paste0('ChromImpute/',fnames[i]),sep="\t",header=T)
    if(nrow(mat) > 1){
        print(paste0('Evaluating ',mark))
        mat = process.dist.mat(mat)
        nam = colnames(mat)
        mdf = t(sapply(nam, function(x){strsplit(x, "_")[[1]]}))
        iset = which(mdf[,2] == 'imp')
        oset = which(mdf[,2] == 'obs')
        # Split up matrices:
        imat = mat[iset, iset]
        omat = mat[oset, oset]
        colnames(imat) = mdf[iset, 1]
        rownames(imat) = mdf[iset, 1]
        colnames(omat) = mdf[oset, 1]
        rownames(omat) = mdf[oset, 1]
        # Add each cell list
        ct.list <- sort(unique(c(ct.list,colnames(imat))))
        obsct.list <- sort(unique(c(obsct.list,colnames(omat))))
        fullct.list<- sort(unique(c(fullct.list,colnames(mat))))
        # Add matrices:
        ll[[mark]] <- imat
        obsll[[mark]] <- omat
        full.ll[[mark]] <- mat
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
        tmp[is.na(tmp)] = 0
        full = full + tmp[ct.list, ct.list]
        tmp = tmp * NA
    }
}
full = full / idm
