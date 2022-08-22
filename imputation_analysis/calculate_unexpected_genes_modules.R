#!/usr/bin/R
# -----------------------------------------------------
# Look at the number of unexpected genes in the modules
# -----------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))

# Arguments for loading data:
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
runall = TRUE # Internal argument for re-running calc.

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only) }
source(paste0(bindir, 'load_gwastree_analysis.R'))
rm(dflist)


# Load in the enhancer sets:
flatset = 'modules'
enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
load(enhsetfile)
NF = length(enhsets)

# Find nearest genes for each enhancer, from the full set of protein coding
gdf = read.delim('Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz', sep="\t", header=F)
names(gdf) = c('chr', 'tss','tss2','name')
gdgr = GRanges(gdf$chr, IRanges(gdf$tss, gdf$tss))
ndf = data.frame(queryHits=1:length(dmgr), subjectHits=nearest(dmgr, gdgr))
ndf$name = as.character(gdf[ndf$subjectHits, 'name'])

gmap = read.delim('Annotation/gencode.gene.name.mapping.tsv', sep="\t", header=F)
names(gmap) = c('name','symbol')
# Annotate genes as known/unknown function:
gmap$unann = 0
gmap$unann[grep("AC[0-9]*.[0-9]", gmap$symbol)] = 1
gann = aggregate(unann ~ name, gmap, max)
rownames(gann) = gann$name
ndf$unann = gann[ndf$name, 'unann']

# Enhancers closest to unannotated.
# TODO: If not in GO annotations.

# Perform Fisher's tests for all modules:
ftests = sapply(1:NF, function(i){ 
                    set = enhmap[enhsets[[i]]]
                    ua = ndf$unann[set]
                    mat = matrix(c(nrow(ndf), length(set),
                                   sum(ndf$unann), sum(ua)), ncol=2, byrow=TRUE)
                    fisher.test(mat, alternative='greater')$p.value })
                    # fisher.test(mat, alternative='less')$p.value })

barplot(-log10(ftests))
ftests[17]

# Specific modules of interest:
mod = c('c117', 'c90','c43','c31','c272','c122','c173',
        'c16', 'c249',
        'c161','c166','c63','c89',
        'c267','c12')

names(ftests) =paste0('c', 1:300 - 1)
round(-log10(ftests[mod]),1)

round(sort(-log10(ftests)),2)
# c27  c156   c95   c10   c28  c190  c169  c168   c61   c51  c278  c200  c133 c9
# 4.53  4.55  4.57  4.64  5.32  5.80  6.36  7.56  7.91  8.12 10.43 13.92 23.82 50.48

# 9, 133, 200, 278 are
# LCL, LCL, cancer, cancer







