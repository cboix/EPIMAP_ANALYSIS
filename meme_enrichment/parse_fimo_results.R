#!/usr/bin/env Rscript
# ------------------------------------------
# Parse raw FIMO outputs
# Calculates pvals using Fisher's exact test
# ------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
library(caTools)
library(gtools)
options(scipen=20)
options(warn=2)

cpref='cls_merge2_wH3K27ac100_raw'
# cpref='cls_merge2_wH3K27ac100_300'
calldir = '/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls'
infofile = paste0(calldir, '/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/', cpref, '/enrichment_infofile.tsv')
resultsdir = paste0(calldir, '/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/', cpref, '/motif_enrichment')

args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied.")
} else {
    infofile = args[1]
    resultsdir = args[2]
}

# NOTE: 0-indexed ids:
chunk = as.integer(as.numeric(Sys.getenv("SGE_TASK_ID"))) - 1

# Directories:
outputdir <- paste(resultsdir, "parsed_fimo_results", sep="/");
dir.create(outputdir, showWarnings=FALSE, recursive=TRUE);

# Read info table:
info = read.delim(infofile, header=F, stringsAsFactors=F) 
names(info) = c('file', 'name')

# Read aggregated files: 
print("[STATUS] Loading alldf")
alldf = read.delim(gzfile(paste0(outputdir, '/all_agg_fimo.tsv.gz')), sep="\t")
print("[STATUS] Loading q10df")
q10df = read.delim(gzfile(paste0(outputdir, '/q10_agg_fimo.tsv.gz')), sep="\t")
print("[STATUS] Parsed results loaded")

# Perform fisher's exact to get pval:
run.fisher <- function(y){
    table <-  matrix(as.numeric(c(y[1], y[2], y[3], y[4])), ncol = 2, byrow = TRUE)
    if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
    return(p)
}

# NOTE: Split p-value calculation:
cols = c('id','motif','db','pval','log2fc')
subdf = alldf[alldf$id == chunk,]
subdf = subdf[subdf$total - subdf$counts >= 0,]
fdf = data.frame(c1 = subdf$counts, n1 = subdf$total - subdf$counts, 
                 c2 = subdf$c.tot, n2 = subdf$t.tot - subdf$c.tot)

print(paste("[STATUS] Computing pvals for", nrow(fdf), "rows in alldf"))
pvals = apply(fdf, 1, run.fisher)
head(sort(pvals), 40)
subdf$pval = pvals
write.table(subdf[,cols], paste0(outputdir, '/all_parsed_pval_',chunk,'_fimo.tsv'), sep="\t", quote=F, row.names=F)

# For q10 matrix
subdf = q10df[q10df$id == chunk,]
subdf = subdf[subdf$total - subdf$counts >= 0,]
fdf = data.frame(c1 = subdf$counts, n1 = subdf$total - subdf$counts, 
                 c2 = subdf$c.tot, n2 = subdf$t.tot - subdf$c.tot)

print(paste("[STATUS] Computing pvals for", nrow(fdf), "rows in q10df"))
pvals = apply(fdf, 1, run.fisher)
head(sort(pvals), 40)
subdf$pval = pvals
write.table(subdf[,cols], paste0(outputdir, '/q10_parsed_pval_',chunk, '_fimo.tsv'), sep="\t", quote=F, row.names=F)

print(paste("[STATUS] Finished all pvalue calculations for chunknum:", chunk))
