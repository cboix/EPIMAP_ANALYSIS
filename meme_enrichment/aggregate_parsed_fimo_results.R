#!/usr/bin/env Rscript
# ----------------------------------
# Aggregate parsed fimo
# Correct pvals + diagnostics
# Plot reordered enrichment matrices
# ----------------------------------
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

cpref='cls_merge2_wH3K27ac100_300'
# cpref='cls_merge2_wH3K27ac100_raw'
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

# Directories:
outputdir <- paste(resultsdir, "parsed_fimo_results", sep="/");
dir.create(outputdir, showWarnings=FALSE, recursive=TRUE);

info = read.delim(infofile, header=F, stringsAsFactors=F) 
names(info) = c('file', 'name')

# Make all plots, first with alldf and then with q10df:
tags = c('All Matches','Matches with FDR < 0.1')
prefs = c('all','q10')

for (i in 1:2){
    tagline = tags[i]
    prefix = prefs[i]
    plotpref = paste0(outputdir, '/', prefix, '_')
    print(paste("[STATUS] Plotting under prefix", plotpref))

    # Aggregate all files:
    alldf = c()
    for (id in info$name){
        fname = paste0(outputdir, '/',prefix,'_parsed_pval_',id,'_fimo.tsv')
        if (file.exists(fname)){
            adf = try(read.delim(fname, sep="\t"))
            cols = c('id','motif','db','pval','log2fc')
            if (class(adf) != 'try-error'){
                alldf = rbind(alldf, adf[,cols])
            }
        } else { print(fname) }
    }

    alldf$pval = as.numeric(as.character(alldf$pval))

    # Plot pval histogram:
    png(paste0(plotpref, 'pval_hist.png'), units='in', res=250, width=7, height=5)
    hist(alldf$pval, 50, border='white', col='darkgrey')
    dev.off()

    # Calculate expectations (number of tests under uniform):
    # Plot qqplot - pvals vs. expected:
    # NOTE: takes a while to plt points - find better way to do this?
    # exp.pvals = (rank(alldf$pval, ties.method="first") + .5) / (length(alldf$pval) + 1)
    # png(paste0(plotpref, 'pval_qqplot.png'), units='in', res=250, width=6, height=5)
    # plot(-log10(exp.pvals), -log10(alldf$pval), asp=1)
    # abline(0, 1)
    # dev.off()

    # Adjust pvalues to q-value, get log10qval
    alldf$padj = p.adjust(alldf$pval)
    png(paste0(plotpref, 'pval_adj_hist.png'), units='in', res=250, width=7, height=5)
    hist(alldf$padj, 50, border='white', col='darkgrey')
    dev.off()

    # Add log10p and log10padj (more compressed):
    alldf$log10p = -log10(alldf$pval)
    alldf$log10pa = -log10(alldf$padj)

    # Theshold:
    print("Thresholding at p.adj < 0.05")
    sigdf = alldf[alldf$padj < 0.05,]
    print(paste(length(unique(alldf$id)), 'to', length(unique(sigdf$id)), 'ids'))
    print(paste(length(unique(alldf$motif)), 'to', length(unique(sigdf$motif)), 'motifs'))

    # TODO: Add reduced plot here:
    # Plot matrix of log2FC and of log10qval
    # Plot the fimo output matrices alone

    # --------------------------
    # Save as RDA and as tsv.gz:
    # --------------------------
    save(sigdf, file=paste(outputdir, "/", prefix, "_sig_merged_fimo_pval.RData", sep=""))
    write.table(sigdf, gzfile(paste(outputdir, "/", prefix, "_sig_merged_fimo_pval.tsv.gz", sep="")), sep="\t", quote=F, row.names=F)

    save(alldf, file=paste(outputdir, "/", prefix, "_merged_fimo_pval.RData", sep=""))
    write.table(alldf, gzfile(paste(outputdir, "/", prefix, "_merged_fimo_pval.tsv.gz", sep="")), sep="\t", quote=F, row.names=F)
}

