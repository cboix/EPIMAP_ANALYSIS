#!/usr/bin/env Rscript
# --------------------------
# Aggregate raw FIMO outputs
# --------------------------
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

# Directories:
outputdir <- paste(resultsdir, "parsed_fimo_results", sep="/");
dir.create(outputdir, showWarnings=FALSE, recursive=TRUE);

info = read.delim(infofile, header=F, stringsAsFactors=F) 
names(info) = c('file', 'name')

# Aggregate all files:
alldf = c()
q10df = c()
for (id in info$name){
    filepref = paste0(resultsdir, '/regions_', id, '_fimo')
    adf = try(read.delim(gzfile(paste0(filepref, '.counts.tsv.gz')), sep="\t", header=F))
    if (class(adf) != 'try-error'){
        names(adf) = c('motif','counts','total','id','db')
        alldf = rbind(alldf, adf)
    }
    qdf = try(read.delim(gzfile(paste0(filepref, '_q10.counts.tsv.gz')), sep="\t", header=F))
    if (class(qdf) != 'try-error'){
        names(qdf) = c('motif','counts','total','id','db')
        q10df = rbind(q10df, qdf)
    }
}

# NOTE: Merge with incomplete to get 0 count areas:
# All regions tested - against all motifs:
motifs = unique(alldf[,c('motif','db')])
print(paste("There are", dim(motifs)[1], "motifs"))

# All regions
iddf = unique(alldf[,c('total','id')])
print(paste("Loaded", dim(iddf)[1], "sets of regions"))

alltests = merge(motifs, iddf)
print(paste("Total", dim(alltests)[1], "tests"))

# iddf = aggregate(total ~ id, alldf, unique)
allm = merge(alldf, alltests, all.y=TRUE)
q10m = merge(q10df, alltests, all.y=TRUE)
allm$counts[is.na(allm$counts)] = 0
q10m$counts[is.na(q10m$counts)] = 0

# Aggregate all counts - background level of motifs:
alltot = aggregate(cbind(counts, total) ~ motif + db, allm, sum)
q10tot = aggregate(cbind(counts, total) ~ motif + db, q10m, sum)
names(alltot)[3:4] = c('c.tot','t.tot')
names(q10tot)[3:4] = c('c.tot','t.tot')

# Remove all with total > counts #
# Specifically a degenerate HOMOCO MAZ motif 
# NOTE: If aggr from DB, only took top hit in region 
# If multiple hits, do not have imbalance issue
alltot = alltot[alltot$c.tot < alltot$t.tot,]
q10tot = q10tot[q10tot$c.tot < q10tot$t.tot,]

# Merge and get log2FC
alldf = merge(allm, alltot)
q10df = merge(q10m, q10tot)
fc = with(alldf, (counts/total) / (c.tot/t.tot))
qfc = with(q10df, (counts/total) / (c.tot/t.tot))
alldf$log2fc = log2(fc)
q10df$log2fc = log2(qfc)

# Write to file:
write.table(alldf, gzfile(paste0(outputdir, '/all_agg_fimo.tsv.gz')), sep="\t", quote=F, row.names=F)
write.table(q10df, gzfile(paste0(outputdir, '/q10_agg_fimo.tsv.gz')), sep="\t", quote=F, row.names=F)
