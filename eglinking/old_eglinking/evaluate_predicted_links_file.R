#!/usr/bin/R
# -------------------------------------------
# Given a file, evaluate the predicted links:
# -------------------------------------------
# 1. Evaluate with: 
# - Best matched GTEx eQTLs
# - Best matched HiC loops
# - GWAS SNP capture (enrichment)

# 2. Compare against:
# - Nearest enh
# - Nearest N (5) enh
# - Nearest all enh penalized by distance
# - Roadmap links
# - ABC model

# For now, evaluate with AUPRC (F1 score) or some calc of sig. enrichment statistics

# --------------------------
gtargs=(commandArgs(TRUE))
print(gtargs)
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
    plotting.trees=TRUE
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
    plotting.trees=TRUE
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))

library(dplyr)
library(cba)
library(GenomicRanges)
library(Matrix)
options(scipen=45) # So we dont get issues writing integers into bedfiles
options(repos='http://cran.rstudio.com/')

# NOTE: Using full dataset needs more than 120G of space, may not even work in R

filename='BSS00001_t1.75.tsv'
filename='BSS00007_t1.75.tsv'
if (length(gtargs) > 0){
    filename = as.character(gtargs[1])
} else {
    print("[WARNING] Need a links file")
}
print(paste("Using", filename))


# -----------------
# Make directories:
# -----------------
lddir = "linking_data/"
evdir = "linking_data/eval/"
lkdir = "linking_data/links/"
imgdir = paste0(img, "links_eval/")
cmd = paste('mkdir -p ', imgdir, lddir, evdir, lkdir)
system(cmd)
basename = sub(paste0(".*/"), "", filename)  # Strip dir for now
setprefix = sub("\\.", "_", sub("\\.tsv","",basename))
imgpref = paste0(imgdir, setprefix, "_")


# -------------------------
# Load the predicted links:
# -------------------------
sampleid = sub("_.*","", basename)
ldf = read.delim(paste0(lkdir, basename), header=F)
names(ldf) = c('chr','start', 'end','gene','score','dist')
ld2 = unique(ldf)
ldf = unique(ldf)




# ---------------------
# Load validation data:
# ---------------------
# - Gene annotation
# - Best matched GTEx eQTLs
# - Best matched HiC loops
# - GWAS SNP capture (enrichment)










