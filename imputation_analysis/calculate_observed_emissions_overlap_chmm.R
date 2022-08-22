#!/usr/bin/R
# -----------------------------------------------------
# qsub -cwd -P compbio_lab -l h_vmem=5G -t 1-833 -l h_rt=2:00:00 -N calc_obs_emission_ovl -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env > /dev/null; R --slave -f $BINDIR/calculate_observed_emissions_overlap_chmm.R"
# Compare the observed + states for a single epigenome:
# -----------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))

chunk = as.integer(as.numeric(Sys.getenv('SGE_TASK_ID')))
sample = cellorder[chunk]

# Load in the 
calldir='ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/'
stdir=paste0(calldir, 'STATEBYLINE/')
binpref = paste0('ChromHMM/files_by_cell/',sample, '/', 
                 sample, '_observed_aux_18_on_mixed_impobs_QCUT_')

# Set up vars, states x marks:
ncall = rep(0,18)
nmark = rep(0,6)
ovmat = matrix(0, nrow=18, ncol=6)
colnames(ovmat) = c('H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3')
 
# Iterate through chromosomes, get state + observed overlap:
chrlist = paste0('chr', c(1:22, 'X'))
for (chr in chrlist){
    print(chr)
    # Read in the chromatin state calls:
    cfile = paste0(stdir, sample, '_18_CALLS_PER_LINE_', chr,'_statebyline.txt.gz')
    x = as.numeric(scan(gzfile(cfile), 'c', sep='\r', skip=2))
    # Read in marks:
    bfile = paste0(binpref, chr, '_binary.txt.gz')
    y = read.delim(gzfile(bfile), header=T, sep="\t", skip=1)
    y = as.matrix(y)
    # Calculate the # per mark occurrences:
    tform = make.tform(x, u=1:18)
    ovmat = ovmat + t(tform) %*% y
    ncall = ncall + apply(tform, 2, sum)
    nmark = nmark + apply(y, 2, sum)
} 

norm.mat = sweep(ovmat, 1, ncall, '/')
normdf = data.frame(norm.mat)
normdf$state = 1:18
normdf = gather(normdf, mark, percent, -state)
normdf$sample = sample

write.table(normdf, file=paste0(calldir, sample, '_observed_emissions_overlap.tsv'), 
            row.names=F, quote=F, sep="\t")
