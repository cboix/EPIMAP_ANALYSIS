#!/usr/bin/R
# ---------------------------------------------------------
# Runner to calculate the permuted pvals + tree enrichments
# To submit: 
#   - For shuffled catalogs:
# qsub -cwd -t 1-10000 -P compbio_lab -l h_vmem=12G -l h_rt=06:00:00 -N gtchunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "conda activate mv_env; R --slave -f $BINDIR/calculate_permuted_gwastree_enrichments.R"
# qsub -cwd -t 110-220 -P compbio_lab -l h_vmem=15G -l h_rt=01:00:00 -N gtchunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env; R --slave -f $BINDIR/calculate_permuted_gwastree_enrichments.R"
#   - For permuted by snp:
# qsub -cwd -t 1-243 -P compbio_lab -l h_vmem=25G -l h_rt=36:00:00 -N gtchunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env; R --slave -f $BINDIR/calculate_permuted_gwastree_enrichments.R"
# ----------------------------------------------------
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
# usetree = 'roadmap'  # For old epi only
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only) }
source(paste0(bindir, 'load_gwastree_analysis.R'))

rm(dflist)

print(paste("Images go to: ", treeimgpref))

# --------------------------------------
# Make sub-directories for permuted data
# --------------------------------------
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, regdir, perdir)
system(cmd)

# ----------------------------------
# Determine which uids we will plot:
# ----------------------------------
task = as.integer(as.numeric(Sys.getenv("SGE_TASK_ID")))
catalog = task %/% 10 + 1
chunk = task %% 10
if (chunk == 0){catalog = catalog - 1; chunk = 10}

# Get the nhits: 
load(gwrdafile)
load(gwintrdafile)

# Logistic regression arguments:
type = 'cons'
against = 'parent'
weighted = FALSE
apref = paste0(type, '_', against)
if (weighted){ 
    weights = sqrt(1 / matmarg[,2])
    apref = paste0(apref, '_weighted')
} else {
    weights = NULL
}

# Load in the permutations:
run.shuffled=TRUE
if (run.shuffled){
    # Run a full shuffled catalog:
    NPERM = 1000
    gwperfile = paste0(perdir, sub(".txt", paste0("_", NPERM, sprintf("_%05d_shuffled_catalog.Rda", catalog)), gwcatfile))
    ipref = paste0(apref, sprintf("_permuted_%05d_shuffled", catalog))
    print(paste(catalog, '----- output with prefix', ipref))
    # Put regressions file in the permuted catalog directory:
    load(gwperfile)  # prdf, prgr, qdf, qnumdf, qallnumdf
    pmatfile = paste0(perpref, ipref, '_', NPERM, '_', chunk, '_lreg_pmat.Rda')
    uidlist = sort(unique(qnumdf$uid))
    print(paste(chunk, '----- chunking to file', pmatfile))
    # Chunklength:
    chunklen = length(uidlist) %/% 10 + 1
    chunkid = (chunklen * (chunk - 1) + 1):min(chunklen * chunk, length(uidlist))
    sub.uidlist = uidlist[chunkid]
    NCAT = length(sub.uidlist) # Here NCAT is number GWAS in a catalog
} else {
    # Separate by number of SNPs, not intersections:
    nsnpdf = aggregate(pValue ~ uid, gwdf, length)
    nsnpdf = aggregate(uid ~ pValue, nsnpdf, length)
    head(nsnpdf)
    NSNPS = nsnpdf$pValue[catalog]
    # Get permutations files:
    NCAT = 1000 # Here NCAT is num catalogs
    gwperfile = paste0(perdir, sub(".txt", paste0("_", NCAT, sprintf("_%06d.Rda", NSNPS)), gwcatfile))
    ipref = paste0(apref, sprintf("_permuted_%06d", NSNPS))
    print(paste(NSNPS, '----- output with prefix', ipref))
    load(gwperfile)  # prdf, prgr, qdf, qnumdf, qallnumdf
    # Put regressions file in the permuted catalog directory:
    pmatfile = paste0(perpref, ipref, '_', NCAT, '_', chunk, '_lreg_pmat.Rda')
    uidlist = paste0('X', 1:NCAT)
    # For chunking:
    chunklen = length(uidlist) %/% 10 + 1
    chunkid = (chunklen * (chunk - 1) + 1):min(chunklen * chunk, length(uidlist))
    sub.uidlist = uidlist[chunkid]
    NCAT = length(sub.uidlist) # Here NCAT is number GWAS in a catalog
}


# Load the permuted values:
pmat = matrix(1, nrow=NCAT, ncol=NN)
if (file.exists(pmatfile)){
    load(pmatfile)
}


# Run NCAT (1000 or number of UIDs) regressions:
for (i in 1:NCAT){
    # Only run if no non-1 values:
    if (mean(pmat[i,]) == 1){
        print(i)
        suid = sub.uidlist[i]
        # Run regression:
        ll = try(get_disjoint_lr(trait=suid, type=type, against=against, qdf=qdf, weights=weights))
        if (class(ll) == 'try-error'){
            ll = 'error'
            # Here, p-values kept at 1 across the board.
        } else {
            # Here, update pmat with p-values from null:
            pmat[i,] = ll$df$pout
        }
        if (i %% 10 == 0){
            save(pmat, file=pmatfile)
        }
    }
}

print(paste("[STATUS] Finished calculating permuted values for run:", chunk))
