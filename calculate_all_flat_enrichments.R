#!/usr/bin/R
# ------------------------------------------------
# Runner to calculate all the flat enrichments (after perm.)
# To submit:
# qsub -cwd -t 1-55 -P compbio_lab -l h_vmem=25G -l h_rt=3:00:00 -hold_jid flat_gtchunked -N flat_calcgtchunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env; R --slave -f $BINDIR/calculate_all_flat_enrichments.R"
# ------------------------------------------------
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

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "gwas_tree_analysis/examples/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)

# Under dbdir:
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
CHUNKSIZE = 100
chunk = as.integer(as.numeric(Sys.getenv("SGE_TASK_ID")))
uids = sort(as.character(unique(gwdf$uid)))
plotrange = ((chunk - 1)* CHUNKSIZE + 1):min(chunk * CHUNKSIZE, length(uids))
plotuids = uids[plotrange]
print("[STATUS] Will calculate the following uids:") 
print(plotuids)

# --------------------------
# Run analysis for each uid:
# --------------------------
# Arguments:
NCAT=1000
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

# Counts of snps:
nsnpdf = aggregate(pValue ~ uid, gwdf, length)

for (flatset in c('epigenomes', 'modules')){
    # Load in the enhancer sets:
    enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
    load(enhsetfile)
    NF = length(enhsets)
    # Evaluate the intersections (pre-hypergeom test)
    dflist = eval.intersections(enhsets, qdf)

    for (plotind in plotrange){
        suid = uids[plotind]
        ipref = paste0(apref, '_', flatset, sprintf("_%05d", plotind))
        print(paste(plotind, '-', suid, '----- output with prefix', ipref))

        # Perform or load hg.test:
        regfile = paste0(regpref, ipref, '_hg.Rda')
        if (!file.exists(regfile) || runall){
            ll = try(get_pvalues_flat(suid, dflist, lenlist, 
                                      qdf, NF=NF, cutp=CUTP))
            if (class(ll) == 'try-error'){ ll = 'error' }
            save(ll, file=regfile)
        } else {
            load(regfile)
        }
        ll.fixed = ll

        # -----------------------
        # Correct the regression:
        # -----------------------
        # Apply appropriate correction per node or set:
        nsnps = nsnpdf$pValue[nsnpdf$uid == suid]
        pmatpref = paste0(apref,'_', flatset , sprintf("_permuted_%06d", nsnps))
        pmatfile = paste0(perpref, pmatpref, '_', NCAT, '_lreg_pmat.Rda')
        load(pmatfile)
        # Adjust at diff levels:
        for (lvl in c(1,5,10)){
            if (class(ll) == 'list'){
                cregfile = paste0(regpref, ipref, '_hg_adj', NCAT, '_', lvl, '.Rda')
                if (!file.exists(cregfile) || runall){
                    # Load in the null p-values matrix:
                    ll = fdr.adjust(ll.fixed, pmat, NPERM=NCAT, NBELOW=lvl)
                    save(ll, file=cregfile)
                }
            }
        }
    }
}

print(paste("[STATUS] Finished calculating all uids for chunk", chunk))
