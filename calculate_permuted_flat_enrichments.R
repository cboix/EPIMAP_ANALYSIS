#!/usr/bin/R
# ---------------------------------------------------------
# Runner to calculate the permuted pvals + FLAT enrichments
# To submit: 
# qsub -cwd -t 1-243 -P compbio_lab -l h_vmem=25G -l h_rt=3:00:00 -N flat_gtchunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env; R --slave -f $BINDIR/calculate_permuted_flat_enrichments.R"
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
chunk = as.integer(as.numeric(Sys.getenv("SGE_TASK_ID")))
# Get the nhits: 
load(gwrdafile)
load(gwintrdafile)
# nhitdf = aggregate(subjectHits ~ uid, qdf, length)
# nhitdf = aggregate(uid ~ subjectHits, nhitdf, length)

# Separate by number of SNPs, not intersections:
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
nsnpdf = aggregate(uid ~ pValue, nsnpdf, length)
head(nsnpdf)
NSNPS = nsnpdf$pValue[chunk]

# ---------------------------
# Run permutations for NSNPS:
# ---------------------------
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
NCAT = 1000

for (flatset in c('epigenomes', 'modules')){
    ipref = paste0(apref, '_', flatset, sprintf("_permuted_%06d", NSNPS))
    print(paste(NSNPS, '----- output with prefix', ipref))
    gwperfile = paste0(perdir, sub(".txt", paste0("_", NCAT, sprintf("_%06d.Rda", NSNPS)), gwcatfile))
    # Put regressions file in the permuted catalog directory:
    pmatfile = paste0(perpref, ipref, '_', NCAT, '_lreg_pmat.Rda')

    # Load in the enhancer sets:
    enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
    if (!file.exists(enhsetfile)){
        if (flatset == 'epigenomes'){
            mat = read.delim(gzfile(enhmatfile), sep="\t", header=T)
            uqepi = sort(unique(mat$col))
            enhsets = sapply(uqepi, function(x){ mat$row[mat$col == x] })
            lenlist = unlist(lapply(enhsets, length))
            save(enhsets, lenlist, file=enhsetfile)
            rm(mat)
        } else if (flatset == 'modules'){
            moddf = read.delim(file='cls_merge2_wH3K27ac100_300_assignments_ENH_052819.bed', sep="\t", header=T)
            uqcls = sort(unique(moddf$cls))
            enhsets = sapply(uqcls, function(x){ moddf$id[moddf$cls == x] })
            lenlist = unlist(lapply(enhsets, length))
            save(enhsets, lenlist, file=enhsetfile)
        }
    } else {
        load(enhsetfile)
        NF = length(enhsets)
    }

    # Load the permuted values:
    load(gwperfile)  # prdf, prgr, qdf, qnumdf, qallnumdf
    # Evaluate the intersections (pre-hypergeom test)
    dflist = eval.intersections(enhsets, qdf)
    pmat = matrix(1, nrow=NCAT, ncol=NF)
    if (file.exists(pmatfile)){
        load(pmatfile)
    }

    # Run NCAT hgtest:
    for (i in 1:NCAT){
        if (mean(pmat[i,]) == 1){
            print(i)
            suid = paste0('X', i)
            # Run regression:
            suid = paste0('X', i)
            ll = try(get_pvalues_flat(suid, dflist, lenlist, qdf, NF=NF, cutp=CUTP))
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
}

print(paste("[STATUS] Finished calculating permuted values for run:", chunk))
