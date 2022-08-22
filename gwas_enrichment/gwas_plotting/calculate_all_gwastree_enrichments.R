#!/usr/bin/R
# ----------------------------------------------------
# Runner to plot for tree + gene tree for NCHUNK files
# To submit: 
# qsub -cwd -t 1-556 -P compbio_lab -l h_vmem=25G -l h_rt=2:00:00 -hold_jid gtchunked -N calcgtchunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "source activate mv_env; R --slave -f $BINDIR/calculate_all_gwastree_enrichments.R"
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

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "gwas_tree_analysis/examples/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
treeimgpref = paste0(imgdir, usetree, '_e', tol, '_')

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
CHUNKSIZE = 10
chunk = as.integer(as.numeric(Sys.getenv("SGE_TASK_ID")))
plotrange = ((chunk - 1)* CHUNKSIZE + 1):(chunk * CHUNKSIZE)
uids = sort(as.character(unique(gwdf$uid)))
plotuids = uids[plotrange]
print("[STATUS] Will calculate + plot the following uids:") 
print(plotuids)

# ------------------------------
# Load RNA-seq data for linking:
# ------------------------------
# Read in the RNA-seq (availability):
rnadir = 'RNA-seq/files/RNA-seq/tsv/'
availfile = paste0(rnadir, 'avail_rna.txt')
pcid = scan(availfile, 'c', sep="\n")
rna.avail = 1 * (labels(dend) %in% pcid)

# Load in TSS annotation:
tssfile = 'Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz'
tssdf = read.delim(gzfile(tssfile), header=F, sep="\t", stringsAsFactors=F)
names(tssdf) = c('chr','tss','t2','gene')

gmfile = 'Annotation/gencode.gene.name.mapping.tsv'
gmdf = read.delim(gmfile, header=F, sep="\t", stringsAsFactors=F)
names(gmdf) = c('gene', 'symbol')

# -----------------
# RNA-seq consensus
# -----------------
rnallfile = paste0(rnadir, 'rna_consensus_object_', usetree, '_070819.Rdata')
print("[STATUS] Loading rnall from file")
load(rnallfile)

NCAT=1000
# Plotting function:
plot.treebasic = function(file, ll, dd, suid, lab, gg=NULL, legend=FALSE){
    pdf(file,width=14.5,height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    if (is.null(gg)){
        circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE, hit.track=dd$hitdf, leaf.pval=all.leafmat[suid,])
    } else {
        circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE, hit.track=dd$hitdf, leaf.pval=all.leafmat[suid,], nodedf=gg$nodedf, rna.avail=rna.avail)
    }
    upViewport()
    circos.clear()
    if (legend){
        draw(pd.legend.ext, x = circle_size, just = "left")
    }
    title(paste(suid, '-', ll$title))
    dev.off()
}


# --------------------------
# Run analysis for each uid:
# --------------------------
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

# Counts of snps:
nsnpdf = aggregate(pValue ~ uid, gwdf, length)

for (plotind in plotrange){
    ipref = paste0(apref, sprintf("_%05d", plotind))
    suid = uids[plotind]
    print(paste(suid, '----- output with prefix', ipref))
    nsnps = nsnpdf$pValue[nsnpdf$uid == suid]
    npref = paste0(apref, sprintf("_permuted_%06d", nsnps))

    # Perform or load regression:
    regfile = paste0(regpref, ipref, '_lreg.Rda')
    if (!file.exists(regfile)){
        ll = try(get_disjoint_lr(trait=suid, qdf=qdf, type=type, against=against, weights=weights))
        if (class(ll) == 'try-error'){
            ll = 'error'
        }
        save(ll, file=regfile)
    } else {
        load(regfile)
    }
    ll.fixed = ll

    # -----------------------
    # Correct the regression:
    # -----------------------
    # Adjust at diff levels:
    for (lvl in c(1,5,10)){
        if (ll != 'error'){
            cregfile = paste0(regpref, ipref, '_lreg_adj', NCAT, '_', lvl, '.Rda')
            if (!file.exists(cregfile)){
                # Load in the null p-values matrix:
                pmatfile = paste0(perpref, npref, '_', NCAT, '_lreg_pmat.Rda')
                load(pmatfile)
                ll = fdr.adjust(ll.fixed, pmat, NPERM=NCAT, NBELOW=lvl)
                save(ll, file=cregfile)
            } else {
                load(cregfile)
            }
        }
    }

    # TODO: Plot with and without the correction?
    if (ll != 'error'){
        # Plot if there are any significant hits:
        for (lvl in c(1,5,10)){
            ltail = paste0('_adj', NCAT, '_', lvl, '.Rda')
            cregfile = paste0(regpref, ipref, '_lreg', ltail)
            treefile = paste0(treeimgpref, ipref, '_lrtree', ltail) 
            load(cregfile)
            if (sum(ll$log10p > 0) && (!file.exists(treefile))){
                # Plot without linking genes:
                dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=3)
                plot.treebasic(treefile, ll, dd, suid, lab=lab, gg=NULL)
            }
        }

        plot.extra=FALSE
        if (plot.extra){
            # Link genes:
            gg = try(get.linked.genes(ll, suid, allgenes=FALSE, minp=3))
            if (class(gg) != 'try-error'){
                ggfile = paste0(regpref, ipref, '_linked_genes.tsv')
                if ((!file.exists(ggfile)) && (nrow(gg$cgdf) > 0)){
                    # TODO: Add the nearest gene as well.
                    write.table(gg$cdf, file=ggfile, quote=F, row.names=F, sep="\t")
                }
                # Plot dendrogram - with genes on it
                linkedtreefile = paste0(treeimgpref, ipref, '_lrtree_genes.pdf')
                if (!file.exists(linkedtreefile)){
                    dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=3, altline=0)
                    plot.treebasic(linkedtreefile, ll, dd, suid, lab=lab, gg=gg)
                }
            }

            # Link genes:
            gg = try(get.linked.genes(ll2, suid, allgenes=FALSE, minp=2))
            if (class(gg) != 'try-error'){
                ggfile = paste0(regpref, ipref, '_linked_genes2.tsv')
                if ((!file.exists(ggfile)) && (nrow(gg$cgdf) > 0)){
                    # TODO: Add the nearest gene as well.
                    write.table(gg$cdf, file=ggfile, quote=F, row.names=F, sep="\t")
                }
                # Plot dendrogram - with genes on it
                linkedtreefile = paste0(treeimgpref, ipref, '_lrtree_genes2.pdf')
                if (!file.exists(linkedtreefile)){
                    dd  = setup_dendrogram(dend3, ll2, udf, declist=declist, bcutoff=3, altline=0)
                    plot.treebasic(linkedtreefile, ll2, dd, suid, lab=lab, gg=gg)
                }
            }

            # Add ALL lower threshold genes + their association.
            # NOTE: won't do this 2x.
            gg = try(get.linked.genes(ll2, suid, allgenes=TRUE, minp=2))
            if (class(gg) != 'try-error'){
                ggfile = paste0(regpref, ipref, '_linked_genes_all.tsv')
                if ((!file.exists(ggfile)) && (nrow(gg$cgdf) > 0)){
                    # TODO: Add the nearest gene as well.
                    write.table(gg$cdf, file=ggfile, quote=F, row.names=F, sep="\t")
                }
                # Plot dendrogram - with genes on it
                linkedtreefile = paste0(treeimgpref, ipref, '_lrtree_genes_all.pdf')
                if (!file.exists(linkedtreefile)){
                    dd  = setup_dendrogram(dend3, ll2, udf, declist=declist, bcutoff=3, altline=0)
                    plot.treebasic(linkedtreefile, ll2, dd, suid, lab=lab, gg=gg)
                }
            }
        }
    }
}

print(paste("[STATUS] Finished plotting all uids for chunk", chunk))
