#!/usr/bin/R
# ------------------------------------------------------------------
# Calculate the tree-based enrichments, updated to SNP-centric view:
# ------------------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'auxiliary_gwastree_functions.R'))
source(paste0(bindir, 'load_metadata.R'))

library(GenomicRanges)
library(dplyr)
library(cba)
library(ggplot2)
library(argparser)


# --------------------------
# Load in the gwas datasets:
# -------------------------------------------------
# Only keep GWAS with 10+ lead SNPs (after pruning)
# and with at least 10k individuals.
# -------------------------------------------------
gwcatfile = 'gwascatalog_may03_2019_noquotes.txt'
gwrdafile = sub(".txt", ".Rda", gwcatfile)
load(gwrdafile)

tol = 2500  # Plus/minus distance - window for enhancer overlaps

# Re-do pruning, this time to 1Mb:
run.pruning = FALSE
if (run.pruning){
    ut = as.character(unique(gwdf$uid))
    kept.snps = ldply(ut, df=gwdf, dist=1e6, quiet=FALSE, prune.snps)
    print(dim(gwdf))
    gwdf = merge(kept.snps, gwdf)
    print(dim(gwdf))
}

# Subset to high-quality GWAS:
subsetgwas = TRUE
if (subsetgwas){
    NIND = 10000
    NSNP = 10
    nsnpdf = aggregate(pValue ~ uid, gwdf, length)
    gwssdf = aggregate(sampsize ~ uid, gwssdf, max)
    gwssdf = merge(gwssdf, nsnpdf)
    keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]
    kept.uids =  sort(unique(as.character(keptgw$uid)))
} else {
    kept.uids =  sort(unique(as.character(gwdf$uid)))
}
NUID = length(kept.uids)
print(paste("Number of kept uids:", NUID))
gwdf = gwdf[gwdf$uid %in% kept.uids,]
gwgr = GRanges(paste0('chr', gwdf$chrom), IRanges(gwdf$chromStart, gwdf$chromEnd))

TOTSNP = nrow(gwdf)
NCAT = 1000
print(paste0("Pruned to: ", TOTSNP, " snps in ", NUID, " GWAS."))

# ------------------------------------------------
# Load all of the relevant enhancer-tree datasets:
# ------------------------------------------------
print("[STATUS] Loading DHS list")
ddir = 'DHS_Index_WM201902/'
dpref = 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19'
dmlrdafile = paste0(ddir, dpref, '_enh.core.srt.Rda')
# Load indices (needed for mapping, etc.)
enhind = as.numeric(scan('Enhancer_indices.txt', 'c')) + 1
load(dmlrdafile)
dmgr = GRanges(enhdf$chr, IRanges(enhdf$start - tol, enhdf$end + tol), name=enhdf$name)

# Intersections with GWAS:
gwintrdafile = sub(".txt", "_intersections.Rda", gwcatfile)
load(gwintrdafile)

# -----------------
# Load in the tree:
# -----------------
treerdafile = 'Enhancer_jaccard_tree.Rda'
if (!file.exists(treerdafile)){
    emat = read.delim(gzfile('Enhancer_jaccard.tsv.gz'), sep="\t", header=F)
    matnames = scan('Enhancer_matrix_names.txt', "c")
    rownames(emat) = matnames
    colnames(emat) = matnames
    dt <- as.dist(emat)
    method = 'complete'
    ht <- hclust(dt, method=method)
    ht$order <- order.optimal(dt, ht$merge)$order
    dend = as.dendrogram(ht)
    lab = labels(dend)
    info = meta[lab, 'infoline']
    col = meta[lab, 'COLOR']
    group = meta[lab, 'GROUP']
    dend2 = set(dend, "labels", info)
    labels_colors(dend2) <- col
    NCLUST=20
    colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(NCLUST)
    dend3 <- color_branches(dend2, k=NCLUST, col=colpair)
    dend3 = set(dend3, "labels_cex", .18)
    NL = length(labels(dend3))
    memb = get_nodes_attr(dend3, 'member')
    NN = length(memb)
    # For relabeling:
    names(lab) = NULL
    # NOTE: Will only work with enh. may need to fix 
    labmapping = sapply(lab, function(x){which(matnames == x)})
    # Unique labels:
    labels_dend <- labels(dend3)
    if (as.logical(anyDuplicated(labels_dend))) {
        labels(dend3) <- paste0(seq_along(labels_dend), "_", labels_dend)
        labels_dend <- labels(dend3) }
    save(dend, dend3, NL, NN, lab, labmapping, file=treerdafile)
} else {
    load(treerdafile)
}

# ------------------------------------------
# Load large CDLL objects and tree metadata:
# ------------------------------------------
usetree = 'enhancers'
cdllfile = paste0('consensus_object_', usetree, '_062819.Rdata')
load(cdllfile)

cdlenfile = paste0('consensus_object_lengths_', usetree, '_062819.Rdata')
load(cdlenfile)

nbpfile = paste0('consensus_object_nbp_', usetree, '_062819.Rdata')
load(nbpfile)

# Enhancer index map to matrix (intersection on enhancers only):
midpref = paste0('_', usetree, '_e', tol, '_all_')
enhmap = rep(0, max(enhind))
enhmap[enhind] = 1:length(enhind)

# Tree + descendant metadata:
ntmeta.rda = paste0('enhancer_tree_metadata.Rda')
load(ntmeta.rda)

NN = nrow(nodetissue)
leafrep = sapply(1:NN, function(x, max.terms=3){
                     blacklist = c('GLAND', 'TISSUE', 'CELL')
                     x = sub("[0-9]*_", "", declist$dec[[x]])
                     # x = unique(x)
                     nx = length(x)
                     if (nx > 1){
                         term_words <- strsplit(x, "[ _,.]");
                         tab_all <- sort(table(unlist(term_words)));
                         tab_all = tab_all[tab_all > 1]
                         tab_all = tab_all[!(names(tab_all) %in% blacklist)]
                         x = paste0(tolower(head(names(sort(tab_all, decreasing=T)), max.terms)), collapse=', ')
                     } else { x = tolower(x) }
                     x = capitalize(x)
                     return(x)})
lind = which(leafrep == '')
nodetissue = nodetissue[order(nodetissue$node),]
leafrep[lind] = nodetissue$GROUP[lind]


# -------------------------------------------
# Turn cdll$diff and cdll$cons into matrices:
# -------------------------------------------
qdf = suppressWarnings(data.frame(findOverlaps(gwgr, dmgr)))
qdf$uid = gwdf$uid[qdf$queryHits]
NMAPPED = length(unique(qdf$queryHits))
names(qdf)[2] = 'id'
int.enh = sort(unique(qdf$id)) # Only intersected enhancers
print(paste("Number SNPs in any enhancer region:", NMAPPED))

NSNP = nrow(gwdf)
NE = max(int.enh)
NC = length(cdll$cons)

# Matrix: SNPs to Enhancers:
qmat = sparseMatrix(i = qdf[,1], j = qdf[,2], x = rep(1, nrow(qdf)), dims = c(NSNP,NE))

# Matrix: SNPs to UIDs.
NUID = length(kept.uids)
gwdf$uid = factor(gwdf$uid, levels=kept.uids)
umat = sparseMatrix(i = as.numeric(gwdf$uid), j = 1:nrow(gwdf), x = rep(1, nrow(gwdf)), dims = c(NUID, NSNP))

# Matrix: Sets to Enhancers.
smat = list()
tmat = list()
fullmat = list()
for (set in c('cons','diff')){
    ce.num = unlist(sapply(1:NC, function(x){ rep(x, length(cdll[[set]][[x]]))} ))
    ce.vec = enhmap[unlist(cdll[[set]])]
    ce.ind = ce.vec %in% int.enh
    smat[[set]] = sparseMatrix(i = ce.vec[ce.ind], j = ce.num[ce.ind], x = rep(1, sum(ce.ind)), dims = c(NE, NC))

    # Calculate overlaps of SNPs with CLS:
    tmat[[set]] = qmat %*% smat[[set]]
    fullmat[[set]] = tmat
    tmat[[set]]@x = rep(1, length(tmat[[set]]@x)) # Set to unique SNPs captured only
}
rm(ce.num, ce.vec, ce.ind)
# rm(cdll)
gc()


# ------------------------------------------------
# Calculate the hypergeometric enrichment results:
# ------------------------------------------------
run.enhsnp.hyper = function(umat, tmat, seed=1, rand=FALSE){
    # Calculate numsnp in each GWAS x CLS intersection
    set.seed(seed)
    if (rand){
        # Shuffle SNPs, leaving the same # of SNPs per GWAS:
        sind = sample(1:ncol(umat), ncol(umat), replace=FALSE)
        ns.mat = umat[,sind] %*% tmat
    } else { ns.mat = umat %*% tmat }
    # Number of captured snps:
    totuid = rowSums(umat)
    totcls = colSums(tmat)
    # Data frame of number of SNPs for each calc:
    ctdf = as.data.frame(summary(ns.mat))
    ctdf$j = ctdf$j - 1
    names(ctdf) = c('i','cls','nsnp')
    # Alternatively, if we want to keep + test the 0 intersections:
    # matdf = data.frame(as.matrix(ns.mat))
    # matdf$i = 1:nrow(matdf)
    # ctdf = gather(matdf, cls, nsnp, -i)
    # ctdf$cls = as.numeric(sub("^X","",ctdf$cls)) - 1
    ctdf$uid = kept.uids[ctdf$i]
    ctdf$totcls = totcls[ctdf$cls + 1]
    ctdf$totuid = totuid[ctdf$i]
    ctdf$totsnp = length(gwgr)
    # Run SNP-centric hypergeometric test:
    pvdf = ctdf[,c('nsnp','totuid','totcls','totsnp')]
    pout <- apply(pvdf, 1, run.hyper)
    ctdf$p = pout
    ctdf$lp = -log10(pout)
    ctdf = ctdf[order(ctdf$p),]
    ctdf$padj = p.adjust(ctdf$p,'BH') # Correct by BH
    return(ctdf)
} 

hgenr.rda = 'tree_snpcentric_hg_enrichments.tsv'
if (!file.exists(hgenr.rda)){
    # If we run hypergeom on cons, get 334 / 262 / 232 (very similar)
    t1 = proc.time()
    resdf = run.enhsnp.hyper(umat, tmat[['cons']], seed=1, rand=FALSE)
    elt = proc.time() - t1
    cat(paste0(round(elt[3], 2), 's elapsed\n'))
    resdf$lbl = leafrep[resdf$cls + 1]

    dim(resdf[resdf$padj < 0.01,])
    length(unique(resdf$uid[resdf$padj < 0.01]))
    length(unique(resdf$uid[resdf$padj < 0.02]))
    length(unique(resdf$uid[resdf$padj < 0.05]))

    # If we run hypergeom on diff, get 324 / 256 / 213 (also very similar, and overall far fewer combinations)
    t1 = proc.time()
    resdiff.df = run.enhsnp.hyper(umat, tmat[['diff']], seed=1, rand=FALSE)
    elt = proc.time() - t1
    cat(paste0(round(elt[3], 2), 's elapsed\n'))
    resdiff.df$lbl = leafrep[resdiff.df$cls + 1]

    dim(resdiff.df[resdiff.df$padj < 0.01,])
    length(unique(resdiff.df$uid[resdiff.df$padj < 0.01]))
    length(unique(resdiff.df$uid[resdiff.df$padj < 0.02]))
    length(unique(resdiff.df$uid[resdiff.df$padj < 0.05]))

    # Save each set of enrichments:
    write.table(resdf, file='tree_cons_snpcentric_enrichments.tsv', quote=F, row.names=F, sep="\t")
    write.table(resdiff.df, file='tree_diff_snpcentric_enrichments.tsv', quote=F, row.names=F, sep="\t")
    save(resdf, resdiff.df, file=hgenr.rda)
} else { load(hgenr.rda) }

# Look at specific traits, such as CAD:
trait = kept.uids[481] # CAD
subdf = resdf[resdf$uid == trait & resdf$padj < 0.01,]
# subdf = resdf[resdf$uid == trait,]
subdiff.df = resdiff.df[resdiff.df$uid == trait & resdiff.df$padj < 0.01,]
subdiff.df = resdiff.df[resdiff.df$uid == trait,]

# Note that liver is not enriched in CAD at a high level, despite strong evidence at the locus level.

# --------------------------------------------------------------------
# Count the total number of bases covered and perform a binomial test:
# --------------------------------------------------------------------
dmgr.red = reduce(dmgr)
dmwidth = sum(dmgr.red@ranges@width)
enh.widths = sapply(1:NN, function(x){
                        ind = enhmap[cdll$cons[[x]]]
                        dmgr.red = reduce(dmgr[ind])
                        return(sum(dmgr.red@ranges@width)) })

enh.widths = sapply(1:NN, function(x){
                        ind = enhmap[cdll$cons[[x]]]
                        dmgr.red = reduce(dmgr[ind])
                        out = c(sum(dmgr.red@ranges@width),
                                sum(dmgr[ind]@ranges@width))
                        return(out)})

cv = GRanges(coverage(dmgr))
cv = cv[cv$score > 0]

# Look at which reduce more or less - is this an issue?
rwdf = resdf
rwdf = resdiff.df
rwdf$width = enh.widths[1, rwdf$cls + 1]
rwdf$fullwidth = enh.widths[2, rwdf$cls + 1]
rwdf$totwidth = dmwidth

rwdf$est.p = rwdf$totuid / rwdf$totwidth # Overall estimated rate
# Estimate of rate accounting for variability in epigenomes:
rwdf$est.p2 = rwdf$totcls / rwdf$width * rwdf$totuid / rwdf$totsnp 
# hist(log2(rwdf$est.p / rwdf$est.p2))

tdf = rwdf[,c('nsnp', 'width','est.p2')]
poiss.pvals = apply(tdf, 1, function(x){
                        x = as.numeric(x)
                        # NOTE: two.sided is very slow
                        poisson.test(x[1], T=x[2], r=x[3], 't')$p.value
                        })
gc()

rwdf$poiss.pval = poiss.pvals
rwdf$poiss.padj = p.adjust(rwdf$poiss.pval, 'fdr')
rwdf$expsnp = ceiling(rwdf$totuid / rwdf$totwidth * rwdf$width)
rwdf$logFC = log2(rwdf$nsnp / rwdf$expsnp)
rwdf = rwdf[order(rwdf$poiss.padj),]
rwdf = rwdf[order(-rwdf$logFC),]

sum(rwdf$poiss.padj < 0.01)
length(unique(rwdf$uid[rwdf$poiss.padj < 0.01])) # Fewer with est.p2 and slightly more with est.p
length(unique(rwdf$uid[rwdf$poiss.padj < 0.01 & (rwdf$expsnp < rwdf$nsnp)])) # Fewer with est.p2 and slightly more with est.p
length(unique(rwdf$uid[rwdf$padj < 0.01]))

subdf = rwdf[rwdf$uid == trait & rwdf$poiss.padj < 0.05 & rwdf$logFC > 0,]
subdf = subdf[order(subdf$poiss.padj),]

png(paste0(img, 'poissontest_pvalhist.png'), res=450, units='in', width=10, height=4.5)
par(mfrow=c(1,2), yaxs='i', xaxs='i')
hist(rwdf$poiss.pval, main='Poisson Test p-values')
hist(rwdf$poiss.padj, main='Poisson Test p-adj')
dev.off()


# Binomial test:
# --------------
binom.pvals = apply(tdf, 1, function(x){
                        x = as.numeric(x)
                        # NOTE: two.sided is very slow
                        binom.test(x=x[1], n=x[2], p=x[3], alternative='greater')$p.value
                        })
gc()

hist(binom.pvals)

rwdf$binom.pval = binom.pvals
rwdf$binom.padj = p.adjust(rwdf$binom.pval, 'fdr')
rwdf = rwdf[order(rwdf$binom.pval),]

sum(rwdf$binom.padj < 0.01)
length(unique(rwdf$uid[rwdf$binom.padj < 0.01])) # Fewer with est.p2 and slightly more with est.p
length(unique(rwdf$uid[rwdf$padj < 0.01]))

png(paste0(img, 'binomialtest_pvalhist.png'), res=450, units='in', width=10, height=4.5)
par(mfrow=c(1,2), yaxs='i', xaxs='i')
hist(rwdf$binom.pval, main='Binomial Test p-values')
hist(rwdf$binom.padj, main='Binomial Test p-adj')
dev.off()

# ------------------------------------
# Alternatively, try rate ratio tests:
# ------------------------------------
# Epigenome vs. all
library(rateratio.test)
tdf = rwdf[,c('nsnp', 'width','totuid', 'totwidth')]
rr1.pval = apply(tdf, 1, function(x){
                    x = as.numeric(x)
                    rateratio.test(x=c(x[1], x[3]), n=c(x[2], x[4]), 
                                   alternative='greater')$p.value })

rwdf$rr1.pval = rr1.pval
rwdf$rr1.padj = p.adjust(rwdf$rr1.pval, 'fdr')
rwdf = rwdf[order(rwdf$rr1.pval),]

sum(rwdf$rr1.padj < 0.01)
length(unique(rwdf$uid[rwdf$rr1.padj < 0.01]))
length(unique(rwdf$uid[rwdf$padj < 0.01]))

png(paste0(img, 'rr1_test_pvalhist.png'), res=450, units='in', width=10, height=4.5)
par(mfrow=c(1,2), yaxs='i', xaxs='i')
hist(rwdf$rr1.pval, main='Rate Ratio Test p-values')
hist(rwdf$rr1.padj, main='Rate Ratio Test p-adj')
dev.off()


# Rate-ratio test is not powered for the hypergeometric test scenario:
# --------------------------------------------------------------------
tdf = rwdf[,c('nsnp', 'totcls','totuid', 'totsnp')]
rr0.pval = apply(tdf, 1, function(x){
                    x = as.numeric(x)
                    rateratio.test(x=c(x[1], x[3]), n=c(x[2], x[4]), 
                                   alternative='two.sided')$p.value })

rwdf$rr0.pval = rr0.pval
rwdf$rr0.padj = p.adjust(rwdf$rr0.pval, 'fdr')
rwdf = rwdf[order(rwdf$rr0.pval),]

sum(rwdf$rr0.padj < 0.01)
length(unique(rwdf$uid[rwdf$rr0.padj < 0.01]))
length(unique(rwdf$uid[rwdf$padj < 0.01]))

# Rate-ratio test against expectation given properties of dataset:
# ----------------------------------------------------------------
rwdf$expsnp = ceiling(rwdf$totuid * rwdf$totcls / rwdf$totsnp)

plot(rwdf$nsnp, rwdf$expsnp, pch='.')
abline(0,1)

# tdf = rwdf[,c('nsnp', 'totcls','totuid', 'totsnp')]
tdf = rwdf[,c('nsnp', 'width','expsnp', 'width')]

rr0.pval = apply(tdf, 1, function(x){
                    x = as.numeric(x)
                    rateratio.test(x=c(x[1], x[3]), n=c(x[2], x[4]), 
                                   alternative='two.sided')$p.value })

rwdf$rr0.pval = rr0.pval
rwdf$rr0.padj = p.adjust(rwdf$rr0.pval, 'fdr')
rwdf = rwdf[order(rwdf$rr0.pval),]

sum(rwdf$rr0.padj < 0.01)
length(unique(rwdf$uid[rwdf$rr0.padj < 0.01]))
length(unique(rwdf$uid[rwdf$padj < 0.01]))

aa = rwdf[rwdf$expsnp < rwdf$nsnp,]
head(aa[,c('nsnp','width','expsnp','totwidth','rr0.padj', 'uid','lbl')])

rwdf$ratio = rwdf$nsnp / rwdf$width

hist(rwdf$ratio[rwdf$cls == 8], 200)
hist(rwdf$ratio[rwdf$cls == 408], 200)


# Epigenome vs. non-epigenome for each epigenome: *******
# -----------------------------------------------
tdf = rwdf[,c('nsnp', 'width','totuid', 'totwidth')]
rr2.pval = apply(tdf, 1, function(x){ x = as.numeric(x)
                 rateratio.test(x=c(x[1], x[3] - x[1]), n=c(x[2], x[4] - x[2]), 
                                alternative='t')$p.value })

rwdf$rr2.pval = rr2.pval
rwdf$rr2.padj = p.adjust(rwdf$rr2.pval, 'fdr')
rwdf = rwdf[order(rwdf$rr2.pval),]

sum(rwdf$rr2.padj < 0.01)
length(unique(rwdf$uid[rwdf$rr2.padj < 0.01])) # Slightly more, not by much.
length(unique(rwdf$uid[rwdf$padj < 0.01]))

rwdf$logFC = log2((rwdf$nsnp / rwdf$width) / (rwdf$totuid / rwdf$totwidth))
length(unique(rwdf$uid[rwdf$rr2.padj < 0.01 & rwdf$logFC > 0])) # Slightly more, not by much.
subdf = rwdf[rwdf$uid == trait & rwdf$rr2.padj < 0.01 & rwdf$logFC > 0,]

png(paste0(img, 'rr2_test_pvalhist.png'), res=450, units='in', width=10, height=4.5)
par(mfrow=c(1,2), yaxs='i', xaxs='i')
hist(rwdf$rr2.pval, main='Rate Ratio Test (In/out) p-values')
hist(rwdf$rr2.padj, main='Rate Ratio Test p-adj')
dev.off()



# Prop.test on three proportions (child, all, parent):
# ----------------------------------------------------
# TODO: include parent + children with 0 SNPs
rwdf$parent = declist$parent[rwdf$cls + 1]
pdf = rwdf[,c('cls','width','uid','nsnp')]
names(pdf) = c('parent','pwidth','uid','pnsnp')
pdf$parent = pdf$parent + 1
rwdf = merge(rwdf, pdf) # TODO: Missed rwdf??

# Doesn't really address issue...
tdf = rwdf[,c('nsnp', 'width','totuid', 'totwidth', 'pnsnp','pwidth')]
# tdf = rwdf[,c('nsnp', 'width','pnsnp','pwidth')]
rr2.pval = apply(tdf, 1, function(x){ 
                     x = as.numeric(x)
                     prop.test(x=c(x[1], x[3] - x[1],x[5]), n=c(x[2], x[4] - x[2], x[6]), 
                     # prop.test(x=c(x[1], x[3]), n=c(x[2], x[4]), 
                               alternative='t')$p.value })

rwdf$rr2.pval = rr2.pval
rwdf$rr2.padj = p.adjust(rwdf$rr2.pval, 'fdr')
rwdf = rwdf[order(rwdf$rr2.pval),]

sum(rwdf$rr2.padj < 0.01)
length(unique(rwdf$uid[rwdf$rr2.padj < 0.01])) # Slightly more, not by much.
length(unique(rwdf$uid[rwdf$padj < 0.01]))

# rwdf$logFC = log2((rwdf$nsnp / rwdf$width) / (rwdf$totuid / rwdf$totwidth))
rwdf$logFC = log2((rwdf$nsnp / rwdf$width) / (rwdf$pnsnp / rwdf$pwidth))
length(unique(rwdf$uid[rwdf$rr2.padj < 0.01 & rwdf$logFC > 0])) # Slightly more, not by much.
subdf = rwdf[rwdf$uid == trait & rwdf$rr2.padj < 0.01 & rwdf$logFC > 0,]






# Stratify by number of SNPs, adjust separately:
# ----------------------------------------------
rwdf$lt = log(rwdf$totuid + 1)
bks = seq(min(rwdf$lt), max(rwdf$lt), length.out=10)
rwdf$loguid.bins = cut(rwdf$lt, bks, include.lowest=TRUE)

gplot = ggplot(rwdf, aes(rr2.pval)) + 
    facet_wrap(~loguid.bins) + 
    geom_histogram() + 
    labs(x='Rate-ratio p-value', y='Number of tests') + 
    theme_bw()
ggsave(paste0(img, 'rr2_test_pvalhist_strat.png'), gplot, dpi=450, units='in', width=8, height=5)

# Epigenome vs. non-epigenome for each epigenome:
# -----------------------------------------------
tdf = rwdf[,c('nsnp', 'width','totuid', 'totwidth')]
rr3.pval = apply(tdf, 1, function(x){ x = as.numeric(x)
                 rateratio.test(x=c(x[1], x[3] - x[1]), n=c(x[2], x[4] - x[2]), 
                                alternative='two.sided')$p.value})

# Remove 1s to properly actually estimate the pi_0
rwdf$rr3.pval = rr3.pval
ind = which(rwdf$rr3.pval < 1)
rwdf$rr3.padj = 1
rwdf$rr3.padj[ind] = p.adjust(rwdf$rr3.pval[ind], 'fdr') 
# rwdf$rr3.padj = p.adjust(rwdf$rr3.pval, 'fdr')
rwdf = rwdf[order(rwdf$rr3.pval),]

ind = rwdf$rr3.padj < 0.01 & rwdf$log2fc >= 0
sum(ind)
length(unique(rwdf$uid[ind])) # Slightly more, not by much.
length(unique(rwdf$uid[rwdf$padj < 0.01]))

png(paste0(img, 'rr3_test_pvalhist_twosided.png'), res=450, units='in', width=10, height=4.5)
par(mfrow=c(1,2), yaxs='i', xaxs='i')
hist(rwdf$rr3.pval, main='Rate Ratio Test (two-sided) p-values')
hist(rwdf$rr3.padj, main='Rate Ratio Test p-adj')
dev.off()

rwdf$log2fc = with(rwdf, log2((nsnp / width) /
                              (totuid/ totwidth)))

par(mfrow=c(1,2), yaxs='i', xaxs='i')
hist(rwdf$rr3.pval[rwdf$log2fc >= 0], main='Rate Ratio Test (two-sided) p-values')
hist(rwdf$rr3.padj[rwdf$log2fc >= 0], main='Rate Ratio Test p-adj')

# TODO: Calibrate FDR - is it OK? 
trait = kept.uids[481] # CAD
subdf = rwdf[rwdf$uid == trait & rwdf$rr3.padj < 0.01 & rwdf$log2fc > 0,]
# subdf = resdf[resdf$uid == trait,]
subdiff.df = resdiff.df[resdiff.df$uid == trait & resdiff.df$padj < 0.01,]
subdiff.df = resdiff.df[resdiff.df$uid == trait,]


# TODO: Also compare rateratio.test with prop.test






# ------------------------------------------------------------------------------
# TODO: ****
# Permute and test that the FDR is properly controlled with the rate ratio test:
# ------------------------------------------------------------------------------















# ---------------------------------------------------
# Try random with this protocol (prob. very inflated)
# ---------------------------------------------------
set = 'diff'
sind = sample(1:ncol(umat), ncol(umat), replace=FALSE)
ns.mat = umat[,sind] %*% tmat[[set]]
# Number of captured snps:
kept.snps = which(rowSums(tmat[['diff']]) > 0)
totuid = rowSums(umat[,kept.snps])
totcls = colSums(tmat[[set]])
# Data frame of number of SNPs for each calc:
ctdf = as.data.frame(summary(ns.mat))
ctdf$j = ctdf$j - 1
names(ctdf) = c('i','cls','nsnp')
ctdf$uid = kept.uids[ctdf$i]
ctdf$totuid = totuid[ctdf$i]
ctdf$width = enh.widths[ctdf$cls + 1]
ctdf$totwidth = dmwidth

pvdf = ctdf[,c('nsnp','totuid','width','totwidth')]
pout <- apply(pvdf, 1, run.hyper)
ctdf$p = pout
ctdf$lbl = leafrep[ctdf$cls + 1]
ctdf = ctdf[ctdf$width != 0,]
ctdf = ctdf[order(ctdf$p),]
ctdf$padj = p.adjust(ctdf$p, 'BH')

# Very poorly separated - very poor FDR control
length(unique(ctdf$uid[ctdf$padj < 0.01]))
length(unique(ctdf$uid[ctdf$padj < 0.05]))

pcut = 10^c(-seq(1, 8, .5))
nreal = sapply(pcut, function(x){length(rwdf$uid[rwdf$padj < x])})
nperm = sapply(pcut, function(x){length(ctdf$uid[ctdf$padj < x])})
rpr = nperm / nreal 

nuid = sapply(pcut, function(x){length(unique(rwdf$uid[rwdf$padj < x]))})


# -----------------------------------
# Try logistic regression approaches:
# -----------------------------------
trait = kept.uids[740] # T2D (no enr.)
trait = kept.uids[792] # ASD (no enr.)
trait = kept.uids[557] # SCZ (weird in all)
trait = kept.uids[481] # CAD
trait = kept.uids[682] # SCZ (worse than hypergeom - back muscle??)
allsnps = 1 * (gwdf$uid == trait)
# kept.uids[grep("Schizophrenia", kept.uids)]

# Numeric overlap - number of enhancers nearby/underlying a SNP:
nmat = list()
for (set in c('cons','diff')){
    nmat[[set]] = qmat %*% smat[[set]]
}


# Direct logistic regression - enhancer neighborhood to SNPs (looks good as well...)
library(Rfast)
ptm = proc.time()
pvals = rep(1, NC)
for (i in 3:NC){
    print(i) 
    x = log(nmat[['cons']][,i])
    fit = try(glm_logistic(x=x, y=allsnps, full=TRUE, maxiters=100), silent=TRUE)
    if (class(fit) != 'try-error'){ pvals[i] = fit$info['x','p-value'] }
}
ptmend = proc.time() - ptm
print(ptmend)

df = data.frame(node=1:NC, p=pvals, leaf=leafrep)
df = df[order(df$p),]
df$padj = p.adjust(df$p,'BH') # Correct by BH


# About 2x slower to run with parent / descendant 
library(Rfast)
ptm = proc.time()
x = matrix(0, ncol=2, nrow=nrow(nmat[['cons']]))
pvals1 = rep(1, NC)
pvals2 = rep(1, NC)
for (i in 3:NC){
    print(i) 
    p = declist$parent[i]
    x[,1] = 1 * (nmat[['cons']][,p] > 0)
    x[,2] = 1 * ((nmat[['cons']][,i] - nmat[['cons']][,p]) > 0)
    fit = try(glm_logistic(x=x, y=allsnps, full=TRUE, maxiters=100),silent=TRUE)
    if (class(fit) != 'try-error'){ pvals1[i] = fit$info['X1','p-value'] }
    if (class(fit) != 'try-error'){ pvals2[i] = fit$info['X2','p-value'] }
}
ptmend = proc.time() - ptm
print(ptmend)

df2 = data.frame(node=1:NC, p1=pvals1, p2=pvals2, leaf=leafrep)
df2 = df2[order(df2$p2),]
df2$padj = p.adjust(df2$p2,'BH') # Correct by BH



# About 2x slower to run with parent / descendant 
library(Rfast)
ptm = proc.time()
x = matrix(0, ncol=2, nrow=nrow(nmat[['cons']]))
pvals1 = rep(1, NC)
pvals2 = rep(1, NC)
for (i in 3:NC){
    print(i) 
    p = declist$parent[i]
    x[,1] = nmat[['cons']][,p]
    x[,2] = nmat[['cons']][,i] - nmat[['cons']][,p]
    fit = try(glm_logistic(x=x, y=allsnps, full=TRUE, maxiters=100),silent=TRUE)
    if (class(fit) != 'try-error'){ pvals1[i] = fit$info['X1','p-value'] }
    if (class(fit) != 'try-error'){ pvals2[i] = fit$info['X2','p-value'] }
}
ptmend = proc.time() - ptm
print(ptmend)

df2 = data.frame(node=1:NC, p1=pvals1, p2=pvals2, leaf=leafrep)
df2 = df2[order(df2$p2),]
df2$padj = p.adjust(df2$p2,'BH') # Correct by BH


# Looking instead at difference:
library(Rfast)
ptm = proc.time()
pvals3 = rep(1, NC)
for (i in 3:NC){
    print(i) 
    fit = try(glm_logistic(x=nmat[['diff']][,i], y=allsnps, full=TRUE, maxiters=100), silent=TRUE)
    if (class(fit) != 'try-error'){ pvals3[i] = fit$info['x','p-value'] }
}
ptmend = proc.time() - ptm
print(ptmend)

df3 = data.frame(node=1:NC, p=pvals3, leaf=leafrep)
df3 = df3[order(df3$p),]
df3$padj = p.adjust(df3$p,'BH') # Correct by BH






library(Rfast)
ptm = proc.time()
x = matrix(0, ncol=2, nrow=nrow(tmat))
pvals2 = rep(1, NC)
for (i in 3:NC){
    print(i) 
    p = declist$parent[i]
    x[,1] = tmat[,p]
    x[,2] = tmat[,i] - tmat[,p]
    # fitp = try(glm_logistic(x=x[,1,drop=F], y=allsnps, full=TRUE, maxiters=100),silent=TRUE)
    fit = try(glm_logistic(x=x, y=allsnps, full=TRUE, maxiters=100),silent=TRUE)
    if (class(fit) != 'try-error'){ pvals2[i] = fit$info['X2','p-value'] }
}
ptmend = proc.time() - ptm
print(ptmend)



# ----------------------
# Run for example trait:
# ----------------------
grep("^2921", kept.uids)
trait = kept.uids[481] # CAD
trait = kept.uids[482]
allsnps = 1 * (gwdf$uid == trait)

ptm = proc.time()
pvals_llk = rep(1, NN)
pvals_alone = rep(1, NN)
pvals = rep(1, NN)
pllk = rep(NA, NN)
dllk = rep(NA, NN)
diffllk = rep(NA, NN)
isnp = rep(0, NN)
verbose = FALSE
for (i in 1:NC){
    print(i)
    p = declist$parent[i]
    nann = tmat[,i]
    pann = tmat[,p]
    dann = nann - pann
    isnp[i] = sum(nann * allsnps)
    gamma = 0.0001
    if (sum(pann) == 0 || sum(nann) == 0 || sum(nann) == sum(pann)){
        pvals[i] = 1
    } else {
        # Run the parent regression:
        pll.condition = (is.na(pllk[p]) || pllk[p] > 0)
        while (is.na(pllk[p]) || pllk[p] > 0){ 
            pmodel = calc_disjoint_lr(y=allsnps, x1=pann, gamma=gamma, calc.pvals=FALSE)
            pllk[p] = pmodel$loglikelihood
            pll.condition = (is.na(pllk[p]) || pllk[p] > 0)
            if (pll.condition){
                if (verbose) { print("NA parent LLK or bad regression, reduced gamma") } 
                gamma = gamma * 0.1
            }
            if (gamma < 10^-10) break
        }
        # Run the descendant regression:
        dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                         || dllk[i] > 0 || abs(diffllk[i]) > 5000)
        while (dll.condition){ 
            dmodel = try(calc_disjoint_lr(y=allsnps, x1=pann, x2=dann, gamma=gamma), silent=TRUE)
            nmodel = calc_disjoint_lr(y=allsnps, x1=nann, gamma=gamma, calc.pvals=T)
            if (class(dmodel) != 'try-error'){
                dllk[i] = dmodel$loglikelihood 
                diffllk[i] = 2 * (dllk[i] - pllk[p])  # NOTE: the LRT uses 2 x log(R1/R0)
                # Check condition:
                dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                                 || dllk[i] > 0 || abs(diffllk[i]) > 1000)
                if (dll.condition){
                    if (verbose) {print("NA diff LLK or bad regression, reduced gamma"); }
                    gamma = gamma * 0.1
                    if (verbose) {print(paste(sapply(c(i, pllk[p], dllk[i], diffllk[i]), digits=2, round)))}
                }
                if (gamma < 10^-10) break
                # Calculate p-values:
                diffllk[i] = 2 * (dllk[i] - pllk[p])  # NOTE: the LRT uses 2 x log(R1/R0)
                pvals_llk[i] = pchisq(diffllk[i], df=1, lower.tail=FALSE)
                pvals_alone[i] = nmodel$pvals[1]
                pvals[i] = dmodel$pvals[2]
            } else {
                # If super problematic run binomial directly:
                x = matrix(0, ncol=2, nrow=nrow(tmat))
                x[,1] = pann
                x[,2] = dann
                fit = try(glm_logistic(x=x, y=allsnps, full=TRUE, maxiters=100), silent=TRUE)
                if (class(fit) != 'try-error'){ pvals[i] = fit$info['X2','p-value'] }
                break
            }
        }
    }
}
ptmend = proc.time() - ptm
print(ptmend) # about 2 min per run - prob can speed up by cutting parent reg.

df = data.frame(node=1:NN, parent=declist$parent[1:NN], pout=pvals, isnp=isnp)
df$pll = pvals_llk
df$palone = pvals_alone
df = merge(df, nodetissue)

# These match pretty closely if there are more than 10 SNPs
plot(-log10(pvals_llk), -log10(pvals)) 
plot(-log10(pvals_llk), -log10(pvals_alone)) 

# Which ranking makes most sense?
head(df[order(df$pout),], 15)
head(df[order(df$pll),], 10) # p-llk may be slightly better on # SNPs?
head(df[order(df$palone),], 20) # p-llk may be slightly better on # SNPs?

padj = p.adjust(pvals_alone, 'BH')

# ll = pvalsdf.tolist(df, minp, cutp)
# ll$title= title
ll$isnp = isnp
if (return.coeff){ ll$coeff = coeffmat }



# If we want to evaluate concordance with Rfast / other method
x <- model.matrix(y ~ ., data.frame(x))
mod <- .Call(Rfast:::Rfast_glm_logistic, x, y, tol=1e-8, maxiters=100)
names(mod$be) <- colnames(x)
res <- list(be = mod$be, devi = mod$deviance)

library(Rfast)
ptm = proc.time()
x = matrix(0, ncol=2, nrow=nrow(tmat))
pvals2 = rep(1, NC)
for (i in 3:NC){
    print(i) 
    p = declist$parent[i]
    x[,1] = tmat[,p]
    x[,2] = tmat[,i] - tmat[,p]
    # fitp = try(glm_logistic(x=x[,1,drop=F], y=allsnps, full=TRUE, maxiters=100),silent=TRUE)
    fit = try(glm_logistic(x=x, y=allsnps, full=TRUE, maxiters=100),silent=TRUE)
    if (class(fit) != 'try-error'){ pvals2[i] = fit$info['X2','p-value'] }
}
ptmend = proc.time() - ptm
print(ptmend)

plot(-log10(pvals), -log10(pvals2))

# library(RcppNumerical)
trait = kept.uids[787]
allsnps = 1 * (gwdf$uid == trait)
df = data.frame(y=allsnps)
pvals3 = rep(1,NC)
for (i in 3:NC){
    print(i) 
    p = declist$parent[i]
    df$par = tmat[,p]
    df$dec = tmat[,i] - tmat[,p]
    # TODO: Speed up!
    fit = try(glm(y ~ 1 + par + dec, df, family=binomial()))
    if (class(fit) != 'try-error'){
        cfit = coefficients(summary(fit))
        if (nrow(cfit) == 3){ pvals3[i] = cfit['dec','Pr(>|z|)'] }
    }
}

plot(-log10(pvals), -log10(pvals2))

i = 280
plot(-log10(pvals)[1:i], -log10(pvals3)[1:i])

# Rfast p-vals look ok, issue with pvals orig.
plot(-log10(pvals2)[1:i], -log10(pvals3)[1:i])




dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                 || dllk[i] > 0 || abs(diffllk[i]) > 5000)

while (dll.condition){ 
    dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=nann, gamma=gamma, weights=weights)
    dllk[i] = dmodel$loglikelihood 
    diffllk[i] = 2 * (dllk[i] - pllk[p])  # NOTE: the LRT uses 2 x log(R1/R0)
    # Check condition:
    dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                     || dllk[i] > 0 || abs(diffllk[i]) > 1000)
    if (dll.condition){
        if (verbose) {print("NA diff LLK or bad regression, reduced gamma"); }
        gamma = gamma * 0.1
        if (verbose) {print(paste(sapply(c(i, pllk[p], dllk[i], diffllk[i]), digits=2, round)))}
    }
    if (gamma < 10^-10) break
}






dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                 || dllk[i] > 0 || abs(diffllk[i]) > 5000)
while (dll.condition){ 
    # Different purposes:
    if (type == 'union' || against == 'enhancers'){
        if (alt.lr){
            dmodel = alt_lr(y=allsnps, x1=pann, x2=nann)
        } else {
            dmodel = calc_subset_lr(y=allsnps, x1=pann, x2=nann, gamma=gamma, weights=weights)
        }
    } else if (type == 'cons') {
        if (alt.lr){
            dmodel = alt_lr(y=allsnps, x1=pann, x2=dann)
        } else {
            dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=dann, gamma=gamma, weights=weights, max.iter=max.iter)
        }
    } else {
        if (alt.lr){
            dmodel = alt_lr(y=allsnps, x1=pann, x2=nann)
        } else {
            dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=nann, gamma=gamma, weights=weights)
        }
    }
    dllk[i] = dmodel$loglikelihood 
    diffllk[i] = 2 * (dllk[i] - pllk[p])  # NOTE: the LRT uses 2 x log(R1/R0)
    # Check condition:
    dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                     || dllk[i] > 0 || abs(diffllk[i]) > 1000)
    if (dll.condition){
        if (verbose) {print("NA diff LLK or bad regression, reduced gamma"); }
        gamma = gamma * 0.1
        if (verbose) {print(paste(sapply(c(i, pllk[p], dllk[i], diffllk[i]), digits=2, round)))}
    }
    if (gamma < 10^-10) break




pll.condition = (is.na(pllk[p]) || pllk[p] > 0)
while (is.na(pllk[p]) || pllk[p] > 0){ 
    if (alt.lr){
        pmodel = alt_lr(y=allsnps, x1=pann)
    } else {
        pmodel = calc_disjoint_lr(y=allsnps, x1=pann, gamma=gamma, weights=weights, max.iter=max.iter)
    }
    pcoeff[[i]] = pmodel$coefficients
    pllk[p] = pmodel$loglikelihood
    pll.condition = (is.na(pllk[p]) || pllk[p] > 0)
    if (pll.condition){
        if (verbose) { print("NA parent LLK or bad regression, reduced gamma") } 
        gamma = gamma * 0.1
    }
    if (gamma < 10^-10) break
}
gamma = 0.0001 
if (type == 'union' || (!is.null(weights))){ gamma = gamma * 0.01 } # Tends to have high NA issues
dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                 || dllk[i] > 0 || abs(diffllk[i]) > 5000)
while (dll.condition){ 
    # Different purposes:
    if (type == 'union' || against == 'enhancers'){
        if (alt.lr){
            dmodel = alt_lr(y=allsnps, x1=pann, x2=nann)
        } else {
            dmodel = calc_subset_lr(y=allsnps, x1=pann, x2=nann, gamma=gamma, weights=weights)
        }
    } else if (type == 'cons') {
        if (alt.lr){
            dmodel = alt_lr(y=allsnps, x1=pann, x2=dann)
        } else {
            dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=dann, gamma=gamma, weights=weights, max.iter=max.iter)
        }
    } else {
        if (alt.lr){
            dmodel = alt_lr(y=allsnps, x1=pann, x2=nann)
        } else {
            dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=nann, gamma=gamma, weights=weights)
        }
    }
    dllk[i] = dmodel$loglikelihood 
    diffllk[i] = 2 * (dllk[i] - pllk[p])  # NOTE: the LRT uses 2 x log(R1/R0)
    # Check condition:
    dll.condition = (is.na(dllk[i]) || is.infinite(dllk[i]) 
                     || dllk[i] > 0 || abs(diffllk[i]) > 1000)
    if (dll.condition){
        if (verbose) {print("NA diff LLK or bad regression, reduced gamma"); }
        gamma = gamma * 0.1
        if (verbose) {print(paste(sapply(c(i, pllk[p], dllk[i], diffllk[i]), digits=2, round)))}
    }
    if (gamma < 10^-10) break
}
pvals[i] = pchisq(diffllk[i], df=1, lower.tail=FALSE)
if (diffllk[i] > 1000){
    print(paste0("Very large difference in log likelihood reported: Check: ", trait))
    print(paste(sapply(c(i, pllk[p], dllk[i], diffllk[i]), digits=2, round)))
}
if(return.coeff){ coeffmat[i,] = dmodel$coeff }
}
}
close(pb)
ptmend = proc.time() - ptm
print(ptmend)
# Make into dataframe:
df = data.frame(node=1:NN, parent=declist$parent[1:NN], pout = pvals)
ll = pvalsdf.tolist(df, minp, cutp)
ll$title= title
ll$isnp = isnp
if (return.coeff){ ll$coeff = coeffmat }
return(ll)




dmodel = alt_lr(y=allsnps, x1=pann, x2=nann)
} else {
    dmodel = calc_disjoint_lr(y=allsnps, x1=pann, x2=nann, gamma=gamma, weights=weights)
}
}
dllk[i] = dmodel$loglikelihood 
diffllk[i] = 2 * (dllk[i] - pllk[p])  # NOTE: the LRT uses 2 x log(R1/R0)


trait = kept.uids[787]

NE = nrow(enhdf)
pvals = rep(1, NN)
pllk = rep(NA, NN)
dllk = rep(NA, NN)
diffllk = rep(NA, NN)
isnp = rep(0, NN)
# if (return.coeff){ coeffmat = matrix(rep(0, NN * 3), nrow=NN, ncol=3) }
pcoeff = sapply(rep(NA, NN), list)
allsnpset = qdf[qdf$uid == suid, 'subjectHits']
allsnps = rep(0, NE)
allsnps[allsnpset] = 1
title = paste('Node', capitalize(type),'vs.')

nann = rep(0, NE)
pann = rep(0, NE)

if (against == 'parent'){
    p = declist$parent[i]
    pset = enhmap[cdll[[type]][[p]]]
} else if (against == 'sibling'){
    parent = declist$parent[i]
    d2 = which(declist$parent == parent)
    p = d2[d2 != i]
    # TODO: write get sibling
    # TODO: will need different regression logic
} else if (against == 'enhancers'){
    pset = 1:NE
    p = 1 # Put the output here - same parent regression
    # TODO: Will different regression - for difference
    # All enhancers vs. all enhancers + node enhancers (separately)
    # - need to update calc_disjoint_lr to do NONDISJOINT.
}
nset = enhmap[cdll[[type]][[i]]]





ipref = paste0(apref, sprintf("_%05d", plotind))
suid = uids[plotind]
print(paste(suid, '----- output with prefix', ipref))
nsnps = nsnpdf$pValue[nsnpdf$uid == suid]
npref = paste0(apref, sprintf("_permuted_%06d", nsnps))

type = 'cons'
against = 'parent'
weighted = FALSE



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
    }
}

print(paste("[STATUS] Finished plotting all uids for chunk", chunk))
