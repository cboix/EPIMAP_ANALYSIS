#!/usr/bin/R
# -----------------------------------------------------
# Plot small circle figures for tree-based enrichments:
# For updating website, updated 07/26/21
# -----------------------------------------------------
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


# -------------------------------
# Load in the enrichment results:
# -------------------------------
# Plotting directories:
eximgdir = paste0(img, 'gwas_tree_analysis/examples/snpcentric/')

usehg = TRUE
useset = 'cons'
if (usehg){
    tol = 2500
    resdf = read.delim(paste0('tree_', useset, '_snpcentric_enrichments.tsv'), header=T)
    eximgpref = paste0(eximgdir, 'enhancers_e', tol, '_hg_', useset)
} else {
    break  # TODO: Add the non-hypergeometric comparison results.
}

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
    agg.gwssdf = aggregate(sampsize ~ uid, gwssdf, max)
    agg.gwssdf = merge(agg.gwssdf, nsnpdf)
    keptgw = agg.gwssdf[agg.gwssdf$sampsize >= NIND & agg.gwssdf$pValue >= NSNP,]
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
    save(dend3, NL, NN, lab, labmapping, file=treerdafile)
} else {
    load(treerdafile)
}

nodedf = data.frame(get_nodes_xy(dend3, type = 'rectangle'))
nodedf$node = 1:NN

# -------------------
# Load tree metadata:
# -------------------
usetree = 'enhancers'

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

cdlenfile = paste0('consensus_object_lengths_', usetree, '_062819.Rdata')
load(cdlenfile)

nl = which(declist$isleaf == 1)
dlen = cdlenlist[['diff']]
clen = cdlenlist[['cons']]
fracuq = dlen[nl] / clen[nl]
udf = data.frame(uq=dlen[nl], total = clen[nl], frac=fracuq)
udf$col='grey65'

# --------------------------------------------
# Create a matrix of enrichments for plotting:
# --------------------------------------------
df = resdf[,c('cls','uid','padj', 'nsnp')]
df$cls = df$cls + 1
opt = expand.grid(cls=1:NN, uid=kept.uids)
df = merge(opt, df, all.x=TRUE)
df$padj[is.na(df$padj)] = 1
df$nsnp[is.na(df$nsnp)] = 0
gc()

pwide = spread(df[,c('cls','uid','padj')], cls, padj)
regmat = as.matrix(pwide[,-1])
rownames(regmat) = pwide$uid

nwide = spread(df[,c('cls','uid','nsnp')], cls, nsnp)
nsnp.regmat = as.matrix(nwide[,-1])
rownames(nsnp.regmat) = nwide$uid

# Function to obtain ll object for plotting:
subset.results = function(regmat, nsnp.regmat, suid, cutp=8){
    ll = list()
    ll$df = data.frame(node=1:NN, pout=1)
    # Fill in for specific trait:
    ll$df$pout = regmat[suid,]
    ll$isnp = nsnp.regmat[suid,]
    ll$rawlp = -log10(regmat[suid,])
    ll$df$rawlog10p = -log10(regmat[suid,])
    ll$log10p = ll$rawlp
    ll$log10p[ll$log10p > cutp] = cutp
    ll$df$log10p = ll$log10p 
    return(ll)
}


# ---------------------------------
# Plot all of these as small plots:
# ---------------------------------
# suid = kept.uids[481] # CAD
# suid = kept.uids[557] # SCZ
for (j in 1:length(kept.uids)){
    suid = kept.uids[j]
    cat(j,'\t', suid)
    # Load regression for suid:
    if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
        suid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                      "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
    } 
    strait = unique(gwdf$trait[gwdf$uid == suid])
    spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
    sgwssdf = gwssdf[gwssdf$uid == suid,]
    sinit = sgwssdf$initSample[order(sgwssdf$sampsize, decreasing=T)][1]
    s2init = split.text(sinit, width=90)
    s2trait = split.text(strait, width=50)
    trait = sub(".* - ", "", suid)
    pmid = sub(" - .*", "", suid)
    traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
    traitstrnoparen = gsub("/", "_", gsub("\\)" ,"", gsub("\\(","",traitstr)))

    for (pcut in c(5, 2, 1, 0.1)){
        # Fix with corrected at either threshold:
        if (pcut == 5){
            pcutstr = '_5pct'
        } else if (pcut == 2){
            pcutstr = '_2pct'
        } else if (pcut == 1){
            pcutstr = '_1pct'
        } else {
            pcutstr = '_pt1pct'
        }
        # Truncate results below this adjusted p-value threshold:
        lpcut = -log10(pcut / 100)
        # Subset to trait's results:
        ll = subset.results(regmat, nsnp.regmat, suid)
        ll$log10p[ll$log10p < lpcut] = 0
        if (sum(ll$log10p) > 0){
            cat('\t', pcut)
            # Setup dendrogram:
            dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=lpcut, cutp=8, altline=0)
            dd$dend = try(set(dd$dend, "branches_lwd", .5))
            # Plot small dendrogram
            NTOP = 10
            ldf = ll$df
            ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
            ldf = ldf[1:NTOP,] # TOP N only
            ldf = ldf[ldf$rawlog10p > lpcut,]
            ldf$rank = 1:nrow(ldf)
            ntdf = merge(ldf, nodetissue)
            ntdf$GROUP = paste0(ntdf$rank,'. ', ntdf$GROUP, "\n(", leafrep[ntdf$node], ")")
            ntdf = ntdf[,c('node','GROUP','COLOR')]
            names(ntdf) = c('node','symbol','color')
            ntdf = merge(ntdf, nodedf)
            ntdf$cex=.4

            imgfilename = paste0(eximgpref, '_gwas_small_', traitstrnoparen, '_', pmid, pcutstr, '.png')
            if (!file.exists(imgfilename)){
                cat('...')
                png(imgfilename, res=400, units='in', width=3.5,height=4)
                par(mar=c(1,0,1,0))
                circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
                           plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
                circos.clear()
                if (s2trait != strait){
                    mtext(s2trait, side=3, line=-.5, cex=.5)
                    mtext(spmid, side=3, line=-1, cex=.3)
                } else {
                    mtext(s2trait, side=3, line=0, cex=.7)
                    mtext(spmid, side=3, line=-.75, cex=.5)
                }
                mtext(s2init, side=1, line=0, cex=.3)
                dev.off()
            }
        }
    }
    cat('\n')
}

