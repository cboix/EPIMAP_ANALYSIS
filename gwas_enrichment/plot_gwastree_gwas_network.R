#!/usr/bin/R
# -------------------------------------------------
# Plot Figure 3 network for tree-based enrichments:
# Updated on 07/26/21
# -------------------------------------------------
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
library(igraph)

# todo: make networks from regmat (condense gwas <- gwas)

# -------------------------------
# load in the enrichment results:
# -------------------------------
# plotting directories:
eximgdir = paste0(img, 'gwas_tree_analysis/')

usehg = TRUE
useset = 'cons'
ksuff = 'allgwas'
pcut = 5
if (usehg){
    resdf = read.delim(paste0('tree_', useset, '_snpcentric_enrichments.tsv'), header=T)
    tol = 2500
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

# -------------------
# Load tree metadata:
# -------------------
usetree = 'enhancers'
treerdafile = 'Enhancer_jaccard_tree.Rda'
load(treerdafile)

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

tmat = t(sapply(1:NN, function(x){labels(dend3) %in% declist$dec[[x]]}))
ntmat = 1 * t(sapply(1:NN, function(x){odf$GROUP %in% nodetissue$GROUP[nodetissue$node == x]}))
ltmat = 1 * t(sapply(1:NL, function(x){odf$GROUP %in% meta$GROUP[meta$id == labels(dend)[x]]}))
colnames(tmat) = labels(dend3)
colnames(ntmat) = odf$GROUP
colnames(ltmat) = odf$GROUP

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

# -----------------------------------------------------
# Image of the kept GWAS, clustered by node occurrence:
# -----------------------------------------------------
keep.cols = kept.uids

if (pcut == 0.1){
    pcutstr = '_pt1pct'
} else { pcutstr = paste0('_',pcut, 'pct') }
eximgpref = paste0(eximgpref, '_', ksuff, pcutstr)
extpref = paste0('snpcentric_hg_',useset, pcutstr, '_')

redmat = -log10(regmat[keep.cols,])
lpcut = -log10(pcut / 100)
redmat[redmat < lpcut] = 0
Zmat = 1 * (redmat >= lpcut)
leafmat = redmat %*% tmat
Zmat = 1 * (redmat >= lpcut) %*% ntmat
kuid = rownames(Zmat)[apply(Zmat, 1, max) > 0]
Zmat = Zmat[kuid,]

Zmat = sweep(Zmat, 1, apply(Zmat, 1, sum), '/')
sapply(c(0.5, 0.75, 1), function(x){sum(apply(Zmat, 1, max) >= x)})
# hist(apply(Zmat, 1, max), 20, col='grey')

rn = names(which(apply(Zmat, 1, max) < .50))
irn = names(which(apply(Zmat, 1, max) >= .50))
summary(keptgw$sampsize[keptgw$uid %in% rn])
summary(keptgw$sampsize[keptgw$uid %in% irn])
summary(apply(Zmat[rn,] > 0, 1, sum))

prn = names(which(apply(Zmat > 0, 1, sum) >= 34 * .25))
mean(apply(Zmat[prn,] > 0,1,sum))

mrn = rn[!(rn %in% prn)]
mean(apply(Zmat[mrn,] > 0,1,sum))

# Label mono-factorial, multi-factorial, poly-factorial.
write.table(names(which(apply(Zmat, 1, max) == 1)), paste0(extpref, 'fdrp1_names_only_1.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")
# Above 75%:
write.table(names(which(apply(Zmat, 1, max) >= .75)), paste0(extpref, 'fdrp1_names_75over.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")
# Above 50%
write.table(names(which(apply(Zmat, 1, max) >= .50)),paste0(extpref, 'fdrp1_names_50over.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")
# Above 50% of GROUPS (poly-factorial)
write.table(names(which(apply(Zmat>0, 1, sum) >= 34 * .25)), paste0(extpref, 'fdrp1_names_9overgroup.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")
write.table(names(which(apply(Zmat>0, 1, sum) >= 34 * .5)), paste0(extpref, 'fdrp1_names_17overgroup.tsv'),
            quote=F, row.names=F, col.names=F, sep="\t")

# Turn into tissue matrix - either normalize or reduce to hit:
# Ztmat = Zmat %*% ntmat
Ztmat = redmat[kuid,] %*% ntmat
Ztmat = sweep(Ztmat, 1, apply(Ztmat, 1, sum), '/')

# Statistics of single tissue/multi-tissue:
groupord = odf$GROUP
Ztmat = Ztmat[,groupord] # Reverse for plotting
# # Diagonalize:
dt = dist(Ztmat, 'ejaccard')
dt = dist(Ztmat, 'cosine')
ht <- hclust(dt, method='ward.D')
ht <- hclust(dt, method='complete')
cocl <- order.optimal(dt, ht$merge)$order
rn <- names(cocl)[cocl]
Ztmat = Ztmat[rn,]

ord = order(apply(Ztmat, 1, which.max), decreasing=TRUE)
Ztmat = Ztmat[ord,]
cto = apply(Ztmat, 1, which.max)
diagcto = cto
breaks = calc.breaks.acut(cto)
blabels = colnames(Ztmat)[unique(cto)]
bcol = odf$COLOR[unique(cto)]
ends = (c(0, nrow(Ztmat)) - .5) / (nrow(Ztmat) - 1)
bloc = (c(ends[1], breaks) + c(breaks, ends[2])) / 2
# Break by category diff:
colat = which(diff(as.numeric(odf$category)) == 1)
colbreaks = (colat - .5) / (ncol(Ztmat) - 1)

# Flip for plotting:
Ztmat = t(Ztmat)

png(paste0(eximgpref, 'gwas_red_clust.png'),res=450,units='in',width=5,height=18)
par(mar=c(.5,15,3,.5))
par(yaxs='i')
par(xaxs='i')
image(Ztmat, axes=F, col=col1)
yt = seq(0, 1, length.out=ncol(Ztmat))
xt = seq(0, 1, length.out=nrow(Ztmat))
yaxlab=colnames(Ztmat)
text(y=yt, x=par()$usr[1] - 0.01*(diff(par()$usr[1:2])),
     labels=colnames(Ztmat), srt=0, adj=1, xpd=TRUE,cex=1 / ncol(Ztmat) * 100)
text(x=xt, y=par()$usr[4] + 0.002*(diff(par()$usr[3:4])),
     labels=rownames(Ztmat), srt=90, adj=0, xpd=TRUE, cex=.35)
abline(h=breaks, lty='dotted')
abline(v=colbreaks, lty='dotted', lwd=.5)
box()
par(xpd=TRUE)
text(y=bloc, x=.02, adj=0,
     labels=blabels, col=bcol, cex=.75)
dev.off()

# Replace:
repl = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
              "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
colnames(Ztmat)[colnames(Ztmat) == repl] = "26974007 - Chronic inflammatory diseases (pleiotropy)"

plot_sub_Ztmat = function(pltid){
    # Make breaks, labels:
    szmat = Ztmat[,pltid]
    sbreaks = calc.breaks.acut(cto[pltid])
    sblabels = rownames(szmat)[unique(cto[pltid])]
    sbcol = odf$COLOR[unique(cto[pltid])]
    sends = (c(0, ncol(szmat)) - .5) / (ncol(szmat) - 1)
    sbloc = (c(sends[1], sbreaks) + c(sbreaks, sends[2])) / 2
    # Break by category diff:
    colat = which(diff(as.numeric(odf$category)) == 1)
    colbreaks = (colat - .5) / (nrow(szmat) - 1)
    # 
    par(mar=c(0.25,11,2.5,0.25))
    par(yaxs='i')
    par(xaxs='i')
    image(szmat, axes=F, col=col1, useRaster=TRUE)
    yt = seq(0, 1, length.out=ncol(szmat))
    xt = seq(0, 1, length.out=nrow(szmat))
    yaxlab=colnames(szmat)
    text(y=yt, x=par()$usr[1] - 0.01*(diff(par()$usr[1:2])),
         # labels=colnames(szmat), srt=0, adj=1, xpd=TRUE,cex=.5 / ncol(szmat) * 100)
         labels=colnames(szmat), srt=0, adj=1, xpd=TRUE,cex=.27)
    text(x=xt, y=par()$usr[4] + 0.002*(diff(par()$usr[3:4])),
         labels=rownames(szmat), srt=90, adj=0, xpd=TRUE, cex=.3)
    abline(h=sbreaks, lty='dashed', lwd=.25)
    abline(v=colbreaks, lty='dashed', lwd=.25)
    box(lwd=.5)
    par(xpd=TRUE)
    box.pad = 0.008
    xx = space.1d(sbloc, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    xx = space.1d(xx, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    text(y=xx, x=.02, adj=0, labels=sblabels, col=sbcol, cex=.35)
    par(xpd=FALSE)
}


pdf(paste0(eximgpref, 'gwas_red_clust_1x3.pdf'),width=8,height=5)
layout(matrix(1:3, 1,3))
midSR = round(ncol(Ztmat) / 3)
plot_sub_Ztmat((2*midSR+1):ncol(Ztmat))
plot_sub_Ztmat((midSR+1):(2*midSR))
plot_sub_Ztmat(1:midSR)
dev.off()


# -------------------------------------
# GWAS network: 
# Color is maximal of mult
# Adjacency is cosine distance
# -------------------------------------
# Use the normalized node matrix:
# TODO: do the sharing/colors directly on the node x GWAS matrix?
# binary.mat=TRUE # Old figs not binary.
binary.mat=FALSE # Old figs not binary.
all.nodes=FALSE
metric = 'jaccard'
# metric = 'cosine'
npref = paste('gwas', metric, 'network', sep="_")

redmat = -log10(regmat)
redmat[redmat < lpcut] = 0
if (binary.mat){
    Ztmat = 1 * (redmat[kuid,] > 0) %*% ntmat
    npref = paste0(npref, '_binary')
} else {
    Ztmat = redmat[kuid,] %*% ntmat
}
if (metric == 'jaccard'){ 
    Ztmat[Ztmat > 0.05] = 1 
    # cutoff = 0.49
    cutoff = 0.55
} else {
    cutoff = 0.75
    # cutoff = 0.7
    # cutoff = 0.65
}
npref = paste0(npref, '_', sub("\\.", "_", cutoff))

dt = dist(Ztmat, metric)
dmat = 1 - as.matrix(dt)
rn = rownames(dmat)
dmat = data.frame(dmat)
colnames(dmat) = rn
rownames(dmat) = rn
dmat$T1 = rownames(dmat)
ddf  = gather(dmat, T2, sim, -T1)
ddf = ddf[ddf$sim >= cutoff,]
ddf = ddf[ddf$T1 != ddf$T2,]
if (all.nodes){
    nodes = sort(guid)
    npref = paste0(npref, '_allnodes')
} else {
    nodes = sort(unique(c(ddf$T1, ddf$T2)))
}
# Remove edges in opposite direction
ddf$T1 = factor(ddf$T1, levels=nodes)
ddf$T2 = factor(ddf$T2, levels=nodes)
ddf = ddf[as.numeric(ddf$T1) < as.numeric(ddf$T2),]
# Kept pct:
dim(ddf)[1]
dim(ddf)[1] / (dim(dmat)[1]^2)
m = apply(ddf, 1, function(x){which.max(Ztmat[x[1],] * Ztmat[x[2],])})
ddf$COLOR = odf$COLOR[m]
diagcto = apply(Ztmat, 1, which.max)
ctodf = data.frame(node=names(diagcto), GROUP=odf$GROUP[diagcto], COLOR=odf$COLOR[diagcto])
rownames(ctodf) = ctodf$node

sum(ctodf[as.character(ddf$T1), 'COLOR'] ==
    ctodf[as.character(ddf$T2), 'COLOR']) /nrow(ddf)

# Simple network: just the links/points:
sdf = ddf
net <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
vcol = as.character(ctodf[nodes,'COLOR'])
vcol[is.na(vcol)] = 'black'
ecol = sapply(sdf$COLOR, alpha=0.25, tsp.col)
V(net)$size = 2
V(net)$label = ""
V(net)$color = vcol
V(net)$frame.color <- 'black' # vcol
V(net)$frame.color <- NA
V(net)$pch = 19
E(net)$color = ecol 
elty = rep('dotted', length(sdf$sim))
elty[sdf$sim >= .85] = 'dashed'
elty[sdf$sim >= .95] = 'solid'
E(net)$lty = elty
ewidth = (sdf$sim >= .95) * .4 + (sdf$sim >= .85) * 0.4 + (sdf$sim >= 0.75) * 0.4 + 0.4
E(net)$width = sdf$sim
# E(net)$weight = sdf$sim * .5
E(net)$weight = sdf$sim * .075
# E(net)$weight = sdf$sim * .5
# set.seed(2)
set.seed(8)
l <- layout_with_fr(net, grid='nogrid') # Usually best

png(paste0(eximgpref, npref, '.png'),res=450,units='in',width=6,height=6)
# pdf(paste0(eximgpref, npref, '.pdf'), width=6, height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(net, layout=l, curved=F)
# plot(net, layout=l, edge.curved=seq(-0.5, 0.5, length = ecount(net)))
dev.off()

# Pie chart experimental
pie.values = lapply(nodes, function(x){ Ztmat[x,] })
maxp = apply(redmat[nodes,], 1, max)
bkcut = c(lpcut, 10, 25, 50, 1000)
fmaxp = cut(maxp, breaks=bkcut, right=FALSE)
maxv = (as.numeric(fmaxp) + 1) / 2
maxv = (log10(maxp) + 0.5) * 1.25
csizes = sapply(bkcut[1:4], function(x){ min(maxv[which(maxp >= x)])})

pdf(paste0(eximgpref, npref, '_pie_linetypes_discretesizepts.pdf'),width=6,height=6)
sp = 0.1
par(mar = rep(sp,4))
plot(net, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR),
     vertex.pie.border=vcol,
     vertex.pie.lty=1,
     vertex.size=maxv)
legend('bottomright', 
       legend=c('max -log10p >= 3', 'max -log10p >= 10', 'max -log10p >= 25', 'max -log10p >= 50', 
                'Cosine Sim. >= 0.75', 'Cosine Sim. >= 0.85', 'Cosine Sim. >= 0.95'),
       col='black', lty=c(NA, NA, NA, NA, 'dotted','dashed','solid'),
       cex=.4, pch=c(19, 19,19,19, NA, NA, NA), pt.cex=c(1, 1.5, 2, 2.5, .75, .75, .75) / 4, 
       inset=0, box.col=NA)
legend('bottomleft', legend=odf$GROUP, col=odf$COLOR, ncol=5,
       cex=.4, pch=15, pt.cex=1, inset=0, box.col=NA)
dev.off()

lnet = net
traits = sub("^[0-9][0-9]* - ", "", nodes)
labs = traits
labs = sapply(labs, width=50,split.text)
V(lnet)$label = NA
V(lnet)$label.cex = .25

l2 = l
lrange = apply(l2, 2, range)
lablocs = sweep(l2, 2, lrange[1,], '-')
lablocs = sweep(lablocs, 2, lrange[2,] - lrange[1,], '/') * 2 - 1

png(paste0(eximgpref, npref, '_pie_labels_traits.png'),res=450,units='in',width=6,height=6)
# pdf(paste0(eximgpref, npref, '_pie_labels_traits.pdf'), width=10, height=10)
sp = 0.1
par(mar = rep(sp,4))
plot(lnet, vertex.shape="pie", vertex.pie=pie.values, layout=l,
     vertex.pie.color=list(odf$COLOR),
     vertex.pie.border=vcol,
     vertex.pie.lty=1,
     # vertex.label.color=vlabcol,
     vertex.size=maxv)
text(x=lablocs[,1], y=lablocs[,2], labels=labs, 
     srt=0, adj=0, xpd=TRUE,cex=.5, col=rgb(0,0,0,.5))
dev.off()








