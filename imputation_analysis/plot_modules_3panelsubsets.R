#!/usr/bin/R
# ----------------------------------------------------------------
# Test plotting the three panel subsets from figure 3 for modules:
# ----------------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'great_enrichment/auxfunctions_GREAT.R'))
library(png)
library(stringr)

# Defaults:
# motiffile='merged_enrichments_0_-99.tsv.gz'
filepref = 'cls_merge2_wH3K27ac100_300'
tagline = 'ChromHMM Enhancers'
imgdir = paste0(img, "clusters/") 

# Arguments:
args=(commandArgs(TRUE))
if (length(args) == 0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    motiffile = args[1]
    filepref = args[2]
    tagline = args[3]
    if (length(args) > 3){ 
        imgdir = args[4] 
    }
}

# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){
    c(filepref, tagline, imgdir) }
source(paste0(bindir, 'load_modules_data.R'))

# Prefix:
cmd = paste('mkdir -p', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'clusters_motif_updated_') # don't overwrite old plots:
# imgpref = paste0(imgdir, 'clusters_motif_')
mpngdir = 'motif_data/logos2/' # TODO: need new motif logo files

# For fixing the binomial enrichments (binom conf interval):
calc_binomial_CI = function(ns, n, z, max=T){
    nf = n - ns
    zsq = z^2
    p = (ns + zsq/2) / (n + zsq)
    pm = (z / (n + zsq)) * sqrt((ns * nf) / n + zsq / 4)
    return(p + (2 * max - 1) * pm)}

# --------------------------------------------
# Cluster center variables for plotting order:
# --------------------------------------------
rnorder = rngroup
acutnamed = acutgroup.nam
set = tagline 
lablist = lablist.group

# Diagonalize the motifs matrix:
dist.order = plot.centers(centlist, set, rnorder, counts=counts, cls=acutnamed, calc.only=TRUE)
clsord = dist.order[[1]] # Cluster order
clscuts = dist.order[[3]]

# ---------------------------------
# Add GREAT enrichment information:
# ---------------------------------
load('parsed_results/cls_merge2_wH3K27ac100_300_filtered_go.Rda')

# Diagonalize, ignore all <zmin
tmin = 2.5
tmat = resmat[,clsord]
tmat[tmat < tmin] = 0
ll = diag.mat(t(tmat))
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
dresmat = resmat[diagord,clsord]
vcut.go =  c(cto[cto ==0], clscuts[cto])
hbreaks = calc.breaks.acut(vcut.go)

plot.go = function(mat, zmax=5, zmin=0.5){
    pmat = t(mat)
    pmat[pmat > zmax] <- zmax
    pmat[pmat < zmin] <- 0
    pmat[pmat < 2] = 1
    pmat[pmat < 3 & pmat >= 2] = 2
    pmat[pmat < 4 & pmat >= 3] = 3
    pmat[pmat >= 4] = 4
    rcols = c('white','yellow','orange','red')
    image(pmat, axes=FALSE,col=rcols, zlim=c(zmin, zmax), useRaster=TRUE)
}

# Write out for res:
write.table(resmat, file=gzfile('parsed_results/modules_GO_enrichments_matrix_ordered.tsv.gz'), sep="\t", quote=F, col.names=T, row.names=T)
write.table(dresmat, file=gzfile('parsed_results/modules_GO_enrichments_matrix_ordered.tsv.gz'), sep="\t", quote=F, col.names=T, row.names=T)


# ---------------------------------
# Add motif enrichment information:
# ---------------------------------
motif.resfile = 'motifs/collated.bkgdhs.cls.enrich.tsv.gz'
motif.rdafile = 'motifs/collated.bkgdhs.cls.enrich.processed.Rda'
motif.fulltsv = 'motifs/collated.bkgdhs.cls.enrich.processed.full.tsv.gz'
motif.redtsv = 'motifs/collated.bkgdhs.cls.enrich.processed.reduced.tsv.gz'
motif.fullnam = 'motifs/collated.bkgdhs.cls.enrich.processed.full.fullnames.tsv'
motif.rednam = 'motifs/collated.bkgdhs.cls.enrich.processed.reduced.fullnames.tsv'
cutoff = 1
pcutoff = -log10(0.005)
if (!file.exists(motif.rdafile)){
    # Read in the new results with updated motifs:
    df = read.delim(gzfile(motif.resfile), header=T)
    map = read.delim('motif_clusters_vierstra2020.tsv', header=T)
    colnames(map)[2] = 'motif'
    df = merge(df, map[,c('Cluster_ID','motif')])
    df$pval_raw[is.na(df$pval_raw)] = 0 
    df$pval_ctrl[is.na(df$pval_ctrl)] = 0 
    df = df[!is.na(df$pfg),]

    # Fix such that all enrichment ratios (enrich/deplete) are more conservative:
    df$rc = log2((df$pfg / df$pbg) / (df$cfg / df$cbg))
    pind = which(df$rc >= 0)
    nind = which(df$rc < 0)
    df$pfrac[pind] = sapply(pind, function(i){calc_binomial_CI(df$pfg[i], df$pbg[i], 1.5, max=F)})
    df$cfrac[pind] = sapply(pind, function(i){calc_binomial_CI(df$cfg[i], df$cbg[i], 1.5, max=T)})
    df$pfrac[nind] = sapply(nind, function(i){calc_binomial_CI(df$pfg[i], df$pbg[i], 1.5, max=T)})
    df$cfrac[nind] = sapply(nind, function(i){calc_binomial_CI(df$cfg[i], df$cbg[i], 1.5, max=F)})
    df$pfrac[is.na(df$pfrac)] = 0
    df$cfrac[is.na(df$cfrac)] = 0
    df$ratio_ctrl = log2(df$pfrac / df$cfrac)
    # Re-do the pval to acct. for depletion. ** (two-sided pval at 0.05)
    enrp = apply(df[,c('pfg','pbg','cfg','cbg')], 1, 
                 function(y){phyper(q=y[1]-1, m=y[3], n=y[4] - y[3], k=y[2], lower.tail=FALSE)})
    depp = apply(df[,c('pfg','pbg','cfg','cbg')], 1, 
                 function(y){phyper(q=y[1], m=y[3], n=y[4] - y[3], k=y[2], lower.tail=TRUE)})
    df$pval_ctrl = apply(cbind(-log10(enrp), -log10(depp)),1, max)

    # Keep motifs with p < 0.01 and log2FC > 1.5
    sdf = df[df$pval_ctrl >= pcutoff & abs(df$ratio_ctrl) >= cutoff,]
    top.ctrl = aggregate(ratio_ctrl ~ Cluster_ID, sdf, function(x){x[which.max(abs(x))]}) 
    cdf = aggregate(motif ~ Cluster_ID, merge(sdf, top.ctrl), function(x){x[1]})
    ctrl.motifs = cdf$motif
    length(ctrl.motifs)

    # Compute full and reduced motif matrices:
    rwide = spread(df[df$motif %in% ctrl.motifs, c('set','motif','ratio_ctrl')], set, ratio_ctrl, fill=0)
    pwide = spread(df[df$motif %in% ctrl.motifs, c('set','motif','pval_ctrl')], set, pval_ctrl, fill=0)
    rwide = spread(df[, c('set','motif','ratio_ctrl')], set, ratio_ctrl, fill=0)
    pwide = spread(df[, c('set','motif','pval_ctrl')], set, pval_ctrl, fill=0)
    rmat = as.matrix(rwide[,-1])
    pmat = as.matrix(pwide[,-1])
    rownames(rmat) = rwide[,1]
    rownames(pmat) = pwide[,1]
    pmat[pmat < 2] = 0
    pmat[pmat >= 2] = 1
    rmat[is.na(rmat)] = 0
    rmat = rmat * pmat[rownames(rmat), colnames(rmat)]
    print(dim(rmat))
    full.mmat = t(rmat)
    rownames(full.mmat) = paste0('c', rownames(full.mmat))
    mmat = full.mmat[, sort(ctrl.motifs)]

    mmap = colnames(mmat)
    mmap = sub("_MA[0-9]*\\.[0-9]", "", mmap)
    mmap = sub("_HUMAN\\.H11MO.*", "", mmap)
    mmap = sub("_.*", "", mmap)
    mmap = toupper(mmap)
    mmap = make.unique(mmap)
    names(mmap) = colnames(mmat)
    colnames(mmat) = mmap[colnames(mmat)]

    full.mmap = colnames(full.mmat)
    full.mmap = sub("_MA[0-9]*\\.[0-9]", "", full.mmap)
    full.mmap = sub("_HUMAN\\.H11MO.*", "", full.mmap)
    full.mmap = sub("_.*", "", full.mmap)
    full.mmap = toupper(full.mmap)
    full.mmap = make.unique(full.mmap)
    names(full.mmap) = colnames(full.mmat)
    colnames(full.mmat) = full.mmap[colnames(full.mmat)]

    # Save motif matrices in quickly loadable files:
    save(mmat, full.mmat, mmap, full.mmap, file=motif.rdafile)
    write.table(mmat, file=gzfile(motif.redtsv), quote=F, sep="\t", row.names=F)
    write.table(full.mmat, file=gzfile(motif.fulltsv), quote=F, sep="\t", row.names=F)

    write.table(as.data.frame(full.mmap), file=motif.fullnam, quote=F, sep="\t", col.names=F)
    write.table(as.data.frame(mmap), file=motif.rednam, quote=F, sep="\t", col.names=F)
} else {
    load(motif.rdafile)
}



# -----------------------------------
# Order motifs, prepare for plotting:
# -----------------------------------
ord = clsord[clsord %in% rownames(mmat)] # Make sure overlaps
ll = diag.mat2(abs(mmat[ord,])) # Diagonalize motifs
col.motiford = ll[[2]]
motnam = col.motiford

ll2 = diag.mat2(abs(full.mmat[ord,])) # Diagonalize motifs
full.motnam = ll2[[2]]

# Margin of go:
gomarg = apply(dresmat, 1, sum)


# ---------------
# Plot resources:
# ---------------
plt.mmat = mmat[clsord,motnam, drop=F]
plt.mmat = plt.mmat[,apply(abs(plt.mmat),2, max) >= 1, drop=F]

png(paste0(imgpref, 'motifs_matrixonly_resource.png'), res=600, units='in', width=7, height=4)
ind = which(abs(plt.mmat) > 0, arr.ind=T)
xat = seq(0,1,length.out=nrow(plt.mmat))
yat = seq(0,1,length.out=ncol(plt.mmat))
vals = round(plt.mmat[ind],1)
zlim = 1.5
plt.mmat[plt.mmat < -zlim] = -zlim
plt.mmat[plt.mmat > zlim] = zlim
par(mar=rep(0,4))
image(plt.mmat, zlim=c(-zlim, zlim), col=rev(colrb), axes=F, useRaster=T)
dev.off()


png(paste0(imgpref, 'centers_matrixonly_resource.png'), res=600, units='in', width=7, height=4)
par(mar=rep(0,4))
cent = centlist[[1]][rnorder, clsord]
image(t(cent), col=col1, axes=F, useRaster=T)
dev.off()


# GO TERMS:
load(paste0("resmat_rownames_filt.RData"))
tmin = 2.5
ind = which(apply(resmat >= 2, 1, sum) < 30 & apply(resmat, 1, max) >= 4)
tmat = resmat[ind, clsord]
tmat[tmat < tmin] = 0

ll = diag.mat2(t(tmat), ratio=.1)
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
dresmat = resmat[diagord,clsord]
# Which do we keep?

png(paste0(imgpref, 'great_matrixonly_resource.png'), res=600, units='in', width=7, height=4)
par(mar=rep(0,4))
plt.gmat = dresmat[,clsord,drop=F]
plot.go(plt.gmat)
dev.off()




# -----------------------------------
# TODO: Test plotting only a subset:
# Get only the core of the module and the core motifs
# Plot each modules' breakdown - fraction from diff. tissues
# Names of each of the samples
# Embryonic etc.
# -----------------------------------
group = 'Heart'
group = 'Brain'
# group = 'Eye'
group = 'HSC & B-cell'
group = 'Placenta & EEM'
group = 'Epithelial'
groups = odf$GROUP
groups = groups[groups != 'Multiple']
figwidth = c()
for (group in groups){ 

    groupstr = tolower(gsub("__", "", str_replace_all(group, "[^[:alnum:]]", "_")))
    print(groupstr)

    # Size of the final figure:
    gind = which(meta[rnorder,'GROUP'] == group)
    subrn = rnorder[gind]
    cind = which(dist.order[[3]] == unique(acutnamed[subrn]))
    # GO:
    plt.gmat = dresmat[,clsord[cind],drop=F]
    goind = which(apply(abs(plt.gmat),1, max) >= 4)
    # Reduce the GO terms to more specific ones:
    if (length(goind) > 1){
        gofrac = apply(abs(plt.gmat[goind,]),1, sum) / gomarg[goind]
        rownames(plt.gmat)
        thresh = .1
        if (sum(gofrac > thresh) > 1){
            goind = goind[gofrac > thresh]
        }
    }
    plt.gmat = plt.gmat[goind,, drop=F]




    # Motifs:
    plt.mmat = mmat[clsord[cind],motnam, drop=F]
    plt.mmat = plt.mmat[,apply(abs(plt.mmat),2, max) >= 1, drop=F]
    add.motifs = (!is.null(ncol(plt.mmat)) & ncol(plt.mmat) > 0)
    add.go = (!is.null(ncol(plt.gmat)) & nrow(plt.gmat) > 0)

    # Define plotting area:
    nplots = 2 + add.motifs + add.go
    widths = c(5, length(cind))
    heights = c(5 + 2, length(subrn))
    if (add.motifs){ heights = c(heights, ncol(plt.mmat)) }
    if (add.go){ heights = c(heights, nrow(plt.gmat) / 3) }
    heights = heights / 2
    w = sum(widths) / 5
    h = sum(heights) / 5
    figwidth = rbind(figwidth, data.frame(group=groupstr, w=w))
    zlim = 1.5

    if (length(cind) > 1){
        pdf(paste0(imgpref, 'groupdist_snapshot_group_', groupstr, '.pdf'), width=w,height=h)
        # png(paste0(imgpref, 'groupdist_snapshot_group_', groupstr, '.png'), res=450, units='in', width=w,height=h)
        layout(matrix(c(1:(2 * nplots)),nrow=nplots,ncol=2), heights=heights, widths=widths, TRUE)
        par(yaxs="i")
        par(xaxs="i")
        lsp = 6
        lsp = 1.5
        par(mar=c(0.25, lsp, 0, 0))
        image(t(as.matrix(metamat[subrn,]) == ""), col='white', axes=F, useRaster=T)
        metaclass = sapply(rev(colnames(metamat)), capitalize)
        text(x=seq(0,1, length.out=ncol(metamat)),
             y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), 
             labels=metaclass, srt=90, adj=0, xpd=TRUE, cex=.7)
        par(mar=c(.25, lsp, 0, 0))
        meta.image(metamat[subrn,5:1], colvals=colvals, cex=0, horiz=T, useRaster=T)
        mtext(paste(length(subrn),group, 'Biosamples'), side=2, cex=ifelse(length(subrn) < 10,0.5,1))
        abline(h=par()$usr[3:4],lty=1,lw=0.5)
        abline(v=par()$usr[1:2],lty=1,lw=0.5)
        # Motif names
        if (add.motifs){
            par(mar=c(0.25, lsp, 0, 0))
            image(plt.mmat, col='white', axes=F, useRaster=T)
            text(y=seq(0, 1, length.out=ncol(plt.mmat)),
                 x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
                 labels=colnames(plt.mmat), srt=0, adj=1, xpd=TRUE,cex=.6)
        }
        # GO term names:
        if (add.go) {
            par(mar=c(0.25, lsp, 0, 0))
            image(t(plt.gmat), col='white', axes=F, useRaster=T)
            text(y=seq(0, 1, length.out=nrow(plt.gmat)),
                 x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]), 
                 labels=rownames(plt.gmat), srt=0, adj=1, xpd=TRUE,cex=.4)
        }
        # Add labels
        par(mar=c(0.25, 0.25, 1.5, 0.25))
        meta.image(enrichmat[clsord[cind],1:5], colvals=colvals, cex=0, horiz=F, useRaster=T)
        box(lwd=0.25)
        mtext(paste(length(cind), group, 'Modules'), side=3, cex=ifelse(length(cind) < 10,0.6,1))
        # Plot clusters and counts
        set = tagline
        par(mar=c(0.25, 0.25, 0, 0.25))
        subcent = centlist[[1]][subrn,clsord[cind]]
        image(t(subcent), col=col1, axes=F, useRaster=T)
        text(y=seq(0, 1, length.out=length(subrn)),
             x=par()$usr[1]+0.001*(par()$usr[2]-par()$usr[1]), 
             labels=paste0(meta[subrn,'id'], ': ', meta[subrn,'infoline']), srt=0, adj=0, xpd=TRUE,cex=.6)
             # labels=meta[subrn,'infoline'], srt=0, adj=0, xpd=TRUE,cex=.6)
        # Plot motifs:
        if (add.motifs){
            ind = which(abs(plt.mmat) > 0, arr.ind=T)
            xat = seq(0,1,length.out=nrow(plt.mmat))
            yat = seq(0,1,length.out=ncol(plt.mmat))
            vals = round(plt.mmat[ind],1)
            plt.mmat[plt.mmat < -zlim] = -zlim
            plt.mmat[plt.mmat > zlim] = zlim
            image(plt.mmat, zlim=c(-zlim, zlim), col=rev(colrb), axes=F, useRaster=T)
            text(xat[ind[,1]], yat[ind[,2]], vals, cex=.65, col=ifelse(abs(vals) >1.25, 'white','black') )
        }
        # Plot GO terms:
        if (add.go){ plot.go(plt.gmat) }
        dev.off()
    }

}


# ---------------------------
# Code resources for website:
# ---------------------------
motif.rdafile = 'motifs/collated.bkgdhs.cls.enrich.processed.Rda'
load(motif.rdafile)
load('parsed_results/cls_merge2_wH3K27ac100_300_filtered_go.Rda')

# Diagonalize, ignore all <zmin
tmin = 2.5
tmat = resmat[,clsord]
tmat[tmat < tmin] = 0
ll = diag.mat(t(tmat))
tmp = ll[[1]]
cto = ll[[3]] # For breaks
diagord = colnames(tmp)
dresmat = resmat[diagord,clsord]
vcut.go =  c(cto[cto ==0], clscuts[cto])
hbreaks = calc.breaks.acut(vcut.go)

# Order motifs, prepare for plotting:
ord = clsord[clsord %in% rownames(mmat)] # Make sure overlaps
ll = diag.mat2(abs(mmat[ord,])) # Diagonalize motifs
col.motiford = ll[[2]]
motnam = col.motiford





# Motif - module bipartite graph:
library(igraph)


# ------------------------------
# Motif - motif adjacency graph:
# ------------------------------
imat = mmat
imat = imat * (imat > 1.5)
mtm = t(imat) %*% imat
rn = rownames(mtm)
mtm = data.frame(mtm)
colnames(mtm) = rn
rownames(mtm) = rn
mtm$M1 = rownames(mtm)
ddf  = gather(mtm, M2, adj, -M1)
ddf = ddf[ddf$M1 != ddf$M2,]
print(dim(ddf))
ddf = ddf[abs(ddf$adj) >= 3,]
print(dim(ddf))
# Remove edges in opposite direction
nodes = sort(unique(c(ddf$M1, ddf$M2)))
ddf$M1 = factor(ddf$M1, levels=nodes)
ddf$M2 = factor(ddf$M2, levels=nodes)
ddf = ddf[as.numeric(ddf$M1) < as.numeric(ddf$M2),]
# Kept pct:
dim(ddf)[1]
dim(ddf)[1] / (dim(mtm)[1]^2)

# TODO: should network just be a gt > 2 * 2 match in any??

# Map max enr match to the group:

# Diag colors for 
m = apply(ddf, 1, function(x){which.max(abs(mmat[,x[1]] * mmat[,x[2]]))})
m = rownames(mmat)[m]
m = enrichmat[m, 'group']
m[m == 'unknown/mixed'] = 'Other'
ddf$COLOR = colvals$group[m]


# Simple network: just the links/points:
# vcol = as.character(ctodf[nodes,'COLOR'])
# vcol[is.na(vcol)] = 'black'
sdf = ddf
net <- graph_from_data_frame(d=sdf, vertices=nodes, directed=F) 
vcol = rep('grey85', length(nodes))
ecol = sapply(sdf$COLOR, alpha=0.75, tsp.col)
V(net)$size = 8
V(net)$color = vcol
V(net)$frame.color <- NA
V(net)$label.color <- 'black'
V(net)$pch = 19
E(net)$color = ecol 
V(net)$label.cex=.4
V(net)$label.font=2

set.seed(2)
l <- layout_with_fr(net, grid='nogrid') # Usually best

png('~/test_motifmotif.png', res=450, units='in', width=6, height=6)
par(mar=rep(0.1, 4))
plot(net, layout=l) 
dev.off()


# elty = rep('dotted', length(sdf$adj))
# elty[abs(sdf$adj) >= .85] = 'dashed'
# elty[abs(sdf$adj) >= .95] = 'solid'
# E(net)$lty = elty
# ewidth = (sdf$adj >= .95) * .4 + (abs(sdf$adj) >= .85) * 0.4 + (abs(sdf$adj) >= 0.75) * 0.4 + 0.4
# E(net)$width = abs(sdf$adj)
# E(net)$weight = abs(sdf$adj) * .5


# Colors for clusters:
clsnam = rownames(mmat)
clsgrp = enrichmat[clsnam, 'group']
clsgrp[clsgrp == 'unknown/mixed'] = 'Other'
clscols = colvals$group[clsgrp]
names(clscols) = clsnam

# Max cols:
mxrn = apply(centers[,clsnam], 2, which.max)
mxrn = rownames(centers)[mxrn]
clsmaxcols = meta[mxrn,'COLOR']
names(clsmaxcols) = clsnam

# Single threshold:
imat = mmat
# imat = full.mmat
imat[abs(imat) < 1.5] = 0
imat = imat
color.edges = TRUE
color.max = FALSE

# Reduce to containing motif:
imat = imat[apply(abs(imat) > 0, 1, sum) > 0,]
imat = imat[,apply(abs(imat) > 0, 2, sum) > 0]
dim(imat)

# Make and modify graph:
net2 <- graph_from_incidence_matrix(imat, weighted=TRUE)
ewt = E(net2)$weight
E(net2)$weight = 1
V(net2)$color <- c("steel blue", "orange")[V(net2)$type+1]
V(net2)$shape <- c("square", "circle")[V(net2)$type+1]
V(net2)$frame.color = NA
# Labels + cols:
lbnames = attributes(V(net2))$names
lind = lbnames %in% rownames(imat)
if (color.max){
    V(net2)$color[lind] = clsmaxcols[lbnames[lind]]
} else {
    V(net2)$color[lind] = clscols[lbnames[lind]]
}
V(net2)$label[lind] = ""
V(net2)$label[!lind] = lbnames[!lind]
V(net2)$label.cex=.4
V(net2)$label.font=2
V(net2)$label.family=2

E(net2)$width <- abs(ewt) / 4
ewt[ewt > 3] = 3
ewt[ewt < -3] = -3
col_fun = function(x, pal=rev(colrb[10:90])){
    bin <- cut(x, seq(-3, 3, length.out=length(pal)), include.lowest=T) 
    pal[bin] }
if (color.edges){ E(net2)$color = col_fun(ewt) }
set.seed(1)
l <- layout_with_fr(net2, grid='nogrid') # Usually best

png('~/test.png', res=450, units='in', width=7, height=7)
par(mar=rep(0.1, 4))
plot(net2, layout=l, vertex.label.color="black", vertex.size=(4 - 3 * V(net2)$type)*1) 
dev.off()


# ------------------------------------
# Save mapping and centers attributes:
# ------------------------------------
save(mmat, full.mmat, clscols, clsmaxcols, mmap, full.mmap, motnam, full.motnam, file='motifs/motif_modules_network_data.Rda')
save(centers, enrichmat, rngroup, clsord, file='motifs/centers_enrichmat_data.Rda')



