#!/usr/bin/R
# library(infotheo)
# library(ape)
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
today <- format(Sys.time(), "%m%d%y")
library(ggplot2)
library(viridis)

# Prefix:
imgdir = paste0(img, "clusters/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'clusters_')

# TODO: get the number of elements in each cluster as well.
# --------------
# Load datasets:
# --------------
# Centers: 
centers = list()
centers[['DNase-seq']] = as.matrix(read.delim('dnase_300_centers.tsv',header=F))
centers[['Promoters']] = as.matrix(read.delim('n15_prom_300_centers.tsv', header=F))
centers[['Enhancers']] = as.matrix(read.delim('n15_enh_300_centers.tsv', header=F))
for (set in c('Enhancers','Promoters','DNase-seq')){
    colnames(centers[[set]]) = paste0('c',0:299)
}

counts = list()
counts[['DNase-seq']] = read.delim('DNase-seq_k300_counts.tsv',header=F)
counts[['Promoters']] = read.delim('n15_prom_k300_counts.tsv', header=F)
counts[['Enhancers']] = read.delim('n15_enh_k300_counts.tsv', header=F)
for (set in c('Enhancers','Promoters','DNase-seq')){
    ct = counts[[set]]$V2
    names(ct) = paste0('c',counts[[set]]$V1)
    counts[[set]] = ct
}

# Names:
n15n = scan('n15_names.tsv','c')
dn = scan('DNase_names.tsv', 'c')
rownames(centers[['DNase-seq']]) = dn
rownames(centers[['Promoters']]) = n15n
rownames(centers[['Enhancers']]) = n15n


# Annotation
map <- read.delim('Annotation/all_submitted_released_biosample_mapping_20180924.tsv',header=F)
names(map) <- c('ct','id')
rownames(map) <- map[,2]
# Order:
dford = read.delim('Annotation/bssid_order_20190106.tsv', header=F)
names(dford) <- c('ct','id', 'cls')
idorder = as.character(dford$id)
rn = as.character(dford$ct)
acut = dford$cls
dist.breaks = calc.breaks.acut(acut)
acut.nam = acut
names(acut.nam) = dford$id
# Availability:
mat = read.delim('ChromImpute/sample_mark_table.tsv', header=F)
names(mat) <- c('CellType','Epitope','file')
mat$file=1
wm <- spread(mat, CellType, file, fill=0)
rownames(wm) <- wm$Epitope 
colnames(wm) <- map[colnames(wm),1]
epitopes <- aggregate(file ~ Epitope, mat, length)
main.marks = as.character(head(epitopes$Epitope[order(epitopes$file, decreasing=T)], 9))

# -----------------------------------------------------
# Read in the colors and metadata for roadmap datasets:
# -----------------------------------------------------
rdmp = read.delim('Annotation/roadmap_metadata.tsv', header=T)
rdcol = unique(rdmp[,c('GROUP','COLOR')])
rdcol = rbind(rdcol, data.frame('GROUP'='NONE', 'COLOR' = 'white'))
rd.map = read.delim('Annotation/all_submitted_released_biosample_roadmap_mapping_20190106.tsv',header=F)
names(rd.map) = c('GROUP','ct','id')
rd.map = merge(rd.map, rdcol)
rownames(rd.map) = rd.map$id
rd.map$GROUP = factor(rd.map$GROUP, levels = as.character(rdcol$GROUP))
labels = rd.map[idorder, 'GROUP']
faclabels = as.matrix(as.numeric(labels))
colset = as.character(rdcol$COLOR)

# Label runs (that are not NONE):
rl = rle(as.numeric(faclabels))
cs = cumsum(rl$lengths)
cb = c(0, cs[-length(cs)] + 1)
loc = (cs + cb )/ 2 / length(labels)
idx = which(rl$lengths > 2)
lab = levels(labels)[rl$values[idx]]

# Remove NONE and average consecutive values:
kid = which(lab != 'NONE')
col.lab = as.character(rdcol$COLOR[rl$values[idx[kid]]])
lab = lab[kid]
loci = loc[idx[kid]]
rl = rle(lab)
cs = cumsum(rl$lengths)
cb = c(0, cs[-length(cs)] + 1)
locfinal = c()
lab = rl$values
col.lab = col.lab[cs]
locfinal = sapply(1:length(cs), function(i){mean(loci[cb[i]:cs[i]])})

# ---------------------------
# Plot clusters side by side:
# ---------------------------
plot.centers = function(centers, set, idorder, counts, 
                        subset=FALSE, cls=acut.nam, cls.ord=NULL){
    mat = centers[[set]]
    ct = counts[[set]]
    if (subset) {
        cidx = keep.cls[[set]]
        mat = mat[,cidx]
        ct = ct[cidx]
    }
    idx = as.numeric(which(ct != 0))
    mat = mat[,idx]
    if (is.null(cls.ord)){
        subid = idorder[idorder %in% rownames(mat)]
        tmp = matrix(NA, nrow=length(idorder), ncol=ncol(mat))
        rownames(tmp) = idorder
        # Reorder:
        ll = diag.mat(mat[subid,])
        tmp[subid,]= ll[[1]]
        # Use cto to get breaks:
        cto = ll[[3]]
        vcut =  c(cto[cto ==0], cls[subid][cto])
        vbreaks = calc.breaks.acut(vcut)
    } else {
        subid = idorder[idorder %in% rownames(mat)]
        subclsord = cls.ord[cls.ord %in% colnames(mat)]
        tmp = matrix(NA, nrow=length(idorder), ncol=length(subclsord),
                     dimnames=list(idorder, subclsord))
        # Reorder:
        tmp[subid,subclsord] = mat[subid,subclsord]
        ll = list(subclsord, subclsord)
        vcut = round(1:ncol(tmp) / 20) 
        vbreaks = calc.breaks.acut(vcut)
    }
    image(t(tmp), axes=F, col=col, zlim=c(0,1))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
    mtext(set, side=3, cex=1.3)
    return(list(ll[[2]], vbreaks))
}

plot.counts = function(counts, set, ordll, ylim=NULL, subset=FALSE){
    ct = counts[[set]]
    ord = ordll[[1]]
    vbreaks = ordll[[2]]
    # Using 
    if (subset) { ct = ct[keep.cls[[set]]] }
    ct = ct[ct != 0]
    subord = ord[ord %in% names(ct)]
    ct = ct[subord]
    barplot(-ct, ylim=ylim, col='darkgrey', border=NA, yaxt='n',xaxt='n')
    abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
}

plot.motifs <-  function(motifs, set, ordll, zlim=0.5, motiford=NULL){
    ct = counts[[set]]
    id = names(which(ct[keep.cls[[set]]] != 0))
    mot = motifs[[set]]
    mot = mot[id,]
    ord = ordll[[1]]
    vbreaks = ordll[[2]]
    subord = ord[ord %in% rownames(mot)]
    mot = mot[subord,]
    # Threshold for plotting:
    mot[mot < -zlim] <- -zlim
    mot[mot > zlim] <- zlim
    if (!is.null(motiford)){ mot = mot[,motiford] }
    image(mot, axes=FALSE,col=rev(colrb), zlim=c(-zlim, zlim))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
}

# TODO: ALSO ORDER CELLS BY THEIR ELEMENTS?
# TODO: REDUCE CENTERS.
# TODO: Breaks by hierarchical clusters of distmat?

png(paste0(imgpref,'cluster_centers_orddist.png'),res=300,units='in',width=15,height=10)
layout(matrix(c(1:10),2,5), heights=c(8,1), widths=c(1,4,4,4,1), TRUE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=locfinal,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lab, srt=0, adj=1, xpd=TRUE,cex=.8, col=col.lab)
plot.new()
# Plot clusters and counts
for (set in c('Enhancers','Promoters','DNase-seq')){
    par(mar=c(0.25, 0.25, 2, 0.25))
    dist.order[[set]] = plot.centers(centers, set, idorder, counts=counts)
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, dist.order[[set]])
}
# Availability:
par(mar=c(.25, 0.25, 2, 0.25))
avail = as.matrix(wm[main.marks, rn])
image(avail, axes=F, col=c('white', 'darkgrey'))
par(mar=c(.25, 0.25, 0, 0.25))
image(avail, axes=F, col='white')
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[4]-0.01*(par()$usr[4]-par()$usr[3]), 
     labels=main.marks, srt=90, adj=1, xpd=TRUE,cex=.6)
dev.off()

png(paste0(imgpref,'cluster_centers_orddist_nodnase.png'),res=300,units='in',width=12,height=10)
layout(matrix(c(1:8),2,4), heights=c(8,1), widths=c(1,4,4,1), TRUE)
par(yaxs="i")
par(xaxs="i")
dist.order = list()
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=locfinal,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lab, srt=0, adj=1, xpd=TRUE,cex=.8, col=col.lab)
plot.new()
# Plot clusters and counts
for (set in c('Enhancers','Promoters')){
    par(mar=c(0.25, 0.25, 2, 0.25))
    dist.order[[set]] = plot.centers(centers, set, idorder, counts=counts)
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, dist.order[[set]])
}
# Availability:
par(mar=c(.25, 0.25, 2, 0.25))
avail = as.matrix(wm[main.marks, rn])
image(avail, axes=F, col=c('white', 'darkgrey'))
par(mar=c(.25, 0.25, 0, 0.25))
image(avail, axes=F, col='white')
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[4]-0.01*(par()$usr[4]-par()$usr[3]), 
     labels=main.marks, srt=90, adj=1, xpd=TRUE,cex=.6)
dev.off()


# Sort first by rdmp colors then by clusters
ord2 = order(labels, decreasing=TRUE)
rn2 = rn[ord2]
idorder2 = idorder[ord2]
labels2 = rd.map[idorder2, 'GROUP']
colcut = as.numeric(labels2)
col.breaks = calc.breaks.acut(colcut)
faclabels2 = as.matrix(colcut)
colcut.nam = colcut
names(colcut.nam) = idorder2

# Loc labels. Adjust vals if very close to each other.
ldf = data.frame(label=labels2, quant=1:length(labels2) * 1.0/ length(labels2) )
labeldf = aggregate(quant ~ label, ldf, mean)

# TODO: compute NEW BREAKS:
# acut <- cutree(ht, nclust)[cocl]
# cuts <- cumsum(rle(acut)$lengths)
# step <- head(diff(seq(0,1,length.out=length(cocl))),1)
# cuts <-  cuts[-length(cuts)] -1
# cuts * step + step/2

png(paste0(imgpref,'cluster_centers_BYrdcol.png'),res=300,units='in',width=15,height=10)
layout(matrix(c(1:10),2,5), heights=c(8,1), widths=c(1,4,4,4,1), TRUE)
par(yaxs="i")
par(xaxs="i")
col.order = list()
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels2), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=labeldf$quant,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=labeldf$label, srt=0, adj=1, xpd=TRUE,cex=.8, col=colset)
plot.new()
# Plot clusters and counts
for (set in c('Enhancers','Promoters','DNase-seq')){
    par(mar=c(0.25, 0.25, 2, 0.25))
    col.order[[set]] = plot.centers(centers, set, idorder2, 
                                    counts=counts, cls=colcut.nam)
    abline(h=col.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, col.order[[set]])
}
# Plot data availability
par(mar=c(.25, 0.25, 2, 0.25))
avail = as.matrix(wm[main.marks, rn2])
image(avail, axes=F, col=c('white', 'darkgrey'))
par(mar=c(.25, 0.25, 0, 0.25))
image(avail, axes=F, col='white')
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[4]-0.01*(par()$usr[4]-par()$usr[3]), 
     labels=main.marks, srt=90, adj=1, xpd=TRUE,cex=.6)
dev.off()


# ---------------------------------
# Add motif enrichment information:
# ---------------------------------
# Load in enrichment files:
enrwide = c()
for (filename in c('dnase_300_enrichments.tsv', 'n15_enh_300_enrichments.tsv','n15_prom_300_enrichments.tsv')){
    enrdf = read.delim(filename, header=T)
    print(filename)
    print(length(unique(enrdf$Prefix)))
    enrdf = enrdf[enrdf$compare == '+_.', c('Motif', 'Prefix', 'log2_enrich')]
    wide = spread(enrdf, Motif, log2_enrich)
    mat = as.matrix(wide[,-1])
    rownames(mat) = wide$Prefix
    enrwide = rbind(enrwide, mat)
}
prefixes = c('DNase-seq' = 'DNase-seq_k300_',
             'Enhancers' = 'n15_enh_k300_', 
             'Promoters' = 'n15_prom_k300_')

# Motif information:
motfam = read.delim('Annotation/motifs-clust-names.txt', header=F, stringsAsFactors=F)
motmap = c()
for (i in 1:nrow(motfam)){
    members = strsplit(motfam[i,1], ";")[[1]]
    motmap = rbind(motmap, data.frame(motif = members, num=i,
                                      main = motfam[i,2], fam = motfam[i,3]))
}
rownames(motmap) = motmap$motif

# Select only one motif per family:
# sum(colnames(enrwide) %in% motmap$motif) / ncol(enrwide) * 100
locmax = apply(abs(enrwide),2, max)
maxdf = data.frame(topval=locmax, motif=names(locmax))
maxdf = merge(motmap, maxdf)
keep.motifs = unlist(sapply(unique(motmap$num), function(i){
                         subdf = maxdf[maxdf$num == i,]
                         as.character(subdf[which.max(subdf$topval),'motif']) }))
keep.motifs = keep.motifs[locmax[keep.motifs] >= 1.5]
redenr = enrwide[,keep.motifs]
print(dim(redenr))

# Motif order:
motifs = list()
keep.cls = list()
for (set in c('Enhancers','Promoters','DNase-seq')){
    pref = prefixes[set]
    idx = grep(pref, rownames(redenr))
    num = as.numeric(sub(pref, "", rownames(redenr)[idx]))
    motifs[[set]] = redenr[idx[order(num)],]
    rownames(motifs[[set]]) = paste0('c',num)
    keep.cls[[set]] = rownames(motifs[[set]])
}

# Get color sorted and dist sorted orders:
# By diagonalizingthe motifs matrix.
main = c()
# for (set in c('Enhancers','Promoters','DNase-seq')){
for (set in c('Enhancers')){
    mat = motifs[[set]]
    ord = col.order[[set]][[1]]
    ord = ord[ord %in% rownames(mat)]
    main = rbind(main,mat[ord,])
}
col.motiford = diag.mat(main)[[2]]

# Diagonalize the motifs matrix:
main = c()
# for (set in c('Enhancers','Promoters','DNase-seq')){
for (set in c('Enhancers')){
    mat = motifs[[set]]
    ord = dist.order[[set]][[1]]
    ord = ord[ord %in% rownames(mat)]
    main = rbind(main,mat[ord,])
}
# dist.motiford = diag.mat(main)[[2]]
dist.motiford = rownames(reord(t(main)))

# Ratios:
widths = unlist(lapply(motifs, nrow)) / nrow(redenr) * 10

png(paste0(imgpref,'centers_motifs_orddist.png'),res=300,units='in',width=15,height=12)
layout(matrix(c(1:12),3,4), heights=c(6,3,1), widths=c(1,widths), TRUE)
# layout(matrix(c(1:12),3,4), heights=c(6,3,1), widths=c(1.5, widths), TRUE)
par(yaxs="i")
par(xaxs="i")
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=locfinal,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lab, srt=0, adj=1, xpd=TRUE,cex=.8, col=col.lab)
par(mar=c(0.25, 0.25, 0, 0.25))
image(motifs[['Enhancers']], col='white', axes=F)
text(y=seq(0,1,length.out=length(dist.motiford)), 
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]),
     labels=dist.motiford, srt=0, adj=1, xpd=TRUE,cex=.5)
#
plot.new()
for (set in names(motifs)){
    par(mar=c(0.25, 0.25, 2, 0.25))
    ord = plot.centers(centers, set, idorder, counts=counts,
                       subset=TRUE, cls=acut.nam)
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.motifs(motifs, set, ord, motiford=dist.motiford)
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, ord, subset=TRUE)
}
dev.off()

png(paste0(imgpref,'centers_motifs_orddist_nodnase.png'),res=300,units='in',width=15,height=12)
layout(matrix(c(1:12),3,3), heights=c(6,3,1), widths=c(1,widths[1:2]), TRUE)
# layout(matrix(c(1:12),3,4), heights=c(6,3,1), widths=c(1.5, widths), TRUE)
par(yaxs="i")
par(xaxs="i")
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=locfinal,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lab, srt=0, adj=1, xpd=TRUE,cex=.8, col=col.lab)
par(mar=c(0.25, 0.25, 0, 0.25))
image(motifs[['Enhancers']], col='white', axes=F)
text(y=seq(0,1,length.out=length(dist.motiford)), 
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]),
     labels=dist.motiford, srt=0, adj=1, xpd=TRUE,cex=.5)
#
plot.new()
for (set in names(motifs)[1:2]){
    par(mar=c(0.25, 0.25, 2, 0.25))
    ord = plot.centers(centers, set, idorder, counts=counts,
                       subset=TRUE, cls=acut.nam)
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.motifs(motifs, set, ord, motiford=dist.motiford)
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, ord, subset=TRUE)
}
dev.off()



# By color
png(paste0(imgpref,'centers_motifs_BYrdcol.png'),res=300,units='in',width=15,height=10)
layout(matrix(c(1:12),3,4), heights=c(6,4,1), widths=c(1,widths), TRUE)
par(yaxs="i")
par(xaxs="i")
# Empty square:
col.order = list()
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels2), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=labeldf$quant,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=labeldf$label, srt=0, adj=1, xpd=TRUE,cex=.8, col=colset)
par(mar=c(0.25, 0.25, 0, 0.25))
image(motifs[['Enhancers']], col='white', axes=F)
text(y=seq(0,1,length.out=length(dist.motiford)), 
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]),
     labels=dist.motiford, srt=0, adj=1, xpd=TRUE,cex=.25)
#
plot.new()
for (set in names(motifs)){
    par(mar=c(0.25, 0.25, 2, 0.25))
    col.order[[set]] = plot.centers(centers, set, idorder2, counts=counts, 
                                    subset=TRUE, cls=colcut.nam)
    abline(h=col.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.motifs(motifs, set, col.order[[set]], motiford=col.motiford)
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, col.order[[set]], subset=TRUE)
}
dev.off()



# ---------------------
# Add the GWAS results:
# ---------------------
gwasdf = c()
gwasfiles = c('DNase-seq' = 'DNase-seq_k300_enrich.txt',
              'Promoters' = 'n15_prom_k300_enrich.txt',
              'Enhancers' = 'n15_enh_k300_enrich.txt')
for (set in c('Enhancers','Promoters','DNase-seq')){
    filename = gwasfiles[[set]]
    df = read.delim(filename, header=F)
    names(df) = c('pvalue','cluster','pmid','trait','counthit','countall','fold')
    df$set = set
    df$pmt = paste0(df$pmid, '_', df$trait)
    df$cls = paste0(set, '_', df$cluster)
    gwasdf = rbind(gwasdf, df)
}
gwasdf$logpval = -log10(gwasdf$pvalue)
gwaslong = aggregate(logpval ~ cls + pmt, gwasdf, max)
wide = spread(gwaslong, pmt, logpval, fill=0)
gwasmat = as.matrix(wide[,-1])
rownames(gwasmat) = wide$cls
fullgwas = gwasmat
gwasmat[gwasmat < 1] = 0
# aggregate(logpval ~ set, gwasdf[gwasdf$logpval > 2,], length)

# TODO: either difference or totals.

# Put all matrices in list:
gwas = list()
allgwas = list()
keep.cls.gwas = list()
setlist = c('Enhancers','Promoters','DNase-seq')
for (set in setlist){
    idx = grep(set, rownames(gwasmat))
    mat = gwasmat[idx,]
    rownames(mat) = sub(paste0(set,"_"), "c",rownames(mat))
    gwas[[set]] = mat
    mat = fullgwas[idx,]
    rownames(mat) = sub(paste0(set,"_"), "c",rownames(mat))
    allgwas[[set]] = mat
    keep.cls.gwas[[set]] = rownames(gwas[[set]])
}

SHOWGWAS=100
gwasmarg = sort(apply(gwasmat, 2, sum), decreasing=T)
keep.studies = names(head(gwasmarg, SHOWGWAS))
# Order the top studies:
r2 = reord(t(gwasmat[, keep.studies]), 'cosine')
studyord = rownames(r2)
r2 = reord(gwasmat[, studyord], 'Jaccard')
r3 = gwasmat[, studyord]
r3[r3 < .5] = 0
r3 = reord(r3, 'cosine')
# image(r3, col=colred)
# Rownames:
rngwas = rownames(r3)
ord.gwas = list()
setlist = c('Enhancers','Promoters','DNase-seq')
for (set in setlist){
    idx = grep(set, rngwas)
    ord.gwas[[set]] = sub(paste0(set,"_"), "c",rngwas[idx])
}

plot.gwas <-  function(gwas, set, ordll, zmin=0.5, zmax=2, gwasord=NULL){
    ct = counts[[set]]
    id = names(which(ct[keep.cls.gwas[[set]]] != 0))
    gmat = gwas[[set]]
    gmat = gmat[id,]
    print(dim(gmat))
    ord = ordll[[1]]
    vbreaks = ordll[[2]]
    subord = ord[ord %in% rownames(gmat)]
    gmat = gmat[subord,]
    print(dim(gmat))
    # Threshold for plotting:
    gmat[gmat > zmax] <- zmax
    gmat[gmat < zmin] <- 0
    if (!is.null(gwasord)){ gmat = gmat[,gwasord] }
    image(gmat, axes=FALSE,col=colred, zlim=c(zmin, zmax))
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    abline(v=vbreaks,lty='dotted',lw=1, col='darkgrey')
}


png(paste0(imgpref,'centers_gwas_orddist.png'),res=450,units='in',width=15,height=10)
layout(matrix(c(1:12),3,4), heights=c(4,4,1), widths=c(1.5,4,4,4), TRUE)
par(yaxs="i")
par(xaxs="i")
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=locfinal,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lab, srt=0, adj=1, xpd=TRUE,cex=.8, col=col.lab)
par(mar=c(0.25, 0.25, 0.1, 0.25))
# mtext('GWAS Enrichments', side=2, line=.5)
image(gwas[['Enhancers']], col='white', axes=F)
text(y=seq(par()$usr[3],par()$usr[4], length.out=length(studyord)),
     x=par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),
     labels=studyord, srt=0, adj=1, xpd=TRUE,cex=.35)
##
plot.new()
for (set in setlist){
    par(mar=c(0.25, 0.25, 2, 0.25))
    ord = plot.centers(centers, set, idorder, counts=counts,
                       subset=FALSE, cls=acut.nam, cls.ord=ord.gwas[[set]])
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0.1, 0.25))
    plot.gwas(gwas, set, ord, gwasord=studyord)
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, ord, subset=FALSE)
}
dev.off()

png(paste0(imgpref,'centers_gwas_orddist_nodnase.png'),res=450,units='in',width=11,height=10)
layout(matrix(c(1:9),3,3), heights=c(4,4,1), widths=c(1,4,4), TRUE)
par(yaxs="i")
par(xaxs="i")
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=locfinal,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lab, srt=0, adj=1, xpd=TRUE,cex=.8, col=col.lab)
par(mar=c(0.25, 0.25, 0.1, 0.25))
image(gwas[['Enhancers']], col='white', axes=F)
text(y=seq(par()$usr[3],par()$usr[4], length.out=length(studyord)),
     x=par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),
     labels=studyord, srt=0, adj=1, xpd=TRUE,cex=.35)
mtext('GWAS Enrichments', side=2, line=.5)
##
plot.new()
for (set in setlist[1:2]){
    par(mar=c(0.25, 0.25, 2, 0.25))
    ord = plot.centers(centers, set, idorder, counts=counts,
                       subset=FALSE, cls=acut.nam, cls.ord=ord.gwas[[set]])
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0.1, 0.25))
    plot.gwas(gwas, set, ord, gwasord=studyord)
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, ord, subset=FALSE)
}
dev.off()




png(paste0(imgpref,'centers_gwas_BYrdcol.png'),res=450,units='in',width=15,height=10)
layout(matrix(c(1:12),3,4), heights=c(4,4,1), widths=c(1.5,4,4,4), TRUE)
par(yaxs="i")
par(xaxs="i")
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels2), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=labeldf$quant,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=labeldf$label, srt=0, adj=1, xpd=TRUE,cex=.8, col=colset)
par(mar=c(0.25, 0.25, 0.1, 0.25))
mtext('GWAS Enrichments', side=2, line=.5)
image(gwas[['Enhancers']], col='white', axes=F)
text(y=seq(par()$usr[3],par()$usr[4], length.out=length(studyord)),
     x=par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),
     labels=studyord, srt=0, adj=1, xpd=TRUE,cex=.35)
##
plot.new()
for (set in setlist){
    par(mar=c(0.25, 0.25, 2, 0.25))
    ord = plot.centers(centers, set, idorder2, counts=counts,
                       subset=FALSE, cls=acut.nam, cls.ord=ord.gwas[[set]])
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0.1, 0.25))
    plot.gwas(gwas, set, ord, gwasord=studyord)
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, ord, subset=FALSE)
}
dev.off()


# Enrichment curves + difference:
# Put all matrices in list:
diffgwas = list()
setlist = c('Enhancers','Promoters','DNase-seq')
dnmat = allgwas[['DNase-seq']]
emat = allgwas[['Enhancers']]
ev = apply(emat,2,sum) 
dev = apply(dnmat,2,sum)

madf = data.frame(sum=ev + dev,
                  diff=ev - dev,
                  label=names(ev))
idx = which(((madf$sum / abs(madf$diff)) < 5) & (madf$sum > 50) | madf$sum > 100)

library(ggrepel)
labdf = data.frame(x=c(75,75), y=c(-7,7), label=c('DNase-seq','Enhancers'))
gplot = ggplot(madf, aes(sum, diff)) + geom_point(alpha=0.25) + 
    geom_hline(yintercept=0, linetype='dotted') + 
    geom_text_repel(data=madf[idx,], aes(sum, diff, label=label), size=2.5)  + 
    geom_text(data=labdf, aes(x,y,label=label), color='darkblue',alpha=0.75, size=10) + 
    labs(x='Sum of Enhancer + DNase GWAS Enrichments across modules',
         y='Difference of Enhancer - DNase GWAS Enrichments across modules') + 
    theme_minimal()
ggsave(paste0(imgpref, 'gwas_MAplot.png'), gplot, dpi=300, width=11, height=8,units='in')


# Plot most specific
idx = which(((madf$sum / abs(madf$diff)) < 3 & madf$diff > 0 & madf$sum > 35) | (madf$sum > 50 & madf$diff > 2))
topenh = as.character(madf$label[idx])



png(paste0(imgpref,'centers_gwas_orddist_topenh.png'),res=450,units='in',width=15,height=8)
layout(matrix(c(1:12),3,4), heights=c(3,3,1), widths=c(1.5,4,4,4), TRUE)
par(yaxs="i")
par(xaxs="i")
# Empty square:
par(mar=c(0.5, 6, 2, 0.15))
image(t(faclabels), axes=F, col=colset)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=locfinal,
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lab, srt=0, adj=1, xpd=TRUE,cex=.8, col=col.lab)
par(mar=c(0.25, 0.25, 0.1, 0.25))
mtext('GWAS Enrichments', side=2, line=.5)
image(gwas[['Enhancers']], col='white', axes=F)
text(y=seq(par()$usr[3],par()$usr[4], length.out=length(keep.studies)),   x=par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),
     labels=keep.studies, srt=0, adj=1, xpd=TRUE,cex=.15)
##
plot.new()
for (set in setlist){
    par(mar=c(0.25, 0.25, 2, 0.25))
    ord = plot.centers(centers, set, idorder, counts=counts,
                       subset=FALSE, cls=acut.nam)
    abline(h=dist.breaks,lty='dotted', lw=1, col='darkgrey')
    par(mar=c(0.25, 0.25, 0.1, 0.25))
    plot.gwas(gwas, set, ord, gwasord=studyord)
    par(mar=c(0.25, 0.25, 0, 0.25))
    plot.counts(counts, set, ord, subset=FALSE)
}
dev.off()







plot(ev + dev, ev - dev, 
     pch=19, cex=.5, col=rgb(0,0,0,.25), 
     xlab='Sum of Enhancer + DNase enrichments',
     ylab='Enhancer - DNase enrichments')
text(ev[idx] + dev[idx], ev[idx] - dev[idx], 
     labels=names(ev)[idx],
     adj=c(.5,1), pch=19, cex=.5, col=rgb(0,0,0,.75))


abline(h=0)



for (set in setlist[1:2]){
    idx = grep(set, rownames(gwasmat))

    mat = gwasmat[idx,]
    rownames(mat) = sub(paste0(set,"_"), "c",rownames(mat))
    gwas[[set]] = mat
    keep.cls.gwas[[set]] = rownames(gwas[[set]])
}

SHOWGWAS=200
gwasmarg = sort(apply(gwasmat, 2, max), decreasing=T)
keep.studies = names(head(gwasmarg, SHOWGWAS))




# NOTE: TOP 200 only
png(paste0(imgpref,'only_gwas_orddist.png'),res=300,units='in',width=13,height=12)

layout(matrix(c(1:4),1,4), widths=c(1.5,4,4,4), TRUE)
par(yaxs="i")
par(xaxs="i")
# GWAS LABELS:
par(mar=c(0.25, 5, 2, 0.25))
image(gwas[['Enhancers']], col='white', axes=F)
text(y=seq(0,1,length.out=length(keep.studies)), 
     x=par()$usr[2]-0.001*(par()$usr[2]-par()$usr[1]),
     labels=keep.studies, srt=0, adj=1, xpd=TRUE,cex=.5)
for (set in setlist){
    par(mar=c(0.25, 0.25, 2, 0.25))
    plot.gwas(gwas, set, ord, gwasord=studyord)
}

dev.off()


# DIFFERENTIAL
# ALL 
# SHOW MORE TISSUE SPECIFIC ENRICHMENTS than DNase?
