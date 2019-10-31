#!/usr/bin/R
# ---------------------------------------
# Plot the imputed, observed, and 
# mix of imputed/observed matrices
# Plot MDS
# Validate imputed diff
# Look for consistency in clades
# Plot dendrograms colored by tissue type
# Identify consistent breakpoints
# ---------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
library(ggplot2)
library(viridis)
library(ggrepel)
library(dendextend)
library(uwot)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
options(scipen=45) # So we dont get issues writing integers into bedfiles

# Load specific distance matrices:
fixedstate = FALSE
if (fixedstate){
    source(paste0(bindir, 'load_region_distance_matrices.R'))
    setprefix = paste0('region_',nregions,'_distances_')
} else { 
    source(paste0(bindir, 'load_distance_matrices.R'))
    setprefix = 'distances_'
}

today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "imp_distance/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, setprefix)

# --------------------------------------
# Plot dendrograms with colored branches
# And extended metadata:
# --------------------------------------
# Reorder matrix:
rn = rownames(full)
keep.ct = rn[!is.na(meta[rn, 'infoline'])]
subfull = full[keep.ct, keep.ct]
dt = as.dist(subfull)
# Create the dend:
dend <- hclust(dt, method = "ward.D") %>% as.dendrogram
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

# Save dendrogram, order, and metadata:
save(dend3, file=paste0(setprefix, 'dendrogram_wardD_full.rda'))
write.table(cbind(id=lab, metamat[lab,]), file=paste0(setprefix, 'tree_ordered_metadata.tsv'), 
            row.names=F, col.names=T, quote=F, sep="\t")

# Make figures of dendrogram:
pdf(paste0(imgpref,'full_tree_metadata.pdf'),width=16,height=3)
layout(matrix(c(1,2,3), 3,1), heights=c(8,1,1), TRUE)
par(yaxs="i")
par(xaxs="i")
par(mar=c(2.5,6,3,3))
plot(dend3, horiz=FALSE)
par(mar=c(0.05,6,.25,3))
meta.image(metamat[lab,], colvals=colvals, cex=.6)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
par(mar=c(0.25,6,.05,3))
avail = t(as.matrix(wm[main.marks, lab]))
image(avail, yaxt='n',xaxt='n', col=c('white','darkgrey'))
text(y=seq(0,1, length.out=length(main.marks)),
     x=par()$usr[1]-0.001*(par()$usr[2]-par()$usr[1]), 
     labels=main.marks, srt=0, adj=1, xpd=true,cex=.25)
dev.off()

pdf(paste0(imgpref,'full_tree_metadata_horiz.pdf'),width=3,height=16)
layout(matrix(c(1,2,3), 1,3), widths=c(8,1,1), TRUE)
par(mar=c(4,1,1,2.5))
par(yaxs="i")
par(xaxs="i")
plot(dend3, horiz=TRUE)
par(mar=c(4,.25,1,.05))
meta.image(metamat[lab,5:1], colvals=colvals, cex=.6, horiz=TRUE)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
par(mar=c(4,.05,1,.25))
avail = as.matrix(wm[main.marks, lab])
image(avail, yaxt='n',xaxt='n', col=c('white','darkgrey'))
text(x=seq(0,1, length.out=length(main.marks)),
     y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
     labels=main.marks, srt=90, adj=1, xpd=true,cex=.25)
dev.off()

# --------------------
# Make circular plots:
# --------------------
NLAB = length(labels(dend3))
mll = meta.image(metamat[lab,5:1], colvals=colvals, cex=.6, horiz=TRUE, return.mat=TRUE)
nummat = mll[[1]]
colsmeta = mll[[2]]
colmat = matrix(colsmeta[nummat], nrow=nrow(nummat))
labels_dend <- labels(dend3)
if (as.logical(anyDuplicated(labels_dend))) {
    labels(dend3) <- paste0(seq_along(labels_dend), "_", labels_dend)
    labels_dend <- labels(dend3)
}

# ---------------------------
# Full with names and matrix:
# ---------------------------
pdf(paste0(imgpref,'full_tree_metadata_circ.pdf'),width=12,height=12)
circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
circos.initialize(factors = "a", xlim = c(0, NLAB)) 
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
                 circos.text(1:NLAB-0.5, rep(0, NLAB), labels(dend3), col = labels_colors(dend3), cex=.25,
                             facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) }, bg.border = NA, track.height = 0.1)
circos.track(ylim = c(0, 5), bg.border = NA, panel.fun = function(x, y) {
    m = t(colmat[,5:1])
    nr = nrow(m)
    nc = ncol(m)
    for(i in 1:nr) {
        circos.rect(1:nc - 1, rep(nr - i, nc), 
            1:nc, rep(nr - i + 1, nc), 
            border = NA, col = m[i, ])
    }
}, track.height=0.125) 
max_height = max(attr(dend3, "height"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                 circos.dendrogram(dend3, max_height = max_height) }, track.height = 0.5, bg.border = NA)
circos.clear()
dev.off()

# --------------------------------------
# Plot dendrogram, split by NCLUST cuts:
# --------------------------------------
NCLUST = 20
colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(NCLUST)
dend3 <- color_branches(dend2, k=NCLUST, col=colpair)
dend3 = set(dend3, "labels_cex", .18)
labels_dend <- labels(dend3)
if (as.logical(anyDuplicated(labels_dend))) {
    labels(dend3) <- paste0(seq_along(labels_dend), "_", labels_dend)
    labels_dend <- labels(dend3)
}
# Factors/splits:
kcol = unlist(get_leaves_attr(dend3, 'edgePar'))
letext = c(letters, paste0('a', letters))
factors = letext[as.numeric(factor(kcol, levels = unique(kcol)))]
dend_list = get_subdendrograms(dend3, k=NCLUST)
# Correct the dendrogram order:
dnames = sapply(dend_list, function(x){as.numeric(strsplit(labels(x)[1],'_')[[1]][[1]]) })
names(dend_list) = letext[1:NCLUST][order(order(dnames))]


pdf(paste0(imgpref,'full_tree_metadata_circ_splitcls_',NCLUST,'.pdf'),width=12,height=12)
circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
circos.initialize(factors, xlim = cbind(c(0, 0), table(factors)))
# Text labels:
circos.track(ylim = c(0, 1), 
             panel.fun = function(x, y) {
                 sector.index = CELL_META$sector.index
                 dend = dend_list[[sector.index]]
                 txt = labels(dend)
                 nc = length(txt)
                 circos.text(1:nc-0.5, rep(0, nc), txt, 
                             col = labels_colors(dend), cex=.25,
                             facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) 
             }, bg.border = NA, track.height = 0.1)
# Metadata:
circos.track(ylim = c(0, 5), bg.border = NA, panel.fun = function(x, y) {
    sector.index = CELL_META$sector.index
    m = t(colmat[factors == sector.index,5:1])
    nr = nrow(m)
    nc = ncol(m)
    for(i in 1:nr) {
        circos.rect(1:nc - 1, rep(nr - i, nc), 
            1:nc, rep(nr - i + 1, nc), 
            border = NA, col = m[i, ])
    }
}, track.height=0.125)
# Dendrogram
max_height = max(sapply(dend_list, function(x) attr(x, "height")))
circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.3, 
    panel.fun = function(x, y) {
        sector.index = CELL_META$sector.index
        dend = dend_list[[sector.index]]
        circos.dendrogram(dend, max_height = max_height)
})
dev.off()



# ------------------------------------------
# Identify, plot, and write out breakpoints:
# ------------------------------------------
nl = get_nodes_attr(dend, 'members')
ind = which(nl >= 10)
stlab = partition_leaves(dend)
mlist = sapply(ind,function(i){ meta[stlab[[i]],'GROUP'] })
ilist = sapply(ind,function(i){ meta[stlab[[i]],'infoline'] })
intern = sapply(mlist, function(x){ tab = table(as.character(x)); max(tab)/sum(tab) })
len = lapply(mlist, length)
# plot(1:length(len), len)
# plot(intern, len)
possible = ind[which(intern > 0.7)]
large = ind[which(unlist(len) >= 7)]
large = large[!(large %in% possible)]
pos2 = unique(c(possible, large))
NCLUST = 20
colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(NCLUST)
dend3 <- color_branches(dend2, k=NCLUST, col=colpair)
dend3 = set(dend3, "labels_cex", .18)
labels_dend <- labels(dend3)
if (as.logical(anyDuplicated(labels_dend))) {
    labels(dend3) <- paste0(seq_along(labels_dend), "_", labels_dend)
    labels_dend <- labels(dend3)
}
dend3.save = dend3
nodecex = rep(1, nnodes(dend3))
nodepch = rep(NA, nnodes(dend3))
nodecex[pos2] = 1
# nodepch[possible] = 18
nodepch[possible] = 'o'
nodepch[large] = 'x'
# nodepch[possible] = as.character(1:length(possible))
# nodepch[possible] = c(letters, letters, letters[1:(length(possible)-2*26)] )
dend3 = set(dend3, "nodes_pch", nodepch)
dend3 = set(dend3, "nodes_cex", nodecex) 
nodedf = data.frame(cbind(get_nodes_xy(dend3, type = 'rectangle')[pos2,], node=pos2))

# Factors/splits:
kcol = unlist(get_leaves_attr(dend3, 'edgePar'))
letext = c(letters, paste0('a', letters))
factors = letext[as.numeric(factor(kcol, levels = unique(kcol)))]
dend_list = get_subdendrograms(dend3, k=NCLUST)
# Correct the dendrogram order:
dnames = sapply(dend_list, function(x){as.numeric(strsplit(labels(x)[1],'_')[[1]][[1]]) })
names(dend_list) = letext[1:NCLUST][order(order(dnames))]

pdf(paste0(imgpref,'full_tree_metadata_circ_possible_',NCLUST,'.pdf'),width=12,height=12)
circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1))
circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
circos.initialize(factors, xlim = cbind(c(0, 0), table(factors)))
# Text labels:
circos.track(ylim = c(0, 1), 
             panel.fun = function(x, y) {
                 sector.index = CELL_META$sector.index
                 dend = dend_list[[sector.index]]
                 txt = labels(dend)
                 nc = length(txt)
                 circos.text(1:nc-0.5, rep(0, nc), txt, 
                             col = labels_colors(dend), cex=.25,
                             facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) 
             }, bg.border = NA, track.height = 0.1)
# Metadata:
circos.track(ylim = c(0, 5), bg.border = NA, panel.fun = function(x, y) {
    sector.index = CELL_META$sector.index
    m = t(colmat[factors == sector.index,5:1])
    nr = nrow(m)
    nc = ncol(m)
    for(i in 1:nr) {
        circos.rect(1:nc - 1, rep(nr - i, nc), 
            1:nc, rep(nr - i + 1, nc), 
            border = NA, col = m[i, ])
    }
}, track.height=0.125)
# Dendrogram
max_height = max(sapply(dend_list, function(x) attr(x, "height")))
circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.3, 
    panel.fun = function(x, y) {
        sector.index = CELL_META$sector.index
        dend = dend_list[[sector.index]]
        circos.dendrogram(dend, max_height = max_height)
})
circos.clear()
# Inner dendrogram
par(new = TRUE) # <- allows new fig inside.
max_height = max(attr(dend, "height"))
circos.par("canvas.xlim" = c(-2.25, 2.25), "canvas.ylim" = c(-2.25, 2.25))
circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
circos.initialize(factors = 'single', xlim = c(0, NLAB))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                 circos.dendrogram(dend3.save, max_height = max_height) 
                 circos.text(nodedf$V1, max_height - nodedf$V2, as.character(nodedf$node),
                             col = 'black', cex=.25,
                             facing = "reverse.clockwise", niceFacing=TRUE, adj = c(0, 0.5)) 
}, track.height = 0.5, bg.border = NA)
circos.clear()
dev.off()


# --------------------------------------------
# Plot breakpoints chosen with the full tree: 
# --------------------------------------------
bkdf = read.delim('Annotation/tissue_bkpts.tsv',header=T, stringsAsFactors=F)
dend3 = dend3.save
nodecex = rep(1, nnodes(dend3))
nodepch = rep(NA, nnodes(dend3))
dend3 = set(dend3, "nodes_pch", nodepch) 
nodecex[bkdf$Node] = 1
nodepch[bkdf$Node] = 'o'
dend3 = set(dend3, "nodes_pch", nodepch)
dend3 = set(dend3, "nodes_cex", nodecex) 
nodedf = data.frame(cbind(get_nodes_xy(dend3, type = 'rectangle')[bkdf$Node,]), node=as.character(bkdf$Name))

# Make rectangles if breaks pair contains a chosen cluster:
acut = cutree(dend, 1)
acut = acut[order.dendrogram(dend)]
ll = c()
for (i in 1:nrow(bkdf)){
    nam = unlist(stlab[bkdf$Node[i]])
    acut[nam] = i + 1
    ll = rbind(ll, data.frame(id=nam, Name=bkdf$Name[i], Num=bkdf$Num[i]))
}
write.table(ll, paste0(setprefix, 'tissue_bkpts_subsets.tsv'),quote=F, sep="\t", row.names=F)
breaks = c(0,calc.breaks.acut(acut) * NLAB, NLAB)
boxdf = c()
for (i in 1:(length(breaks)-1)){
    id = which(breaks[i] < nodedf$X1 & breaks[i+1] > nodedf$X1)
    if (length(id) > 0){
        boxdf = rbind(boxdf, data.frame(X1=breaks[i], X2=breaks[i+1]))
    }
}

circleplot = function(with.donors=FALSE, with.boxes=TRUE, with.tree=TRUE, only.diff=FALSE){
    circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
    circos.initialize(factors = "single", xlim = c(0, NLAB)) 
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
                     circos.text(1:NLAB-0.5, rep(0, NLAB), labels(dend3), col = labels_colors(dend3), cex=.25,
                                 facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) }, bg.border = NA, track.height = 0.1)
    circos.track(ylim = c(0, 5), bg.border = NA, panel.fun = function(x, y) {
                     m = t(colmat[,5:1])
                     nr = nrow(m)
                     nc = ncol(m)
                     for(i in 1:nr) {
                         circos.rect(1:nc - 1, rep(nr - i, nc), 
                                     1:nc, rep(nr - i + 1, nc), 
                                     border = NA, col = m[i, ])
                     }
                     }, track.height=0.125) 
    max_height = max(attr(dend3, "height"))
    if (with.tree){
        circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                         circos.dendrogram(dend3, max_height = max_height) },
                         track.height = 0.5, bg.border = NA)
        circos.clear()
    }
    # Inner dendrogram
    if (with.boxes){
        par(new = TRUE) # <- allows new fig inside.
        lim = 1.295
        circos.par("canvas.xlim" = c(-lim, lim), "canvas.ylim" = c(-lim, lim))
        circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
        circos.initialize(factors = 'single', xlim = c(0, NLAB))
        circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                         circos.rect(boxdf$X1, rep(max_height, nrow(boxdf)),
                                     boxdf$X2, rep(.3 * max_height, nrow(boxdf)),
                                     lwd=1, lty='dashed', col=NA,border='darkgrey')
                         circos.text(nodedf$X1, .95 * max_height - nodedf$X2, as.character(nodedf$node),
                                     col = 'black', cex=.75,
                                     facing = "reverse.clockwise", niceFacing=TRUE, adj = c(0, 0.5)) 
                     }, track.height = 0.5, bg.border = NA)
        circos.clear()
    }
    if (with.donors){
        par(new = TRUE) # <- allows new fig inside.
        # lim = 1.295
        lim = .65
        circos.par("canvas.xlim" = c(-lim, lim), "canvas.ylim" = c(-lim, lim))
        circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
        circos.initialize(factors = 'single', xlim = c(0, NLAB))
        circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                         # Order is: lab:
                         allpair = c()
                         for (donor in dlist$Donor){
                             indx = which(meta$Donor == donor)
                             ids = meta$id[indx]
                             groups = meta$GROUP[indx]
                             pts = which(lab %in% ids)
                             # Make links for all pairs
                             if (length(pts) < 3){
                                 pairs = pts
                             } else {
                                 pairs = t(combn(pts, 2))
                             }
                             allpair = rbind(allpair, pairs)
                         }
                         allpair = unique(allpair)
                         n = length(x)
                         for (i in 1:nrow(allpair)){
                             circos.link('single', allpair[i,1], 'single', allpair[i,2], 
                                         h.ratio=.8, col=rgb(0,0,0,.5))
                         }
                     }, track.height = 0.5, bg.border = NA)
    }
    circos.clear()
}



pdf(paste0(imgpref,'full_tree_metadata_circ_bkpts_',NCLUST,'.pdf'),width=14.5,height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot()
upViewport()
draw(pd.legend, x = circle_size, just = "left")
dev.off()


pdf(paste0(imgpref,'full_tree_metadata_circ_donors.pdf'),width=14.5,height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(with.donors=T, with.tree=F, with.boxes=F)
upViewport()
draw(pd.legend, x = circle_size, just = "left")
dev.off()


pdf(paste0(imgpref,'full_tree_metadata_circ_donors_onlydiff.pdf'),width=14.5,height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(with.donors=T, with.tree=F, with.boxes=F, only.diff=TRUE)
upViewport()
draw(pd.legend, x = circle_size, just = "left")
dev.off()




# --------------------------------------------
# Evaluate different methods of building tree:
# --------------------------------------------
# Try different orderings:
hclust_methods <- c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2")
for (method in hclust_methods){
    print(method)
    dend <- hclust(dt, method=method) %>% as.dendrogram
    # Reorder metadata:
    lab = labels(dend)
    metamat = metamat[lab,]
    info = meta[lab, 'infoline']
    col = meta[lab, 'COLOR']
    # Edit file:
    dend2 = set(dend, "labels", info)
    labels_colors(dend2) <- col
    colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(NCLUST)
    dend3 <- color_branches(dend2, k=NCLUST, col=colpair)
    dend3 = set(dend3, "labels_cex", .15)
    # Make plot:
    pdf(paste0(imgpref,'full_tree_metadata_horiz_',method,'.pdf'),width=3,height=15)
    layout(matrix(c(1,2,3), 1,3), widths=c(8,1,1), TRUE)
    par(mar=c(4,1,1,3))
    par(yaxs="i")
    par(xaxs="i")
    plot(dend3, horiz=TRUE)
    par(mar=c(4,.25,1,.05))
    meta.image(metamat[lab,5:1], colvals=colvals, cex=.6, horiz=TRUE)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    par(mar=c(4,.05,1,.25))
    avail = as.matrix(wm[main.marks, lab])
    image(avail, yaxt='n',xaxt='n', col=c('white','darkgrey'))
    text(x=seq(0,1, length.out=length(main.marks)),
         y=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), 
         labels=main.marks, srt=90, adj=1, xpd=true,cex=.25)
    dev.off()
}

full_dendlist <- dendlist()
for(i in seq_along(hclust_methods)) {
   hc_full <- hclust(dt, method = hclust_methods[i])   
   full_dendlist <- dendlist(full_dendlist, as.dendrogram(hc_full))
}
names(full_dendlist) <- hclust_methods

pdf(paste0(imgpref,'full_tree_allmethods.pdf'),width=6,height=15)
par(mfrow = c(4,2))
for(i in 1:8) {
    par(mar=c(2,.25,1,.05))
    dend = full_dendlist[[i]] 
    lab = labels(dend)
    info = meta[lab, 'infoline']
    col = meta[lab, 'COLOR']
    # Edit file:
    dend = set(dend, "labels", info)
    labels_colors(dend) <- col
    dend = set(dend, "labels_cex", .15)
    dend %>% color_branches(k=NCLUST, col=colpair) %>% plot(axes = FALSE, horiz = TRUE)
    title(names(full_dendlist)[i])
}
dev.off()

# -----------------------------------------------------------
# Look at positions of each group:
# -----------------------------------------------------------
grouplist=sort(as.character(unique(meta$GROUP)))
pdf(paste0(imgpref,'full_tree_pergroup_horiz.pdf'),width=7,height=11)
layout(matrix(c(1:35), 5,7, byrow=T), TRUE)
for (group in grouplist){
    par(mar=c(0,2.5,3,1.5))
    print(group)
    idx = which(metamat[lab,'group'] == group)
    gleaf_pch = rep(26,length(col))
    gleaf_pch[idx] = 19
    dend3 %>% set("leaves_pch", gleaf_pch) %>% 
        set("leaves_cex", .75) %>% set("labels", '') %>% plot(horiz=T, axes=F)
    mtext(paste0(group, '\n(', length(idx),' samples)'), side=3, line=.5, 
          cex=.8, col=colvals[['group']][group])
}
dev.off()



# dendlist_cor2 <- cor.dendlist(full_dendlist, method = "common")
# corrplot::corrplot(dendlist_cor2, "pie", "lower")
# -----------------------------------------------------------
# Compare the orderings according to each mark, against full:
# -----------------------------------------------------------

# png(paste0(imgpref,'mds_marks_col.png'),res=300,units='in',width=11,height=6.5)

# Full:
dt = as.dist(full)
dend <- hclust(dt, method = "ward.D") %>% as.dendrogram
lab = labels(dend)
info = meta[lab, 'infoline']
col = meta[lab, 'COLOR']
group = meta[lab, 'GROUP']
full.dend = set(dend, "labels", info)
labels_colors(full.dend) <- col

mark = 'H3K27ac'
# For mark
# mat = obsll[[mark]]
mat = ll[[mark]]
mdt = as.dist(mat)
mark.dend <- hclust(mdt, method = "ward.D") %>% as.dendrogram
mlab = labels(mark.dend)
info = meta[mlab, 'infoline']
col = meta[mlab, 'COLOR']
group = meta[mlab, 'GROUP']
mark.dend = set(mark.dend, "labels", info)
labels_colors(mark.dend) <- col


par(mar=c(0.5, 0.5, 2, 0.5))
mat = obsll[[mark]]
subrn = rn[rn %in% colnames(mat)]
mat2 = tmp
mat2[subrn,subrn] = mat[subrn, subrn]
plot.cov(mat2, clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')
mtext(mark, side=3, cex=1.3)
# Imputed:
par(mar=rep(0.5,4))
mat = ll[[mark]]
mat = clamp.mat(mat)
plot.cov(mat[rn, rn], clamp=TRUE, palette=palette, breaks=breaks, breakscol='darkgrey')

dl = dendlist(full.dend, mark.dend)
dl = ladderize(dl)
# dl = untangle(dl, method = "step1side", k_seq = 2:6) 
# dl = untangle(dl,  order=lab)
# , k_seq = 2:6) 
# dl = set(dl, "branches_k_color", k=6)
# labels_colors(full.dend) <- col

# pdf(paste0(imgpref,'full_tree_allmethods.pdf'),width=6,height=15)
# dl = tanglegram(dl, faster = TRUE) # (common_subtrees_color_branches = TRUE)

# ("labels", "ladderize", "random", "step1side", "step2side", "DendSer")
     


# # Plot comparisons:
# full_dendlist %>% dendlist(which = c(3,4)) %>% ladderize %>% 
#     untangle(method = "step1side", k_seq = 2:6) %>%
#     set("branches_k_color", k=2) %>% 
#     tanglegram(faster = TRUE) # (common_subtrees_color_branches = TRUE)

