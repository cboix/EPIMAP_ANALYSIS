#!/usr/bin/R
# ----------------------------------------
# Plot breakpoints on enhancer tree matrix
# ----------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'auxiliary_gwastree_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
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

imgdir = paste0(img, "gwas_tree_analysis/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'enhancers_')


# ------------------------------------------
# Identify, plot, and write out breakpoints:
# ------------------------------------------
print('[STATUS] Loading enhancers jaccard matrix:')
emat = read.delim(gzfile('Enhancer_jaccard.tsv.gz'), sep="\t", header=F)
matnames = scan('Enhancer_matrix_names.txt', "c")
rownames(emat) = matnames
colnames(emat) = matnames
# Tree from jaccard distance:
dt <- as.dist(emat)
method = 'ward.D'
method = 'ward.D2'
method = 'complete'
ht <- hclust(dt, method=method)
cocl <- order.optimal(dt, ht$merge)$order
# tree = as.phylo(ht)
# dend = as.dendrogram(tree)
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
NLAB = length(lab)

# Size/internal:
nl = get_nodes_attr(dend, 'members')
ind = which(nl >= 10)
stlab = partition_leaves(dend)
mlist = sapply(ind,function(i){ meta[stlab[[i]],'GROUP'] })
ilist = sapply(ind,function(i){ meta[stlab[[i]],'infoline'] })
intern = sapply(mlist, function(x){ tab = table(as.character(x)); max(tab)/sum(tab) })
len = lapply(mlist, length)

# Choose nodes:
possible = ind[which(intern > 0.7)]
large = ind[which(unlist(len) >= 7)]
large = large[!(large %in% possible)]
pos2 = sort(unique(c(possible, large)))

# Lab nodes
nodecex = rep(1, nnodes(dend3))
nodepch = rep(NA, nnodes(dend3))
nodecex[pos2] = 1
nodepch[possible] = 'o'
nodepch[large] = 'x'
dend3 = set(dend3, "nodes_pch", nodepch)
dend3 = set(dend3, "nodes_cex", nodecex) 
nodedf = data.frame(cbind(get_nodes_xy(dend3, type = 'rectangle')[pos2,], node=pos2))
dend <- color_branches(dend, col=rgb(0,0,0,.5))

pdf(paste0(imgpref,'gwastree_possible_',NCLUST,'.pdf'),width=12,height=12)
max_height = attr(dend3, 'height')
circos.initialize(factors = 'single', xlim = c(0, NLAB))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                 circos.dendrogram(dend, max_height = max_height) 
                 circos.text(nodedf$V1, max_height - nodedf$V2, as.character(nodedf$node),
                             col = 'black', cex=.5,
                             facing = "reverse.clockwise", niceFacing=TRUE, adj = c(0, 0.5))
}, bg.border = NA, track.height = .8)
dev.off()


# --------------------------------------------
# Plot breakpoints chosen with the full tree: 
# --------------------------------------------
memb = get_nodes_attr(dend3, 'member')
NN = length(memb)
dend4 = set(dend3, 'nodes_pch', 1:NN)
dend4 = set_labels(dend4, 1:NLAB)
declist = list(dec=sapply(rep(NA, NN), list), isleaf=rep(NA, NN), parent=rep(NA, NN))
declist = get_dec(dend4, declist=declist)
nodeloc = get_nodes_xy(dend3, type = 'rectangle')
max_height = attr(dend3, 'height')


# Get breakpoints:
bkdf = read.delim('Annotation/gwastree_breakpts.tsv',header=T, stringsAsFactors=F)
bkdf = bkdf[bkdf$Name != '',]
# Get positions/cex:
pos3 = bkdf$Node
nodecex[] = 0
nodecex[pos3] = .5
dend3 = set(dend3, "nodes_pch", 19)
dend3 = set(dend3, "nodes_cex", nodecex) 
boxlabdf = data.frame(cbind(get_nodes_xy(dend3, type = 'rectangle')[bkdf$Node,]), node=as.character(bkdf$Name))
# Get breaks:
boxdf = c()
for (i in 1:nrow(bkdf)){
    rnx = declist$dec[[bkdf$Node[i]]]
    boxdf = rbind(boxdf, data.frame(X1=min(rnx)-1, X2=max(rnx), Y1=max_height, 
                                    Y2=max(0, max_height - nodeloc[bkdf$Node[i],2])))
}
boxdf$col= ifelse((1:nrow(boxdf))%%2==1, 'black','grey50')
boxdf$lwd=0.75

# Plot figure:
pdf(paste0(imgpref,'gwastree_annotated_',NCLUST,'.pdf'),width=14.5,height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dend3, lab=lab, boxdf=boxdf, boxlabdf=boxlabdf)
upViewport()
circos.clear()
draw(pd.legend, x=circle_size, just="left")
dev.off()

# Plot figure:
pdf(paste0(imgpref,'gwastree_annotated_',NCLUST,'_horiz.pdf'),width=14.5,height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dend3, lab=lab, boxdf=boxdf, boxlabdf=boxlabdf, boxhoriz=TRUE)
upViewport()
circos.clear()
draw(pd.legend, x=circle_size, just="left")
dev.off()



