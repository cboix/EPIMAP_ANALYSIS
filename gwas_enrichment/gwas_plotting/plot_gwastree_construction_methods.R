#!/usr/bin/R
# ----------------------------------------------
# Evaluates the construction of different trees:
# ----------------------------------------------
gtargs=(commandArgs(TRUE))
print(gtargs)
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))
source(paste0(bindir, 'auxiliary_gwastree_functions.R'))

library(ggplot2)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(GenomicRanges)
library(dplyr)
library(cba)

library(Matrix)
# Tree libraries:
library(ape)
library(seqinr)
library(phangorn)
library(dendextend)
library(phylogram)
options(scipen=45) # So we dont get issues writing integers into bedfiles
options(repos='http://cran.rstudio.com/')

# Arguments:
usetree = 'enhancers'
tol = 2500  # Plus/minus distance - window for enhancer overlaps
plot.trees = FALSE 
setprefix = paste0(usetree, '_treemethods_')

# Make directories:
imgdir = paste0(img, "gwas_tree_analysis/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, setprefix)
treeimgpref = paste0(imgpref, 'e', tol, '_')

# --------------------------------------------
# Load enhancer jaccard - see different trees:
# --------------------------------------------
print('[STATUS] Loading enhancers jaccard matrix:')
emat = read.delim(gzfile('Enhancer_jaccard.tsv.gz'), sep="\t", header=F)
matnames = scan('Enhancer_matrix_names.txt', "c")
rownames(emat) = matnames
colnames(emat) = matnames
if (usetree == 'roadmap'){
    emat = emat[keepbss, keepbss]
}
jacc.dt <- as.dist(emat)

meta$uqinfo = paste0(1:nrow(meta), '_', meta$info)

color.dend = function(dend, cex=.25){
    lab = labels(dend)
    info = meta[lab, 'uqinfo']
    col = meta[lab, 'COLOR']
    group = meta[lab, 'GROUP']
    dend2 = set(dend, "labels", info)
    labels_colors(dend2) <- col
    NCLUST=20
    colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(NCLUST)
    # dend3 <- color_branches(dend2, k=NCLUST, col=colpair)
    dend3 = set(dend2, "labels_cex", cex)
    return(dend3)
}


# Trees from jaccard distance:
methodlist = c('ward.D', 'ward.D2', 'complete', 'single', 'average', 'upgma','nj')
for(method in methodlist){
    print(method)
    if (method == 'upgma'){
        ht = upgma(jacc.dt)
    } else if(method == 'nj'){
        ht = NJ(jacc.dt)
    } else {
        ht <- hclust(jacc.dt, method=method)
    }
    dend = as.dendrogram(ht)
    dend2 = color.dend(dend)
    NCLUST=20
    dend3 <- color_branches(dend2, k=NCLUST, col=colpair)
    NL = length(lab)
    # Plot the tree:
    pdf(paste0(treeimgpref, sub("\\.","_",method), "_link_jacc.pdf"), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    title(paste0(method, '-linkage of Jaccard Similarity of Enhancers'))
    dev.off()
}


# Plot/save NJ tree separately:
ht = NJ(jacc.dt)
write.tree(ht, '~/tree_nj.nwk')
lab = labels(ht)
info = meta[lab, 'infoline']
col = meta[lab, 'COLOR']
labels(ht) = info
dend = as.dendrogram(ht)

pdf('~/test.pdf', width=4, height=20)
par(mar=rep(0,4))
plot.phylo(ht, tip.color=col, cex=.20, edge.width=.25)
dev.off()

ht.nj = NJ(jacc.dt)
ht.jacc <- hclust(jacc.dt, method='complete')

## convert phylo objects to dendrograms
dnd1 <- as.dendrogram(ht.nj)
dnd2 <- as.dendrogram(ht.jacc)
dnd1 = color.dend(dnd1)
dnd2 = color.dend(dnd2)
## rearrange in ladderized fashion
dnd1 <- ladder(dnd1)
dnd2 <- ladder(dnd2)
## plot the tanglegram
dndlist <- dendextend::dendlist(dnd1, dnd2)
pdf('~/test_comp.pdf', width=8, height=20)
dendextend::tanglegram(dndlist, fast = TRUE, margin_inner=2, edge.lwd=.25, lwd=.15)
dev.off()

# ----------------------------------------
# Load matrix and get epigenomes per node:
# ----------------------------------------
fullmatfile = 'Enhancer_H3K27ac_matrix_062619.mtx.gz'
enhmatfile = 'Enhancer_H3K27ac_matrix_enhonly_062619.mtx.gz'
enhmargfile = 'Enhancer_H3K27ac_margins_enhonly_062619.tsv.gz'
dna.file = 'Enhancer_H3K27ac_enhonly_062619.dna'
print("[STATUS] Making interleaved format file")
if (!file.exists(enhmatfile)){
    print("[STATUS] Loading enhancer matrix")
    # Load in the product of enhancer and H3K27ac mtx: 
    mat = read.delim(gzfile(fullmatfile), sep="\t", header=F)
    # Matrix is 0-indexed - turn to 1-indexing
    mat[,1] = mat[,1] + 1
    mat[,2] = mat[,2] + 1
    names(mat) = c('row','col')
    # Keep only enhancers:
    kid = which(mat$row %in% enhind)
    mat = mat[kid,]
    rm(kid)
    # Margin (for weighted regression):
    matmarg = aggregate(col ~ row, mat, length)
    matmarg = matmarg[order(matmarg$row),]
    print("Saving just enhancers. Might take a while to write.")
    write.table(matmarg, gzfile(enhmargfile), quote=F, sep="\t", row.names=F)
    write.table(mat, gzfile(enhmatfile), quote=F, sep="\t", row.names=F)
} else { 
    mat = read.delim(gzfile(enhmatfile), sep="\t", header=T)
    matmarg = read.delim(gzfile(enhmargfile), sep="\t", header=T)
}
matnames = scan('Enhancer_matrix_names.txt', "c")
# Mapping all to just enhancers:
enhind = as.numeric(scan('Enhancer_indices.txt', 'c')) + 1
enhmap = rep(0, max(enhind))
enhmap[enhind] = 1:length(enhind)
mrow = enhmap[mat$row]
mat$row = mrow
matmarg$row = enhmap[matmarg$row]
# Create Matrix, fill in:
mm = Matrix(0, nrow=nrow(matmarg), ncol=max(mat$col))
mm[as.matrix(mat)] = 1
gc()


# ---------------------------------------
# Compute distance metrics on the matrix:
# ---------------------------------------
# Jaccard
NF = ncol(mm)
jacc.mat = matrix(NA, nrow=NF, ncol=NF, dimnames=list(matnames, matnames))
colmarg = colSums(mm)
for (i in 1:NF){
    print(i)
    xi = mm[,i]
    sapply(i:NF, function(x){ 
               xj = mm[,j]
               int = sum(xi * xj) 
               union = (colmarg[i] + colmarg[j] - int)
               return(int / union) })

    for (j in i:NF){
        xj = mm[,j]
        int = sum(xi * xj) 
        union = (colmarg[i] + colmarg[j] - int)
        jacc.mat[i,j] = int / union
        jacc.mat[j,i] = int / union
    }
}


# Copy for tf-idf:
mm2 = mm

# tf-IDF, to dot product/ejaccard?
mm2 = sweep(mm, 2, colmarg, '/')
matmarg$idf = -log(matmarg$col / ncol(mm))
# Remove idpt:


# writeMM(mm, matrix.file)
# jmm = jacc.dist(mm)

# Compute jaccard here
# jaccard(mm)


if (!file.exists(dna.file)){
    write.table(paste0(ncol(mm), '\t', nrow(mm)), file=dna.file, quote=F, sep="\t", row.names=F, col.names=F)
    # Write line by line:
    for (i in 1:ncol(mm)){
        print(i)
        mat.col = as.character(mm[,i])
        mat.col[mat.col == '1'] = 'G'
        mat.col[mat.col == '0'] = 'A'
        mat.col = paste0(mat.col, collapse='')
        write.table(paste0(matnames[i], '\t', mat.col), file=dna.file, quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
    }
    gc()
}

# Make one pseudo alignment per cluster:
phydir = "enh_as_phylip/"
moddir = paste0(phydir, 'modules/')
cmd = paste('mkdir -p', phydir, moddir)
system(cmd)

# Load in modules, write to each file:
load('modules_enhancer_sets.Rda')
NMOD = length(enhsets)
for (k in 1:NMOD){
    print(k)
    ind = enhmap[enhsets[[k]]]
    kpref = paste0(moddir, 'Enhancer_H3K27ac_enhonly_062619_module_', k)
    mod.file = paste0(kpref, '.phylip')
    mod.tree.file = paste0(kpref, '_nj.tree')
    if (!file.exists(mod.file)){
        write.table(paste0(ncol(mm), '\t', length(ind)), file=mod.file, quote=F, sep="\t", row.names=F, col.names=F)
        pb = txtProgressBar(min=0, max=ncol(mm), style = 3)
        for (i in 1:ncol(mm)){
            setTxtProgressBar(pb, i)
            mat.col = as.character(mm[ind,i])
            mat.col = paste0(mat.col, collapse='')
            write.table(paste0(matnames[i], '\t', mat.col), file=mod.file, quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
        }
        close(pb)
    }
    # Also make tree:
    if (!file.exists(mod.tree.file)){
        submat = as.matrix(mm[ind,])
        sub.jacc = jacc.dist(t(submat))
        sub.jacc[is.na(sub.jacc)] = median(sub.jacc, na.rm=T)
        sub.tree = NJ(sub.jacc)
        write.tree(sub.tree, file=mod.tree.file)
    }
    gc()
}

# Load in modules, write to each file:
load('modules_enhancer_sets.Rda')
NMOD = length(enhsets)
for (k in 1:NMOD){
    print(k)
    ind = enhmap[enhsets[[k]]]
    mod.file = paste0('Enhancer_H3K27ac_enhonly_062619_module_', k, '.phylip')
    if (!file.exists(mod.file)){
        write.table(paste0(ncol(mm), '\t', length(ind)), file=mod.file, quote=F, sep="\t", row.names=F, col.names=F)
        pb = txtProgressBar(min=0, max=ncol(mm), style = 3)
        for (i in 1:ncol(mm)){
            setTxtProgressBar(pb, i)
            mat.col = as.character(mm[ind,i])
            mat.col = paste0(mat.col, collapse='')
            write.table(paste0(matnames[i], '\t', mat.col), file=mod.file, quote=F, sep="\t", row.names=F, col.names=F, append=TRUE)
        }
        close(pb)
        gc()
    }
}


# rm(mat, matmarg, mm)
# gc()


# NOTE: CANNOT READ IN FULL DNA FILE
# enh.dna <- read.dna(dna.file, format="interleaved")
# mammals_phyDat <- phyDat(mammals, type = "DNA", levels = NULL)
# mammals10 <- subset(mammals_phyDat, 1:10)
# mammals10_phyDat <- phyDat(mammals10, type = "DNA", levels = NULL)

tree.upgma <- upgma(jacc.dt)
tree.nj  <- NJ(jacc.dt)





# mp10k.tree = read.tree('10k_first.phylip.parstree')
mp10k.tree = read.tree('50k_first.phylip.parstree')
dend = as.dendrogram(mp10k.tree)
lab = labels(dend)
dend2 = color.dend(dend)
NCLUST=20
# dend3 <- color_branches(dend2, k=NCLUST, col=colpair)
dend3 = dend2
NL = length(labels(dend))
# Plot the tree:
pdf(paste0(treeimgpref, sub("\\.","_",'first_10k'), "_mpars.pdf"), width=14.5, height=12, onefile=T)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dend3, lab=lab)
upViewport()
circos.clear()
title('First 10k enhancers - max parsimony on binary activities')
dev.off()



exseq <- read.dna("10k_first.dna", format="interleaved")
exseq_phyDat <- phyDat(exseq, type="DNA", levels=NULL)

ex10 <- subset(exseq_phyDat, 1:10)
ex10_phyDat <- phyDat(ex10, type="DNA", levels=NULL)

# Perform model test? --> Should be binary
mt <- modelTest(ex10)
print(mt)
dm <- dist.ml(exseq, model="JC69")



mat10k = as.matrix(mm[1:1e4,])
colnames(mat10k) = matnames
binseq = as.phyDat(t(mat10k), type="USER", levels = c(0, 1))

# Preliminary distance (to improve):
dbin <- dist.ml(binseq, model="JC69")
dbin[is.na(dbin)] = median(dbin, na.rm=T)
dbin[is.infinite(dbin)] = 2

bin_UPGMA <- upgma(dbin)
bin_NJ  <- NJ(dbin)
plot(bin_UPGMA, main="UPGMA")
plot(bin_NJ, main = "Neighbor Joining")

# Calculate parsimony:
parsimony(bin_UPGMA, binseq)
parsimony(bin_NJ, binseq)
# Starting from NJ (better), optimize parsimony (on 10k ~5min?)
bin_optim <- optim.parsimony(bin_NJ, binseq)

dend = as.dendrogram(bin_optim)
lab = labels(dend)
dend2 = color.dend(dend)
NCLUST=20
# dend3 <- color_branches(dend2, k=NCLUST, col=colpair)
dend3 = dend2
NL = length(labels(dend))
# Plot the tree:
pdf(paste0(treeimgpref, 'nj_optim_10k', "_mpars.pdf"), width=14.5, height=12, onefile=T)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dend3, lab=lab)
upViewport()
circos.clear()
title('First 10k enhancers - max parsimony on binary activities')
dev.off()



# Optimize parsimony using parsimony ratchet (on 10k ~ ....)
bin_pratchet <- pratchet(binseq)
write.tree(bin_pratchet, 'pratchet_10k_mpars.tree')

dend = as.dendrogram(bin_pratchet)
lab = labels(dend)
dend2 = color.dend(dend)
NCLUST=20
dend3 = dend2
NL = length(labels(dend))
# Plot the tree:
pdf(paste0(treeimgpref, 'pratchet_10k', "_mpars.pdf"), width=14.5, height=12, onefile=T)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dend3, lab=lab)
upViewport()
circos.clear()
title('First 10k enhancers - max parsimony on binary activities')
dev.off()










# ------------------------------------
# Get tree from fused distance matrix:
# ------------------------------------
if (usetree == 'correlation'){
    print("MARK CORRELATION")
    print(usetree)
    dt <- as.dist(full)
    method='ward.D'
} else if (usetree == 'enhancers' || usetree == 'roadmap'){
    print("ENH JACCARD")
    dt <- as.dist(emat)
    method = 'complete'
} else if (usetree == 'gwastree'){
    print("GWAS JACCARD")
    dt <- as.dist(emat)
    dt = dist(epimat > 0, method='jaccard')
    method='ward.D'
}
ht <- hclust(dt, method=method)
ht$order <- order.optimal(dt, ht$merge)$order

# -------------------------
# Make the dendextend tree:
# -------------------------
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
NL = length(lab)
memb = get_nodes_attr(dend3, 'member')
NN = length(memb)

# For relabeling:
names(lab) = NULL
# NOTE: Will only work with enh. may need to fix 
labmapping = sapply(lab, function(x){which(matnames == x)})
