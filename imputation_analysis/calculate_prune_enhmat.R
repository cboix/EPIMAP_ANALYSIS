#!/usr/bin/R
# --------------------------
# De-dup the enhancer matrix
# --------------------------
gtargs=(commandArgs(TRUE))
print(gtargs)
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
    plotting.trees=TRUE
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
    plotting.trees=TRUE
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))

library(dplyr)
library(cba)
library(GenomicRanges)

# Tree libraries:
library(Matrix)
library(ape)
library(seqinr)
library(phangorn)
library(dendextend)
library(phylogram)
options(scipen=45) # So we dont get issues writing integers into bedfiles
options(repos='http://cran.rstudio.com/')

# NOTE: Using full dataset needs more than 120G of space, may not even work in R

# Arguments:
usetree = 'enhancers'
tol = 2500  # Plus/minus distance - window for enhancer overlaps
plot.trees = FALSE 
setprefix = paste0(usetree, '_dedupmat_')

# Make directories:
gtdir = "gwas_tree_analysis/"
treedir = paste0(gtdir, 'mp_trees/')
imgdir = paste0(img, "gwas_tree_analysis/")
cmd = paste('mkdir -p ', imgdir, gtdir, treedir)
system(cmd)
imgpref = paste0(imgdir, setprefix)
treeimgpref = paste0(imgpref, 'e', tol, '_')

# Plotting functions:
if (plotting.trees){
    library(circlize)
    library(ComplexHeatmap)
    library(gridBase)
    source(paste0(bindir, 'auxiliary_gwastree_functions.R'))

    NCLUST = 20
    meta$uqinfo = paste0(1:nrow(meta), '_', meta$info)
    color.dend = function(dend, cex=.25, k=NCLUST){
        lab = labels(dend)
        info = meta[lab, 'uqinfo']
        col = meta[lab, 'COLOR']
        group = meta[lab, 'GROUP']
        dend2 = set(dend, "labels", info)
        labels_colors(dend2) <- col
        # colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(k)
        # dend3 <- color_branches(dend2, k=k, col=colpair)
        dend3 = set(dend2, "labels_cex", cex)
        return(dend3)
    }
}

# ----------------------------------------
# Load matrix and get epigenomes per node:
# ----------------------------------------
fullmatfile = 'Enhancer_H3K27ac_matrix_062619.mtx.gz'
enhmatfile = 'Enhancer_H3K27ac_matrix_enhonly_062619.mtx.gz'
enhmargfile = 'Enhancer_H3K27ac_margins_enhonly_062619.tsv.gz'
enhrdafile = 'Enhancer_H3K27ac_matrix_enhonly_062619.Rda'
dna.file = 'Enhancer_H3K27ac_enhonly_062619.dna'
matnames = scan('Enhancer_matrix_names.txt', "c")
if (!file.exists(enhrdafile)){
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
    colnames(mm) = matnames
    save(mm, file=enhrdafile)
    rm(mat, matmarg, mrow)
    gc()
} else {
    load(enhrdafile)
}

# Load in the margins:
enhmargfile = 'Enhancer_H3K27ac_margins_enhonly_062619.tsv.gz'
matmarg = read.delim(gzfile(enhmargfile), sep="\t", header=T)
# 316,316 locations are exactly 1 --> remove these in the overall tree calc. 
# We will re-add them in the branch length calculations
flag1 = which(matmarg$col == 1)
length(flag1)

flag3 = which(matmarg$col == 2) # Potentially, also flag 1-2
length(flag3)


# ---------------------------------------
# Quantify the multiplicity of enhancers:
# ---------------------------------------
dedup.file = paste0(gtdir, 'dedup_matrix_indices.Rda')
if (!file.exists(dedup.file)){
    ddir = 'DHS_Index_WM201902/'
    dpref = 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19'
    dmlfile = paste0(ddir, dpref, '.core.srt.txt')
    dmlnamfile = paste0(ddir, dpref, '_r200_e0_names.core.srt.tsv')
    dmlrdafile = paste0(ddir, dpref, '_enh.core.srt.Rda')
    load(dmlrdafile)

    # -------------------------------------------------
    dmgr = GRanges(enhdf$chr, IRanges(enhdf$start, enhdf$end))
    crossovl = findOverlaps(dmgr, dmgr, minoverlap=50)
    crossovl = data.frame(crossovl)
    crossovl = crossovl[crossovl[,1] < crossovl[,2],]
    keep.mult = which(!(crossovl[,1] %in% flag1 | crossovl[,2] %in% flag1))
    crossovl = crossovl[keep.mult,]
    dim(crossovl)
    # With no min-overlap: There are 545k overlaps, 491k in non-1 regions. 
    # With min-overlap = 50: There are 304k in non-1 regions. 

    # For each pair, calculate jaccard sim:
    m2 = mm[crossovl[,1],] * mm[crossovl[,2],]
    cross.int = apply(m2, 1, sum)
    cross.both = matmarg$col[crossovl[,1]] + matmarg$col[crossovl[,2]]
    gc()
    cross.jacc = cross.int / (cross.both - cross.int)

    png(paste0(imgpref, 'jaccard_of_overlapping_enhancers.png'), width=8, height=6, units='in', res=450)
    hist(cross.jacc, col='darkgrey', border='white', 
         ylab='Number of Pairs', xlab='Jaccard similarity of enhancer activity', main='')
    dev.off()

    # 66k are exactly 1: 
    sum(cross.jacc == 1)
    sum(cross.jacc > .75) # About 138k pairs that are gt than 0.75+
    # Properties of enh:
    cross.out = crossovl[cross.jacc > 0.75,]
    cross.enh = data.frame(enh=c(cross.out[,1], cross.out[,2]), id=1)
    cross.num = aggregate(id ~ enh, cross.enh, sum)
    sum(cross.num$id > 1) # 23.6k used more than once

    # Now, choose which stays:
    cross.out$m1 = matmarg$col[cross.out[,1]] 
    cross.out$m2 = matmarg$col[cross.out[,2]]
    cross.out$gt12 = 1 * (cross.out$m1 > cross.out$m2)
    cross.out$kept = with(cross.out, gt12 * queryHits + (1 - gt12) * subjectHits)
    all.ind = sort(unique(c(cross.out$queryHits, cross.out$subjectHits)))
    kept.ind = sort(unique(c(cross.out$kept)))
    flag2 = all.ind[!(all.ind %in% kept.ind)]
    cross.out$k1 = cross.out[,1] %in% kept.ind
    cross.out$k2 = cross.out[,2] %in% kept.ind
    cross.out[cross.out$k1 + cross.out$k2 == 2,]
    # Kept indices
    # flagged.ind = unique(c(flag1, flag2))
    flagged.ind = unique(c(flag1, flag2, flag3))
    kept.mat.ind = 1:nrow(matmarg)
    kept.mat.ind = kept.mat.ind[!(kept.mat.ind %in% flagged.ind)]
    # Make and save matrix:
    mm2 = mm[kept.mat.ind,]
    print(dim(mm2))
    save(mm2, kept.mat.ind, file=dedup.file)
} else {
    load(dedup.file)
}




# ---------------------------------
# Calculate the Jaccard similarity 
# ---------------------------------
dedup.jacc.file = paste0(gtdir, 'dedup_jacc_matrices.Rda')
if (!(file.exists(dedup.jacc.file))){
    jacc.int = matrix(NA, nrow=ncol(mm2), ncol=ncol(mm2))
    rs = rep(0, ncol(mm2))
    chunk = 40
    for (i in 1:21){
        print(i)
        j = ((i-1) * chunk + 1):min(i * chunk, ncol(mm2))
        vec = t(mm2)[j,] %*% mm2
        jacc.int[j, ] = as.matrix(vec)
        rs[j] = apply(mm2[,j], 2, sum)
    }
    gc()

    NT = length(rs)
    rsm = matrix(rep(rs, NT), nrow=NT,ncol=NT, byrow=T)
    rsm = t(rsm) + rsm

    jacc.union = rsm - jacc.int
    jacc.dt = 1.0 - jacc.int / jacc.union 
    rownames(jacc.dt) = colnames(mm2)
    colnames(jacc.dt) = colnames(mm2)
    save(jacc.int, jacc.union, jacc.dt, file=dedup.jacc.file)
} else {
    load(dedup.jacc.file)
}


# --------------------------------------------------
# Plot the complete-linkage tree of de-dupped matrix
# --------------------------------------------------
dt <- as.dist(jacc.dt)
method = 'complete'
ht <- hclust(dt, method=method)

# Plot the end results:
# NL = length(labels(ht))
# tree = root(ht, outgroup=labels(bin_optim)[NL], resolve.root = TRUE)
dend = as.dendrogram(ht)
lab = labels(dend)
dend2 = color.dend(dend)
NCLUST=20
dend3 = dend2

# Plot the tree:
pdf(paste0(treeimgpref, 'jacc_complete_dedupped.pdf'), width=14.5, height=12, onefile=T)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dend3, lab=lab)
upViewport()
circos.clear()
title(paste0('Complete-linkage clustering on Jaccard of de-dupped matrix'))
dev.off()



# ------------------------------------------
# Calculate + plot the neighbor joining tree 
# ------------------------------------------
tree = NJ(jacc.dt)
NL = length(labels(tree))
tree = root(tree, outgroup=labels(tree)[NL], resolve.root = TRUE)
dend = as.dendrogram(tree)
lab = labels(dend)
dend2 = color.dend(dend)
NCLUST=20
dend3 = dend2

# Plot the tree:
pdf(paste0(treeimgpref, 'nj_tree_dedupped.pdf'), width=14.5, height=12, onefile=T)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend=dend3, lab=lab)
upViewport()
circos.clear()
title(paste0('Neighbor Joining on Jaccard of de-dupped matrix'))
dev.off()





