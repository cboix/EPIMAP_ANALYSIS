#!/usr/bin/R
# ---------------------------------
# Calculate the consensus gwas tree
# on the set of 2.1M enhancers
# ---------------------------------
gtargs=(commandArgs(TRUE))
print(gtargs)
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
    plotting.trees=FALSE
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
    plotting.trees=TRUE
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))

library(dplyr)
library(cba)

# Tree libraries:
library(Matrix)
library(ape)
library(seqinr)
library(phangorn)
library(dendextend)
library(phylogram)
library(phytools) # For tree consensus methods
options(scipen=45) # So we dont get issues writing integers into bedfiles
options(repos='http://cran.rstudio.com/')

# Arguments for chunking enhancer list:
chunksize = 1e5

# Arguments:
usetree = 'enhancers'
tol = 2500  # Plus/minus distance - window for enhancer overlaps
plot.trees = FALSE 
setprefix = paste0(usetree, '_treemethods_')

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



# ---------------------------------------
# Load in the tree files, without matrix:
# ---------------------------------------
chunkpref = paste0(chunksize, "_")
treeimgpref = paste0(treeimgpref, chunkpref)
treepref = paste0(treedir, setprefix, chunkpref)
pr.files = list.files(treedir, pattern=paste0(setprefix, chunkpref, '.*pratchet_mpars.tree'))

emat = read.delim(gzfile('Enhancer_jaccard.tsv.gz'), sep="\t", header=F)
matnames = scan('Enhancer_matrix_names.txt', "c")
rownames(emat) = matnames
colnames(emat) = matnames

pr.trees = list()
internal.trees = list()
for (file in pr.files){
    pr.trees[[file]] = read.tree(paste0(treedir, file))
    # NOTE: Basic bl method (long-term, shouldnt do on full mat):
    pr = pr.trees[[file]]
    internal.trees[[file]] = nnls.phylo(pr, emat)
    # TODO: Load actual scored matrices (by acctran)
}


# Ape consensus (strict + majority rule --> could do scale between these)
cons.strict.tree = ape::consensus(pr.trees, p=1, check.labels=TRUE)
cons.tree = ape::consensus(pr.trees, p=0.5, check.labels=TRUE)


if (plotting.trees){
    # Plot the end results:
    dend = as.dendrogram(cons.tree)
    lab = labels(dend)
    dend2 = color.dend(dend)
    NCLUST=20
    dend3 = dend2
    NL = length(labels(dend))

    # Plot the tree:
    pdf(paste0(treeimgpref, 'mpars_cons_p5.pdf'), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    title(paste0('Majority-rule consensus tree from max parsimony from NJ on binary activities'))
    dev.off()


    # Plot the end results:
    dend = as.dendrogram(cons.strict.tree)
    lab = labels(dend)
    dend2 = color.dend(dend)
    NCLUST=20
    dend3 = dend2
    NL = length(labels(dend))

    # Plot the tree:
    pdf(paste0(treeimgpref, 'mpars_cons_1_strict.pdf'), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    title(paste0('Strict consensus tree from max parsimony from NJ on binary activities'))
    dev.off()

}


# From phytools:
sub.trees = internal.trees[c(1:5)]

rf.tree <- averageTree(sub.trees, method="symmetric.difference")
# rf.tree <- averageTree(internal.trees, method="symmetric.difference")
# symmetric.difference throws error??

# TODO: Make this work with more trees?
qd.tree = averageTree(sub.trees, start=sub.trees[[1]], method="quadratic.path.difference", tol=1e-12, quiet=FALSE)
# qd.tree = averageTree(internal.trees, start=NULL, method="quadratic.path.difference", tol=1e-12, quiet=FALSE)

ls.consensus(sub.trees, start=sub.trees[[1]], tol=1e-3, quiet=FALSE)

if (plotting.trees){

    # TODO: ROOT?
    # Plot the end results:
    dend = as.dendrogram(qd.tree)
    lab = labels(dend)
    dend2 = color.dend(dend)
    NCLUST=20
    dend3 = dend2
    NL = length(labels(dend))

    # Plot the tree:
    pdf(paste0(treeimgpref, 'mpars_qdcons.pdf'), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    title(paste0('Quadratic path difference consensus tree from max parsimony trees'))
    dev.off()

}


# ls.consensus(internal.trees, start=NULL, tol=1e-12, quiet=FALSE, ...)
# minTreeDist(internal.trees[[1]], internal.trees, method="quadratic.path.difference", ...)



library(doMC)
options(cores = 4)
registerDoMC()

trees = internal.trees
D <- Reduce("+", lapply(trees, function(x,t) cophenetic(x)[t,t], t=trees[[1]]$tip.label))/length(trees)
# D = as.matrix(emat)

if(is.null(start)) start<-NJ(D)
if(hasArg(ultrametric)) { 
    ultrametric<-list(...)$ultrametric  ## should the consensus tree be ultrametric
} else { ultrametric<-all(sapply(trees,is.ultrametric))}

if(ultrametric&&!is.rooted(start)) { start<-midpoint(start)}

curr.file = paste0(treepref, 'mpars_chunks_ls_consensus.tree')
if (!file.exists(curr.file)){
    curr <- nnls.tree(D, tree=start, rooted=ultrametric,trace=0)
} else {
    curr = read.tree(curr.file)
}

tol = 1e-3
if(hasArg(optNNI)) { optNNI<-list(...)$optNNI } else { optNNI<-TRUE }
Q = Inf
Qp = attr(curr,"RSS")
if(is.null(Qp)) Qp <- phytools:::rss(D,curr)
ct = 0
curr.file = paste0(treepref, 'mpars_chunks_ls_consensus.tree')
while((Q-Qp)>tol){
    print(ct)
    print(Qp)
    Q <- Qp
    NNIs <- .uncompressTipLabel(nni(curr))
    curr <- list(curr)
    class(curr) <- "multiPhylo"
    NNIs <- c(NNIs,curr)
    # Takes forever here --> Multi core helps a bit.
    print("Computing NNI")
    # Utilize doMC to paralellize
    if (domain == 'broadinstitute.org'){
        # Scales in number of edges + # NNI to calc.
        NNIs <- lapply(NNIs, nnls.tree, dm=D, rooted=ultrametric, trace=0)
    } else {
        NNI_data= foreach(i=1:length(NNIs)) %dopar% {
            nnls.tree(dm=D, tree=NNIs[[i]], rooted=ultrametric, trace=0)
        }
        NNIs = NNI_data
    }

    print("Calculating Qs on all NNIs")
    qs <- sapply(NNIs, phytools:::rss, D=D)
    ii <- which(qs==min(qs))[1]
    if(!quiet) message(paste("Best Q =",qs[ii]))
    Qp <- qs[ii]
    curr <- NNIs[[ii]]
    ct <- ct+1
    write.tree(curr, file=curr.file)

    # Plot the interim results:
    dend = as.dendrogram(curr)
    lab = labels(dend)
    dend2 = color.dend(dend)
    NCLUST=20
    dend3 = dend2
    NL = length(labels(dend))

    # Plot the tree:
    pdf(paste0(treeimgpref, 'mpars_qdcons.pdf'), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    title(paste0('Quadratic path difference consensus tree from max parsimony trees'))
    dev.off()
}
message(paste("Solution found after",ct, "set of nearest neighbor interchanges."))


# OTHER TREES:
# write.tree(bin_NJ, paste0(treepref, 'nj.tree'))
# write.tree(bin_UPGMA, paste0(treepref, 'upgma.tree'))
# write.tree(bin_optim, paste0(treepref, 'nj_optim_mpars.tree'))
# write.tree(bin_stochastic, paste0(treepref, 'pratchet_stochastic_mpars.tree')

# mtree = minTreeDist(NJ(as.matrix(emat)),
#                     internal.trees, method="quadratic.path.difference")


