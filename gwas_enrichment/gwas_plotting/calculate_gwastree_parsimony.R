#!/usr/bin/R
# ----------------------------------------------
# Calculate optimal gwas tree by max parsimony
# on the set of 2.1M enhancers
# 
# Run as:
# source /broad/compbio/cboix/EPIMAP_ANALYSIS/bin/config_ChromImpute.sh
# qsub -cwd -P compbio_lab -t 1-21 -l h_vmem=45G -l h_rt=48:00:00 -N enh_parstree_chunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "conda activate mv_env > /dev/null; R --slave -f $BINDIR/calculate_gwastree_parsimony.R"
# qsub -cwd -P compbio_lab -t 1-21 -l h_vmem=45G -l h_rt=02:00:00 -N enh_parstree_chunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "conda activate mv_env > /dev/null; R --slave -f $BINDIR/calculate_gwastree_parsimony.R"
# Test:
# qsub -cwd -P compbio_lab -t 1-2 -l h_vmem=20G -l h_rt=01:00:00 -N enh_parstree_chunked -j y -b y -V -r y -o $DBDIR/out/ChromImpute "conda activate mv_env > /dev/null; R --slave -f $BINDIR/calculate_gwastree_parsimony.R"
# ----------------------------------------------
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

# Tree libraries:
library(Matrix)
library(ape)
library(seqinr)
library(phangorn)
library(dendextend)
library(phylogram)
library(parallel)
options(scipen=45) # So we dont get issues writing integers into bedfiles
options(repos='http://cran.rstudio.com/')

# Arguments for chunking enhancer list:
print(gtargs)
if (length(gtargs) > 1){
    chunksize = as.numeric(gtargs[1])
} else {
    # chunksize = 1e5
    chunksize = 2.5e5
}
chunk = as.numeric(Sys.getenv("SGE_TASK_ID"))
# NOTE: Using full dataset needs more than 120G of space, may not even work in R

# Arguments:
usetree = 'enhancers'
tol = 2500  # Plus/minus distance - window for enhancer overlaps
plot.trees = FALSE 
setprefix = paste0(usetree, '_treemethods_')
use.dedup = TRUE

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
matnames = scan('Enhancer_matrix_names.txt', "c")
meta$uqinfo = paste0(1:nrow(meta), '_', meta$info)
if (use.dedup){
    dedup.file = paste0(gtdir, 'dedup_matrix_indices.Rda')
    load(dedup.file)
    mm = mm2 
    rm(mm2)
    dedup.jacc.file = paste0(gtdir, 'dedup_jacc_matrices.Rda')
    load(dedup.jacc.file)
    setprefix = paste0(setprefix, 'dedupped_')
} else {
    fullmatfile = 'Enhancer_H3K27ac_matrix_062619.mtx.gz'
    enhmatfile = 'Enhancer_H3K27ac_matrix_enhonly_062619.mtx.gz'
    enhmargfile = 'Enhancer_H3K27ac_margins_enhonly_062619.tsv.gz'
    enhrdafile = 'Enhancer_H3K27ac_matrix_enhonly_062619.Rda'
    dna.file = 'Enhancer_H3K27ac_enhonly_062619.dna'
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
}
gc()


# ---------------------------
# Load in the relevant chunk:
# ---------------------------
print(paste(chunk, chunksize))
NENH = nrow(mm)
ind = (chunksize * (chunk-1) + 1):min((chunksize * chunk), NENH)
print(range(ind))
chunkpref = paste0(chunksize, "_", chunk)
treepref = paste0(treedir, setprefix, chunkpref, '_')
treeimgpref = paste0(treeimgpref, chunkpref, '_')

# Get matrix, clean up full matrix for space:
sub.mat = as.matrix(mm[ind,])
# rm(mm)
gc()
colnames(sub.mat) = matnames
binseq = as.phyDat(t(sub.mat), type="USER", levels = c(0, 1))
# rm(sub.mat)
gc()

# ----------------------------------------
# Compute preliminary and parsimony trees:
# ----------------------------------------
# Preliminary distance (may want to improve):
dbin.file = paste0(treepref, '_distml_matrix.Rda')
if (!file.exists(dbin.file)){
    dbin <- dist.ml(binseq, model="JC69")
    dbin[is.na(dbin)] = median(dbin, na.rm=T)
    dbin[is.infinite(dbin)] = 2
    gc()
    save(dbin, file=dbin.file)
} else {
    load(dbin.file)
}


# Preliminary trees:
bin_UPGMA <- upgma(dbin)
bin_NJ  <- NJ(dbin)
gc()
write.tree(bin_NJ, paste0(treepref, 'nj.tree'))
write.tree(bin_UPGMA, paste0(treepref, 'upgma.tree'))

# Calculate parsimony scores for preliminary trees:
print(paste("UPGMA:", parsimony(bin_UPGMA, binseq)))
print(paste("NJ:", parsimony(bin_NJ, binseq)))


if (plotting.trees){
    # Plot the end results:
    dend = as.dendrogram(bin_NJ)
    lab = labels(dend)
    dend2 = color.dend(dend)
    NCLUST=20
    dend3 = dend2
    NL = length(labels(dend))

    # Plot the tree:
    pdf(paste0(treeimgpref, 'nj_mpars.pdf'), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    title(paste0('Chunk: ', paste0(range(ind), collapse=' to '), ', max parsimony from NJ on binary activities'))
    dev.off()
}

# ---------------------------------------------
# Starting from NJ (better), optimize parsimony
# ---------------------------------------------
if (chunksize < 5e5){
    op.file = paste0(treepref, 'nj_optim_mpars.tree')
    if (!file.exists(op.file)){
        bin_optim <- optim.parsimony(bin_NJ, binseq)
        write.tree(bin_optim, op.file)
    } else {
        bin_optim = read.tree(op.file) 
    }

    if (plotting.trees){
        # Plot the end results:
        dend = as.dendrogram(bin_optim)
        lab = labels(dend)
        dend2 = color.dend(dend)
        NCLUST=20
        dend3 = dend2
        NL = length(labels(dend))

        # Plot the tree:
        pdf(paste0(treeimgpref, 'nj_optim_mpars.pdf'), width=14.5, height=12, onefile=T)
        plot.new()
        circle_size = unit(1, "snpc") # snpc unit gives you a square region
        pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
        par(omi = gridOMI(), new = TRUE)
        circleplot(dend=dend3, lab=lab)
        upViewport()
        circos.clear()
        title(paste0('Chunk: ', paste0(range(ind), collapse=' to '), ', max parsimony from NJ on binary activities'))
        dev.off()
    }
}

# ------------------------------------------
# Optimize parsimony using parsimony ratchet
# ------------------------------------------
pr.file = paste0(treepref, 'pratchet_mpars.tree')
if (!file.exists(pr.file)){
    # bin_pratchet <- pratchet(binseq, start=bin_NJ, maxit=20, k=5)
    bin_pratchet <- mod.pratchet(binseq, start=bin_NJ, maxit=20, k=5, 
                                 savefile=paste0(pr.file, '.current'), mc.cores=8)
    write.tree(bin_pratchet, pr.file)
} else {
    bin_pratchet = read.tree(pr.file)
}


# Get branch lengths:
blpr.file = paste0(treepref, 'bl_pratchet_mpars.tree')
if (!file.exists(blpr.file)){
    bl_pratchet = acctran(bin_pratchet, binseq)
    write.tree(bl_pratchet, blpr.file)
} else {
    bl_pratchet = read.tree(blpr.file)
}


if (plotting.trees){
    # Plot the end results:
    dend = as.dendrogram(bin_pratchet)
    lab = labels(dend)
    dend2 = color.dend(dend)
    NCLUST=20
    dend3 = dend2
    NL = length(labels(dend))

    # Plot the tree:
    pdf(paste0(treeimgpref, 'pratchet_mpars.pdf'), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    title('First 10k enhancers - max parsimony on binary activities')
    dev.off()


    # Plot the end results:
    dend = as.dendrogram(bl_pratchet)
    lab = labels(dend)
    dend2 = color.dend(dend)
    NCLUST=20
    dend3 = dend2
    NL = length(labels(dend))

    # Plot the tree:
    pdf(paste0(treeimgpref, 'bl_pratchet_mpars.pdf'), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    title('Max parsimony on binary activities + branch lengths')
    dev.off()

}


stpr.file = paste0(treepref, 'pratchet_stochastic_mpars.tree')
bin_stochastic <- pratchet(binseq, start=bin_NJ, maxit=20, k=5)
write.tree(bin_stochastic, stpr.file)


# # TODO: Try max-likelihood methods:
# fit <- pml(bin_NJ, binseq)
# print(fit)

# # Optimize by ML, rearranging stochastically:
# fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
# logLik(fitJC)

# bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))


# plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")





# Look at modules:
analyze.modules=FALSE
if (analyze.modules){
    # ---------------------------
    # Load in the relevant chunk:
    # ---------------------------
    load('modules_enhancer_sets.Rda')
    enhind = as.numeric(scan('Enhancer_indices.txt', 'c')) + 1
    enhmap = rep(0, max(enhind))
    enhmap[enhind] = 1:length(enhind)
    for (chunk in 4:300){
        print(paste(chunk))
        ind = enhmap[enhsets[[chunk]]]

        chunkpref = paste0("modules_", chunk)
        treepref = paste0(treedir, setprefix, chunkpref, '_')
        treeimgpref = paste0(imgpref, 'e', tol, '_', chunkpref, '_')

        # Get matrix, clean up full matrix for space:
        sub.mat = as.matrix(mm[ind,])
        colnames(sub.mat) = matnames
        binseq = as.phyDat(t(sub.mat), type="USER", levels = c(0, 1))
        gc()

        # ----------------------------------------
        # Compute preliminary and parsimony trees:
        # ----------------------------------------
        # Preliminary distance (may want to improve):
        dbin <- dist.ml(binseq, model="JC69")
        dbin[is.na(dbin)] = median(dbin, na.rm=T)
        dbin[is.infinite(dbin)] = 2

        # Preliminary trees:
        bin_UPGMA <- upgma(dbin)
        bin_NJ  <- NJ(dbin)
        write.tree(bin_NJ, paste0(treepref, 'nj.tree'))
        write.tree(bin_UPGMA, paste0(treepref, 'upgma.tree'))

        # Calculate parsimony scores for preliminary trees:
        print(paste("UPGMA:", parsimony(bin_UPGMA, binseq)))
        print(paste("NJ:", parsimony(bin_NJ, binseq)))

        # ---------------------------------------------
        # Starting from NJ (better), optimize parsimony
        # ---------------------------------------------
        bin_optim <- optim.parsimony(bin_NJ, binseq)
        write.tree(bin_optim, paste0(treepref, 'nj_optim_mpars.tree'))

        if (plotting.trees){
            # Plot the end results:
            NL = length(labels(bin_optim))
            tree = root(bin_optim, outgroup=labels(bin_optim)[NL], 
                        resolve.root = TRUE)
            dend = as.dendrogram(tree)
            lab = labels(dend)
            dend2 = color.dend(dend)
            NCLUST=20
            dend3 = dend2

            # Plot the tree:
            pdf(paste0(treeimgpref, 'nj_optim_mpars.pdf'), width=14.5, height=12, onefile=T)
            plot.new()
            circle_size = unit(1, "snpc") # snpc unit gives you a square region
            pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
            par(omi = gridOMI(), new = TRUE)
            circleplot(dend=dend3, lab=lab)
            upViewport()
            circos.clear()
            title(paste0('Module ', chunk, ', max parsimony from NJ on binary activities'))
            dev.off()
        }
    }

    # smat = sub.mat[,labels(dend)]
    # image(smat)

    #     dt <- dist(smat, 'eJaccard')
#     ht <- hclust(dt, method='ward.D')
#     cocl <- order.optimal(dt, ht$merge)$order
#     reord <- names(cocl)[cocl]
#     ss = reord(smat)
#     image(ss)


}

