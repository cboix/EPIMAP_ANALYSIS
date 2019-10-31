#!/usr/bin/R
# ---------------------------------------------------
# Preliminary analysis for tree GWAS - 
# Counts on a fixed matrix from epigenomic similarity
# ---------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
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
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(GenomicRanges)
library(ape)
library(dplyr)
# Fast libraries for glm:
library(RcppNumerical)  # for fastLR (fastest)
library(speedglm)  # for speedglm
options(scipen=45) # So we dont get issues writing integers into bedfiles
# options(repos='http://cran.rstudio.com/')

# Load specific distance matrices:
fixedstate = FALSE
if (fixedstate){
    source(paste0(bindir, 'load_region_distance_matrices.R'))
    setprefix = paste0('region_',nregions,'_distances_')
} else { 
    source(paste0(bindir, 'load_distance_matrices.R'))
    setprefix = 'distances_'
}

usetree = 'enhancers'
# usetree = 'gwastree'
# usetree = 'correlation'
tol = 2500
singlematch = FALSE # Use 1-to-1 SNP enhancer mapping only?

today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "gwas_tree_analysis/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, setprefix)
treeimgpref = paste0(imgdir, usetree, '_e', tol, '_')

# ----------------------------------------
# Load in GWAS and overlap with enhancers:
# ----------------------------------------
# Load GWAS Catalog, reduce to loc, pmid, pval:
gwdf = read.delim('gwascatalog_may03_2019_noquotes.txt', header=T, stringsAsFactors=F, sep="\t")
gwdf = gwdf[,c('chrom', 'chromStart', 'chromEnd', 'pubMedID','trait', 'pValue')]
gwdf$uid = paste0(gwdf$pubMedID, ' - ', gwdf$trait)
# Collapse multi-counted SNPs:
gwdf = aggregate(pValue ~ chrom + chromStart + chromEnd + uid + trait + pubMedID, gwdf, min)
gwgr = GRanges(paste0('chr', gwdf$chrom), IRanges(gwdf$chromStart, gwdf$chromEnd))

# Totals:
aggw = aggregate(pValue ~ uid + trait, gwdf, length)

# Load in the enhancer coordinates:
ddir = 'DHS_Index_WM201902/'
dmlfile = paste0(ddir, 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt')
dmlnamfile = paste0(ddir, 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv')
dmldf = read.table(dmlfile, header=F, stringsAsFactors=F, sep="\t")
names(dmldf) = c('chr','start','end','name')
# Reorder dml:
dmlnam = read.table(dmlnamfile, header=T, stringsAsFactors=F, sep="\t")
dmldf = merge(dmldf, dmlnam)
dmldf = dmldf[order(dmldf$cls), ]
rm(dmlnam)

# Enhancer indices are 0-indexed - turn to 1-index:
enhind = as.numeric(scan('Enhancer_indices.txt', 'c')) + 1
enhdf = data.frame(cls = enhind)
enhdf = merge(enhdf, dmldf)
rm(dmldf)

# Overlap SNPs with enhancers, with some tolerance:
# -------------------------------------------------
dmgr = GRanges(enhdf$chr, IRanges(enhdf$start - tol, enhdf$end + tol))
qdf = suppressWarnings(data.frame(findOverlaps(gwgr, dmgr)))
qdf$uid = gwdf$uid[qdf$queryHits]
# Non-unique and unique counts:
qnumdf = aggregate(queryHits ~ uid,qdf, function(x){length(unique(x))})
qallnumdf = aggregate(queryHits ~ uid,qdf, function(x){length(x)})
# Almost all SNPs are within an enhancer +/- tolerance:
print(round(length(unique(qdf$queryHits)) / length(gwgr) * 100, 2))

# 1-1 assign each SNP to nearest ENHANCER:
# ----------------------------------------
dmgr2 = GRanges(enhdf$chr, IRanges(enhdf$start, enhdf$end))
nl = nearest(gwgr, dmgr2)
ndf = data.frame(queryHits=1:length(gwgr), subjectHits=nl)
# Filter to keep ones in tolerance:
ndf = merge(ndf, qdf)
ndf = ndf[order(ndf$queryHits),]
ndf$uid = gwdf$uid[ndf$queryHits]
nnumdf = aggregate(queryHits ~ uid, ndf, function(x){length(unique(x))})
nallnumdf = aggregate(queryHits ~ uid, ndf, function(x){length(x)})
# Keep ~ 99k of 121k SNPs
print(round(length(unique(ndf$queryHits)) / length(gwgr) * 100, 2))


# --------------------------------------------
# Load enhancer jaccard - see different trees:
# --------------------------------------------
emat = read.delim(gzfile('Enhancer_jaccard.tsv.gz'), sep="\t", header=F)
matnames = scan('Enhancer_matrix_names.txt', "c")
rownames(emat) = matnames
colnames(emat) = matnames

dt <- as.dist(emat)
# ht <- hclust(dt, method='ward.D')
method = 'ward.D'
method = 'ward.D2'
ht <- hclust(dt, method=method)
cocl <- order.optimal(dt, ht$merge)$order
tree = as.phylo(ht)

dend = as.dendrogram(tree)
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

# NOTE: NEED CIRCLEPLOT FUN.
plot.trees=FALSE
if (plot.trees){
    NCLUST=20
    pdf(paste0(imgpref,sub("\\.","_",method),"_link_jacc.pdf"), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x=circle_size, just="left")
    title(paste0(method, '-linkage of Jaccard Similarity of Enhancers'))
    dev.off()
}

# ----------------------------------
# Load gwas matrix for tree creation
# ----------------------------------
gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_5000_enrich.tsv'
filepref = 'cls_merge2_wH3K27ac100_raw'
gwlindf = read.delim(gwasfile, header=F)
names(gwlindf) = c('pvalue','cluster','pmid','trait',
              'counthit','countall','fold')
namesfile = paste0(filepref, '_names.tsv')
epinames = scan(namesfile,'c')
gwlindf$pmt = paste0(gwlindf$pmid, '_', gwlindf$trait)
gwlindf$cls = paste0('c', gwlindf$cluster)
gwlindf$logpval = -log10(gwlindf$pvalue)
gwlong = aggregate(logpval ~ cls + pmt, gwlindf, max)
wide = spread(gwlong, pmt, logpval, fill=0)
gwmat = as.matrix(wide[,-1])
rownames(gwmat) = wide$cls
gwmat[gwmat < 1] = 0

# Threshold for plotting:
zmax=12
zmin=2
gwmat[gwmat > zmax] <- zmax
gwmat[gwmat < zmin] <- 0

clsn = paste0('c',1:length(epinames) - 1)
names(epinames) = clsn
epimat = gwmat
rownames(epimat) = epinames[rownames(gwmat)]

metric='euclidean'
dt = dist(epimat, method=metric)
if (metric == 'jaccard') { dt = dist(epimat > 0, method=metric) }
# method = 'complete'
method = 'ward.D'
ht <- hclust(dt, method=method)
ht$order <- order.optimal(dt, ht$merge)$order
tree = as.phylo(ht)

dend = as.dendrogram(tree)
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

if (plot.trees){
    NCLUST=20
    pdf(paste0(imgpref,metric, "_", sub("\\.","_",method),"_link_gwmat.pdf"), width=14.5, height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab)
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x=circle_size, just="left")
    title(paste0(method, '-linkage on ', metric, ' distance of GWAS enrichments'))
    dev.off()
}

# ------------------------------------
# Get tree from fused distance matrix:
# ------------------------------------
if (usetree == 'correlation'){
    print("MARK CORRELATION")
    print(usetree)
    dt <- as.dist(full)
    method='ward.D'
} else if (usetree == 'enhancers'){
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
tree = as.phylo(ht)

# ----------------------------------------
# Load matrix and get epigenomes per node:
# ----------------------------------------
# Load in the product of enhancer and H3K27ac mtx: 
mat = read.delim(gzfile('Enhancer_H3K27ac_matrix_062619.mtx.gz'), sep="\t", header=F)
# Matrix is 0-indexed - turn to 1-indexing
mat[,1] = mat[,1] + 1
mat[,2] = mat[,2] + 1
names(mat) = c('row','col')
# Keep only enhancers:
kid = which(mat$row %in% enhind)
mat = mat[kid,]
rm(kid)

# -------------------------
# Make the dendextend tree:
# -------------------------
dend = as.dendrogram(tree)
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
labmapping = sapply(lab, function(x){which(matnames == x)})

# ------------------------------------
# Recursively fill out the tree nodes:
# ------------------------------------
cdllfile = paste0('consensus_object_', usetree, '_062819.Rdata')
if (!file.exists(cdllfile)){
    # Recursive function to get consensus and merge:
    get_consensus <- function(subdend, cdll){
        node = attributes(unclass(subdend))$nodePar$pch
        # GET ID OF NODE:
        print(node)
        nset = cdll$cons[[node]]
        print(head(nset))
        # if (length(nset) == 0){ cdll$cons[[node]] = NA }
        if (length(nset) == 0 || is.na(nset)){
            if (length(subdend) == 2){
                # If internal node, get consensus by merging:
                # Get subtrees and their node #s:
                dend1 = subdend[[1]]
                dend2 = subdend[[2]]
                d1 = attributes(unclass(dend1))$nodePar$pch
                d2 = attributes(unclass(dend2))$nodePar$pch
                # Update each node:
                cdll = get_consensus(dend1, cdll)
                cdll = get_consensus(dend2, cdll)
                # Merge nodes for consensus (intersect) or union:
                cdll$cons[[node]] = intersect(cdll$cons[[d1]], cdll$cons[[d2]])
                cdll$union[[node]] = union(cdll$union[[d1]], cdll$union[[d2]])
                # Update each of the descendants - diff is in node, not parent
                cdll$diff[[d1]] = setdiff(cdll$cons[[d1]], cdll$cons[[node]])
                cdll$diff[[d2]] = setdiff(cdll$cons[[d2]], cdll$cons[[node]])
                # Update: novel is in parent, not node
                cdll$novel[[d1]] = setdiff(cdll$union[[node]], cdll$union[[d1]])
                cdll$novel[[d2]] = setdiff(cdll$union[[node]], cdll$union[[d1]])
            } else {
                print('leaf')
                # If leaf, get epigenome from matrix:
                id = labmapping[labels(subdend)]
                cdll$cons[[node]] = mat[mat[,2] == id, 1]
                cdll$union[[node]] = mat[mat[,2] == id, 1]
            }
        } else { print("Already have consensus at this node") }
        print(length(cdll$cons[[node]]))
        return(cdll)
    }
    # Building from the bottom, fill in all:
    clist = sapply(rep(NA, NN), list)
    cdll = list(cons=sapply(rep(NA, NN), list),
                diff=sapply(rep(NA, NN), list),
                union=sapply(rep(NA, NN), list),
                novel=sapply(rep(NA, NN), list))
    # Store node id as pch (use to ID where we are)
    set(dend, 'nodes_pch', 1:NN) -> dend
    # Run consensus function to update clist:
    cdll = get_consensus(dend, cdll)
    save(cdll, file=cdllfile)
} else {
    print("[STATUS] Loading cdll from file")
    load(cdllfile)
}

# Set top node difference to top node consensus:
if (is.na(cdll$diff[[1]])){
    print("toplevel is na")
    cdll$diff[[1]] = cdll$cons[[1]]
}

if (is.na(cdll$novel[[1]])){
    print("toplevel is na")
    cdll$novel[[1]] = cdll$union[[1]]
}

# Look at lengths - if no 1 (NA), ok.
cdlenlist = lapply(cdll, function(x){sapply(x,length)})
lapply(cdlenlist, function(x){head(x, 10)})
# dlen = sapply(cdll$diff, length)
# clen = sapply(cdll$cons, length)
# ulen = sapply(cdll$union, length)
# nlen = sapply(cdll$novel, length)

# Plot set sizes (vs. height in tree?)
# Only a 8 internal nodes with 0, and more than 600 have 10k+ elements
# NOTE: Numbers are slightly inflated - contain ~ 200k non-enhancer elements (for if we want to look at promoters??
# NOTE: Why some of following leaves with VERY few enhancers??
# c2 = sort(sapply(clist, length)[1:NL], decreasing=T)

# Enhancer index map to matrix:
# (did intersection on enhancers only)
enhmap = rep(0, max(enhind))
enhmap[enhind] = 1:length(enhind)

print(max(qdf$subjectHits))
print(max(ndf$subjectHits))

# ----------------------------
# From intersection, get SNPs:
# ----------------------------
if (singlematch){
    midpref = paste0('_', usetree, '_e', tol, '_single_')
} else {
    midpref = paste0('_', usetree, '_e', tol, '_all_')
}
treefile = paste0('treedf', midpref, 'diff_snpint.tsv')
dflist = list()
types = c('diff','cons','union','novel')
if (!file.exists(treefile)){
    for (type in types){ dflist[[type]] = c() }
    for (i in 1:NN){
        print(i)
        for (type in types){
            enhset = cdll[[type]][[i]]
            # NOTE: Map enhset to enhind
            enhset = enhmap[enhset]
            if (singlematch){
                enhsnp = merge(ndf, data.frame(subjectHits = enhset))$queryHits
            } else {
                enhsnp = merge(qdf, data.frame(subjectHits = enhset))$queryHits
            }
            # Old version: unique snps (seems to artificially shift pvalues to edges of trees due to overlaps)
            # enhsnp = unique(enhsnp)
            if (length(enhsnp) > 0){
                dfsnp = aggregate(pValue ~ uid, gwdf[enhsnp, ], length)
                names(dfsnp) = c('uid','nsnp')
                dfsnp$node = i
                dflist[[type]] = rbind(dflist[[type]], dfsnp)
            }
        }
    }
    # Write out:
    for (type in types){ 
        print(paste("[STATUS] Writing for",type))
        tfile = paste0('treedf', midpref, type, '_snpint.tsv')
        write.table(dflist[[type]], tfile, , quote=F, sep="\t", row.names=F)
    }
} else { 
    for (type in types){ 
        print(paste("[STATUS] Writing for",type))
        tfile = paste0('treedf', midpref, type, '_snpint.tsv')
        dflist[[type]] = read.delim(treefile, sep="\t", header=T)
    }
}

# Remove mat and list for space:
if (nrow(dflist[['diff']])  > 0){
    gc()
    rm(mat)
    rm(cdll)
    gc()
}

# =======================================
# Plot gwas intersection results on tree:
# =======================================
# ------------------------
# Tree plotting functions:
# ------------------------
# Reduced circleplot FN:
circleplot = function(dend, lab, with.tree=TRUE, nodedf=NULL, udf=NULL,
                      # For plotting metadata:
                      add.metadata=TRUE,
                      # For plotting tracks around:
                      fractional=TRUE, fulltrack=TRUE){
    # Set up dendrogram:
    NLAB = length(labels(dend))
    mll = meta.image(metamat[lab,5:1], colvals=colvals, cex=.6, horiz=TRUE, return.mat=TRUE)
    nummat = mll[[1]]
    colsmeta = mll[[2]]
    colmat = matrix(colsmeta[nummat], nrow=nrow(nummat))
    labels_dend <- labels(dend)
    if (as.logical(anyDuplicated(labels_dend))) {
        labels(dend) <- paste0(seq_along(labels_dend), "_", labels_dend)
        labels_dend <- labels(dend) }
    # Plot figure:
    circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
    circos.initialize(factors = "single", xlim = c(0, NLAB)) 
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
                     circos.text(1:NLAB-0.5, rep(0, NLAB), labels(dend), col = labels_colors(dend), cex=.25,
                                 facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) }, bg.border = NA, track.height = 0.1)
    if (!is.null(udf)){
        # Add barplots:
        nr = nrow(udf)
        udf$col = as.character(udf$col)
        if (fulltrack){
            topy = max(udf$total)
            circos.track(ylim = c(0, topy), bg.border = NA, panel.fun = function(x, y) {
                             for(i in 1:nr) {
                                 circos.rect(rep(i - 1, 2),c(0, udf$uq[i]),
                                             rep(i, 2),c(udf$uq[i], udf$total[i]),
                                             col=c(udf$col[i], 'grey85'), border=NA)
                     } }, track.height=0.06) 
        }
        if (fractional){
            circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
                             for(i in 1:nr) {
                                 circos.rect(rep(i - 1, 2),c(0, udf$frac[i]),
                                             rep(i, 2),c(udf$frac[i], 1),
                                             col=c(udf$col[i], 'grey85'), border=NA)
                     } }, track.height=0.06) 
        } 
    }
    if (add.metadata){
        circos.track(ylim = c(0, 5), bg.border = NA, panel.fun = function(x, y) {
                         m = t(colmat[,5:1])
                         nr = nrow(m)
                         nc = ncol(m)
                         for(i in 1:nr) {
                             circos.rect(1:nc - 1, rep(nr - i, nc), 
                                         1:nc, rep(nr - i + 1, nc), 
                                         border = NA, col = m[i, ])
                     } }, track.height=0.075) 
    }
    max_height = max(attr(dend, "height"))
    if (with.tree){
        circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                         circos.dendrogram(dend, max_height = max_height)
                         if (!is.null(nodedf)){
                             circos.text(nodedf$V1, max_height - nodedf$V2, as.character(nodedf$node),
                                         col = 'black', cex=.25,
                                         facing = "reverse.clockwise", niceFacing=TRUE, adj = c(0, 0.5)) 
                         } }, track.height = 0.5, bg.border = NA)
        circos.clear()
    }
    circos.clear()
}

# Same as above, but split:
circleplot_split = function(dend, lab, NCLUST, with.tree=TRUE, only.diff=FALSE, nodedf=NULL){
    # Set up dendrogram:
    NLAB = length(labels(dend))
    mll = meta.image(metamat[lab,5:1], colvals=colvals, cex=.6, horiz=TRUE, return.mat=TRUE)
    nummat = mll[[1]]
    colsmeta = mll[[2]]
    colmat = matrix(colsmeta[nummat], nrow=nrow(nummat))
    labels_dend <- labels(dend)
    if (as.logical(anyDuplicated(labels_dend))) {
        labels(dend) <- paste0(seq_along(labels_dend), "_", labels_dend)
        labels_dend <- labels(dend) }
    # Plot figure:
    # Make split dendrogram:
    colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(NCLUST)
    dend <- color_branches(dend, k=NCLUST, col=colpair)
    dend = set(dend, "labels_cex", .18)
    labels_dend <- labels(dend)
    if (as.logical(anyDuplicated(labels_dend))) {
        labels(dend) <- paste0(seq_along(labels_dend), "_", labels_dend)
        labels_dend <- labels(dend)
    }
    # Factors/splits:
    kcol = unlist(get_leaves_attr(dend, 'edgePar'))
    letext = c(letters, paste0('a', letters))
    factors = letext[as.numeric(factor(kcol, levels = unique(kcol)))]
    dend_list = get_subdendrograms(dend, k=NCLUST)
    # Correct the dendrogram order:
    dnames = sapply(dend_list, function(x){as.numeric(strsplit(labels(x)[1],'_')[[1]][[1]]) })
    names(dend_list) = letext[1:NCLUST][order(order(dnames))]
    # ----------------
    # Plot dendrogram:
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
                 } }, track.height=0.125)
    # Dendrogram
    max_height = max(sapply(dend_list, function(x) attr(x, "height")))
    circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.3, 
                 panel.fun = function(x, y) {
                     sector.index = CELL_META$sector.index
                     dend = dend_list[[sector.index]]
                     circos.dendrogram(dend, max_height = max_height)
                     if (!is.null(nodedf)){
                         circos.text(nodedf$V1, max_height - nodedf$V2, as.character(nodedf$node),
                                     col = 'black', cex=.25,
                                     facing = "reverse.clockwise", niceFacing=TRUE, adj = c(0, 0.5)) 
                     }
                 }
    )
}

# ----------------
# Get descendants:
# ----------------
get_dec <- function(subdend, declist, parent=0){
    node = attributes(unclass(subdend))$nodePar$pch
    # GET ID OF NODE:
    nset = declist$dec[[node]]
    # if (length(nset) == 0){ declist$dec[[node]] = NA }
    if (length(nset) == 0 || is.na(nset)){
        if (length(subdend) == 2){
            # Get subtrees and their node #s:
            dend1 = subdend[[1]]
            dend2 = subdend[[2]]
            d1 = attributes(unclass(dend1))$nodePar$pch
            d2 = attributes(unclass(dend2))$nodePar$pch
            # Update each node:
            declist = get_dec(dend1, declist, parent=node)
            declist = get_dec(dend2, declist, parent=node)
            # Merge each side's descendants:
            declist$dec[[node]] = c(declist$dec[[d1]], declist$dec[[d2]])
            declist$isleaf[node] = 0
        } else {
            # If leaf, get epigenome from matrix:
            slab = labels(subdend)
            declist$dec[[node]] = slab
            declist$isleaf[node] = 1
        }
        declist$parent[node] = parent
    } else { print("Already have consensus at this node") }
    return(declist)
}


labels_dend <- labels(dend3)
if (as.logical(anyDuplicated(labels_dend))) {
    labels(dend3) <- paste0(seq_along(labels_dend), "_", labels_dend)
    labels_dend <- labels(dend3) }

# For counting:
dend4 = set(dend3, 'nodes_pch', 1:NN)
declist = list(dec=sapply(rep(NA, NN), list), isleaf=rep(NA, NN), parent=rep(NA, NN))
declist = get_dec(dend4, declist=declist)
declist$parent[1] = 1
pdf = data.frame(node=1:length(declist$parent), parent=declist$parent)

# Mapping: dendlist
mapdl = data.frame(id=labels(dend), lab=labels(dend3))
mapnode = data.frame(node=which(declist$isleaf == 1), lab=sapply(which(declist$isleaf == 1), function(x){declist$dec[[x]]}))
mapmat = data.frame(id=names(labmapping), matid=labmapping)
mapnode = merge(merge(mapdl, mapnode), mapmat)

# Get number of unique and total enh:
nl = which(declist$isleaf == 1)
dlen = cdlenlist[['diff']]
clen = cdlenlist[['cons']]
fracuq = dlen[nl] / clen[nl]
udf = data.frame(uq=dlen[nl], total = clen[nl], frac=fracuq)
udf$col='grey65'

# Plot basic tree - verify working:
pdf(paste0(treeimgpref, 'basic_tree.pdf'), width=14.5, height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend3, lab=lab, udf=udf)
upViewport()
draw(pd.legend, x = circle_size, just = "left")
dev.off()

# -----------------------------------------
# Plot the NUMBER of regions at each point:
# -----------------------------------------
# DIFFERENCE:
# Color by amount (viridis?)
palette = colryb
dr = dlen / 1000
bins <- cut(dr, seq(0, max(dr), length.out=length(palette)), include.lowest=T) 
dc = palette[bins]
dend3 = set(dend3, 'nodes_pch', 19)
dend3 = set(dend3, 'nodes_cex', .5)
dend3 = set(dend3, 'nodes_col', dc)

pdf(paste0(treeimgpref, 'numdiff_regions.pdf'), width=14.5,height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend3, lab=lab)
upViewport()
title("Number of sites differing from parent")
draw(pd.legend, x = circle_size, just = "left")
dev.off()

# Color the number of CONSENSUS:
palette = colryb
cr = clen / 1000
bins <- cut(cr, seq(0, max(cr), length.out=length(palette)), include.lowest=T) 
cc = palette[bins]
dend3 = set(dend3, 'nodes_pch', 19)
dend3 = set(dend3, 'nodes_cex', .5)
dend3 = set(dend3, 'nodes_col', cc)

pdf(paste0(treeimgpref,'numcons_regions.pdf'),width=14.5,height=12)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circleplot(dend3, lab=lab)
upViewport()
title("Number of consensus sites")
draw(pd.legend, x = circle_size, just = "left")
dev.off()


# ----------------------
# Look at specific GWAS:
# ----------------------
run.hyper <- function(y){phyper(q=y[1]-1, m=y[3], n=y[4] - y[3], k=y[2], lower.tail=FALSE)} 
run.fisher <- function(y){
    table <-  matrix(as.numeric(c(y[1], y[2], y[3], y[4])), ncol = 2, byrow = TRUE)
    if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
    return(p)
}

# Update palette and legend for pvalues:
palette = col3
CUTP = 12
col_fun = colorRamp2(c(0, CUTP), c(palette[1], palette[length(palette)]))
pval.legend = Legend(at = seq(3, CUTP, 3), col_fun=col_fun, title_position = "topleft", title = "p-value")
pd.legend.ext = packLegend(col.list$group, col.list.h$project, col.list$sex, col.list$type, col.list$lifestage, pval.legend)

# Q: Bottom-up: Is diff signif vs. consensus?
# # Pvalues for difference node-parent vs. consensus node:
get_pvalues = function(strait, dflist, cdlenlist, pdf=pdf, cutp=CUTP, 
                       ingroup='diff', outgroup='cons', against.parent=TRUE){
    # Find the GWAS uid with most snps:
    sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == strait,], length)
    suid = sgw[order(sgw$pValue, decreasing=T), 'uid'][1]
    # Measure numbers of ingroup SNPs (in any group):
    isubdf = filter(dflist[[ingroup]], uid == suid)
    isnp = rep(0, NN)
    isnp[isubdf$node] = isubdf$nsnp
    ilen = cdlenlist[[ingroup]]
    # Make title:
    title = paste('Node', capitalize(ingroup),'vs.')
    # Measure numbers of outgroup SNPs (in any group + enh):
    if (outgroup == 'enh'){
        if (singlematch){
            # osnp = filter(nnumdf, uid == suid)$queryHits
            osnp = filter(nallnumdf, uid == suid)$queryHits
        } else {
            # osnp = filter(qnumdf, uid == suid)$queryHits
            osnp = filter(qallnumdf, uid == suid)$queryHits
        }
        olen = length(enhind)
        title = paste(title, 'All Enhancers')
    } else { 
        osubdf = filter(dflist[[outgroup]], uid == suid)  # SNPS in novel 
        osnp = rep(0, NN)
        osnp[osubdf$node] = osubdf$nsnp
        olen = cdlenlist[[outgroup]]
        if (against.parent){
            # Reorder snps and lengths by parent:
            osnp = osnp[declist$parent] 
            olen = olen[declist$parent]
            title = paste(title, 'Parent', capitalize(outgroup))
        } else {
            title = paste(title, 'Node', capitalize(outgroup))
        }
    }
    # Hyper-geometric testing:
    df = cbind(q=isnp, draw=ilen,
               m=osnp, N=olen)
    pout <- apply(df, 1, run.hyper)
    df = data.frame(df)
    df$pout = pout
    df$padj = p.adjust(df$pout)
    df$padj[df$padj == 0] = min(df$padj[df$padj != 0])
    df$rawlog10p = -log10(df$padj)
    # Cap and threshold the adjusted log10p:
    df$log10p = df$rawlog10p
    df$log10p[df$log10p < 3] = 0
    df$log10p[df$log10p > cutp] = cutp
    return(list(log10p=df$log10p, isnp=isnp, rawlp=df$rawlog10p, df=df, title=title))
}


setup_dendrogram <- function(dend, ll, udf, declist=declist, altline=3){
    log10p = ll$log10p
    # Color for strength of p-value:
    bins <- cut(log10p, seq(0, CUTP, length.out=length(palette)), include.lowest=T) 
    cc = palette[bins]
    # Highlight the top hits:
    sord = order(ll$rawlp, ll$isnp, decreasing=T)
    tophits = head(sord, 4)
    cexlist = .5 * (log10p > 0)
    cexlist[tophits] = 1
    pchlist = rep(19,NN)
    pchlist[tophits] = 21
    dend = set(dend, 'nodes_pch', pchlist)
    dend = set(dend, 'nodes_cex', cexlist)
    dend = set(dend, 'nodes_col', cc)
    # For counting:
    dend4 = set(dend, 'nodes_pch', 1:NN)
    declist = list(dec=sapply(rep(NA, NN), list), isleaf=rep(NA, NN), parent=rep(NA, NN))
    declist = get_dec(dend4, declist=declist)
    # Get all descendant branches:
    cutoff = 6
    keeplp = which((log10p >= cutoff) * (1 - declist$isleaf) == 1) # Remove if "only" leaf
    leafhits = (log10p >= cutoff)
    hits = unique(unlist(sapply(keeplp, function(x){declist$dec[[x]]})))
    # NOTE: Could set altline to either 3 or 0 (blank - faster plotting)
    dend = set(dend, "by_labels_branches_lty", value=hits, TF_values = c(1, altline))
    dend = set(dend, "by_labels_branches_lwd", value=hits, TF_values = c(1, 1))
    # Edit colors:
    nl = which(declist$isleaf == 1)
    udf$col = 'grey65'
    udf$col[which(nl %in% which(leafhits == 1))] = 'red4'
    return(list(dend=dend, hits=leafhits, udf=udf))
}

# ---------------------
# Plot basic dendrogram
# ---------------------
# List of traits:
ut = unique(gwdf$trait)
# Choose one gwas: 
# strait = 'Age-related macular degeneration'
# strait = 'Schizophrenia'
# strait = "Alzheimer's disease"
# strait = 'Type 2 diabetes'
strait = 'Glomerular filtration rate'
strait = 'HDL cholesterol'
strait = 'LDL cholesterol'
strait = "Crohn's disease"
# strait = "Alzheimer's disease"

ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                 ingroup='diff', outgroup='enh', against.parent=FALSE)
ldf = ll$df[order(ll$df$pout),]
head(ldf,50)

dd  = setup_dendrogram(dend3, ll, udf, declist=declist)

# Plot basic dendrogram
if(usetree == 'enhancers'){

    testlist = list(c('union','union',TRUE),
                    # c('novel','union',TRUE),
                    c('union','enh',FALSE),
                    # c('novel','enh',FALSE),
                    # c('diff','enh', FALSE),
                    c('diff','cons',FALSE),
                    c('cons','enh', FALSE))

    # Pick four semi-obvious GWAS:
    testtraits = c('Macular thickness', "QT interval", 'LDL cholesterol', "Crohn's disease")
    pdf(paste0(treeimgpref,'full_test_types.pdf'),width=13,height=12, onefile=T)
    for (strait in testtraits) {
        for (testset in testlist){
            print(testset)
            ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                             ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
            dd  = setup_dendrogram(dend3, ll, udf, declist=declist, altline=0)
            plot.new()
            circle_size = unit(1, "snpc") # snpc unit gives you a square region
            pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                                  just = c("left", "center")))
            par(omi = gridOMI(), new = TRUE)
            circleplot(dend=dd$dend, lab=lab, udf=dd$udf, add.metadata=FALSE)
            upViewport()
            circos.clear()
            # draw(pd.legend.ext, x = circle_size, just = "left")
            title(paste(strait, '-', ll$title))
        }
    }
    dev.off()


    pdf(paste0(treeimgpref,'full_test.pdf'),width=14.5,height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dd$dend, lab=lab, udf=dd$udf)
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x = circle_size, just = "left")
    title(paste(strait, '-', ll$title))
    dev.off()


    # Plot consensus:
    # NOTE: Can try both treedf and consdf:
    # ll = get_pvalues(strait, treedf)
    ll = get_pvalues(strait, treedf=dflist[['diff']], consdf=dflist[['cons']],
                     against.enh=FALSE, against.parent=FALSE)
    dd  = setup_dendrogram(dend3, ll)
    dend3 = dd$dend
    # Edit colors:
    nl = which(declist$isleaf == 1)
    udf$col = 'grey65'
    udf$col[which(nl %in% which(dd$hits == 1))] = 'red4'
    # Plot:
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab, udf=udf)
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x = circle_size, just = "left")
    title(paste0(strait, ' - Branch Difference vs. Subtree Consensus'))
    # PARENT AGAINST NODE:
    ll = get_pvalues(strait, treedf=dflist[['diff']], consdf=dflist[['cons']], against.enh=FALSE, against.parent=TRUE)
    dd  = setup_dendrogram(dend3, ll)
    dend3 = dd$dend
    # Edit colors:
    nl = which(declist$isleaf == 1)
    udf$col = 'grey65'
    udf$col[which(nl %in% which(dd$hits == 1))] = 'red4'
    # Plot:
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dend3, lab=lab, udf=udf)
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x = circle_size, just = "left")
    title(paste0(strait, ' - Parent Consensus vs. Child Consensus'))
    dev.off()
} else {
    NCLUST=20
    pdf(paste0(treeimgpref,'split_k', NCLUST, '_test.pdf'),width=14.5,height=12, onefile=T)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot_split(dend=dend3, lab=lab, NCLUST=NCLUST)
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x=circle_size, just="left")
    title(strait)
    dev.off()
}


# ---------------------------
# Plot many different traits:
# ---------------------------
# List of traits to test:
tlist = c(# Liver:
          'HDL cholesterol', 'LDL cholesterol', 'Triglycerides',
          # Mostly immune:
          'Type 1 diabetes', 'Type 2 diabetes', "Crohn's disease", 'Rheumatoid arthritis',
          'Multiple sclerosis', 'Systemic lupus erythematosus', "Inflammatory bowel disease",
          # Brain:
          'Neuroticism', "Alzheimer's disease", "Parkinson's disease",
          'Schizophrenia', 'Age-related macular degeneration',
          'Macular thickness', 'Glaucoma',
          # Heart:
          "QT interval", "Resting heart rate", 
          # Other:
          'Migraine', 'Proinsulin levels', "Ulcerative colitis", 
          'Glomerular filtration rate', "Pulmonary function")


expandedlist = c(# Lung:
                 "Post bronchodilator FEV1/FVC ratio", "Lung function (FEV1/FVC)", "Post bronchodilator FEV1", "Pulmonary function",
                 # Bone/body
                 "Heel bone mineral density", "Total body bone mineral density", "Adolescent idiopathic scoliosis", "Height",
                 "Hair color", "Body mass index", "Waist-hip ratio", "Waist-to-hip ratio adjusted for body mass index",
                 "Fat-free mass", "Obesity-related traits", "Total cholesterol levels",
                 # Heart:
                 "QT interval", "Resting heart rate", "Coronary artery disease", "Cardiovascular disease", "Atrial fibrillation",
                 # Brain:
                 "General cognitive ability", "Educational attainment (years of education)", "Educational attainment (MTAG)",
                 "Intelligence (MTAG)", "Intelligence", "Depressive symptoms", "Well-being spectrum (multivariate analysis)",
                 "General risk tolerance (MTAG)", "Cognitive performance (MTAG)", "Highest math class taken (MTAG)", "Self-reported math ability (MTAG)",
                 "Self-reported math ability", "Smoking initiation (ever regular vs never regular) (MTAG)", "Itch intensity from mosquito bite adjusted by bite size",
                 "Reaction time", "Chronotype", "Autism spectrum disorder or schizophrenia", "Alzheimer's disease or family history of Alzheimer's disease",
                 # Brain:
                 'Neuroticism', "Alzheimer's disease", "Parkinson's disease", 'Schizophrenia', 
                 # Eye: 
                 'Age-related macular degeneration', 'Macular thickness', 'Glaucoma', "Intraocular pressure",
                 # Liver:
                 'HDL cholesterol', 'LDL cholesterol', 'Triglycerides',
                 "Blood metabolite levels", "Serum metabolite ratios in chronic kidney disease", 
                 # Sex/age covariate related?
                 "Breast cancer", "Menarche (age at onset)", "Male-pattern baldness", "Balding type 1",
                 "DNA methylation variation (age effect)",
                 # Mostly immune:
                 'Type 1 diabetes', 'Type 2 diabetes', "Crohn's disease", 'Rheumatoid arthritis', "Ulcerative colitis",
                 'Multiple sclerosis', 'Systemic lupus erythematosus', "Inflammatory bowel disease",
                 # Blood traits: See how refined we get?
                 "Blood protein levels", "Red blood cell count", "White blood cell count",
                 "Monocyte count", "Eosinophil counts", "Platelet count",
                 "Plateletcrit", "Mean corpuscular volume", "Mean corpuscular hemoglobin",
                 "Pulse pressure", "Systolic blood pressure", "Diastolic blood pressure", "IgG glycosylation", 
                 # Other:
                 'Migraine', 'Proinsulin levels', 'Glomerular filtration rate')


# EVAL MANY, reduced size:
pdf(paste0(treeimgpref,'full_examples_gwaslist.pdf'),width=13,height=12, onefile=T)
for (strait in tlist) {
    testset = c('cons','enh', FALSE)
    print(strait)
    ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                     ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
    dd  = setup_dendrogram(dend3, ll, udf, declist=declist, altline=0)  # Plot with none-lines
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dd$dend, lab=lab, udf=dd$udf, add.metadata=FALSE)  # Plot without metadata!
    upViewport()
    circos.clear()
    title(paste(strait, '-', ll$title))
}
dev.off()


NCLUST=10
# ONLY TEST diff/enh
pdf(paste0(treeimgpref,"split_k",NCLUST,"_examples.pdf"), width=14.5, height=12, onefile=T)
for (strait in expandedlist){
    print(strait)
    testlist = list(c('diff','enh', FALSE))
    for (testset in testlist){
        ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                         ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist)
        plot.new()
        circle_size = unit(1, "snpc") # snpc unit gives you a square region
        pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                              just = c("left", "center")))
        par(omi = gridOMI(), new = TRUE)
        circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fulltrack=FALSE)
        upViewport()
        circos.clear()
        # draw(pd.legend.ext, x = circle_size, just = "left")
        title(paste(strait, '-', ll$title))
    }
}
dev.off()


# Test methods:
testlist = list(c('union','union',TRUE), 
                c('diff','enh', FALSE))
# c('novel','union',TRUE)
# c('diff','cons',FALSE)
# c('cons','enh', FALSE)
# Chunk traits:
chunksize=20
nchunk = ceiling(length(expandedlist) / chunksize)
for (i in 1:nchunk){
    pdf(paste0(treeimgpref,"examples_comp_methods", sprintf("%03d",i), ".pdf"), width=14.5, height=12, onefile=T)
    chunklist = expandedlist[((i-1) * chunksize + 1):min(i * chunksize, length(expandedlist))]
    print(i)
    print(chunklist)
    for (strait in chunklist){
        print(strait)
        for (testset in testlist){
            ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                             ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
            dd  = setup_dendrogram(dend3, ll, udf, declist=declist)
            plot.new()
            circle_size = unit(1, "snpc") # snpc unit gives you a square region
            pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                                  just = c("left", "center")))
            par(omi = gridOMI(), new = TRUE)
            circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE)
            upViewport()
            circos.clear()
            draw(pd.legend.ext, x = circle_size, just = "left")
            title(paste(strait, '-', ll$title))
        }
    }
    dev.off()
}



# Look at Triglycerides - do we need to prune, and how:
strait = 'Triglycerides'
sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == strait,], length)
suid = sgw[order(sgw$pValue, decreasing=T), 'uid'][1]
subdf = gwdf[gwdf$uid == suid,]
dim(subdf)
subdf = aggregate(pValue ~ chrom + chromStart + chromEnd + uid + trait, subdf, min)
dim(subdf)

subdf = subdf[order(subdf$chrom, subdf$chromStart),]

# Want to prune - take top snp in +/-1Mb at all times

pdf(paste0(treeimgpref,'full_test_types.pdf'),width=13,height=12, onefile=T)
for (testset in testlist){
    print(testset)
    ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                     ingroup=testset[1], outgroup=testset[2], against.parent=testset[3])
    dd  = setup_dendrogram(dend3, ll, udf, declist=declist)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dd$dend, lab=lab, udf=dd$udf)
    upViewport()
    circos.clear()
    # draw(pd.legend.ext, x = circle_size, just = "left")
    title(paste(strait, '-', ll$title))
}
dev.off()


snps = gwdf[gwdf$uid == suid,]



# ---------------------------------------
# Pilot the logistic regression approach:
# ---------------------------------------
# For exposition, easy trait + top hit:
# Will show pvals with all vs with reduced 1-to-1
strait = 'LDL cholesterol'
# strait = "Crohn's disease"
type = 'cons'

# Trait and example place:
sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == strait,], length)
suid = sgw[order(sgw$pValue, decreasing=T), 'uid'][1]
ll = get_pvalues(strait, dflist, cdlenlist=cdlenlist, pdf=pdf,
                 ingroup='diff', outgroup='enh', against.parent=FALSE)
ldf = ll$df[order(ll$df$pout),]
ldf$node = rownames(ldf)
head(ldf, 5)

i = as.numeric(ldf[ldf$draw > 10000, 'node'][1])
print(paste("Is leaf:", declist$isleaf[i]))
p = declist$parent[i]

# Create vectors:
NE = nrow(enhdf)
snps = rep(0, NE)
allsnps = rep(0, NE)
ann = rep(0, NE)
pann = rep(0, NE)
snpset = ndf[ndf$uid == suid, 'subjectHits']
allsnpset = qdf[qdf$uid == suid, 'subjectHits']
nset = enhmap[cdll[[type]][[i]]]
pset = enhmap[cdll[[type]][[p]]]

ann[nset] = 1
pann[pset] = 1
snps[snpset] = 1
allsnps[allsnpset] = 1
print(paste(sum(ann), sum(pann), sum(snps), sum(allsnps)))

testdf = data.frame(snp=snps, parent=pann, node=ann, allsnp=allsnps, intercept=1)
testdf$diff = testdf$node - testdf$parent

# Compare options by speed: 
# Regression using 1-to-1 snps only:
nmod = speedglm(snp ~ parent + 1, data=testdf, family=binomial())
dmod = speedglm(snp ~ parent + diff + 1, data=testdf, family=binomial())
lrd = summary(dmod)$log - summary(nmod)$log
print(paste(sapply(c(summary(nmod)$log, summary(dmod)$log, lrd), digits=2, round)))
pchisq(lrd, df=1, lower.tail=FALSE)

# Log reg for each (~ 14.36 seconds)
ptm <- proc.time()
nmod = speedglm(allsnp ~ parent + 1, data=testdf, family=binomial())
dmod = speedglm(allsnp ~ parent + diff + 1, data=testdf, family=binomial())
lrd = summary(dmod)$log - summary(nmod)$log
print(paste(sapply(c(summary(nmod)$log, summary(dmod)$log, lrd), digits=2, round)))
pchisq(lrd, df=1, lower.tail=FALSE)
print(proc.time() - ptm)

# Log reg for each (~ 3.5 seconds)
ptm <- proc.time()
nmat = as.matrix(testdf[,c('parent', 'intercept')])
dmat = as.matrix(testdf[,c('parent', 'diff', 'intercept')])
nmodflr = fastLR(x=nmat, y=testdf$allsnp)
dmodflr = fastLR(x=dmat, y=testdf$allsnp)
lrdflr = dmodflr$loglikelihood - nmodflr$loglikelihood
print(paste(sapply(c(nmodflr$loglikelihood, dmodflr$loglikelihood, lrdflr), digits=2, round)))
pchisq(lrdflr, df=1, lower.tail=FALSE)
print(proc.time() - ptm)

# Time comparison for just functions:
system.time(dmod <- speedglm(allsnp ~ parent + diff + 1, data=testdf, family=binomial()))
system.time(dmod <- fastLR(x=dmat, y=testdf$allsnp))



# Function using fastLR:
for (i in 1:NN){
    print(i)
    for (type in types){
        enhset = cdll[[type]][[i]]
        if (singlematch){
            enhsnp = unique(merge(ndf, data.frame(subjectHits = enhset))$queryHits)
        } else {
            enhsnp = unique(merge(qdf, data.frame(subjectHits = enhset))$queryHits)
        }
        if (length(enhsnp) > 0){
            dfsnp = aggregate(pValue ~ uid, gwdf[enhsnp, ], length)
            names(dfsnp) = c('uid','nsnp')
            dfsnp$node = i
            dflist[[type]] = rbind(dflist[[type]], dfsnp)
        }
    }
}



# ---------------------------------------
# Pilot the logistic regression approach:
# ---------------------------------------
# For exposition, easy trait + top hit:
# Will show pvals with all vs with reduced 1-to-1
strait = 'LDL cholesterol'
# strait = "Crohn's disease"
type = 'cons'


# Is slow. can we speed up by not evaluating certain ones?
get_lr_pvalues = function(strait, dflist, cdlenlist, pdf=pdf, cutp=CUTP, 
                       ingroup='diff', outgroup='cons', against.parent=TRUE){
    # Trait and example place:
    sgw = aggregate(pValue ~ uid + trait, gwdf[gwdf$trait == strait,], length)
    suid = sgw[order(sgw$pValue, decreasing=T), 'uid'][1]
    NE = nrow(enhdf)
    pvals = rep(1, NN)
    for (i in 1:NN){
        print(paste(i, "is leaf?", declist$isleaf[i]))
        p = declist$parent[i]
        # Create vectors:
        allsnps = rep(0, NE)
        dann = rep(0, NE)
        pann = rep(0, NE)
        # Sets:
        if (singlematch){
            allsnpset = ndf[ndf$uid == suid, 'subjectHits']
        } else {
            allsnpset = qdf[qdf$uid == suid, 'subjectHits']
        }
        pset = enhmap[cdll[['cons']][[p]]]
        dset = enhmap[cdll[['diff']][[i]]]
        # Fill:
        dann[dset] = 1
        pann[pset] = 1
        allsnps[allsnpset] = 1
        # print(paste(sum(pann), sum(dann), sum(allsnps)))
        testdf = data.frame(snp=snps, parent=pann, node=ann, allsnp=allsnps, intercept=1)
        testdf$diff = testdf$node - testdf$parent
        # Log reg for each (~ 3.5 seconds)
        ptm <- proc.time()
        nmat = cbind(pann, 1)
        dmat = cbind(pann, dann, 1)
        nmodflr = fastLR(x=nmat, y=testdf$allsnp)
        dmodflr = fastLR(x=dmat, y=testdf$allsnp)
        lrdflr = dmodflr$loglikelihood - nmodflr$loglikelihood
        print(paste(sapply(c(nmodflr$loglikelihood, dmodflr$loglikelihood, lrdflr), digits=2, round)))
        pvals[i] = pchisq(lrdflr, df=1, lower.tail=FALSE)
        print(proc.time() - ptm)
    }
    # Make into dataframe:
    df = data.frame(node=1:NN, parent=declist$parent[1:NN])
    df$pout = pvals
    df$padj = p.adjust(df$pout)
    df$padj[df$padj == 0] = min(df$padj[df$padj != 0])
    df$rawlog10p = -log10(df$padj)
    # Cap and threshold the adjusted log10p:
    df$log10p = df$rawlog10p
    df$log10p[df$log10p < 3] = 0
    df$log10p[df$log10p > cutp] = cutp
    title = ''
    return(list(log10p=df$log10p, isnp=isnp, rawlp=df$rawlog10p, df=df, title=title))
}




