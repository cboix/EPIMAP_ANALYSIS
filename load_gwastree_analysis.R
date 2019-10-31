#!/usr/bin/R
# --------------------------------------------------
# Loads primary variables for the GWAS tree analyses
# --------------------------------------------------
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
library(dendextend)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(GenomicRanges)
library(dplyr)
library(cba)
# Fast libraries for glm:
options(scipen=45) # So we dont get issues writing integers into bedfiles
options(repos='http://cran.rstudio.com/')

# Arguments:
usetree = 'enhancers'
# usetree = 'gwastree'
# usetree = 'correlation'
# usetree = 'roadmap'
tol = 2500  # Plus/minus distance - window for enhancer overlaps
singlematch = FALSE # Use 1-to-1 SNP enhancer mapping only?
plotting.only = TRUE # Load data for plotting only?
plot.trees = FALSE 
if (length(gtargs)==0) {
    print("Using default arguments. Only loading what is needed for plotting")
} else {        
    usetree = gtargs[1]
    tol = as.integer(gtargs[2])
    singlematch = as.logical(gtargs[3])
    if (length(gtargs) > 3){ 
        plotting.only = as.logical(gtargs[4])
    }
    if (length(gtargs) > 4){ 
        plot.trees = as.logical(gtargs[5])
    }
}

if (singlematch){
    midpref = paste0('_', usetree, '_e', tol, '_single_')
} else {
    midpref = paste0('_', usetree, '_e', tol, '_all_')
}

# Load specific distance matrices:
if (usetree == 'correlation'){
    print("[STATUS] Loading distance matrices")
    source(paste0(bindir, 'load_distance_matrices.R'))
    setprefix = 'distances_'
} else {
    setprefix = paste0(usetree, '_')
}

# Make directories:
today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "gwas_tree_analysis/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, setprefix)
treeimgpref = paste0(imgdir, usetree, '_e', tol, '_')

if (usetree == 'roadmap'){
    keepbss = meta$id[meta$Project %in% c('ENCODE 2012','Roadmap 2015')]
}

# ----------------------------------------
# Load in GWAS and overlap with enhancers:
# ----------------------------------------
print("[STATUS] Loading GWAS catalog")
gwcatfile = 'gwascatalog_may03_2019_noquotes.txt'
gwrdafile = sub(".txt", ".Rda", gwcatfile)
gwrsidfile = sub(".txt", "_rsid.txt", gwcatfile)
gwintrdafile = sub(".txt", "_intersections.Rda", gwcatfile)
if (!file.exists(gwrdafile)){
    # Load GWAS Catalog, reduce to loc, pmid, pval:
    gwdf = read.delim(gwcatfile, header=T, stringsAsFactors=F, sep="\t")
    gwdf$uid = paste0(gwdf$pubMedID, ' - ', gwdf$trait)
    write.table(gwdf[, c('chrom','chromStart','uid','name')], gwrsidfile, sep="\t", col.names=T, row.names=F, quote=F)
    gwssdf = gwdf[, c('pubMedID','trait','initSample', 'replSample','uid', 'pubDate')]
    gwssdf = unique(gwssdf)
    gwssdf$old.sampsize = sapply(gwssdf$initSample, function(x){sum(munge.nos(x))})
    gwssdf$old.rep.size = sapply(gwssdf$replSample, function(x){sum(munge.nos(x))})
    gwssdf$sampsize = sapply(gwssdf$initSample, prune.cases)
    gwssdf$rep.size = sapply(gwssdf$replSample, prune.cases)
    # sum(gwssdf$sampsize > 50000)
    # sum(gwssdf$sampsize > 100000)
    gwdf = gwdf[,c('chrom', 'chromStart', 'chromEnd', 'pubMedID','trait', 'pValue', 'uid')]
    # Collapse multi-counted SNPs:
    gwdf = aggregate(pValue ~ chrom + chromStart + chromEnd + uid + trait + pubMedID, gwdf, min)
    # Filter out chrY: 
    gwdf = gwdf[gwdf$chrom != 'Y',]
    # Prune the HLA region
    # chr6: 29691116-33054976
    gwdf = gwdf[!(gwdf$chrom == '6' & gwdf$chromEnd > 29691116 & gwdf$chromEnd < 33054976),]
    # sum(gwdf$chrom == '6' & gwdf$chromStart > 29691116 & gwdf$chromEnd < 33054976)
    # ------------------------------------------------------------------------
    # Prune SNPs according to Roadmap pruning procedure:
    # The pruning procedure considered each SNP in ranked order of P value 
    # with the most significant coming first, and we retained a SNP 
    # if there was no already retained SNP on the same chromosome within 1 Mb.
    # ------------------------------------------------------------------------
    # Order catalog by SNP significance:
    print("[STATUS] Pruning GWAS catalog - keeping most significant w/in 1 Mb")
    gwdf = gwdf[order(gwdf$pValue),]
    # Prune, on a per-uid basis:
    prune.snps = function(suid, df=gwdf, dist=1e6, quiet=TRUE){
        subdf = df[df$uid == suid,]
        keptdf = c()
        chrlist = c(as.character(1:22), 'X')
        keptlist = sapply(rep(0, 23), function(x){c()})
        names(keptlist) = chrlist
        for (i in 1:nrow(subdf)){
            chrom = as.character(subdf$chrom[i])
            loc = subdf$chromStart[i]
            # Check if any within range/add:
            if (is.null(keptlist[[chrom]])){
                keptlist[[chrom]] = loc
            } else {
                nclose = sum(abs(loc - keptlist[[chrom]]) < dist)
                if (nclose ==  0){
                    keptlist[[chrom]] = c(keptlist[[chrom]], loc)
                }
            }
        }
        # Report out - chrom, loc, uid is sufficient:
        kdf = c()
        for (chrom in chrlist){
            if (!is.null(keptlist[[chrom]])){
                kdf = rbind(kdf, data.frame(chrom=chrom, 
                                            chromStart=keptlist[[chrom]],
                                            uid=suid))
            }
        }
        if (!quiet){ print(paste(suid, ': Kept', nrow(kdf), 'of',nrow(subdf))) }
        return(kdf)
    }
    # Prune:
    ut = unique(gwdf$uid)
    kept.snps = ldply(ut, df=gwdf, dist=5e3, prune.snps)
    print(dim(gwdf))
    gwdf = merge(kept.snps, gwdf)
    print(dim(gwdf))
    # Ranges object:
    gwgr = GRanges(paste0('chr', gwdf$chrom), IRanges(gwdf$chromStart, gwdf$chromEnd))
    # Save dataframe and the gr object:
    save(gwdf, gwgr, gwssdf, file=gwrdafile)
} else {
    load(gwrdafile)
}

if (!plotting.only){
    # ---------------------------------
    # Load in the enhancer coordinates:
    # ---------------------------------
    print("[STATUS] Loading DHS list")
    ddir = 'DHS_Index_WM201902/'
    dpref = 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19'
    dmlfile = paste0(ddir, dpref, '.core.srt.txt')
    dmlnamfile = paste0(ddir, dpref, '_r200_e0_names.core.srt.tsv')
    dmlrdafile = paste0(ddir, dpref, '_enh.core.srt.Rda')
    # Load indices (needed for mapping, etc.)
    enhind = as.numeric(scan('Enhancer_indices.txt', 'c')) + 1
    if (!file.exists(dmlrdafile)){
        dmldf = read.table(dmlfile, header=F, stringsAsFactors=F, sep="\t")
        names(dmldf) = c('chr','start','end','name')
        # Reorder dml:
        dmlnam = read.table(dmlnamfile, header=T, stringsAsFactors=F, sep="\t")
        dmldf = merge(dmldf, dmlnam)
        dmldf = dmldf[order(dmldf$cls), ]
        # Enhancer indices are 0-indexed - turn to 1-index:
        enhdf = data.frame(cls = enhind)
        enhdf = merge(enhdf, dmldf)
        rm(dmldf, dmlnam)
        save(enhdf, file=dmlrdafile)
    } else {
        load(dmlrdafile)
    }

    # -------------------------------------------------
    # Overlap SNPs with enhancers, with some tolerance:
    # -------------------------------------------------
    dmgr = GRanges(enhdf$chr, IRanges(enhdf$start - tol, enhdf$end + tol))
    if (!file.exists(gwintrdafile)){
        print('[STATUS] Overlapping SNPs with enhancers')
        qdf = suppressWarnings(data.frame(findOverlaps(gwgr, dmgr)))
        qdf$uid = gwdf$uid[qdf$queryHits]
        # Non-unique and unique counts:
        qnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(unique(x))})
        qallnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(x)})
        # Almost all SNPs are within an enhancer +/- tolerance:
        print(round(length(unique(qdf$queryHits)) / length(gwgr) * 100, 2))
        # Alternatively: 1-1 assign each SNP to nearest enhancer:
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
        # Keep ~ 93k of 113k pruned + no HLA SNPs
        # Keep ~ 72k of 88k very aggresively pruned SNPs
        print(round(length(unique(ndf$queryHits)) / length(gwgr) * 100, 2))
        save(ndf, qdf, nnumdf, nallnumdf, qnumdf, qallnumdf, file=gwintrdafile)
    } else {
        load(gwintrdafile)
    }
}

# --------------------------------------
# Calculate the span of the DHS dataset:
# --------------------------------------
calc.span = FALSE
if (calc.span){
    # For all:
    dmlgr = GRanges(dmldf$chr, IRanges(dmldf$start, dmldf$end))
    dmlgr2 = reduce(dmlgr)
    gen.width = 3036303846
    sum(width(dmlgr2)) / gen.width
    # For reduced:
    enhgr = GRanges(enhdf$chr, IRanges(enhdf$start, enhdf$end))
    enhgr2 = reduce(enhgr)
    sum(width(enhgr2)) / gen.width
}

# --------------------------------------------
# Load enhancer jaccard - see different trees:
# --------------------------------------------
if (usetree == 'enhancers' || usetree == 'roadmap'){
    print('[STATUS] Loading enhancers jaccard matrix:')
    emat = read.delim(gzfile('Enhancer_jaccard.tsv.gz'), sep="\t", header=F)
    matnames = scan('Enhancer_matrix_names.txt', "c")
    rownames(emat) = matnames
    colnames(emat) = matnames
    if (usetree == 'roadmap'){
        emat = emat[keepbss, keepbss]
    }
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
    NL = length(lab)

    if (plot.trees){
        pdf(paste0(imgpref, sub("\\.","_",method),"_link_jacc.pdf"), width=14.5, height=12, onefile=T)
        plot.new()
        circle_size = unit(1, "snpc") # snpc unit gives you a square region
        pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                              just = c("left", "center")))
        par(omi = gridOMI(), new = TRUE)
        circleplot(dend=dend3, lab=lab)
        upViewport()
        circos.clear()
        draw(pd.legend, x=circle_size, just="left")
        # title(paste0(method, '-linkage of Jaccard Similarity of Enhancers'))
        dev.off()

        # Plot the matrix as well:
        emat.reord = as.matrix(emat[lab,lab])
        labels = meta[lab, 'GROUP']
        faclabels = as.matrix(as.numeric(labels))
        lablist = label.runs(faclabels, labels, rdcol)
        # Update palette and legend for jaccard:
        jacc.col_fun = function(x, pal=rev(viridis(100))){
            palette = rev(pal)
            bin <- cut(x, seq(0, 100, length.out=length(palette)), include.lowest=T) 
            palette[bin] }
        jacc.legend = Legend(at = seq(0, 100, 25), 
                             labels = c('0%','25%','50%','75%','100%'),
                             labels_gp = gpar(fontsize=5),
                             title_gp = gpar(fontsize=5, fontface='bold'),
                             col_fun=jacc.col_fun, title_position = "topleft", 
                             title='Jaccard Similarity', direction = 'vertical')
        plegend = packLegend(jacc.legend)
        # Figure:
        pdf(paste0(imgpref, sub("\\.","_",method),"_link_jacc_matrix.pdf"), width=9.5, height=8, onefile=T)
        # image(emat.reord, col=rev(col1), zlim=c(0,1), axes=F)
        sp=0.1
        layout(matrix(c(1:2),1,2), widths=c(1.5,8), heights=c(8), TRUE)
        par(yaxs="i")
        par(xaxs="i")
        # Metadata matrix:
        par(mar=c(sp,5,sp,sp))
        meta.image(metamat[lab,5:1], colvals=colvals, cex=0, horiz=T, useRaster=TRUE)
        box.pad = 0.011
        xx = space.1d(lablist[[1]], box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
        xx = space.1d(xx, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
        rx = c(0.01, 0.06, 0.2, 0.25, 0.27)
        x = par()$usr[1]-rx*(diff(par()$usr[1:2]))
        text(y=xx, x=x[5], labels=lablist[[2]],
             srt=0, adj=1, xpd=TRUE, cex=.5, col=lablist[[3]])
        par(xpd=TRUE)
        segments(x0=x[1], y0=lablist[[1]], x1=x[2], y1=lablist[[1]], col=lablist[[3]])
        segments(x0=x[3], y0=xx, x1=x[2], y1=lablist[[1]], col=lablist[[3]])
        segments(x0=x[3], y0=xx, x1=x[4], y1=xx, col=lablist[[3]])
        par(xpd=FALSE)
        draw(plegend, x = unit(0.4,'in'), y=unit(7.5,'in'), just = "top")
        abline(h=par()$usr[3:4],lty=1,lw=0.5)
        abline(v=par()$usr[1:2],lty=1,lw=0.5)
        par(mar=rep(sp,4))
        image(emat.reord, col=rev(viridis(100)), zlim=c(0,1), axes=F, useRaster=TRUE)
        dev.off()
    }
}

# ----------------------------------
# Load gwas matrix for tree creation
# ----------------------------------
if (usetree == 'gwastree'){
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
    # Rename:
    clsn = paste0('c',1:length(epinames) - 1)
    names(epinames) = clsn
    epimat = gwmat
    rownames(epimat) = epinames[rownames(gwmat)]
    # Tree from GWAS enrichments:
    metric='euclidean'
    dt = dist(epimat, method=metric)
    if (metric == 'jaccard') { dt = dist(epimat > 0, method=metric) }
    # method = 'complete'
    method = 'ward.D'
    ht <- hclust(dt, method=method)
    ht$order <- order.optimal(dt, ht$merge)$order
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
}

# ----------------------------------------
# Load matrix and get epigenomes per node:
# ----------------------------------------
if (!plotting.only){
    print("[STATUS] Loading enhancer matrix")
    fullmatfile = 'Enhancer_H3K27ac_matrix_062619.mtx.gz'
    enhmatfile = 'Enhancer_H3K27ac_matrix_enhonly_062619.mtx.gz'
    enhmargfile = 'Enhancer_H3K27ac_margins_enhonly_062619.tsv.gz'
    if (!file.exists(enhmatfile)){
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
}


# Calculate the average genomic coverage from enhmat:
if (calc.span){
    spans = c()
    uqc = sort(unique(mat$col))
    for (i in uqc){
        print(i)
        urow = mat$row[mat$col == i]
        u.dmlgr = reduce(dmlgr[urow,])
        spans = c(spans,sum(width(u.dmlgr)) / gen.width)
    }
    ms = mean(spans) * 100
    print(ms)
    png(paste0(img, 'clusters/hist_genome_coverage_epi.png'), res=450, units='in', width=6, height=4)
    par(mar=c(4,4,1,1), yaxs='i', xaxs='i')
    hist(spans * 100, 40, col='darkgrey', border='white', ylim=c(0,150),yaxt='n',
         ylab='Number of Epigenomes', xlab='% of Genome Covered', main='')
    axis(2,las=1)
    abline(v=ms, col='red', lty='dashed')
    text(x=ms, y=125, label=paste0('Mean = ', round(ms, 2), '%'), col='red', adj=-.1) 
    box()
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

# ------------------------------------
# Recursively fill out the tree nodes:
# ------------------------------------
if (!plotting.only){
    cdllfile = paste0('consensus_object_', usetree, '_062819.Rdata')
    if (!file.exists(cdllfile)){
        # Recursive function to get consensus and merge:
        get_consensus <- function(subdend, cdll){
            node = attributes(unclass(subdend))$nodePar$pch
            # Get id of node:
            print(node)
            nset = cdll$cons[[node]]
            print(head(nset))
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
        # Set top node difference to top node consensus:
        if (is.na(cdll$diff[[1]])){
            print("toplevel is na")
            cdll$diff[[1]] = cdll$cons[[1]]
        }
        if (is.na(cdll$novel[[1]])){
            print("toplevel is na")
            cdll$novel[[1]] = cdll$union[[1]]
        }
    } else {
        print("[STATUS] Loading cdll from file")
        load(cdllfile)
    }
}

# Need for outer rim of plotting:
cdlenfile = paste0('consensus_object_lengths_', usetree, '_062819.Rdata')
if (!file.exists(cdlenfile)){
    cdlenlist = lapply(cdll, function(x){sapply(x,length)})
    # Look at lengths - if no 1 (NA), ok.
    lapply(cdlenlist, function(x){head(x, 10)})
    save(cdlenlist, file=cdlenfile)
} else {
    load(cdlenfile)
}

nbpfile = paste0('consensus_object_nbp_', usetree, '_062819.Rdata')
if (!file.exists(nbpfile)){
    enhdf$size = enhdf$end - enhdf$start
    nbplist = list(cons=sapply(rep(0, NN), list),
                   diff=sapply(rep(0, NN), list))
    for (i in 1:NN){
        print(i)
        for (type in c('cons','diff')){
            cl = cdll[[type]][[i]]
            if (length(cl) > 0){
                sl = enhdf$size[enhdf$cls %in% cl]
                nbplist[[type]][[i]] = sum(sl)
            }
        }
    }
    save(nbplist, file=nbpfile)
} else {
    load(nbpfile)
}


if (!plotting.only){
    # ----------------------------
    # From intersection, get SNPs:
    # ----------------------------
    # Enhancer index map to matrix:
    # (did intersection on enhancers only)
    enhmap = rep(0, max(enhind))
    enhmap[enhind] = 1:length(enhind)
    print(max(qdf$subjectHits))
    print(max(ndf$subjectHits))

    # NOTE: Don't load this if doing only logistic regression:
    treefile = paste0('treedf', midpref, 'diff_snpint.tsv')
    dflist = list()
    types = c('diff','cons','union','novel')
    print("[STATUS] Loading dflist - intersection of sets and overlaps")
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
            print(paste("[STATUS] Writing for:",type))
            tfile = paste0('treedf', midpref, type, '_snpint.tsv')
            write.table(dflist[[type]], tfile, , quote=F, sep="\t", row.names=F)
        }
    } else { 
        for (type in types){ 
            print(paste("[STATUS] Reading in:",type))
            tfile = paste0('treedf', midpref, type, '_snpint.tsv')
            dflist[[type]] = read.delim(treefile, sep="\t", header=T)
        }
    }

    # Remove mat and list for space:
    if (nrow(dflist[['diff']])  > 0){
        gc()
        rm(mat)
        gc()
    }
}

# --------------------------------------
# Make tree and descendants information:
# --------------------------------------
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

# Get the component cell/tissue groups for the large nodes
leafmeta = data.frame(label=labels(dend3), id=lab, GROUP=meta[lab,'GROUP'])
leafmeta$label = as.character(leafmeta$label)
nodemeta = ldply(1:NN, function(i){ 
                     x = declist$dec[[i]]
                     df = data.frame(node=i, GROUP=leafmeta[leafmeta$label %in% x,'GROUP'], count=1) 
                     df = aggregate(count ~ node + GROUP, df, length)
                     df$total = sum(df$count)
                     df$frac = df$count/df$total
                     # NOTE: Deal with 50-50 cases (neither max) with avg:
                     df$rank = rank(-df$frac, ties.method='average')
                     # Maximal if at least some fraction and max:
                     fcut = 0.5
                     j = which((df$rank == 1) * (df$frac > fcut) == 1)
                     if (length(j) == 0){ df$maxgroup = 'Multiple' 
                     } else { df$maxgroup = df$GROUP[j] }
                     return(df)
                    })

nodetissue = unique(nodemeta[, c('node','total','maxgroup')])
names(nodetissue)[3] = 'GROUP'
nodetissue =merge(nodetissue, rdcol, all.x=TRUE)
nodetissue$COLOR[is.na(nodetissue$COLOR)] = 'grey80'  # Add color for Multiple
nodetissue$category[is.na(nodetissue$category)] = 'Other'  # Add color for Multiple
nodetissue.stats = aggregate(node ~ GROUP, nodetissue, length) # Only 191 assigned to 1

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

# Update palette and legend for pvalues:
palette = col3
CUTP = 12
col_fun = colorRamp2(c(0, CUTP), c(palette[1], palette[length(palette)]))
pval.legend = Legend(at = seq(3, CUTP, 3), col_fun=col_fun, title_position = "topleft", title = "p-value")
pd.legend.ext = packLegend(col.list$group, col.list.h$project, col.list$sex, col.list$type, col.list$lifestage, pval.legend)

# For labeling nodes:
nodedf = data.frame(get_nodes_xy(dend3, type = 'rectangle'))
nodedf$node = 1:NN

# UID list for chunks:
uids = sort(as.character(unique(gwdf$uid)))
if (!file.exists('full_uidlist.txt')){
    write.table(uids, 'full_uidlist.txt', quote=F, row.names=F, col.names=F, sep="\t")
}
# UIDs with NSNPs:
uids = sort(as.character(unique(gwdf$uid)))
countgw = aggregate(pValue ~ uid, gwdf, length)
countgw$tail = ''
countgw$tail[countgw$pValue > 1] = 's'
snpuids = paste0(countgw$uid, ' (', countgw$pValue, ' snp', countgw$tail, ')')
snpuids = sort(snpuids)
if (!file.exists('full_snpuidlist.txt')){
    write.table(snpuids, 'full_snpuidlist.txt', quote=F, row.names=F, col.names=F, sep="\t")
}

# List of 24 traits for testing:
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


# Expanded list of 80 traits:
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

if (!file.exists('original_traitlist.txt')){
    write.table(tlist, 'original_traitlist.txt', quote=F, row.names=F, col.names=F, sep="\t")
}

if (!file.exists('expanded_traitlist.txt')){
    write.table(expandedlist, 'expanded_traitlist.txt', quote=F, row.names=F, col.names=F, sep="\t")
}

# ------------------------------
# Get the leaf pvals by HG test:
# ------------------------------
leaffile = paste0('logreg', midpref, 'leaf_hgpvals_elist.Rda')
if (!file.exists(leaffile)){
    print("[STATUS] Calculating pvals for leaves by HG test")
    leafmat = matrix(0, nrow=length(expandedlist), ncol=NL)
    rownames(leafmat) = expandedlist
    for (strait in expandedlist){
        leafdf = get_pvalues_leaves(strait, dflist=dflist, 
                                    cdlenlist=cdlenlist, declist=declist)
        leafmat[strait,] = leafdf$log10p
    }
    save(leafmat, file=leaffile)
} else {
    load(leaffile)
}

# ------------------------------
# Get the leaf pvals by HG test:
# ------------------------------
all.leaffile = paste0('logreg', midpref, 'leaf_hgpvals_alltraits.Rda')
if (!file.exists(all.leaffile)){
    ut = as.character(unique(gwdf$uid))
    print("[STATUS] Calculating pvals for leaves by HG test")
    all.leafmat = matrix(0, nrow=length(ut), ncol=NL)
    rownames(all.leafmat) = ut
    for (i in 1:length(ut)){
        print(i)
        suid = ut[i]
        leafdf = get_pvalues_leaves(suid, dflist=dflist, 
                                    cdlenlist=cdlenlist, declist=declist)
        all.leafmat[suid,] = leafdf$log10p
    }
    save(all.leafmat, file=all.leaffile)
} else {
    load(all.leaffile)
}



# Plot example/metadata trees:
if (plot.trees){
    # ----------------
    # Plot basic tree:
    # ----------------
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
}

