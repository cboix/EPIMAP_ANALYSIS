#!/usr/bin/R
# ------------------------------------------------
# Test a multiview local NMF approximation for genes
# - Get local region - H3K27ac + corrs + data, look at 
# - May not work with many marks
# - Compare to std corr...
# ------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))

library(rhdf5)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(Matrix)
library(scales)
library(plotrix)

#' Calculate position relative to par()$usr 
#'
#' @param axis 1 is x; 2 is y;
#' @param shift percent of axis to shift over
#' @return position shifted from start of x or y axis
#' @export
parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}

# For plotting multiple arrows on a segment:
multarrows <- function(x0, y0, x1, y1, n_arr, ...) {
    x <- seq(x0, x1, length=n_arr + 1)
    y <- seq(y0, y1, length=n_arr + 1)
    arrows(x[-length(x)], y[-length(y)], x[-1], y[-1], ...) }

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "linking/")
genedir = paste0(img, "linking/genes/")
gwasdir = paste0(img, "linking/gwas/")
cmd = paste('mkdir -p', imgdir, genedir, gwasdir)
system(cmd)
imgpref = paste0(imgdir, 'corrtest_')

# ---------------------
# Load expression data:
# ---------------------
lddir = 'linking_data/'
tssfile = 'Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz'
matfile = paste0(lddir, 'merged_log2fpkm.pc.mtx.gz')
generdafile = paste0(lddir, 'gene_mat_tss.Rda')
gmapfile = 'Annotation/gencode.gene.name.mapping.tsv'

if (!file.exists(generdafile)){
    tssdf = read.delim(gzfile(tssfile), header=F, sep="\t", stringsAsFactors=F)
    names(tssdf) = c('chr','tss','tss2','gene')
    tssdf = tssdf[!(tssdf$chr %in%  c('chrY', 'chrM')),]

    # Matrix to wide:
    matlong = read.delim(gzfile(matfile), header=T)
    matwide = spread(matlong, id, log2fpkm)
    gmat = as.matrix(matwide[,-1])
    rownames(gmat) = matwide[,1]
    rm(matwide, matlong)

    # Subset to shared genes:
    kept.genes = sort(tssdf$gene[tssdf$gene %in% rownames(gmat)])
    tssdf = tssdf[tssdf$gene %in% kept.genes,]
    tssdf$gene = factor(tssdf$gene, levels=kept.genes)
    tssdf = tssdf[order(tssdf$gene), ]
    gmat = gmat[kept.genes,]

    gc()
    save(gmat, tssdf, file=generdafile)
} else {
    load(generdafile)
}

genemap = read.delim(gmapfile, header=F)
names(genemap) = c('gene','symbol')
anno = merge(tssdf, genemap) # Note: not same as tssdf file.


# --------------------------------------
# Read in the dhs locations and indices:
# --------------------------------------
lddir = 'linking_data/'
mapfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv'
locfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt'
indfile = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_masterlist_indices.tsv'
marknamesfile = paste0(lddir, 'mark_matrix_names.txt')
dhsrdafile = paste0(lddir, 'dhs_loc.Rda')

# Load data:
if (!file.exists(dhsrdafile)){
    # Locations:
    locmap = read.delim(mapfile, header=T, stringsAsFactors=F)
    locdf = read.delim(locfile, header=F, stringsAsFactors=F)
    colnames(locdf) = c('chr','start','end','name')
    head(locmap)
    head(locdf)
    # Merge mapping:
    locdf = merge(locmap, locdf)
    locdf = locdf[order(locdf$cls),]
    locdf$mid = (locdf$start + locdf$end) / 2
    rm(locmap)
    gc()
    # Matrix names + indices:
    mnames = scan(marknamesfile,'c')
    matind = as.numeric(scan(indfile, 'c')) + 1
    locdf = locdf[matind,]
    # Save for loading:
    save(mnames, locdf, matind, file=dhsrdafile)
} else {
    load(dhsrdafile)
}
gc()


# -------------------------------------------
# Read and process the GTF for gene plotting:
# -------------------------------------------
parse.gtf.info = function(x){
    x = strsplit(x, "; ")[[1]]
    x = unlist(strsplit(x, " "))
    out = x[seq(2,length(x), by=2)]
    names(out) = x[seq(1,length(x), by=2)] 
    return(out) }

gtf.file = 'Annotation/gencode.v30lift37.basic.annotation.gtf.gz'
gtf.rdafile= 'Annotation/gencode.v30lift37.basic.annotation.Rda'
if (!file.exists(gtf.rdafile)){
    tab = read.delim(gtf.file, sep="\t", skip=5, header=F, stringsAsFactors=F)
    # Find longest tx per gene:
    txdf = tab[tab[,3] == 'transcript',c(1,4,5,7,9)]
    names(txdf) = c('chr','start','end','strand','info')
    tlist = sapply(as.character(txdf$info), parse.gtf.info)
    txdf$txid = sapply(tlist, function(x){x[['transcript_id']]})
    txdf$symbol = sapply(tlist, function(x){x[['gene_name']]})
    txdf$txtype = sapply(tlist, function(x){x[['transcript_type']]})
    txdf$length = txdf$end - txdf$start
    txdf = txdf[txdf$txtype == 'protein_coding',]
    # Reduce to longest transcript only:
    txdf = merge(txdf, aggregate(length ~ symbol, txdf, max))
    txdf$tss = with(txdf,start * (strand == '+') + end * (strand == '-'))
    # Match exons to these transcripts:
    exdf = tab[tab[,3] == 'exon',c(1,4,5,7,9)]
    names(exdf) = c('chr','start','end','strand','info')
    exdf = exdf[grep('protein_coding', exdf$info),]
    elist = sapply(as.character(exdf$info), parse.gtf.info)
    exdf$txid = sapply(elist, function(x){x[['transcript_id']]})
    exdf$exnum = sapply(elist, function(x){x[['exon_number']]})
    exdf$exnum = as.numeric(exdf$exnum)
    exdf = exdf[exdf$txid %in% txdf$txid,]
    exdf = merge(exdf, txdf[,c('txid','symbol')])
    # Match UTRs to these transcripts:
    uxdf = tab[tab[,3] == 'UTR',c(1,4,5,7,9)]
    names(uxdf) = c('chr','start','end','strand','info')
    uxdf = uxdf[grep('protein_coding', uxdf$info),]
    ulist = sapply(as.character(uxdf$info), parse.gtf.info)
    uxdf$txid = sapply(ulist, function(x){x[['transcript_id']]})
    uxdf = uxdf[uxdf$txid %in% txdf$txid,]
    uxdf = merge(uxdf, txdf[,c('txid','symbol')])
    save(exdf, txdf, uxdf, file=gtf.rdafile)
} else {
    load(gtf.rdafile)
}


# -------------------------------------
# Read in the precomputed correlations:
# -------------------------------------
corrfile = paste0(lddir, 'mark_precomputed_corr.Rda')  # For fast reload
if (!file.exists(corrfile)){
    h1 = read.delim(gzfile(paste0(lddir, 'H3K27ac_precomputed_corr.tsv.gz')), header=F)
    h1 = as.numeric(h1$V1)
    h2 = read.delim(gzfile(paste0(lddir, 'H3K4me1_precomputed_corr.tsv.gz')), header=F)
    h2 = as.numeric(h2$V1)
    save(h1, h2, file=corrfile)
} else {
    load(corrfile)
}

# Data files:
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'


# -------------------------------------------------------
# Choose the gwas of interest and find all relevant loci:
# -------------------------------------------------------
TRAIT = 'Coronary artery disease'
# TRAIT = 'QT interval'
# TRAIT = 'Atrial fibrillation'
# TRAIT = "Crohn's disease"
# TRAIT = "Alzheimer's disease"
# TRAIT = "Waist-to-hip ratio"
traitstr = gsub("'","_", gsub(" ", "_", tolower(TRAIT)))

# Load in the gwas file, pruned - load CAD:
gwrdafile = '~/EPIMAP_ANALYSIS/db/gwascatalog_may03_2019_noquotes.Rda'
load(gwrdafile)
NIND = 10000 # Remove all the low GWAS
NSNP = 10
nsnpdf = aggregate(pValue ~ uid, gwdf[gwdf$pValue < 1e-8,], length)
# NOTE: Some duplicates exist.
gwssdf2 = aggregate(sampsize ~ uid, gwssdf, max)
gwssdf = merge(gwssdf, gwssdf2)
gwssdf = unique(merge(gwssdf, nsnpdf))
keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]

# Get the kept gwas with most pvals for trait in question:
sub.gwssdf = gwssdf[grep(TRAIT, gwssdf$trait),]
sub.gwssdf = sub.gwssdf[order(sub.gwssdf$pValue, decreasing=T),]
guid = sub.gwssdf$pubMedID[1]
sub.gwdf = unique(gwdf[gwdf$pubMedID == guid,c('chrom','chromStart','pValue')])
sub.gwdf = sub.gwdf[order(sub.gwdf$chrom),]
sub.gwdf = sub.gwdf[sub.gwdf$pValue < 1e-8,]
sub.gwdf$chr = paste0('chr', sub.gwdf$chrom)

# Get the nearest genes to the GWAS loci:
locgr = with(sub.gwdf, GRanges(paste0('chr', chrom), IRanges(chromStart, chromStart)))
tssgr = GRanges(seqnames=tssdf$chr, IRanges(tssdf$tss, tssdf$tss), name=tssdf$gene)
sub.gwdf$gene = as.character(tssgr$name[nearest(locgr, tssgr)])
sub.gwdf = sub.gwdf[order(sub.gwdf$chr),]
center.genes = unique(sub.gwdf$gene)
center.symbols = sapply(center.genes, function(ensg){as.character(anno$symbol[anno$gene == ensg])})

# Make UID specific directory for images:
traitdir = paste0(gwasdir, traitstr, '_', guid, '/')
cmd = paste('mkdir -p', traitdir)
system(cmd)

# -----------------------------------
# For each gene, make a linking plot:
# -----------------------------------
# TODO: later test if JCAD looks exactly same with this plotting routine.
# gnames = as.character(sapply(center.genes, function(x){genemap$symbol[genemap$gene == x]}))
window = 1e6
cpairsdf = NULL
for (k in 1:length(center.genes)){
    ensg = center.genes[k]
    gene = as.character(anno$symbol[anno$gene == ensg])
    print(paste(k, gene))
    infoline = anno[anno$symbol == gene,]
    chrom = infoline$chr

    # Get overlapping genes:
    genetss <- with(infoline, GRanges(seqnames=chr, ranges=IRanges(start=tss - window/2, end=tss + window/2)))
    tssgr = GRanges(seqnames=tssdf$chr, IRanges(tssdf$tss, tssdf$tss), name=tssdf$gene)
    anno.ovl = data.frame(findOverlaps(genetss, tssgr))
    kept.gene = as.character(tssgr$name[anno.ovl$subjectHits])
    gloc = tssdf$tss[anno.ovl$subjectHits]
    gnames = as.character(sapply(kept.gene, function(x){genemap$symbol[genemap$gene == x]}))
    ord = order(gloc) 
    gnames = gnames[ord]
    kept.gene = kept.gene[ord]
    gloc = gloc[ord]
    names(gloc) = kept.gene

    # Load tested pairs:
    if (is.null(cpairsdf) || cpairsdf$chrom[1] != chrom){
        print(paste("Loading links for", chrom))
        cfile = paste0(lddir, chrom, '_pairs_ENH_ovl_df.tsv.gz')
        # genes = sub.gwdf$gene[sub.gwdf$chr == chrom]
        cpairsdf = read.delim(cfile, header=T, sep="\t", stringsAsFactors=F)
        cpairsdf$mind = cpairsdf$mind + 1  # Was 0-indexed
        cpairsdf$chrom = chrom
    }
    sub.ovldf = cpairsdf[cpairsdf$gene %in% kept.gene,]
    sub.ovldf$chr = chrom
    sub.ovldf$corr = h1[sub.ovldf$ind]
    sub.ovldf$corr_k4 = h2[sub.ovldf$ind]

    # Get all links:
    sub.ovldf$name = locdf$name[sub.ovldf$mind]
    sub.ovldf$mid = locdf$mid[sub.ovldf$mind]
    sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]
    sub.ovldf$dist = sub.ovldf$mid - infoline$tss 
    sub.ovldf = sub.ovldf[abs(sub.ovldf$dist) < (window / 2),]
    sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]
    # Kept locations:
    kept.locdf = unique(sub.ovldf[,c('mid','mind', 'name')])
    kept.dhs = kept.locdf$mind
    kept.mid = kept.locdf$mid
    rownames(kept.locdf) = kept.locdf$name
    kept.locdf$id = 1:nrow(kept.locdf)

    # Slice locations from hdf5:
    h5f = H5Fopen(H3K27ac.file)
    h5d = h5f&"matrix"
    sub.mmat = h5d[,matind[kept.dhs]]
    H5Dclose(h5d)
    H5Fclose(h5f)
    sub.mmat = t(sub.mmat)
    colnames(sub.mmat) = mnames

    # Make average dhs tracks:
    tform = make.tform(meta[mnames,'GROUP'], norm=TRUE, u=odf$GROUP)
    sub.avg.mmat = sub.mmat %*% tform 

    # Average expression:
    tform = make.tform(meta[colnames(gmat),'GROUP'], norm=TRUE, u=odf$GROUP)
    sub.gmat = gmat[kept.gene,]
    sub.avg.expr = sub.gmat %*% tform
    sub.avg.expr[is.na(sub.avg.expr)] = 0


    # Quick multi-view factorization here:
    library(IntNMF)
    datalist = list(expr=t(sub.gmat), mark=log2(t(sub.mmat[, colnames(sub.gmat)]) + 1))

    k = 2
    nmf.mnnals(dat = datalist, k = k, maxiter = 200, st.count = 20, n.ini = 30, ini.nndsvd = TRUE,
               seed = TRUE)

    # Alternatively with liger:
    library(liger)
    rownames(sub.mmat) = as.character(1:nrow(sub.mmat))
    datalist = list(expr=t(sub.gmat), mark=t(log2(sub.mmat[, colnames(sub.gmat)] + 1)))
    # datalist = list(expr=sub.gmat, mark=log2(sub.mmat[, colnames(sub.gmat)] + 1))
    dl <- createLiger(datalist, make.sparse=T, remove.missing=F)

    # sub.mmat = 
    # Norm, select genes with RNA, scale data:
    dl <- liger::normalize(dl)
    dl@scale.data = lapply(dl@norm.data, function(x){as.matrix(t(x))})
    # dl <- selectGenes(dl , datasets.use = 2)
    # dl <- scaleNotCenter(dl)
    # dl@norm.data[[i]][dl@var.genes,]
   # image(dl$H[[1]])
   # dev.new()
   # image(dl$H[[2]])
    # ------------------------------------
    # Perform the multiview factorization:
    # + plot results
    # ------------------------------------
    NFACT = 5
    # datalist
    # datalist = list(expr=sub.gmat, mark=log2(sub.mmat[, colnames(sub.gmat)] + 1))
    # dl <- optimizeALS(datalist, k=NFACT)
    # dl@scale.data = datalist
    dl <- optimizeALS(dl, k=NFACT)

    # Plot the t-SNE without normalization:
    dl <- liger::runTSNE(dl, use.raw=T)
    p1 <- plotByDatasetAndCluster(dl, return.plots=T)
    print(p1[[1]])

    # Normalization
    dl <- quantile_norm(dl)

    # t-SNE with normalization:
    dl <- liger::runTSNE(dl)
    plots <- plotByDatasetAndCluster(dl, return.plots=T) # , clusters=as.character(rna.df$broad.cell.type) )
    print(plots[[1]])

    gplot = p1[[1]] + plots[[1]]
    ggsave(paste0(traitdir, 'gene_NMF_tSNE_', gene, '_w', window, '.png'), gplot, dpi=350,units='in', width=8, height=4)


    # Color each enhancer (will be to link) by only H3K27ac
    group.max = apply(sub.avg.mmat, 1, which.max)
    group.col = odf$COLOR[group.max]
    group.nam = odf$GROUP[group.max]

    # TODO: Experiment with H3K27ac or H3K4me1 or combined scores
    THRESHOLD = 0.4
    # plot(sub.ovldf$corr, sub.ovldf$corr_k4, pch=19,
    #      col=ifelse(sub.ovldf$corr < THRESHOLD, 'grey50','indianred'))
    # abline(h=0); abline(v=0)
    linkdf = sub.ovldf[sub.ovldf$corr > THRESHOLD, c('gene','name', 'mid', 'corr')]
    linkdf$gene = as.character(linkdf$gene)
    linkdf$tss = gloc[linkdf$gene]
    # linkdf$mid = kept.locdf[linkdf$name, 'mid']
    linkdf$id = kept.locdf[linkdf$name, 'id']
    linkdf$color = group.col[linkdf$id]
    linkdf$nam = group.nam[linkdf$id]
    linkdf$symbol = as.character(sapply(linkdf$gene, function(x){genemap$symbol[genemap$gene == x]}))
    write.table(linkdf, paste0('linking_data/corr_links/gene_links_rawcorr_H3K27ac_', gene, '_w', window,'.tsv'), sep="\t", quote=F, row.names=F)

    plot.red = TRUE
    if (plot.red){
        add.groups = c('Heart','Sm. Muscle','Endothelial','HSC & B-cell')
        linkdf = linkdf[linkdf$nam %in% add.groups,]
    }

    # Calc link parameters:
    link.center = with(linkdf, (tss + mid) / 2)
    link.radius = abs(link.center - linkdf$tss)
    # GWAS SNPs:
    gwloc = sub.gwdf$chromStart[sub.gwdf$chr == chrom]
    # TODO: Find linked SNPs - for each gene (using gwdf)
    
    # TODO: FINE MAP each SNP to the closest RELEVANT ENHANCER 

    # Set up the data limits:
    xlim = with(infoline, c(tss - window /2 , tss + window / 2))
    ylim.atac = c(0, max(sub.avg.mmat))

    # To convolve signal slightly
    # yconv = c(0.01, .05, .15, .30, .4, .30, .15, .05, 0.01) / .4
    # xconv = c(-200, -150, -100, -50, 0, 50, 100, 150, 200)
    # For axis:
    step = ((window * 2) / 5)
    mx = (floor(xlim[1] / step) + 1) * step
    atpar = seq(mx, xlim[2], by=step)
    lbls = comma_format()(atpar)

    add.genes = TRUE
    if (add.genes){
        genesuf = '_with_tx'
        sub.txdf = aggregate(txid ~ ., txdf[txdf$symbol %in% gnames,c('chr','start','end','strand','txid', 'symbol'),], function(x){x[1]})
        # Pad for name or other
        sub.txdf$pad = 40000
        sub.txdf = within(sub.txdf, pad.start <- start - pad * (strand == '+'))
        sub.txdf = within(sub.txdf, pad.end <- end + pad * (strand == '-'))
        sub.txdf = sub.txdf[order(sub.txdf$pad.start),]
        # Figure out how many lines we need:
        stgr = with(sub.txdf, GRanges(seqnames=chr, IRanges(pad.start, pad.end)))
        cvg = coverage(stgr)
        nlines = max(cvg)
        genes.height = .5 * nlines + .1
        # Order of genes:
        tx.ovl = data.frame(findOverlaps(stgr, stgr))
        tx.ovl = tx.ovl[tx.ovl$queryHits != tx.ovl$subjectHits,]
        # tx.ovl$qn = sub.txdf$symbol[tx.ovl$queryHits]
        NGENES = nrow(sub.txdf)
        sub.txdf$line = rep(0, NGENES)
        poss.lines = rev(1:nlines)
        for (i in 1:NGENES){
            if (i %in% tx.ovl$queryHits){
                pdf = tx.ovl[tx.ovl$queryHits == i,]
                comp = sort(pdf[,2])
                notline = c()
                for (j in comp){
                    if (i > j){
                        notline = c(notline, sub.txdf$line[j])
                    }
                }
                if (length(notline) > 0){
                    l = poss.lines[!(poss.lines %in% notline)][1]
                } else { l = nlines }
                sub.txdf$line[i] = l
            } else {
                sub.txdf$line[i] = 1
            }
        }
        sub.exdf = merge(exdf, sub.txdf[,c('txid','symbol', 'line')])
        sub.uxdf = merge(uxdf, sub.txdf[,c('txid','symbol', 'line')])
    } else { genesuf = ''; genes.height = 0}

    auto.addgroups = (gene != 'JCAD')
    if (auto.addgroups){
        # mtdf = aggregate(corr ~ nam, linkdf, max)
        if (nrow(linkdf) > 0){
            mtdf = aggregate(corr ~ nam, linkdf, length)
            mtdf = mtdf[order(mtdf$corr, decreasing=T),]
            add.groups = mtdf$nam[1:8]
            add.groups = add.groups[!is.na(add.groups)]
        } else { add.groups = c() }
    } else {
        add.groups = c('Heart','Sm. Muscle','Muscle','Endothelial','Stromal', 'Kidney', 'Liver', 'HSC & B-cell', 'Brain')
        add.groups = c('Heart','Sm. Muscle','Muscle','Endothelial','Kidney', 'Liver', 'HSC & B-cell')
    }
    if (length(add.groups) < 8){
        add.groups = unique(c(add.groups, 'HSC & B-cell','Muscle','Liver','Brain','Kidney','Pancreas','Heart', 'Stromal'))[1:8]
    }


    if (plot.red){
        add.groups = c('Heart','Sm. Muscle','Endothelial','HSC & B-cell')
        genesuf = paste0(genesuf, '_reduced')
    }
    NCT = length(add.groups)

    top.arc = TRUE
    plot.width=8
    plot.height=1 + 3 * (NCT / 8) + 1 * add.genes * genes.height / 1.5
    plot.height=2 / 1.2
    png(paste0(traitdir, 'gene_link_NMF_vis_', gene, '_w', window, genesuf, '.png'), units='in', res=450, width=plot.width, height=plot.height)
    # arc.sect=3
    arc.sect=2.5
    if (top.arc){
        lmat = matrix(c(NCT+2 + 1 * add.genes, 1:(NCT + 1 + add.genes * 1)), nrow=NCT + 2 + 1 * add.genes, ncol=1)
        heights = c(arc.sect, rep(1, NCT + 1))
    } else {
        lmat = matrix(1:(NCT + 1), nrow=NCT + 1 + 1 * add.genes, ncol=1)
        heights = c(rep(1, NCT), arc.sect)
    }
    if (add.genes) { heights = c(heights, genes.height) }
    layout(lmat, heights=heights)
    par(xaxs='i')
    sp=0
    par(mar=rep(sp, 4))
    for (i in 1:NCT){ 
        sgroup = add.groups[i]
        x = kept.locdf$mid
        y = c(sub.avg.mmat[,sgroup])
        ylim = ylim.atac
        plot(x, y, type='n', axes=F, ylab='', xlab='', ylim=ylim, xlim=xlim)
        # Add highlights:
        if (i == NCT){ 
            axis(1, at=atpar, labels=rep('', length(atpar)), lwd=.25)
            text(x=atpar, y=parpos(2, .7 - .2 * plot.red), labels=lbls, xpd=NA, cex=.7)
        }
        # Add the highlighted linked locations:
        tol = 500
        if (nrow(linkdf) > 0){
            link.atac = linkdf$mid
            rect(xleft=link.atac - tol, xright=link.atac + tol, 
                 ybottom=ylim.atac[1], ytop=ylim.atac[2], col='grey90', border=NA)
        }
        # abline(v=gwloc, col='slateblue', lwd=.5, lty='dashed')
        diffpar = seq(-10, 10, by=2.5) * 1e5 
        abline(v=infoline$tss + diffpar, lwd=.5, col='grey50', lty='dashed')
        abline(v=gloc, lwd=.5, col='red', lty='dashed')
        # Plot actual track:
        rect(xleft=x-tol, ybottom=0, xright=x+tol, ytop=y, col=colvals$group[sgroup], border=NA)
        text(x=xlim[1] + 0.028 * diff(xlim), y=mean(ylim.atac), sgroup, 
             font=2, cex=1, col=colvals$group[sgroup], adj=0)
        box(lwd=.25)
        arrowcol='grey25'
        if (i == 1){
            if (!(add.genes)) {
                text(gloc + -1 * diff(xlim) * 1e-3, parpos(2, -.4),
                     labels=gnames, adj=1, col=arrowcol, font=2, cex=.9)
            }
            for (dd in diffpar){
                tstr = paste(dd / 1e6,'Mb')
                adj = 1.2 * (dd > 0) - .1
                if (dd == 0){ tstr = 'TSS'; adj=-.05 }
                text(x=infoline$tss + dd, y = ylim[1] + .85 * (diff(ylim.atac)),
                     labels=tstr, cex=.5, adj=adj, col='grey30')
            }
        }
    }
    # TADs and SNPs on the bottom:
    par(mar=rep(sp, 4))
    plot(0, 1, xlim = xlim, ylim=c(0,1), type='n', axes=F)
    # Plot snps:
    gwind = which(gwloc >= infoline$tss - window/2 & gwloc <= infoline$tss + window/2)
    if (length(gwind) >  0){
        kept.gwloc = gwloc[gwind]
    } else {
        kept.gwloc = gwloc[gwloc >= infoline$tss - window]
    }
    rect(xleft=kept.gwloc-tol * 2, xright=kept.gwloc + tol * 2, 
         ybottom=0.8, ytop=1, border=NA, col='slateblue', lwd=1)
    tcex = 0.6
    tlab = paste(TRAIT, 'GWAS lead SNPs')
    twidth = strwidth(tlab, cex=tcex)
    tloc = xlim[1] + .01 * diff(xlim)
    text(tloc, y=.15,labels=tlab,
         font=2, cex=tcex, col='slateblue', adj=0)
    min.snp = min(kept.gwloc) # First possible snp.
    max.snp = max(kept.gwloc) # First possible snp.
    # Bottom:
    yloc = c(0.1, 0.75)
    segments(x0=tloc + twidth + .008 *  diff(xlim),
             x1=max.snp - 0.005 *  diff(xlim),
             y0=yloc[1],y1=yloc[1], col='slateblue', lwd=.5, xpd=TRUE)
    segments(x0=kept.gwloc - 0.005 *  diff(xlim), x1=kept.gwloc,
             y0=yloc[1],y1=yloc[2], col='slateblue', lwd=.5, xpd=TRUE)
    segments(x0=kept.gwloc - 0.001 * diff(xlim),
             x1=kept.gwloc + 0.001 * diff(xlim),
             y0=yloc[2],y1=yloc[2], col='slateblue', lwd=.5, xpd=TRUE)
    if (add.genes){
        par(mar=rep(sp, 4))
        plot(0, 1, xlim = xlim, ylim=c(.5,nlines + .5), type='n', axes=F)
        # box(lwd=.5)
        genecol = 'midnightblue'
        tpad = 1000
        for (i in 1:nrow(sub.txdf)){
            with(sub.txdf[i,], text(x=(start - tpad) * (strand == '+') + (end + tpad) * (strand == '-'),
                                y=line, labels=symbol, cex=.6, xpd=TRUE, adj=(strand == '+') * 1))
            # segments(x0=sub.txdf$start[i], x1=sub.txdf$end[i],
            #          y0=sub.txdf$line[i], y1=sub.txdf$line[i], 
            #          lwd=.25, col='slateblue')
            start = with(sub.txdf[i,], start * (strand == '+') + end * (strand == '-'))
            end = with(sub.txdf[i,], start * (strand == '-') + end * (strand == '+'))
            multarrows(x0=start, x1=end,
                       y0=sub.txdf$line[i], y1=sub.txdf$line[i], 
                       n_arr = round((sub.txdf$end[i] - sub.txdf$start[i]) / 5000),
                       length=.02,
                       lwd=.25, col='midnightblue')
        }
        lpad =0.2
        # Add exons:
        rect(xleft=sub.exdf$start, xright=sub.exdf$end,
             ybottom=sub.exdf$line - lpad, ytop=sub.exdf$line + lpad, 
             lwd=.25, col='midnightblue', border=NA)
        rect(xleft=sub.uxdf$start, xright=sub.uxdf$end,
             ybottom=sub.uxdf$line - lpad /2, ytop=sub.uxdf$line + lpad /2, 
             lwd=.25, col='midnightblue', border=NA)
    }
    # TODO: Add rsid
    if (nrow(linkdf) > 0){
        if (top.arc){ 
            par(mar=rep(sp, 4))
            plot(0, 1, xlim = xlim, ylim=c(0,1), type='n', axes=F)
        }
        asp.ratio = plot.width / (plot.height * (arc.sect / (NCT + add.genes * genes.height + top.arc + arc.sect)))
        mod.link.radius = abs(link.radius) / cos(pi/4)
        yvar = mod.link.radius / diff(xlim) * asp.ratio / 1.03
        # yvar = mod.link.radius / diff(xlim) * asp.ratio / (1.03 + 0.04 * add.genes * (genes.height / 1.5)**2)
        if (top.arc){ ang1 = pi/4 }
        draw.arc(x=link.center, y= 1 * (!top.arc) + (2 * (!top.arc) - 1) * (yvar) * (pi/4), xpd=TRUE, 
                 lwd=(linkdf$gene == ensg) * .5 + .5,
                 radius=mod.link.radius, angle1=pi * (!top.arc) + pi/4,angle2=pi * (!top.arc)+ pi - pi/4, col=linkdf$col)
    }
    dev.off()

}








