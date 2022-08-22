#!/usr/bin/R
# -----------------------------------------------
# Parallel to test_calculate_correlations.py
# Evaluate/debug data:
# -----------------------------------------------
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
cmd = paste('mkdir -p', imgdir, genedir)
system(cmd)
imgpref = paste0(imgdir, 'corrtest_')

# --------------------------
# Read in data mark matrix:
# TODO: Could read all marks
# --------------------------
lddir = 'linking_data/'
mark = 'H3K27ac'
chrom = 'chr10'
use.dense = TRUE
if (use.dense){
    runpref = paste0('dense_', chrom, '_')
    imgpref = paste0(imgpref, runpref)
    # Using chr10 for now:
    markfile = paste0(lddir, mark, '_all_bin_dense_on_mixed_impobs_chr10_r25_e100_merged.mtx.gz')
    mapfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv'
    locfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt'
    indfile = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_masterlist_indices.tsv'
} else {
    markfile = paste0(lddir, mark, '_all_bin_nonovl_on_mixed_impobs_r25_e100_allchr_merged.mtx.gz')
    mapfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_nonovl_any_coords_hg19_r200_e0_names.core.srt.tsv'
    locfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_nonovl_any_coords_hg19.core.srt.txt'
    indfile = 'masterlist_matindices/matindices_nonovl_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_masterlist_indices.tsv'
}
marknamesfile = paste0(lddir, 'mark_matrix_names.txt')
if (use.dense){
    markrdafile = paste0(lddir, mark, '_', runpref, 'mat_loc.Rda')
} else {
    markrdafile = paste0(lddir, mark, '_mat_loc.Rda')
}

# Load data:
if (!file.exists(markrdafile)){
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
    if (use.dense) {
        locdf = locdf[locdf$chr == 'chr10',]
        chrind = locdf$cls
        locdf$cls.old = locdf$cls
        locdf$cls = locdf$cls - min(locdf$cls.old) + 1
    }
    # Mark matrix:
    mmat = readMM(gzfile(markfile))
    mnames = scan(marknamesfile,'c')
    # Indices:
    matind = as.numeric(scan(indfile, 'c')) + 1
    if (use.dense) {
        # Subset + re-number matind for specific chromosome:
        matind = matind[matind %in% locdf$cls.old]
        matind = matind - min(locdf$cls.old) + 1
        mmat = as.matrix(mmat[matind,])
    } else {
        mmat = mmat[matind,]
    }
    gc()
    locdf = locdf[matind,]
    # Save for loading:
    save(mmat, mnames, locdf, matind, file=markrdafile)
} else {
    load(markrdafile)
}
gc()


# ---------------------
# Load expression data:
# ---------------------
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

# ------------------------------------------
# Find the samples in both + slice matrices:
# ------------------------------------------
# Make the overall matrix with ALL samples:
tform = make.tform(meta[mnames,'GROUP'], norm=TRUE, u=odf$GROUP)
avg.mmat = mmat %*% tform 

# Reduce to shared samples:
kept.samples = mnames[mnames %in% colnames(gmat)]
gmat = gmat[,kept.samples]
kept.idx = as.numeric(sapply(kept.samples, function(x){which(mnames == x)}))
mmat = mmat[,kept.idx]
mnames = mnames[kept.idx]
gc()

# Pre-compute the means for the matrices:
gmean = apply(gmat,1, mean)
mmean = apply(mmat,1, mean)

# ------------------------
# Make all testable pairs:
# ------------------------
window = 1e6
tssgr = GRanges(seqnames=tssdf$chr, IRanges(tssdf$tss - window, tssdf$tss + window), name=tssdf$gene)
locgr = GRanges(seqnames=locdf$chr, IRanges(locdf$start, locdf$end), name=locdf$name)

if (use.dense){
    ovlrdafile = paste0(lddir, mark, '_dense_chr10', '_', window, '_ovl_corr_df.Rda')
} else {
    ovlrdafile = paste0(lddir, mark, '_', window, '_ovl_corr_df.Rda')
}
if (!file.exists(ovlrdafile)){
    ovldf = data.frame(findOverlaps(tssgr, locgr))
    ovldf$dist = with(ovldf, locdf$mid[subjectHits] - tssdf$tss[queryHits])
    ovldf$corr = 0
    ovldf$dist.cut = cut(ovldf$dist, breaks=seq(-window, window, length.out=100), include.lowest=T)

    # Calculate the correlation
    chunksize = 1e5
    NT = nrow(ovldf)
    nchunks = floor(NT / chunksize) + 1
    for (i in 1:nchunks){
        print(i)
        ind = ((i-1) * chunksize + 1):min((i * chunksize), NT)
        gind = ovldf[ind,'queryHits']
        mind = ovldf[ind,'subjectHits']
        gsub = sweep(gmat[gind,], 1, gmean[gind],'-')
        msub = sweep(as.matrix(mmat[mind,]), 1, mmean[mind],'-')
        ssG = apply(gsub^2, 1, sum)
        ssM = apply(msub^2, 1, sum)
        out = apply(gsub * msub,1,sum) / sqrt(ssG * ssM)
        out[is.na(out)] = 0
        ovldf$corr[ind] = out
    }
    save(ovldf, file=ovlrdafile)
    gc()
} else {
    load(ovlrdafile)
}

ovldf$name = locgr$name[ovldf$subjectHits]
ovldf$gene = as.character(tssgr$name[ovldf$queryHits])


# Overall histogram of correlations:
png(paste0(imgpref, 'corr_hist.png'), res=450, units='in', width=6, height=3)
par(mar=c(2,2,2,0))
hist(ovldf$corr, 50, col='darkgrey',border='white')
dev.off()

# Plot boxplots of correlation against distance...
gplot = ggplot(ovldf, aes(dist.cut, corr)) + 
    geom_boxplot(fill='grey75', outlier.size=.5) + 
    labs(x='Distance from TSS', y='Correlation') + 
    theme_minimal() + 
    theme(axis.text.x=element_text(angle=90, size=5))
ggsave(paste0(imgpref, 'dist_v_corr.png'), gplot, dpi=450, units='in', width=8, height=6)


# -------------------------
# Load the validation data:
# -------------------------
compare.pchic = FALSE
if (compare.pchic){
    # Some validation:
    pcfile = paste0(lddir, 'ActivePromoterEnhancerLinks_JavierrePChIC.tsv.gz')
    pchdf = read.delim(gzfile(pcfile), header=T)
    btdf = unique(pchdf[,c('baitChr','baitSt','baitEnd','baitID')])
    oedf = unique(pchdf[,c('oeChr','oeSt','oeEnd','oeID')])

    # TSS with no window + PCHiC GR objects
    tssgr = GRanges(seqnames=tssdf$chr, IRanges(tssdf$tss, tssdf$tss), name=tssdf$gene)
    btgr = GRanges(seqnames=btdf$baitChr, IRanges(btdf$baitSt, btdf$baitEnd), name=btdf$baitID)
    oegr = GRanges(seqnames=oedf$oeChr, IRanges(oedf$oeSt, oedf$oeEnd), name=oedf$oeID)

    # Overlaps for mapping both names and genes to pchic data:
    btovl = data.frame(findOverlaps(tssgr, btgr))
    oeovl = data.frame(findOverlaps(locgr, oegr))
    btovl$gene = tssgr$name[btovl$queryHits] 
    btovl$baitID = btgr$name[btovl$subjectHits] 
    oeovl$name = locgr$name[oeovl$queryHits] 
    oeovl$oeID = oegr$name[oeovl$subjectHits] 
    pchdf = merge(merge(pchdf, btovl[,c('gene','baitID')]), oeovl[,c('name','oeID')])

    red.pchdf = unique(pchdf[,c('gene','name')])
    red.pchdf$is.pch = 1
    red.pchdf$gene = as.character(red.pchdf$gene)
    red.pchdf$name = as.character(red.pchdf$name)

    # Merge into ovldf (only the genes w/ bait):
    sub.ovldf = ovldf[ovldf$gene %in% red.pchdf$gene,]
    red.pchdf$ng = with(red.pchdf, paste0(gene,"_", name))
    sub.ovldf$ng = with(sub.ovldf, paste0(gene,"_", name))

    # Much much faster than merging:
    sub.ovldf$is.pch = 0
    sub.ovldf$is.pch[sub.ovldf$ng %in% red.pchdf$ng] = 1
    print(paste(sum(sub.ovldf$is.pch), nrow(sub.ovldf), 
                round(sum(sub.ovldf$is.pch) / nrow(sub.ovldf), 4)))

    gplot = ggplot(sub.ovldf, aes(corr)) + 
        facet_wrap(~is.pch, scales='free_y') + 
        geom_histogram() + 
        labs(x='Correlation', y='Count') + 
        theme_minimal()
    ggsave(paste0(imgpref, 'corr_by_ispchic.png'), gplot, dpi=450, units='in', width=6, height=3)


    gplot = ggplot(sub.ovldf, aes(corr, fill=factor(is.pch))) + 
        # facet_wrap(~is.pch, scales='free_y') + 
        geom_density(alpha=.5) + 
        labs(x='Correlation' ) + 
        theme_minimal()
    ggsave(paste0(imgpref, 'corr_by_ispchic_density.png'), gplot, dpi=450, units='in', width=6, height=3)

    # Plot best - shows that you can't validate with PCHiC because oe end is way too long
    # r2pch = merge(sub.ovldf[,c('gene','name','corr')], pchdf)
    # r2agg = aggregate(corr ~ oeID, r2pch, max)
}

# ---------------------------------------
# Look at neighborhood of specific genes:
# ---------------------------------------
genemap = read.delim(gmapfile, header=F)
names(genemap) = c('gene','symbol')
anno = merge(tssdf, genemap) # Note: not same as tssdf file.

gene = 'JCAD'
ensg = as.character(anno$gene[anno$symbol == gene])
infoline = anno[anno$symbol == gene,]
sub.gmat = gmat[ensg,]

# Matrix:
sub.ovldf = ovldf[ovldf$gene == ensg,]
sub.ovldf = merge(sub.ovldf, locdf[,c('name','mid')])
sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]
sub.mmat = as.matrix(mmat[sub.ovldf$subjectHits,])
sub.avg.mmat = avg.mmat[sub.ovldf$subjectHits,]
colnames(sub.mmat) = mnames

tform = make.tform(meta[mnames,'GROUP'], norm=TRUE, u=odf$GROUP)
sub.avg.expr = sub.gmat %*% tform
sub.avg.expr[is.na(sub.avg.expr)] = 0
ord = order(sub.avg.expr)
t(sub.avg.expr)[ord,]


# sub.mmat[sub.mmat > 10] = 10

# TODO: Should we quantile norm?
# Figure out which is top epigenome?
# tform = make.tform(meta[colnames(sub.mmat),'GROUP'], norm=TRUE, u=odf$GROUP)
# sub.avg.mmat = sub.mmat %*% tform 
# sub.avg.mmat[is.na(sub.avg.mmat)] = 0
group.max = apply(sub.avg.mmat, 1, which.max)
group.col = odf$COLOR[group.max]
group.nam = odf$GROUP[group.max]

# plot(sub.ovldf$corr, at=sub.ovldf$mid, pch=19)
# plot(sub.ovldf$mid, sub.ovldf$corr, type='l', xlab='', ylab ='', lwd=.5)

png(paste0(imgpref, 'gene_corr_', gene, '.png'), res=450, units='in', width=8, height=3)
par(xaxs='i')
# par(yaxs='i')
sp=.1
par(mar=c(3,3,sp,sp))
# plot(sub.ovldf$mid, sub.ovldf$corr, type='l', xlab='', ylab ='', lwd=.5)
x = sub.ovldf$mid
y = sub.ovldf$corr
plot(x, y, type='n')
barW = 500
rect(xleft=x-barW, ybottom=0, xright=x+barW, ytop=y, col='gray25', border=NA)
abline(v=infoline$tss, lwd=1, col='red')
ind = sub.ovldf$corr > 0.3
points(sub.ovldf$mid[ind], sub.ovldf$corr[ind], col=group.col[ind], pch=19, cex=.5)
ind2 = which(ind)[group.nam[ind] == 'Heart']
points(sub.ovldf$mid[ind2], sub.ovldf$corr[ind2], col='black', lwd=.5, bg=group.col[ind2], pch=21, cex=.5)
ind3 = sub.ovldf$corr > 0.4
text(sub.ovldf$mid[ind3] + barW * 20, sub.ovldf$corr[ind3], labels=group.nam[ind3], adj=0, cex=.25)
abline(h=0, lwd=.5, col='grey50', lty='dashed')
text(x=infoline$tss - 1e4, y=-.2,labels='TSS\n(JCAD)', adj=1, col='red')
mtext('Correlation', side=2, line=2)
mtext('chr10 Location (bp)', side=1, line=2)
box()
dev.off()



# TODO: ADD GROUPS AND ADD GWAS SIGNAL

add.groups = c('Heart','Sm. Muscle','Muscle','Endothelial','Stromal', 'Kidney', 'Liver', 'HSC & B-cell', 'Brain')
NG = length(add.groups)

w = 8
h = (3 + NG) * .3
apply(sub.avg.mmat[,add.groups], 2, max)
print(range(y))

# Load in the gwas file, pruned - load CAD:
gwrdafile = '~/EPIMAP_ANALYSIS/db/gwascatalog_may03_2019_noquotes.Rda'
load(gwrdafile)
gwssdf = gwssdf[grep("Coronary artery disease", gwssdf$trait),]
gwssdf = gwssdf[gwssdf$sampsize > 10000,]
gwssdf[gwssdf$pubMedID == '29212778',]
gwdf = unique(merge(gwdf, gwssdf)[,c('chrom','chromStart','pValue')])
gwdf = gwdf[order(gwdf$chrom),]
gwdf = gwdf[gwdf$pValue < 1e-8,]
gwdf = gwdf[gwdf$chrom == '10',]


png(paste0(imgpref, 'gene_corr_', gene, '_withtracks.png'), res=450, units='in', width=w, height=h)
par(xaxs='i')
layout(matrix(1:(NG + 2), nrow=(NG + 2)), heights=c(3, rep(1, NG), .5), TRUE)
par(yaxs='i')
sp=.1; rsp = 3
barW = 500
par(mar=c(sp, rsp,sp,sp))
x = sub.ovldf$mid
y = sub.ovldf$corr
ylim = range(y) * 1.05
plot(x, y, type='n', axes=F, ylim=ylim)
step = 0.2
atpar = c(rev(seq(0, min(ylim), by=-step)[-1]), c(seq(0, max(ylim), by=step)))
axis(2, at=atpar, lwd=.5, labels=FALSE)
text(x=par()$usr[1] - .01 * diff(par()$usr[1:2]),
     y=atpar, labels=atpar, srt=0, adj=1, xpd=TRUE, cex=.8)
ind = sub.ovldf$corr > 0.3
ind2 = which(ind)[group.nam[ind] == 'Heart']
ind3 = sub.ovldf$corr > 0.4
# Highlights
gwloc = gwdf$chromStart
rect(xleft=gwloc - barW * 5, ybottom=ylim[1], xright=gwloc+barW * 5, ytop=ylim[2], col=tsp.col('slateblue'), border=NA)
rect(xleft=x[ind3]-barW * 5, ybottom=ylim[1], xright=x[ind3]+barW, ytop=ylim[2], col=tsp.col('gray80'), border=NA)
rect(xleft=x-barW, ybottom=0, xright=x+barW, ytop=y, col='gray25', border=NA)
abline(v=infoline$tss + seq(-10, 10, by=2.5) * 1e5, lwd=.25, col='grey50', lty='dashed')
abline(v=infoline$tss, lwd=1, col='red')
points(sub.ovldf$mid[ind], sub.ovldf$corr[ind], col=group.col[ind], pch=19, cex=.5)
points(sub.ovldf$mid[ind2], sub.ovldf$corr[ind2], col='black', lwd=.5, bg=group.col[ind2], pch=21, cex=.5)
text(sub.ovldf$mid[ind3], sub.ovldf$corr[ind3], labels=group.nam[ind3], adj=-.25, cex=.25)
abline(h=0, lwd=.5, col='grey50', lty='dashed')
text(x=infoline$tss - 1e4, y=-.2,labels='TSS\n(JCAD)', adj=1, col='red', cex=.9)
mtext('Correlation', side=2, line=2, cex=.7)
# mtext('chr10 Location (bp)', side=1, line=2)
box(lwd=.5)
for (i in 1:NG){
    sgroup = add.groups[i]
    y = c(sub.avg.mmat[,sgroup])
    ylim = c(0, max(y))
    plot(x, y, type='n', axes=F, ylab='', xlab='', ylim=ylim)
    rect(xleft=x[ind3]-barW * 5, ybottom=ylim[1], xright=x[ind3]+barW, ytop=ylim[2], col=tsp.col('gray80'), border=NA)
    ind4 = which(ind3)[group.nam[ind3] == sgroup]
    if (length(ind4) > 0){
        rect(xleft=x[ind4]-barW * 5, ybottom=ylim[1], xright=x[ind4]+barW, ytop=ylim[2], col=tsp.col('yellow'), border=NA)
    }
    rect(xleft=x-barW, ybottom=0, xright=x+barW, ytop=y, col=colvals$group[sgroup], border=NA)
    text(x=range(x)[1] + 0.05 * diff(range(x)), y=mean(range(y)), sgroup, 
         font=2, cex=1, col=colvals$group[sgroup])
    abline(v=infoline$tss + seq(-10, 10, by=2.5) * 1e5, lwd=.25, col='grey50', lty='dashed')
    abline(v=infoline$tss, lwd=.25, col='red', lty='dashed')
    box(lwd=.25)
}
plot(x, rep(1, length(x)), type='n', axes=F, ylab='', xlab='', ylim=c(0,1))
rect(xleft=gwloc - barW * 5, ybottom=ylim[1], xright=gwloc+barW * 5, ytop=ylim[2], col=tsp.col('slateblue'), border=NA)
dev.off()



# -----------------------------------------
# Load in the gwas file, pruned - load CAD:
# -----------------------------------------
gwrdafile = '~/EPIMAP_ANALYSIS/db/gwascatalog_may03_2019_noquotes.Rda'
load(gwrdafile)
gwssdf = gwssdf[gwssdf$sampsize > 10000,]
sub.gwssdf = gwssdf[grep("Coronary artery disease", gwssdf$trait),]
sub.gwssdf[sub.gwssdf$pubMedID == '29212778',]
sub.gwdf = unique(merge(gwdf, sub.gwssdf)[,c('chrom','chromStart','pValue')])
sub.gwdf = sub.gwdf[order(sub.gwdf$chrom),]
sub.gwdf = sub.gwdf[sub.gwdf$pValue < 1e-8,]
sub.gwdf = sub.gwdf[sub.gwdf$chrom == '10',]
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

# ----------------------------------------------------------------
# Make nicer vis for the gene, adding other genes in neighborhood:
# ----------------------------------------------------------------
genemap = read.delim(gmapfile, header=F)
names(genemap) = c('gene','symbol')
anno = merge(tssdf, genemap) # Note: not same as tssdf file.

# genes:
chrom.genes = as.character(anno$symbol[anno$chr == chrom])
gene = 'JCAD'
for (k in 1:length(chrom.genes)){
    gene = chrom.genes[k]
    print(paste(k, gene))
    ensg = as.character(anno$gene[anno$symbol == gene])
    infoline = anno[anno$symbol == gene,]

    window = 1e6
    # Get ovl genes:
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

    # Get all links:
    sub.ovldf = ovldf[ovldf$gene %in% kept.gene,]
    sub.ovldf = merge(sub.ovldf, locdf[,c('name','mid')])
    sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]
    sub.ovldf$dist = sub.ovldf$mid - infoline$tss 
    sub.ovldf = sub.ovldf[abs(sub.ovldf$dist) < (window / 2),]
    sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]
    kept.locdf = unique(sub.ovldf[,c('mid','subjectHits', 'name')])
    kept.dhs = kept.locdf$subjectHits
    kept.mid = kept.locdf$mid
    rownames(kept.locdf) = kept.locdf$name
    kept.locdf$id = 1:nrow(kept.locdf)

    # Get locations:
    sub.mmat = as.matrix(mmat[kept.dhs,])
    sub.avg.mmat = avg.mmat[kept.dhs,]
    colnames(sub.mmat) = mnames

    # Average expression:
    tform = make.tform(meta[mnames,'GROUP'], norm=TRUE, u=odf$GROUP)
    sub.gmat = gmat[kept.gene,]
    sub.avg.expr = sub.gmat %*% tform
    sub.avg.expr[is.na(sub.avg.expr)] = 0

    # Top by only H3K27ac
    group.max = apply(sub.avg.mmat, 1, which.max)
    group.col = odf$COLOR[group.max]
    group.nam = odf$GROUP[group.max]

    THRESHOLD = 0.4
    linkdf = sub.ovldf[sub.ovldf$corr > THRESHOLD, c('gene','name', 'mid', 'corr')]
    linkdf$tss = gloc[linkdf$gene]
    # linkdf$mid = kept.locdf[linkdf$name, 'mid']
    linkdf$id = kept.locdf[linkdf$name, 'id']
    linkdf$color = group.col[linkdf$id]
    linkdf$nam = group.nam[linkdf$id]
    linkdf$symbol = as.character(sapply(linkdf$gene, function(x){genemap$symbol[genemap$gene == x]}))
    write.table(linkdf, paste0('linking_data/corr_links/gene_links_rawcorr_H3K27ac_', gene, '_w', window,'.tsv'), sep="\t", quote=F, row.names=F)
    # Calc link parameters:
    link.center = with(linkdf, (tss + mid) / 2)
    link.radius = abs(link.center - linkdf$tss)
    gwloc = gwdf$chromStart

    # TODO: Find linked SNPs - for each gene (using gwdf)

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

    # GWAS SNPs:
    gwloc = sub.gwdf$chromStart

    add.genes = TRUE
    if (add.genes){
        genesuf = '_with_tx'
        sub.txdf = aggregate(txid ~ ., txdf[txdf$symbol %in% gnames,c('chr','start','end','strand','txid', 'symbol'),], function(x){x[1]})
        # TODO: GET UTR??
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
    NCT = length(add.groups)

    top.arc = TRUE
    plot.width=8
    plot.height=1 + 3 * (NCT / 8) + 1 * add.genes * genes.height / 1.5
    plot.height=2
    png(paste0(genedir, 'gene_link_vis_', gene, '_w', window, genesuf, '.png'), units='in', res=450, width=plot.width, height=plot.height)
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
            text(x=atpar, y=parpos(2, .7), labels=lbls, xpd=NA, cex=.7)
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
    # Add snps:
    if (gene == 'JCAD'){
        rect(xleft=gwloc-tol * 2, xright=gwloc + tol * 2, 
             ybottom=0.8, ytop=1, border=NA, col='slateblue', lwd=1)
        text(xlim[1] + .18 *  diff(xlim), 
             y=.15,labels='Coronary Artery Disease (CAD) GWAS lead SNPs',
             font=2, cex=.6, col='slateblue', adj=0)
    segments(x0=xlim[1] + .374 *  diff(xlim),
             x1=xlim[1] + .41 *  diff(xlim),
             y0=.15,y1=.15, col='slateblue', lwd=.5, xpd=TRUE)
    segments(x0=xlim[1] + .41 *  diff(xlim),
             x1=xlim[1] + .415 *  diff(xlim),
             y0=.15,y1=.75, col='slateblue', lwd=.5, xpd=TRUE)
    segments(x0=xlim[1] + .412 *  diff(xlim),
             x1=xlim[1] + .420 *  diff(xlim),
             y0=.75,y1=.75, col='slateblue', lwd=.5, xpd=TRUE)
    }
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









# ---------------------
# CODE FOR LINKING OLD:
# ---------------------
target.ind = which(genelocdf$symbol == gene)
geneinfo = genelocdf[target.ind,]
# Find overlapping peaks:
genetss <- with(geneinfo, GRanges(seqnames=chr, ranges=IRanges(start=tss - window, end=tss + window)))
tad.ovl = data.frame(findOverlaps(genetss, tad.gr)) # tads for plotting
atac.ovl = data.frame(findOverlaps(genetss, atac.gr))
geneloc.gr = with(genelocdf, GRanges(seqnames=chr, ranges=IRanges(start=tss, end=tss+1)))
anno.ovl = data.frame(findOverlaps(genetss, geneloc.gr))
kept.gene = which(gene.ind %in% anno.ovl$subjectHits)
kept.atac = which(atac.ind %in% atac.ovl$subjectHits)

sub.geneloc = genelocdf[anno.ovl$subjectHits,]
gloc = sub.geneloc$tss
gstrand = as.character(sub.geneloc$strand)
aloc = apply(atac.df[atac.ovl$subjectHits, c('start','end')], 1, mean)
sub.gmat = as.matrix(gmat[kept.gene, ])
sub.amat = as.matrix(amat[kept.atac, ])
tform = make.tform(sapply(colnames(gmat), function(x){sub("\\..*","", x)}), norm=TRUE)
if (ncol(sub.gmat) == 1){
    ms.gmat = t(sub.gmat) %*% tform 
} else { 
    ms.gmat = sub.gmat %*% tform 
}
ms.amat = sub.amat %*% tform 

# Get the links:
lmat = as.matrix(chrfit[kept.atac, kept.gene])
link.ind = which(lmat > 0, arr.ind=T)
link.str = lmat[link.ind]
link.center = (aloc[link.ind[,1]] + gloc[link.ind[,2]]) /2
link.radius = abs(link.center - gloc[link.ind[,2]])

sub.genes = as.character(sub.geneloc$symbol)
sub.qtl = merge(qtldf, sub.geneloc)
if (nrow(sub.qtl) > 0){
    qtl.center = (sub.qtl$loc + sub.qtl$tss)/2
    qtl.radius = abs(qtl.center - sub.qtl$tss)
}


loop.df$mid1 = apply(loop.df[,c('start1','end1')], 1, mean)
loop.df$mid2 = apply(loop.df[,c('start2','end2')], 1, mean)
sub.loop.df = loop.df[(loop.df$mid1 <= xlim[2] & loop.df$mid1 >= xlim[1]) | (loop.df$mid2 <= xlim[2] & loop.df$mid2 >= xlim[1]),]
sub.loop.df = sub.loop.df[sub.loop.df$chr == chr,]
if (nrow(sub.loop.df) > 0){
    print(gene)
    loop.center = (sub.loop.df$mid1 + sub.loop.df$mid2)/2
    loop.radius = abs(loop.center - sub.loop.df$mid1)
}

# Set up the data limits:
xlim = with(geneinfo, c(tss - window, tss + window))
ylim.atac = c(0, max(ms.amat))
ylim.rna = c(0, max(ms.gmat) * 1.25)

# To convolve signal slightly
yconv = c(0.01, .05, .15, .30, .4, .30, .15, .05, 0.01) / .4
xconv = c(-200, -150, -100, -50, 0, 50, 100, 150, 200)

# For axis:
step = ((window * 2) / 5)
mx = (floor(xlim[1] / step) + 1) * step
atpar = seq(mx, xlim[2], by=step)
lbls = comma_format()(atpar)

top.arc = TRUE
NCT = 8
plot.width=8
plot.height=4
png(paste0('~/linking_visualization_',gene,'_w',window,'.png'), units='in', res=450, width=plot.width, height=plot.height)
arc.sect=3
if (top.arc){
    lmat = matrix(c(NCT+2, 1:(NCT + 1)), nrow=NCT + 2, ncol=1)
    heights = c(arc.sect, rep(1, NCT + 1))
} else {
    lmat = matrix(1:(NCT + 1), nrow=NCT + 1, ncol=1)
    heights = c(rep(1, NCT), arc.sect)
}
layout(lmat, heights=heights)
par(xaxs='i')
sp=0
par(mar=rep(sp, 4))
for (ct in short.nam[c(1:6,8:9)]){
    ctcol = map.scell[ct,'col']
    plot(0, 1, xlim = xlim, ylim=ylim.rna, type='n', axes=F)
    text(xlim[1] + .005 *  diff(xlim), y=ylim.rna[1] * .5 + ylim.rna[2] * .5, 
         labels=map.scell[ct,'t3cell'],
         font=2, col=ctcol, cex=1.4, adj=c(0,.5))
    abline(h=0, lwd=.5)
    if (ct == short.nam[9]){ 
        axis(1, at=atpar, labels=rep('', length(atpar)), lwd=.5)
        text(x= atpar, y=parpos(2, .5), labels=lbls, xpd=NA)
    }
    # Add the highlighted linked locations:
    tol = 500
    if (nrow(link.ind) > 0){
        link.atac = aloc[link.ind[,1]]
        rect(xleft=link.atac - tol, xright=link.atac + tol, 
             ybottom=ylim.rna[1], ytop=ylim.rna[2], col='grey90', border=NA)
    }
    rect(xleft=gloc - tol, xright=gloc + tol, 
         ybottom=ylim.rna[1], ytop=ylim.rna[2], col=tsp.col('goldenrod1'), border=NA)
    # Plot atac and RNA:
    yscaled = ms.amat[,ct] * (ylim.rna[2] / ylim.atac[2])
    xmat = matrix(rep(aloc, length(xconv)), ncol=length(xconv))
    ymat = matrix(rep(yscaled, length(xconv)), ncol=length(xconv))
    xmat = sweep(xmat, 2, xconv, '+')
    ymat = sweep(ymat, 2, yconv, '*')
    xvec = unlist(xmat)
    yvec = unlist(ymat)
    ord = order(xvec)
    xvec = xvec[ord]
    yvec = yvec[ord]
    polygon(c(xlim[1], xvec, xlim[2]), c(0, yvec, 0), 
            col=tsp.col(ctcol), border=ctcol, xpd=TRUE, lwd=1)
    # RNA:
    ashift = ms.gmat[,ct]
    ashift[ashift < 15] = 15
    arrowshift = (2 * (gstrand == '+') - 1) * diff(xlim) * 0.001 * ashift
    arrowcol='grey25'
    p.arrows2(x1=gloc, x2=gloc + arrowshift, y1=ms.gmat[,ct], y2=ms.gmat[,ct], 
              pch=19, col=arrowcol, fill=arrowcol, border=arrowcol, lwd=1, size=.4)
    segments(x0=gloc, x1=gloc, y0=0, y1= ms.gmat[,ct], pch=19, col=arrowcol, lwd=1)
    if (ct == short.nam[1]){
        abline(h=ylim.rna[2], col='grey25', lwd=.5)
        nind = which(gstrand == '-')
        pind = which(gstrand == '+')
        if (length(nind) > 0){
            text(gloc[nind] + -1 * diff(xlim) * 5e-3, mean(ylim.rna), labels=sub.geneloc$symbol[nind],
                 adj =1, col=arrowcol, font=2, cex=1)
        }
        if (length(pind) > 0){
            text(gloc[pind] + 1 * diff(xlim) * 5e-3, mean(ylim.rna), labels=sub.geneloc$symbol[pind],
                 adj =0, col=arrowcol, font=2, cex=1)
        }
    }
}
plot(0, 1, xlim = xlim, ylim=c(0,1), type='n', axes=F)
kept.tad = tad.df[tad.ovl$subjectHits,]
if (nrow(kept.tad) > 0){
    rect(xleft=kept.tad$start, xright=kept.tad$end, 
         ybottom=.05, ytop=.25, 
         border=NA, col='grey80', lwd=1)
}
# Add snps:
rect(xleft=sub.gwdf$chromStart-tol, 
     xright=sub.gwdf$chromStart+tol, 
     ybottom=00, ytop=.30, 
     border=NA, col='slateblue', lwd=1)
text(xlim[1] + .005 *  diff(xlim), 
     y=.15,labels='Schizophrenia PGC lead SNPs:',
     font=2, cex=.7, col='slateblue', adj=0)
if (nrow(link.ind) > 0){
    if (top.arc){ plot(0, 1, xlim = xlim, ylim=c(0,1), type='n', axes=F)}
    if (nrow(sub.qtl) > 0 & plot.eqtl){
        mod.qtl.radius = qtl.radius / cos(pi/4)
        qtl.yvar = mod.qtl.radius / diff(xlim) * asp.ratio / 1.03
        draw.arc(x=qtl.center, y= 1 * (!top.arc) + (2 * (!top.arc) - 1) * (qtl.yvar) * (pi/4), xpd=TRUE,
                 radius=mod.qtl.radius, angle1=pi * (!top.arc) + pi/4,angle2=pi * (!top.arc)+ pi - pi/4, col=tsp.col('grey80'))
        legend('right', legend=c('Raw GTEx v7 eQTL - Frontal Cortex', 'Predicted Enhancer-Gene Links'),
               lwd=2,,bty='n', col=c('grey80','slateblue'))
    }
    asp.ratio = plot.width / (plot.height * (arc.sect / (NCT + top.arc + arc.sect)))
    mod.link.radius = abs(link.center - gloc[link.ind[,2]]) / cos(pi/4)
    yvar = mod.link.radius / diff(xlim) * asp.ratio / 1.03
    if (top.arc){ 
        ang1 = pi/4
    }
    draw.arc(x=link.center, y= 1 * (!top.arc) + (2 * (!top.arc) - 1) * (yvar) * (pi/4), xpd=TRUE,
             radius=mod.link.radius, angle1=pi * (!top.arc) + pi/4,angle2=pi * (!top.arc)+ pi - pi/4, col='slateblue')
    if (nrow(sub.loop.df) > 0){
        print(gene)
        mod.loop.radius = loop.radius / cos(pi/4)
        loop.yvar = mod.loop.radius / diff(xlim) * asp.ratio / 1.03
        draw.arc(x=loop.center, y= 1 * (!top.arc) + (2 * (!top.arc) - 1) * (loop.yvar) * (pi/4), xpd=TRUE,
                 radius=mod.loop.radius, angle1=pi * (!top.arc) + pi/4,angle2=pi * (!top.arc)+ pi - pi/4, col='indianred')
    }
}
dev.off()









# ---------------------------------
# Testing the gene expression data:
# ---------------------------------
# Run PCA on the log2 FPKM matrix:
norm.gmat = sweep(gmat, 2, colMeans(gmat),'-')
pc = prcomp(t(norm.gmat), rank=10)
# pc = prcomp(t(gmat), rank=10)
dim(pc$x)
bcol = meta[colnames(gmat), 'COLOR']


library(uwot)

pc = prcomp(t(norm.gmat), rank=50)
ux = umap(pc$x, n_neighbors=40)
plot(ux[,1], ux[,2], col = bcol, main = "UMAP", xlab = "UMAP 1", ylab = "UMAP 3", pch=19)

dev.new()
pc = prcomp(t(gmat), rank=50)
ux = umap(pc$x, n_neighbors=40)

png(paste0(imgpref, 'samples_expr_umap.png'), res=450, units='in', width=5, height=5)
par(mar=rep(1,4))
# plot(ux[,1], ux[,2], col = gcol, main = "", xlab="", ylab="", # pch=19, cex=.2, axes=F)
plot(ux[,1], ux[,2], col = bcol, main = "", xlab = "", ylab = "", pch=19, axes=F)
mtext("UMAP 1",side=1)
mtext("UMAP 2",side=2)
mtext("UMAP on Gene Expression - Samples", side=3, line=0, cex=1)
dev.off()



# Plot genes:
# norm.gmat = sweep(gmat, 1, rowMeans(gmat),'-')
ux = umap(gmat, n_neighbors=50)
gcol = bcol[apply(norm.gmat,1, which.max)]

png(paste0(imgpref, 'genes_umap.png'), res=450, units='in', width=6, height=6)
par(mar=rep(1,4))
plot(ux[,1], ux[,2], col = gcol, main = "", xlab="", ylab="", 
     pch=19, cex=.2, axes=F)
mtext("UMAP 1",side=1)
mtext("UMAP 2",side=2)
mtext("UMAP on Gene Expression - Genes", side=3, line=0, cex=1)
dev.off()





# ve = pc$sdev^2 / sum(pc$sdev^2)
# plot(1:10, ve[1:10] * 100, type='b', xlab='PC', ylab='% variance explained', ylim=c(0,100))
# plot(1:50, cumsum(ve[1:50] * 100), type='b', xlab='PC', ylab='Cummulative % variance explained', ylim=c(0,100))

plot(pc$x[, 1], pc$x[, 2], col = bcol, main = "PCA", xlab = "PC1", ylab = "PC2", pch=19)
plot(pc$x[, 1], pc$x[, 3], col = bcol, main = "PCA", xlab = "PC1", ylab = "PC3", pch=19)
plot(pc$x[, 2], pc$x[, 3], col = bcol, main = "PCA", xlab = "PC2", ylab = "PC3", pch=19)


# pcg$

plot(sub.ovldf$dist, sub.ovldf$corr, type='l')
abline(0,0)





