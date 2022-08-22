#!/usr/bin/R
# ---------------------------------------
# Make correlation based FDR predictions:
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
library(PRROC)
library(huxtable)

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
    locmids = locdf$mid  # Keep mids for distances
    rm(locmap)
    gc()
    # Matrix names + indices:
    mnames = scan(marknamesfile,'c')
    matind = as.numeric(scan(indfile, 'c')) + 1
    # Mapping 3.6M to 2M
    enhmap = rep(0,nrow(locdf))
    enhmap[matind] = 1:length(matind)
    locdf = locdf[matind,]
    # Save for loading:
    save(mnames, locdf, matind, enhmap, locmids, file=dhsrdafile)
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


# -------------------------------------------------------
# Load in validation datasets into genome ranges objects:
# -------------------------------------------------------
vdir = paste0(lddir, 'validation/')
bnames = list.files(path=vdir, pattern='.*Benchmark.*bed')
blist = list()
for (bfile in bnames){
    btype = sub("GM12878.","",sub("-Benchmark.*","", bfile))
    btab = read.delim(paste0(vdir, bfile), header=F)
    names(btab) = c('chr','start','end','gene','score')
    blist[[btype]] = with(btab, GRanges(seqnames=chr, IRanges(start=start, end=end), name=gene, score=score))
    rm(btab)
}

# -------------------------------------
# Read in the precomputed correlations:
# -------------------------------------
# Data files:
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
fnames = list.files(path=lddir, pattern='*_precomputed_corr.hdf5')

# ------------------------------------------------
# Create some predictions for all chr for BSS00439
# ------------------------------------------------
sample = 'BSS00439'
sampind = which(mnames == sample)
# sampind = 200
meta[mnames[sampind],]
# 1. FDR corr?
# 2. Distance metrics, etc.
# 3. Some mixture

# Get strength for this sample:
# H3K27ac
h5f = H5Fopen(H3K27ac.file)
h5d = h5f&"matrix"
h27samp = h5d[sampind,]
H5Dclose(h5d)
H5Fclose(h5f)

# H3K4me1
h5f = H5Fopen(H3K4me1.file)
h5d = h5f&"matrix"
h4samp = h5d[sampind,]
H5Dclose(h5d)
H5Fclose(h5f)


# Read in the roadmap links:
rddf = c()
for (state in c('6','7','12')){
    rddf = rbind(rddf, read.delim(paste0(lddir, 'roadmap_links/links_E116_', state,'_2.5.txt'), header=F))
}
names(rddf) = c('chrom','start','end','gene','score', 'bin')

# XGB alone, with dist, with mark, (with both)
xgl = list()
xgl[['xgb']] = read.delim(paste0(lddir, 'predictions/collated/BSS00439_collated_pred.tsv.gz'), sep="\t", header=F, stringsAsFactors=F)
xgl[['xgb_d50k']] = read.delim(paste0(lddir, 'predictions_d50000/BSS00439_collated_pred_d50000.tsv.gz'), sep="\t", header=F, stringsAsFactors=F)
xgl[['xgb_mark']] = read.delim(paste0(lddir, 'predictions_mark/BSS00439_collated_pred_mark.tsv.gz'), sep="\t", header=F, stringsAsFactors=F)
xgl[['xgb_mark_d50k']] = read.delim(paste0(lddir, 'predictions_mark_d50000/collated/BSS00439_collated_pred_mark_d50000.tsv.gz'), sep="\t", header=F, stringsAsFactors=F)
for (stype in names(xgl)){
    names(xgl[[stype]]) = c('chrom','start','end','gene','proba', 'state')
    xgl[[stype]]$dhsloc = with(xgl[[stype]], (end + start) / 2)
}

# # Alternatively: read some (full) xgboost predictions:
# xgdf = c()
# for (state in c('E7','E8', 'E9','E10', 'E11', 'E15')){
#     etab = read.delim(paste0(lddir, sample, '_',state,'_0.txt'), sep="\t", header=F, stringsAsFactors=F)
#     xgdf = rbind(xgdf, etab)
# }
# names(xgdf) = c('chrom','start','end','gene','proba')
# xgdf$dhsloc = (xgdf$end + xgdf$start) / 2

chrlist = paste0('chr', c(1:23,'X'))
for (chrom in chrlist){
    # Load links:
    print(paste("Loading links for", chrom))
    cfile = paste0(lddir, chrom, '_pairs_ENH_ovl_df.tsv.gz')
    # genes = sub.gwdf$gene[sub.gwdf$chr == chrom]
    cpairsdf = read.delim(cfile, header=T, sep="\t", stringsAsFactors=F)
    cpairsdf$mind = cpairsdf$mind + 1  # Was 0-indexed
    cpairsdf$chrom = chrom
    cind = cpairsdf$ind

    # Get matrix of correlations:
    mat = matrix(0, ncol=length(fnames), nrow=length(cind))
    marks = c()
    for (i in 1:length(fnames)){
        marks = c(marks, sub("_precomputed.*","",fnames[i]))
        print(marks)
        hfile = paste0(lddir, fnames[i])
        h5f = H5Fopen(hfile)
        mat[, i]= h5f$matrix[cind]
        H5Fclose(h5f)
    }
    colnames(mat) = marks

    mat[is.na(mat)] = 0 
    mcoeff = c(.5,1,1,-.5,-.5,.5)
    cscore = as.numeric(mat %*% as.matrix(mcoeff))

    # Raw values from matrix:
    c27vals = h27samp[matind[cpairsdf$mind]]
    c4vals = h4samp[matind[cpairsdf$mind]]
    # Get strengths, locdf, etc.
    cpairsdf$dhsloc = locdf$mid[cpairsdf$mind]
    cpairsdf$dist = abs(cpairsdf$dhsloc - cpairsdf$tss)

    tol = 100
    tol = 250
    cgr = with(cpairsdf, GRanges(seqnames=chrom, IRanges(start=dhsloc-tol, end=dhsloc+tol), name=gene))
    cpairsdf$acd7 = c27vals * cpairsdf$dist^-.7
    tot.acd7 = aggregate(acd7 ~ gene, cpairsdf, sum)
    names(tot.acd7)[2] = 'tot'
    rownames(tot.acd7) = tot.acd7$gene
    cp2df = cpairsdf[,c('ind','acd7','gene')]
    cp2df$tot = tot.acd7[cp2df$gene, 'tot']
    cp2df$score = cp2df$acd7 / cp2df$tot

    # Parse all the xgboost results:
    sgl = list()
    print("Parsing xgboost results")
    for (stype in names(xgl)){
        sgdf = xgl[[stype]][xgl[[stype]]$chrom == chrom,]
        sgdf = aggregate(proba ~ dhsloc + gene, sgdf, max)
        sgdf = merge(cpairsdf[,c('gene','ind', 'dhsloc')], sgdf, all.x=TRUE)
        sgdf = sgdf[order(sgdf$ind),]
        ind = which(is.na(sgdf$proba))
        sgdf$proba[ind] = 0
        print(paste(stype, sum(sgdf$proba > 0)))  # Output number of links
        sgl[[stype]] = sgdf$proba
    }

    # Roadmap scores (usually small subset)
    srddf = rddf[rddf$chrom == chrom,]
    srgr = with(srddf, GRanges(seqnames=chrom, IRanges(start=start, end=end), name=gene, score=score))
    srddf = as.data.frame(findOverlaps(srgr, cgr))
    srddf$rgene = srgr$name[srddf$queryHits]
    srddf$score = srgr$score[srddf$queryHits]
    srddf$gene = cgr$name[srddf$subjectHits]
    srddf$ind = cpairsdf$ind[srddf$subjectHits]
    srddf = srddf[srddf$rgene == srddf$gene,]
    srddf = aggregate(score ~ ind, srddf, max)
    srddf = merge(cpairsdf,srddf, all.x=TRUE)
    srddf = srddf[order(srddf$ind),]
    srddf$score[is.na(srddf$score)] = 0

    # Merge xgboost predictions into cpairsdf:
    # SCORES:
    # Distance + distance * strength:
    slist = list()
    # slist[['d7']] = cpairsdf$dist^-.7
    slist[['d1']] = cpairsdf$dist^-1
    slist[['rdmp']] = srddf$score
    # Add xgboost:
    for (stype in names(sgl)){ slist[[stype]] = sgl[[stype]] }

    # slist[['rdmp + d1']] = srddf$score + cpairsdf$dist^-1
    # slist[['xgboost x d1']] = (sgdf$proba / (1 - sgdf$proba)) * cpairsdf$dist^-1
    # xg2 = sgdf$proba
    # xg2[xg2 < 0.714] = 0 
    # slist[['xgboost + d1']] = (xg2 / (1 - xg2)) + cpairsdf$dist^-1
    # # slist[['xgboost + ac + d1']] = (xg2 / (1 - xg2)) + cpairsdf$dist^-1 + c27vals / max(c27vals)
    # xm2 = smdf$proba
    # xm2[xm2 < 0.714] = 0 
    # slist[['xgboost2 + d1']] = (xm2 / (1 - xm2)) + cpairsdf$dist^-1

    # dind = (cpairsdf$dist < 50000)
    # slist[['rdmp mixed']] = dind * cpairsdf$dist^-1 + (!dind) * srddf$score
    # slist[['xgboost mixed']] = dind * cpairsdf$dist^-1 + (!dind) * sgdf$proba
    # slist[['ac_d7']] = c27vals * cpairsdf$dist^-.7
    # slist[['ac_d1']] = c27vals * cpairsdf$dist^-1
    # slist[['me_d7']] = c4vals * cpairsdf$dist^-.7
    # slist[['me_d1']] = c4vals * cpairsdf$dist^-1
    # slist[['all_d7']] = (c27vals + c4vals) * (cpairsdf$dist^-.7)
    # slist[['all2_d7']] = (c27vals + c4vals)^2 * (cpairsdf$dist^-.7)
    # slist[['all_log_d7']] = log10(c27vals + c4vals + 1) * (cpairsdf$dist^-.7)
    # slist[['all_d1']] = (c27vals + c4vals) * (cpairsdf$dist^-1)
    # slist[['cscore_d1']] = cscore * cpairsdf$dist^-1
    # slist[['ac_d7_norm']] = cp2df$score
    # slist[['ac_alone']] = c27vals
    # slist[['me_alone']] = c4vals
    # slist[['ac_corr']] = mat[,'H3K27ac'] 
    # slist[['me_corr']] = mat[,'H3K4me1']
    # slist[['cscore']] = cscore 
    # slist[['ac_corr_d7']] = mat[,'H3K27ac'] * cpairsdf$dist^-.7
    # slist[['me_corr_d7']] = mat[,'H3K4me1'] * cpairsdf$dist^-.7
    # slist[['ac_corr_d1']] = mat[,'H3K27ac'] * cpairsdf$dist^-1
    # slist[['me_corr_d1']] = mat[,'H3K4me1'] * cpairsdf$dist^-1
    # slist[['ac2_corr_d1']] = mat[,'H3K27ac'] * cpairsdf$dist^-1 * c27vals
    # slist[['me2_corr_d1']] = mat[,'H3K4me1'] * cpairsdf$dist^-1 * c4vals

    sapply(slist, length)

    # Colors:
    snap.cols = scan('Annotation/snap_colors.tsv', 'c')

    # MERGE WITH VALIDATIONS
    aucdf = c()
    glist = list()
    gpr = list()
    for (btype in names(blist)){
        print(btype)
        bgr = blist[[btype]]
        bovl = as.data.frame(findOverlaps(bgr, cgr))
        ind = which(bgr$name[bovl$queryHits] == cgr$name[bovl$subjectHits])
        bovl = bovl[ind,]
        # BG and FG:
        fg = bovl$subjectHits[bgr$score[bovl$queryHits] == 1]
        bg = bovl$subjectHits[bgr$score[bovl$queryHits] == 0]
        print(paste0(length(fg), ', ', length(bg)))

        if (length(fg) > 0){
            curvedf = c()
            for (stype in names(slist)){
                score = slist[[stype]]
                score[is.na(score)] = 0
                roc = roc.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
                pr = pr.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
                print(paste0(stype, ': ', round(roc$auc,3), ", ", round(pr$auc.integral,3)))
                aucdf = rbind(aucdf, data.frame(auroc=round(roc$auc,4), auprc=round(pr$auc.integral,3),stype=stype,valid=btype))
                srdf = data.frame(x=roc$curve[,1], y=roc$curve[,2], stype=stype, curve='ROC')
                spdf = data.frame(x=pr$curve[,1], y=pr$curve[,2], stype=stype, curve='PR')
                curvedf = rbind(curvedf, rbind(srdf, spdf))
            }
        }

        gplot = ggplot(curvedf, aes(x,y,col=stype)) + 
            facet_wrap(~curve) + labs(title=btype) + 
            geom_line() + theme_pubr() + 
            scale_color_brewer(palette='Spectral', name='Method') + 
            theme(legend.position='right')
        ggsave(paste0(imgpref, btype, '_roc_pr_comp_', chrom, '.png'), gplot, units='in', dpi=450, width=9, height=5)
        glist[[btype]] = gplot

        gplot = ggplot(curvedf[curvedf$curve == 'PR',], aes(x,y,col=stype)) + 
            labs(title=btype, x='Recall',y='Precision') + 
            geom_line() + theme_pubr() + 
            # scale_color_brewer(palette='Spectral', name='Method') + 
            scale_color_brewer(palette='Paired', name='Method') + 
            theme(legend.position='right')
        ggsave(paste0(imgpref, btype, '_pr_comp_', chrom, '.png'), gplot, units='in', dpi=450, width=9, height=5)
        gpr[[btype]] = gplot
    }
    prtab = spread(aucdf[,c('auprc','stype','valid')], valid, auprc)
    roctab = spread(aucdf[,c('auroc','stype','valid')], valid, auroc)

    gplot = patchwork::wrap_plots(glist)
    ggsave(paste0(imgpref, 'allmetric_roc_pr_comp_', chrom, '.png'), gplot, units='in', dpi=450, width=18, height=8)

    # gplot = patchwork::wrap_plots(gpr)
    gplot = ggarrange(gpr[[1]], gpr[[2]], gpr[[3]],
                      gpr[[4]], gpr[[5]], gpr[[6]], ncol=3, nrow=2,
                      align='hv', common.legend = TRUE, legend = "bottom")
    ggsave(paste0(imgpref, 'allmetric_pr_comp_', chrom, '.png'), gplot, units='in', dpi=450, width=10, height=8)

    prtab = prtab[order(prtab$CHiC, decreasing=T),]
    hp = as_hux(prtab, add_colnames=TRUE)
    position(hp) = "left"
    top_border(hp)[1,]=2
    bottom_border(hp)[1,]=1
    bottom_border(hp)[nrow(hp),]=2
    width(hp) = .75
    quick_pdf(hp, file=paste0(imgpref, "allmetric_auprc_", chrom, ".pdf"))

    roctab = roctab[order(roctab$CHiC, decreasing=T),]
    hp = as_hux(roctab, add_colnames=TRUE)
    position(hp) = "left"
    top_border(hp)[1,]=2
    bottom_border(hp)[1,]=1
    bottom_border(hp)[nrow(hp),]=2
    width(hp) = .75
    quick_pdf(hp, file=paste0(imgpref, "allmetric_auroc_", chrom, ".pdf"))



    # Compare to another validation source (PCHiC)
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
    btovl$gene = tssgr$name[btovl$queryHits] 
    # Map gene onto pchdf:
    btovl$baitID = btgr$name[btovl$subjectHits] 
    pchdf = merge(pchdf, btovl[,c('gene','baitID')])
    # pchdf$oeMid = (pchdf$oeSt + pchdf$oeEnd) / 2

    # Make all possible links:
    geneuq = unique(pchdf[,c('baitChr','gene', 'baitMid')])
    names(geneuq) = c('chr','gene')
    geneuq = merge(geneuq, anno[,c('gene','tss')])
    genegr = GRanges

    window = 1e6
    gqgr = GRanges(seqnames=geneuq$chr, IRanges(geneuq$tss-window, geneuq$tss+window), name=geneuq$gene)
    govl = data.frame(findOverlaps(gqgr, oegr))
    possdf = data.frame(gene=gqgr$name[govl$queryHits], oeID=oegr$name[govl$subjectHits])
    possdf = merge(possdf, oedf)
    red.pchdf = unique(pchdf[,c('gene','oeID')])
    red.pchdf$hit = 1
    possdf = merge(possdf, red.pchdf, all.x=TRUE)
    possdf$hit[is.na(possdf$hit)] = 0

    btype='PCHiC'
    bgr = with(possdf, GRanges(seqnames=oeChr, IRanges(start=oeSt, end=oeEnd), name=gene, score=hit))
    bovl = as.data.frame(findOverlaps(bgr, cgr))
    ind = which(bgr$name[bovl$queryHits] == cgr$name[bovl$subjectHits])
    bovl = bovl[ind,]
    # BG and FG:
    fg = bovl$subjectHits[bgr$score[bovl$queryHits] == 1]
    bg = bovl$subjectHits[bgr$score[bovl$queryHits] == 0]
    print(paste0(length(fg), ', ', length(bg)))

    if (length(fg) > 0){
        curvedf = c()
        for (stype in names(slist)){
            score = slist[[stype]]
            score[is.na(score)] = 0
            roc = roc.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
            pr = pr.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
            print(paste0(stype, ': ', round(roc$auc,3), ", ", round(pr$auc.integral,3)))
            aucdf = rbind(aucdf, data.frame(auroc=round(roc$auc,4), auprc=round(pr$auc.integral,3),stype=stype,valid=btype))
            srdf = data.frame(x=roc$curve[,1], y=roc$curve[,2], stype=stype, curve='ROC')
            spdf = data.frame(x=pr$curve[,1], y=pr$curve[,2], stype=stype, curve='PR')
            curvedf = rbind(curvedf, rbind(srdf, spdf))
        }
    }

    stypes = unique(curvedf$stype)
    scoldf = data.frame(stype=unique(curvedf$stype), 
                        distalone=(stypes %in% c('d1','d7')),
                        ac=(1:length(stypes) %in% grep('ac',stypes)),
                        me=(1:length(stypes) %in% grep('me',stypes)),
                        corr=(1:length(stypes) %in% grep('corr',stypes)),
                        dist=(1:length(stypes) %in% grep('d',stypes)))

    gplot = ggplot(curvedf, aes(x,y,col=stype)) + 
        facet_wrap(~curve) + labs(title=btype) + 
        geom_line() + theme_pubr() + 
        theme(legend.position='right')
    ggsave(paste0(imgpref, btype, '_roc_pr_comp_', chrom, '.png'), gplot, units='in', dpi=450, width=9, height=5)
    glist[[btype]] = gplot


    # -----------------------------------------------------------
    # TEST ALL SAMPLES to check for indexing or validation issue.
    # -----------------------------------------------------------
    # Get strength for this sample:
    # H3K27ac

    btype='CHiC'
    # MERGE WITH VALIDATIONS
    aucdf = c()
    glist = list()
    print(btype)
    bgr = blist[[btype]]
    bovl = as.data.frame(findOverlaps(bgr, cgr))
    ind = which(bgr$name[bovl$queryHits] == cgr$name[bovl$subjectHits])
    bovl = bovl[ind,]
    # BG and FG:
    fg = bovl$subjectHits[bgr$score[bovl$queryHits] == 1]
    bg = bovl$subjectHits[bgr$score[bovl$queryHits] == 0]
    print(paste0(length(fg), ', ', length(bg)))

    for (i in 1:833){
        h5f = H5Fopen(H3K27ac.file)
        h5d = h5f&"matrix"
        ac_mat= h5d[i,]
        H5Dclose(h5d)
        H5Fclose(h5f)
        score = ac_mat[matind[cpairsdf$mind]]
        score[is.na(score)] = 0
        roc = roc.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
        pr = pr.curve(scores.class0 = score[fg], scores.class1 = score[bg], curve=TRUE)
        print(paste0(i, ': ', round(roc$auc,3), ", ", round(pr$auc.integral,3)))
        aucdf = rbind(aucdf, data.frame(auroc=round(roc$auc,4), auprc=round(pr$auc.integral,3),stype=i,valid=btype))
        srdf = data.frame(x=roc$curve[,1], y=roc$curve[,2], stype=i, curve='ROC')
        spdf = data.frame(x=pr$curve[,1], y=pr$curve[,2], stype=i, curve='PR')
        # curvedf = rbind(curvedf, rbind(srdf, spdf))
    }

    # H3K4me1
    h5f = H5Fopen(H3K4me1.file)
    h5d = h5f&"matrix"
    h4samp = h5d[sampind,]
    H5Dclose(h5d)
    H5Fclose(h5f)



}

