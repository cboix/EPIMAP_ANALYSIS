#!/usr/bin/R
# ----------------------------------------------------------------
# Compare the imputed and observed data within enhancers of chr19:
# TODO: compare within enhancers.
# ----------------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
today <- format(Sys.time(), "%m%d%y")

library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(rhdf5)
library(ggpubr)
library(PRROC)

# Arguments:
today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "variance/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'nondhs_recovery_')

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

# ---------------------
# Load expression data:
# ---------------------
H3K27ac.file = 'state_regionlists/H3K27ac_all_bin_on_matched_chromhmm_enh_chr19_r25_e100_matrix.hdf5'
H3K27ac.nam = 'linking_data/H3K27ac_all_bin_dense_on_matched_r25_e100_allchr_merged_names.tsv'

# Load colnames as metadata:
nam.ac = scan(H3K27ac.nam, 'c', sep="\n")
nam.ac = sub(".*name=","", nam.ac)
acdf = data.frame(ind=1:length(nam.ac), full=nam.ac, id=sub("_.*","", nam.ac), type=sub(".*_", "", nam.ac))

# Matrix check: 
h5ls(H3K27ac.file)

# To remove flagged datasets:
abdf = read.delim('Annotation/flagged_potential_abswap.tsv')
smdf = read.delim('Annotation/flagged_potential_sampleswap.tsv')
lwdf = read.delim('Annotation/flagged_low_agreement.tsv')
flagdf = rbind(abdf[,c(1,2)], smdf[,c(1,2)], lwdf)
flagdf$is.flagged = 1

# ------------------------------------------------------------
# Compare recovery of quantiles as a method of estimating TPR:
# ------------------------------------------------------------
nondhs.quantilerec.rda = 'nondhs_quantile_recovery_metrics_validation.Rda'
if (!file.exists(nondhs.quantilerec.rda)){
    h5f = H5Fopen(H3K27ac.file)
    h5d = h5f&"matrix"
    hmat = h5d[]
    H5Dclose(h5d)
    H5Fclose(h5f)
    gc()
    hmat = t(hmat)
    # Separate out obs/imp:
    oind = seq(1, ncol(hmat), by=2)
    iind = seq(2, ncol(hmat), by=2)
    oid = as.character(acdf[oind,'id'])
    omat = hmat[,oind] # Keep raw values for observed
    thresh = apply(omat, 2, function(x){quantile(x, c(.9,.95, .975,.99, .999))})
    imat = 1 * (hmat[, iind] >= 2)
    imat2 = 1 * (hmat[, iind] >= 4)
    imat10 = 1 * (hmat[, iind] >= 10)
    bmat = omat * imat
    bmat2 = omat * imat2
    bmat10 = omat * imat10
    rm(hmat, omat)
    gc()

    # Imputed column sums:
    csi = apply(imat, 2, sum)
    csi2 = apply(imat2, 2, sum)
    csi10 = apply(imat10, 2, sum)
    rm(imat, imat2, imat10)
    # Count number of calls passing each of the observed quantiles: 
    # NOTE: can compare this to closest?? -> precision/recall of quantiles?
    # Makes sense because we used the signal = 2 cutoff for H3K27ac --> gives some estimate of the TPR.
    tdf = c()
    for (i in 1:nrow(thresh)){
        tt = thresh[i,]
        cut = rownames(thresh)[i]
        cat(i, '\t', cut, '\n')
        # Count for each of cutoffs:
        bmat.tt = sweep(bmat, 2, tt, '/')
        bmat.tt2 = sweep(bmat2, 2, tt, '/')
        bmat.tt10 = sweep(bmat10, 2, tt, '/')
        cs = apply(bmat.tt >= 1, 2, sum)
        cs2 = apply(bmat.tt2 >= 1, 2, sum)
        cs10 = apply(bmat.tt10 >= 1, 2, sum)
        df = data.frame(id=oid, cut=cut, thresh=tt, tot=csi, 
                        pass=cs, frac=cs / csi, icut=2)
        df2 = data.frame(id=oid, cut=cut, thresh=tt, tot=csi2, 
                        pass=cs2, frac=cs2 / csi2, icut=4)
        df10 = data.frame(id=oid, cut=cut, thresh=tt, tot=csi10, 
                        pass=cs10, frac=cs10 / csi10, icut=10)
        tdf = rbind(tdf, rbind(rbind(df, df2), df10))
        gc()
    }
    save(tdf, file=nondhs.quantilerec.rda)
} else {
    load(nondhs.quantilerec.rda)
}




# -------------------------------
# Near non-dhs quantile recovery:
# -------------------------------
# Which observed can we match:
rcdir = 'ChromImpute/recovery_metrics/'
rcdf = read.delim(gzfile(paste0(rcdir, 'recovery_metrics_H3K27ac_chr19_stats.tsv.gz')))
# sum(rcdf$near_id %in% acdf$id)
rcdf = rcdf[(rcdf$near_id %in% acdf$id),]
mainid = as.character(rcdf$id)
nearid = as.character(rcdf$near_id)
# sum(mainid == nearid)
oind = seq(1, nrow(acdf), by=2)
oid = as.character(acdf[oind,'id'])
mainind = oind[sapply(mainid, function(x){which(oid == x)})]
nearind = oind[sapply(nearid, function(x){which(oid == x)})]

near.nondhs.quantilerec.rda = 'nondhs_npeak_nearest_quantile_recovery_metrics_validation.Rda'
if (!file.exists(near.nondhs.quantilerec.rda)){
    h5f = H5Fopen(H3K27ac.file)
    h5d = h5f&"matrix"
    hmat = h5d[]
    H5Dclose(h5d)
    H5Fclose(h5f)
    gc()
    hmat = t(hmat)
    # Separate out obs/imp:
    oind = seq(1, ncol(hmat), by=2)
    iind = seq(2, ncol(hmat), by=2)
    oid = as.character(acdf[mainind,'id'])
    omat = hmat[,mainind] # Keep raw values for observed
    thresh = apply(omat, 2, function(x){quantile(x, c(.9,.95, .975,.99, .999))})
    imat = 1 * (hmat[, nearind] >= 2)
    imat2 = 1 * (hmat[, nearind] >= 4)
    imat10 = 1 * (hmat[, nearind] >= 10)
    bmat = omat * imat
    bmat2 = omat * imat2
    bmat10 = omat * imat10
    rm(hmat, omat)
    gc()

    # Imputed column sums:
    csi = apply(imat, 2, sum)
    csi2 = apply(imat2, 2, sum)
    csi10 = apply(imat10, 2, sum)
    rm(imat, imat2, imat10)
    # Count number of calls passing each of the observed quantiles: 
    # NOTE: can compare this to closest?? -> precision/recall of quantiles?
    # Makes sense because we used the signal = 2 cutoff for H3K27ac --> gives some estimate of the TPR.
    near.tdf = c()
    for (i in 1:nrow(thresh)){
        tt = thresh[i,]
        cut = rownames(thresh)[i]
        cat(i, '\t', cut, '\n')
        # Count for each of cutoffs:
        bmat.tt = sweep(bmat, 2, tt, '/')
        bmat.tt2 = sweep(bmat2, 2, tt, '/')
        bmat.tt10 = sweep(bmat10, 2, tt, '/')
        cs = apply(bmat.tt >= 1, 2, sum)
        cs2 = apply(bmat.tt2 >= 1, 2, sum)
        cs10 = apply(bmat.tt10 >= 1, 2, sum)
        df = data.frame(id=oid, cut=cut, thresh=tt, tot=csi, 
                        pass=cs, frac=cs / csi, icut=2)
        df2 = data.frame(id=oid, cut=cut, thresh=tt, tot=csi2, 
                        pass=cs2, frac=cs2 / csi2, icut=4)
        df10 = data.frame(id=oid, cut=cut, thresh=tt, tot=csi10, 
                        pass=cs10, frac=cs10 / csi10, icut=10)
        near.tdf = rbind(near.tdf, rbind(rbind(df, df2), df10))
        gc()
    }
    save(near.tdf, file=near.nondhs.quantilerec.rda)
} else {
    load(near.nondhs.quantilerec.rda)
}



# -----------------------------------------------------------------
# Compare the two marks within the DHS regions only (+ in ENH only)
# -----------------------------------------------------------------
comp.io = function(ox, ix, thresh=2, calc.top1=FALSE){
    # Sort scores for faster roc, etc:
    ord = order(ix)
    ix = ix[ord]
    ox = ox[ord]
    oo = 1 * (ox >= thresh)
    # Perform correlation + prediction comparisons:
    gwcorr = cor(ox, ix)
    cat('CORR:', '\t', round(gwcorr,3), '\t')
    roc = roc.curve(scores.class0=ix[oo == 1],
                    scores.class1=ix[oo == 0], curve=TRUE, sorted=TRUE)
    pr = pr.curve(scores.class0=ix[oo == 1],
                  scores.class1=ix[oo == 0], curve=TRUE, sorted=TRUE)
    cat('AUC:', '\t', round(roc$auc,3), '\t', round(pr$auc.integral,3), '\t')
    # ROC and AUC on the top 1%:
    if (calc.top1){
        threshq = quantile(ox, .99)
        oo2 = 1 * (ox >= threshq)
        roc2 = roc.curve(scores.class0=ix[oo2 == 1],
                         scores.class1=ix[oo2 == 0], curve=TRUE, sorted=TRUE)
        pr2 = pr.curve(scores.class0=ix[oo2 == 1],
                       scores.class1=ix[oo2 == 0], curve=TRUE, sorted=TRUE)
        cat('AUC1:', '\t', round(roc2$auc,3), '\t', round(pr2$auc.integral,3), '\t')
    }
    # Calculate the catch metrics:
    ot1 = 1 * (ox >= quantile(ox, .99))
    ot5 = 1 * (ox >= quantile(ox, .95))
    it1 = 1 * (ix >= quantile(ix, .99))
    it5 = 1 * (ix >= quantile(ix, .95))
    # Catch top 1 % from top 5%
    co1 = sum(ot1 * it5) / sum(ot1)
    ci1 = sum(it1 * ot5) / sum(it1)
    mi1 = sum(it1 * ot1) / sum(it1)
    mo1 = sum(it1 * ot1) / sum(ot1)
    # Aggregate the statistics:
    sdf = data.frame(gwcorr=gwcorr,
                     auroc=roc$auc, auprc=pr$auc.integral,
                     catchobs=co1, catchimp=ci1,
                     matchobs=mo1, matchimp=mi1)
    if (calc.top1){
        sdf$auroc1=roc2$auc
        sdf$auprc1=pr2$auc.integral
    }
    cat("\n")
    return(sdf)
}

nondhs.recovery.rda = 'nondhs_recovery_metrics_validation.Rda'
if (!file.exists(nondhs.recovery.rda)){
    # H3K27ac:
    h5f = H5Fopen(H3K27ac.file)
    h5d = h5f&"matrix"
    hmat = h5d[]
    H5Dclose(h5d)
    H5Fclose(h5f)
    hmat = t(hmat)
    # 
    hids = as.character(unique(acdf$id))
    statsdf = c()
    for (i in 1:length(hids)){
        oind = acdf[acdf$id == hids[i] & acdf$type == 'observed', 'ind']
        iind = acdf[acdf$id == hids[i] & acdf$type == 'imputed', 'ind']
        # Compare obs with imputed
        cat(i,"\t", hids[i],'\t')
        sdf = comp.io(hmat[,oind], hmat[,iind], calc.top1=TRUE)
        sdf$id = hids[i]
        sdf$set = 'Imputed'
        statsdf = rbind(statsdf, sdf)
        if (oind %in% mainind){
            cat(i,"n\t", hids[i],'\t')
            nind = nearind[which(mainind == oind)]
            sdf = comp.io(hmat[,oind], hmat[,nind], calc.top1=TRUE)
            sdf$id = hids[i]
            sdf$set = 'Nearest'
            statsdf = rbind(statsdf, sdf)
        }
    }
    statsdf$mark = 'H3K27ac'
    save(statsdf, file=nondhs.recovery.rda)
} else {
    load(nondhs.recovery.rda)
}




# -----------------------------------------
# Plot these statistics in per-mark manner:
# -----------------------------------------
# Keep only all DHS eval + remove flagged:
sdf = merge(statsdf, flagdf, all.x=TRUE)
sdf = sdf[is.na(sdf$is.flagged),]
sdf$is.flagged = NULL
slong = gather(sdf, metric, value, -id, -set, -mark)

col.paired = brewer.pal(12, 'Paired')
ggplot(slong, aes(metric, value, fill=set)) + 
    geom_boxplot() + 
    scale_fill_manual(values=col.paired[c(1,7)], name='Mark:') + 
    scale_y_continuous(labels=scales::percent) + 
    labs(y='Metric Value', x='Metric (DHSs only)') + 
    theme_pubr()


ggsave(paste0(imgpref, 'metrics_markcomp.png'), dpi=400, units='in', width=6, height=6)
ggsave(paste0(imgpref, 'metrics_markcomp.pdf'), width=6, height=6)


# Look at nearest over auprc:
swide = spread(slong, set, value)
swide = swide[!is.na(swide$Nearest),]

swide$over = swide$Imputed > swide$Nearest
aggregate(over ~ metric, swide, mean)


# Evaluate restriction of the enhancer regions > if we use only H3K4me1

# Remove flagged:
tdf$mark = 'H3K27ac'
tdf$set = 'Imputed'
near.tdf$mark = 'H3K27ac'
near.tdf$set = 'Nearest'
both.tdf = rbind(tdf, near.tdf)

both.tdf = merge(both.tdf, flagdf, all.x=TRUE)
both.tdf = both.tdf[is.na(both.tdf$is.flagged),]
both.tdf$is.flagged = NULL

both.tdf$set.cut = paste0(both.tdf$set, '.', both.tdf$icut)

# col.paired = colorRampPalette(brewer.pal(9, 'Blues'))(100)

col.blues = colorRampPalette(brewer.pal(9, 'Blues'))(100)
col.reds = colorRampPalette(brewer.pal(9, 'Reds'))(100)
rbcols = c(col.blues[c(25,50,75)], col.reds[c(25,50, 75)])
names(rbcols) = c('Imputed.2','Imputed.4','Imputed.10','Nearest.2','Nearest.4','Nearest.10')
both.tdf$set.cut = factor(both.tdf$set.cut, levels=names(rbcols))


gplot = ggplot(both.tdf, aes(cut, frac, fill=set.cut)) + 
    geom_boxplot() + 
    scale_fill_manual(values=rbcols, name='Imputed signal above:') + 
    scale_y_continuous(labels=scales::percent) + 
    # labs(y='Percent of DHS sites with imputed H3K27ac above threshold\nin the top quantile of matched observed data', 
    labs(y='% of DHS sites with high imp. H3K27ac\nin obs. H3K27ac DHS above quantile', 
         x='Observed data quantile (DHSs only)') + 
    theme_pubr()

ggsave(paste0(imgpref, 'pct_in_quantile.png'), gplot, dpi=400, units='in', width=5, height=6)
ggsave(paste0(imgpref, 'pct_in_quantile.pdf'), gplot, width=5, height=6)

ggsave(paste0(imgpref, 'pct_in_quantile_flip.pdf'), gplot + coord_flip(), width=6, height=4.5)


dhs.quantilerec.rda = 'dhs_quantile_recovery_metrics_validation.Rda'
load(dhs.quantilerec.rda)
dhs.tdf = tdf
dhs.tdf$set = paste0('DHS.',dhs.tdf$icut)
nondhs.tdf = both.tdf[both.tdf$set == 'Imputed',]
nondhs.tdf$set = paste0('ENH.',nondhs.tdf$icut)

nd.tdf = rbind(dhs.tdf[,c('id','cut','thresh', 'frac','icut', 'set')],
               nondhs.tdf[,c('id','cut','thresh', 'frac','icut', 'set')])

col.blues = colorRampPalette(brewer.pal(9, 'Blues'))(100)
col.ors = colorRampPalette(brewer.pal(9, 'Oranges'))(100)
rbcols = c(col.blues[c(25,50,75)], col.ors[c(25,50, 75)])
names(rbcols) = c('DHS.2','DHS.4','DHS.10','ENH.2','ENH.4','ENH.10')
rbcols = rbcols[c(1,4,2,5,3,6)]
nd.tdf$set = factor(nd.tdf$set, levels=names(rbcols))
nd.tdf = nd.tdf[nd.tdf$cut != '90%',]

gplot = ggplot(nd.tdf, aes(cut, frac, fill=set)) + 
    geom_boxplot(outlier.cex=.5) + 
    scale_fill_manual(values=rbcols, name='Imputed signal above:') + 
    scale_y_continuous(labels=scales::percent) + 
    # labs(y='Percent of DHS sites with imputed H3K27ac above threshold\nin the top quantile of matched observed data', 
    labs(y='% of DHS or ENH sites with high imp. H3K27ac\nin obs. H3K27ac DHS or ENH above quantile', 
         x='Observed data quantile') + 
    theme_pubr()
ggsave(paste0(imgpref, 'pct_in_quantile_dhsenh_comp.png'), gplot, dpi=400, units='in', width=5, height=6)
ggsave(paste0(imgpref, 'pct_in_quantile_dhsenh_comp.pdf'), gplot, width=5, height=6)
ggsave(paste0(imgpref, 'pct_in_quantile_dhsenh_comp_flip.pdf'), gplot + coord_flip(), width=6, height=4)


# Quantile cutoffs:
gplot = ggplot(nd.tdf, aes(cut, thresh, fill=set)) + 
    geom_boxplot(outlier.size=.5) + 
    scale_fill_manual(values=rbcols, name='Imputed signal above:') + 
    scale_y_continuous(labels=scales::percent) + 
    labs(y='% of DHS or ENH sites with high imp. H3K27ac\nin obs. H3K27ac DHS or ENH above quantile', 
         x='Observed data quantile') + 
    theme_pubr()
gplot + coord_flip()


