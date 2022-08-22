#!/usr/bin/R
# -----------------------------------------------------------
# Compare the imputed and observed data within the DHS sites:
# TODO: compare within enhancers.
# -----------------------------------------------------------
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
imgpref = paste0(imgdir, 'dhs_recovery_')

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
mpref = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_'
enh.file = paste0(mpref, 'matrix.hdf5')
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_matched_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K4me1_all_bin_dense_on_matched_r25_e100_allchr_merged.hdf5'
H3K27ac.nam = 'linking_data/H3K27ac_all_bin_dense_on_matched_r25_e100_allchr_merged.hdf5'
H3K4me1.nam = 'linking_data/H3K4me1_all_bin_dense_on_matched_r25_e100_allchr_merged_names.tsv'

# Load colnames as metadata:
nam.ac = scan('linking_data/H3K27ac_all_bin_dense_on_matched_r25_e100_allchr_merged_names.tsv', 'c', sep="\n")
nam.me = scan('linking_data/H3K4me1_all_bin_dense_on_matched_r25_e100_allchr_merged_names.tsv', 'c', sep="\n")
nam.ac = sub(".*name=","", nam.ac)
nam.me = sub(".*name=","", nam.me)
acdf = data.frame(ind=1:length(nam.ac), full=nam.ac, id=sub("_.*","", nam.ac), type=sub(".*_", "", nam.ac))
medf = data.frame(ind=1:length(nam.me), full=nam.me, id=sub("_.*","", nam.me), type=sub(".*_", "", nam.me))

# Load enhancer indices:
indfile = paste0(mpref, 'masterlist_indices.tsv')
# 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_masterlist_indices.tsv'
matind = as.numeric(scan(indfile, 'c')) + 1
# enhmap = rep(0,nrow(locdf))
# enhmap[matind] = 1:length(matind)

# Matrix check: 
h5ls(enh.file)
h5ls(H3K27ac.file)
h5ls(H3K4me1.file)

# Names for enhancer matrix vs. for other matrices:
enam = scan('Enhancer_matrix_names.txt', "c")
dnam = sort(enam) # Full mark matrices:

abdf = read.delim('Annotation/flagged_potential_abswap.tsv')
smdf = read.delim('Annotation/flagged_potential_sampleswap.tsv')
lwdf = read.delim('Annotation/flagged_low_agreement.tsv')
flagdf = rbind(abdf[,c(1,2)], smdf[,c(1,2)], lwdf)
flagdf$is.flagged = 1

# Remove flagged datasets?
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
    cat('CORR:', '\t', round(gwcorr,3), '\n')
    roc = roc.curve(scores.class0=ix[oo == 1],
                    scores.class1=ix[oo == 0], curve=TRUE, sorted=TRUE)
    pr = pr.curve(scores.class0=ix[oo == 1],
                  scores.class1=ix[oo == 0], curve=TRUE, sorted=TRUE)
    cat('AUC:', '\t', round(roc$auc,3), '\t', round(pr$auc.integral,3), '\n')
    # ROC and AUC on the top 1%:
    if (calc.top1){
    threshq = quantile(ox, .99)
    oo2 = 1 * (ox >= threshq)
    roc2 = roc.curve(scores.class0=ix[oo2 == 1],
                     scores.class1=ix[oo2 == 0], curve=TRUE, sorted=TRUE)
    pr2 = pr.curve(scores.class0=ix[oo2 == 1],
                   scores.class1=ix[oo2 == 0], curve=TRUE, sorted=TRUE)
    cat('AUC1:', '\t', round(roc2$auc,3), '\t', round(pr2$auc.integral,3), '\n')
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
    return(sdf)
}

dhs.recovery.rda = 'dhs_recovery_metrics_validation.Rda'
if (!file.exists(dhs.recovery.rda)){
    # H3K27ac:
    hids = unique(acdf$id)
    statsdf = c()
    for (i in 1:length(hids)){
        cat(i,"\t", hids[i],'\n')
        oind = acdf[acdf$id == hids[i] & acdf$type == 'observed', 'ind']
        iind = acdf[acdf$id == hids[i] & acdf$type == 'imputed', 'ind']
        # Read in the datasets:
        h5f = H5Fopen(H3K27ac.file)
        h5d = h5f&"matrix"
        hmat = h5d[c(oind, iind),]
        H5Dclose(h5d)
        H5Fclose(h5f)
        hmat = t(hmat)
        # Calculate statistics on all the DHS
        sdf = comp.io(hmat[,1], hmat[,2])
        sdf$id = hids[i]
        sdf$set = 'All'
        statsdf = rbind(statsdf, sdf)
        # Calculate statistics on the enhancers only
        sdf = comp.io(hmat[matind,1], hmat[matind,2])
        sdf$id = hids[i]
        sdf$set = 'Enhancers' 
        statsdf = rbind(statsdf, sdf)
    }
    statsdf$mark = 'H3K27ac'

    # H3K4me1:
    hids = unique(medf$id)
    for (i in 1:length(hids)){
        cat(i,"\t", hids[i],'\n')
        oind = medf[medf$id == hids[i] & medf$type == 'observed', 'ind']
        iind = medf[medf$id == hids[i] & medf$type == 'imputed', 'ind']
        # Read in the datasets:
        h5f = H5Fopen(H3K4me1.file)
        h5d = h5f&"matrix"
        hmat = h5d[c(oind, iind),]
        H5Dclose(h5d)
        H5Fclose(h5f)
        hmat = t(hmat)
        # Calculate statistics on all the DHS
        sdf = comp.io(hmat[,1], hmat[,2])
        sdf$id = hids[i]
        sdf$set = 'All'
        sdf$mark = 'H3K4me1'
        statsdf = rbind(statsdf, sdf)
        # Calculate statistics on the enhancers only
        sdf = comp.io(hmat[matind,1], hmat[matind,2])
        sdf$id = hids[i]
        sdf$set = 'Enhancers' 
        sdf$mark = 'H3K4me1'
        statsdf = rbind(statsdf, sdf)
    }
    save(statsdf, file=dhs.recovery.rda)
} else {
    load(dhs.recovery.rda)
}


# -----------------------------------------
# Plot these statistics in per-mark manner:
# -----------------------------------------
# Keep only all DHS eval + remove flagged:
sdf = statsdf[statsdf$set == 'All',]
sdf = merge(sdf, flagdf, all.x=TRUE)
sdf = sdf[is.na(sdf$is.flagged),]
sdf$is.flagged = NULL

swide = gather(sdf, metric, value, -id, -set, -mark)

col.paired = brewer.pal(12, 'Paired')
ggplot(swide, aes(metric, value, fill=mark)) + 
    geom_boxplot() + 
    scale_fill_manual(values=col.paired[c(1,7)], name='Mark:') + 
    scale_y_continuous(labels=scales::percent) + 
    labs(y='Metric Value', x='Metric (DHSs only)') + 
    theme_pubr()
ggsave(paste0(imgpref, 'metrics_markcomp.png'), dpi=400, units='in', width=6, height=6)
ggsave(paste0(imgpref, 'metrics_markcomp.pdf'), width=6, height=6)


# -------------------------------------------------------
# Calculate FDR estimates on all of these simultaneously:
# -------------------------------------------------------
# h5f = H5Fopen(H3K4me1.file)
h5f = H5Fopen(H3K27ac.file)
h5d = h5f&"matrix"
hmat = h5d[]
H5Dclose(h5d)
H5Fclose(h5f)
gc()
hmat = t(hmat)

oind = seq(1, ncol(hmat), by=2)
iind = seq(2, ncol(hmat), by=2)

oid = as.character(medf[oind,'id'])
iid = as.character(medf[iind,'id'])

omat = 1 * (hmat[matind, oind] >= 2)
imat = 1 * (hmat[matind, iind] >= 2)
bmat = omat * imat

# Best marg:
amat = t(omat) %*% omat
amat2 = sweep(amat, 2, diag(amat), '/')
diag(amat2) = 0

# Shows a super high range:
mdf = data.frame(id=medf$id[oind], mark='H3K27ac')
mdf$obs = colSums(omat)
mdf$imp = colSums(imat)
mdf$both = colSums(bmat)
mdf$rec = mdf$both / mdf$imp
mdf$near = apply(amat2, 1, max, na.rm=T)

# Calculate F1 instead:
mdf$precision = mdf$both / mdf$imp
mdf$recall = mdf$both / mdf$obs
mdf$f1 = 2 * (mdf$precision * mdf$recall) / (mdf$precision + mdf$recall)

pmat = sweep(amat, 1, diag(amat), '/')
rmat = sweep(amat, 2, diag(amat), '/')
fmat = 2 * (pmat * rmat) / (pmat + rmat)
diag(fmat) = 0
mdf$near_f1 = apply(fmat,1, max)
mean(mdf$f1 > mdf$near_f1) # Not good, only 50%, 58% when above 2
mean(mdf$f1 > mdf$near_f1) # Not good, only 50%, 58% when above 2

cind = which(mdf$obs < 20000) # For H3K4me1, there are problematic tracks.

# Probably not actually comparable:
mean(mdf$rec > mdf$near) # Highly variable depending on what cutoffs we use...
mean((mdf$rec > mdf$near)[-cind]) 
summary(mdf$rec)
summary(mdf$rec[-cind])
# sort(mdf$rec[-cind])


# Compare 

# Evaluate: AUPRC and AUROC, CATCHonOBS vs CATCHonIMP for each (not top %)


# -------------------------------------------------------
# Estimate recovery and false positive rate in enhancers:
# -------------------------------------------------------
# Reduce to only datasets where BOTH marks present:
# Evaluate agreement of H3K4me1 * H3K27ac IN the enhancer regions


# Evaluate restriction of the enhancer regions > if we use only H3K4me1

# ------------------------------------------------------------
# Compare recovery of quantiles as a method of estimating TPR:
# ------------------------------------------------------------
dhs.quantilerec.rda = 'dhs_quantile_recovery_metrics_validation.Rda'
if (!file.exists(dhs.quantilerec.rda)){
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
    omat = hmat[matind,oind] # Keep raw values for observed
    thresh = apply(omat, 2, function(x){quantile(x, c(.9,.95, .975,.99, .999))})
    imat = 1 * (hmat[matind, iind] >= 2)
    imat2 = 1 * (hmat[matind, iind] >= 4)
    imat10 = 1 * (hmat[matind, iind] >= 10)
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
    save(tdf, file=dhs.quantilerec.rda)
} else {
    load(dhs.quantilerec.rda)
}

# Remove flagged:
tdf$mark = 'H3K27ac'
tdf = merge(tdf, flagdf, all.x=TRUE)
tdf = tdf[is.na(tdf$is.flagged),]
tdf$is.flagged = NULL

col.paired = colorRampPalette(brewer.pal(9, 'Blues'))(100)

gplot = ggplot(tdf, aes(cut, frac, fill=factor(icut))) + 
    geom_boxplot() + 
    scale_fill_manual(values=col.paired[c(25,50,75)], name='Imputed signal above:') + 
    scale_y_continuous(labels=scales::percent) + 
    # labs(y='Percent of DHS sites with imputed H3K27ac above threshold\nin the top quantile of matched observed data', 
    labs(y='% of DHS sites with high imp. H3K27ac\nin obs. H3K27ac DHS above quantile', 
         x='Observed data quantile (DHSs only)') + 
    theme_pubr()
ggsave(paste0(imgpref, 'pct_in_quantile.png'), gplot, dpi=400, units='in', width=5, height=6)
ggsave(paste0(imgpref, 'pct_in_quantile.pdf'), gplot, width=5, height=6)

ggsave(paste0(imgpref, 'pct_in_quantile_flip.pdf'), gplot + coord_flip(), width=6, height=4.5)

# Which is the really bad sample:
tdf[tdf$frac < .5 & tdf$icut == 10 & tdf$cut=='90%',]
fid = 'BSS01366' # Neural deriv doesn't really match
tdf[tdf$id == fid,]

tstats = spread(aggregate(frac ~ cut + icut, tdf, mean), icut, frac)
# Output:
#     cut          2          4        10
# 1   90% 0.69314722 0.84991336 0.9455563
# 2   95% 0.51365005 0.74208339 0.9118000
# 3 97.5% 0.32513438 0.56908472 0.8437754
# 4   99% 0.15147658 0.32624866 0.6692870
# 5 99.9% 0.01678249 0.04590063 0.1882259
# On avg, 30-50% not in obs with cutoff of 2, goes up to 15-25% and then to 5-10% (depending on cutoff).

# -----------------------------------------------------------------
# Compare recovery of quantiles conditioned on the npeak statistic:
# -----------------------------------------------------------------
npeak.quantilerec.rda = 'dhs_npeak_quantile_recovery_metrics_validation.Rda'
if (!file.exists(npeak.quantilerec.rda)){
    # H3K27ac matrix:
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
    omat = hmat[matind,oind] # Keep raw values for observed
    thresh = apply(omat, 2, function(x){quantile(x, c(.9,.95, .975,.99, .999))})
    imat = 1 * (hmat[matind, iind] >= 2)
    imat2 = 1 * (hmat[matind, iind] >= 4)
    imat10 = 1 * (hmat[matind, iind] >= 10)
    bmat = omat * imat
    bmat2 = omat * imat2
    bmat10 = omat * imat10
    rm(hmat, omat)
    gc()

    # Peak definitions:
    lddir = 'linking_data/'
    suff = "_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged_stats_df.tsv.gz"
    markdf = read.delim(gzfile(paste0(lddir, 'H3K27ac', suff)), header=T)
    # Transform matrix for the DHS npeak stats:
    npeak = markdf$nsamp[matind]
    uvals = sort(unique(npeak))
    npdf = aggregate(Mean ~ nsamp, markdf[matind,], length)
    names(npdf)[2] = 'nloc'
    tform = make.tform(x=npeak, u=uvals, norm=FALSE)
    tform = t(tform)

    # Imputed column sums:
    csi = tform %*% imat
    csi2 = tform %*% imat2
    csi10 = tform %*% imat10
    rm(imat, imat2, imat10)
    # Count avg frac of calls passing each of the observed quantiles per nsamp cutoff
    # Gives some estimate of the TPR.
    # NOTE: can compare this to closest?? -> precision/recall of quantiles?
    np.tdf = c()
    for (i in 1:nrow(thresh)){
        tt = thresh[i,]
        cut = rownames(thresh)[i]
        cat(i, '\t', cut, '\n')
        # Count for each of cutoffs:
        bmat.tt = sweep(bmat, 2, tt, '/')
        bmat.tt2 = sweep(bmat2, 2, tt, '/')
        bmat.tt10 = sweep(bmat10, 2, tt, '/')
        # Get the intersections per DHS type:
        cs = tform %*% (bmat.tt >= 1)
        cs2 = tform %*% (bmat.tt2 >= 1)
        cs10 = tform %*% (bmat.tt10 >= 1)
        # Margins for plotting:
        df = rbind(data.frame(nsamp=uvals, frac=apply(cs/csi, 1, mean, na.rm=T), frac.sd=apply(cs/csi, 1, sd, na.rm=T), icut=2, cut=cut),
                   data.frame(nsamp=uvals, frac=apply(cs2/csi2, 1, mean, na.rm=T), frac.sd=apply(cs2/csi2, 1, sd, na.rm=T), icut=4, cut=cut),
                   data.frame(nsamp=uvals, frac=apply(cs10/csi10, 1, mean, na.rm=T), frac.sd=apply(cs10/csi10, 1, sd, na.rm=T), icut=10, cut=cut))
        df = merge(df, npdf, all.x=TRUE)
        np.tdf = rbind(np.tdf, df)
        gc()
    }
    save(np.tdf, file=npeak.quantilerec.rda)
} else {
    load(npeak.quantilerec.rda)
}


# Remove flagged:
col.blues = colorRampPalette(brewer.pal(9, 'Blues'))(100)

gplot = ggplot(np.tdf, aes(nsamp, frac, color=factor(icut), fill=factor(icut))) + 
    facet_wrap(. ~ cut, nrow=1) + 
    geom_ribbon(aes(ymin=frac - 2 * frac.sd / sqrt(227), ymax = frac + 2 * frac.sd / sqrt(227)), alpha=0.5, color=NA) + 
    geom_line() + 
    scale_color_manual(values=col.blues[c(25,50,75)], name='Imputed signal above:') + 
    scale_fill_manual(values=col.blues[c(25,50,75)], name='Imputed signal above:') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    labs(y='% of DHS sites with high imp. H3K27ac\nin obs. H3K27ac DHS above quantile', 
         x='Number of samples for which the DHS is labeled as an active enhancer') + 
    theme_pubr()
ggsave(paste0(imgpref, 'pct_in_quantile_npeak.png'), gplot, dpi=400, units='in', width=9, height=3.5)
ggsave(paste0(imgpref, 'pct_in_quantile_npeak.pdf'), gplot, width=9, height=3.5)


ggsave(paste0(imgpref, 'pct_in_quantile_npeak_xlog10.png'), gplot + scale_x_log10(), dpi=400, units='in', width=9, height=3.5)
ggsave(paste0(imgpref, 'pct_in_quantile_npeak_xlog10.pdf'), gplot + scale_x_log10(), width=9, height=3.5)



# ---------------------------------------------------------------------
# Compare recovery of quantiles conditioned on the npeak statistic:
# For the NEAREST sample rather than the appropriately matched imputed.
# ---------------------------------------------------------------------
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

near.quantilerec.rda = 'dhs_npeak_nearest_quantile_recovery_metrics_validation.Rda'
if (!file.exists(near.quantilerec.rda)){
    # H3K27ac matrix:
    h5f = H5Fopen(H3K27ac.file)
    h5d = h5f&"matrix"
    hmat = h5d[]
    H5Dclose(h5d)
    H5Fclose(h5f)
    gc()
    hmat = t(hmat)

    # Separate out obs/imp:
    omat = hmat[matind,mainind] # Keep raw values for observed
    thresh = apply(omat, 2, function(x){quantile(x, c(.9,.95, .975,.99, .999))})
    # Treat the nearest as imputed:
    imat = 1 * (hmat[matind, nearind] >= 2)
    imat2 = 1 * (hmat[matind, nearind] >= 4)
    imat10 = 1 * (hmat[matind, nearind] >= 10)
    bmat = omat * imat
    bmat2 = omat * imat2
    bmat10 = omat * imat10
    rm(hmat, omat)
    gc()

    # Peak definitions:
    lddir = 'linking_data/'
    suff = "_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged_stats_df.tsv.gz"
    markdf = read.delim(gzfile(paste0(lddir, 'H3K27ac', suff)), header=T)
    # Transform matrix for the DHS npeak stats:
    npeak = markdf$nsamp[matind]
    uvals = sort(unique(npeak))
    npdf = aggregate(Mean ~ nsamp, markdf[matind,], length)
    names(npdf)[2] = 'nloc'
    tform = make.tform(x=npeak, u=uvals, norm=FALSE)
    tform = t(tform)

    # Imputed column sums:
    csi = tform %*% imat
    csi2 = tform %*% imat2
    csi10 = tform %*% imat10
    rm(imat, imat2, imat10)
    near.tdf = c()
    for (i in 1:nrow(thresh)){
        tt = thresh[i,]
        cut = rownames(thresh)[i]
        cat(i, '\t', cut, '\n')
        # Count for each of cutoffs:
        bmat.tt = sweep(bmat, 2, tt, '/')
        bmat.tt2 = sweep(bmat2, 2, tt, '/')
        bmat.tt10 = sweep(bmat10, 2, tt, '/')
        # Get the intersections per DHS type:
        cs = tform %*% (bmat.tt >= 1)
        cs2 = tform %*% (bmat.tt2 >= 1)
        cs10 = tform %*% (bmat.tt10 >= 1)
        # Margins for plotting:
        df = rbind(data.frame(nsamp=uvals, frac=apply(cs/csi, 1, mean, na.rm=T), frac.sd=apply(cs/csi, 1, sd, na.rm=T), icut=2, cut=cut),
                   data.frame(nsamp=uvals, frac=apply(cs2/csi2, 1, mean, na.rm=T), frac.sd=apply(cs2/csi2, 1, sd, na.rm=T), icut=4, cut=cut),
                   data.frame(nsamp=uvals, frac=apply(cs10/csi10, 1, mean, na.rm=T), frac.sd=apply(cs10/csi10, 1, sd, na.rm=T), icut=10, cut=cut))
        df = merge(df, npdf, all.x=TRUE)
        near.tdf = rbind(near.tdf, df)
        gc()
    }
    save(near.tdf, file=near.quantilerec.rda)
} else {
    load(near.quantilerec.rda)
}


# TODO: Plot this
col.reds = colorRampPalette(brewer.pal(9, 'Reds'))(100)

gplot = ggplot(near.tdf, aes(nsamp, frac, color=factor(icut), fill=factor(icut))) + 
    facet_wrap(. ~ cut, nrow=1) + 
    geom_ribbon(aes(ymin=frac - 2 * frac.sd / sqrt(227), ymax = frac + 2 * frac.sd / sqrt(227)), alpha=0.5, color=NA) + 
    geom_line() + 
    scale_color_manual(values=col.reds[c(25,50,75)], name='Nearest signal above:') + 
    scale_fill_manual(values=col.reds[c(25,50,75)], name='Nearest signal above:') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    labs(y='% of DHS sites with high imp. H3K27ac\nin obs. H3K27ac DHS above quantile', 
         x='Number of samples for which the DHS is labeled as an active enhancer') + 
    theme_pubr()
ggsave(paste0(imgpref, 'pct_in_quantile_npeak_near.png'), gplot, dpi=400, units='in', width=9, height=3.5)
ggsave(paste0(imgpref, 'pct_in_quantile_npeak_near.pdf'), gplot, width=9, height=3.5)


np.tdf$set = paste0('I', np.tdf$icut)
near.tdf$set = paste0('N', near.tdf$icut)
both.tdf = rbind(near.tdf, np.tdf)
rbcols = c(col.blues[c(25,50,75)], col.reds[c(25,50, 75)])
names(rbcols) = c('I2','I4','I10','N2','N4','N10')
both.tdf$set = factor(both.tdf$set, levels=c('I2','I4','I10','N2','N4','N10'))


gplot = ggplot(both.tdf, aes(nsamp, frac, color=set, fill=set)) + 
    facet_wrap(. ~ cut, nrow=1) + 
    geom_ribbon(aes(ymin=frac - 2 * frac.sd / sqrt(227), ymax = frac + 2 * frac.sd / sqrt(227)), alpha=0.5, color=NA) + 
    geom_line() + 
    scale_color_manual(values=rbcols, name='Signal above:') + 
    scale_fill_manual(values=rbcols, name='Signal above:') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    labs(y='% of DHS sites with high imp. H3K27ac\nin obs. H3K27ac DHS above quantile', 
         x='Number of samples for which the DHS is labeled as an active enhancer') + 
    theme_pubr()
ggsave(paste0(imgpref, 'pct_in_quantile_npeak_both.png'), gplot, dpi=400, units='in', width=9, height=3.5)
ggsave(paste0(imgpref, 'pct_in_quantile_npeak_both.pdf'), gplot, width=9, height=3.5)

ggsave(paste0(imgpref, 'pct_in_quantile_npeak_both_xlog10.png'), gplot + scale_x_log10(), dpi=400, units='in', width=9, height=3.5)
ggsave(paste0(imgpref, 'pct_in_quantile_npeak_both_xlog10.pdf'), gplot + scale_x_log10(), width=9, height=3.5)


mean(meta[nearid,'infoline'] == meta[mainid,'infoline'])
mean(meta[nearid,'infoline'] == meta[mainid,'infoline'])


# Average:
enhtot = sum(unique(np.tdf[,c('nsamp', 'nloc')])$nloc)
both.tdf$wfrac = both.tdf$frac * both.tdf$nloc / enhtot
aggregate(wfrac ~ set + cut, both.tdf, sum)

aggregate(np.tdf$set





# -------------------------------------------------------
# Calculate FDR estimates on both marks together:
# -------------------------------------------------------
mids = unique(as.character(medf$id))
cids = unique(as.character(acdf$id))
bids = mids[mids %in% cids] # Overlapping datasets:

mind = c(sapply(bids, function(x){medf$ind[medf$id == x]}))
cind = c(sapply(bids, function(x){acdf$ind[acdf$id == x]}))

# Load in both matrices:
h5f = H5Fopen(H3K4me1.file)
h5d = h5f&"matrix"
mmat = h5d[mind,]
H5Dclose(h5d)
H5Fclose(h5f)
gc()
mmat = t(mmat)
# H3K27ac:
h5f = H5Fopen(H3K27ac.file)
h5d = h5f&"matrix"
cmat = h5d[cind,]
H5Dclose(h5d)
H5Fclose(h5f)
gc()
cmat = t(cmat)

# --------------------------------------
# Calculations of f1, recovery, and fdr:
# --------------------------------------
oind = seq(1, ncol(mmat), by=2)
iind = seq(2, ncol(mmat), by=2)
omat = 1 * (mmat[matind, oind] >= 2) * (cmat[matind, oind] >= 2)
imat = 1 * (mmat[matind, iind] >= 2) * (cmat[matind, iind] >= 2)
bmat = omat * imat

# Try to match the best marginal:
amat = t(omat) %*% omat
aind = which(diag(amat) < 200) # 6 samp
amat2 = sweep(amat, 2, diag(amat), '/')
diag(amat2) = 0
amat2[aind,] = 0
amat2[,aind] = 0

# Shows a super high range:
mdf = data.frame(id=medf$id[mind[oind]], mark='H3K27ac')
mdf$obs = colSums(omat)
mdf$imp = colSums(imat)
mdf$both = colSums(bmat)
mdf$rec = mdf$both / mdf$imp
mdf$near = apply(amat2, 1, max, na.rm=T)
mean((mdf$rec > mdf$near)[-aind]) # Precision super high for near if very small datasets.

# Calculate F1 as a balanced score:
mdf$precision = mdf$both / mdf$imp
mdf$recall = mdf$both / mdf$obs
mdf$f1 = 2 * (mdf$precision * mdf$recall) / (mdf$precision + mdf$recall)
# Get precision across nearest:
rmat = sweep(amat, 1, diag(amat), '/') # max over rows is both / obs (recall
pmat = sweep(amat, 2, diag(amat), '/')
fmat = 2 * (pmat * rmat) / (pmat + rmat)
diag(fmat) = 0
# Remove problematic locations (17 samples)
fmat[aind,] = 0
fmat[,aind] = 0
mdf$near_f1 = apply(fmat,1, max)
mean((mdf$f1 > mdf$near_f1)[-aind]) # 88.4% of the time f1 is higher 
mean((mdf$f1 / mdf$near_f1)[-aind]) # 24.7% relative increase
mean((mdf$f1 - mdf$near_f1)[-aind]) # 0.086 flat increase


# --------------------------------------------------------------------
# Compare relative quantiles, because it is hard to evaluate the rest.
# --------------------------------------------------------------------
# TODO: Load the enhancer matrix --> How many of the "Active enh are FP?"

# Q: how many of calls we make of >= 2 are in the top 1% or top 5% - because not well calibrated.
oind = seq(1, ncol(mmat), by=2)
iind = seq(2, ncol(mmat), by=2)
omat = 1 * (mmat[matind, oind] >= 2) * (cmat[matind, oind] >= 2)
imat = 1 * (mmat[matind, iind] >= 2) * (cmat[matind, iind] >= 2)
bmat = omat * imat



# --------------------------------------
# Flag datasets by odds ratio + overlap:
# --------------------------------------
# BSS00095 looks odd (DNase-seq is bad - plot all but 4 flagged DNase-seq samples.)
# TODO: Make sure the DNase-seq is flagged.
# Also look at BSS005XX (MPP)


