#!/usr/bin/R
# ---------------------------------------
# Plot the modules validation
# Make correlation based FDR predictions:
# ------------------------------------------------
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
# library(GenomicRanges)
# library(Matrix)
# library(scales)
# library(plotrix)
library(PRROC)
# library(huxtable)
# library(DescTools) # For Gini

# Arguments:
today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "variance/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'dnase_enh_')

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
enh.file = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_matrix.hdf5'
dnase.file = 'linking_data/DNase-seq_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K4me1_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'

h5ls(enh.file)
h5ls(dnase.file)

# Names for enhancer matrix vs. for other matrices:
enam = scan('Enhancer_matrix_names.txt', "c")
dnam = sort(enam)

# --------------------------
# Get statistics on overlap:
# --------------------------
dnovl.rda = 'DNase-seq_enhancer_overlap_statistics.Rda'
if (!file.exists(dnovl.rda)){
    statsdf = c()
    chunksize = 10
    cind = c()
    for (i in 1:length(dnam)){
        sample = dnam[i]
        j = which(enam == sample)
        cat(i,'\t', sample,'\t')
        # DNase-seq - read as chunks:
        if (!(i %in% cind)){
            cind = i:(i + chunksize - 1)
            h5f = H5Fopen(dnase.file)
            h5d = h5f&"matrix"
            dmat = h5d[cind,]
            H5Dclose(h5d)
            H5Fclose(h5f)
        }
        dX = dmat[cind == i,]
        # Enhancers:
        h5f = H5Fopen(enh.file)
        h5d = h5f&"matrix"
        eX = h5d[j,]
        H5Dclose(h5d)
        H5Fclose(h5f)
        # Calculate the roc:
        roc = roc.curve(scores.class0=dX[eX == 1],
                        scores.class1=dX[eX == 0], curve=TRUE)
        pr = pr.curve(scores.class0=dX[eX == 1],
                      scores.class1=dX[eX == 0], curve=TRUE)
        cat(round(roc$auc,3), '\t', round(pr$auc.integral,3), '\t')
        # Margin calculations:
        sdf = data.frame(i=i, id=sample,
                         auroc=roc$auc, auprc=pr$auc.integral,
                         tot.enh=sum(eX == 1), tot.pk=sum(dX >= 2), tot.high.pk=sum(dX >= 5),
                         pct.inenh=mean(dX[eX == 1] >= 2), pct.outenh=mean(dX[eX == 0] >= 2),
                         pct.high.inenh=mean(dX[eX == 1] >= 5), pct.high.outenh=mean(dX[eX == 0] >= 5),
                         mean.inenh = mean(dX[eX == 1]), mean.outenh=mean(dX[eX == 0]),
                         pct.inpk = mean(eX[dX >= 2]), pct.outpk = mean(eX[dX < 5]))
        # Enrichment fraction:
        mat = matrix(c(sum(dX[eX == 1] >= 2), sdf$tot.enh, sdf$tot.pk, length(eX)), nrow=2, byrow=T)
        fit = fisher.test(mat)
        sdf$or.est = fit$estimate # Odds ratio estimate
        cat(round(fit$estimate,2), '\n')
        statsdf = rbind(statsdf, sdf)
    }
    save(statsdf, file=dnovl.rda)
} else { 
    load(dnovl.rda)
}

flagdn.df = statsdf[statsdf$tot.pk >= 3e6,]
statsdf = statsdf[statsdf$tot.pk < 3e6,]

# ---------------------------
# Plot statistics on overlap:
# ---------------------------
g1 = ggplot(statsdf, aes(tot.enh, pct.inenh)) + 
    geom_point() + 
    labs(x='Total # Enhancers', y='% Enhancers with DNase-seq peak') +
    theme_pubr()

g2 = ggplot(statsdf, aes(tot.enh, auroc)) + 
    geom_point() + 
    labs(x='Total # Enhancers', y='AUROC of predicting ENH from DNase alone') +
    theme_pubr()

g3 = ggplot(statsdf, aes(pct.inpk, pct.inenh)) + 
    geom_point() + 
    labs(y='% Enh. with DNase-seq peak', x='% DNase-seq peaks that are Enh.') +
    theme_pubr()

g4 = ggplot(statsdf, aes(tot.enh, tot.enh * pct.inenh)) + 
    geom_point() + 
    labs(x='Total # Enhancers', y='Total # Enhancers with DNase-seq peaks') +
    theme_pubr()

garr = ggarrange(g1,g2, g3, g4)
ggsave(paste0(imgpref, 'stats_on_int_enh.png'), garr, dpi=400, units='in', width=9, height=9)

g5 = ggplot(statsdf, aes(pct.inenh, or.est)) + 
    geom_point() + 
    theme_pubr()

png(paste0(imgpref, 'stats_on_int_enh_hists.png'), res=400, units='in', width=9, height=4)
layout(matrix(1:2, ncol=2))
par(mar=c(4,4,.25,.25))
hist(statsdf$auroc, col='grey75', border='white', main='', xlab='Percent of enhancers with DNase-seq signal > 2')
par(mar=c(4,4,.25,.25))
hist(statsdf$pct.inenh, col='grey75', border='white', main='', xlab='AUROC for predicting enhancers with DNase-seq')
dev.off()


hist(statsdf$or.est, col='grey75', border='white')

hist(statsdf$mean.inenh, col='grey75', border='white')
hist(statsdf$mean.outenh, col='grey75', border='white')

hist(statsdf$pct.outenh, col='grey75', border='white')

ggplot(statsdf, aes(auroc, auprc)) + 
    geom_point() + 
    theme_pubr()




hist(statsdf$auroc, col)
hist(statsdf$auprc)


# --------------------------------------
# Flag datasets by odds ratio + overlap:
# --------------------------------------
# BSS00095 looks odd (DNase-seq is bad - plot all but 4 flagged DNase-seq samples.)
# TODO: Make sure the DNase-seq is flagged.
# Also look at BSS005XX (MPP)


