#!/usr/bin/R
# ----------------------------------------------------
# Plot chromImpute metrics for imputing existing files
# Metrics (from ChromImpute manual @ Jason Ernst):
# (1) the fraction of the observed top percent1 locations in the imputed top percent1 locations
# (2) the fraction of the imputed top percent1 in the observed top percent2
# (3) the fraction of the observed top percent1 in the imputed top percent2
# (4) the correlation between the observed and imputed data
# (5) the area under the ROC for predicting the top percent1 imputed signal with the full range of observed signal
# (6) the area under the ROC for predicting the top percent1 observed signal with the full range of imputed signal
# ----------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(ggrepel)
library(ggpubr)
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    # recent = '20190107-0201'
    # recent = '20190224-1305'
    recent = '20190905-1037'
    qcfile=paste0('all_qc_',recent,'.tsv')
    evalfile=paste0('all_eval_',recent,'.tsv')
    cat(paste0("[STATUS] No arguments supplied: Need evalfile and qcfile filenames.\n", 
                 "[STATUS] DEFAULTING TO RUNS FROM ", recent, "\n"))
} else {        
    evalfile = args[1]
    qcfile = args[2]
}

# Directories:
imgpref = paste0(img, "metrics/")
system(paste('mkdir -p', imgpref))

# Read in metrics and get ranks:
metricsdf = read.delim(evalfile)
names(metricsdf)[1] = 'id'  # BSSID
metrics = colnames(metricsdf)[!(colnames(metricsdf) %in% c('id','mark'))]
ranks = 0
for (metric in metrics){
    ranks = ranks + rank(metricsdf[[metric]])
}
metricsdf$avg.rank = ranks / length(metricsdf)
metricsdf$cor.rank = rank(metricsdf$Correlation)
pm = c('ATAC-seq','DNase-seq','H2AFZ','H3K27ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac')
bm = c('H3K27me3','H3K36me3','H3K79me2','H3K9me3','H4K20me1')
fc = c('CTCF','EP300','POLR2A','SMC3','RAD21')
markdef = rbind(data.frame(mark=pm, type='punctate'), 
                data.frame(mark=bm, type='broad'), 
                data.frame(mark=fc, type='factor'))
metricsdf = merge(metricsdf, markdef)
aggregate(mark ~ type, metricsdf, length)

# write.table(metricsdf[,c('id','mark')], file='impobs_comb_tocompare.tsv', 
#             quote=F, sep="\t", row.names=F, col.names=F)

# View the metrics as a heatmaps:
logit = function(x, scale=1){ x = x / scale; log(x / (1 - x)) }
sigm = function(x, scale=1){ scale / (1 + exp(-x)) }

# Load QC file
qcdf = read.delim(qcfile, header=F)
names(qcdf) = c('file','nreads','estFragLen', 'corr_estFragLen','PhantomPeak',
                'corr_phantomPeak','argmin_corr','min_corr','NSC','RSC','QualityTag')
mat = t(sapply(qcdf$file, function(x){
                   x = sub(".sub.*","",sub("FINAL_","",x)); 
                   strsplit(x,"_")[[1]]}))
qcdf$id = mat[,2]
qcdf$mark = mat[,1]

# ------------------------------------------------------
# 2. Color according to elbow changes in the correlation
# GOAL: Flag files for review
# ------------------------------------------------------
metricsdf$status = 'MEDIUM'
metricsdf = metricsdf[order(metricsdf$cor.rank, decreasing=T),]
for(mark in c(pm,bm)){
    idx = which(metricsdf$mark == mark)
    subdf = metricsdf[idx,]
    jumps = c(diff(subdf$Correlation),0) / subdf$Correlation
    if (mark %in% c('ATAC-seq', 'H4K20me1', 'H2AFZ')){
        cut = head(which(jumps < -0.07),1)
    } else {
        cut = head(which(jumps < -0.05),1)
    }
    status = c(rep(0, cut), rep(1, nrow(subdf) - cut))
    metricsdf$status[idx[(cut + 1):length(idx)]] <- 'LOW'
    metricsdf$status[idx[1:5]] <- 'TOP'
} 
ovals = c('TOP' = 'forestgreen','MEDIUM' = 'lightgrey','LOW' = 'indianred')
vals = c('TOP' = 'lightgrey','MEDIUM' = 'lightgrey','LOW' = 'indianred')
mapmetric = data.frame(metric=c("BOTH_1.0", "IMPUTE_1.0_OBSERVED_5.0", "OBSERVED_1.0_IMPUTE_5.0","Correlation", "IMPUTE_1.0_AUC_PREDICT_OBSERVE","OBSERVED_1.0_AUC_PREDICT_IMPUTE", 'nreads','RSC','NSC'), 
                       short.metric=c('Match1','Catch1imp','Catch1obs', 'GWCorr', 'AucImp1','AucObs1','# Reads','RSC','NSC'))

# Write out the low metric files:
lowmdf = metricsdf[metricsdf$status == 'LOW', c('id','mark', 'Correlation', 'cor.rank')]
write.table(lowmdf, file='low_qc_metrics_files.tsv', row.names=F, quote=F, sep="\t")
qcwide = merge(qcdf, metricsdf)

# -----------------------------------------------------------------------------------------
# 1b. For the correlation metric, model the maximal GWCorr of observed tracks vs. imputed. 
# See which have better correlates.
# Also add QC
# -----------------------------------------------------------------------------------------
# Load matrices
source(paste0(bindir, 'load_distance_matrices.R'))

# -----------------------------
# Add in the sample comparisons
# -----------------------------
mdf = c()
for (mark in names(full.ll)){
    mat = full.ll[[mark]]
    oidx = grep("obs", colnames(mat))
    iidx = grep("imp", colnames(mat))
    mx = 1 - apply(mat[oidx, iidx], 1, min)
    wmx = apply(mat[oidx, iidx], 1, which.min)
    ct = sapply(names(mx), function(x){strsplit(x,"_")[[1]][1]})
    closeto = sapply(colnames(mat[,iidx])[wmx], function(x){strsplit(x, "_")[[1]][1]})
    ct = sapply(names(mx), function(x){strsplit(x,"_")[[1]][1]})
    mdf = rbind(mdf, data.frame(mark=mark, id=ct, maxCorr=mx, msamp=closeto))
}
rownames(mdf) = NULL
qcwide = merge(qcwide, mdf, all.x=TRUE)

# ---------------------------
# Add in the mark comparisons
# ---------------------------
markdf = read.delim(gzfile('mark_withinsamp_dist.tsv.gz'), header=F, sep="\t")
names(markdf) = c('id','mark','against','dist')
markdf = markdf[grep('imp',markdf$against),]
wmark = spread(markdf, against, dist)

# Track where it came from
wsdf = aggregate(dist ~ mark + id, markdf, min)
wsmdf = aggregate(against ~ mark + id, markdf, function(x){strsplit(as.character(x[1]), "_")[[1]][1]})
wsdf = merge(wsdf, wsmdf)

qcwide = merge(qcwide, wsdf, all.x=TRUE)
qcwide$dist[is.na(qcwide$dist)] = 1 - qcwide$Correlation[is.na(qcwide$dist)]

cors = cbind(qcwide$maxCorr, 1 - qcwide$dist, NA)
assign = cbind(as.character(qcwide$msamp), as.character(qcwide$against), NA)
lid = apply(cors, 1, function(x){ id = which.max(x); if(length(id) == 0){id = 3}; id})

qcwide$maxCorr2 = cors[cbind(1:length(lid),lid)]
qcwide$assign = assign[cbind(1:length(lid),lid)]
qcwide$closeby = factor(lid)
qcwide$QCpass = factor((qcwide$RSC > 1) + (qcwide$NSC > 1.05))
qcwide$difMax = qcwide$maxCorr2 - qcwide$Correlation

# NOTE: Can't order just by id, must order by both 
qcwide = qcwide[order(qcwide$cor.rank, decreasing=TRUE),]
qcwide$cm = paste0(qcwide$id, "_", qcwide$mark)
qcwide$cm = factor(qcwide$cm, levels=unique(qcwide$cm))

filter(qcwide, mark == 'ATAC-seq' & id == 'BSS01543')
aggregate(status ~ mark + id, qcwide, length)
dim(unique(qcwide[,c('mark','id')]))

ggplot(qcwide, aes(cm, Correlation, fill=status)) +
    facet_wrap(~ mark, scales='free', nrow=3) + 
    labs(x=paste0('Observed Sample + Mark'), y='Correlation') +
    geom_bar(stat='identity', width=1) +
    scale_fill_manual(values=vals) +
    ylim(0,1) +
    theme_minimal() + 
    theme(legend.position='none') + 
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'corr_all18assays.pdf'), dpi=300, width=8, height=4,units='in')
ggsave(paste0(imgpref, 'corr_all18assays.png'), dpi=300, width=8, height=4,units='in')


ocols = c("id", "mark", "avg.rank", "cor.rank", "status")
emetrics = c("nreads", "NSC", "RSC", "BOTH_1.0", "IMPUTE_1.0_OBSERVED_5.0", 
"OBSERVED_1.0_IMPUTE_5.0", "Correlation", "IMPUTE_1.0_AUC_PREDICT_OBSERVE", "OBSERVED_1.0_AUC_PREDICT_IMPUTE") 
qcwide$NSC[qcwide$NSC > 4] = 4
qcwide$RSC[qcwide$RSC > 4] = 4

qclong = gather(qcwide[,c(ocols, emetrics)], metric, value, -id, -mark, -avg.rank, -cor.rank, -status)
qclong = merge(qclong, mapmetric)
# NOTE: Can't order just by id, must order by both id + mark:
qclong = qclong[order(qclong$avg.rank, decreasing=TRUE),]
qclong$cm = paste0(qclong$id, "_", qclong$mark)
qclong$cm = paste0(qclong$mark, "_", qclong$id)
qclong$cm = factor(qclong$cm, levels= unique(qclong$cm))

qclong = qclong[order(qclong$cor.rank, decreasing=TRUE),]
qclong$cm = factor(qclong$cm, levels= unique(qclong$cm))

gplot = ggplot(qclong, aes(cm, value, fill=status)) +
    facet_grid(short.metric ~ mark, scales='free') + 
    labs(x=paste0('Observed Sample + Mark'), y='QC Metric Value') +
    geom_bar(stat='identity', size=.1) +
    scale_fill_manual(values=vals) +
    theme_minimal() + 
    theme(axis.text.x = element_blank(), legend.position='none') + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'impute_qc_metrics_cor_qc_all',today,'.png'), gplot, dpi=300, width=14, height=7,units='in')
ggsave(paste0(imgpref, 'impute_qc_metrics_cor_qc_all',today,'.pdf'), gplot, dpi=300, width=14, height=7,units='in')


qcwide$NSC[qcwide$NSC > 4] = 4
qcwide$RSC[qcwide$RSC > 4] = 4

qcwide$Nlow = qcwide$NSC < 1.05
qcwide$Rlow = qcwide$RSC < 1.00

# Quartiles (text)
qcwide$qnn = qcwide$Nlow & qcwide$Rlow
qcwide$qpr = qcwide$Nlow & !(qcwide$Rlow)
qcwide$qpn = !(qcwide$Nlow) & qcwide$Rlow
qcwide$qpnr = !(qcwide$Nlow) & !(qcwide$Rlow)
qcwide$all = 1
qcnr = aggregate(cbind(qnn, qpr, qpn, qpnr, all) ~ status, qcwide, sum)
qcnr = gather(qcnr, quartile, count, -all, -status)
qcnr$pct = round(qcnr$count / qcnr$all * 100,2)
qcnr = merge(qcnr, data.frame(quartile=c('qnn','qpr','qpn','qpnr'), 
                              NSC = c(.5, .5, 2.5, 2.5), 
                              RSC = c(.5, 2.5, .5, 2.5)))
qcnr = merge(qcnr, data.frame(status=c('LOW','MEDIUM','TOP'), 
                              Radj=c(-.15,0,.15)))

gplot = ggplot(qcwide, aes(NSC, RSC, color=status)) + 
    geom_point(cex=.3, alpha=0.75) + 
    geom_text(data=qcnr, aes(NSC, RSC + Radj, label=paste0(pct, '%'), color=status)) +
    scale_color_manual(values=ovals) +
    theme_pubclean() + 
    xlim(0, max(qcwide$NSC)) + 
    ylim(0, max(qcwide$RSC)) + 
    geom_hline(yintercept=1) + geom_vline(xintercept=1.05) + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'RSC_NSC_scatter_all_',today,'.png'), gplot, dpi=300, width=5, height=5,units='in')
ggsave(paste0(imgpref, 'RSC_NSC_scatter_all_',today,'.pdf'), gplot, dpi=300, width=5, height=5,units='in')


# ---------------------
# Look at the low qual:
# ---------------------
aggregate(mark ~ id, lowmdf, function(x){paste(x, collapse=', ')})

nmid = aggregate(mark ~ id, lowmdf, length)
names(nmid)[2] = 'count'
lowqc = merge(qcwide, lowmdf)
lowqc = merge(nmid, lowqc)


png(paste0(imgpref, 'lowqual_vmaxcorr_',today,'.png'),res=300,units='in',width=12,height=7)
layout(matrix(c(1:2),1,2))
par(xaxs='i')
par(yaxs='i')
par(mar=c(4,4,1,1))
plot(lowqc$Correlation, 1 - lowqc$dist, pch=19, 
     xlim=c(0,1), ylim=c(0,1), col=rgb(0,0,0,.5), cex=.5,
     xlab='Correlation Obs-Imp', ylab='Max Correlation Obs-Any mark imputed')
lab = paste0(lowqc$cm, '\n', lowqc$against)
lab[(lowqc$Correlation >= (1 - lowqc$dist))] = ''
lab[lowqc$count == 1] = ''
text(lowqc$Correlation, 1 - lowqc$dist, lab, cex=.4, xpd=TRUE)
abline(0,1, lty='dashed')
par(mar=c(4,4,1,1))
plot(lowqc$Correlation, lowqc$maxCorr, pch=19,
     xlim=c(0,1), ylim=c(0,1), col=rgb(0,0,0,.5), cex=.5,
     xlab='Correlation Obs-Imp', ylab='Max Correlation Obs-Any sample imputed')
lab = paste0(lowqc$cm, '\n', lowqc$msamp)
lab[(lowqc$Correlation >= (lowqc$maxCorr))] = ''
lab[lowqc$count == 1] = ''
text(lowqc$Correlation, lowqc$maxCorr, lab, cex=.4, xpd=TRUE)
abline(0,1, lty='dashed')
dev.off()






# ------------
# ALT FIGURES:
# ------------
rqcdf = qcwide[,c('id','mark','Correlation','maxCorr','dist', 'against','msamp')]
rqcdf$dist = 1 - rqcdf$dist
names(rqcdf) = c('id','mark','corr','corr.msamp','corr.mmark', 'mmark','msamp')
rqcdf = merge(rqcdf, markdef, all.x=T)
rqcdf$corr.mmax = apply(rqcdf[,c('corr.msamp','corr.mmark')], 1,max)
rqcdf = merge(rqcdf, metricsdf[,c('id','mark','status')])
# rqcdf = rqcdf[rqcdf$mark %in% main.marks,]
sum(rqcdf$status == 'LOW')
lowdf = aggregate(status ~ mark, rqcdf, function(x){sum(x == 'LOW')})




rqcdf$lab = paste(rqcdf$mark, '-',rqcdf$mmark)
rqcdf$lab[rqcdf$corr.mmark - rqcdf$corr < .1] = ''

tvals = vals
tvals[2] = 'darkgrey'
# ggplot(rqcdf, aes(mark, corr.mmax - corr, color=status)) +
# ggplot(rqcdf, aes(mark, corr.mmark - corr, label=lab, color=status)) +
gp = ggplot(rqcdf, aes(mark, corr.msamp - corr, label=lab, color=status)) +
    scale_color_manual(values=tvals) +
    geom_jitter(cex=.1) + 
    geom_text(cex=1) + 
    theme_pubr() + xlab('') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(legend.position='none')
ggsave('~/test.png', gp, dpi=450, units='in', width=4, height=5)


ggplot(rqcdf, aes(mark, corr.mmark - corr, fill=type)) +
    geom_boxplot() + 
    theme_pubr()



ggplot(rqcdf, aes(mark, corr.msamp - corr, fill=type)) +
    geom_boxplot() + 
    theme_pubr()


# --------------
# Make heatmaps:
# --------------
# Correlation:
zlim = c(0,1 + 1e-6)
wide = spread(aggregate(Correlation ~ id + mark, qcwide, max), id, Correlation)
mat = as.matrix(wide[,-1])
rownames(mat) = wide$mark
subeporder = eporder[eporder %in% rownames(mat)]
subcellorder = cellorder[cellorder %in% colnames(mat)]

# Maxcorr:
mwide = spread(aggregate(maxCorr2 ~ id + mark, qcwide, max), id, maxCorr2)
mmat = as.matrix(mwide[,-1])
rownames(mmat) = mwide$mark
diffmat = mmat - mat
diffmat[diffmat > 0.4] = 0.5
dzlim = range(diffmat, na.rm=T)

qcwide$sampdist = (qcwide$maxCorr - qcwide$Correlation)
qcwide$logFCsamp = log2(qcwide$maxCorr / qcwide$Correlation)
tmpm = aggregate(sampdist ~ mark, qcwide, mean)
tmps = aggregate(sampdist ~ mark, qcwide, sd)
names(tmpm)[2] = 'sdmean'
names(tmps)[2] = 'sdsd'
qcwide = merge(merge(qcwide, tmpm), tmps)
qcwide$sampz = (qcwide$sampdist - qcwide$sdmean) / qcwide$sdsd
# NOTE: doesnt work very well, because some marks have no swaps - everything is well modeled.

qcwide$markdist = ((1 - qcwide$dist) - qcwide$Correlation)

# Origin:
dcut = 0.2
# qcwide$origin = (qcwide$sampz > 4) + 2 * (qcwide$markdist > dcut)
# qcwide$origin = (qcwide$sampdist > dcut) + 2 * (qcwide$markdist > dcut)
qcwide$origin = (qcwide$sampdist > dcut) + 2 * (qcwide$markdist > dcut)
owide = spread(aggregate(origin ~ id + mark, qcwide, min), id, origin)
omat = as.matrix(owide[,-1])
rownames(omat) = owide$mark

subqcwide = filter(qcwide, origin > 0)

# png(paste0(imgpref, txt, '_heatmap_vsmax_d',dcut,'.png'),res=300,units='in',width=7,height=3)
pdf(paste0(imgpref, txt, '_heatmap_vsmax_d',dcut,'.pdf'),width=7,height=3)
layout(matrix(c(1:3),3,1), TRUE)
par(xaxs='i')
par(yaxs='i')
par(mar=c(.15,2.75,1.5,.15))
image(t(mat[subeporder, subcellorder]),axes=F, col=rev(viridis_pal()(100)), zlim=zlim)
text(y=seq(0,1, length.out=nrow(mat)),
     x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]), 
     labels=subeporder, srt=0, adj=1, xpd=TRUE,cex=.4)
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
mtext('Imputed vs. Observed Correlation', side=3, cex=1)
par(mar=c(.15,2.75,1.5,.15))
image(t(diffmat[subeporder, subcellorder]),axes=F, col=col, zlim=dzlim)
text(y=seq(0,1, length.out=nrow(mat)),
     x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]), 
     labels=subeporder, srt=0, adj=1, xpd=TRUE,cex=.4)
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
mtext('Difference with maximum within-sample or within-mark correlation', side=3, cex=1)
par(mar=c(.15,2.75,1.5,.15))
image(t(omat[subeporder, subcellorder]),axes=F, 
      col=c('lightgrey','indianred','royalblue', 'mediumpurple4'), zlim=c(0,3))
text(y=seq(0,1, length.out=nrow(mat)),
     x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]), 
     labels=subeporder, srt=0, adj=1, xpd=TRUE,cex=.4)
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
mtext('Flagged samples (red/blue=Sample/Ab swap)', side=3, cex=1)
dev.off()


# With flagged samples, eval imputation quality:
lmodel = lm(lgtCorr ~ mark, data=qcwide)

# Can evaluate the feeding datasets by the similarity of datasets with feeding information.

labwide = qcwide[qcwide$maxCorr2/qcwide$Correlation > 1.40,]
ggplot(qcwide, aes(Correlation, difMax, color=closeby)) +
    geom_point(alpha=0.3) +
    scale_color_manual(values=c('red','black', 'black')) + 
    geom_text_repel(data=labwide, aes(Correlation, difMax, label=paste0(id, ' - ', mark, '\n(to ',assign,')')), cex=1.75) +
    theme_minimal()
ggsave(paste0(imgpref, 'corr_maxcorr_bothswap_',today,'.pdf'), dpi=300, width=11, height=8,units='in')

labwide = filter(qcwide, maxCorr2 - Correlation  > 0.05 & closeby == 2)
ggplot(qcwide, aes(Correlation, difMax, color=closeby)) +
    geom_point(alpha=0.3) +
    scale_color_manual(values=c('grey','black', 'black')) + 
    geom_text_repel(data=labwide, aes(Correlation, difMax, label=paste0(id, ' - ', mark, '\n(to ',assign,')')), cex=1.75) +
    theme_minimal()
ggsave(paste0(imgpref, 'corr_maxcorr_markswap_',today,'.pdf'), dpi=300, width=11, height=8,units='in')

# Automatically find (easy to assign) lines:
segdf = c()
for (i in 1:nrow(qcwide)){
    os = qcwide$id[i]
    ss = qcwide$msamp[i]
    om = qcwide$mark[i]
    sm = qcwide$against[i]
    pt1 = c(qcwide$Correlation[i], qcwide$difMax[i])
    # 1) Potentially swapped samples:
    pts = filter(qcwide, id == ss & msamp == os & mark == om)
    if (nrow(pts) == 1){
        pt2 = c(pts$Correlation[1], pts$difMax[1])
        segdf = rbind(segdf, c(pt1, pt2, 1))
    }
    # 2) Potentially wapped marks:
    # NOTE: H3K4me3 strongly corr with H3K9ac.
    ptm = filter(qcwide, mark == sm & against == om & id == os)
    if (nrow(ptm) == 1){
        pt3 = c(ptm$Correlation[1], ptm$difMax[1])
        segdf = rbind(segdf, c(pt1, pt3, 2))
    }
}
segdf = data.frame(segdf)
names(segdf) = c('x1','y1','x2','y2', 'closeby')
segdf$closeby = factor(segdf$closeby, levels=levels(qcwide$closeby))
segdf = filter(segdf, x1 != x2 & y1 !=y2)

labwide = filter(qcwide, maxCorr2 - Correlation  > 0.05 & closeby == 2)
ggplot(qcwide, aes(Correlation, difMax, color=closeby)) +
    scale_color_manual(values=c('grey','black', 'black')) + 
    geom_segment(data=segdf, aes(x=x1, xend=x2, y=y1, yend=y2)) +
    geom_point(alpha=0.3, cex=2) +
    geom_text_repel(data=labwide, aes(Correlation, difMax, label=paste0(id, ' - ', mark, '\n(to ',assign,')')), cex=1.75) +
    theme_minimal()

# Figure out which it is closest to.


# # Plot boxplots of corr + flags.
# ggplot(corrdf, aes(mark, Correlation - pmain, color=badcls)) +
#     geom_violin() +
#     geom_jitter() +
#     theme_minimal()


# -------------------------------------
# Plot metrics with the flagged tracks:
# -------------------------------------
mlong = gather(metricsdf, metric, value, -id, -mark, -avg.rank, -cor.rank, -status)
mlong = merge(mlong, map)
mlong = merge(mlong, mapmetric)

# Punctate vs. broad:
pm = c('ATAC-seq','DNase-seq','H2AFZ','H3K27ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac')
bm = c('H3K27me3','H3K36me3','H3K79me2','H3K9me3','H4K20me1')
fc = c('CTCF','EP300','POLR2A','SMC3','RAD21')
markdef = rbind(data.frame(mark=pm, type='punctate'), 
                data.frame(mark=bm, type='broad'), 
                data.frame(mark=fc, type='factor'))
mlong = merge(mlong, markdef, all.x=T)
# Averages:
mavgdf = aggregate(value ~ short.metric + mark, mlong, mean)
mavgdf = mavgdf[order(mavgdf$short.metric),]

tavgdf = aggregate(value ~ short.metric + type, mlong, mean)
tavgdf = tavgdf[order(tavgdf$short.metric),]
print(tavgdf)

# NOTE: Can't order just by id, must order by both 
mlong = mlong[order(mlong$avg.rank, decreasing=TRUE),]
# mlong$id = factor(mlong$id, levels=unique(mlong$id)) 
mlong$cm = paste0(mlong$id, "_", mlong$mark)
mlong$cm = paste0(mlong$mark, "_", mlong$id)
mlong$cm = factor(mlong$cm, levels= unique(mlong$cm))
NSAMP = length(unique(mlong$id))
NMARK = length(unique(metricsdf$mark))

# Write out the status:
outdf = unique(metricsdf[,c('id','mark','status')])
outdf = outdf[order(outdf$id, outdf$mark),]
write.table(outdf, paste0('best_worst_impobs_agreement_', today, '.tsv'), sep="\t", col.names=F, row.names=F, quote=F)


mlong = mlong[mlong$mark %in% c('DNase-seq','H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9me3'),]

# Plot without colors:
gplot = ggplot(mlong, aes(cm, value)) +
    facet_grid(short.metric ~ mark, scales='free') + 
    labs(x=paste0('Observed Sample + Mark'), y='AUC/% Recovery/Correlation') +
    geom_bar(stat='identity', color='darkgrey') +
    theme_minimal() + 
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'impute_qc_metrics_',today,'.png'), gplot, dpi=300, width=11, height=8,units='in')
ggsave(paste0(imgpref, 'impute_qc_metrics_',today,'.pdf'), gplot, dpi=300, width=11, height=8,units='in')

# Plot with labels:
mlong = mlong[order(mlong$cor.rank, decreasing=TRUE),]
mlong$cm = factor(mlong$cm, levels= unique(mlong$cm))
gplot = ggplot(mlong, aes(cm, value, fill=status, color=status)) +
    facet_grid(short.metric ~ mark, scales='free') + 
    labs(x=paste0('Observed Sample + Mark'), y='AUC/% Recovery/Correlation') +
    geom_bar(stat='identity', size=.1) +
    scale_fill_manual(values=vals) +
    scale_color_manual(values=vals) +
    theme_minimal() + 
    theme(axis.text.x = element_blank(), legend.position='none') + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'impute_qc_metrics_corrank_',today,'.png'), gplot, dpi=300, width=11, height=5,units='in')
ggsave(paste0(imgpref, 'impute_qc_metrics_corrank_',today,'.pdf'), gplot, dpi=300, width=11, height=5,units='in')

mlong$value = as.numeric(mlong$value)
mlong = mlong[order(mlong$cor.rank, decreasing=TRUE),]
mlong$cm = factor(mlong$cm, levels= unique(mlong$cm))
sublong = mlong[mlong$metric == 'Correlation',]
gplot = ggplot(sublong, aes(cm, value, fill=status, color=status)) +
    facet_wrap(~ mark, scales='free') + 
    labs(x=paste0('Observed Sample + Mark'), y='Correlation (Observed to Imputed)') +
    geom_bar(stat='identity', size=.1) +
    lims(y=c(0,1)) +
    scale_fill_manual(values=vals) +
    scale_color_manual(values=vals) +
    theme_minimal() + 
    theme(axis.text.x = element_blank(), legend.position='none') + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'impute_qc_metrics_onlycor_',today,'.png'), gplot, dpi=300, width=7, height=2.5,units='in')
ggsave(paste0(imgpref, 'impute_qc_metrics_onlycor_',today,'.pdf'), gplot, dpi=300, width=7, height=2.5,units='in')

mlong = mlong[order(mlong$cor.rank, decreasing=TRUE),]
mlong$cm = factor(mlong$cm, levels= unique(mlong$cm))
sublong = mlong[mlong$metric == 'Correlation',]
gplot = ggplot(sublong, aes(cm, value, fill=status, color=status)) +
    facet_wrap(~ mark, scales='free',ncol=1) + 
    labs(x=paste0('Observed Sample + Mark'), y='Correlation (Observed to Imputed)') +
    geom_bar(stat='identity', size=.1) +
    lims(y=c(0,1)) +
    scale_fill_manual(values=vals) +
    scale_color_manual(values=vals) +
    theme_minimal() + 
    theme(axis.text.x = element_blank(), legend.position='none') + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'impute_qc_metrics_onlycor_tall_',today,'.png'), gplot, dpi=300, width=2.5, height=7,units='in')
ggsave(paste0(imgpref, 'impute_qc_metrics_onlycor_tall_',today,'.pdf'), gplot, dpi=300, width=2.5, height=7,units='in')

# -------------------------------------------------------------------------------------------
# Evaluate other QC metrics (# reads, NSC, RSC) for the poor quality tracks:
# -------------------------------------------------------------------------------------------
# Read in the QC file:
qcdf = read.delim(qcfile, header=F)
names(qcdf) = c('file','nreads','estFragLen', 'corr_estFragLen','PhantomPeak','corr_phantomPeak','argmin_corr','min_corr','NSC','RSC','QualityTag')
mat = t(sapply(qcdf$file, function(x){
                 x = sub(".sub.*","",sub("FINAL_","",x)); 
                 strsplit(x,"_")[[1]]}))
qcdf$id=mat[,2]
qcdf$mark=mat[,1]
qcwide = merge(qcdf[,c('mark','id','nreads','NSC','RSC', 'QualityTag')], metricsdf)

qclong = gather(qcwide, metric, value, -id, -mark, -avg.rank, -cor.rank, -status, -QualityTag)
qclong = merge(qclong, map)
qclong = merge(qclong, mapmetric)
# NOTE: Can't order just by id, must order by both 
qclong = qclong[order(qclong$avg.rank, decreasing=TRUE),]
# qclong$id = factor(qclong$id, levels=unique(qclong$id)) 
qclong$cm = paste0(qclong$id, "_", qclong$mark)
qclong$cm = paste0(qclong$mark, "_", qclong$id)
qclong$cm = factor(qclong$cm, levels= unique(qclong$cm))
# sublong = mlong[mlong$metric == 'Correlation',]

qclong = qclong[order(qclong$cor.rank, decreasing=TRUE),]
qclong$cm = factor(qclong$cm, levels= unique(qclong$cm))
sublong = qclong[qclong$metric %in%  c('Correlation', 'nreads','NSC','RSC'),]

gplot = ggplot(sublong, aes(cm, value, fill=status, color=status)) +
    facet_grid(short.metric ~ mark, scales='free') + 
    labs(x=paste0('Observed Sample + Mark'), y='QC Metric Value') +
    geom_bar(stat='identity', size=.1) +
    scale_fill_manual(values=vals) +
    scale_color_manual(values=vals) +
    theme_minimal() + 
    theme(axis.text.x = element_blank(), legend.position='none') + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'impute_qc_metrics_cor_qc_',today,'.png'), gplot, dpi=300, width=10, height=4,units='in')
ggsave(paste0(imgpref, 'impute_qc_metrics_cor_qc_',today,'.pdf'), gplot, dpi=300, width=10, height=4,units='in')


gplot = ggplot(qclong, aes(cm, value, fill=status, color=status)) +
    facet_grid(short.metric ~ mark, scales='free') + 
    labs(x=paste0('Observed Sample + Mark'), y='QC Metric Value') +
    geom_bar(stat='identity', size=.1) +
    scale_fill_manual(values=vals) +
    scale_color_manual(values=vals) +
    theme_minimal() + 
    theme(axis.text.x = element_blank(), legend.position='none') + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'impute_qc_metrics_all_',today,'.png'), gplot, dpi=300, width=10, height=6,units='in')
ggsave(paste0(imgpref, 'impute_qc_metrics_all_',today,'.pdf'), gplot, dpi=300, width=10, height=6,units='in')

# aka rsc vs. nsc, colored by corr.

qcwide$status = factor(qcwide$status, levels=c('MEDIUM','TOP','LOW'))
qcwide = qcwide[order(qcwide$status),]
qcwide$col = vals[as.character(qcwide$status)]


gplot = ggplot(qcwide, aes(NSC, RSC, color=status)) + 
    facet_wrap(~ mark) + geom_point() + 
    scale_x_log10() +
    scale_color_manual(values=vals) +
    theme_minimal() + 
    geom_hline(yintercept=1) + geom_vline(xintercept=1.05) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(paste0(imgpref, 'RSC_NSC_scatter_all_',today,'.png'), gplot, dpi=300, width=10, height=6,units='in')

png(paste0(imgpref, 'RSC_NSC_scatter.png'),res=300,units='in',width=7,height=6)
par(mar=c(3,3,1,1))
plot(qcwide$NSC, qcwide$RSC, col=qcwide$col, 
     pch=19, ylab='',xlab='', cex=.5)
abline(h=1, lty='dashed')
abline(v=1.05, lty='dashed')
mtext('NSC',side=1, line=2, cex=1.2)
mtext('RSC',side=2, line=2, cex=1.2)
dev.off()


qcwide$lab = paste0(qcwide$id,' - ',qcwide$mark)
qcwide$lab[qcwide$NSC < 1.05] = '' 
qcwide$lab[qcwide$RSC < 1] = ''
qcwide$lab[qcwide$status != 'LOW'] = ''

png(paste0(imgpref, 'RSC_NSC_scatter_lab.png'),res=300,units='in',width=7,height=6)
par(mar=c(3,3,1,1))
plot(qcwide$NSC, qcwide$RSC, col=qcwide$col, 
     pch=19, ylab='',xlab='', cex=.5)
text(qcwide$NSC, qcwide$RSC, labels=qcwide$lab, 
     col=qcwide$col, adj=-.1, cex=.75, font=2)
abline(h=1, lty='dashed')
abline(v=1.05, lty='dashed')
mtext('NSC',side=1, line=2, cex=1.2)
mtext('RSC',side=2, line=2, cex=1.2)
dev.off()




# Show what quality metrics look like.
# Show bad tracks identified even when metrics are bad.


# # Can plot as heatmap as well (separate by histonemark?)
# metricsdf = metricsdf[order(metricsdf$mark, metricsdf$avg.rank, decreasing=TRUE),]
# mat = as.matrix(metricsdf[,metrics])
# norm = rep(100, length(metrics))
# norm[apply(mat, 2, max) <= 1] = 1
# mat = sweep(mat, 2, norm, '/')
# image(mat, axes=F, col=col, zlim=c(0,1))
# mtext(paste(unique(metricsdf$mark), collapse="\t"),3)



# samplesheet <- read.delim('ChromImpute/sample_mark_table.tsv',header=F)
# names(samplesheet) <- c('BSSID','epitope','file')
# samplesheet <- merge(samplesheet, map)
# samplesheet$file = 1
# widesample = spread(samplesheet, epitope, file)
# widesample[is.na(widesample)] <- 0
# rownames(widesample) = widesample$BSSID
