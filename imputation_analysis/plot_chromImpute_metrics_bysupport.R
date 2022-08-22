#!/usr/bin/R
# ----------------------------------------------------
# Plot chromImpute metrics for imputing existing files
# Against the amount of supporting experiments
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
markdef = rbind(data.frame(mark=pm, type='Punctate'), 
                data.frame(mark=bm, type='Broad'), 
                data.frame(mark=fc, type='factor'))
metricsdf = merge(metricsdf, markdef)
aggregate(mark ~ type, metricsdf, length)

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

# Make wide + long qc tables: 
qcwide = merge(qcdf, metricsdf)
ocols = c("id", "mark", "avg.rank", "cor.rank", "status")
emetrics = c("BOTH_1.0", "IMPUTE_1.0_OBSERVED_5.0", 
             "OBSERVED_1.0_IMPUTE_5.0", "Correlation",
             "IMPUTE_1.0_AUC_PREDICT_OBSERVE", "OBSERVED_1.0_AUC_PREDICT_IMPUTE") 
qcwide$NSC[qcwide$NSC > 4] = 4
qcwide$RSC[qcwide$RSC > 4] = 4
qclong = gather(qcwide[,c(ocols, emetrics)], metric, value, -id, -mark, -avg.rank, -cor.rank, -status)
qclong = merge(qclong, mapmetric)

# Quantify number of supporting datasets, per cell + assay:
wmat = as.matrix(wm[,-1])
marg.assay = apply(wmat, 1, sum)
marg.cell = apply(wmat, 2, sum)
marg.assay = data.frame(mark=names(marg.assay), marg.assay = marg.assay)
marg.cell = data.frame(id=names(marg.cell), marg.cell = marg.cell)


qcwide = merge(merge(qcwide, marg.assay), marg.cell)
qclong = merge(merge(qclong, marg.assay), marg.cell)
qclong = merge(qclong, markdef, all.x=TRUE)
qclong$marg.cell.short = as.character(qclong$marg.cell - 1)
qclong$marg.cell.short[(qclong$marg.cell - 1) >= 10] = '10+'
qclong$marg.cell.short = factor(qclong$marg.cell.short, levels=c(as.character(1:9), '10+'))
qclong$type = factor(qclong$type, levels=c('Punctate','Broad'))

gplot = ggplot(qclong, aes(factor(marg.cell.short), value, color=type)) + 
    facet_wrap(~short.metric, scales='free_y', nrow=3) + 
    geom_boxplot(position=position_dodge(.85), width=.6, fill=NA, outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.45, dodge.width=.85), cex=.1) + 
    theme_pubr() + 
    scale_color_manual(values=c('Broad'='indianred', 'Punctate'='slateblue'), name='Assay/Mark Type') + 
    labs(x='Number of supporting datasets from the same sample', y='QC Metric Value') + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'metrics_by_support_type.pdf'), gplot, dpi=300, width=6, height=5,units='in')
ggsave(paste0(imgpref, 'metrics_by_support_type.png'), gplot, dpi=300, width=6, height=5,units='in')


qclong$marg.assay.short = qclong$marg.assay - 1
qclong$mark.wmarg = paste0(qclong$mark, ' (', qclong$marg.assay - 1, ')')
qclong = qclong[order(qclong$marg.assay),]
qclong = qclong[order(qclong$type),]
qclong$mark.wmarg = factor(qclong$mark.wmarg, levels = unique(qclong$mark.wmarg))

gplot = ggplot(qclong, aes(mark.wmarg, value, color=type)) + 
    facet_wrap(~short.metric, scales='free_y') + 
    geom_boxplot(position=position_dodge(.85), width=.6, fill=NA, outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.45, dodge.width=.85), cex=.1) + 
    theme_pubr() + 
    scale_color_manual(values=c('slateblue', 'indianred'), name='Assay/Mark Type') + 
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    labs(x='Number of supporting datasets from the same assay type', y='QC Metric Value') + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'metrics_by_support_assay_type.pdf'), gplot, dpi=300, width=8, height=4,units='in')
ggsave(paste0(imgpref, 'metrics_by_support_assay_type.png'), gplot, dpi=300, width=8, height=4,units='in')


# -------------------------------------------------------------------------------
# For each point - find the avg. sim of the closest 10 cells to the tested track:
# -------------------------------------------------------------------------------
# Load matrices
source(paste0(bindir, 'load_distance_matrices.R'))

# Calculate nearest neighbors:
all.m10df = NULL
for (mark in names(full.ll)){
    mat = full.ll[[mark]]
    oidx = grep("obs", colnames(mat))
    # For each in tested:
    comp.id = qcwide$id[qcwide$mark == mark]
    obs.mat = 1 - mat[oidx,oidx]
    # Mean of 10 closest cells:
    cmx = sapply(paste0(comp.id,'_obs'), function(x){
                     x = head(sort(obs.mat[x,], decreasing=T),10)
                     mean(x) })
    m10df = data.frame(id=comp.id, mark=mark, m10close=cmx)
    all.m10df = rbind(all.m10df, m10df)
}

qcwide = merge(qcwide, all.m10df)
qclong = merge(qclong, all.m10df)


gplot = ggplot(qclong, aes(m10close, value, color=type)) + 
    facet_wrap(~short.metric, scales='free_y', nrow=3) + 
    # geom_point(cex=.5, alpha=.5) + 
    geom_point(cex=.25) + 
    theme_pubr() + 
    scale_color_manual(values=c('slateblue', 'indianred'), name='Assay/Mark Type') + 
    labs(x='Average correlation of 10 closest datasets with observed dataset', y='QC Metric Value') + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'metrics_10closest_type_scatter.pdf'), gplot, dpi=300, width=6, height=5,units='in')
ggsave(paste0(imgpref, 'metrics_10closest_type_scatter.png'), gplot, dpi=300, width=6, height=5,units='in')

cordf = qclong[qclong$short.metric == 'GWCorr',]
mean(cordf$value > cordf$m10close) # 95.3% of time correlation is greater with imputed data.
mean(cordf$value / cordf$m10close) # 36.3% higher on average


# Plot this against the # of cells:
qclong$mcell.bins = cut(qclong$marg.cell - 1, breaks=c(1,2,4,6,8,10,28), include.lowest=T)
colb.full <- colorRampPalette(brewer.pal(n=9,name="Blues"))(100) 
colb.10 = colb.full[seq(40, 100, length.out = length(levels(qclong$mcell.bins)))]

gplot = ggplot(qclong, aes(m10close, value, color=mcell.bins)) + 
    facet_wrap(~short.metric, scales='free_y',nrow=3) + 
    geom_point(cex=.25, alpha=1) + 
    theme_pubr() + 
    scale_color_manual(values=colb.10, name='# Supporting\nSamples:') + 
    labs(x='Average correlation of 10 closest datasets with observed dataset', y='QC Metric Value') + 
    guides(color = guide_legend(override.aes=list(size=4, pch=15)))
ggsave(paste0(imgpref, 'metrics_by_support_10closest_type_scatter.pdf'), gplot, dpi=300, width=6, height=5,units='in')
ggsave(paste0(imgpref, 'metrics_by_support_10closest_type_scatter.png'), gplot, dpi=300, width=6, height=5,units='in')
ggsave(paste0(imgpref, 'metrics_by_support_10closest_type_scatter_nolegend.png'), gplot + theme(legend.position='none'),
       dpi=300, width=6, height=5,units='in')


qclong$m10.bins = cut(qclong$m10close, breaks=seq(0,1,by=.25), include.lowest=F)
gplot = ggplot(qclong, aes(m10.bins, value, color=mcell.bins)) + 
    facet_wrap(~short.metric, scales='free_y',nrow=3) + 
    geom_boxplot(position=position_dodge(.85), width=.7, fill=NA, outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.35, dodge.width=.85), cex=.1) + 
    theme_pubr() + 
    scale_color_manual(values=colb.10, name='# Supporting\ndatasets in\nsame sample') + 
    labs(x='Average correlation of 10 closest datasets with observed dataset', y='QC Metric Value') + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    theme(legend.position='none')
ggsave(paste0(imgpref, 'metrics_by_support_10closest_type_boxplot.pdf'), gplot, dpi=300, width=6, height=5,units='in')
ggsave(paste0(imgpref, 'metrics_by_support_10closest_type_boxplot.png'), gplot, dpi=300, width=6, height=5,units='in')

# ----------------------------------------------------------------
# Show how many are based on one expt. and their level of support:
# ----------------------------------------------------------------
single.id = as.character(marg.cell$id[marg.cell$marg.cell == 1])
single.id = single.id[single.id %in% cellorder]
print(length(single.id))
obs.id = paste0(single.id,'_obs')

# Get the sample closeness metric:
# - closeness for single samples
# - closeness for other samples
# - closeness for single samples to non-single samples
# - closensess for other sample to non-single samples
single.m10df = NULL
for (mark in names(full.ll)){
    mat = full.ll[[mark]]
    oidx = grep("obs", colnames(mat))
    # For each in tested:
    obs.mat = 1 - mat[oidx,oidx]
    onam = colnames(obs.mat)
    ind = which(onam %in% obs.id)
    nind = which(!(onam %in% obs.id))
    # Mean of 10 closest cells:
    if (length(ind) > 0){
        ci = sapply(onam[ind], function(x){ mean(head(sort(obs.mat[x,], decreasing=T),10)) })
        cn = sapply(onam[nind], function(x){ mean(head(sort(obs.mat[x,], decreasing=T),10)) })
        cin = sapply(onam[ind], function(x){ mean(head(sort(obs.mat[x,nind], decreasing=T),10)) })
        cnn = sapply(onam[nind], function(x){ mean(head(sort(obs.mat[x,nind], decreasing=T),10)) })
        vals = c(ci, cn, cin, cnn)
        ids = sub("_obs",'', names(vals))
        sups = c(rep('single', length(ind)), 
                 rep('non-single', length(nind)), 
                 rep('single', length(ind)), 
                 rep('non-single', length(nind)))
        types = c(rep('all', length(ind)), 
                  rep('all', length(nind)), 
                  rep('non-single', length(ind)), 
                  rep('non-single', length(nind)))
        m10df = data.frame(id=ids, mark=mark, m10close=vals, 
                           support=sups, against=types)
        single.m10df = rbind(single.m10df, m10df)
    }
}
rownames(single.m10df) = NULL
single.m10df = merge(single.m10df, markdef)
single.m10df$comp = with(single.m10df, paste0(support, ' vs. ', against))

sub.m10df = single.m10df[single.m10df$against == 'non-single',]
sub.m10df$status = 'Single Expt.'
sub.m10df$status[sub.m10df$support == 'non-single'] = 'Multiple Expt.'

gplot = ggplot(sub.m10df, aes(m10close, mark, color=status)) + 
    geom_boxplot(position=position_dodge(.75), width=.6, outlier.shape=NA) + 
    geom_jitter(position=position_jitterdodge(jitter.width=.1, dodge.width=.75), cex=.25) + 
    theme_pubr() + lims(x=c(0,1)) + 
    scale_color_manual(values=c(colpair[10], colpair[60]), name='Sample:') + 
    labs(y='Assay/Histone Mark', x='Average correlation of 10 closest datasets\nfrom samples with multiple experiments') + 
    theme(legend.position='right')
ggsave(paste0(imgpref, 'singledata_comparisontypes_avg10closest_boxplot.pdf'), gplot, dpi=300, width=5.5, height=3,units='in')
ggsave(paste0(imgpref, 'singledata_comparisontypes_avg10closest_boxplot.png'), gplot, dpi=300, width=5.5, height=3,units='in')

single.m10df = sub.m10df[sub.m10df$support == 'single' & sub.m10df$against == 'non-single',]

aggregate(m10close ~ mark, single.m10df, function(x){mean(x >= .5)})
mean(single.m10df$m10close >= .5)

# TODO: Color by actual type of data:

# TODO: See how close the DNase-seq that is unobserved is to others
# TODO: Show that if DNase-seq is close, others are close.


