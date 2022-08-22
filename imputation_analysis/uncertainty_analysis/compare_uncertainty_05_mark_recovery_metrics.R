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
library(ggrepel)
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
imgpref = paste0(imgdir, 'mark_recovery_')

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

#' Aggregate and rename the aggregated column
#'
#' @param formula as in aggregate
#' @param data as in aggregate
#' @param FUN as in aggregate
#' @param name name for the aggregated column e.g. total
#' @return aggregate result with column renamed
#' @export
agg.rename = function(formula, data, FUN, name){
    agg.df = aggregate(formula, data, FUN)
    names(agg.df)[ncol(agg.df)] = name
    return(agg.df) }


# -----------
# Load files:
# -----------
rcdir = 'ChromImpute/recovery_metrics/'
fnames = list.files(path=rcdir, pattern='recovery_metrics_.*_chr.*_stats.tsv.gz')

statsdf = c()
for (file in fnames){
    mark = sub('_chr.*', '', sub('recovery_metrics_', '', file))
    chrom = sub('_stats.*', '', sub('recovery_metrics_.*chr', 'chr', file))
    df = read.delim(paste0(rcdir, file))
    df$mark = mark
    df$chrom = chrom
    statsdf = rbind(statsdf, df)
}

# ---------------------------------
# Read in flagged samples + remove:
# ---------------------------------
abdf = read.delim('Annotation/flagged_potential_abswap.tsv')
smdf = read.delim('Annotation/flagged_potential_sampleswap.tsv')
lwdf = read.delim('Annotation/flagged_low_agreement.tsv')
flagdf = rbind(abdf[,c(1,2)], smdf[,c(1,2)], lwdf)
flagdf$is.flagged = 1

statsdf = merge(statsdf, flagdf, all.x=TRUE)
statsdf = statsdf[is.na(statsdf$is.flagged),]
statsdf = merge(statsdf, meta[,c('id','GROUP','COLOR', 'infoline')], all.x=TRUE)

infomap = meta[,c('id','infoline')]
names(infomap) = c('near_id','near_infoline')
statsdf = merge(statsdf, infomap, all.x=TRUE)
statsdf = statsdf[statsdf$id %in% cellorder,]

# -------------------------------------
# Plot the statistics on these samples:
# -------------------------------------
gplot = ggplot(statsdf, aes(mark, aurocq)) + 
    geom_boxplot(fill='grey75') + 
    labs(x='Mark or Assay',y='AUROC') + 
   theme_pubr() + 
   theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'auroc_overall_boxplots.png'), gplot, dpi=400, units='in', width=5, height=7)

gplot = ggplot(statsdf, aes(mark, apq)) + 
    geom_boxplot(fill='grey75') + 
    labs(x='Mark or Assay',y='Average Precision') + 
    theme_pubr() + 
    theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'auprc_overall_boxplots.png'), gplot, dpi=400, units='in', width=5, height=7)


# ------------------------------------------------
# Plot a scatter of average precision comparisons:
# ------------------------------------------------
statsdf$near.lab = paste0(statsdf$infoline, '\n', statsdf$near_infoline)
statsdf$near.lab[statsdf$apq * 1.25 > statsdf$near_apq] = ''
pct.over1 = mean(statsdf$apq > statsdf$near_apq)
g1 = ggplot(statsdf, aes(near_apq, apq, color=GROUP)) + 
    geom_abline(intercept=0,slope=1) + 
    geom_point(cex=.5) + 
    geom_text(aes(x=.15, y=.85, label=paste0(round(100 * pct.over1, 1), '%\nof datasets')), col='grey65') +
    geom_text(aes(x=.85, y=.1, label=paste0(round(100 * (1 - pct.over1) , 1), '%\nof datasets')), col='grey65') +
    geom_text_repel(aes(label=near.lab), cex=1.25, segment.size=.25, box.padding=.2, label.padding=2, max.iter=4000) +
    scale_color_manual(values=colvals$group) + 
    labs(x='AP of Nearest Observed Dataset',y='AP of Imputed Track') + 
    theme_pubr() + theme(legend.position='none')
    # theme(axis.text.x=element_text(angle=90, hjust=1))

statsdf$mean.lab = statsdf$infoline
statsdf$mean.lab[statsdf$apq * 1.1 > statsdf$mean_apq] = ''
pct.over2 = mean(statsdf$apq > statsdf$mean_apq)
g2 = ggplot(statsdf, aes(mean_apq, apq, color=GROUP)) + 
    geom_abline(intercept=0,slope=1) + 
    geom_point(cex=.5) + 
    geom_text(aes(x=.15, y=.85, label=paste0(round(100 * pct.over2, 1), '%\nof datasets')), col='grey65') +
    geom_text(aes(x=.85, y=.1, label=paste0(round(100 * (1 - pct.over2) , 1), '%\nof datasets')), col='grey65') +
    geom_text_repel(aes(label=mean.lab), cex=1.5, lwd=.5) +
    scale_color_manual(values=colvals$group) + 
    labs(x='AP of Mean Observed Data',y='AP of Imputed Track') + 
    theme_pubr() + theme(legend.position='none')

garr = ggarrange(g1, g2)#   + rremove('y.text') + rremove('y.axis') + rremove('ylab') + rremove('y.ticks'))
ggsave(paste0(imgpref, 'auprc_scatter_comparison.png'), garr, dpi=400, units='in', width=8, height=4)

statsdf[statsdf$apq < statsdf$mean_apq,c('infoline','near_infoline')]

mean(statsdf$apq / statsdf$mean_apq)
mean(statsdf$apq / statsdf$near_apq)

statsdf$near_ratio = statsdf$apq / statsdf$near_apq
statsdf$mean_ratio = statsdf$apq / statsdf$mean_apq
rdf = aggregate(cbind(near_ratio, mean_ratio) ~ GROUP, statsdf, mean)
rdf = rdf[order(rdf$near_ratio, decreasing=T),]

statsdf$beats_mean = statsdf$apq > statsdf$mean_apq
statsdf$beats_near = statsdf$apq > statsdf$near_apq
rdf = aggregate(cbind(beats_near, beats_mean) ~ GROUP, statsdf, mean)

mean(statsdf$apq / statsdf$mean_apq)
mean(statsdf$apq / statsdf$near_apq)

tab = statsdf[statsdf$apq < statsdf$near_apq,]
aggregate(id ~ mark, tab, length)
stab = tab[!is.na(tab$near_infoline),]
mean(stab$infoline == stab$near_infoline)
# In 49% of cases, infoline is exactly the same - near replicates.
stab = stab[stab$mark == 'DNase-seq',]
mean(stab$infoline == stab$near_infoline)

head(stab[,c('infoline','near_infoline')],50)


# ----------------------------
# Gather and plot as boxplots:
# ----------------------------
slong = gather(statsdf[,c('mark','id','apq', 'mean_apq','near_apq')], type, value, - mark, -id)
along = gather(statsdf[,c('mark','id','aurocq', 'mean_aurocq','near_aurocq')], type, value, - mark, -id)

gplot = ggplot(along, aes(mark, value, fill=type)) + 
    geom_boxplot() + 
    labs(x='Mark or Assay',y='AUROC') + 
    theme_pubr() + 
    theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'auroc_comparison_boxplots.png'), gplot, dpi=400, units='in', width=9, height=9)


gplot = ggplot(slong, aes(mark, value, fill=type)) + 
    geom_boxplot() + 
    labs(x='Mark or Assay',y='Average Precision') + 
    theme_pubr() + 
    theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'ap_comparison_boxplots.png'), gplot, dpi=400, units='in', width=9, height=9)


# --------------------------------------
# Plot as bar charts, both AP and AUROC:
# --------------------------------------
mdf = agg.rename(value ~ type + mark, slong, mean, 'mean')
sdf = agg.rename(value ~ type + mark, slong, sd, 'sd')
ndf = agg.rename(value ~ type + mark, slong, length, 'n')

mdf = merge(merge(mdf, sdf), ndf)
mdf$se = mdf$sd / sqrt(mdf$n)
mdf$mark = factor(mdf$mark, levels=t12marks)

maptype = data.frame(type=c('apq','mean_apq','near_apq'),
                     type.nam=c('Imputed','Mean of Observed','Nearest Observed Dataset'))
mdf = merge(mdf, maptype)
# TODO: Order the marks (tier 1 vs tier 2?)
col.paired = brewer.pal(12, 'Paired')

g1 = ggplot(mdf, aes(mark, mean, fill=type.nam)) + 
    geom_hline(yintercept=seq(.2, .8, by=.2), lty='dotted', lwd=.5, col='grey75') +
    geom_bar(stat='identity', position='dodge') + 
    geom_errorbar(aes(ymin=mean-2 * se, ymax=mean + 2 * se), width=.2,
                  position=position_dodge(.9))  + 
    labs(x='Mark or Assay',y='Average Precision on Observed') + 
    theme_pubr() + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    # scale_fill_manual(values=c('indianred','grey70','grey85'), name='Predictor') + 
    scale_fill_manual(values=c(col.paired[5],'grey85',col.paired[1]), name='Predictor') + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(paste0(imgpref, 'ap_comparison_bars.png'), g1, dpi=400, units='in', width=6, height=4)


mdf = agg.rename(value ~ type + mark, along, mean, 'mean')
sdf = agg.rename(value ~ type + mark, along, sd, 'sd')
ndf = agg.rename(value ~ type + mark, along, length, 'n')

mdf = merge(merge(mdf, sdf), ndf)
mdf$se = mdf$sd / sqrt(mdf$n)
mdf$mark = factor(mdf$mark, levels=t12marks)

maptype = data.frame(type=c('aurocq','mean_aurocq','near_aurocq'),
                     type.nam=c('Imputed','Mean of Observed','Nearest Observed Dataset'))
mdf = merge(mdf, maptype)
# TODO: Order the marks (tier 1 vs tier 2?)
col.paired = brewer.pal(12, 'Paired')

g2 = ggplot(mdf, aes(mark, mean, fill=type.nam)) + 
    geom_hline(yintercept=seq(.5, 1, by=.1), lty='dotted', lwd=.5, col='grey75') +
    geom_bar(stat='identity', position='dodge') + 
    geom_errorbar(aes(ymin=mean-2 * se, ymax=mean + 2 * se), width=.2,
                  position=position_dodge(.9))  + 
    labs(x='Mark or Assay',y='AUROC on Observed Dataset') + 
    theme_pubr() + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    coord_cartesian(ylim=c(.5,1)) + 
    # scale_fill_manual(values=c('indianred','grey70','grey85'), name='Predictor') + 
    scale_fill_manual(values=c(col.paired[5],'grey85',col.paired[1]), name='Predictor') + 
    theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'auroc_comparison_bars.png'), g2, dpi=400, units='in', width=6, height=4)


garr = ggarrange(g1, g2, common.legend=TRUE)#   + rremove('y.text') + rremove('y.axis') + rremove('ylab') + rremove('y.ticks'))
ggsave(paste0(imgpref, 'ap_and_auroc_comparison_bars.png'), garr, dpi=400, units='in', width=10, height=4)



# -------------------------------------
# Load in the DHS recovery data + plot:
# -------------------------------------




# TODO: Remove 
# TODO: Plot scatter.
# TODO: Plot bar chart (fill hist) - how often each is best.






