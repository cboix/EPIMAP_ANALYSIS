#!/usr/bin/R
# -----------------------------------------------------------
# Plot the imputed agreement metrics for chr1 on the IC data:
# -----------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggbeeswarm)

library(viridis)
library(gridBase)
library(tidyr)
options(scipen=45)

imgdir = paste0(img, 'validation/')
system(paste0('mkdir -p ', imgdir))
imgpref = paste0(imgdir, 'external_')

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

# -----------------
# Read in the data:
# -----------------
evdir = 'external_validation/'
evfile = paste0(evdir, 'allic.chr1.stats.tsv')
df = read.delim(evfile)
df$ic.id = sub("M.*", "", df$prefix)
df$ic.mark = sub(".*M", "M", df$prefix)

# Read in mappings:
nmap = read.delim('Annotation/ic_namespace_mapping.tsv',header=T)
mmap = read.delim('Annotation/ic_mark_mapping.tsv',header=F)
names(mmap) = c('ic.mark', 'mark')
df = merge(merge(df, mmap, all.x=TRUE), nmap, all.x=TRUE)
metrics = c( 'gwcorr', 'auroc', 'auprc', 'auroc1', 'auprc1', 'catchobs', 'catchimp', 'matchobs', 'matchimp')

# Split + reduce to best sample, (by each metric):
idf = df[df$comp.id == df$prefix,c(metrics, 'prefix','mark', 'step')]
mdf = df[df$comp.id == 'Mean',c(metrics, 'prefix','mark', 'step')]
odf = df[df$comp.id != 'Mean' & df$comp.id != df$prefix,c(metrics, 'prefix','mark', 'step')]

o2df = aggregate(cbind(gwcorr, auroc, auprc, auroc1, auprc1, catchobs, catchimp, matchobs, matchimp) ~ step + mark + prefix, odf, max)
olong = gather(o2df[, c(metrics, 'mark','prefix', 'step')], metric, value, -prefix, -mark, -step)
mlong = gather(mdf[, c(metrics, 'mark','prefix', 'step')], metric, value, -prefix, -mark, -step)
ilong = gather(idf[, c(metrics, 'mark','prefix', 'step')], metric, value, -prefix, -mark, -step)
olong$type = 'Nearest'
mlong$type = 'Mean'
ilong$type = 'Imputed'
all.long = rbind(mlong, rbind(olong, ilong))
all.long = merge(all.long, unique(df[,c('prefix','ic.infoline','ic.name')]))
slong = all.long[all.long$step == 200, c('metric', 'value', 'type', 'mark','prefix')]

# How often does imputed beat the others:
tab = spread(slong, type, value)
tab$im = tab$Imputed > tab$Mean
tab$inr = tab$Imputed > tab$Nearest
tab$imnr = tab$Imputed > tab$Nearest & tab$Imputed > tab$Mean

unique(aggregate(im ~ metric + mark, tab, length)[, c('im','mark')])
colnames(markdef)[2] = 'mark.type'
tab = merge(tab, markdef, all.x=TRUE)
t1df = aggregate(cbind(im, inr, imnr) ~ metric, tab, mean)
t2df = aggregate(cbind(im, inr, imnr) ~ metric + mark, tab, mean)
t3df = aggregate(cbind(im, inr, imnr) ~ metric + mark.type, tab, mean)
t2df = gather(t2df, comparison, value, -mark, -metric)
# merge(t2df$
# t2df$comparison

aa = unique(tab[,c('prefix','mark.type')])

aggregate(cbind(im, inr, imnr) ~ metric, tab[tab$metric == 'auprc1' & tab$mark %in% c('H3K27me3','H3K9me3'),], mean)
aggregate(cbind(im, inr, imnr) ~ metric, tab[tab$metric == 'auprc1' & tab$mark %in% c('H3K27me3','H3K9me3'),], sum)

gplot = ggplot(t2df, aes(mark, value, fill=comparison)) + 
    facet_wrap(~metric) + 
    geom_hline(yintercept=seq(.2, .8, by=.2), lty='dotted', lwd=.5, col='grey75') +
    geom_bar(stat='identity', position='dodge') + 
    labs(x='Mark or Assay',y='Average Precision on Observed Dataset') + 
    theme_pubr() + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    scale_fill_manual(values=c('grey95','grey85','grey75'), name='Predictor') + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(paste0(imgpref, 'allmetrics_permark_200bp_frac.png'), gplot, dpi=400, units='in', width=12, height=8)


# --------------------------------------
# Plot as bar charts, both AP and AUROC:
# --------------------------------------
mdf = agg.rename(value ~ type + mark + metric, slong, mean, 'mean')
sdf = agg.rename(value ~ type + mark + metric, slong, sd, 'sd')
ndf = agg.rename(value ~ type + mark + metric, slong, length, 'n')

mdf = merge(merge(mdf, sdf), ndf)
mdf$se = mdf$sd / sqrt(mdf$n)
mdf$mark = factor(mdf$mark, levels=t12marks)

col.paired = brewer.pal(12, 'Paired')

gplot = ggplot(mdf, aes(mark, mean, fill=type)) + 
    facet_wrap(~metric) + 
    geom_hline(yintercept=seq(.2, .8, by=.2), lty='dotted', lwd=.5, col='grey75') +
    geom_bar(stat='identity', position='dodge') + 
    geom_errorbar(aes(ymin=mean-2 * se, ymax=mean + 2 * se), width=.2,
                  position=position_dodge(.9))  + 
    labs(x='Mark or Assay',y='Average Precision on Observed Dataset') + 
    theme_pubr() + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    scale_fill_manual(values=c(col.paired[5],'grey85',col.paired[1]), name='Predictor') + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(paste0(imgpref, 'allmetrics_permark_200bp.png'), gplot, dpi=400, units='in', width=12, height=8)

# Plot reduced:
# metric.map = data.frame(metric=c('auroc1','auprc1','gwcorr'),
#                         Metric=c('AUROC','AUPRC','Correlation'))
metric.map = data.frame(metric=c('auroc1','auprc1','catchobs','catchimp'),
                        Metric=c('AUROC','AUPRC','Catch1Obs','Catch1Imp'))
mdf2 = merge(mdf, metric.map)
mdf2$mark = factor(mdf2$mark, levels=c('H3K27ac','H3K4me1','H3K4me3','ATAC-seq','DNase-seq','H3K36me3','H3K27me3','H3K9me3'))

gplot = ggplot(mdf2, aes(mark, mean, fill=type)) + 
    facet_wrap(~Metric) + 
    geom_hline(yintercept=seq(.25, 1, by=.25), lty='dotted', lwd=.5, col='grey75') +
    geom_bar(stat='identity', position='dodge') + 
    geom_errorbar(aes(ymin=mean-2 * se, ymax=mean + 2 * se), width=.2,
                  position=position_dodge(.9))  + 
    labs(x='Mark or Assay',y='AUC or % Recovery') + 
    theme_pubr() + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    scale_fill_manual(values=c(col.paired[5],'grey85',col.paired[1]), name='Predictor') + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))
ggsave(paste0(imgpref, 'mainmetrics_permark_200bp.png'), gplot, dpi=400, units='in', width=6.5, height=5)
ggsave(paste0(imgpref, 'mainmetrics_permark_200bp.pdf'), gplot, width=6.5, height=5)


mwide = spread(mdf2[,c('Metric','type','mean', 'mark')], type, mean)
mwide$ImpNear = mwide$Imputed / mwide$Nearest
mwide$ImpMean = mwide$Imputed / mwide$Mean

# Statistics for the text:
mwide[mwide$Metric == 'AUPRC',]
t3df = t2df[t2df$metric == 'auprc1',]
t3df$value = paste0(round(100 * t3df$value,1), '%')
t3df = spread(t3df, comparison, value)




# See where the H3K4me3 is beaten in auroc:
swide = spread(slong, type, value)
subwide = swide[swide$mark == 'H3K4me3' & swide$Imputed < swide$Nearest & swide$metric == 'auroc',]

all.long[all.long$prefix %in% subwide$prefix & all.long$metric == 'auroc',]

df[df$auroc > 0.8 & df$prefix %in% subwide$prefix,]

df[df$auroc > 0.8 & df$prefix == 'C40M22',]
tdf = df[df$prefix == 'C40M22',]
idf[idf$prefix %in% subwide$prefix,]



# ----------------
# Plot as scatter:
# ----------------
swide = spread(slong, type, value)
sub.wide = swide[swide$metric == 'auprc1',]
sub.wide = merge(sub.wide, markdef)
sub.wide = merge(sub.wide, unique(df[, c('prefix','mark','ic.name')]))
indp = sub.wide$mark.type == 'punctate'
indb = sub.wide$mark.type == 'broad'

pct.over1p = mean(sub.wide$Imputed[indp] >= sub.wide$Nearest[indp])
pct.over1b = mean(sub.wide$Imputed[indb] >= sub.wide$Nearest[indb])
labdf = data.frame(rbind(c(x=.15, y=.9, label=paste0(round(100 * pct.over1p, 1), '%'), mark.type='punctate'),
                 c(x=.8, y=.15, label=paste0(round(100 * (1 - pct.over1p) , 1), '%'), mark.type='punctate'),
                 c(x=.15, y=.8, label=paste0(round(100 * pct.over1b, 1), '%'), mark.type='broad'),
                 c(x=.8, y=.05, label=paste0(round(100 * (1 - pct.over1b) , 1), '%'), mark.type='broad')))
labdf$x = as.numeric(labdf$x)
labdf$y = as.numeric(labdf$y)
g1 = ggplot(sub.wide, aes(Nearest, Imputed, color=mark.type)) + 
    geom_abline(intercept=0,slope=1) + 
    geom_point(cex=.75) + 
    geom_text(data=labdf,aes(x=x, y=y, label=label, col=mark.type), show.legend = FALSE) +
    scale_color_manual(values=col.paired[c(6,2)], name='Assay type:') + 
    labs(x='AP of Nearest Observed Dataset',y='AP of Imputed Track') + 
    theme_pubr() # + theme(legend.position='none')


pct.over2p = mean(sub.wide$Imputed[indp] >= sub.wide$Mean[indp])
pct.over2b = mean(sub.wide$Imputed[indb] >= sub.wide$Mean[indb])
labdf = data.frame(rbind(c(x=.15, y=.9, label=paste0(round(100 * pct.over2p, 1), '%'), mark.type='punctate'),
                 c(x=.8, y=.15, label=paste0(round(100 * (1 - pct.over2p) , 1), '%'), mark.type='punctate'),
                 c(x=.15, y=.8, label=paste0(round(100 * pct.over2b, 1), '%'), mark.type='broad'),
                 c(x=.8, y=.05, label=paste0(round(100 * (1 - pct.over2b) , 1), '%'), mark.type='broad')))
labdf$x = as.numeric(labdf$x)
labdf$y = as.numeric(labdf$y)
g2 = ggplot(sub.wide, aes(Mean, Imputed, color=mark.type)) + 
    geom_abline(intercept=0,slope=1) + 
    geom_point(cex=.75) + 
    geom_text(data=labdf,aes(x=x, y=y, label=label, col=mark.type), show.legend = FALSE) +
    scale_color_manual(values=col.paired[c(6,2)], name='Assay type:') + 
    labs(x='AP of Mean Observed Data',y='AP of Imputed Track') + 
    theme_pubr() # + theme(legend.position='none')

garr = ggarrange(g1, g2, common.legend=TRUE)#   + rremove('y.text') + rremove('y.axis') + rremove('ylab') + rremove('y.ticks'))
ggsave(paste0(imgpref, 'validation_auprc_scatter_comparison.png'), garr, dpi=400, units='in', width=6, height=3)

mean((sub.wide$Imputed / sub.wide$Mean)[indp])
mean((sub.wide$Imputed / sub.wide$Mean)[indb])
mean((sub.wide$Imputed / sub.wide$Nearest)[indp])
mean((sub.wide$Imputed / sub.wide$Nearest)[indb])




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


plot(sub.wide$Nearest, sub.wide$Imputed)
abline(0,1)




# -------------
# Availability:
# -------------
amarg = colSums(avail)
wmarg = apply(avail, 2, function(x){paste0(names(which(x > 0)), collapse=', ')})
df$nmarg = amarg[df$id]
df$wmarg = wmarg[df$id]

# TODO: Image of tested data + blind data (for reviews)
ids = sort(unique(df$id))
smat = avail[,ids]
# Add in the blind data:


# Gather and plot all of the different metrics:
mlong = gather(df, metric, value, -ic.id, -ic.mark, -prefix, -mark, -step, -id, -ic.name, -ic.infoline,-nmarg)
ggplot(mlong[mlong$step == 200,], aes(mark, value)) + 
    facet_wrap(~metric, scales='free_y')+
    geom_boxplot(color='grey50') + 
    geom_quasirandom(cex=.5) + 
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'metrics_permark_200bp.png'), dpi=400, units='in', width=8, height=6)


metric.map = data.frame(metric=c('auroc','auprc','catchobs','catchimp'),
                        Metric=c('AUROC','AUPRC','Catch1Obs','Catch1Imp'))
m2long = merge(mlong, metric.map)

gplot = ggplot(m2long, aes(mark, value)) + 
    facet_wrap(~Metric, nrow=1) + # , scales='free', nrow=1)+
    scale_y_continuous(labels=scales::percent) + 
    geom_boxplot(fill='grey75') + 
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'redmetrics_permark_200bp.png'), gplot, dpi=400, units='in', width=8, height=4)



ggplot(mlong[mlong$step == 200,], aes(mark, value, color=ic.infoline)) + 
    facet_wrap(~metric, scales='free')+
    geom_boxplot(color='grey50') + 
    geom_quasirandom() + 
    theme_pubr()


ggplot(mlong, aes(mark, value, color=ic.infoline, fill=factor(step))) + 
    facet_wrap(~metric, scales='free')+
    geom_boxplot(color='grey50') + 
    geom_quasirandom() + 
    theme_pubr()



# Against




# TODO: See how much each is based on.




