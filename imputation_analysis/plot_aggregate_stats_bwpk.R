#!/usr/bin/R
# -----------------------------------------------
# Given intersection from BW and PK file:
# 1. Collect aggregate statistics
# 2. Plot each epigenome
# 3. Save as matrix
# -----------------------------------------------
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
library(ggpubr)
library(scales)

options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")

# Args
filename='all_bwpk_statistics.tsv'
tagline='H3K27ac in DNase-seq peaks'
plotprefix='all_bwpk_statistics'
infofile='all_pkbw_statistics_infofile.tsv'

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied. Need:
          filename tagline plotprefix infofile")
} else {        
    filename = args[1]
    tagline = args[2]
    plotprefix = args[3]
    infofile = args[4]
    # TODO: DETECT OR FEED IN EXTENSION (pdf, png?)
}

# Read in intersection:
df = read.delim(filename, header=T, sep="\t")
metrics = as.character(unique(df$metric))
mmap = c('frac2' = '% Peaks with H3K27ac > 2', 
        'gt2' = '# Peaks with H3K27ac > 2',
        'gt3' = '# Peaks with H3K27ac > 3',
        'gt4' = '# Peaks with H3K27ac > 4',
        'gt5' = '# Peaks with H3K27ac > 5',
        'tot.dhs' = 'Total # DNase-seq peaks')
df = merge(df, data.frame(ext.metric=mmap, metric=names(mmap)))
keptid = scan('Annotation/kept_bssid_20190322.txt','c')
df = df[df$id %in% keptid,]

info = read.delim(infofile, header=F, sep="\t")
names(info) = c('id','pk','bw')
info$pkimp = "Observed Peaks"
info$pkimp[grep("impute", info$pk)] = "Imputed Peaks"
info$bwimp = "Observed H3K27ac"
info$bwimp[grep("impute", info$bw)] = "Imputed H3K27ac"

df = merge(df, info)

# --------------------
# Plot all statistics:
# --------------------
# All statistics - pdf: 
gp = ggplot(df, aes(value)) +
    geom_histogram(color='white', fill='darkgrey', bins=30) + 
    theme_pubclean() + 
    scale_x_continuous(labels = scales::comma) +
    labs(x='Number or Percent of Peaks', 
         y='Number of Epigenomes', 
         title=paste('Statistics for', tagline)) +
    facet_wrap(~ext.metric, scales='free')
ggsave(paste0(plotprefix, '_allstat_pdf.png'), gp, dpi=300, width=10, height=6, units='in')


# All statistics as cdf: 
gp = ggplot(df, aes(value)) +
    stat_ecdf(aes(ymin=0,ymax=..y..),geom="ribbon", fill='lightgrey', alpha=0.5) +
    stat_ecdf(geom="step") +
    theme_pubclean() + 
    scale_x_continuous(labels = scales::comma) +
    labs(x='Number or Percent of Peaks', 
         y='Cumulative % of Epigenomes', 
         title=paste('Statistics for', tagline)) +
         facet_wrap(~ext.metric, scales='free')
ggsave(paste0(plotprefix, '_allstat_cdf.png'), gp, dpi=300, width=10, height=6, units='in')


# ----------------------
# Plot count statistics:
# ----------------------
# Stats counts only:
ctmetrics = c('gt5','gt4','gt3','gt2','tot.dhs')
cdf = df[df$metric %in% ctmetrics,]
cvals = col[round(seq(25, 100, length.out = 5))]
names(cvals) = mmap[ctmetrics]
cdf$ext.metric = factor(as.character(cdf$ext.metric), levels=rev(mmap[ctmetrics]))

gp = ggplot(cdf, aes(value, fill=ext.metric)) +
    geom_histogram(color='white', bins=30) + 
    theme_pubclean() + 
    scale_fill_manual(values = cvals, name='Metric') + 
    scale_x_continuous(labels = scales::comma) +
    theme(legend.position='right') + 
    labs(x='Number of Peaks', 
         y='Number of Epigenomes', 
         title=paste('Count Statistics for', tagline)) +
    facet_wrap(~ext.metric, scales='free')
ggsave(paste0(plotprefix, '_countstat_pdf.png'), gp, dpi=300, width=12.5, height=6, units='in')


gp = ggplot(cdf, aes(value, fill=ext.metric)) +
    geom_histogram(color='white', bins=30) + 
    theme_pubclean() + 
    scale_fill_manual(values = cvals, name='Metric') + 
    scale_x_continuous(labels = scales::comma) +
    theme(legend.position='right') + 
    labs(x='Number of Peaks', 
         y='Number of Epigenomes', 
         title=paste('Count Statistics for', tagline)) +
    facet_grid(pkimp~ext.metric, scales='free')
ggsave(paste0(plotprefix, '_countstat_pkimp_pdf.png'), gp, dpi=300, width=15, height=6, units='in')


gp = ggplot(cdf, aes(value, fill=ext.metric)) +
    geom_histogram(color='white', bins=30) + 
    theme_pubclean() + 
    scale_fill_manual(values = cvals, name='Metric') + 
    scale_x_continuous(labels = scales::comma) +
    theme(legend.position='right') + 
    labs(x='Number of Peaks', 
         y='Number of Epigenomes', 
         title=paste('Count Statistics for', tagline)) +
    facet_grid(bwimp~ext.metric, scales='free')
ggsave(paste0(plotprefix, '_countstat_bwimp_pdf.png'), gp, dpi=300, width=15, height=6, units='in')


gp = ggplot(cdf, aes(value, fill=ext.metric)) +
    geom_histogram(color='white', bins=30) + 
    theme_pubclean() + 
    scale_fill_manual(values = cvals, name='Metric') + 
    scale_x_continuous(labels = scales::comma) +
    theme(legend.position='right') + 
    labs(x='Number of Peaks', 
         y='Number of Epigenomes', 
         title=paste('Count Statistics for', tagline)) +
    facet_grid(bwimp +pkimp~ext.metric, scales='free')
ggsave(paste0(plotprefix, '_countstat_bothimp_pdf.png'), gp, dpi=300, width=18, height=6, units='in')


# All counts as cdf: 
gp = ggplot(cdf, aes(value, color=ext.metric)) +
    stat_ecdf(geom="step", pad=FALSE) +
    theme_pubclean() + 
    scale_color_manual(values=cvals, name='Metric') + 
    scale_x_continuous(labels = scales::comma) +
    theme(legend.position='right') + 
    labs(x='Number or Percent of Peaks', 
         y='Cumulative % of Epigenomes', 
         title=paste('Count Statistics for', tagline))
ggsave(paste0(plotprefix, '_countstat_cdf.png'), gp, dpi=300, width=8, height=5, units='in')

# ----------------------
# Plot frac statistics:
# ----------------------
# Get the fracs for others + plot:
ctmetrics = c('gt5','gt4','gt3','gt2')
cdf = df[df$metric %in% ctmetrics,]
tdf =df[df$metric == 'tot.dhs',c('value','id')]
names(tdf)[1] = 'total'
cdf = merge(cdf, tdf)
cdf$frac.metric = sub("# ", "% ", as.character(cdf$ext.metric))
lvs =sub("#","%", mmap[ctmetrics])
cdf$frac.metric = factor(as.character(cdf$frac.metric), levels=rev(lvs))
cvals = col3[round(seq(25, 100, length.out = length(ctmetrics)))]
names(cvals) = lvs

gp = ggplot(cdf, aes(value/total, fill=frac.metric)) +
    geom_histogram(color='white', bins=30) + 
    theme_pubclean() + 
    scale_fill_manual(values = cvals, name='Metric') + 
    scale_x_continuous(labels = scales::percent, limits=c(0,1)) +
    theme(legend.position='right') + 
    labs(x='Percent of Peaks', 
         y='Number of Epigenomes', 
         title=paste('Percent Statistics for', tagline)) +
    facet_wrap(~frac.metric, scales='free')
ggsave(paste0(plotprefix, '_fracstat_pdf.png'), gp, dpi=300, width=10, height=6, units='in')


gp = ggplot(cdf, aes(value/total, fill=frac.metric)) +
    geom_histogram(color='white', bins=30) + 
    theme_pubclean() + 
    scale_fill_manual(values = cvals, name='Metric') + 
    scale_x_continuous(labels = scales::percent, limits=c(0,1)) +
    theme(legend.position='right') + 
    labs(x='Percent of Peaks', 
         y='Number of Epigenomes', 
         title=paste('Percent Statistics for', tagline)) +
    facet_grid(bwimp+pkimp ~ frac.metric, scales='free')
ggsave(paste0(plotprefix, '_fracstat_bothimp_pdf.png'), gp, dpi=300, width=15, height=6, units='in')

# --------------------------------
# Plot tot.dhs against count vals:
# --------------------------------
gp = ggplot(cdf, aes(value, total, color=frac.metric)) +
    geom_point(alpha=0.5) + 
    theme_pubclean() + 
    scale_color_manual(values = cvals, name='Metric') + 
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    theme(legend.position='right') + 
    labs(x='Number of Kept Peaks', 
         y='Number of Total Peaks)', 
         title=paste('Peak Inclusion for', tagline)) +
    facet_wrap(~frac.metric, scales='free')
ggsave(paste0(plotprefix, '_countstat_scatter.png'), gp, dpi=300, width=10, height=6, units='in')

gp = ggplot(cdf, aes(value, total, color=frac.metric)) +
    geom_point(alpha=0.5) + 
    theme_pubclean() + 
    scale_color_manual(values = cvals, name='Metric') + 
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    theme(legend.position='right') + 
    labs(x='Number of Kept Peaks', 
         y='Number of Total Peaks)', 
         title=paste('Peak Inclusion for', tagline)) +
    facet_grid(bwimp + pkimp~frac.metric, scales='free')
ggsave(paste0(plotprefix, '_countstat_bothimp_scatter.png'), gp, dpi=300, width=18, height=6, units='in')

print("Plotting finished.")
