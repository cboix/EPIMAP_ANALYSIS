#!/usr/bin/R
library(tidyr) 
library(ggplot2)

# Read in data, reshape to long format:
# -------------------------------------
tab = read.delim('pval_distribution_092721.tsv', header=T)
names(tab)[1] = 'log10p'
df = gather(tab, 'track','nbin', -log10p)
df$mark_index = sub(".*M","M", df$track)
df$cell_index = sub("M.*","", df$track)

# Merge in cell / mark
# --------------------
cmap = read.delim('IC_cell_mapping.tsv', header=T, sep=" ")
mmap = read.delim('IC_mark_mapping.tsv', header=T, sep=" ")
df = merge(merge(df, cmap), mmap)

# Plot on log-log scale:
# ----------------------
gplot = ggplot(df, aes(log10p, nbin, color=cell_type)) + 
    facet_wrap(~mark) + 
    geom_hline(yintercept=1184, lty='dashed') + 
    scale_x_log10() + 
    scale_y_log10() + 
    geom_line() + theme_bw()
       
# gplot

ggsave('pval_distribution_log_log_plot.png', gplot, dpi=300, units='in', width=10, height=8)


