#!/usr/bin/R
# ----------------------------------------------------
# Plot the gwas hg stats for quality of enrich with DNase etc.
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

# Directories:
imgpref = paste0(img, "gwas_hg/")
system(paste('mkdir -p', imgpref))

# Load all metrics files:
stdir = 'gwas_hg_stats/'
fnames = list.files(path=stdir, pattern='.*.table.tsv')
fnames = fnames[grep('nonovl', fnames, invert=TRUE)]

fnames = list.files(path=stdir, pattern='.*merge2_raw.*.table.tsv')
fnames = c(fnames, list.files(path=stdir, pattern='.*100_raw.*.table.tsv'))

df = c()
for (fn in fnames){
    tab = read.delim(paste0(stdir, fn), sep="\t", header=T)
    elem = sub("_.*","", fn)
    tab$elem = elem
    if (length(grep("wH3K27ac", fn)) > 0){
        tab$merge = 'Annotation with H3K27ac'
    } else {
        tab$merge = 'Annotation alone'
    } 
    df = rbind(df, tab)
}

# Plot stats:
sdf = df[df$qt %in% c('99%'),]
sdf = sdf[sdf$ct.type == 'Single cutoff',]
# sdf = sdf[!(sdf$elem == 'H3K9ac'),]
sdf = sdf[sdf$metric %in% c('ngwas', 'ntotal'),]
subdf = sdf[sdf$metric == 'ngwas' & sdf$merge == 'Annotation with H3K27ac',]
eord = subdf$elem[order(subdf$count)]
sdf$elem = factor(sdf$elem, levels=eord)
sdf = sdf[!(sdf$elem == 'H3K27ac' & sdf$merge == 'Annotation with H3K27ac'),]
sdf$ann.type = 1
sdf$ann.type[sdf$elem == 'DNase-seq'] = 2
sdf$ann.type[sdf$elem == 'ENH'] = 3

col.paired = brewer.pal(12, 'Paired')

gplot = ggplot(sdf, aes(elem, count)) +
    geom_bar(stat='identity', fill='grey75') + 
    geom_text(aes(elem, count + space, label=count)) + 
    theme_pubr() + 
    labs(x='Element Type',y='GWAS Enrichment Metric') + 
    facet_grid(metric~merge, scales='free')  + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) 


g1 = ggplot(sdf[sdf$metric == 'ngwas',], aes(elem, count, fill=factor(ann.type))) +
    facet_grid(~merge, scales='free')  + 
    geom_bar(stat='identity') + 
    geom_text(aes(elem, count + space, label=count, color=factor(ann.type))) + 
    scale_y_continuous(expand=c(0,0))+ 
    scale_fill_manual(values=c('grey75',col.paired[7], col.paired[8])) + 
    scale_color_manual(values=c('grey50',col.paired[7], col.paired[8])) + 
    theme_pubr() + 
    labs(x='DHS annotation used to test for GWAS enrichments',y='Number of enriched\nGWAS (of 803 total)') + 
    # theme(axis.text.x = element_text(angle=45, hjust=1)) 
    theme(axis.text.x = element_text(angle=45, hjust=1))  + 
    theme(legend.position='none')

ggsave(paste0(imgpref, 'gwas_hg_stats_ngwas.pdf'), g1, dpi=300, width=5, height=3,units='in')
ggsave(paste0(imgpref, 'gwas_hg_stats_ngwas.png'), g1, dpi=300, width=5, height=3,units='in')


g2 = ggplot(sdf[sdf$metric == 'ntotal',], aes(elem, count, fill=factor(ann.type))) +
    facet_grid(~merge, scales='free')  + 
    geom_bar(stat='identity') + 
    geom_text(aes(elem, count + space, label=count, color=factor(ann.type))) + 
    scale_y_continuous(expand=c(0,0))+ 
    scale_fill_manual(values=c('grey75',col.paired[7], col.paired[8])) + 
    scale_color_manual(values=c('grey50',col.paired[7], col.paired[8])) + 
    theme_pubr() + 
    labs(x='DHS annotation used to test for GWAS enrichments',y='Number of enrichments\n(of 669k total)') + 
    # theme(axis.text.x = element_text(angle=45, hjust=1)) 
    theme(axis.text.x = element_text(angle=45, hjust=1))  + 
    theme(legend.position='none')


ggsave(paste0(imgpref, 'gwas_hg_stats_ntotal.pdf'), g2, dpi=300, width=5, height=3,units='in')
ggsave(paste0(imgpref, 'gwas_hg_stats_ntotal.png'), g2, dpi=300, width=5, height=3,units='in')

garr = ggarrange(g1, g2, ncol=1)

ggsave(paste0(imgpref, 'gwas_hg_stats.pdf'), garr, dpi=300, width=6, height=6,units='in')
ggsave(paste0(imgpref, 'gwas_hg_stats.png'), garr, dpi=300, width=6, height=6,units='in')


