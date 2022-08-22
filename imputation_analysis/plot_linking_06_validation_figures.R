#!/usr/bin/R
# -----------------------------------
# Plot validation metrics for linking
# -----------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
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
library(plotrix)
library(huxtable)
library(tidyr)
library(kableExtra)

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

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "linking/")
cmd = paste('mkdir -p', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'validation_')

# ---------------------
# Load in metrics data:
# ---------------------
valdir = 'linking_validation/'
f1files = list.files(pattern='f1.*.tsv', path=valdir)
cells = sub("f1_", "", sub(".tsv","", f1files))

ovdf = data.frame()
f1df = data.frame()
apdf = data.frame()
for (cell in cells){
    # F1
    df = read.delim(paste0(valdir, 'f1_', cell, '.tsv'), header=T)
    df$method = rownames(df)
    df = gather(df, metric, value, -method)
    df$cell = cell
    f1df = rbind(f1df, df)
    # AUPRC
    df = read.delim(paste0(valdir, 'auprc_', cell, '.tsv'), header=T)
    df$method = rownames(df)
    df = gather(df, metric, value, -method)
    df$cell = cell
    apdf = rbind(apdf, df)
    # OVERALL DISTANCE STATS:
    df = read.delim(paste0(valdir, 'ovr_distance_stats_', cell, '.tsv'), header=T)
    df$cell = cell
    ovdf = rbind(ovdf, df)
}


# Remap names of methods
mmap = data.frame(method=c("distance", "roadmap",
                           "H3K27ac_m_d1", "H3K27ac_c_d1", "H3K27ac_cm_d1",
                           "epimap-xgb", "epimap-xgb-mark",
                           "epimap-xgb-d50k", "epimap-xgb-mark-d50k"),
                  mname=c('Distance','Roadmap','MarkbyDist','CorrbyDist','CorrMarkbyDist',
                          'Epimap Corr','Epimap CorrMark','Epimap Corr + Dist', 'Epimap CorrMark + Dist'))
mmap$mname = factor(mmap$mname, levels=mmap$mname)

f1df$mname = NULL
apdf$mname = NULL
ovdf$mname = NULL
f1df = merge(f1df, mmap)
apdf = merge(apdf, mmap)
ovdf = merge(ovdf, mmap)

# -------------------
# Plot basic metrics:
# -------------------
col.paired = brewer.pal(12, 'Paired')
link.cols = c('grey', col.paired[4],
              colorRampPalette(col.paired[7:8])(3),
              colorRampPalette(col.paired[1:2])(4))

ggplot(f1df, aes(metric, value, fill=mname)) + 
    geom_bar(stat='identity', position='dodge', color='white') + 
    facet_wrap(~cell, scales='free_x', nrow=1) + 
    labs(x='Validation Link Dataset',y='F1 Score') + 
    scale_fill_manual(values=link.cols, name='Method') + 
    theme_pubr()


ggplot(f1df, aes(metric, value, color=mname)) + 
    # geom_bar(stat='identity', position='dodge', color='white') + 
    geom_point() + 
    coord_flip() + 
    facet_wrap(~cell, scales='free_x', nrow=1) + 
    labs(x='Validation Link Dataset',y='F1 Score') + 
    scale_color_manual(values=link.cols, name='Method') + 
    theme_pubr()


barplot(value ~ metric + mname + cell, f1df)

# Plot the BENGI gold standards first:
f1df = f1df[order(f1df$metric),]
sub.f1df = f1df[!(f1df$metric %in% c('GREAT.withnearest', 'GREAT.withoutnearest')),]
g1 = ggplot(sub.f1df, aes(cell, value, fill=mname)) + 
    geom_bar(stat='identity', position='dodge', color='white') + 
    facet_grid(~metric, space='free_x', scales='free_x') + 
    labs(x='',y='F1 Score') + 
    scale_y_continuous(expand=c(0,0), limits=c(0,1)) + 
    scale_fill_manual(values=link.cols, name='Method') + 
    theme_pubr()

apdf = apdf[order(apdf$metric),]
sub.apdf = apdf[!(apdf$metric %in% c('GREAT.withnearest', 'GREAT.withoutnearest')),]
g2 = ggplot(sub.apdf, aes(cell, value, fill=mname)) + 
    geom_bar(stat='identity', position='dodge', color='white') + 
    facet_grid(~metric, space='free_x', scales='free_x') + 
    labs(x='',y='AUPRC') + 
    scale_y_continuous(expand=c(0,0), limits=c(0,.75)) + 
    scale_fill_manual(values=link.cols, name='Method') + 
    theme_pubr()

garr = ggarrange(g1, g2, ncol=1, common.legend=TRUE)
ggsave(paste0(imgpref, 'auprc_f1_BENGI_datasets.png'), garr, dpi=450, units='in', width=14, height=7)
ggsave(paste0(imgpref, 'auprc_f1_BENGI_datasets.pdf'), garr, width=14, height=7)


# Plot the GREAT enrichment results:
f1df = f1df[order(f1df$metric),]
sub.f1df = f1df[f1df$metric %in% c('GREAT.withnearest', 'GREAT.withoutnearest'),]
sub.f1df$metric[sub.f1df$metric == 'GREAT.withnearest'] = 'GO (with nearest gene)'
sub.f1df$metric[sub.f1df$metric == 'GREAT.withoutnearest'] = 'GO (without nearest gene)'
g3 = ggplot(sub.f1df, aes(cell, value, fill=mname)) + 
    geom_bar(stat='identity', position='dodge', color='white') + 
    facet_grid(~metric, space='free_x', scales='free_x') + 
    # labs(x='Validation Link Dataset',y='F1 Score') + 
    labs(x='',y='F1 Score') + 
    scale_y_continuous(expand=c(0,0), limits=c(0,.33)) + 
    scale_fill_manual(values=link.cols, name='Method') + 
    theme_pubr()

apdf = apdf[order(apdf$metric),]
sub.apdf = apdf[apdf$metric %in% c('GREAT.withnearest', 'GREAT.withoutnearest'),]
sub.apdf$metric[sub.apdf$metric == 'GREAT.withnearest'] = 'GO (with nearest)'
sub.apdf$metric[sub.apdf$metric == 'GREAT.withoutnearest'] = 'GO (without nearest)'
g4 = ggplot(sub.apdf, aes(cell, value, fill=mname)) + 
    geom_bar(stat='identity', position='dodge', color='white') + 
    facet_grid(~metric, space='free_x', scales='free_x') + 
    # labs(x='Validation Link Dataset',y='AUPRC') + 
    labs(x='',y='AUPRC') + 
    scale_y_continuous(expand=c(0,0), limits=c(0,.33)) + 
    scale_fill_manual(values=link.cols, name='Method') + 
    theme_pubr()

garr = ggarrange(g3, g4, ncol=1, common.legend=TRUE)
ggsave(paste0(imgpref, 'auprc_f1_GREAT_benchmark.png'), garr, dpi=450, units='in', width=10, height=6)
ggsave(paste0(imgpref, 'auprc_f1_GREAT_benchmark.pdf'), garr, width=10, height=6)







# Distance f1 + auprc:
sub.ovdf = ovdf[ovdf$offset <= 1e6,]
ggplot(sub.ovdf, aes(offset, f1, color=mname)) + 
    geom_line() + 
    # geom_point() + 
    labs(x='Distance between enhancer center and - gene TSS', y='F1 Score') + 
    facet_grid(cell ~ dataset) + 
    scale_color_manual(values=link.cols, name='Method') + 
    theme_pubr()

# Distance f1 + auprc:
sub.ovdf = ovdf[ovdf$offset <= 1e6,]
ggplot(sub.ovdf, aes(offset, auprc, color=mname)) + 
    geom_line() + 
    # geom_point() + 
    labs(x='Distance between enhancer center and - gene TSS', y='AUPRC') + 
    facet_grid(cell ~ dataset) + 
    scale_color_manual(values=link.cols, name='Method') + 
    theme_pubr()



kbt = kable(f1df, 'latex', booktabs=T)
save_kable(kbt, file='~/test.pdf')


fwide = f1df
fwide$value = round(fwide$value, 2)
fwide = spread(fwide[,c('mname','cell','value','metric')], mname, value)
fwide = fwide[,-1]
colnames(fwide)[1] = 'Metric'
library(dplyr)
kbt = fwide %>% mutate_if(is.numeric, function(x) {
                              cell_spec(x, "latex", bold = T,
                                        background = spec_color(x, option = "A", direction = -1)) }) %>% 
    kable('latex', booktabs=T) %>% 
    pack_rows("CD34", 1, 3) %>%
    pack_rows("GM12878", 4, 11) %>%
    pack_rows("HeLa", 12, 16) %>%
    pack_rows("K562", 17, 20) %>%
    row_spec(c(1,4, 12,17)-1, extra_latex_after = "\\rowcolor{gray!6}") %>%
    collapse_rows(1, latex_hline = "none")
save_kable(kbt, file='~/test.pdf')





# Change with cell_spec:
cell_spec(mpg, "latex", color =ifelse(mpg > 20, "red", "blue")),
cell_spec(cyl, "latex", color = "white", align = "c", angle = 45,background =factor(cyl,c(4, 6, 8),c("#666666", "#999999", "#BBBBBB")))







