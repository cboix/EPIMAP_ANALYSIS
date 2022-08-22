#!/usr/bin/R
# ------------------------------------
# Compare DHSs to screen dELs and pELs
# ------------------------------------
gtargs=(commandArgs(TRUE))
print(gtargs)
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
    plotting.trees=TRUE
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
    plotting.trees=TRUE
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))

library(dplyr)
library(cba)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(scales)
library(huxtable)
options(scipen=45) # So we dont get issues writing integers into bedfiles
options(repos='http://cran.rstudio.com/')

parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}

agg.rename = function(formula, data, FUN, name){
    agg.df = aggregate(formula, data, FUN)
    names(agg.df)[ncol(agg.df)] = name
    return(agg.df) }

today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "validation/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'screen_comparison_')

# ----------------------------------------
# Load in dELS, pELS, and Index locations:
# ----------------------------------------
ddf = read.delim('DHS_Index_WM201902/GRCh38-ccREs.dELS.bed', header=F)
pdf = read.delim('DHS_Index_WM201902/GRCh38-ccREs.pELS.bed', header=F)
adf = read.delim('DHS_Index_WM201902/GRCh38-ccREs.bed', header=F)
idf = read.delim('DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords.txt', header=F)
names(ddf) = c('chr','start','end','id1','id2','type')
names(pdf) = c('chr','start','end','id1','id2','type')
names(adf) = c('chr','start','end','id1','id2','type')
names(idf) = c('chr','start','end','name')

# Size:
cat(nrow(ddf), nrow(pdf), nrow(idf),'\n')

# Make ranges:
dr = GRanges(seqnames=ddf$chr, IRanges(start=ddf$start, end=ddf$end), id=ddf$id1)
pr = GRanges(seqnames=pdf$chr, IRanges(start=pdf$start, end=pdf$end), id=pdf$id1)
ar = GRanges(seqnames=adf$chr, IRanges(start=adf$start, end=adf$end), id=adf$id1, type=adf$type)
ir = GRanges(seqnames=idf$chr, IRanges(start=idf$start, end=idf$end), id=idf$name)

# Map all overlaps:
dovl = as.data.frame(findOverlaps(ir, dr, minoverlap=100))
povl = as.data.frame(findOverlaps(ir, pr, minoverlap=100))
aovl = as.data.frame(findOverlaps(ir, ar, minoverlap=100))
dovl$name = ir$id[dovl$queryHits]
povl$name = ir$id[povl$queryHits]
aovl$name = ir$id[aovl$queryHits]
aovl$type = ar$type[aovl$subjectHits]
# Matched 97%-98% of the locations to a DHS:
round(100 * length(unique(dovl$subjectHits)) / nrow(ddf),2) 
round(100 * length(unique(povl$subjectHits)) / nrow(pdf),2)
round(100 * length(unique(aovl$subjectHits)) / nrow(adf),2)


# ------------------------------------
# Load in the DHS list/enhancers data:
# ------------------------------------
ddir = 'DHS_Index_WM201902/'
dpref = 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19'
dmlfile = paste0(ddir, dpref, '.core.srt.txt')
dmlnamfile = paste0(ddir, dpref, '_r200_e0_names.core.srt.tsv')
dmlrdafile = paste0(ddir, dpref, '_enh.core.srt.Rda')
# enhind = as.numeric(scan('Enhancer_indices.txt', 'c')) + 1 # Old
mdir = 'masterlist_matindices/'
mpref = paste0(mdir, 'matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_')
enhind = as.numeric(scan(paste0(mpref, 'ENH','_masterlist_indices.tsv'),'c')) + 1
promind = as.numeric(scan(paste0(mpref, 'PROM','_masterlist_indices.tsv'),'c')) + 1
dyadind = as.numeric(scan(paste0(mpref, 'DYADIC','_masterlist_indices.tsv'),'c')) + 1
load(dmlrdafile)
dmldf = read.table(dmlfile, header=F, stringsAsFactors=F, sep="\t")
names(dmldf) = c('chr','start','end','name')

# Overall overlaps:
enhnam = dmldf$name[enhind]
promnam = dmldf$name[promind]
dyadnam = dmldf$name[dyadind]
dovl$in.enh = 1 * (dovl$name %in% enhnam)
povl$in.enh = 1 * (povl$name %in% enhnam)
povl$in.prom = 1 * (povl$name %in% promnam)
aovl$in.enh = 1 * (aovl$name %in% enhnam)
aovl$in.prom = 1 * (aovl$name %in% promnam)
aovl$in.dyad = 1 * (aovl$name %in% dyadnam)
dagg = aggregate(in.enh ~ subjectHits, dovl, max)
pagg = aggregate(in.enh ~ subjectHits, povl, max)
aagg = aggregate(in.enh ~ subjectHits + type, aovl, max)

paste0('dELS in enh: ', round(100 * mean(dagg$in.enh), 2),'% of ', nrow(dagg), ' mapped elements')
paste0('pELS in enh: ', round(100 * mean(pagg$in.enh), 2),'% of ', nrow(pagg), ' mapped elements')
paste0('all in enh: ', round(100 * mean(aagg$in.enh), 2),'% of ', nrow(aagg), ' mapped elements')
# dELS in enh: 84.53% of 654,481 mapped elements
# pELS in enh: 40.97% of 137,984 mapped elements
# all in enh: 71.28% of 906,360 mapped elements"
asum = aggregate(in.enh ~ type, aagg, mean)
asum$frac = paste0(round(100 * asum$in.enh, 2),'%')
asum = asum[order(asum$in.enh, decreasing=T),]
#                      dELS 0.8519808  85.2%
#           dELS,CTCF-bound 0.8316043 83.16%
#                      pELS 0.5576247 55.76%
#  DNase-H3K4me3,CTCF-bound 0.4833294 48.33%
#                       PLS 0.4475796 44.76%
#             DNase-H3K4me3 0.3735147 37.35%
#      CTCF-only,CTCF-bound 0.3638384 36.38%
#           pELS,CTCF-bound 0.2859823  28.6%
#            PLS,CTCF-bound 0.1036407 10.36%
agg = aggregate(cbind(in.enh, in.prom, in.dyad) ~ subjectHits + type, aovl, sum)

# Take for each the top:
agg$set = 'NA'
agg$sum = agg$in.enh + agg$in.prom + agg$in.dyad
agg$set[(agg$in.enh > agg$in.prom) & (agg$in.enh > agg$in.dyad)] = 'Enhancer'
agg$set[agg$in.prom > agg$in.enh & agg$in.prom > agg$in.dyad] = 'Promoter'
agg$set[agg$in.dyad > agg$in.enh & agg$in.dyad > agg$in.prom] = 'Dyadic'
agg$set[(agg$in.enh == agg$in.dyad) & agg$in.enh > 0] = 'Enhancer'
agg$set[(agg$in.prom == agg$in.dyad) & agg$in.prom > 0] = 'Promoter'
agg$set[(agg$in.enh == agg$in.prom) & agg$in.enh > 0] = 'Enhancer'
agg$set[agg$sum == 0] = 'None'

totsum = aggregate(subjectHits ~ type + set, agg, length)
totsum$set = factor(totsum$set, levels=rev(c('Promoter','Enhancer','Dyadic','None')))
totsum$main = sapply(totsum$type, function(x){strsplit(x,",")[[1]][1]})
tt = aggregate(subjectHits ~ set + main, totsum, sum)
twide = spread(tt, set, subjectHits)
tmat = as.matrix(twide[,-1])
rownames(tmat) = twide[,1]
tmat = tmat / rowSums(tmat)
round(tmat, 4)

gp = ggplot(totsum, aes(type, subjectHits, fill=set)) + 
    geom_bar(stat='identity', position='fill', width=.98) + 
    scale_fill_manual(values=c('Promoter'=rgb(255,0,0, max=255), 'Enhancer'=rgb(255,195,77, max=255),
                               'Dyadic'='royalblue', 'None'='grey75'), name='Epimap:') + 
    coord_flip() + 
    labs(y='Percent of SCREEN elements recovered', x='SCREEN element set') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    theme_pubr()
ggsave(paste0(imgpref, 'element_recovery_pct.png'), gp, dpi=450, units='in', width=7, height=4)
ggsave(paste0(imgpref, 'element_recovery_pct.pdf'), gp, dpi=450, units='in', width=7, height=4)

gp2 = ggplot(totsum, aes(type, subjectHits, fill=set)) + 
    geom_bar(stat='identity', position='stack', width=.98) + 
    scale_fill_manual(values=c('Promoter'=rgb(255,0,0, max=255), 'Enhancer'=rgb(255,195,77, max=255),
                               'Dyadic'='royalblue', 'None'='grey75'), name='Epimap:') + 
    coord_flip() + 
    labs(y='Number of SCREEN elements recovered', x='SCREEN element set') + 
    scale_y_continuous(labels=scales::comma, expand=c(0,0)) + 
    theme_pubr()
ggsave(paste0(imgpref, 'element_recovery_num.png'), gp2, dpi=450, units='in', width=7, height=4)
ggsave(paste0(imgpref, 'element_recovery_num.pdf'), gp2, dpi=450, units='in', width=7, height=4)


# Reverse set intersection:
NE = nrow(enhdf)
d2agg = aggregate(in.enh ~ queryHits, dovl, max)
p2agg = aggregate(in.enh ~ queryHits, povl, max)
a2agg = aggregate(in.enh ~ queryHits, aovl, max)
paste0('Enh in dELS: ', round(100 * sum(d2agg$in.enh) / NE, 2),'% of ', NE, ' elements')
paste0('Enh in pELS: ', round(100 * sum(p2agg$in.enh) / NE, 2),'% of ', NE, ' elements')
paste0('Enh in all: ', round(100 * sum(a2agg$in.enh) / NE, 2),'% of ', NE, ' elements')
# Enh in dELS: 32.49% of 2,069,090 elements
# Enh in pELS: 3.33% of 2,069,090 elements
# Enh in all: 37.81% of 2,069,090 elements
a2agg = aggregate(in.enh ~ queryHits + type, aovl, max)
a2sum = aggregate(in.enh ~ type, a2agg, sum)
a2sum$frac = paste0(round(100 * a2sum$in.enh / NE, 2),'%')
a2sum = a2sum[order(a2sum$in.enh, decreasing=T),]
#                      dELS 457005 22.09%
#           dELS,CTCF-bound 221144 10.69%
#                      pELS  42925  2.07%
#           pELS,CTCF-bound  26658  1.29%
#      CTCF-only,CTCF-bound  22581  1.09%
#             DNase-H3K4me3   6979  0.34%
#  DNase-H3K4me3,CTCF-bound   4830  0.23%
#                       PLS   4310  0.21%
#            PLS,CTCF-bound   3590  0.17%

agg2 = aggregate(cbind(in.enh, in.prom, in.dyad) ~ queryHits + type, aovl, sum)
agg2$sum = agg2$in.enh + agg2$in.prom + agg2$in.dyad
agg2 = agg2[agg2$sum > 0,]
agg2$set[(agg2$in.enh > agg2$in.prom) & (agg2$in.enh > agg2$in.dyad)] = 'Enhancer'
agg2$set[agg2$in.prom > agg2$in.enh & agg2$in.prom > agg2$in.dyad] = 'Promoter'
agg2$set[agg2$in.dyad > agg2$in.enh & agg2$in.dyad > agg2$in.prom] = 'Dyadic'
agg2$set[(agg2$in.enh == agg2$in.dyad) & agg2$in.enh > 0] = 'Enhancer'
agg2$set[(agg2$in.prom == agg2$in.dyad) & agg2$in.prom > 0] = 'Promoter'
agg2$set[(agg2$in.enh == agg2$in.prom) & agg2$in.enh > 0] = 'Enhancer'

totsum2 = aggregate(queryHits ~ type + set, agg2, length)
totsum2$set = factor(totsum2$set, levels=rev(c('Promoter','Enhancer','Dyadic')))

fdf = aggregate(queryHits ~ set, totsum2, sum)
totsum2 = rbind(totsum2, data.frame(type='None',set='Enhancer', queryHits=length(enhind) - fdf$queryHits[fdf$set == 'Enhancer'])) 
totsum2 = rbind(totsum2, data.frame(type='None',set='Promoter', queryHits=length(promind) - fdf$queryHits[fdf$set == 'Promoter']))
totsum2 = rbind(totsum2, data.frame(type='None',set='Dyadic', queryHits=length(dyadind) - fdf$queryHits[fdf$set == 'Dyadic']))


sc.cols = c('PLS' = '#ff0000ff', 'pELS' = '#ffa700ff',
            'dELS' = '#ffcd00ff', 'CTCF-only' = '#00b0f0ff',
            'DNase-H3K4me3' = '#ffaaaaff', 'None' ='grey75')

totsum2$main = sapply(totsum2$type, function(x){strsplit(x,",")[[1]][1]})
totsum2$main = factor(totsum2$main, levels=rev(names(sc.cols)))
tt2 = aggregate(queryHits ~ set + main, totsum2, sum)
twide2 = spread(tt2, main, queryHits)
tmat2 = as.matrix(twide2[,-1])
rownames(tmat2) = twide2[,1]
tmat2 = tmat2 / rowSums(tmat2)
round(tmat2, 4)


gp = ggplot(totsum2, aes(set, queryHits, fill=main)) + 
    geom_bar(stat='identity', position='fill', width=.98) + 
    scale_fill_manual(values=sc.cols, name='SCREEN:') + 
    coord_flip() + 
    labs(y='Percent of Epimap elements recovered', x='Epimap element set') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    theme_pubr()
ggsave(paste0(imgpref, 'reverse_element_recovery_pct.png'), gp, dpi=450, units='in', width=7, height=2.5)
ggsave(paste0(imgpref, 'reverse_element_recovery_pct.pdf'), gp, dpi=450, units='in', width=7, height=2.5)

gp2 = ggplot(totsum2, aes(set, queryHits, fill=main)) + 
    geom_bar(stat='identity', position='stack', width=.98) + 
    scale_fill_manual(values=sc.cols, name='SCREEN:') + 
    coord_flip() + 
    labs(y='Number of Epimap elements recovered', x='Epimap element set') + 
    scale_y_continuous(labels=scales::comma, expand=c(0,0)) + 
    theme_pubr()
ggsave(paste0(imgpref, 'reverse_element_recovery_num.png'), gp2, dpi=450, units='in', width=7, height=2.5)
ggsave(paste0(imgpref, 'reverse_element_recovery_num.pdf'), gp2, dpi=450, units='in', width=7, height=2.5)

# Write out the enhancer chunk names for epimap-only/screen-only/both
# -------------------------------------------------------------------
agg2$name = ir$id[agg2$queryHits]
agg2$main = sapply(agg2$type, function(x){strsplit(x,",")[[1]][1]})
denh = unique(agg2$name[agg2$set == 'Enhancer' & agg2$main == 'dELS'])
nenh = unique(agg2$name[agg2$set == 'Enhancer' & agg2$main != 'dELS'])
# Check fully in:
sum(denh %in% enhnam) / length(denh)
sum(nenh %in% enhnam) / length(nenh)
# Epimap-only
eenh = enhnam[!(enhnam %in% c(denh, nenh))]

# Write each out:
write(denh, 'epimap_screen_dELS.txt')
write(nenh, 'epimap_screen_non_dELS.txt')
write(eenh, 'epimap_non_screen_elements.txt')


# TODO: Make plots from some of these:

# Eval quality of each (only epimap - both - only screen) 
# - how many of the GWAS loci do they not recover? 
# - Enhancer modules
# - partition into ours/theirs -> Motifs/GWAS more/less enriched??

# Partition modules into the three <-> compare enrichments for three sets...

# -------------------------------
# Enhancers per epigenome/module:
# -------------------------------
# Enhancers with matches: 
enhdf$in.dELS = enhdf$name %in% dovl$name
enhdf$in.pELS = enhdf$name %in% povl$name
enhdf$in.all = enhdf$name %in% aovl$name

aodf = c()
for (flatset in c('epigenomes','modules')){
    enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
    load(enhsetfile)
    NF = length(enhsets)

    # Which have high epigenome occ: 
    occdf = c()
    for (i in 1:NF){
        print(i)
        x = enhsets[[i]]
        nx = length(x)
        subdf = enhdf[enhdf$cls %in% x,]
        sdf = data.frame(i=i, nx=nx, 
                         in.dELS=sum(subdf$in.dELS) / nx,
                         in.pELS=sum(subdf$in.pELS) / nx,
                         in.all=sum(subdf$in.all) / nx)
        occdf = rbind(occdf, sdf)
    }

    if (flatset == 'modules'){
        load('motifs/motif_modules_network_data.Rda')
        occdf$COLOR = clscols[occdf$i]
        occdf$id = paste0('id',occdf$i)
        occdf = merge(occdf, odf)
    } else if (flatset == 'epigenomes'){
        matnames = scan('Enhancer_matrix_names.txt', "c")
        occdf$id = matnames[occdf$i]
        occdf = merge(occdf, meta[,c('id','GROUP','COLOR')])
    }
    occdf$set = flatset
    cols = c('i','id','GROUP','COLOR','set','nx','in.dELS','in.pELS','in.all')
    aodf = rbind(aodf, occdf[,cols])

}
aodf$GROUP = factor(aodf$GROUP, levels=odf$GROUP)

# Look at distributions:
gp = ggplot(aodf, aes(set, in.dELS, color=GROUP)) + 
    geom_violin(color='black', scale='width') + 
    geom_boxplot(width=.25, color='black') + 
    geom_quasirandom() + 
    scale_color_manual(values=colvals$group) + 
    coord_flip() + theme_pubr() + 
    labs(x='Set',y='Percent of elements in SCREEN dELS') + 
    theme(legend.position='none')

ggsave(paste0(imgpref, 'modepi_recovery.png'), gp, dpi=450, units='in', width=6, height=4)
ggsave(paste0(imgpref, 'modepi_recovery.pdf'), gp, dpi=450, units='in', width=6, height=4)

# Plot the margin + recovery:
enhmargfile = 'Enhancer_H3K27ac_margins_enhonly_062619.tsv.gz'
matmarg = read.delim(gzfile(enhmargfile), sep="\t", header=T)
names(matmarg) = c('cls','margin')
enhdf = merge(enhdf, matmarg)
# Look at distributions:
mdf = aggregate(in.dELS ~ margin, enhdf, mean)

gp = ggplot(mdf, aes(margin, in.dELS)) + 
    geom_line() + theme_pubr() + 
    geom_smooth() + 
    labs(x='Number of epigenomes containing element',y='Percent of elements in SCREEN dELS') + 
    scale_y_continuous(labels=scales::percent) + 
    scale_x_log10() + theme(legend.position='none')

ggsave(paste0(imgpref, 'recov_rarity.png'), gp, dpi=450, units='in', width=4, height=4)
ggsave(paste0(imgpref, 'recov_rarity.pdf'), gp, dpi=450, units='in', width=4, height=4)




