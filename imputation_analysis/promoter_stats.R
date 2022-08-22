#!/usr/bin/R
# Calculate promoter statistics vs. current gencode
library(GenomicRanges)

# Read in called promoters (all files pulled from web):
promind = as.numeric(scan('PROM_masterlist_indices_0indexed.tsv', 'c')) + 1
dhsdf = read.delim('masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt',header=F)
promdf = dhsdf[promind,]
rm(dhsdf); gc()
names(promdf) = c('chr','start','end','chunk')
promdf$mid = (promdf$start + promdf$end) / 2

# Read in gencode annotation:
gtf = read.delim('gene.gencode.v38lift37.annotation.tsv.gz',header=F)
names(gtf) = c('chr','src','class','start','end','v1','strand','v2','ensg','type','symbol')

gtf$tss = gtf$start
nind = gtf$strand == '-'
gtf$tss[nind] = gtf$end[nind]
gtf = gtf[!is.na(gtf$tss),]

# Make GenomicRanges objects:
subgtf = gtf[gtf$type == 'lncRNA',]
usegtf = subgtf
usegtf = gtf
pr = with(promdf,GRanges(seqnames=chr, IRanges(start, end)))
gr = with(usegtf,GRanges(seqnames=chr, IRanges(tss, tss+1)))

# Find nearest:
novl1 = nearest(pr, gr)
dist.ptog1 = promdf$mid - usegtf$tss[novl1]
type1 = usegtf$type[novl1]

dcut = 10000
dist.ptog1[dist.ptog1 > dcut] = dcut
dist.ptog1[dist.ptog1 < -dcut] = -dcut


novl2 = nearest(gr, pr)
dist.ptog2 = promdf$mid[novl2] - usegtf$tss
type2 = usegtf$type

dcut = 10000
dist.ptog2[dist.ptog2 > dcut] = dcut
dist.ptog2[dist.ptog2 < -dcut] = -dcut

par(xaxs='i')
par(yaxs='i')
layout(matrix(1:2, ncol=1))
hist(dist.ptog1, 50, main='Histogram of distance to nearest GENCODE TSS for each annotated promoter (GRCh37)', border='white', xlab='Distance (bp)')
hist(dist.ptog2, 50, main='Histogram of distance to nearest Epimap promoter for each GENCODE TSS (GRCh37)', border='white', xlab='Distance (bp)')

df1 = data.frame(dist=dist.ptog1, type=type1)
df2 = data.frame(dist=dist.ptog2, type=type2)

library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggrepel)

snapcols = scan('snap_colors.tsv','c')
pcols = brewer.pal(12, 'Paired')

types = unique(gtf$type)
j = 2
tcols = snapcols[(j+1):(j + length(types))]
names(tcols) = types

gp1 = ggplot(df1, aes(dist, fill=type)) + 
    geom_histogram(color=NA, position='stack') + 
    scale_fill_manual(values=tcols) + 
    scale_y_continuous(expand=c(0,0), labels=scales::comma) + 
    scale_x_continuous(expand=c(0,0), labels=scales::comma) + 
    labs(title='Histogram of distance to nearest GENCODE TSS for each annotated promoter (GRCh37)', x='Distance (bp)', y='Number of elements') + 
    theme_pubr()
ggsave('promstats_dist_to_gencode.pdf',gp1,dpi=300, units='in',width=8, height=7)

gp2 = ggplot(df2, aes(dist, fill=type)) + 
    geom_histogram(color=NA, position='stack') + 
    scale_fill_manual(values=tcols) + 
    scale_y_continuous(expand=c(0,0), labels=scales::comma) + 
    scale_x_continuous(expand=c(0,0), labels=scales::comma) + 
    labs(title='Histogram of distance to nearest EpiMap promoter for each GENCODE TSS (GRCh37)', x='Distance (bp)', y='Number of elements') + 
    theme_pubr()
ggsave('promstats_dist_to_epimap.pdf',gp2,dpi=300, units='in',width=8, height=7)


# Total (# promoters per gene):







# Find overlaps:



# Stats on promoters:


