#!/usr/bin/R
# Plot the data origins/timeline
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
library(ggrepel)

# Locations:
imgdir = paste0(img, "metadata/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'meta_')

# Read extended files metadata (for histones + DNase):
accmap = read.delim('accession_bssid_mark.tsv', header=F, stringsAsFactors=F)
names(accmap) = c('Accession','id','Epitope')
filemeta = read.delim('Annotation/file_extended_metadata_20180128.tsv', header=T, stringsAsFactors=F)
accmap = merge(accmap, unique(filemeta[,c('Accession','Project', 'replicates.library.biosample.donor.accession', 'Lab', 'Month.released')]), all.x=T)

print(aggregate(Accession ~ Lab, accmap, length))

accmap$month = sub(", .*", "", accmap$Month.released)
accmap$month = factor(accmap$month, levels=month.name)
accmap$year = sub(".*,", "", accmap$Month.released)
accmap = accmap[order(accmap$month, decreasing=F),]
accmap = accmap[order(accmap$year, decreasing=F),]
accmap$Month.released = factor(accmap$Month.released, 
                               levels=unique(accmap$Month.released))

lab.cols = brewer.pal(10, 'Paired')

# All labs
gplot = ggplot(accmap, aes(Lab, fill=Lab)) + 
    geom_bar() + theme_bw() + 
    ylab('Number of assays') + 
    theme(legend.position = 'none')+ 
    scale_fill_manual(values=lab.cols) + 
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'assays_by_lab.png'),gplot, units='in', dpi=300, width=5, height=7)
ggsave(paste0(imgpref, 'assays_by_lab.pdf'),gplot, units='in', dpi=300, width=5, height=7)


# All labs, by month 
gplot = ggplot(accmap, aes(Month.released, fill=Lab)) + 
    geom_bar() +
    theme_bw() + 
    ylab('Number of assays') + 
    theme(legend.position = 'top')+ 
    scale_fill_manual(values=lab.cols) + 
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'assays_by_lab_month.png'),gplot, units='in', dpi=300, width=9, height=6)
ggsave(paste0(imgpref, 'assays_by_lab_month.pdf'),gplot, units='in', dpi=300, width=9, height=6)

# Facetted by lab:
gplot = ggplot(accmap, aes(Month.released, fill=Lab)) + 
    facet_wrap(~Lab, scales='free_y') + 
    geom_bar() + theme_bw() + 
    ylab('Number of assays') + 
    theme(legend.position = 'top')+ 
    scale_fill_manual(values=lab.cols) + 
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'assays_by_lab_month_scales_bylab.png'),gplot, units='in', dpi=300, width=15, height=8)
ggsave(paste0(imgpref, 'assays_by_lab_month_scales_bylab.pdf'),gplot, units='in', dpi=300, width=15, height=8)

# Facetted by type of assay:
gplot = ggplot(accmap, aes(Month.released, fill=Lab)) + 
    facet_wrap(~Epitope, scales='free_y') + 
    geom_bar() + theme_bw() + 
    ylab('Number of assays') + 
    theme(legend.position = 'top')+ 
    scale_fill_manual(values=lab.cols) + 
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'assays_by_lab_month_byepitope.png'),gplot, units='in', dpi=300, width=25, height=10)
ggsave(paste0(imgpref, 'assays_by_lab_month_byepitope.pdf'),gplot, units='in', dpi=300, width=25, height=10)


# Finally, plot a cumulative plot of assays:
acclen = aggregate(Accession ~ Lab + Month.released, accmap, length)
all.pos = expand.grid(unique(acclen$Lab), unique(acclen$Month.released))
names(all.pos) = c('Lab','Month.released')
acclen = merge(acclen, all.pos, all.y=TRUE)
acclen$Accession[is.na(acclen$Accession)] = 0
acclen = acclen[order(acclen$Month.released, decreasing=F),]
acclen$cs = 0

for (lab in unique(acclen$Lab)){
    idx = which(acclen$Lab == lab)
    acclen$cs[idx] = cumsum(acclen$Accession[idx])
}

# All labs, by month 
gplot = ggplot(acclen, aes(Month.released, cs, fill=Lab)) + 
    geom_bar(stat='identity') + theme_bw() + 
    ylab('Number of assays') + 
    theme(legend.position = 'top')+ 
    scale_fill_manual(values=lab.cols) + 
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'assays_by_lab_month_cumulative.png'),gplot, units='in', dpi=300, width=9, height=6)
ggsave(paste0(imgpref, 'assays_by_lab_month_cumulative.pdf'),gplot, units='in', dpi=300, width=9, height=6)


accmap2 = merge(accmap, meta[,c('id','Project')])


# All labs, by month 
gplot = ggplot(accmap, aes(Month.released, fill=Lab)) + 
    facet_wrap(~Project, nrow=3, scales='free_y') + 
    geom_bar() + theme_bw() + 
    ylab('Number of assays') + 
    theme(legend.position = 'top')+ 
    scale_fill_manual(values=lab.cols) + 
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'assays_by_lab_month_byproj.png'),gplot, units='in', dpi=300, width=9, height=9)
ggsave(paste0(imgpref, 'assays_by_lab_month_byproj.pdf'),gplot, units='in', dpi=300, width=9, height=9)


all.pos2 = expand.grid(unique(acclen$Lab), month.name, unique(accmap$year))
names(all.pos2) = c('Lab','month', 'year')
all.pos2$Month.released = paste0(all.pos2$month, ',', all.pos2$year)
mlvl = unique(all.pos2$Month.released)

acclen = aggregate(Accession ~ Lab + Month.released, accmap, length)
acclen$Month.released = as.character(acclen$Month.released)
acclen = merge(acclen, all.pos2, all.y=TRUE)
acclen$Accession[is.na(acclen$Accession)] = 0
acclen$Month.released = factor(acclen$Month.released, levels=mlvl)

# All labs, by month 
gplot = ggplot(acclen, aes(Month.released, Accession, fill=Lab)) + 
    geom_bar(stat='identity') + theme_bw() + 
    ylab('Number of assays') + 
    theme(legend.position = 'top')+ 
    scale_fill_manual(values=lab.cols) + 
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(paste0(imgpref, 'assays_by_lab_month_allmonth.png'),gplot, units='in', dpi=300, width=12, height=6)
ggsave(paste0(imgpref, 'assays_by_lab_month_allmonth.pdf'),gplot, units='in', dpi=300, width=12, height=6)












