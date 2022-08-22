#!/usr/bin/R
# -------------------------------------------
# Use distances obs v. all to id sample swaps
# -------------------------------------------
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
library(scales)
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
markdef = rbind(data.frame(mark=pm, type='punctate'), 
                data.frame(mark=bm, type='broad'), 
                data.frame(mark=fc, type='factor'))
metricsdf = merge(metricsdf, markdef)
aggregate(mark ~ type, metricsdf, length)

# write.table(metricsdf[,c('id','mark')], file='impobs_comb_tocompare.tsv', 
#             quote=F, sep="\t", row.names=F, col.names=F)

# View the metrics as a heatmaps:
logit = function(x, scale=1){ x = x / scale; log(x / (1 - x)) }
sigm = function(x, scale=1){ scale / (1 + exp(-x)) }

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
qcwide = merge(qcdf, metricsdf)


# ------------------------------
# Load the distance comparisons:
# ------------------------------
compdf = read.delim(gzfile('ChromImpute/all_obscomp_t10.tsv.gz'), header=F, sep="\t")
names(compdf) = c('id','mark','against','avg.corr')
compdf$avg.corr = 1 - compdf$avg.corr
cwide = spread(compdf, against, avg.corr)

filter(compdf, id == 'BSS00132')
filter(compdf, id == 'BSS00132' & against == 'H3K4me1')
filter(compdf, id == 'BSS00132' & mark == 'H3K4me1')

# -------------------------------------
# Plot each mark, flag specific tracks:
# -------------------------------------
# kmark = 'H3K27ac'
all.df = c()
all.labdf = c()
for (kmark in unique(compdf$mark)){
    print(kmark)
    subdf = filter(cwide, mark == kmark)
    subdf = gather(subdf, against, avg.corr, -id, -mark, -kmark)
    names(subdf)[3] = 'kmark'
    subdf = merge(subdf, metricsdf[,c('mark','id','status')])
    # Find the standard deviation in either direction,
    # Needs to be out in both directions to be labeled:
    subdf = merge(subdf, setNames(aggregate(avg.corr ~ against, subdf, mean), 
                                  c('against','m.corr')))
    subdf = merge(subdf, setNames(aggregate(avg.corr ~ against, subdf, sd), 
                                  c('against','sd.corr')))
    subdf$m.own = mean(unique(subdf[,c('mark','id','kmark')])$kmark)
    subdf$sd.own = sd(unique(subdf[,c('mark','id','kmark')])$kmark)
    subdf$z.corr = (subdf$avg.corr - subdf$m.corr) / subdf$sd.corr
    subdf$z.own = (subdf$kmark - subdf$m.own) / subdf$sd.own
    # NOTE: better would be to regress out and get residuals
    model = lm(avg.corr ~ against + kmark, subdf)
    subdf$pred.lm = predict(model, subdf)
    # Not really accurate match, but works decently:
    subdf$z.pred = (subdf$avg.corr - subdf$pred.lm) / subdf$sd.corr
    # Save the flagged tracks.
    cutoff = head(sort(subdf$z.pred, decreasing=TRUE), 5)[5]
    labdf = filter(subdf, z.pred > min(10, cutoff))
    # Save the data.frames:
    all.df = rbind(all.df, subdf)
    all.labdf = rbind(all.labdf, labdf)

    gp = ggplot(subdf, aes(kmark, avg.corr)) +
        geom_point(aes(color=status), cex=.5) +
        geom_text_repel(data=labdf, aes(kmark, avg.corr, label=id), cex=2) +
        geom_smooth(method='lm') + 
        scale_color_manual(values=ovals) + 
        facet_wrap(~against) +
        ylim(0,1) + xlim(0,1) + 
        xlab(paste0('Avg. Corr. with top 10 closest ', kmark, ' tracks')) +
        ylab(paste0('Avg. Corr. of other mark')) +
        theme_pubr() + theme(legend.position='none') + 
        theme(axis.text.x= element_text(size=8)) +
        theme(axis.text.y= element_text(size=8))
    ggsave(paste0(imgpref, 'flag_swaps_', kmark, '.pdf'), gp,  width=9, height=6, units='in')
}



# TODO: No space between boxes
# TODO fix the y axis
# Plot all comparisons:
gp = ggplot(all.df, aes(kmark, avg.corr)) +
    facet_grid(mark ~ against) +
    geom_point(aes(color=status), cex=.25) +
    geom_point(data=all.labdf, aes(kmark, avg.corr, color=status), cex=.5, pch=21, fill='white') +
    geom_text_repel(data=all.labdf, aes(kmark, avg.corr, label=id), cex=1, segment.size=.25) +
    geom_smooth(method='lm', lwd=.5) + 
    scale_color_manual(values=ovals) + 
    scale_fill_manual(values=ovals) + 
    scale_x_continuous(breaks=c(0,.5,1), labels=c(0,.5,1), lim=c(0,1)) + 
    scale_y_continuous(breaks=c(0,.5,1), labels=c(0,.5,1), lim=c(0,1)) + 
    xlab(paste0('Avg. Corr. with top 10 closest tracks of putative mark')) +
    ylab(paste0('Avg. Corr. with top 10 closest tracks of comparison mark')) +
    theme_bw() + theme(legend.position='none') + 
    theme(axis.text.x= element_text(size=6)) +
    theme(axis.text.y= element_text(size=6)) + 
    theme(panel.spacing = unit(0, "lines"))
ggsave(paste0(imgpref, 'flag_swaps_all.pdf'), gp,  width=10, height=10, units='in')


# ----------------------
# Evaluate sample swaps:
# ----------------------
# Samples that aren't ab-swapped, but show much lower corr with own than max:
samedf = compdf[compdf$mark == compdf$against,]
samedf = merge(samedf, metricsdf[,c('mark','id','status','Correlation')])
samedf = merge(samedf, setNames(aggregate(avg.corr ~ against, samedf, mean), 
                                c('against','m.corr')))
samedf = merge(samedf, setNames(aggregate(avg.corr ~ against, samedf, sd), 
                                c('against','sd.corr')))
# Remove ab swaps:
amat = matrix(c('BSS00333', 'H3K36me3', 'BSS00333', 'H3K27ac',
                'BSS01857', 'H3K9ac',   'BSS01857', 'H3K79me2',
                'BSS00141', 'H3K9me3',  'BSS00132', 'H3K4me1',
                'BSS01420', 'H3K4me1'), ncol=2, byrow=T)
adf = data.frame(amat)
names(adf) = c('id','mark')
adf$abswap = 1
samedf = merge(samedf, adf, all.x=TRUE)
samedf = samedf[is.na(samedf$abswap),]
model = lm(Correlation ~ avg.corr + against, samedf)
samedf$pred.lm = predict(model, samedf)
samedf = merge(samedf, setNames(aggregate((pred.lm - Correlation) ~ against, samedf, sd), 
                                c('against','sd.pred')))
# Not really accurate match, but works decently:
samedf$z.pred = (samedf$pred.lm - samedf$Correlation) / samedf$sd.pred
# Save the flagged tracks.
cutoff = head(sort(samedf$z.pred, decreasing=TRUE), 5)[5]
labdf = filter(samedf, z.pred > min(3, cutoff))

gp = ggplot(samedf, aes(avg.corr, Correlation)) +
    geom_point(aes(color=status), cex=.5) +
    geom_text_repel(data=labdf, aes(avg.corr, Correlation, label=id), 
                    cex=2, nudge_y=-.025, segment.size=.25) +
    scale_color_manual(values=ovals) + 
    geom_abline(lwd=.5, alpha=.5) + 
    facet_wrap(~against, nrow=3) +
    ylim(0,1) + xlim(0,1) + 
    xlab(paste0('Avg. Corr. with top 10 closest same-mark tracks')) +
    ylab(paste0('Correlation with imputed track')) +
    theme_bw() + theme(legend.position='none') + 
    theme(axis.text.x= element_text(size=8)) +
    theme(axis.text.y= element_text(size=8))
ggsave(paste0(imgpref, 'flag_sampswaps_allmarks.pdf'), gp,  width=12, height=8, units='in')


# ---------------------------------------------------
# Figure out if potential for swap was in the top 25:
# Use to try to match sample swaps.
# ---------------------------------------------------
rdf = read.delim('ChromImpute/all_obscomp_ranked_compare.txt', sep="\t", header=F)
names(rdf) = c('id','mark','vid','corr','rank')
subrdf = rdf[rdf$rank == 1,]
# Want to find places where both in/out have high z.pred
scdf = samedf[,c('id','mark','z.pred')]
rdf = merge(rdf, scdf)
names(scdf) = c('vid', 'mark', 'vz.pred')
rdf = merge(rdf, scdf)
rdf = rdf[as.character(rdf$vid) != as.character(rdf$id),]
pdf = merge(rdf, labdf)
pdf$sumz = pdf$z.pred + pdf$vz.pred
# Prune to best possible matches:
pdf = merge(pdf, aggregate(sumz ~ mark + id, pdf, max))

# Add best match:
labdf2 = merge(labdf, subrdf[,c('id','mark','vid')])
labdf2 = labdf2[order(labdf2$mark), c('id','mark','vid')]
names(labdf2)[3] = 'potential.sampswap'

# -----------------------------------
# Write flagged samples/datasets out:
# -----------------------------------
write.table(adf[,c('id', 'mark')], 'Annotation/flagged_potential_abswap.tsv', quote=F, col.names=T, row.names=F, sep="\t")
write.table(labdf2[,c('id','mark', 'potential.sampswap')], 'Annotation/flagged_potential_sampleswap.tsv', quote=F, col.names=T, row.names=F, sep="\t")
lowdf = metricsdf[metricsdf$status == 'LOW', ]
write.table(lowdf[,c('id','mark')], 'Annotation/flagged_low_agreement.tsv', quote=F, col.names=T, row.names=F, sep="\t")

# ----------------------------------------------
# Count number of tracks not flagged and in 833:
# ----------------------------------------------
# Read track lists:
itdf = read.table('data_to_release/all_imputed_tracks.tsv', sep="\t", header=F, stringsAsFactors=F)
otdf = read.table('data_to_release/all_observed_tracks.tsv', sep="\t", header=F, stringsAsFactors=F)
ocols = c('id','mark','file')
names(itdf) = ocols
names(otdf) = ocols

# Subset to 833:
itdf = itdf[itdf$id %in% cellorder,]
otdf = otdf[otdf$id %in% cellorder,]

# Remove the flagged samples:
fldf = lowdf[,c('id','mark', 'status')]
otdf = merge(otdf, fldf, all.x=TRUE)
otdf = otdf[is.na(otdf$status),]

print(dim(itdf))
print(dim(otdf))
print(dim(itdf) + dim(otdf))

write.table(itdf[,ocols], 'data_to_release/postqc_imputed_tracks.tsv', sep="\t", row.names=F, col.names=F, quote=F)
write.table(otdf[,ocols], 'data_to_release/postqc_observed_tracks.tsv', sep="\t", row.names=F, col.names=F, quote=F)



