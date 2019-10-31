#!/usr/bin/R
# --------------------------------------------------------
# Evaluate the deltas between tracks for technical issues.
# - Load distances of every delta to all within-sample tracks
# - Load distances of every delta to all within-mark tracks
# Plot the distribution of M1-M2 for normal cases, find outliers.
# (similarly for S1-S2, but modeling the type of mark as well - gives distribution size)
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
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")

# Load the distance matrices (for within-mark):
source(paste0(bindir, 'load_distance_matrices.R'))
# Full is the distance matrix:
diag(full) = 0
cn = colnames(full)

recent='20190310-1846'
wsfile = paste0('all_diffdist_',recent,'.tsv.gz')
wmfile = paste0('all_sampdiffdist_',recent,'.tsv.gz')
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need within-sample and within-mark filenames.")
} else {        
    wsfile = args[1]
    wmfile = args[2]
}

# Directories:
imgpref = paste0(img, "diff_qc/")
system(paste('mkdir -p', imgpref))

# ------------------------------------------------------
# Within-sample analysis: 
# NOTE: Looking for 2ndary antibodies
# or a mixture of two Ab due to ab swap in one expt.
# ------------------------------------------------------
# Read in ws table
wsdf = read.delim(wsfile, stringsAsFactors=F)
wsdf$cor = 1 - wsdf$dist

# Histograms of each against each other:
marks = sort(unique(as.character(wsdf$mark)))
NMARK = length(marks)
r = max(-min(wsdf$cor), max(wsdf$cor))

pdf(paste0(imgpref, 'withinsample_mark_comparison_boxplots.pdf'),width=11,height=5)
layout(matrix(c(1:((NMARK + 1)* NMARK)),NMARK+1,NMARK), TRUE)
par(yaxs='i')
par(xaxs='i')
for (smark in marks){
    par(mar=rep(0.0,4))
    plot(1,1,type='n',axes=F)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]),
         smark, col=rgb(0,0,0,.45))
    abline(v=par()$usr[1:2])
    abline(h=par()$usr[3:4])
    for (amark in marks){
        par(mar=rep(0.00,4))
        subdf = filter(wsdf, mark == smark & against == amark)
        # hist(subdf$cor, xlab='',ylab='', main='', axes=F, xlim=c(-.75,.75))
        boxplot(subdf$cor, xlab='',ylab='', main='', axes=F, 
                ylim=c(-r,r), horizontal=T, pch=19, lwd=.5, cex=.25)
        abline(v=0, lty='dotted', col=rgb(0,0,1,.5), lwd=1)
        abline(v=par()$usr[1:2], lty='dotted',col='grey',lwd=.5)
        abline(h=par()$usr[3:4], lty='dotted',col='grey',lwd=.5)
        if(amark == smark){
            abline(v=par()$usr[1:2])
            abline(h=par()$usr[3:4])
        }
    }
}
dev.off()


pdf(paste0(imgpref, 'withinsample_mark_comparison_histograms.pdf'),width=11,height=5)
layout(matrix(c(1:((NMARK + 1)* NMARK)),NMARK+1,NMARK), TRUE)
par(yaxs='i')
par(xaxs='i')
for (smark in marks){
    par(mar=rep(0.0,4))
    plot(1,1,type='n',axes=F)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]),
         smark, col=rgb(0,0,0,.45))
    abline(v=par()$usr[1:2])
    abline(h=par()$usr[3:4])
    for (amark in marks){
        par(mar=rep(0.00,4))
        subdf = filter(wsdf, mark == smark & against == amark)
        hist(subdf$cor, xlab='',ylab='', main='', axes=F, xlim=c(-r,r))
        # boxplot(subdf$cor, xlab='',ylab='', main='', axes=F, 
        #         ylim=c(-r,r), horizontal=T, pch=19, lwd=.5, cex=.25)
        abline(v=0, lty='dotted', col=rgb(0,0,1,.5), lwd=1)
        abline(v=par()$usr[1:2], lty='dotted',col='grey',lwd=.5)
        abline(h=par()$usr[3:4], lty='dotted',col='grey',lwd=.5)
        if(amark == smark){
            abline(v=par()$usr[1:2])
            abline(h=par()$usr[3:4])
        }
    }
}
dev.off()

# Assuming approximate normality, calculate the z-score for each:
mcor = aggregate(cor ~ mark + against, wsdf, mean)
scor = aggregate(cor ~ mark + against, wsdf, sd)
names(mcor)[3] = 'mean'
names(scor)[3] = 'sd'
wsdf = merge(merge(wsdf, mcor), scor)
wsdf$zscore = (wsdf$cor - wsdf$mean) / wsdf$sd

# Assign unique ID to each track (delta sample):
uqdata = unique(wsdf[,c('sample', 'mark')])
NTRACK = dim(uqdata)[1]
print(paste("[STATUS] There are", NTRACK, "unique tracks"))
uqdf = data.frame(uqdata, uid=1:NTRACK)
wsdf = merge(wsdf, uqdf)
# Merge with broad/punctate definitions
wsdf = merge(wsdf, markdef)
amarkdef = markdef
names(amarkdef) = c('against','against.type')
wsdf = merge(wsdf, amarkdef)

lsamp = c(aggregate(sample ~ mark, uqdf, length)$sample)

# Plot zscores by mark, highlight/count up the ones that are consistently discrepant from -3 to 3
width = 11
height = 3
# Buffer is for y-axis.
buffer = 120
widths = c(buffer, lsamp) * 11 / (buffer + sum(lsamp))
yr = max(-min(wsdf$zscore), max(wsdf$zscore)) * 1.1
colmarktype = c('punctate' = 'royalblue', 'broad' = 'indianred')

# pdf(paste0(imgpref, 'withinsample_zscore_scatters.pdf'),width=width, height=height)
png(paste0(imgpref, 'withinsample_zscore_scatters.png'),res=300, units='in',width=width, height=height)
layout(matrix(c(1:(NMARK+1)), 1, NMARK + 1), widths=widths, TRUE)
par(xaxs="i")
par(yaxs="i")
par(mar=c(0,4,0,0))
plot(1,1, ylim=c(-yr,yr), type='n', xaxt='n', bty='n',
     ylab='Z-score for cross-mark, within-sample comparison')
for (smark in marks){
    par(mar=c(0,.1,0,.4))
    subdf = filter(wsdf, mark == smark & against != smark)
    hdf = filter(subdf, abs(zscore) > 2)
    ldf = filter(subdf, abs(zscore) < 2)
    xlim = range(subdf$uid)
    plot(ldf$uid, ldf$zscore, ylim=c(-yr,yr), axes=F, xlim=xlim,
         pch=19, cex=.35, col='grey')
    points(hdf$uid, hdf$zscore, ylim=c(-yr,yr), axes=F, xlim=xlim,
           pch=19, cex=.35, col=colmarktype[hdf$against.type])
    abline(h=c(-2,0,2), lty='dotted',col=rgb(0,0,0,.5),lwd=1)
    text(mean(par()$usr[1:2]), 
         par()$usr[3] + 0.05 * diff(par()$usr[3:4]),
         smark, col=colmarktype[subdf$type[1]], srt=90, adj=0)
}
dev.off()

# ------------------------------------------------------
# Within-mark analysis: 
# NOTE: Looking for potential sample swaps 
# or mixtures of two samples due to two diff expt merge
# or mixtures of two cell types due to sample collection
# ------------------------------------------------------
# Read in within-mark table
wmdf = read.delim(gzfile(wmfile), stringsAsFactors=F, header=F)
names(wmdf) = c('sample','mark','against','dist')
samples = sort(unique(as.character(wmdf$sample)))
NSAMP = length(samples)
wmdf$cor = 1 - wmdf$dist
wmdf = merge(wmdf, uqdf)

# First look at max cor (compared to other samp comparisons within same type):
# NOTE: Can id some poor samples but they are likely Ab issues.
mxdf = aggregate(cor ~ uid + sample + mark, wmdf, max)
ggplot(mxdf, aes(mark, cor)) + geom_violin()
# ggplot(mxdf, aes(mark, cor)) + geom_boxplot()

# Evaluate both spearman and pearson correlation:
cordf = c()
pb = txtProgressBar(min = 0, max = nrow(uqdf), style = 3)
for (i in 1:nrow(uqdf)){
    setTxtProgressBar(pb, i)
    subdf = filter(wmdf, uid == i)
    samp = subdf$sample[1]
    maxsamp = subdf$against[which.max(subdf$cor)]
    agg = aggregate(cor ~ against, subdf, mean)
    x = agg$cor
    nam = agg$against
    y = 1 - full[samp, nam]
    scor = cor(x, y, method='spearman')
    pcor = cor(x, y, method='pearson')
    cordf = rbind(cordf, data.frame(uid=i, pearson=pcor, spearman=scor, maxsamp=maxsamp))
}
close(pb)
cordf = merge(uqdf, cordf)
cldf = gather(cordf, cor.type, Correlation, -uid, -sample, -mark, -maxsamp)

# Use zscores to label
mcor = aggregate(Correlation ~ mark + cor.type, cldf, mean)
scor = aggregate(Correlation ~ mark + cor.type, cldf, sd)
names(mcor)[3] = 'mean'
names(scor)[3] = 'sd'
cldf = merge(merge(cldf, mcor), scor)
cldf$zscore = (cldf$Correlation - cldf$mean) / cldf$sd

# Look at range of agreements:
labdf = filter(cldf, zscore < -2)
ggplot(cldf, aes(cor.type, Correlation, fill=cor.type)) + 
    facet_wrap(~mark, nrow=1) + 
    geom_boxplot(alpha=0.5) + 
    geom_text_repel(data=labdf, aes(cor.type, Correlation, label=sample), alpha=0.75, size=2) + 
    theme_bw() + 
    xlab('Correlation Type') + 
    ylab('Correlation of Delta Correlations\n with Full Sample Correlation Matrix') + 
    theme(legend.position='none') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))
ggsave(paste0(imgpref, 'withinmark_sample_agreement.png'),dpi=350, units='in', width=11, height=4)

# Add the info lines to say what is mixed (rather than BSSID)
mo = meta[,c('id','infoline')]
names(mo)[1] = c('sample')
cldf = merge(cldf, mo)
# Merge with maxsamp too:
names(mo) = c('maxsamp', 'maxinfo')
cldf = merge(cldf, mo)

# Look at range of agreements:
subdf = filter(cldf, cor.type=='pearson')
labdf = filter(subdf, zscore < -2)
ggplot(subdf, aes(mark, Correlation)) + 
    geom_boxplot(col='darkgrey', outlier.size=1) + 
    geom_text_repel(data=labdf, aes(mark, Correlation, label=sample), 
                    size=1.5, segment.color='grey') + 
    lims(y=c(-1,1)) +
    theme_bw() + 
    xlab('Histone Mark/Assay') + 
    ylab('Pearson Correlation of Delta Correlations\n with Full Sample Correlation Matrix') + 
    theme(legend.position='none') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7))
ggsave(paste0(imgpref, 'withinmark_sample_agreement_pearson.png'),dpi=350, units='in', width=5, height=4)

# Look at range of agreements:
subdf = filter(cldf, cor.type=='pearson')
labdf = filter(subdf, zscore < -2)
ggplot(subdf, aes(mark, Correlation)) + 
    geom_boxplot(col='darkgrey', outlier.size=1) + 
    geom_text_repel(data=labdf, aes(mark, Correlation, label=paste0(infoline, "\n(to ",maxinfo, ")")), 
                    size=1, segment.color='grey') + 
    lims(y=c(-1,1)) +
    theme_bw() + 
    xlab('Histone Mark/Assay') + 
    ylab('Pearson Correlation of Delta Correlations\n with Full Sample Correlation Matrix') + 
    theme(legend.position='none') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7))
ggsave(paste0(imgpref, 'withinmark_sample_agreement_pearson_infoline.png'),dpi=350, units='in', width=5, height=4)


# From here, output the potentially poor quality ones:
outdf = labdf[order(labdf$zscore), c('sample', 'mark', 'Correlation', 'zscore', 'maxsamp', 'infoline','maxinfo')]

print(paste("Returning",dim(outdf)[1],"potential 2ndary samples."))
write.table(outdf, paste0('flagged_second_samp_', recent, '.tsv'), quote=F, row.names=F, sep="\t")

# NOTE: Also could look at distribution, zscore may not be best measure for flagging data.
# ggplot(cldf, aes(Correlation, fill=cor.type)) + 
#     geom_histogram(alpha=0.5) +
#     geom_vline(xintercept=0, lty='dashed', color='grey') + 
#     facet_grid(cor.type~mark) + 
#     theme_minimal()
