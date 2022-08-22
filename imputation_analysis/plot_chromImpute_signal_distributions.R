#!/usr/bin/R
# -------------------------------------------
# Plot the signal distributions of
# imputed and observed tracks under each mark
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
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")
library(dplyr)
library(ggplot2)
library(scales)

# Functions:
perc.rank <- function(x) trunc(rank(x))/length(x)

# Directories:
imgpref = paste0(img, "signal_distr/")
signaldir = "ChromImpute/signal_distribution/"
system(paste('mkdir -p', imgpref))

# Read in track dataset:
tracktab = read.delim('ChromImpute/all_impobs_tracks_uniq_table.tsv', header=F, stringsAsFactors=F)
names(tracktab) = c('uqsample','Epitope','prefix')
tracktab$sample = sapply(tracktab$uqsample, function(x){strsplit(x,"_")[[1]][1]})
tracktab$type = sapply(tracktab$uqsample, function(x){
                           a = strsplit(x,"_")[[1]][2]
                           if (a == 'imp'){ a = "Imputed" } else {a = "Observed" } })
marks = sort(unique(tracktab$Epitope))
samples = sort(unique(tracktab$sample))

# For each mark, plot the distribution (line per sample) 
# Split by imputed and observed.
mark = 'H3K27ac'
cutdf = c()
for (mark in t12marks[8:13]){
# for (mark in mainmarks){
    print(mark)
    subtrack = filter(tracktab, Epitope == mark)
    markfile = paste0('ChromImpute/signal_distr_',mark,'.tsv.gz')
    # Read in the distributions:
    if (file.exists(markfile)){
        markdf = read.delim(gzfile(markfile), stringsAsFactors=F, header=T)
    } else {
        markdf = c()
        pb = txtProgressBar(min = 0, max = nrow(subtrack), style = 3)
        for (i in 1:nrow(subtrack)){
            setTxtProgressBar(pb, i)
            # df = read.delim(paste0(signaldir, subtrack$prefix[i], '_signalcounts.tsv'), header=F, stringsAsFactors=F)
            filename = paste0(signaldir, subtrack$prefix[i], '_signal_bystate.tsv')
            if (file.exists(filename)){
                df = read.delim(filename, header=F, stringsAsFactors=F)
                names(df) = c('value','count', 'state')
                df$sample = subtrack$sample[i]
                df$type = subtrack$type[i]
                # Round to nearest tenth:
                df$rdval = round(df$value,1)
                df = aggregate(count ~ rdval + sample + type + state, df, sum)
                # Split by instate:
                subdf = filter(df, state == 'instate')
                subdf$cs = cumsum(subdf$count)
                subdf$total = subdf$cs[length(subdf$cs)]
                # and out of state:
                outdf = filter(df, state == 'outside')
                outdf$cs = cumsum(outdf$count)
                outdf$total = outdf$cs[length(outdf$cs)]
                df = rbind(subdf, outdf)
                markdf = rbind(markdf, df)
            } # else {print(filename)}
        }
        close(pb)
        print(dim(markdf))
        # Write out the mark table:
        write.table(markdf, gzfile(markfile), quote=F, sep="\t", row.names=F)
    }

    if (!is.null(markdf)){
        markdf = filter(markdf, sample %in% cellorder)
        gplot = ggplot(markdf, aes(rdval, count / total, color=sample)) + 
            geom_line(alpha=.2) + 
            facet_grid(state~type) + 
            geom_vline(xintercept=2, color='blue',lty='dashed') + 
            theme_bw() + 
            scale_y_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) + 
            labs(y='Percent of Bins', x='Signal Value in Bin', title=paste(mark, 'Signal Distribution')) +
            theme(legend.position='none')
        ggsave(paste0(imgpref, 'distr_xraw_', mark, '.png'), gplot, width=8, height=5, dpi=250, units='in')

        gplot = ggplot(markdf, aes(rdval, count / total, color=sample)) + 
            geom_line(alpha=.2) + 
            facet_grid(state~type) + 
            scale_x_log10() + 
            geom_vline(xintercept=2, color='blue',lty='dashed') + 
            theme_bw() + 
            scale_y_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) + 
            labs(y='Percent of Bins', x='Signal Value in Bin', title=paste(mark, 'Signal Distribution')) +
            theme(legend.position='none')
        ggsave(paste0(imgpref, 'distr_xlog_', mark, '.png'), gplot, width=8, height=5, dpi=250, units='in')

        gplot = ggplot(markdf, aes(rdval, count, lty=sample, color=state)) + 
            geom_line(alpha=.2) + 
            facet_wrap(~type, nrow=2) + 
            scale_x_log10() + 
            geom_vline(xintercept=2, color='blue',lty='dashed') + 
            theme_bw() + 
            scale_y_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) + 
            labs(y='Number of Bins', x='Signal Value in Bin', title=paste(mark, 'Signal Distribution')) +
            theme(legend.position='none')
        ggsave(paste0(imgpref, 'counts_xlog_', mark, '.png'), gplot, width=8, height=5, dpi=250, units='in')

        # Plot cumulative.
        gplot = ggplot(markdf, aes(rdval, cs / total, color=sample)) + 
            geom_line(alpha=.2) + 
            facet_wrap(state~type) + 
            scale_x_log10() + 
            geom_vline(xintercept=2, color='blue',lty='dashed') + 
            theme_bw() + 
            labs(y='Percent of Bins', x='Signal Value in Bin', title=paste(mark, 'Signal Distribution')) +
            theme(legend.position='none')
        ggsave(paste0(imgpref, 'cmdistr_xlog_', mark, '.png'), gplot, width=8, height=5, dpi=250, units='in')

        # Add in the aggregate track overall:
        aggdf = aggregate(count ~ rdval + type + sample, markdf, sum)
        atdf = aggregate(count ~ type + sample, aggdf, sum)
        names(atdf)[3] = 'total'
        aggdf = merge(atdf, aggdf)

        # Complete aggregate:
        totdf = aggregate(count ~ rdval + type, markdf, sum)
        ttdf = aggregate(count ~ type , totdf, sum)
        names(ttdf)[2] = 'total'
        totdf = merge(ttdf, totdf)
        totdf$sample = 'total'
        subdf = filter(totdf, type == 'Imputed')
        subdf$cs = cumsum(subdf$count)
        # and out of state:
        outdf = filter(totdf, type == 'Observed')
        outdf$cs = cumsum(outdf$count)
        totdf = rbind(subdf, outdf)
        totdf$cdf = totdf$cs / totdf$total

        # Quantile of aggregate/average observed track at cutoff
        qcut = totdf$cdf[totdf$rdval == 2.0 & totdf$type == 'Observed']
        subdf = totdf[totdf$type == 'Imputed',]
        rn = c(tail(which(subdf$cdf <= qcut),1), which(subdf$cdf >= qcut)[1])
        # NOTE: May want to revisit, with binsignal from ChromHMM?
        print(subdf[rn,])
        # Take in conservative fashion:
        cutoff = subdf$rdval[rn[2]]
        cutdf = rbind(cutdf, data.frame(mark = mark, cutoff = cutoff))

        # Look at range of cutoffs if we were to cut each off individually:
        # sdf = c()
        # for (samp in unique(aggdf$sample)){
        #     subdf = filter(aggdf, type == 'Imputed' & sample == samp)
        #     subdf$cs = cumsum(subdf$count)
        #     subdf$cdf = subdf$cs / subdf$total
        #     idx = tail(which(subdf$cdf <= qcut),1)
        #     cutoff = subdf$rdval[idx]
        #     sdf = rbind(sdf, data.frame(sample = samp, cutoff = cutoff))
        # }

        gplot = ggplot(aggdf, aes(rdval, count / total, color=sample)) + 
            geom_line(alpha=.2) + 
            geom_line(data=totdf, aes(rdval, count / total, color=sample), lwd=1, color='black') + 
            facet_wrap(~type, nrow=2) + 
            scale_x_log10() + 
            geom_vline(xintercept=2, color='blue',lty='dashed') + 
            theme_bw() + 
            scale_y_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) + 
            labs(y='Percent of Bins', x='Signal Value in Bin', title=paste(mark, 'Signal Distribution')) +
            theme(legend.position='none')
        ggsave(paste0(imgpref, 'distragg_xlog_', mark, '.png'), gplot, width=8, height=5, dpi=250, units='in')
        # cant just pick 10^-4 because that is a raw top intensity regions cutoff.

        # Only observed:
        gplot = ggplot(filter(aggdf, type=='Observed'), aes(rdval, count / total, color=sample)) + 
            geom_line(alpha=.2) + 
            geom_line(data=filter(totdf, type=='Observed'), aes(rdval, count / total, color=sample), lwd=1, color='black') + 
            scale_x_log10() + 
            geom_vline(xintercept=2, color='blue',lty='dashed') + 
            theme_bw() + 
            scale_y_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x))) + 
            labs(y='Percent of Bins', x='Signal Value in Bin', title=paste(mark, 'Signal Distribution')) +
            theme(legend.position='none')
        ggsave(paste0(imgpref, 'obs_distragg_xlog_', mark, '.png'), gplot, width=8, height=5, dpi=250, units='in')
    }
}
print(cutdf)

write.table(cutdf, 'ChromImpute/binarization_cutoffs_match_avg_distr.tsv', quote=F, sep="\t", row.names=F)
