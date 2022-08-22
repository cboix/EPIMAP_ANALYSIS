#!/usr/bin/R
# ------------------------------
# Compare the variability tracks
# ------------------------------
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
imgdir = paste0(img, "variance/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'cov_')

# -----------------------
# Read in variance files:
# -----------------------
cvdir = 'ChromImpute/coeffv/'
# Mark-specific:
mark = 'H3K27me3'
v1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_coeffv.tsv.gz')), header=F)[,1]
vm1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_mean.tsv.gz')), header=F)[,1]
vs1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark ,'_chr19_std.tsv.gz')), header=F)[,1]
vn1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_npeak.tsv.gz')), header=F)[,1]
# (Imputed) mark-specific:
i1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_coeffv.tsv.gz')), header=F)[,1]
im1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_mean.tsv.gz')), header=F)[,1]
is1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_std.tsv.gz')), header=F)[,1]
in1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_npeak.tsv.gz')), header=F)[,1]
# Sample-specific:
m1 = read.delim(gzfile(paste0(cvdir, 'chr19_impute_BSS00001_H3K27me3.wig.gz')), skip=2, header=F)[,1]
c1 = read.delim(gzfile(paste0(cvdir, 'chr19_impute_BSS00001_H3K27me3_coeffv.wig.gz')), skip=2, header=F)[,1]

v1[is.na(v1)] = 0
i1[is.na(i1)] = 0
c1[is.na(c1)] = 0

cor(v1, c1)
cor(v1, m1)
cor(vm1, c1)
cor(vm1, m1)
cor(vs1, c1)
cor(vs1, m1)

x = 1:length(v1)
ind = 1e4:6e4

png(paste0(imgpref, 'chr19_H3K27me3_example_tracks.png'), units='in', res=400, width=8, height=4)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=c(sp,sp,sp,sp))
layout(matrix(1:8, ncol=1), heights=rep(1,8))
plot(x[ind], m1[ind], type='l', bty='n', col='black', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Sample Track Mean', font=2, adj=0)
# lines(x[ind], v1[ind], type='l', bty='n', col='red')
plot(x[ind], c1[ind], type='l', bty='n', col='red', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Sample Track Coefficient of Variation', col='red', font=2, adj=0)
# Imputed
plot(x[ind], im1[ind], type='l', bty='n', col='grey50', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Mark Mean (imputed)', col='grey50', font=2, adj=0)
plot(x[ind], i1[ind], type='l', bty='n', col='indianred', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Mark CV (imputed)', col='indianred',font=2, adj=0)
plot(x[ind], is1[ind], type='l', bty='n', col='orange', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Mark std. (imputed)', col='orange', font=2, adj=0)
# Observed
plot(x[ind], vm1[ind], type='l', bty='n', col='grey50', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Mark Mean (observed)', col='grey50', font=2, adj=0)
plot(x[ind], v1[ind], type='l', bty='n', col='royalblue', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Mark CV (observed)', col='royalblue',font=2, adj=0)
plot(x[ind], vs1[ind], type='l', bty='n', col='slateblue', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Mark std. (observed)', col='slateblue', font=2, adj=0)
dev.off()


png(paste0(imgpref, 'chr19_H3K27me3_example_tracks_cvonly.png'), units='in', res=400, width=8, height=2)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=c(sp,sp,sp,sp))
layout(matrix(1:4, ncol=1), heights=rep(1,4))
plot(x[ind], m1[ind], type='l', bty='n', col='black', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Sample Track Mean', font=2, adj=0)
plot(x[ind], c1[ind], type='l', bty='n', col='red', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Sample Track Coefficient of Variation', col='red', font=2, adj=0)
# Imputed
plot(x[ind], i1[ind], type='l', bty='n', col='royalblue', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Mark CV (imputed)', col='indianred',font=2, adj=0)
# Observed:
plot(x[ind], v1[ind], type='l', bty='n', col='indianred', axes=F)
text(x=parpos(1,-.05), y=parpos(2,-.85), 'Mark CV (observed)', col='royalblue',font=2, adj=0)
dev.off()


# Separate high calls from low in mean --> distr of cv
# Then stratify high calls by the cv --> compare to the cv of the mark overall.
mind = m1 > 2

png(paste0(imgpref, 'chr19_H3K27me3_example_distributions.png'), units='in', res=400, width=6, height=6)
par(xaxs='i')
par(yaxs='i')
sp = 0.1
par(mar=c(2,2,sp,sp))
layout(matrix(1:6, ncol=2, byrow=T), heights=rep(1,3))
hist(c1[mind])
hist(c1[!mind])
hist(v1[mind])
hist(v1[!mind])
hist(i1[mind])
hist(i1[!mind])
dev.off()

cvdf = rbind(data.frame(cv=c1, call=mind, op=vm1 > 2, ip = vm1 > 2, type='Sample'),
             data.frame(cv=v1, call=mind, op=vm1 > 2, ip = vm1 > 2, type='Observed'),
             data.frame(cv=i1, call=mind, op=vm1 > 2, ip = vm1 > 2, type='Imputed'))

cvdf$peak = 'No'
cvdf$peak[cvdf$call == TRUE] = 'Yes'

gplot = ggplot(cvdf, aes(type, cv, fill=peak)) + 
    labs(x='Dataset', y='Coefficient of Variation')+ 
    geom_boxplot() + theme_pubr() + scale_fill_manual(values=c('grey75','grey25'))
ggsave(paste0(imgpref, 'chr19_H3K27me3_example_boxplots.png'), gplot, units='in', dpi=400, width=4, height=6)

sub.cvdf = cvdf[cvdf$call==TRUE,]
gplot = ggplot(sub.cvdf, aes(type, cv, fill=op)) + 
    labs(x='Dataset', y='Coefficient of Variation', title='Peaks matching peaks in the observed')+ 
    geom_boxplot() + theme_pubr() + scale_fill_manual(values=c('grey75','grey25'))
ggsave(paste0(imgpref, 'chr19_H3K27me3_example_boxplots_op.png'), gplot, units='in', dpi=400, width=4, height=6)


# Show presence of locations with high coeff var in obs + imp in mark and present as peaks, but not in mean?
vind = v1 > quantile(v1, .99)
iind = i1 > quantile(i1, .99)
sum(iind * vind) / sum(iind)

png(paste0(imgpref, 'top_cov.png'), units='in', res=400, width=6, height=4)
cind = which(c1 > 1.25)
hist(v1[cind])
dev.off()




# -----------------------------------------
# Make scatterplots for pairs of variables:
# -----------------------------------------
png(paste0(imgpref, 'mean_cv_scatter.png'), units='in', res=400, width=6, height=6)
sp = 0.1
par(mar=c(4,4,sp,sp))
plot(m1, c1, pch=1, cex=.07, col=tsp.col('black', .1),
     bty='n', xlab='Imputed Track (Mean)', ylab='Imputation Coefficient of Variation')
dev.off()

png(paste0(imgpref, 'obs_mean_cv_scatter.png'), units='in', res=400, width=6, height=6)
par(mar=c(4,4,sp,sp))
plot(vm1, v1, pch=1, cex=.07, col=tsp.col('black', .1),
     bty='n', xlab='Mean (Observed Tracks)', ylab='Coefficient of Variation (Observed Tracks)')
dev.off()

png(paste0(imgpref, 'imp_mean_cv_scatter.png'), units='in', res=400, width=6, height=6)
par(mar=c(4,4,sp,sp))
plot(im1, i1, pch=1, cex=.07, col=tsp.col('black', .1),
     bty='n', xlab='Mean (Imputed Tracks)', ylab='Coefficient of Variation (Imputed Tracks)')
dev.off()

png(paste0(imgpref, 'obs_mean_imp_mean_scatter.png'), units='in', res=400, width=6, height=6)
par(mar=c(4,4,sp,sp))
plot(vm1, im1, pch=1, cex=.07, col=tsp.col('black', .1),
     bty='n', xlab='Mean (Observed Tracks)', ylab='Mean (Imputed Tracks)')
dev.off()


png(paste0(imgpref, 'obs_cv_imp_cv_scatter.png'), units='in', res=400, width=6, height=6)
par(mar=c(4,4,sp,sp))
plot(v1, i1, pch=1, cex=.07, col=tsp.col('black', .1),
     bty='n', xlab='Coefficient of Variation (Observed Tracks)', ylab='Coefficient of Variation (Imputed Tracks)')
dev.off()

png(paste0(imgpref, 'obs_np_imp_np_scatter.png'), units='in', res=400, width=6, height=6)
par(mar=c(4,4,sp,sp))
plot(vn1, in1, pch=1, cex=.07, col=tsp.col('black', .1),
     bty='n', xlab='# Sample with Peaks (Observed Tracks)', ylab='# Sample with Peaks (Imputed Tracks)')
dev.off()


df = data.frame(obs.mean=vm1, imp.mean=im1)
fit = lm(obs.mean ~ imp.mean, df)
summary(fit)

df = data.frame(obs.cv=v1, imp.cv=i1)
fit = lm(obs.cv ~ imp.cv, df)
summary(fit)

cor(v1, c1)

df = data.frame(obs.np=vn1, imp.np=in1)
fit = lm(obs.np ~ imp.np, df)
summary(fit)



# ----------------
# Load in H3K27ac: 
# TODO: Collapse later to run all marks.
# ----------------
fnames = list.files(cvdir, pattern='obs.*coeffv.tsv.gz')
marks = sub("_chr.*tsv.gz", "", sub("obs_","", fnames)) 

rdf = c()
all.fdf = c()
cutfrac.rda.file = paste0(cvdir, 'intext_fracrecovery_stats.Rda')
if (!file.exists(cutfrac.rda.file)){
    for (mark in marks){
        print(mark)
        v1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_coeffv.tsv.gz')), header=F)[,1]
        vm1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_mean.tsv.gz')), header=F)[,1]
        vn1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_npeak.tsv.gz')), header=F)[,1]
        i1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_coeffv.tsv.gz')), header=F)[,1]
        im1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_mean.tsv.gz')), header=F)[,1]
        in1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_npeak.tsv.gz')), header=F)[,1]
        v1[is.na(v1)] = 0
        i1[is.na(i1)] = 0

        plot.scatter = FALSE
        if (plot.scatter) {
            png(paste0(imgpref, 'obs_mean_imp_mean_scatter_', mark, '.png'), units='in', res=400, width=6, height=6)
            par(mar=c(4,4,sp,sp))
            plot(vm1, im1, pch=1, cex=.07, col=tsp.col('black', .1),
                 bty='n', xlab='Mean (Observed Tracks)', ylab='Mean (Imputed Tracks)')
            mtext(mark, line=-1)
            dev.off()

            png(paste0(imgpref, 'obs_cv_imp_cv_scatter_', mark, '.png'), units='in', res=400, width=6, height=6)
            par(mar=c(4,4,sp,sp))
            plot(v1, i1, pch=1, cex=.07, col=tsp.col('black', .1),
                 bty='n', xlab='Coefficient of Variation (Observed Tracks)', ylab='Coefficient of Variation (Imputed Tracks)')
            mtext(mark, line=-1)
            dev.off()

            png(paste0(imgpref, 'obs_np_imp_np_scatter_', mark, '.png'), units='in', res=400, width=6, height=6)
            par(mar=c(4,4,sp,sp))
            plot(vn1, in1, pch=1, cex=.07, col=tsp.col('black', .1),
                 bty='n', xlab='# Sample with Peaks (Observed Tracks)', ylab='# Sample with Peaks (Imputed Tracks)')
            mtext(mark, line=-1)
            dev.off()

            png(paste0(imgpref, 'obs_np_obs_cv_scatter_', mark, '.png'), units='in', res=400, width=6, height=6)
            par(mar=c(4,4,sp,sp))
            plot(vn1, v1, pch=1, cex=.07, col=tsp.col('black', .1),
                 bty='n', xlab='# Sample with Peaks (Observed Tracks)', ylab='Coefficient of variation (Observed Tracks)')
            mtext(mark, line=-1)
            dev.off()
        }

        df = data.frame(obs.mean=vm1, imp.mean=im1)
        fit = lm(obs.mean ~ imp.mean, df)
        rdf = rbind(rdf,  data.frame(mark=mark, adjr=summary(fit)$adj.r.squared, type='mean'))

        df = data.frame(obs.cv=v1, imp.cv=i1)
        fit = lm(obs.cv ~ imp.cv, df)
        rdf = rbind(rdf, data.frame(mark=mark, adjr=summary(fit)$adj.r.squared, type='cv'))

        df = data.frame(obs.np=vn1, imp.np=in1)
        fit = lm(obs.np ~ imp.np, df)
        rdf = rbind(rdf, data.frame(mark=mark, adjr=summary(fit)$adj.r.squared, type='np'))

        totdf = agg.rename(imp.np ~ obs.np, df, length, 'tot.loc')
        # fracdf = agg.rename(imp.np ~ obs.np, df, function(x){sum(x > 1)}, 'nimp')
        # fracdf = agg.rename(imp.np ~ obs.np, df, function(x){sum(x > 10)}, 'nimp')
        fracdf = agg.rename(imp.np ~ obs.np, df, function(x){sum(x > 4)}, 'nimp')
        fracdf = merge(fracdf, totdf)
        fracdf$frac = fracdf$nimp / fracdf$tot.loc
        fracdf$mark = mark
        fracdf$scaled.obs = fracdf$obs.np / max(fracdf$obs.np)
        # plot(fracdf$obs.np, fracdf$frac, type='l')
        all.fdf = rbind(all.fdf, fracdf)
    }
    save(all.fdf, rdf, file=cutfrac.rda.file)
} else { 
    load(cutfrac.rda.file)
}

rwdf = spread(rdf, type, adjr)
rwdf = rwdf[order(rwdf$cv, decreasing=T),]

colnames(markdef)[2] = 'mark.type'
rdf = merge(rdf, markdef, all.x=TRUE)
rdf$mark = factor(rdf$mark, levels=rev(rwdf$mark))
rdf$Type = 'Mean Signal'
rdf$Type[rdf$type == 'cv'] = 'Coeff. of Variation'
rdf$Type[rdf$type == 'np'] = 'Number of samples with peak'
rdf$Type = factor(rdf$Type, levels=c('Mean Signal','Number of samples with peak','Coeff. of Variation'))

gplot = ggplot(rdf, aes(mark, adjr)) + 
    facet_wrap(~Type, nrow=1) + 
    geom_bar(stat='identity', alpha=0.5, fill='grey70') + 
    geom_text(data=rdf, aes(y=.025, x=mark, label=paste0(round(100 * adjr,1), '%')), hjust=0) + 
    lims(y=c(0,1)) + 
    labs(x='Assay/Mark', y='Adjusted R^2 for predicting observed statistic from imputed statistic') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    theme_pubr() + coord_flip()
ggsave(paste0(imgpref, 'chr19_adjr_reg_statistics.png'), gplot, units='in', dpi=400, width=7, height=3)
ggsave(paste0(imgpref, 'chr19_adjr_reg_statistics.pdf'), gplot, units='in', dpi=400, width=7, height=3)

aggregate(adjr ~ Type + mark.type, rdf, mean)

# Make a latex table pdf:
hp = as_hux(rwdf, add_colnames=TRUE)
position(hp) = "left"
top_border(hp)[1,]=2
bottom_border(hp)[1,]=1
bottom_border(hp)[nrow(hp),]=2
width(hp) = .75
quick_pdf(hp, file=paste0(imgpref, "rsq_impobs_chr19_statistics.pdf"))



rwdf = spread(rdf, type, adjr)
rwdf = rwdf[order(rwdf$np, decreasing=T),]
colnames(rwdf) = c('Mark','CoV','Mean','NPeak')

# Make a latex table pdf:
hp = as_hux(rwdf, add_colnames=TRUE)
position(hp) = "left"
top_border(hp)[1,]=2
bottom_border(hp)[1,]=1
bottom_border(hp)[nrow(hp),]=2
width(hp) = .75
quick_pdf(hp, file=paste0(imgpref, "rsq_impobs_chr19_statistics_bynp.pdf"))


# --------------------------------
# Quantification of peak recovery:
# --------------------------------
gplot = ggplot(all.fdf, aes(scaled.obs, frac, color=mark)) + 
    geom_line() + theme_pubr() + 
    scale_x_continuous(labels=scales::percent, expand=c(0,0)) + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    labs(x='Percent of Observed Samples where location is a peak', 
         y='Percent of locations that are peaks (val > 2)\nin more than one imputed samples') + 
    scale_color_manual(values=brewer.pal(12,'Paired')) + 
    theme(legend.position=c(.5,.9))
ggsave(paste0(imgpref, 'chr19_peak_recovery_quantification.png'), gplot, units='in', dpi=400, width=8, height=4)
ggsave(paste0(imgpref, 'chr19_peak_recovery_quantification.pdf'), gplot, units='in', dpi=400, width=8, height=4)


# Point at which they cross 75% recovery
thresh = 0.75
crossdf = merge(all.fdf, aggregate(frac ~ mark, all.fdf, function(x){x[x > thresh][1] }))
crossdf = crossdf[,c('mark','frac','scaled.obs')]
crossdf = crossdf[order(crossdf$scaled.obs, decreasing=F),]
crossdf$pct.obs = paste0(round(crossdf$scaled.obs * 100,1), '%')
crossdf$mark = factor(crossdf$mark, levels=crossdf$mark)

gplot = ggplot(crossdf, aes(mark, scaled.obs)) + 
    geom_bar(stat='identity', alpha=0.5, fill='grey70') + 
    geom_text(data=crossdf, aes(y=.005, x=mark, label=pct.obs), hjust=0) + 
    labs(x='Assay/Mark', y='Percent of observed samples that must have peak to\nhave 75% peak recovery (in more than one imputed)') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    theme_pubr() + coord_flip()
ggsave(paste0(imgpref, 'chr19_pkrecovery_cutoff_statistics.png'), gplot, units='in', dpi=400, width=5.5, height=3)
ggsave(paste0(imgpref, 'chr19_pkrecovery_cutoff_statistics.pdf'), gplot, width=5.5, height=3)


# Make a latex table pdf:
names(crossdf) = c('Mark','Frac. Imp','Pct. Obs')
hp = as_hux(crossdf, add_colnames=TRUE)
position(hp) = "left"
top_border(hp)[1,]=2
bottom_border(hp)[1,]=1
bottom_border(hp)[nrow(hp),]=2
width(hp) = .75
quick_pdf(hp, file=paste0(imgpref, "peak_recovery_impobs_chr19_statistics.pdf"))


# --------------------------------------------------------
# Quantification of peak recovery only for punctate marks:
# --------------------------------------------------------
gplot = ggplot(all.fdf, aes(scaled.obs, frac, color=mark)) + 
    geom_line() + theme_pubr() + 
    scale_x_continuous(labels=scales::percent, expand=c(0,0)) + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    labs(x='Percent of Observed Samples where location is a peak', 
         y='Percent of locations that are peaks (val > 2)\nin more than one imputed samples') + 
    scale_color_manual(values=brewer.pal(12,'Paired')) + 
    theme(legend.position=c(.5,.9))
ggsave(paste0(imgpref, 'chr19_peak_recovery_quantification.png'), gplot, units='in', dpi=400, width=8, height=4)
ggsave(paste0(imgpref, 'chr19_peak_recovery_quantification.pdf'), gplot, units='in', dpi=400, width=8, height=4)


# Point at which they cross 75% recovery
thresh = 0.75
crossdf = merge(all.fdf, aggregate(frac ~ mark, all.fdf, function(x){x[x > thresh][1] }))
crossdf = crossdf[,c('mark','frac','scaled.obs')]
crossdf = crossdf[order(crossdf$scaled.obs, decreasing=F),]
crossdf$pct.obs = paste0(round(crossdf$scaled.obs * 100,1), '%')
crossdf$mark = factor(crossdf$mark, levels=crossdf$mark)

gplot = ggplot(crossdf, aes(mark, scaled.obs)) + 
    geom_bar(stat='identity', alpha=0.5, fill='grey70') + 
    geom_text(data=crossdf, aes(y=.005, x=mark, label=pct.obs), hjust=0) + 
    labs(x='Assay/Mark', y='Percent of observed samples that must have peak to\nhave 75% peak recovery (in more than one imputed)') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    theme_pubr() + coord_flip()
ggsave(paste0(imgpref, 'chr19_pkrecovery_cutoff_statistics.png'), gplot, units='in', dpi=400, width=5.5, height=3)
ggsave(paste0(imgpref, 'chr19_pkrecovery_cutoff_statistics.pdf'), gplot, width=5.5, height=3)


# Make a latex table pdf:
names(crossdf) = c('Mark','Frac. Imp','Pct. Obs')
hp = as_hux(crossdf, add_colnames=TRUE)
position(hp) = "left"
top_border(hp)[1,]=2
bottom_border(hp)[1,]=1
bottom_border(hp)[nrow(hp),]=2
width(hp) = .75
quick_pdf(hp, file=paste0(imgpref, "peak_recovery_impobs_chr19_statistics.pdf"))







