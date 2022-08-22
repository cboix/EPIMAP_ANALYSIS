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
fnames = list.files(cvdir, pattern='chr19_impute_.*coeffv.wig.gz')
prefs = sub("_coeffv.wig.gz","", sub("chr19_impute_","", fnames))

# Info table:
fdf = data.frame(file=fnames, prefix=prefs, id=sub("_.*","",prefs), mark=sub(".*_","",prefs))
marks = unique(fdf$mark)
print(aggregate(id ~ mark, fdf, length))
# Kept samp only:
fdf = fdf[fdf$id %in% cellorder,]

# For plotting tracks:
plt.track = function(y, ind, nam='', col='black', ylim=NULL){
    x = 1:length(y)
    plot(x[ind], y[ind], type='l', bty='n', col=col, axes=F, ylim=ylim)
    text(x=parpos(1,-.05), y=parpos(2,-.85), nam, font=2, adj=0, col=col) 
    text(x=parpos(1,-.000), y=c(parpos(2,-0.15), parpos(2,-.85)), cex=.75,ylim, font=1, adj=0, col='black') }

# Per mark, read in NTOP samples:
intext.rda.file = paste0(cvdir, 'coeffv_reg_stats.Rda')
if (!file.exists(intext.rda.file)){
    NTOP = 25
    all.rdf = c()
    for (mark in marks){
        print(mark)
        # Read in overall mark information:
        v1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_coeffv.tsv.gz')), header=F)[,1]
        vm1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_mean.tsv.gz')), header=F)[,1]
        vn1 = read.delim(gzfile(paste0(cvdir, 'obs_', mark, '_chr19_npeak.tsv.gz')), header=F)[,1]
        i1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_coeffv.tsv.gz')), header=F)[,1]
        im1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_mean.tsv.gz')), header=F)[,1]
        in1 = read.delim(gzfile(paste0(cvdir, 'imp_', mark, '_chr19_npeak.tsv.gz')), header=F)[,1]
        v1[is.na(v1)] = 0
        i1[is.na(i1)] = 0

        # Read in NTOP of each mark:
        sfdf = fdf[fdf$mark == mark,]
        set.seed(1)
        sind = sample(1:nrow(sfdf), NTOP)
        sfdf = sfdf[sind,]
        meta[sfdf$id,]
        sfiles = sfdf$file
        cmat = matrix(0, ncol=NTOP, nrow=length(i1))
        snam = c()
        for (i in 1:length(sfiles)){
            fname = sfiles[i]
            spref = sfdf[sfdf$file == fname,'prefix']
            cmat[,i] = read.delim(gzfile(paste0(cvdir, fname)), header=F, skip=2)[,1]
            snam = c(snam, spref)
        }
        colnames(cmat) = snam

        # Compare tracks to the overall CoV:
        ci = cor(i1, cmat)
        co = cor(v1, cmat)
        cim = cor(im1, cmat)
        com = cor(vm1, cmat)
        cin = cor(in1, cmat)
        con = cor(vn1, cmat)

        cdf = data.frame(cmat)
        cdf$obs.cov = v1
        cdf$obs.mean = vm1
        cdf$obs.npeak = vn1
        cdf$imp.cov = i1
        cdf$imp.mean = im1
        cdf$imp.npeak = in1

        rdf = c()
        for (sample in snam){
            print(sample)
            obs.fit = lm(as.formula(paste0(sample, '~ obs.cov + obs.mean + obs.npeak')), cdf)
            imp.fit = lm(as.formula(paste0(sample, '~ imp.cov + imp.mean + imp.npeak')), cdf)
            rdf = rbind(rdf, data.frame(samp=sample, obs.rsq=summary(obs.fit)$adj.r, imp.rsq=summary(imp.fit)$adj.r))
        }
        print(rdf)
        rdf$mark = mark
        all.rdf = rbind(all.rdf, rdf)

        # Model only in locations with peaks:
        # Actually performs worse.
        rsdf = c()
        oind = which(cdf$obs.npeak > 5)
        iind = which(cdf$imp.npeak > 5)
        for (sample in snam){
            print(sample)
            obs.fit = lm(as.formula(paste0(sample, '~ obs.cov + obs.mean + obs.npeak')), cdf[oind,])
            imp.fit = lm(as.formula(paste0(sample, '~ imp.cov + imp.mean + imp.npeak')), cdf[iind,])
            rsdf = rbind(rsdf, data.frame(samp=sample, obs.rsq=summary(obs.fit)$adj.r, imp.rsq=summary(imp.fit)$adj.r))
        }
        print(rsdf)

        plot.tracks = FALSE
        if (plot.tracks){
            # Plot each of these:
            print("Plotting tracks")
            ind = 1e4:6e4
            tm = 3
            hm = which(im1 > tm & vm1 > tm)
            hm = hm[hm %in% ind]
            tc = 2
            cm = which(i1 > tc & v1 > tc)
            cm = cm[cm %in% ind]

            ylim = c(0, max(cmat[ind,]))
            png(paste0(imgpref, 'chr19_',mark, '_example_coeffv_tracks.png'), units='in', res=400, width=8, height=5)
            par(xaxs='i')
            par(yaxs='i')
            sp = 0.1
            par(mar=c(sp,sp,sp,sp))
            layout(matrix(1:(NTOP + 4), ncol=1), heights=rep(1,(NTOP + 4)))
            for (i in 1:ncol(cmat)){
                spref = colnames(cmat)[i]
                prefstr = sub("_", " ",spref)
                plt.track(cmat[,spref], ind, nam=prefstr, col=ifelse(i %% 2 == 1, 'grey40', 'grey60'), ylim=ylim)
                pad = 0.5
                rect(xleft=hm - pad, xright=hm + pad, ybottom=0, ytop=ylim[2], col=tsp.col('goldenrod1'), border=NA)
                if (length(cm) > 0){
                    cpad = 5
                    rect(xleft=cm - cpad, xright=cm + cpad, ybottom=0, ytop=ylim[2], col=tsp.col('indianred'), border=NA)
                }
            }
            # Imputed
            plt.track(vm1, ind, nam='Overall Mean (Observed)', col='slateblue')
            plt.track(v1, ind, nam='Overall CoV (Observed)', col='royalblue')
            plt.track(im1, ind, nam='Overall Mean (Imputed)', col='firebrick')
            plt.track(i1, ind, nam='Overall CoV (Imputed)', col='indianred')
            dev.off()
        }
    }
    save(all.rdf, file=intext.rda.file)
} else {
    load(intext.rda.file)
}



# --------------------------------
# Quantification of peak recovery:
# --------------------------------
rwdf = gather(all.rdf, type, adjr, -samp, -mark)
rwdf$Type = 'Observed'
rwdf$Type[rwdf$type == 'imp.rsq'] = 'Imputed'
aggregate(adjr ~ mark, rwdf, mean)

gplot = ggplot(rwdf, aes(mark, adjr, fill=Type)) + 
    geom_boxplot() + theme_pubr() + 
    labs(x='Imputed Histone Mark', y='Adjusted R-squared', 
         title='Adj. Rsq for using data-wide variation metrics\nto model individual coefficient of variation tracks.') + 
    scale_fill_manual(values=brewer.pal(12,'Paired')[c(5,1)], name='') + 
    theme(legend.position=c(0.11,.95))
ggsave(paste0(imgpref, 'internal_cov_model.png'), gplot, units='in', dpi=400, width=6, height=4.5)


# ------------------------------------------------
# Predicting how well external and internal agree:
# ------------------------------------------------
rwdf$id = sub("_.*","", rwdf$samp)
avail.marg = colSums(avail)
rwdf$is.avail = apply(rwdf, 1, function(x){avail[x['mark'], x['id']]})
rwdf$navail = avail.marg[rwdf$id] - rwdf$is.avail
rwdf$is.alone = (rwdf$navail == 1)
# TODO: need to sum actual
rwdf$navail.cut = cut(rwdf$navail-1, breaks=c(0,1,3,10), include.lowest=TRUE, right=FALSE)
rwdf$n27me3 = avail['H3K27me3', rwdf$id] ==1
rwdf$n27ac = avail['H3K27ac', rwdf$id] == 1
rwdf$ndnase = avail['DNase-seq', rwdf$id] == 1

g1 = ggplot(rwdf[rwdf$Type == 'Observed',], aes(navail, adjr)) + 
    facet_wrap(~mark) + 
    geom_point() + theme_pubr() + 
    geom_smooth(method='lm') + 
    scale_y_continuous(expand=c(0,0)) + 
    labs(x='Number of same-sample datasets', y='Adjusted R-squared', 
         title='Adj. Rsq vs. number of same-sample, different-mark datasets') + 
    scale_fill_manual(values=brewer.pal(12,'Paired')[c(5,1)], name='') + 
    theme(legend.position=c(0.11,.95))
ggsave(paste0(imgpref, 'internal_cov_model_bysamesample.png'), g1, units='in', dpi=400, width=6, height=4.5)


gplot = ggplot(rwdf[rwdf$Type == 'Observed',], aes(mark, adjr, fill=is.alone)) + 
    geom_boxplot() + theme_pubr() + 
    labs(x='Imputed Histone Mark', y='Adjusted R-squared', 
         title='Adj. Rsq vs. whether sample has only one informant dataset') + 
    scale_fill_manual(values=brewer.pal(12,'Paired')[c(1,2)], name='Single-expt.\nin sample') + 
    stat_compare_means(hide.ns=TRUE, label='p.format') + 
    theme(legend.position=c(0.11,.75))
ggsave(paste0(imgpref, 'internal_cov_model_alone.png'), gplot, units='in', dpi=400, width=6, height=4.5)

srwdf = rwdf[rwdf$Type == 'Observed' & rwdf$mark %in% c('H3K27me3','H3K36me3'),]
srwdf = srwdf[srwdf$is.alone,]
srwdf = srwdf[order(srwdf$adjr),]
wmarg = apply(avail, 2, function(x){paste(names(which(x == 1)), collapse=', ')})
srwdf$which.avail = wmarg[srwdf$id]

rwdf = cbind(rwdf, data.frame(t(avail[marks,rwdf$id] == 1)))
rlong = gather(rwdf[rwdf$Type == 'Observed', c('mark','adjr','id', marks)],
               amark, avail, -mark, -adjr, -id)

gplot = ggplot(rlong, aes(mark, adjr, fill=avail)) + 
    facet_wrap(~amark) + 
    geom_boxplot() + theme_pubr() + 
    labs(x='Imputed Histone Mark', y='Adjusted R-squared for external-internal variation', 
         title='Adj. Rsq vs. whether sample has observed experiments for each mark') + 
    scale_fill_manual(values=brewer.pal(12,'Paired')[c(1,2)], name='Sample has ChIP-seq experiment for the mark') + 
    theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(paste0(imgpref, 'internal_cov_model_bymark.png'), gplot, units='in', dpi=400, width=10, height=6)





print(rdf)

rwdf = spread(rdf, type, adjr)
rwdf = rwdf[order(rwdf$cv, decreasing=T),]

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


# # --------------------------------
# # Quantification of peak recovery:
# # --------------------------------
# gplot = ggplot(all.fdf, aes(scaled.obs, frac, color=mark)) + 
#     geom_line() + theme_pubr() + 
#     scale_x_continuous(labels=scales::percent, expand=c(0,0)) + 
#     scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
#     labs(x='Percent of Observed Samples where location is a peak', 
#          y='Percent of locations that are peaks (val > 2)\nin more than one imputed samples') + 
#     scale_color_manual(values=brewer.pal(12,'Paired'))
# ggsave(paste0(imgpref, 'chr19_peak_recovery_quantification.png'), gplot, units='in', dpi=400, width=8, height=6)

# # Point at which they cross 75% recovery
# crossdf = merge(all.fdf, aggregate(frac ~ mark, all.fdf, function(x){x[x > 0.75][1] }))
# crossdf = crossdf[,c('mark','frac','scaled.obs')]
# crossdf = crossdf[order(crossdf$scaled.obs, decreasing=F),]
# crossdf$scaled.obs = paste0(round(crossdf$scaled.obs * 100,1), '%')
# names(crossdf) = c('Mark','Frac. Imp','Pct. Obs')

# # Make a latex table pdf:
# hp = as_hux(crossdf, add_colnames=TRUE)
# position(hp) = "left"
# top_border(hp)[1,]=2
# bottom_border(hp)[1,]=1
# bottom_border(hp)[nrow(hp),]=2
# width(hp) = .75
# quick_pdf(hp, file=paste0(imgpref, "peak_recovery_impobs_chr19_statistics.pdf"))






