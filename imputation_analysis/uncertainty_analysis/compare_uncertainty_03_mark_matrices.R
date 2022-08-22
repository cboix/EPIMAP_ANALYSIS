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
library(viridis)
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
lddir = 'linking_data/'
suff = "_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged_stats_df.tsv.gz"
fnames = list.files(lddir, pattern=paste0('.*', suff))
marks = sub(suff,"", fnames)

# Read in overall mark information:
for (mark in marks){
    print(mark)
    # Read in 3.6M * attr dataframe:
    df = read.delim(gzfile(paste0(lddir, mark, suff)), header=T)
    ind = 1:5e5
    pc = 0.1

    g1 = ggplot(df[ind,], aes(nsamp, CoV, fill=log10(..count..))) + 
        geom_hex() + theme_pubr() + 
        scale_x_log10() + scale_fill_viridis() + 
        scale_y_log10() + 
        labs(x='Number of Samples with Peak', y='Coefficient of Variation')

    g2 = ggplot(df[ind,], aes(Max + pc, CoV, fill=log10(..count..))) + 
        geom_hex() + theme_pubr() + 
        scale_x_log10() + scale_fill_viridis() + 
        scale_y_log10() + 
        labs(x='Maximum value', y='Coefficient of Variation')

    g3 = ggplot(df[ind,], aes(Mean + pc, CoV, fill=log10(..count..))) + 
        geom_hex() + theme_pubr() + 
        scale_x_log10() + scale_fill_viridis() + 
        scale_y_log10() + 
        labs(x='Mean value', y='Coefficient of Variation')

    g4 = ggplot(df[ind,], aes(Mean + pc, nsamp, fill=log10(..count..))) + 
        geom_hex() + theme_pubr() + 
        scale_x_log10() + scale_fill_viridis() + 
        labs(x='Mean value', y='Number of Samples with Peak')
        
    g5 = ggplot(df[ind,], aes(Max + pc, nsamp, fill=log10(..count..))) + 
        geom_hex() + theme_pubr() + 
        scale_x_log10() + scale_fill_viridis() + 
        labs(x='Max value', y='Number of Samples with Peak')
        

    garr = ggarrange(g3 + labs(title=mark), g2, g1,
              g4,g5,ncol=3, nrow=2, common.legend=TRUE)
    ggsave(paste0(imgpref, 'markmatrix_stats_', mark, '.png'), garr, units='in', dpi=400, width=12, height=7)


    rn = log10(range(df$nsamp + 1))
    breaks = seq(rn[1], rn[2], length.out = 8)
    breaks = 10^breaks - 1
    df$nsamp.cut = cut(df$nsamp, breaks=breaks, include.lowest=TRUE)

    qdf = aggregate(CoV ~ nsamp.cut, df, function(x){quantile(x, probs=c(.25, .75))})
    qdf$b1 = breaks[-length(breaks)]
    qdf$b2 = breaks[-1]
    qdf$c1 = qdf$CoV[,1]
    qdf$c2 = qdf$CoV[,2]

    g1 = ggplot(df[ind,], aes(nsamp, CoV)) + 
        geom_hex(aes(fill=log10(..count..))) + theme_pubr() + 
        geom_rect(data=qdf, aes(x=b1, y=c1, xmin=b1, xmax=b2, ymin=c1, ymax=c2), fill=NA, col='red') + 
        scale_fill_viridis() + 
        scale_x_log10() +
        scale_y_log10() + 
        labs(x='Number of Samples with Peak', y='Coefficient of Variation')
    ggsave(paste0(imgpref, 'markmatrix_covrange_', mark, '.png'), g1, units='in', dpi=400, width=4, height=4)
}




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
NTOP = 10
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




    # png(paste0(imgpref, 'obs_np_obs_cv_scatter_', mark, '.png'), units='in', res=400, width=6, height=6)
    # par(mar=c(4,4,sp,sp))
    # plot(vn1, v1, pch=1, cex=.07, col=tsp.col('black', .1),
    #      bty='n', xlab='# Sample with Peaks (Observed Tracks)', ylab='Coefficient of variation (Observed Tracks)')
    # mtext(mark, line=-1)
    # dev.off()
}



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


# --------------------------------
# Quantification of peak recovery:
# --------------------------------
gplot = ggplot(all.fdf, aes(scaled.obs, frac, color=mark)) + 
    geom_line() + theme_pubr() + 
    scale_x_continuous(labels=scales::percent, expand=c(0,0)) + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
    labs(x='Percent of Observed Samples where location is a peak', 
         y='Percent of locations that are peaks (val > 2)\nin more than one imputed samples') + 
    scale_color_manual(values=brewer.pal(12,'Paired'))
ggsave(paste0(imgpref, 'chr19_peak_recovery_quantification.png'), gplot, units='in', dpi=400, width=8, height=6)

# Point at which they cross 75% recovery
crossdf = merge(all.fdf, aggregate(frac ~ mark, all.fdf, function(x){x[x > 0.75][1] }))
crossdf = crossdf[,c('mark','frac','scaled.obs')]
crossdf = crossdf[order(crossdf$scaled.obs, decreasing=F),]
crossdf$scaled.obs = paste0(round(crossdf$scaled.obs * 100,1), '%')
names(crossdf) = c('Mark','Frac. Imp','Pct. Obs')

# Make a latex table pdf:
hp = as_hux(crossdf, add_colnames=TRUE)
position(hp) = "left"
top_border(hp)[1,]=2
bottom_border(hp)[1,]=1
bottom_border(hp)[nrow(hp),]=2
width(hp) = .75
quick_pdf(hp, file=paste0(imgpref, "peak_recovery_impobs_chr19_statistics.pdf"))






