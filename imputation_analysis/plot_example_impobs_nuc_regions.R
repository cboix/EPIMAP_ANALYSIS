#!/usr/bin/R
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
library(tidyr)
library(dplyr)
library(ggridges)
library(ggplot2)
library(viridis)
library(NNLM)
library(WaveletComp)

smooth.wavelet = function(x){
    data = data.frame(val=x)
    my.w <- analyze.wavelet(data, "val", loess.span = 0, 
                            dt = 25, dj = 1/50, lowerPeriod = 50, upperPeriod = 1000,
                            make.pval=F, n.sim = 10, verbose=F)
    my.rec <- reconstruct(my.w, plot.rec=FALSE, verbose=F)
    rec <- my.rec$series$val.r
    return(rec)
}

localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) { y <- y[-1] }
    return(y)
}

mkmat = function(mark){
    df = read.delim(paste0('example_regions/concat_tmp_', mark, '.mtx'), header=F, stringsAsFactors=F) 
    names(df) = c('pos','val','cell','set')
    odf = filter(df, set == 'o')
    idf = filter(df, set == 'i')
    df2 = df
    qm = quantile(df$val, 1 - 5e-4)
    df2$val[df2$val > qm]  = qm
    cells = sort(unique(df$cell))
    pos = sort(unique(df$pos))
    NC = length(cells)
    tmp = matrix(NA, nrow=length(pos), ncol=NC, dimnames=list(pos, cells))
    # Observed data:
    spread(odf[,c('pos','val','cell')], cell, val) -> wide
    mat = as.matrix(wide[,-1])
    rownames(mat) = wide$pos
    qm = quantile(mat, 1 - 5e-4)
    mat[mat > qm] = qm
    omat = tmp
    omat[rownames(mat), colnames(mat)] = mat
    # Imputed data:
    spread(idf[,c('pos','val','cell')], cell, val) -> wide
    mat = as.matrix(wide[,-1])
    rownames(mat) = wide$pos
    qm = quantile(mat, 1 - 5e-4)
    mat[mat > qm] = qm
    imat = tmp
    imat[rownames(mat), colnames(mat)] = mat
    return(list(imat, omat, df2)) 
}


# mark = 'H3K27ac'
mark = 'DNase-seq'
ll = mkmat(mark)
imat = ll[[1]]
omat = ll[[2]]
df2 = ll[[3]]
imat[is.na(imat)] = -.5
omat[is.na(omat)] = -.5


# mark = 'H3K27ac'
mdf = c()
marks = c('H3K27ac', 'H3K4me1','H3K4me3', 'H3K36me3','H3K9me3', 'DNase-seq')
for (m in marks){
    ll = mkmat(m)
    df = ll[[3]]
    mdf = rbind(mdf, data.frame(df, mark=m))
}

cells = sort(unique(df$cell))
pos = sort(unique(df$pos))
tmp = matrix(NA, nrow=length(pos), ncol=length(cells),
             dimnames=list(pos, cells))

# Aggregating doesnt help us see if nucleosomes agree.
amdf = aggregate(val ~ pos + cell + set, mdf, mean)
wide = spread(filter(amdf, set == 'o')[,c('pos','val','cell')], cell, val)
mat = as.matrix(wide[,-1])
rownames(mat) = wide$pos
qm = quantile(mat, 1 - 5e-4)
mat[mat > qm] = qm
momat = tmp
momat[rownames(mat), colnames(mat)] = mat


# For each track, estimate pks:
NP = 15
png(paste0(img, 'impobs_all_n',NP, '.png'), res=300, units='in', width=10, height=NP)
layout(matrix(c(1:NP), NP,1), T)
par(xaxs='i')
par(mar = rep(0,4))
subdf = filter(mdf, set == 'o')
counts = aggregate(mark ~ cell, subdf, function(x){length(unique(x))})
ocell = counts$cell[counts$mark > 4]
for (i in round(seq(from=1, to=length(ocell), length.out=NP))){
    cdf = filter(subdf, cell == ocell[i])
    cmarks = unique(cdf$mark)
    par(xaxs="i")
    par(mar=rep(.15,4))
    ylim = c(0,max(cdf$val) + length(marks) - 1)
    xlim = c(1,max(cdf$pos))
    plot(1,1, xlim=xlim, ylim=ylim, yaxt='n', xaxt='n', col='white', ylab='', xlab='')
    for (j in 1:length(marks)){
        text(x=10, y=j-1, labels=m, cex=1, col=rgb(0,0,0,.5))
        if (j %% 2 == 1){col = 'darkblue'} else {col = 'royalblue'}
        if (m == 'DNase-seq'){col = 'forestgreen' }
        m = marks[j]
        if (m %in% cmarks){
            x = filter(cdf, mark == m)$val
            s = smooth.wavelet(x)
            mx = localMaxima(s)
            my = s[mx]
            lines(x + j - 1, type='s', yaxt='n', xaxt='n', lwd=2, col=col)
            # lines(s + j - 1, type='l', col='red', lwd=1)
            # points(mx, my + j - 1)
        } else {
            abline(h=j - 1, lwd=2, col=col)
        }
    }
}
dev.off()


layout(matrix(c(1,2), 1,2), T)
par(mar=c(.5,.5,2,.5))
image(momat, axes=F, col=c('black',viridis_pal(option='C')(100)))
mtext(paste0('Imputed', mark), side=3, line=.5)


par(mar=c(.5,.5,2,.5))
image(omat[rn,], axes=F, col=c('black',viridis_pal(option='C')(100)))
mtext(paste0('Observed', mark), side=3, line=.5)




# Cluster order:
cmat = rbind(imat, omat)
cmat[is.na(cmat)] = 0
dt = dist(t(cmat), 'euclidean')
ht = hclust(dt, 'ward.D')
cn = colnames(cmat)[ht$order]
rn = 1:100
rn = 1:nrow(imat)



png(paste0('~/EPIMAP_ANALYSIS/img/impobs_signals_region_',mark, '.png'), res=300, units='in', width=4, height=6)
layout(matrix(c(1,2), 1,2), T)
par(mar=c(.5,.5,2,.5))
image(imat[rn,], axes=F, col=c('black',viridis_pal(option='C')(100)))
mtext(paste0('Imputed', mark), side=3, line=.5)
par(mar=c(.5,.5,2,.5))
image(omat[rn,], axes=F, col=c('black',viridis_pal(option='C')(100)))
mtext(paste0('Observed', mark), side=3, line=.5)
dev.off()

png(paste0('~/EPIMAP_ANALYSIS/img/impobs_signals_region_reord_',mark, '.png'), res=300, units='in', width=4, height=6)
par(xaxs='i')
layout(matrix(c(1:4), 2,2), heights=c(1,5), widths=c(2,2), F)
par(mar=c(.25,.25,2,.25))
avgtrack = aggregate(val ~ pos, filter(df2, set == 'i'), mean)
plot(avgtrack, type='s', yaxt='n', xaxt='n')
mtext(paste('Imputed', mark), side=3, line=.5)
par(mar=rep(.25,4))
image(imat[rn,cn], axes=F, col=c('black',viridis_pal(option='C')(100)))
par(mar=c(.25,.25,2,.25))
avgtrack = aggregate(val ~ pos, filter(df2, set == 'o'), mean)
plot(avgtrack, type='s', yaxt='n', xaxt='n')
mtext(paste('Observed', mark), side=3, line=.5)
par(mar=rep(.25,4))
image(omat[rn,cn], axes=F, col=c('black',viridis_pal(option='C')(100)))
dev.off()


# Instead of avg track, explore with NMF:
k = 5
ocell = sort(unique(filter(df2, set == 'o')$cell))
onm = nnmf(omat[,ocell], k, method='scd', loss = 'mkl')
icell = sort(unique(filter(df2, set == 'i')$cell))
inm = nnmf(imat[,icell], k, method='scd', loss = 'mkl')

# matplot(onm$W, type = 'l', col = 1:5, lty='solid')

png(paste0('~/EPIMAP_ANALYSIS/img/impobs_signals_region_nmf_', k, '_',mark, '.png'), res=300, units='in', width=4, height=6)
par(xaxs='i')
layout(matrix(c(1:(2 * k + 4)), k + 2,2), heights=c(1, rep(.5,k),6), widths=c(3,3), F)
# Imputed:
par(mar=c(.15,.15,2,.15))
avgtrack = aggregate(val ~ pos, filter(df2, set == 'i'), mean)
plot(avgtrack, type='s', yaxt='n', xaxt='n', lwd=2)
mtext(paste('Imputed', mark), side=3, line=.5)
for (i in 1:k) {
    par(mar=rep(.15,4))
    plot(inm$W[,i], type='s', yaxt='n', xaxt='n', col='darkblue')
}
par(mar=rep(.15,4))
image(imat[rn,cn], axes=F, col=c('black',viridis_pal(option='C')(100)))
# Observed:
par(mar=c(.15,.15,2,.15))
avgtrack = aggregate(val ~ pos, filter(df2, set == 'o'), mean)
plot(avgtrack, type='s', yaxt='n', xaxt='n', lwd=2)
mtext(paste('Observed', mark), side=3, line=.5)
for (i in 1:k) {
    par(mar=rep(.15,4))
    plot(onm$W[,i], type='s', yaxt='n', xaxt='n', col='darkblue')
}
par(mar=rep(.15,4))
image(omat[rn,cn], axes=F, col=c('black',viridis_pal(option='C')(100)))
dev.off()


# Wavelet smoothing followed by pk calling:
smat = apply(omat[,ocell], 2, smooth.wavelet)

# Look at smoothed matrix.
layout(matrix(c(1:2), 2,1))
par(xaxs="i")
par(mar=c(.15,.15,.25,.15))
image(omat[,ocell], axes=F, col=c('black',viridis_pal(option='C')(100)))
par(mar=c(.15,.15,.25,.15))
image(smat, axes=F, col=c(viridis_pal(option='C')(100)))

x.rec = smooth.wavelet(avgtrack$val)
max.x = localMaxima(x.rec)
max.y = x.rec[max.x]
avg.track = apply(omat[,ocell], 1, mean)
avg.sm = apply(smat, 1, mean)
maxA.x = localMaxima(avg.sm)
maxA.y = avg.sm[maxA.x]

png(paste0(img, 'impobs_region_wavelet_pk_',mark, '.png'), res=300, units='in', width=10, height=2)
layout(matrix(c(1:2), 2,1))
par(xaxs="i")
par(mar=c(.15,.15,1.5,.15))
plot(avgtrack$val, type='s', yaxt='n', xaxt='n', lwd=2)
lines(x.rec, type='l', yaxt='n', xaxt='n', col='red', lwd=1)
points(max.x, max.y)
abline(v=max.x, lty='dotted', lwd=.5)
mtext(paste0(mark, ' - Wavelet Smoothing the Average Track'), line=.25)
par(mar=c(.15,.15,1.5,.15))
plot(avg.track, type='s', yaxt='n', xaxt='n', lwd=2)
lines(avg.sm, type='l', yaxt='n', xaxt='n', col='red', lwd=1)
points(maxA.x, maxA.y)
abline(v=maxA.x, lty='dotted', lwd=.5)
mtext(paste0(mark, ' - Averaging all Wavelet Smoothed Tracks'), line=.25)
dev.off()

# For each track, estimate pks:
NP =20
png(paste0(img, 'impobs_region_wavelet_pk_',mark, '.png'), res=300, units='in', width=10, height=NP / 3)
layout(matrix(c(1:NP), NP,1), T)
par(xaxs='i')
par(mar = rep(0,4))
for (i in round(seq(from=1, to=length(ocell), length.out=NP))){
    x = omat[,ocell[i]]
    s = smat[,i]
    mx = localMaxima(s)
    my = s[mx]
    plot(x, type='s', yaxt='n', xaxt='n', lwd=2)
    lines(s, type='l', axes=F, col='red', lwd=1)
    points(mx, my)
    abline(v=max.x, lwd=5, col=rgb(0,0,0,.15))
    # lines coinciding or close:
    dm = sapply(mx, function(x){min(abs(x  - max.x))})
    idx = which(dm < 2)
    abline(v=mx[idx], lwd=1, col='forestgreen', lty='dashed')
    abline(v=mx[-idx], lwd=1, col='red', lty='dashed')
}
dev.off()


# avg.track = apply(omat[,ocell], 1, mean)
# avg.sm = apply(smat, 1, mean)
# maxA.x = localMaxima(avg.sm)
# maxA.y = avg.sm[maxA.x]


# Dataframe of smoothed + cluster:
dfs = gather(data.frame(smat, pos=1:nrow(smat)), cell, val, -pos)
dt = dist(t(smat), 'euclidean')
ht = hclust(dt, 'ward.D')
scell = colnames(smat)
cn = scell[ht$order]
dfs$cell = factor(dfs$cell, levels=cn)
# NOTE: Lines inacurate representation. Needs to be bins.
# ggplot(filter(odf, pos < 200), aes(pos, cell, height=val, group=cell, fill=val)) +

ggplot(dfs, aes(pos, cell, height=val, group=cell, fill=val)) +
    # geom_vline(xintercept=max.x, lty='dotted', lwd=.5) + 
    geom_vline(xintercept=max.x, lwd=.5) + 
    geom_density_ridges_gradient(stat='identity', scale=5, color=NA) +
    scale_fill_viridis(name = "-log10p", option = "C") + 
    theme_minimal() + 
    theme(axis.line=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position='none',
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave(paste0(img, 'impobs_smoothedridges_pkcall_', mark, '.png'), dpi=250, width=6, height=8, units='in')


# Extract peaks:



# LATER: Integrate both the histone mark data and the DNase-seq




# ggplot(df2, aes(pos, cell, height=val, group=cell, fill=val)) +
#     facet_wrap(~ set, ncol=2) +
#     geom_density_ridges_gradient(stat='identity', scale=8, color=NA) +
#     scale_fill_viridis(name = "-log10p", option = "C") + 
#     theme_minimal() + 
#     theme(axis.line=element_blank(),
#           axis.text.y=element_blank(), axis.ticks.y=element_blank(),
#           axis.title.x=element_blank(), axis.title.y=element_blank(),
#           panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#           panel.grid.minor=element_blank(),plot.background=element_blank())

# Cluster order:
cells = sort(unique(filter(df2, set == 'o')$cell))
dt = dist(t(omat[,cells]), 'euclidean')
ht = hclust(dt, 'ward.D')
cn = cells[ht$order]
odf = filter(df2, set == 'o')
odf$cell = factor(odf$cell, levels=cn)

avgtrack = aggregate(val ~ pos, odf, mean)
# NOTE: Lines inacurate representation. Needs to be bins.
# Look at df
# ggplot(filter(odf, pos < 200), aes(pos, cell, height=val, group=cell, fill=val)) +
ggplot(odf, aes(pos, cell, height=val, group=cell, fill=val)) +
    geom_density_ridges_gradient(stat='identity', scale=5, color=NA) +
    scale_fill_viridis(name = "-log10p", option = "C") + 
    theme_minimal() + 
    theme(axis.line=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
ggsave(paste0(img, 'example_', mark, '_region_ridges.png'), dpi=250, width=4, height=6, units='in')

