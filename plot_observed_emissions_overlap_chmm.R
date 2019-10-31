#!/usr/bin/R
# ------------------------------------------------------------
# Plot the comparison of observed + states for all epigenomes:
# ------------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))


imgpref = paste0(img, "model_diagnostics/")

# Load model parameters:
model = 'observed_aux_18_on_mixed_impobs_QCUT'
modelpref = 'observed_aux_18'
dcol <- read.table(paste0('CHMM_', modelpref, '_colors.tsv'),header=F)
names(dcol) <- c('state','name','color')
dcol$col <- sapply(dcol$color,function(x){
                       s <- unlist(strsplit(as.character(x),split=','))
                       s <- as.numeric(s)
                       rgb(s[1],s[2],s[3],max=255) })
# Get emissions:
emfile = 'emissions_18_core_K27ac.txt'
emdf = read.delim(emfile)
emat = as.matrix(emdf[,-1])
rownames(emat) = emdf[,1]
emdf$state = 1:18
emdf = emdf[,-1]
emlong = gather(emdf, mark, pct, -state)

# Get mnemonics:
mfile = 'CHMM_observed_aux_18_mnemonic.tsv'
mndf = read.delim(mfile, header=F)
names(mndf) = c('state','name','string','colstr','color')
dcol = merge(dcol, mndf)
dcol = dcol[order(dcol$state),]
rownames(dcol) = NULL

# Orders:
pltord = rev(1:18)
markord = colnames(emat)
epiord = cellorder

# Load observed/emission comparison file:
ovfile = 'ChromHMM/calls/all_observed_emissions_overlap.tsv'
ovdf = read.delim(ovfile, header=T, sep="\t")
ovwide = spread(ovdf, sample, percent)

# Order state first, mark second:
ovwide$mark = factor(ovwide$mark, levels=markord)
ovwide = ovwide[order(ovwide$mark, decreasing=T),]
ovwide = ovwide[order(ovwide$state, decreasing=T),]
ovmat = as.matrix(ovwide[,cellorder])

# Order emissions too:
emlong$mark = factor(emlong$mark, levels=markord)
emlong = emlong[order(emlong$mark, decreasing=T),]
emlong = emlong[order(emlong$state, decreasing=T),]

# Counts overall:
freqdir = 'ChromHMM/freqs/'
modeldir = paste0(freqdir, model,'/')

tabdf = read.delim(paste0(modeldir, 'state_counts_092219.tsv'), sep="\t", header=T)
mat = as.matrix(tabdf[,-1])
ord = order(mat[nrow(mat),], decreasing=T)
mat = mat[,ord]
frac = sweep(mat,2,apply(mat,2,sum), '/')
bssids = colnames(frac)
bssids = bssids[bssids %in% cellorder]
frac = frac[,bssids]
dfracs = round(apply(frac, 1, mean) * 100, 2)

pdf(paste0(imgpref, model, '_observed_vs_emissions.pdf'), width=8, height=5.25)
# States:
sp=0.1
bsp=2.5
layout(matrix(1:(19 * 3), 19,3, byrow=TRUE), heights=rep(.2, 19), widths=c(1, .5 ,4.25,2), TRUE)
par(mar=c(sp,sp,sp,sp))
plot(1,1, type='n', axes=F)
text(1,1, 'State')
par(mar=c(sp,sp,sp,sp))
plot(1,1, type='n', axes=F)
text(1,1, 'Emissions')
par(mar=c(sp,sp,sp,sp))
plot(1,1, type='n', axes=F)
text(1,1, '% Mark Occurrence per Epigenome')
for (state in 1:18){
    sub.emlong= emlong[emlong$state == state,]
    ind = which(ovwide$state == state)
    par(mar=c(sp,sp,sp,sp))
    image(t(as.matrix(state)), col=dcol$col[state], axes=F, useRaster=TRUE)
    text(-.95, 0, paste0(state, '  ', as.character(dcol$name)[state], '  (', dfracs[state], '%)'), adj=0)
    box(lwd=.5)
    # Real emissions:
    par(mar=c(sp,bsp,sp,sp))
    image(t(as.matrix(sub.emlong$pct)), col=col1, axes=F, useRaster=TRUE, zlim=c(0,1))
    yat = seq(0, 1, length.out=6)
    text(par()$usr[1] - .05 * diff(par()$usr[1:2]), yat, 
         as.character(sub.emlong$mark), adj=1, xpd=TRUE, cex=.35)
    box(lwd=.5)
    par(mar=c(sp,sp,sp,sp))
    image(t(ovmat[ind,]), col=col1, axes=F, useRaster=TRUE, zlim=c(0,1))
    box(lwd=.5)
}
dev.off()

# Count what fraction of epigenomes look off, in Biv states and EnhG1:
ind = c(which(ovwide$state==14 & ovwide$mark == 'H3K27ac'),
        which(ovwide$state==7 & ovwide$mark == 'H3K27ac'),
        which(ovwide$state==15 & ovwide$mark == 'H3K27ac'))
mean(apply(ovwide[ind,cellorder] > .5, 1, sum)  / 833 )







