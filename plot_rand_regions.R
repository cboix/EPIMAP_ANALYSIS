#!/usr/bin/R
# ------------------------------------
# Load all rand regions and plot them:
# ------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    source('~/data/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
} else {
    source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
}
source('~/EPIMAP_ANALYSIS/bin/auxiliary_chromImpute_functions.R')
library(dplyr)
library(viridis)
library(cba)
options(scipen=45)

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need evalfile and qcfile filenames.")
} else {        
    selfile = args[1]
    regprefix = args[2]
    imgprefix = args[3]
}

randpref = 'rand_regionfig/region_tracks_15344/'
selfile = paste0(randpref, '/track_selection_allmark_allfiles.tsv')
regprefix = paste0(randpref, '/regions/raw_regions')
imgprefix = paste0(randpref, '/rand_region_plot')

# Read in sample tables:
seldf = read.delim(selfile, sep="\t", header=F)
names(seldf) = c('id','mark')

# Temporary files (if rerun):
impfile = paste0(regprefix, '_allimp.tsv')
obsfile = paste0(regprefix, '_allobs.tsv')

if (file.exists(impfile) & file.exists(obsfile)){
    imat = as.matrix(read.delim(impfile, header=F, stringsAsFactors=F, sep="\t"))
    print(paste0("Read imat: ", paste0(dim(imat), collapse=", ")))
    omat = as.matrix(read.delim(obsfile, header=F, stringsAsFactors=F, sep="\t"))
    print(paste0("Read omat: ", paste0(dim(omat), collapse=", ")))
} else {
    print("Loading regions")
    pb = txtProgressBar(min = 0, max = nrow(seldf), style = 3)
    for (i in 1:nrow(seldf)){
        setTxtProgressBar(pb, i)
        # TRY OR PRINT ERROR.
        imp = as.numeric(scan(paste0(regprefix, '_imp_', i, '.tsv'), 'c', quiet=T))
        obs = as.numeric(scan(paste0(regprefix, '_obs_', i, '.tsv'), 'c', quiet=T))
        if (i == 1){
            imat = matrix(0, nrow=nrow(seldf), ncol=length(imp))
            omat = matrix(0, nrow=nrow(seldf), ncol=length(imp))
        }
        imat[i,] = imp
        omat[i,] = obs
    }
    close(pb)
    # Temporarily write tables:
    write.table(imat, impfile, row.names=F, quote=F, sep="\t", col.names=F)
    write.table(omat, obsfile, row.names=F, quote=F, sep="\t", col.names=F)
}

palette = viridis_pal()(100)
palette = cividis(100)

png(paste0(imgprefix, '_both_raw.png'), res=250, units='in', width=7,height=8)
par(yaxs='i')
par(xaxs='i')
layout(matrix(c(1:2),2,1), TRUE)
par(mar=rep(0.5,4))
image(t(omat), axes=F, col=palette, zlim=c(0,15))
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
par(mar=rep(0.5,4))
image(t(imat), axes=F, col=palette, zlim=c(0,15))
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
dev.off()

# Load in the per mark binarization cutoffs
ctdf = read.delim('ChromImpute/binarization_cutoffs_match_avg_distr.tsv', header=T)
ctdf = rbind(ctdf, data.frame(mark=c('ATAC-seq','CTCF','SMC3','RAD21','EP300','POLR2A'), cutoff=2))

# Hierarchically cluster the N random regions in each mark for observed only: 
method ='cosine'
ordfile = paste0(regprefix, method, '_markord.rda')
marks = sort(unique(as.character(seldf$mark)))
if (file.exists(ordfile)){
    load(ordfile)
} else {
    lr = list()
    lc = list()
    for (mark in marks){
        print(mark)
        ind = which(seldf$mark == mark)
        submat = omat[ind,]
        LIM = 4
        submat[submat > LIM] = LIM
        # Get order of columns and of rows:
        dt <- dist(t(submat), method)
        ht <- hclust(dt)
        lc[[mark]] <- order.optimal(dt, ht$merge)$order
        dt <- dist(submat, method)
        ht <- hclust(dt)
        lr[[mark]] <- order.optimal(dt, ht$merge)$order
    }
    # Save ordering of the marks:
    save(lr, lc, file=ordfile)
}

# Make figure:
height = 11
width = 5
ln = unlist(lapply(lr, length))
ln = c(sapply(ln, times=2, rep))
heights = ln / sum(ln) * height
sp=0.05

# png(paste0(imgprefix, '_both_cls.png'), res=450, units='in', width=width,height=height)
pdf(paste0(imgprefix, '_both_cls.pdf'), width=width,height=height)
par(yaxs='i')
par(xaxs='i')
layout(matrix(c(1:length(ln)),length(ln),1), heights=heights, TRUE)
for (mark in marks){
    # Get submatrix:
    ind = which(seldf$mark == mark)
    subomat = omat[ind,]
    subimat = imat[ind,]
    oLIM = 4
    iLIM = 2 * ctdf$cutoff[ctdf$mark == mark]
    subomat[subomat > oLIM] = oLIM
    subimat[subimat > iLIM] = iLIM
    ro = lr[[mark]]
    co = lc[[mark]]
    ur = TRUE
    # Make images
    par(mar=c(sp/2,sp,sp,sp))
    image(t(subomat[ro, co]), axes=F, col=palette, zlim=c(0,oLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), paste0(mark, " (N = ", nrow(subimat),")"), 
         cex=.65, col=rgb(1,1,1,.5))
    par(mar=c(sp,sp,0,sp))
    image(t(subimat[ro,co]), axes=F, col=palette, zlim=c(0,iLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), '(Imputed)', 
         cex=.45, col=rgb(1,1,1,.5))
}
dev.off()


png(paste0(imgprefix, '_both_cls_thresh.png'), res=450, units='in', width=width,height=height)
# pdf(paste0(imgprefix, '_both_cls_thresh.pdf'), width=width,height=height)
par(yaxs='i')
par(xaxs='i')
layout(matrix(c(1:length(ln)),length(ln),1), heights=heights, TRUE)
for (mark in marks){
    # Get submatrix:
    ind = which(seldf$mark == mark)
    oLIM = 2
    iLIM = ctdf$cutoff[ctdf$mark == mark]
    subomat = 1 * (omat[ind,] > oLIM)
    subimat = 1 * (imat[ind,] > iLIM)
    ro = lr[[mark]]
    co = lc[[mark]]
    # Make images
    par(mar=c(sp/2,sp,sp,sp))
    image(t(subomat[ro, co]), axes=F, col=c('darkblue','goldenrod1'))
    abline(v=par()$usr[1:2], lwd=.5)
    abline(h=par()$usr[3], lwd=.5)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), paste0(mark, " (N = ", nrow(subimat),")"), 
         cex=.85, col=rgb(1,1,1,.5))
    par(mar=c(sp,sp,0,sp))
    image(t(subimat[ro,co]), axes=F, col=c('darkblue', 'goldenrod1'))
    box(lwd=.25)
}
dev.off()

# Image (clustered by all??)

# Make figure:
height = 6
width = 3
t12marks = marks[!(marks %in% c('CTCF','SMC3','RAD21','EP300','POLR2A'))]
ln = sapply(lr, length)[t12marks]
ln = c(sapply(ln, times=2, rep))
heights = ln / sum(ln) * height
sp=0.05

png(paste0(imgprefix, '_both_cls_t12.png'), res=450, units='in', width=width,height=height)
# pdf(paste0(imgprefix, '_both_cls_t12.pdf'), width=width,height=height)
par(yaxs='i')
par(xaxs='i')
layout(matrix(c(1:length(ln)),length(ln),1), heights=heights, TRUE)
for (mark in t12marks){
    # Get submatrix:
    ind = which(seldf$mark == mark)
    subomat = omat[ind,]
    subimat = imat[ind,]
    oLIM = 4
    iLIM = 2 * ctdf$cutoff[ctdf$mark == mark]
    subomat[subomat > oLIM] = oLIM
    subimat[subimat > iLIM] = iLIM
    ro = lr[[mark]]
    co = lc[[mark]]
    ur = TRUE
    # Make images
    par(mar=c(sp/2,sp,sp,sp))
    image(t(subomat[ro, co]), axes=F, col=palette, zlim=c(0,oLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), paste0(mark, " (N = ", nrow(subimat),")"), 
         cex=.65, col=rgb(1,1,1,.5))
    par(mar=c(sp,sp,0,sp))
    image(t(subimat[ro,co]), axes=F, col=palette, zlim=c(0,iLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), '(Imputed)', 
         cex=.45, col=rgb(1,1,1,.5))
}
dev.off()

# png(paste0(imgprefix, '_both_cls_t12_flip.png'), res=450, units='in', width=height, height=2)
pdf(paste0(imgprefix, '_both_cls_t12_flip.pdf'), width=height,height=2)
par(yaxs='i')
par(xaxs='i')
layout(matrix(c(1:length(ln)),2, length(ln)), widths=heights, TRUE)
for (mark in t12marks){
    # Get submatrix:
    ind = which(seldf$mark == mark)
    subomat = omat[ind,]
    subimat = imat[ind,]
    oLIM = 4
    iLIM = 2 * ctdf$cutoff[ctdf$mark == mark]
    subomat[subomat > oLIM] = oLIM
    subimat[subimat > iLIM] = iLIM
    ro = lr[[mark]]
    co = lc[[mark]]
    ur = TRUE
    # Make images
    par(mar=c(sp,sp,sp,sp/2))
    image(subomat[ro, co], axes=F, col=palette, zlim=c(0,oLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), paste0(mark, " (N = ", nrow(subimat),")"), 
         srt=90, cex=.65, col=rgb(1,1,1,.5))
    par(mar=c(sp,0,sp,sp))
    image(subimat[ro,co], axes=F, col=palette, zlim=c(0,iLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), '(Imputed)', 
         srt=90, cex=.45, col=rgb(1,1,1,.5))
}
dev.off()


height = 3
pm = sort(c('ATAC-seq','DNase-seq','H2AFZ','H3K27ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac'))
ln = sapply(lr, length)[pm]
heights = ln / sum(ln) * height
# JUST THE Punctate marks
# "peak analysis"

png(paste0(imgprefix, '_both_cls_pm_sidebyside.png'), res=450, units='in', width=4,height=height)
# pdf(paste0(imgprefix, '_both_cls_pm_sidebyside.pdf'), width=4,height=height)
par(yaxs='i')
par(xaxs='i')
sp = 0.025
layout(matrix(c(1:(2*length(ln))),length(ln),2, byrow=TRUE), heights=heights, TRUE)
for (mark in pm){
    # Get submatrix:
    ind = which(seldf$mark == mark)
    subomat = omat[ind,]
    subimat = imat[ind,]
    oLIM = 4
    iLIM = 2 * ctdf$cutoff[ctdf$mark == mark]
    subomat[subomat > oLIM] = oLIM
    subimat[subimat > iLIM] = iLIM
    ro = lr[[mark]]
    co = lc[[mark]]
    ur = TRUE
    # Make images
    par(mar=c(sp,sp,sp,sp))
    image(t(subomat[ro, co]), axes=F, col=palette, zlim=c(0,oLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), paste0(mark, " (N = ", nrow(subimat),")"), 
         cex=.65, col=rgb(1,1,1,.5))
    par(mar=c(sp,sp,sp,sp))
    image(t(subimat[ro,co]), axes=F, col=palette, zlim=c(0,iLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), 'Imputed', 
         cex=.5, col=rgb(1,1,1,.5))
}
dev.off()


# png(paste0(imgprefix, '_both_cls_pm_sidebyside_diag.png'), res=450, units='in', width=4,height=height)
pdf(paste0(imgprefix, '_both_cls_pm_sidebyside_diag.pdf'), width=4,height=height)
par(yaxs='i')
par(xaxs='i')
sp = 0.025
layout(matrix(c(1:(2*length(ln))),length(ln),2, byrow=TRUE), heights=heights, TRUE)
for (mark in pm){
    # Get submatrix:
    ind = which(seldf$mark == mark)
    ro = lr[[mark]]
    co = lc[[mark]]
    subomat = omat[ind,][ro, co]
    subimat = imat[ind,][ro, co]
    tco = names(sort(apply(subomat,2,sum), decreasing=TRUE))[1:2000]
    # dco = diag.mat2(subomat[,tco]/5, cutoff=.5, ratio=.25)[[2]]
    dco = diag.mat2(subimat[,tco]/5, cutoff=.5, ratio=.25)[[2]]
    subomat = subomat[, dco]
    subimat = subimat[, dco]
    oLIM = 4
    iLIM = 2 * ctdf$cutoff[ctdf$mark == mark]
    subomat[subomat > oLIM] = oLIM
    subimat[subimat > iLIM] = iLIM
    ur = TRUE
    # Make images
    par(mar=c(sp,sp,sp,sp))
    image(t(subomat), axes=F, col=palette, zlim=c(0,oLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), paste0(mark, " (N = ", nrow(subimat),")"), 
         cex=.65, col=rgb(1,1,1,.5))
    par(mar=c(sp,sp,sp,sp))
    image(t(subimat), axes=F, col=palette, zlim=c(0,iLIM), useRaster=ur)
    text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), 'Imputed', 
         cex=.5, col=rgb(1,1,1,.5))
}
dev.off()


