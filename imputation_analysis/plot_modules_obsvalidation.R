#!/usr/bin/R
# ---------------------------------------
# Plot the modules validation
# Make correlation based FDR predictions:
# ------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
today <- format(Sys.time(), "%m%d%y")
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

library(rhdf5)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(Matrix)
library(scales)
library(plotrix)
library(PRROC)
library(huxtable)
library(DescTools) # For Gini

# Arguments:
filepref = 'cls_merge2_wH3K27ac100_300'

tagline = 'ChromHMM Enhancers'
imgdir = paste0(img, "clusters/") 
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need clusters filename")
} else {        
    filepref = args[1]
    tagline = args[2]
    if (length(args) > 2){ 
        imgdir = args[3] 
    } else { 
        imgdir = paste0(img, "clusters/") 
    }
}


# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){
    c(filepref, tagline, imgdir) }
source(paste0(bindir, 'load_modules_data.R'))

# Prefix:
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'clusters_')


#' Calculate position relative to par()$usr 
#'
#' @param axis 1 is x; 2 is y;
#' @param shift percent of axis to shift over
#' @return position shifted from start of x or y axis
#' @export
parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}

# ---------------------
# Load expression data:
# ---------------------
lddir = 'linking_data/'
# Data files:
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K4me1_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
fnames = list.files(path=lddir, pattern='*_precomputed_corr.hdf5')

# Load actual ENH + H3K27ac matrix:
fullmatfile = 'Enhancer_H3K27ac_matrix_062619.mtx.gz'
enhmatfile = 'Enhancer_H3K27ac_matrix_enhonly_062619.mtx.gz'
enhmargfile = 'Enhancer_H3K27ac_margins_enhonly_062619.tsv.gz'
mat = read.delim(gzfile(enhmatfile), sep="\t", header=T)
matmarg = read.delim(gzfile(enhmargfile), sep="\t", header=T)
# Enhancer mapping:
matnames = scan('Enhancer_matrix_names.txt', "c")
enhind = as.numeric(scan('Enhancer_indices.txt', 'c')) + 1
enhmap = rep(0, max(enhind))
enhmap[enhind] = 1:length(enhind)

# Module mapping:
enhsetfile = 'modules_enhancer_sets.Rda'
load(enhsetfile)

# Read in which are imputed, observed:
iodf = read.delim('ChromImpute/mixobs_table_fordist.tsv', header=F)
colnames(iodf) = c('id','mark','prefix')
iodf$imputed = FALSE
iodf$imputed[grep('impute', iodf[,3])] = TRUE

# Take 235 samples to recapitulate this:
siodf = iodf[iodf$mark == 'H3K27ac' & iodf$imputed == FALSE,]


# --------------------------------------------------------------
# Create a clusters centers matrix using only the observed data:
# --------------------------------------------------------------
keepcol = which(matnames %in% siodf$id)
smat = mat[mat$col %in% keepcol,]

fracdf = c()
for (i in 1:300){
    print(i)
    rlist = enhsets[[i]]
    nr = length(rlist)
    rmat = smat[smat$row %in% rlist,]
    # Aggregte:
    rdf = aggregate(row ~ col, rmat, length)
    rdf$frac = rdf$row / nr
    rdf$set = paste0('c', i-1)
    fracdf = rbind(fracdf, rdf)
    print(max(rdf$frac))
}
fwide = spread(fracdf[,c('col','set','frac')], set, frac, fill=0)
ocenters = as.matrix(fwide[,-1])
rownames(ocenters) = matnames[fwide$col]
ocenters = ocenters[, colnames(centers)]



# ------------------------------------------------------
# See specificity of module centers (for OBSERVED ONLY):
# ------------------------------------------------------
# TODO: CHANGE to ocenters
Z = 1 * (centers > 0.25)
marg = apply(Z, 2, sum)
oZ = 1 * (ocenters > 0.25)
omarg = apply(oZ, 2, sum)

# Numbers for manuscript
ind = which(marg >= 150)
rind = which(marg < 150)
mean(marg[rind]) / 833 * 100
mean(marg[ind]) / 833 * 100

odf = data.frame(allmarg=marg, obsmarg=omarg, type='Specific')
odf$type[marg >=150] = 'Broad'

odf$allfrac = odf$allmarg / 833
odf$obsfrac = odf$obsmarg / 235

odf$allgini = apply(Z, 2, Gini)
odf$obsgini = apply(oZ, 2, Gini)

g1 = ggplot(odf, aes(allmarg, obsmarg, color=type)) +
    geom_vline(xintercept=150, col='indianred', lty='dashed') + 
    geom_smooth(method='lm', color='black') + 
    geom_point(alpha=0.75) + theme_pubr() +
    scale_y_log10() + 
    scale_x_log10() +
    scale_color_manual(values=c('indianred','royalblue'), name='') + 
    theme(legend.position = c(0.15, 0.95)) + 
    labs(x='# Active samples in module (of 833 total)',
         y='# Active samples in module (of 235 observed)')
ggsave(paste0(imgpref, "observed_validation_margins.png"), g1, width=5,height=5, units='in', dpi=300)
ggsave(paste0(imgpref, "observed_validation_margins.pdf"), g1, width=5,height=5)


g2 = ggplot(odf, aes(allgini, obsgini, color=type)) +
    geom_smooth(method='lm', color='black') + 
    geom_point(alpha=0.75) + theme_pubr() +
    scale_color_manual(values=c('indianred','royalblue'), name='') + 
    theme(legend.position ='none') + 
    labs(x='Gini Coefficient (across 833 samples)',
         y='Gini coefficient (across 235 observed)')
ggsave(paste0(imgpref, "observed_validation_gini.png"), g2, width=5,height=5, units='in', dpi=300)
ggsave(paste0(imgpref, "observed_validation_gini.pdf"), g2, width=5,height=5)

garr = ggarrange(g1,g2, ncol=2, nrow=1)
ggsave(paste0(imgpref, "observed_validation_both.png"), garr, width=8,height=4, units='in', dpi=300)
ggsave(paste0(imgpref, "observed_validation_both.pdf"), garr, width=8,height=4)

bdf = odf[odf$allmarg >=150,]
mean(bdf$obsfrac)

odf[order(odf$obsgini),]

fit = lm(obsmarg ~ allmarg, odf)
summary(fit)


fit = lm(obsgini~ allgini, odf)
summary(fit)


sum(odf$obsmarg == 0)
sum(odf$obsmarg == 1)
sum(odf$obsmarg > 1)

aggregate(obsfrac ~ type, odf, mean)

odf$obstype = 'os'
odf$obstype[odf$obsfrac > .18] = 'ob'

znam = rownames(odf)[odf$obsmarg == 0]

sum(counts[[1]][znam])


# --------------------------------------------
# Make complete modules dataset for other use:
# --------------------------------------------
# ENH locations:
ddir = 'DHS_Index_WM201902/'
dpref = 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19'
dmlrdafile = paste0(ddir, dpref, '_enh.core.srt.Rda')
load(dmlrdafile)

# Data:
module.types = list(broad=rownames(odf)[odf$allmarg >= 150],
                    specific=rownames(odf)[odf$allmarg < 150])
module.assign = enhsets
names(module.assign) = paste0('c', 1:300 -1)
module.centers = centers

save(module.types, module.centers,
     module.assign, enhdf,
     file='modules_types_and_data.Rda')


# -------------------------------------------
# TODO: Also validate with raw observed data?
# -------------------------------------------
# sample = 'BSS00439'
mnames = sort(matnames)
mcol = which(mnames %in% matnames[keepcol])
# sampind = which(mnames == sample)
# meta[mnames[sampind],]

h5ls(H3K27ac.file)

# --------------------------------
# Get mapping between old and new:
# --------------------------------
ndf = read.delim('DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv', header=T)
colnames(ndf) = c('name','chr','cls.new')
new.names = ndf$name
old.names = enhdf$name

# Merge + map over:
bothdf = merge(ndf[,c('name','cls.new')], enhdf[,c('name','cls')])
newmap = rep(0, max(bothdf$cls))
newmap[bothdf$cls] = bothdf$cls.new
new.enhsets = lapply(enhsets, function(x){newmap[x]})

# Write the fix onto the cluster assignments:
clsdf = read.delim('cls_300_assign_withindexingerror.bed')
names(clsdf) = c('cls','id.old')
clsdf$id = newmap[clsdf$id.old]
clsdf = clsdf[order(clsdf$id),c('cls','id')]
write.table(clsdf,'cls_300_assign_new_sorted_coords.bed', quote=F, row.names=F, sep="\t")

# -----------------------------------------------------------
# Get average strengths for each sample + module combination:
# -----------------------------------------------------------
NTOT = 833
mmat = matrix(0, nrow=NTOT, ncol=length(new.enhsets))
for (i in 1:(NTOT %/% 10 + 1)){
    print(i)
    ind = ((i - 1) * 10 + 1):min(i * 10, NTOT)
    h5f = H5Fopen(H3K27ac.file)
    h5d = h5f&"matrix"
    hmat = h5d[ind,]
    H5Dclose(h5d)
    H5Fclose(h5f)
    mmat[ind,] = sapply(new.enhsets, function(x){ apply(hmat[,x],1, mean) })
}
rownames(mmat) = mnames
colnames(mmat) = paste0('c', 1:300 - 1)

# Get average strengths for each sample (H3K4me1)
emat = matrix(0, nrow=NTOT, ncol=length(new.enhsets))
for (i in 1:(NTOT %/% 10 + 1)){
    print(i)
    ind = ((i - 1) * 10 + 1):min(i * 10, NTOT)
    h5f = H5Fopen(H3K4me1.file)
    h5d = h5f&"matrix"
    hmat = h5d[ind,]
    H5Dclose(h5d)
    H5Fclose(h5f)
    emat[ind,] = sapply(new.enhsets, function(x){ apply(hmat[,x],1, mean) })
}
rownames(emat) = mnames
colnames(emat) = paste0('c', 1:300 - 1)

# Plot these against the coeff.
cmat = centers[mnames,]

plot(mmat, cmat)

plot(emat * mmat, cmat, xlim=c(0,50))

plot(emat, mmat)

c2 = t(diag.mat(t(cmat))[[1]])
c2 = diag.mat(c2)[[1]]

c2 = diag.mat(cmat)[[1]]
c2 = t(diag.mat(t(c2))[[1]])


# TODO: Assignments not matching --> Need to reconcile old + new sorted matrices!!!
png(paste0(imgpref, "mark_avg_all.png"), res=400, units='in', width=7, height=3)
layout(matrix(1:4, ncol=4))
sp=0.25
par(mar=c(sp,sp,1.5,sp))
image(t(c2), axes=F, useRaster=T)
box(lwd=.25)
mtext('Cluster Centers')
plt.mmat = mmat
plt.mmat[plt.mmat > 5] = 5
image(t(plt.mmat[rownames(c2), colnames(c2)]), axes=F, useRaster=T)
box(lwd=.25)
mtext('H3K27ac')
plt.emat = emat
plt.emat[plt.emat > 5] = 5
image(t(plt.emat[rownames(c2), colnames(c2)]), axes=F, useRaster=T)
box(lwd=.25)
mtext('H3K4me1')
image(t(plt.emat[rownames(c2), colnames(c2)] * plt.mmat[rownames(c2), colnames(c2)]), axes=F, useRaster=T)
box(lwd=.25)
mtext('H3K27ac x H3K4me1')
dev.off()


# Only observed:
rowobs = rownames(c2)[rownames(c2) %in% siodf$id]
o2 = c2[rowobs,]

# TODO: Assignments not matching --> Need to reconcile old + new sorted matrices!!!
png(paste0(imgpref, "mark_avg_obs.png"), res=400, units='in', width=7, height=2)
layout(matrix(1:4, ncol=4))
sp=0.25
par(mar=c(sp,sp,1.5,sp))
image(t(o2), axes=F, useRaster=T)
box(lwd=.25)
mtext('Cluster Centers')
plt.mmat = mmat
plt.mmat[plt.mmat > 5] = 5
image(t(plt.mmat[rownames(o2), colnames(o2)]), axes=F, useRaster=T)
box(lwd=.25)
mtext('H3K27ac')
plt.emat = emat
plt.emat[plt.emat > 5] = 5
image(t(plt.emat[rownames(o2), colnames(o2)]), axes=F, useRaster=T)
box(lwd=.25)
mtext('H3K4me1')
image(t(plt.emat[rownames(o2), colnames(o2)] * plt.mmat[rownames(o2), colnames(o2)]), axes=F, useRaster=T)
box(lwd=.25)
mtext('H3K27ac x H3K4me1')
dev.off()


# Add margins too:
odf$memarg = apply(emat > 2,2,sum)
odf$obs.memarg = apply(emat[rowobs,] > 2,2,sum)
odf$acmarg = apply(mmat > 2,2,sum)
odf$obs.acmarg = apply(mmat[rowobs,] > 2,2,sum)
odf$obs.bothmarg = apply((emat[rowobs,] > 2) * (mmat[rowobs,] * 2),2,sum)

odf$megini = apply(emat, 2, Gini)
odf$obs.megini = apply(emat[rowobs,], 2, Gini)
odf$acgini = apply(mmat, 2, Gini)
odf$obs.acgini = apply(mmat[rowobs,], 2, Gini)
odf$obs.bothgini = apply((emat * mmat) > 2,2,Gini)
odf$obs.bothgini2 = apply((emat[rowobs,] > 2) * (mmat[rowobs,] * 2),2,Gini)

g3 = ggplot(odf, aes(allmarg, obs.bothmarg, color=type)) +
    geom_vline(xintercept=150, col='indianred', lty='dashed') + 
    geom_smooth(method='lm', color='black') + 
    geom_point(alpha=0.75) + theme_pubr() +
    scale_y_log10() + 
    scale_x_log10() +
    scale_color_manual(values=c('indianred','royalblue'), name='') + 
    theme(legend.position = c(0.15, 0.95)) + 
    labs(x='# Active samples in module (of 833 total)',
         y='# samples with average signal > 2 in both H3K27ac and H3K4me1\n(of 235 observed)')
ggsave(paste0(imgpref, "observed_validation_mark_margins.png"), g3, width=5,height=5, units='in', dpi=300)

g4 = ggplot(odf, aes(allgini, obs.bothgini, color=type)) +
    geom_smooth(method='lm', color='black') + 
    geom_point(alpha=0.75) + theme_pubr() +
    scale_color_manual(values=c('indianred','royalblue'), name='') + 
    theme(legend.position ='none') + 
    labs(x='Gini Coefficient (across 833 samples)',
         y='Gini coefficient in H3K27ac and H3K4me1 (across 235 observed)')
ggsave(paste0(imgpref, "observed_validation_mark_gini.png"), g4, width=5,height=5, units='in', dpi=300)

garr = ggarrange(g3,g4, ncol=2, nrow=1)
ggsave(paste0(imgpref, "observed_validation_mark_both.png"), garr, width=8,height=4, units='in', dpi=300)
ggsave(paste0(imgpref, "observed_validation_mark_both.pdf"), garr, width=8,height=4)






