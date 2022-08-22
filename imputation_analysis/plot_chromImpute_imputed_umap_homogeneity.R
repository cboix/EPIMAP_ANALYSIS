#!/usr/bin/R
# ---------------------------------------
# Plot the imputed + observed UMAP
# Plot similar MDS
# ---------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'auxiliary_function_general_repel.R'))
library(ggplot2)
library(viridis)
library(ggrepel)
library(dendextend)
library(uwot)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(ggpubr)
library(igraph)
library(PRROC)
options(scipen=45) # So we dont get issues writing integers into bedfiles

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    fixedstate = FALSE # Look at the gwcorr matrices.
    nregions = 0 
    dataset = 'spearman'
    print(paste("No arguments supplied.", 
                "Defaulting to loading nregions =", 
                nregions, 'and dataset', dataset))
} else {        
    fixedstate = as.logical(args[1])
    nregions = as.numeric(args[2])
    dataset = args[3]
}

# Load specific distance matrices:
if (fixedstate){
    commandArgs <- function(trailingOnly=TRUE){ c(nregions, dataset) }
    source(paste0(bindir, 'load_region_distance_matrices.R'))
    setprefix = paste('region',nregions, dataset,'distances_', sep="_")
} else { 
    source(paste0(bindir, 'load_distance_matrices.R'))
    setprefix = 'distances_'
    metric = 'genome-wide correlation'
}

today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "imp_distance/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, setprefix)

parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}

# ------------------------------------------------------------
# Calculate simple clustering-coefficient/homogeneity metrics:
# ------------------------------------------------------------
# 1. Clustering coefficient (using the N nearest neighbors)
# 2. Within-group distance vs. Out-group distance?
# Find samples with similar properties:
llist = sapply(cellorder, function(x){
                   cind = which((meta[cellorder,'infoline'] == meta[x,'infoline']) * 
                                (meta[cellorder,'lifestage'] == meta[x,'lifestage']) *
                                (meta[cellorder,'GROUP'] == meta[x,'GROUP']) == 1)
                   cellorder[cind] } )
llist = llist[lapply(llist, length) > 2]
llist = unique(llist)

rlist = lapply(llist, function(x){ sample(cellorder, length(x))})
meta$GROUP = factor(meta$GROUP, levels=odf$GROUP)

# Plot ordered by group:
png(paste0(imgpref, 'matrices_grpord_allmarks.png'),res=300,units='in',width=13,height=2.2)
layout(matrix(1:(13 * 3), nrow=3, byrow=F), heights=c(.2, 1,1), widths=rep(1,13), TRUE)
for (i in 1:(length(marks) -1)){
    # Load both observed and imputed datasets:
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells:
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    omat = omat[cn, cn]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    imat = imat[icn, icn]
    # Also calculate the same exact restricted matrix:
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    bomat = omat[bcn, bcn]
    bimat = imat[bcn, bcn]
    # bicn = rownames(reord(bimat))
    bocn = rownames(reord(bomat))
    # Plot, ordered by group:
    grpord = bocn[order(meta[bocn,'GROUP'])]
    # Plot 2x2 ordering with the colors on the side...
    sp=0.1
    par(mar=c(sp, sp,1,sp))
    meta.image(metamat[grpord,5,drop=F], colvals=colvals, cex=0, horiz=F, useRaster=TRUE)
    mtext(paste(mark), cex=.75)
    par(mar=c(sp, sp,sp,sp))
    image(bomat[grpord, grpord], useRaster=T, axes=F, col=colspec)
    image(bimat[grpord, grpord], useRaster=T, axes=F, col=colspec)
}
dev.off()


# Init:
ccdf = c()
all.wdf = c()
all.samedf = c()
all.randdf = c()
all.wodf = c()
# For dist:
all.stdf = c()
all.bdf = c()

for (i in 1:(length(marks) -1)){
    # Load both observed and imputed datasets:
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells:
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    omat = omat[cn, cn]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    imat = imat[icn, icn]
    # Also calculate the same exact restricted matrix:
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    bomat = omat[bcn, bcn]
    bimat = imat[bcn, bcn]
    bocn = rownames(reord(bomat))
    bicn = rownames(reord(bimat))

    # Within/between:
    bodf = data.frame(bomat)
    bodf$id = rownames(bodf)
    bodf = gather(bodf, id2, dist, -id)
    bodf = bodf[!is.na(bodf$dist),]
    bodf$group1 = meta[bodf$id, 'GROUP']
    bodf$group2 = meta[bodf$id2, 'GROUP']
    bodf$type = 'Observed'
    bidf = data.frame(bimat)
    bidf$id = rownames(bidf)
    bidf = gather(bidf, id2, dist, -id)
    bidf = bidf[!is.na(bidf$dist),]
    bidf$group1 = meta[bidf$id, 'GROUP']
    bidf$group2 = meta[bidf$id2, 'GROUP']
    bidf$type = 'Imputed'
    bdf = rbind(bodf, bidf)

    # Aggregate all the between/within:
    bdf$corr = 1 - bdf$dist
    bdf$set = 'Between'
    bdf$set[bdf$group1 == bdf$group2] = 'Within'
    stdf = aggregate(corr ~ group1 + type + set, bdf, mean)
    stdf = spread(stdf, set, corr)
    stdf = stdf[!is.na(stdf$Within),]
    if (nrow(stdf) > 0){
        stdf$ratio = stdf$Within / stdf$Between
        stdf$mark = mark
        all.stdf = rbind(all.stdf, stdf)
    }
    bdf$mark = mark
    all.bdf = rbind(all.bdf, bdf)


    # Calculate homogeneity statistics:
    # On same set, calc avg. distance between members of the same set
    tform = make.tform(metamat[bcn,'group'], norm=TRUE, u=odf$GROUP)
    # TODO: Mean distance and variance in distances?
    bomat.norm = bomat - mean(bomat)
    bimat.norm = bimat - mean(bimat)
    # gomat = t(tform) %*% (bomat.norm %*% tform)
    # gimat = t(tform) %*% (bimat.norm %*% tform)
    gomat = t(tform) %*% (bomat %*% tform)
    gimat = t(tform) %*% (bimat %*% tform)

    # Only within dist --> (TODO) show fit:
    wdf = data.frame(obs=diag(gomat), imp=diag(gimat), mark=mark)
    wdf$GROUP = rownames(wdf)
    wdf =merge(wdf, odf)
    wdf = rbind(all.wdf, wdf)
    fit = lm(imp ~ obs, wdf)
    print(summary(fit))

    png(paste0(imgpref, 'avg_group_dist_impobs_', mark,'.png'),res=300,units='in',width=8,height=5)
    layout(matrix(1:2, ncol=2))
    par(mar=c(4,4,2,1))
    plot(wdf$obs,wdf$imp, pch=19, xlab='Avg. within dist (Observed)', ylab='Avg. within dist (Imputed)',
         col=odf$COLOR, xlim=c(0,1), ylim=c(0,1), bty='n')
    mtext(paste(mark, '(Only within-group distances)'))
    # Multiple comparisons of between/within dist:
    par(mar=c(4,4,2,1))
    plot(gomat, gimat, pch=19, xlab='Avg. dist (Observed)', ylab='Avg. dist (Imputed)', 
         xlim=c(0,1), ylim=c(0,1), bty='n')
    abline(0,1)
    mtext(paste(mark, '(All cross-group distances)'))
    dev.off()

    # All comparisons:
    bodf = data.frame(gomat)
    bodf$GROUP = rownames(bodf)
    bodf = gather(bodf, group2, dist, -GROUP)
    bodf = merge(bodf, odf)
    bodf$type = 'Observed'
    # 
    bidf = data.frame(gimat)
    bidf$GROUP = rownames(bidf)
    bidf = gather(bidf, group2, dist, -GROUP)
    bidf = merge(bidf, odf)
    bidf$type = 'Imputed'
    bdf = rbind(bodf, bidf)
    bdf = bdf[!is.na(bdf$dist),]

    gplot = ggplot(bdf, aes(group2, dist, color=group2)) + 
        geom_point() + 
        facet_grid(type ~ GROUP, scales='free_y') + 
        scale_color_manual(values=colvals$group) + 
        labs('Average cross-cluster distance') +
        theme_pubclean()  + 
        theme(axis.text.x  = element_blank(),
              axis.ticks.x  = element_blank(),
              panel.background = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              # axis.line = element_line(colour = "black"),
              panel.border = element_rect(colour = "black", fill=NA, size=1))
    ggsave(paste0(imgpref, 'avgdist_crossgroup', mark, '.png'), gplot, dpi=450, units='in', width=15, height=5)


    # Calc some within-group vs out-group metrics:
    wometrics = sapply(odf$GROUP, function(x){
                           samp = bcn[metamat[bcn, 'group'] == x]
                           nosamp = bcn[metamat[bcn, 'group'] != x]
                           # Eval in-group/out-group in both types:
                           somat = bomat[samp, samp]
                           simat = bimat[samp, samp]
                           with.obs = mean(somat[lower.tri(somat)])
                           with.imp = mean(simat[lower.tri(simat)])
                           wv.obs = var(somat[lower.tri(somat)])
                           wv.imp = var(simat[lower.tri(simat)])
                           nomat = bomat[samp, nosamp]
                           nimat = bimat[samp, nosamp]
                           out.obs = mean(nomat)
                           out.imp = mean(nimat)
                           ov.obs = var(c(nomat))
                           ov.imp = var(c(nimat))
                           return(c(with.obs / out.obs, with.imp / out.imp,
                                    wv.obs / ov.obs, wv.imp / ov.imp))
              })


    png(paste0(imgpref, 'with_out_mean_impobs_', mark,'.png'),res=300,units='in',width=5.5,height=5)
    par(mar=c(4,4,1.5,1))
    rn = max(wometrics[1:2,], na.rm=T) * 1.1
    plot(wometrics[1,], wometrics[2,], pch=19, col=odf$COLOR, 
         xlab='Avg. dist in-group / out-group (Observed)', 
         ylab='Avg. dist in-group / out-group (Imputed)', 
         xlim=c(0,rn), ylim=c(0,rn), bty='n')
    abline(0,1)
    mtext(paste(mark,'(in/out-group mean distance)'))
    dev.off()

    png(paste0(imgpref, 'with_out_var_impobs_', mark,'.png'),res=300,units='in',width=5.5,height=5)
    rn = max(wometrics[3:4,], na.rm=T) * 1.1
    par(mar=c(4,4,1.5,1))
    plot(wometrics[3,], wometrics[4,], pch=19, col=odf$COLOR, 
         xlab='Avg. var in-group / out-group (Observed)', 
         ylab='Avg. var in-group / out-group (Imputed)', 
         xlim=c(0,rn), ylim=c(0,rn), bty='n')
    abline(0,1)
    axis(1)
    axis(2)
    mtext(paste(mark,'(in/out-group variance in distance)'))
    dev.off()

    wodf = data.frame(t(wometrics), mark=mark)
    names(wodf)[1:4] = c('wm','om','wv','ov')
    wodf$GROUP = rownames(wodf)
    wodf = merge(wodf, odf)
    all.wodf = rbind(all.wodf, wodf)

}



# Calculate homogeneity statistics:
# Init:
# NTOP = 50
ccdf = c()
all.samedf = c()
all.randdf = c()
for (i in 1:(length(marks) -1)){
    # Load both observed and imputed datasets:
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells:
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    omat = omat[cn, cn]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    imat = imat[icn, icn]
    # Also calculate the same exact restricted matrix:
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    bomat = omat[bcn, bcn]
    bimat = imat[bcn, bcn]
    bocn = rownames(reord(bomat))
    bicn = rownames(reord(bimat))

    NTOP = round(length(bocn)/ 10)
    # NTOP = 25

    # Make adjacency:
    # oadj = apply(bonorm, 1, function(x){
    oadj = apply(bomat, 1, function(x){
                     thresh = sort(x)[NTOP]
                     1 * (x <= thresh) })
    # iadj = apply(binorm, 1, function(x){
    iadj = apply(bimat, 1, function(x){
                     thresh = sort(x)[NTOP]
                     1 * (x <= thresh) })
    og = graph_from_adjacency_matrix(oadj)
    ig = graph_from_adjacency_matrix(iadj)
    occ = transitivity(og)
    icc = transitivity(ig)
    ccdf = rbind(ccdf, data.frame(obs=occ, imp=icc, mark=mark))

    # locc = transitivity(og, type='local')
    # licc = transitivity(ig, type='local')
    # plot(locc, licc, pch=19, col=meta[vertex_attr(og)$name,'COLOR'])

    # Look at similar samples:
    slist = sapply(1:length(llist), function(i){
                       x = llist[[i]]
                       samp = bcn[bcn %in% x]
                       nosamp = bcn[!(bcn %in% samp)]
                       if (length(samp) > 1){
                           somat = bomat[samp, samp]
                           simat = bimat[samp, samp]
                           with.obs = mean(somat[lower.tri(somat)])
                           with.imp = mean(simat[lower.tri(simat)])
                           nomat = bomat[samp, nosamp]
                           nimat = bimat[samp, nosamp]
                           out.obs = mean(nomat)
                           out.imp = mean(nimat)
                           return(c(i,meta[x[1], c('infoline')],
                                    with.obs / out.obs, with.imp / out.imp)) }})
    samedf = data.frame(do.call(rbind, slist))
    samedf$mark = mark
    all.samedf = rbind(all.samedf, samedf)

    rslist = sapply(1:length(rlist), function(i){
                       x = rlist[[i]]
                       samp = bcn[bcn %in% x]
                       nosamp = bcn[!(bcn %in% samp)]
                       if (length(samp) > 1){
                           somat = bomat[samp, samp]
                           simat = bimat[samp, samp]
                           with.obs = mean(somat[lower.tri(somat)])
                           with.imp = mean(simat[lower.tri(simat)])
                           nomat = bomat[samp, nosamp]
                           nimat = bimat[samp, nosamp]
                           out.obs = mean(nomat)
                           out.imp = mean(nimat)
                           return(c(i,meta[x[1], c('infoline')],
                                    with.obs / out.obs, with.imp / out.imp)) }})
    randdf = data.frame(do.call(rbind, rslist))
    randdf$mark = mark
    all.randdf = rbind(all.randdf, randdf)
}


# -----------------------------
# Clustering coefficient plots:
# -----------------------------
gplot = ggplot(ccdf, aes(obs, imp, label=mark)) + 
    geom_point() + 
    geom_text_repel() + 
    xlim(0,1) + ylim(0,1) + 
    geom_abline(intercept=0,slope=1) + 
    labs(x='Observed', y='Imputed', title='Clustering Coefficient (10% closest neighbors)') + 
    theme_pubr()
ggsave(paste0(imgpref, 'transitivity_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)

ccdf = ccdf[order(ccdf$imp/ccdf$obs),]
ccdf$mark = factor(ccdf$mark, levels=ccdf$mark)

gplot = ggplot(ccdf, aes(mark, imp/obs)) + 
    geom_bar(stat='identity', fill='grey75') + 
    geom_hline(yintercept=1) + 
    labs(x='Mark/Assay', y='Imputed/Observed Coefficient', title='Ratio of Clust. Coeff. (10% closest neighbors)') + 
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1)) 
ggsave(paste0(imgpref, 'transitivity_ratio_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)


# Plot the within mark/without etc:
gplot = ggplot(all.wodf, aes(wm, om, color=GROUP)) + 
    geom_point() + 
    geom_abline() + 
    scale_color_manual(values=colvals$group) +
    labs(x='Avg. dist in-group / out-group (Observed)', y='Avg. dist in-group / out-group (Imputed)') +
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'inoutgroup_mean.png'), gplot,dpi=300,units='in',width=5,height=5)


# Plot the within mark/without etc:
gplot = ggplot(all.wodf, aes(wv, ov, color=GROUP)) + 
    geom_point() + 
    geom_abline() + 
    scale_color_manual(values=colvals$group) +
    labs(x='Var. dist. in-group / out-group (Observed)', y='Var. dist. in-group / out-group (Imputed)') +
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'inoutgroup_var.png'), gplot,dpi=300,units='in',width=5,height=5)


names(all.samedf) = c('num','infoline','obs','imp','mark')
all.samedf = merge(all.samedf, unique(meta[,c('GROUP','infoline')])) 
all.samedf = merge(all.samedf, odf) 
all.samedf$obs = as.numeric(as.character(all.samedf$obs))
all.samedf$imp = as.numeric(as.character(all.samedf$imp))

names(all.randdf) = c('num','infoline','obs','imp','mark')
all.randdf = merge(all.randdf, unique(meta[,c('GROUP','infoline')])) 
all.randdf = merge(all.randdf, odf) 
all.randdf$obs = as.numeric(as.character(all.randdf$obs))
all.randdf$imp = as.numeric(as.character(all.randdf$imp))


# Plot close biological samples (TODO: compare to rand samples or samples within same group)
mx = max(c(all.samedf$obs, all.samedf$imp))
gplot = ggplot(all.samedf, aes(obs, imp, color=GROUP)) + 
    geom_point() + 
    geom_abline() + 
    xlim(0, mx) + ylim(0,mx)  +
    scale_color_manual(values=colvals$group) +
    labs(x='Avg. dist in/out-group (Observed)', y='Avg. dist in/out-group (Imputed)') +
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'inout_sametype_mean.png'), gplot,dpi=300,units='in',width=5,height=5)

# Comparison if we do random samples of same size:
mx = max(c(all.randdf$obs, all.randdf$imp))
gplot = ggplot(all.randdf, aes(obs, imp)) + 
    geom_point() + 
    geom_abline() + 
    xlim(0, mx) + ylim(0,mx)  +
    scale_color_manual(values=colvals$group) +
    labs(x='Avg. dist in/out-group (Observed)', y='Avg. dist in/out-group (Imputed)') +
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'inout_randtype_mean.png'), gplot,dpi=300,units='in',width=5,height=5)


# -------------------------------------
# Plots for between/within correlation:
# -------------------------------------
# swdf = spread(all.stdf[,c('type','ratio','mark','group1')], type, ratio)
# rlim = range(all.stdf$ratio)
all.stdf$l2fc = log2(all.stdf$ratio)
swdf = spread(all.stdf[,c('type','l2fc','mark','group1')], type, l2fc)
rlim = range(all.stdf$l2fc)
cor(swdf$Imputed, swdf$Observed, method='pearson')
cor(swdf$Imputed, swdf$Observed, method='spearman')

gplot = ggplot(swdf, aes(Observed, Imputed, color=group1)) + 
    geom_point() + 
    xlim(rlim[1],rlim[2]) + ylim(rlim[1],rlim[2]) + 
    scale_color_manual(values=colvals$group) +
    geom_hline(yintercept=0, lty='dotted') + 
    geom_vline(xintercept=0, lty='dotted') + 
    geom_abline(intercept=0,slope=1) + 
    geom_smooth(method='lm', color='blue') + 
    stat_cor(label.y.npc = .87, color='darkgrey') +
    stat_regline_equation(label.y.npc = .95, color='darkgrey') + 
    labs(x='Observed (log2 ratio)', y='Imputed (log2 ratio)', title='Within / Between dataset distance\nfor each sample group in each mark') + 
    theme_pubr() + theme(legend.position = 'none')
ggsave(paste0(imgpref, 'withinbetween_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)

head(swdf[order(swdf$Observed, decreasing=T),],10)
head(swdf[order(swdf$Observed, decreasing=F),],10)


# Plot connecting:
stlong = gather(all.stdf[,c('group1','type','mark','Between','Within')],  set, value, -group1, -type, -mark)
stlong$full = paste0(stlong$type, '.', stlong$set)
stwide = spread(stlong[,c('group1','full','value','mark')], full, value)

gplot = ggplot(stwide, aes(x=Observed.Within, y=Observed.Between, xend=Imputed.Within, yend=Imputed.Between, color=group1)) + 
    geom_segment(lineend = "round", linejoin = "round",
                 size = .25, arrow = arrow(length = unit(0.1, "inches"))) +
    theme_pubr() + theme(legend.position = 'none') + 
    scale_color_manual(values=colvals$group) +
    labs(x='Within-group average dataset distance', 
         y='Between-group average dataset distance') + 
    geom_abline(intercept=0,slope=1)
ggsave(paste0(imgpref, 'withinbetween_arrows_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)


# gplot = ggplot(stwide, aes(x=1-Observed.Within, y=1-Observed.Between, xend=1-Imputed.Within, yend=1-Imputed.Between, color=group1)) + 
gplot = ggplot(stwide, aes(x=Observed.Within, y=Observed.Between, xend=Imputed.Within, yend=Imputed.Between, color=group1)) + 
    geom_segment(lineend = "round", linejoin = "round", size = .25, arrow = arrow(length = unit(0.1, "inches"))) +
    theme_pubr() + theme(legend.position = 'none') + 
    scale_color_manual(values=colvals$group) +
    labs(x='Within-group average dataset distance', y='Between-group average dataset distance') + 
    scale_x_log10() + scale_y_log10() + 
    geom_abline(intercept=0,slope=1)
ggsave(paste0(imgpref, 'withinbetween_arrows_logscale_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)


# gplot = ggplot(all.stdf, aes(x=Within, y=Between, pch=type, color=group1)) + 
gplot = ggplot(all.stdf, aes(x=1-Within, y=1-Between, color=type)) + 
    geom_point() + 
    theme_pubr() + theme(legend.position=c(0.15, .9)) + 
    scale_color_manual(values=c('indianred','royalblue'), name='') +
    guides(color = guide_legend(override.aes=list(size=4))) + 
    labs(x='Within-group average dataset distance', 
         y='Between-group average dataset distance') + 
    geom_abline(intercept=0,slope=1)
ggsave(paste0(imgpref, 'withinbetween_pts_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)



# gplot = ggplot(all.stdf, aes(x=Within, y=Between, pch=type, color=group1)) + 
gplot = ggplot(all.stdf, aes(x=1-Within, y=1-Between, color=type)) + 
    geom_point() + 
    theme_pubr() + theme(legend.position=c(0.15, .9)) + 
    scale_color_manual(values=c('indianred','royalblue'), name='') +
    guides(color = guide_legend(override.aes=list(size=4))) + 
    labs(x='Within-group average dataset distance', 
         y='Between-group average dataset distance') + 
    scale_x_log10() + scale_y_log10()  + 
    geom_abline(intercept=0,slope=1)
ggsave(paste0(imgpref, 'withinbetween_pts_logscale_impobs_allmarks.png'), gplot,dpi=300,units='in',width=5,height=5)



# Same-group classification task:
# -------------------------------
# Exclude other/cancer:
sub.bdf = all.bdf[!(all.bdf$group1 %in% c('Cancer','Other')),]
sub.bdf = sub.bdf[!(sub.bdf$group2 %in% c('Cancer','Other')),]

aucdf = c()
for (mark in unique(all.bdf$mark)){
    print(mark)
    sodf = sub.bdf[sub.bdf$mark == mark & sub.bdf$type == 'Observed',]
    sidf = sub.bdf[sub.bdf$mark == mark & sub.bdf$type == 'Imputed',]
    # Calculate the AUROCs:
    oroc = roc.curve(scores.class0=sodf$corr[sodf$set == 'Within'],
                     scores.class1=sodf$corr[sodf$set == 'Between'], curve=TRUE)
    iroc = roc.curve(scores.class0=sidf$corr[sidf$set == 'Within'],
                     scores.class1=sidf$corr[sidf$set == 'Between'], curve=TRUE)
    print(oroc$auc)
    print(iroc$auc)
    df = data.frame(mark=mark, auc=c(oroc$auc, iroc$auc), type=c('Observed','Imputed'))
    aucdf = rbind(aucdf, df)
}

gplot = ggplot(aucdf, aes(mark, auc, fill=type)) + 
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('indianred','royalblue'), name='') + 
    labs(x='',y='AUC') + 
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1)) 
ggsave(paste0(imgpref, 'aucclassify_impobs_allmarks_samesamples.png'), gplot,dpi=300,units='in',width=5,height=5)


# -------------------------------------
# Repeat the AUC analysis for ALL data:
# -------------------------------------
# For dist:
all.pairsdf = c()
for (i in 1:(length(marks) -1)){
    # Load both observed and imputed datasets:
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Within/between:
    bodf = data.frame(omat)
    bodf$id = rownames(bodf)
    bodf = gather(bodf, id2, dist, -id)
    bodf = bodf[!is.na(bodf$dist),]
    bodf$group1 = meta[bodf$id, 'GROUP']
    bodf$group2 = meta[bodf$id2, 'GROUP']
    bodf$type = 'Observed'
    bidf = data.frame(imat)
    bidf$id = rownames(bidf)
    bidf = gather(bidf, id2, dist, -id)
    bidf = bidf[!is.na(bidf$dist),]
    bidf$group1 = meta[bidf$id, 'GROUP']
    bidf$group2 = meta[bidf$id2, 'GROUP']
    bidf$type = 'Imputed'
    bdf = rbind(bodf, bidf)
    # Aggregate all the between/within:
    bdf$corr = 1 - bdf$dist
    bdf$set = 'Between'
    bdf$set[bdf$group1 == bdf$group2] = 'Within'
    bdf$mark = mark
    all.pairsdf = rbind(all.pairsdf, bdf)
}

# Exclude other/cancer:
sub.bdf = all.pairsdf[!(all.pairsdf$group1 %in% c('Cancer','Other')),]
sub.bdf = sub.bdf[!(sub.bdf$group2 %in% c('Cancer','Other')),]

aucdf = c()
for (mark in unique(all.bdf$mark)){
    print(mark)
    sodf = sub.bdf[sub.bdf$mark == mark & sub.bdf$type == 'Observed',]
    sidf = sub.bdf[sub.bdf$mark == mark & sub.bdf$type == 'Imputed',]
    # Calculate the AUROCs:
    oroc = roc.curve(scores.class0=sodf$corr[sodf$set == 'Within'],
                     scores.class1=sodf$corr[sodf$set == 'Between'], curve=TRUE)
    iroc = roc.curve(scores.class0=sidf$corr[sidf$set == 'Within'],
                     scores.class1=sidf$corr[sidf$set == 'Between'], curve=TRUE)
    print(oroc$auc)
    print(iroc$auc)
    df = data.frame(mark=mark, auc=c(oroc$auc, iroc$auc), type=c('Observed','Imputed'), set='All Samples')
    aucdf = rbind(aucdf, df)
}

gplot = ggplot(aucdf, aes(mark, auc, fill=type)) + 
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('indianred','royalblue'), name='') + 
    labs(x='',y='AUC') + 
    theme_pubr() + theme(axis.text.x=element_text(angle=90, hjust=1)) 
ggsave(paste0(imgpref, 'aucclassify_impobs_allmarks_allsamples.png'), gplot,dpi=300,units='in',width=5,height=5)








# ---------------------------------------
# Add the read length and other metadata:
# ---------------------------------------
nam = 'all_submitted_released'
outfile <- paste0('Annotation/',nam,'_metadata.tsv')
metadata <- read.delim(outfile,sep="\t",header=T)
metadata$UID <- paste(metadata$Experiment.accession,
                      metadata$Biological.replicate.s.,
                      metadata$Technical.replicate,
                      sep=".")

fqmeta = metadata[metadata$File.format == 'fastq',]
# Year
fqmeta = unique(fqmeta[,c('Experiment.accession', 'Experiment.date.released', 'Lab',
                          'Run.type', 'Read.length', 'Platform', 'File.Status', 'Audit.NOT_COMPLIANT')])
names(fqmeta)[1] = 'Accession'
fqmeta = merge(accmap, fqmeta)
fqmeta$Year = sub("-.*","", fqmeta$Experiment.date.released)

rldf = aggregate(Read.length ~ Accession + id + Epitope, fqmeta, max)

# ------------------------------------------------------------------------------------
# Show that modeling read length, lab, and other covariates 
# is barely accounting for variability when also accounting for lifestage, group, etc.
# ------------------------------------------------------------------------------------
# Prune fqmeta into a unique metadata table:
fqmeta = fqmeta[fqmeta$File.Status %in% c('released','submitted'),]
fqdf = aggregate(cbind(Read.length, Year) ~ Accession + Lab + Epitope + id + Run.type, fqmeta, max)
yrdf = aggregate(Year ~ Accession + Epitope + id, fqmeta, max)
lbdf = aggregate(Lab ~ Accession + Epitope + id, fqmeta, function(x){paste(sort(unique(x)), collapse=', ')})
rtdf = aggregate(Run.type ~ Accession + Epitope + id, fqmeta, function(x){paste(sort(unique(x)), collapse=', ')})
rownames(rtdf) = rtdf$Accession
rownames(rldf) = rldf$Accession
rownames(yrdf) = yrdf$Accession
rownames(lbdf) = lbdf$Accession

fdf = fqmeta[,c('Lab','Read.length','Year','Run.type')]

# ----------------------------
# Plot descriptive statistics:
# ----------------------------
g1 = ggplot(fdf, aes(Read.length)) + 
    labs(x='Read Length (bp)', y='Number of datasets') + 
    geom_histogram() + theme_pubr()
g2 = ggplot(fdf, aes(Run.type)) + 
    labs(x='Run Type', y='Number of datasets') + 
    geom_bar() + theme_pubr() +
    theme(axis.text.x=element_text(angle=90, hjust=1)) 
g3 = ggplot(fdf, aes(Lab)) + 
    labs(x='Lab', y='Number of datasets') + 
    geom_bar() + theme_pubr() +
    theme(axis.text.x=element_text(angle=90, hjust=1)) 
ggsave(paste0(imgpref, 'batches_stats_hist_rlen.png'), g1, dpi=300,units='in',width=3.5,height=2.5)
ggsave(paste0(imgpref, 'batches_stats_hist_runtype.png'), g2, dpi=300,units='in',width=1.5,height=3)
ggsave(paste0(imgpref, 'batches_stats_hist_lab.png'), g3, dpi=300,units='in',width=3.5,height=5)


# # TODO: Move this all to a new script!
# fmat = cbind(rl=fqmeta$Read.length > 36, 
#              yr=fqmeta$Year > 2013,
#              pr=fqmeta$Project == 'Roadmap',
#              rt=fqmeta$Run.type == 'single-ended')
# fmat = cbind(rl=mdf$rl> 36, 
#              yr=mdf$yr > 2013,
#              pr=mdf$rl == 'Roadmap',
#              # rt=mdf$rt == 'single-ended',
#              tp=mdf$type == 'tissue')
# dist(t(fmat))

# Attributes for datasets:
mdf = accmap[,c('id','Epitope','Accession')]
mdf$group = meta[mdf$id,'GROUP']
mdf$ls = meta[mdf$id,'lifestage']
mdf$type = meta[mdf$id,'type']
mdf$proj = meta[mdf$id,'Project']
mdf$rl = rldf[mdf$Accession,'Read.length']
mdf$yr = yrdf[mdf$Accession,'Year']
mdf$lb = lbdf[mdf$Accession,'Lab']
mdf$rt = rtdf[mdf$Accession,'Run.type']

rmdf = aggregate(rl ~ group, mdf, mean)
rmdf = rmdf[order(rmdf$rl),]
mdf$group = factor(mdf$group, levels=rmdf$group)

ggplot(mdf, aes(group, rl, color=group)) + 
    geom_boxplot() + 
    geom_jitter() + 
    theme_pubr() + 
    scale_color_manual(values=colvals$group) + 
    theme(axis.text.x = element_text(angle=90, hjust=1), legend.position='none') 

mdf$rl.cut = cut(mdf$rl,breaks = c(0,36,50,75,90,101), right=FALSE)

# TODO: Enrichments by lab - all 

# -------------------------------
# Gather all distances as a table
# -------------------------------
aodf = c()
aidf = c()
for (i in 1:(length(marks) -1)){
    # Load both observed and imputed datasets:
    mark = marks[i]
    print(mark)
    omat = obsll[[mark]]
    imat = ll[[mark]]
    # Restrict to chosen cells:
    cn = rownames(omat)[rownames(omat) %in% cellorder]
    omat = omat[cn, cn]
    icn = rownames(imat)[rownames(imat) %in% cellorder]
    imat = imat[icn, icn]
    # Also calculate the same exact restricted matrix:
    bcn = icn[icn %in% cn]
    names(bcn) = NULL
    bomat = omat[bcn, bcn]
    bimat = imat[bcn, bcn]
    bodf = data.frame(bomat * upper.tri(bomat))
    bodf$id = rownames(bodf)
    bodf = gather(bodf, id2, dist, -id)
    bodf$Epitope = mark
    bidf = data.frame(bimat * upper.tri(bimat))
    bidf$id = rownames(bidf)
    bidf = gather(bidf, id2, dist, -id)
    bidf$Epitope = mark
    # Add: 
    aodf = rbind(aodf, bodf)
    aidf = rbind(aidf, bidf)
}
aodf = aodf[aodf$id != aodf$id2,]
aidf = aidf[aidf$id != aidf$id2,]
aodf = aodf[aidf$dist != 0,] # Remove double counting.
aidf = aidf[aidf$dist != 0,]
aodf$type = 'Observed'
aidf$type = 'Imputed'

# Need to model distance by is same of each of x attributes.
accmap1 = accmap[,c('id','Epitope','Accession')]
accmap2 = accmap1
names(accmap2) = c('id2', 'Epitope','Accession2')
aodf = merge(merge(aodf, accmap1), accmap2)
aidf = merge(merge(aidf, accmap1), accmap2)

# Add observations:
aodf$same.group = meta[aodf$id,'GROUP'] == meta[aodf$id2, 'GROUP']
aodf$same.ls = meta[aodf$id,'lifestage'] == meta[aodf$id2, 'lifestage']
aodf$same.type = meta[aodf$id,'type'] == meta[aodf$id2, 'type']
aodf$same.proj = meta[aodf$id,'Project'] == meta[aodf$id2, 'Project']
aodf$same.rl = rldf[aodf$Accession,'Read.length'] == rldf[aodf$Accession2, 'Read.length']
aodf$same.yr = yrdf[aodf$Accession,'Year'] == yrdf[aodf$Accession2, 'Year']
aodf$same.lb = lbdf[aodf$Accession,'Lab'] == lbdf[aodf$Accession2, 'Lab']
aodf$same.rt = rtdf[aodf$Accession,'Run.type'] == rtdf[aodf$Accession2, 'Run.type']

# Add observations:
aidf$same.group = meta[aidf$id,'GROUP'] == meta[aidf$id2, 'GROUP']
aidf$same.ls = meta[aidf$id,'lifestage'] == meta[aidf$id2, 'lifestage']
aidf$same.type = meta[aidf$id,'type'] == meta[aidf$id2, 'type']
aidf$same.proj = meta[aidf$id,'Project'] == meta[aidf$id2, 'Project']
aidf$same.rl = rldf[aidf$Accession,'Read.length'] == rldf[aidf$Accession2, 'Read.length']
aidf$same.yr = yrdf[aidf$Accession,'Year'] == yrdf[aidf$Accession2, 'Year']
aidf$same.lb = lbdf[aidf$Accession,'Lab'] == lbdf[aidf$Accession2, 'Lab']
aidf$same.rt = rtdf[aidf$Accession,'Run.type'] == rtdf[aidf$Accession2, 'Run.type']


# ---------------------------------------------------------------------
# Z-score the per-mark correlations (per mark, in observed and imputed)
# ---------------------------------------------------------------------
bdf = rbind(aodf, aidf)
bdf$corr = 1 - bdf$dist
mbdf = aggregate(corr ~ Epitope + type, bdf, mean)
sbdf = aggregate(corr ~ Epitope + type, bdf, sd)
names(mbdf)[3] = 'mean'
names(sbdf)[3] = 'sd'
bdf = merge(bdf, merge(mbdf, sbdf))
bdf$zscore = (bdf$corr - bdf$mean) / bdf$sd

# Regr:
i1 = lm(zscore ~ 0 + type*(same.group + same.ls + same.type + same.proj + same.rl + same.yr + same.lb), bdf)
cfdf = data.frame(coefficients(summary(i1)))
colnames(cfdf) = c('Est', 'SE','tvalue', 'pvalue')
cfdf$term = rownames(cfdf)
cfdf = cfdf[grep('same', rownames(cfdf)),]

i1 = lm(zscore ~ 1 + type*(same.group + same.ls + same.type + same.proj + same.rl + same.yr + same.lb), bdf)
cfdf = data.frame(coefficients(summary(i1)))

ggplot(bdf, aes(type, zscore, fill=same.rl)) + 
    # scale_fill_manual(values=c("#00AFBB", "#E7B800")) + 
    scale_fill_manual(values=colpair[c(1,10)]) + 
    geom_boxplot() + theme_pubr()

ggplot(bdf, aes(type, zscore, fill=same.rt)) + 
    scale_fill_manual(values=colpair[c(1,10)]) + 
    geom_boxplot() + theme_pubr()

ggplot(bdf, aes(type, zscore, fill=same.rl)) + 
    scale_fill_manual(values=colpair[c(1,10)]) + 
    stat_compare_means() + 
    geom_violin() + 
    geom_boxplot(width=.25) +
    theme_pubr()

ggplot(bdf, aes(type, 1 - dist, fill=same.rl)) + 
    facet_wrap(~Epitope) + 
    scale_fill_manual(values=colpair[c(1,10)]) + 
    geom_boxplot() + theme_pubr()


# ----------------------------------------
# Look at it as ratio diff in correlation:
# ----------------------------------------
# bwdf = spread(bdf[,c('Epitope','id','id2','Accession','Accession2', 'type', 'dist', 'same.group', 
#                      'same.ls', 'same.type', 'same.proj', 'same.rl','same.yr','same.lb')], type, dist)
# bwdf$ratio = (1 - bwdf$Imputed) / (1 - bwdf$Observed)
# i1 = lm(ratio ~ same.group + same.ls + same.type + same.proj + same.rl + same.yr + same.lb, bwdf)


# Compare by z-score:
bwdf = spread(bdf[,c('Epitope','id','id2','Accession','Accession2', 'type', 'zscore', 'same.group', 
                     'same.ls', 'same.type', 'same.proj', 'same.rl','same.yr','same.lb', 'same.rt')], type, zscore)
bwdf$diff = bwdf$Imputed - bwdf$Observed

i1 = lm(diff ~ same.group + same.ls + same.type + same.proj + 
        same.rl + same.yr + same.lb + same.rt, bwdf)
cfdf = data.frame(coefficients(summary(i1)))
colnames(cfdf) = c('Estimate', 'SE','tvalue', 'pvalue')
cfdf$term = rownames(cfdf)
cfdf = cfdf[grep('same', rownames(cfdf)),]

cov.map = c('Year','Sample Type','Read Length','Project','Lifestage','Lab','Bio. group', 'Read Type')
names(cov.map) = c('yr','type','rl','proj','ls','lb','group', 'rt')
rownames(cfdf) = sub("TRUE","", rownames(cfdf))
rownames(cfdf) = sub("same.","", rownames(cfdf))
cfdf$Attribute = cov.map[rownames(cfdf)]

library(huxtable)
hp = as_hux(cfdf[,c(6,1:4)], add_colnames=TRUE)
position(hp) = "left"
top_border(hp)[1,]=2
bottom_border(hp)[1,]=1
bottom_border(hp)[nrow(hp),]=2
width(hp) = .75
quick_pdf(hp, file="coeff_comparison_batches.pdf")

cfdf$ymax = cfdf$Estimate + cfdf$SE * 2
cfdf$ymin = cfdf$Estimate - cfdf$SE * 2

# Basic plot (unfinished)
bt = barplot(cfdf$Estimate, col=ifelse(cfdf$Estimate < 0, colrb[90],colrb[10]), 
             border=NA, ylim=range(c(cfdf$ymax, cfdf$ymin)))
arrows(bt, cfdf$ymax, bt,cfdf$ymin, angle=90, code=3, length=0.1)


# Separate calculations:
i1 = lm(zscore ~ same.group + same.ls + same.type + same.proj + 
        same.rl + same.yr + same.lb + same.rt, bdf[bdf$type == 'Imputed',])
o1 = lm(zscore ~ same.group + same.ls + same.type + same.proj + 
        same.rl + same.yr + same.lb + same.rt, bdf[bdf$type == 'Observed',])

cfdf = rbind(data.frame(coefficients(summary(i1)), type='Imputed'), 
             data.frame(coefficients(summary(o1)), type='Observed'))
colnames(cfdf) = c('Estimate', 'SE','tvalue', 'pvalue', 'type')
cfdf$term = rownames(cfdf)
cfdf = cfdf[grep('same', rownames(cfdf)),]
cfdf$term = sub("TRUE.*$","", cfdf$term)
cfdf$term = sub("same.","", cfdf$term)
cfdf$Attribute = cov.map[cfdf$term]
cfdf$ymax = cfdf$Estimate + cfdf$SE * 2
cfdf$ymin = cfdf$Estimate - cfdf$SE * 2
cfdf$term = factor(cfdf$term, levels=c('group','ls','type','lb','proj','rl','rt','yr'))
cfdf$type = factor(cfdf$type, levels=rev(c('Observed','Imputed')))
cfdf = cfdf[order(cfdf$type),]
cfdf = cfdf[order(cfdf$term),]

xlim = range(c(cfdf$ymax, cfdf$ymin))
xlim[1] = floor(10 * xlim[1]) / 10
xlim[2] = ceiling(10 * xlim[2]) / 10

png(paste0(imgpref, 'coeff_batcheffects_indpt_reg.png'), res=300,units='in',width=6,height=3)
par(mar=c(2.5,7,.25,.75))
bt = barplot(Estimate ~ type + term, cfdf, beside=T, col=rev(c('royalblue','indianred')),
             space=rep(c(.5,0),nrow(cfdf)/2),
             border=NA, xlim=xlim, horiz=T,ylab=NA, xlab=NA, axisnames=F,xaxt='n')
yat = apply(bt, 2, mean)
text(x=parpos(1,.0), y=yat, cov.map[levels(cfdf$term)], xpd=TRUE, adj=1)
xat = seq(xlim[1], xlim[2], by=.1)
axis(1, at=xat, tck=-.02, labels=F)
text(x=xat, y=parpos(2,0.06), xat, xpd=TRUE)
mtext('Estimated effect size on z-scored sample similarity', side=1, line=1.3)
abline(v=0)
arrows(cfdf$ymax, bt, cfdf$ymin, bt, angle=90, code=3, length=0.035)
legend('topright',c('Observed','Imputed'), pch=15, pt.cex=2, cex=1, col=c('royalblue','indianred'), bty='n')
ypad = 0.75
segments(x0=parpos(1,0.235), y0=c(yat[1], yat[4]) - ypad,
         x1=parpos(1,0.235), y1=c(yat[3], yat[8]) + ypad, lwd=1, xpd=TRUE)
text(x=parpos(1,0.27), y=c(mean(yat[c(1,3)]), mean(yat[c(4,8)])),
     labels=c('Biological', 'Batch/Technical'), srt=90,adj=.5, xpd=TRUE)
dev.off()


# Separate calculations:
bio.i1 = lm(zscore ~ same.group + same.ls + same.type, bdf[bdf$type == 'Imputed',])
bio.o1 = lm(zscore ~ same.group + same.ls + same.type, bdf[bdf$type == 'Observed',])
summary(bio.i1)$r.sq
summary(bio.o1)$r.sq

# Separate calculations:
tech.i1 = lm(zscore ~ same.proj + same.rl + same.yr + same.lb + same.rt, bdf[bdf$type == 'Imputed',])
tech.o1 = lm(zscore ~ same.proj + same.rl + same.yr + same.lb + same.rt, bdf[bdf$type == 'Observed',])
summary(tech.i1)$r.sq
summary(tech.o1)$r.sq











# # Lab does increase correlation:
# i0 = lm(diff ~ same.group + same.ls + same.type + same.proj + same.rl + same.yr, bwdf)
# anova(i0, i1)



