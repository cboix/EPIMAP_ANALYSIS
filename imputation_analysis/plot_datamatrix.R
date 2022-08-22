#!/usr/bin/R
# Plot the processing status report for ENCODE3/ROADMAP DATA
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
library(ggrepel)

# Locations:
imgdir = paste0(img, "metadata/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'meta_')

# -----------------------
# Plot the data matrices:
# -----------------------
png(paste0(imgpref,'matrix_alone.png'),res=300,units='in',width=10,height=2)
par(mar=c(.25,3,.25,.25))
image(t(avail[eporder,cellorder]), col=c('white','darkblue'), axes=F, useRaster=TRUE)
# mtext('Histone ChIP + DNase-seq + ATAC-seq', side=3)
text(y=seq(0,1, length.out=nrow(avail)),
     x=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=eporder, srt=0, adj=1, xpd=TRUE,cex=.3)
dev.off()

# Keep 833 cells
co833 = cellorder

# Order by # assays and by group
rownames(meta) = meta$id
cells = merge(cells, meta)
print(dim(cells))
# cells$group = as.character(meta[as.character(cells$id), 'GROUP'])
# cells$group = factor(cells$group, levels=rdcol$GROUP)
cells = cells[order(-cells$file, cells$GROUP), ]
cellorder = as.character(cells$id,1)

# Write/clean up data for release:
write.table(sort(cellorder), '~/EPIMAP_ANALYSIS/db/data_to_release/samplelist_859_original.tsv', quote=F, row.names=F, col.names=F, sep="\t")
write.table(sort(co833), '~/EPIMAP_ANALYSIS/db/data_to_release/samplelist_833_final.tsv', quote=F, row.names=F, col.names=F, sep="\t")
mt859 = meta[sort(cellorder),]
write.table(mt859, '~/EPIMAP_ANALYSIS/db/data_to_release/metadata_main.tsv', quote=F, row.names=F, sep="\t")


# Make labels:
labels = meta[cellorder, 'GROUP']
faclabels = as.matrix(as.numeric(labels))
colset = as.character(rdcol$COLOR)
lablist = label.runs(faclabels, labels, rdcol)
breaks = calc.breaks.acut(c(rep(0,7), rep(1,6), rep(2,22)))

png(paste0(imgpref,'matrix_col.png'),res=300,units='in',width=10,height=3)
par(mar=c(.25,3,.25,.25))
layout(matrix(c(1,2),2,1), heights=c(3,1.5), TRUE)
image(t(avail[eporder,cellorder]), col=c('white','darkblue'), axes=F, useRaster=TRUE)
text(y=seq(0,1, length.out=nrow(avail)),
     x=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=eporder, srt=0, adj=1, xpd=TRUE,cex=.3)
text(x=mean(par()$usr[1:2]), y = mean(par()$usr[3:4]), 
     labels = paste0(length(eporder), ' Assays (Marks + Accessibility) by\n', length(cellorder), ' Cell Types/States & Tissues'),
     col='darkgrey', srt=0, adj=.5, xpd=TRUE,cex=1.5)
abline(h=breaks, lty='dotted')
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
par(mar=c(3,3,.25,.25))
# image(faclabels, axes=F, col=colset)
meta.image(metamat[cellorder,5:1], colvals=colvals, cex=.5, horiz=F, useRaster=TRUE)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
text(x=lablist[[1]],
     y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), 
     labels=lablist[[2]], srt=90, adj=1, xpd=TRUE,cex=.4, col=lablist[[3]])
dev.off()

# Add impute vs. observed vs. both vs. none:
impmat = avail * 0
impmat[eporder[1:13],] = 2
o1 = names(which(apply(avail,2,sum) == 1))
impmat[eporder, o1] = impmat[eporder, o1] - 2 * avail[eporder, o1]
finalmat = avail + impmat
# 0 = neither, 1 = obs only, 2 = imp only, 3 = both.

sum(finalmat[,cellorder] %in%  c(1,3))
sum(finalmat[,cellorder] %in% c(2,3))
amarg = apply(avail[,cellorder], 1, sum)

sum(finalmat[,co833] %in%  c(1,3))
sum(efinalmat[,co833] %in%  c(1,3))
sum(finalmat[,co833] %in% c(2,3))
sum(efinalmat[,co833] %in%  c(2,3))
sum(finalmat[,co833] %in% c(2,3))
amarg = apply(avail[,co833], 1, sum)
sum(amarg)
# Remove places where had only one

col.list.h = list()
for (attr in names(colvals)){
    # Add number of values for each name
    acol = metamat[cellorder, attr]
    nv = sapply(names(colvals[[attr]]), function(x, mat=acol){ paste0(x, ' (', sum(mat == x), ')') })
    col.list.h[[attr]] = Legend(labels = nv, legend_gp = gpar(fill = colvals[[attr]], color='black'), title = capitalize(attr), direction='horizontal', nrow=8, 
                                border='black', grid_width=unit(.25, 'cm'), grid_height=unit(.25, 'cm'),
                                title_gp = gpar(fontsize = 8, fontface='bold'), labels_gp = gpar(fontsize = 6))
}
# pd.legend.horiz = packLegend(col.list.h$group, col.list.h$sex, col.list.h$type, col.list.h$lifestage, direction='horizontal', max_height=unit(1,'in'))
pd.legend.horiz = packLegend(col.list.h$project, col.list.h$group, col.list.h$sex, col.list.h$type, col.list.h$lifestage, direction='horizontal', max_height=unit(1,'in'))

ev = sapply(eporder, function(x){ paste0(x, ' (', sum(avail[x,cellorder]), ')') })

png(paste0(imgpref,'matrix_withimp_col_matrixonly_resource.png'),res=600,units='in',width=10,height=3)
par(mar=rep(0,4))
image(t(finalmat[eporder,cellorder]), col=c('white','royalblue', 'lightgrey', 'darkblue'), axes=F, useRaster=TRUE)
dev.off()

png(paste0(imgpref,'matrix_withimp_col_metadata_resource.png'),res=600,units='in',width=10,height=1)
par(mar=rep(0,4))
meta.image(metamat[cellorder,5:1], colvals=colvals, cex=.5, horiz=F, useRaster=TRUE)
dev.off()

png(paste0(imgpref,'matrix_withimp_col.png'),res=300,units='in',width=10,height=3)
par(mar=c(.25,3,.25,.25))
layout(matrix(c(1,2),2,1), heights=c(3,1.5), TRUE)
image(t(finalmat[eporder,cellorder]), col=c('white','royalblue', 'lightgrey', 'darkblue'), axes=F, useRaster=TRUE)
text(y=seq(0,1, length.out=nrow(avail)),
     x=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=ev, srt=0, adj=1, xpd=TRUE,cex=.3)
abline(h=breaks, lty='dotted')
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
par(mar=c(3,3,.25,.25))
meta.image(metamat[cellorder,5:1], colvals=colvals, cex=.5, horiz=F, useRaster=TRUE)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
text(x=lablist[[1]],
     y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), 
     labels=lablist[[2]], srt=90, adj=1, xpd=TRUE,cex=.4, col=lablist[[3]])
draw(pd.legend.horiz, x = unit(5.5,'in'), y=unit(2.925,'in'), just = "top")
dev.off()

pdf(paste0(imgpref,'matrix_withimp_col.pdf'),width=10,height=3)
par(mar=c(.25,3,.25,.25))
layout(matrix(c(1,2),2,1), heights=c(3,1.5), TRUE)
image(t(finalmat[eporder,cellorder]), col=c('white','royalblue', 'lightgrey', 'darkblue'), axes=F, useRaster=TRUE)
text(y=seq(0,1, length.out=nrow(avail)),
     x=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=ev, srt=0, adj=1, xpd=TRUE,cex=.3)
abline(h=breaks, lty='dotted')
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
par(mar=c(3,3,.25,.25))
meta.image(metamat[cellorder,5:1], colvals=colvals, cex=.5, horiz=F, useRaster=TRUE)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
text(x=lablist[[1]],
     y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), 
     labels=lablist[[2]], srt=90, adj=1, xpd=TRUE,cex=.4, col=lablist[[3]])
draw(pd.legend.horiz, x = unit(6,'in'), y=unit(2.9,'in'), just = "top")
dev.off()


# ADDING EXTRA EPITOPES:
# Add impute vs. observed vs. both vs. none:
eimpmat = extavail * 0
eimpmat[exteporder[1:18],] = 2
o1 = names(which(apply(avail,2,sum) == 1))
eimpmat[eporder, o1] = eimpmat[eporder, o1] - 2 * avail[eporder, o1]
efinalmat = extavail + eimpmat
ebreaks = calc.breaks.acut(c(rep(0,7), rep(1,6), rep(3,5), rep(2,22)))
# 0 = neither, 1 = obs only, 2 = imp only, 3 = both.
eev = sapply(exteporder, function(x){ paste0(x, ' (', sum(extavail[x,cellorder]), ')') })


# STATS:
table(efinalmat[,cellorder])
round(table(efinalmat[,cellorder]) / sum(table(efinalmat[,cellorder])) * 100,2)

png(paste0(imgpref,'matrix_withimp_acc_col_matrixonly_resource.png'),res=600, units='in', width=10,height=3)
par(mar=rep(0,4))
image(t(efinalmat[exteporder,cellorder]), col=c('white','royalblue', 'lightgrey', 'darkblue'), axes=F, useRaster=TRUE)
dev.off()

pdf(paste0(imgpref,'matrix_withimp_acc_col.pdf'),width=10,height=3)
# png(paste0(imgpref,'matrix_withimp_acc_col.png'),res=300, units='in', width=10,height=3)
sp=0.15
par(mar=c(sp,3,sp,sp))
layout(matrix(c(1,2),2,1), heights=c(3,1.6), TRUE)
image(t(efinalmat[exteporder,cellorder]), col=c('white','royalblue', 'lightgrey', 'darkblue'), axes=F, useRaster=TRUE)
text(y=seq(0,1, length.out=nrow(extavail)),
     x=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), 
     labels=eev, srt=0, adj=1, xpd=TRUE,cex=.3)
abline(h=ebreaks, lty='dotted')
abline(v=par()$usr[1:2])
abline(h=par()$usr[3:4])
par(mar=c(3.5,3,0,sp))
meta.image(metamat[cellorder,5:1], colvals=colvals, cex=.45, horiz=F, useRaster=TRUE)
abline(v=par()$usr[1:2],lty=1,lw=0.5)
abline(h=par()$usr[3:4],lty=1,lw=0.5)
box.pad = 0.0075
xx = space.1d(lablist[[1]], box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
xx = space.1d(xx, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
rx = c(0.05, 0.125, 0.3, 0.375, 0.425)
x = par()$usr[3]-rx*(diff(par()$usr[3:4]))
text(x=xx, y=x[5], labels=lablist[[2]],
     srt=90, adj=1, xpd=TRUE, 
     cex=0.4, col=lablist[[3]])
par(xpd=TRUE)
segments(y0=x[4], x0=xx, y1=x[3], x1=xx, col=lablist[[3]], lwd=.5)
segments(y0=x[3], x0=xx, y1=x[2], x1=lablist[[1]], col=lablist[[3]], lwd=.5)
segments(y0=x[1], x0=lablist[[1]], y1=x[2], x1=lablist[[1]], col=lablist[[3]], lwd=.5)
par(xpd=FALSE)
draw(pd.legend.horiz, x = unit(5.5,'in'), y=unit(2.9,'in'), just = "top")
dev.off()


# REMOVE: Cell-types that only have ATAC-seq:
ared = avail[,avail['ATAC-seq',] > 0]
atacIDs = names(which(apply(ared,2,sum) == 1))
write.table(atacIDs, 'remove_IDs_ATAC.tsv',
            quote=F, row.names=F, col.names=F, sep='\t')
