#!/usr/bin/R
# -----------------------------------------
# Plot the gwas tree enrichment statistics:
# Tissue-centric analyses/plots:
# -----------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
library(qvalue)

# Arguments for loading data:
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
use.adj = TRUE
use.strict = TRUE

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict) }
source(paste0(bindir, 'load_statistics_gwastree_enrichments.R'))

rm(dflist, cdll)

# --------------------------------------------
# Barplots/boxplots of tissue characteristics:
# --------------------------------------------
# Reduce regdf to GWAS above threshold:
regdf = regdf[regdf$uid %in% keptgwas,]
# Number of enrichments:
tdf = aggregate(log10p ~ COLOR + GROUP + category, regdf, length)
tdf = tdf[order(tdf$GROUP),]
tdf$space = c(0, diff(as.numeric(tdf$category))) * 1 + 0.1
tdf$space[1] = 0
cs = cumsum(tdf$space + 1)
cs = cs - cs[1]

png(paste0(imgpref, 'nenrich_group.png'),res=450,units='in',width=8,height=4)
par(mar=c(5,3,1,.5))
par(yaxs='i')
par(xaxs='i')
# barplot(tdf$log10p, space=tdf$space, width=.5, col=tdf$COLOR, # names.arg=tdf$Project, 
minor = 250
barplot(tdf$log10p, space=tdf$space, width=1, col=tdf$COLOR, # names.arg=tdf$Project, 
        ylim = c(0, (round(max(tdf$log10p) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
abline(h=seq(0,max(tdf$log10p), minor), col='white', lwd=.5, lty='dotted')
abline(h=seq(0,max(tdf$log10p), minor*2), col='white', lwd=2)
txtoffset = diff(par()$usr[1:2]) / (2 * max(cs))
txtloc = (cs /max(cs) * (diff(par()$usr[1:2]) - txtoffset * 2)) + txtoffset
axis(2, cex.axis=.75)
mtext("Total number of enrichments", side=2, cex=.75, line=2)
par(xpd=TRUE)
text(txtloc, tdf$log10p + minor / 10, tdf$log10p, cex=.5, srt=90,adj=0)
text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     tdf$GROUP, adj=1, srt=90,cex=.75)
dev.off()


# Strength of enrichments:
png(paste0(imgpref, 'l10p_enrich_group.png'),res=450,units='in',width=8,height=4)
par(mar=c(5,3,1,.5))
par(yaxs='i')
par(xaxs='i')
ylim = c(0,quantile(regdf$log10p,.999))
atpar = cumsum(c(0, diff(as.numeric(odf$category))) + 1)
boxplot(log10p ~ GROUP, data=regdf,
        at = atpar, width=rep(.01,nrow(odf)), pch=19, cex=.25,
        ylim=ylim, axes=FALSE, col=odf$COLOR, xaxs=FALSE, xlab='')
print(par()$usr)
axis(2, cex.axis=.75)
abline(h=seq(ylim[1],ylim[2],10), col='grey', lwd=.5, lty='dotted')
par(xpd=TRUE)
mtext("-log10 p-value of enrichment", side=2, cex=.75, line=2)
text(atpar, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     odf$GROUP, adj=1, srt=90,cex=.75)
dev.off()


# Number by nodes of enrichments:
ndf = aggregate(log10p ~ node + GROUP + type, regdf, length)
png(paste0(imgpref, 'nenrich_node_group.png'),res=450,units='in',width=8,height=4)
par(mar=c(5,3,1,.5))
par(yaxs='i')
par(xaxs='i')
ylim = c(0,quantile(ndf$log10p,1))
atpar = cumsum(c(0, diff(as.numeric(odf$category))) + 1)
boxplot(log10p ~ GROUP, data=ndf,
        at = atpar, width=rep(.01,nrow(odf)), pch=19, cex=.25,
        ylim=ylim, axes=FALSE, col=odf$COLOR, xaxs=FALSE, xlab='')
print(par()$usr)
axis(2, cex.axis=.75)
abline(h=seq(ylim[1],ylim[2], 20), col='grey', lwd=.5, lty='dotted')
boxplot(log10p ~ GROUP, data=ndf, add=TRUE,
        at = atpar, width=rep(.01,nrow(odf)), pch=19, cex=.25,
        ylim=ylim, axes=FALSE, col=odf$COLOR, xaxs=FALSE, xlab='')
par(xpd=TRUE)
mtext("Number of GWAS per node", side=2, cex=.75, line=2)
text(atpar, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     odf$GROUP, adj=1, srt=90,cex=.75)
dev.off()


# --------------------------------------------------------------------------------------------
# How specific are enrichmets: (total # GWAS hit) and how many nodes at a time: enrichments + 
# --------------------------------------------------------------------------------------------
# Number of unique GWAS hit:
tdf = aggregate(uid ~ COLOR + GROUP + category, regdf, function(x){length(unique(x))})
# tdf = aggregate(log10p ~ COLOR + GROUP + category, regdf, length)
tdf = tdf[order(tdf$GROUP),]
tdf$space = c(0, diff(as.numeric(tdf$category))) * 1 + 0.1
tdf$space[1] = 0
cs = cumsum(tdf$space + 1)
cs = cs - cs[1]

png(paste0(imgpref, 'ngwas_group.png'),res=450,units='in',width=8,height=4)
par(mar=c(5,3,1,.5))
par(yaxs='i')
par(xaxs='i')
minor = 50
barplot(tdf$uid, space=tdf$space, width=1, col=tdf$COLOR,
        ylim = c(0, (round(max(tdf$uid) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
abline(h=seq(0,max(tdf$uid),minor), col='white', lwd=.5, lty='dotted')
abline(h=seq(0,max(tdf$uid),minor *2), col='white', lwd=2)
txtoffset = diff(par()$usr[1:2]) / (2 * max(cs))
txtloc = (cs /max(cs) * (diff(par()$usr[1:2]) - txtoffset * 2)) + txtoffset
axis(2, cex.axis=.75)
mtext("Total number of unique GWAS", side=2, cex=.75, line=2)
par(xpd=TRUE)
text(txtloc, tdf$uid + minor / 10, tdf$uid, cex=.5, srt=90,adj=0)
text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     tdf$GROUP, adj=1, srt=90,cex=.75)
dev.off()


# Number of unique GWAS hit (with assignment to majority)
# Find max signal per tissue - or summed max?
tsdf = aggregate(log10p ~ uid + COLOR + GROUP + category, regdf, sum)
mtsdf = aggregate(log10p ~ uid, tsdf, max)
tsdf = merge(mtsdf, tsdf)
dim(tsdf)
print(length(unique(tsdf$uid)))
ctsdf = aggregate(uid ~ COLOR + GROUP + category, tsdf, function(x){length(unique(x))})
ctsdf = merge(ctsdf, tdf[, c('GROUP','COLOR','category','space')], all.y=T)
ctsdf = ctsdf[order(ctsdf$GROUP),]
pct = round(ctsdf$uid / tdf$uid * 100, 1)
pct[is.na(pct)] = 0

# png(paste0(imgpref, 'ngwas_summajority_group.png'),res=450,units='in',width=8,height=4)
png(paste0(imgpref, 'ngwas_summajority_group.pdf'), width=8,height=4)
par(mar=c(5,3,1,.5))
par(yaxs='i')
par(xaxs='i')
# minor = 50
minor = 25
tcol = sapply(tdf$COLOR, tsp.col)
barplot(tdf$uid, space=tdf$space, width=1, col=tcol,
        ylim = c(0, (round(max(tdf$uid) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
par(new=TRUE)
barplot(ctsdf$uid, space=tdf$space, width=1, col=tdf$COLOR,
        ylim = c(0, (round(max(tdf$uid) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
abline(h=seq(0,max(tdf$uid),minor), col='white', lwd=.5, lty='dotted')
abline(h=seq(0,max(tdf$uid),minor *2), col='white', lwd=2)
txtoffset = diff(par()$usr[1:2]) / (2 * max(cs))
txtloc = (cs /max(cs) * (diff(par()$usr[1:2]) - txtoffset * 2)) + txtoffset
axis(2, cex.axis=.75)
mtext("Total number of unique GWAS", side=2, cex=.75, line=2)
par(xpd=TRUE)
text(txtloc, tdf$uid + minor / 10, paste0(tdf$uid, ' (', pct,'%)'), cex=.5, srt=90,adj=0)
# text(txtloc, tdf$uid + minor / 10, tdf$uid, cex=.5, srt=90,adj=0)
text(txtloc, ctsdf$uid + minor / 10, ctsdf$uid, cex=.5, srt=90,adj=0)
text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     tdf$GROUP, adj=1, srt=90,cex=.75)
dev.off()


# Number of unique GWAS hit (with assignment to majority)
# Find max signal per tissue - or summed max?
tsdf = aggregate(log10p ~ uid + COLOR + GROUP + category, regdf, max)
mtsdf = aggregate(log10p ~ uid, tsdf, max)
tsdf = merge(mtsdf, tsdf)
dim(tsdf)
print(length(unique(tsdf$uid)))
ctsdf = aggregate(uid ~ COLOR + GROUP + category, tsdf, function(x){length(unique(x))})
ctsdf = merge(ctsdf, tdf[, c('GROUP','COLOR','category','space')], all.y=T)
ctsdf = ctsdf[order(ctsdf$GROUP),]
pct = round(ctsdf$uid / tdf$uid * 100, 1)
pct[is.na(pct)] = 0

# png(paste0(imgpref, 'ngwas_majority_group.png'),res=450,units='in',width=8,height=4)
pdf(paste0(imgpref, 'ngwas_majority_group.pdf'),width=8,height=4)
par(mar=c(5,3,1,.5))
par(yaxs='i')
par(xaxs='i')
minor = 50
tcol = sapply(tdf$COLOR, tsp.col)
barplot(tdf$uid, space=tdf$space, width=1, col=tcol,
        ylim = c(0, (round(max(tdf$uid) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
par(new=TRUE)
barplot(ctsdf$uid, space=tdf$space, width=1, col=tdf$COLOR,
        ylim = c(0, (round(max(tdf$uid) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
abline(h=seq(0,max(tdf$uid),minor), col='white', lwd=.5, lty='dotted')
abline(h=seq(0,max(tdf$uid),minor *2), col='white', lwd=2)
txtoffset = diff(par()$usr[1:2]) / (2 * max(cs))
txtloc = (cs /max(cs) * (diff(par()$usr[1:2]) - txtoffset * 2)) + txtoffset
axis(2, cex.axis=.75)
mtext("Total number of unique GWAS", side=2, cex=.75, line=2)
par(xpd=TRUE)
text(txtloc, tdf$uid + minor / 10, paste0(tdf$uid, ' (', pct,'%)'), cex=.5, srt=90,adj=0)
# text(txtloc, tdf$uid + minor / 10, tdf$uid, cex=.5, srt=90,adj=0)
text(txtloc, ctsdf$uid + minor / 10, ctsdf$uid, cex=.5, srt=90,adj=0)
text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     tdf$GROUP, adj=1, srt=90,cex=.75)
dev.off()


# Order majority by ratio:
pord = order(pct, decreasing=TRUE)
gord = as.character(tdf$GROUP[pord])
pdf(paste0(imgpref, 'ngwas_majority_group_reord.pdf'),width=8.5,height=4)
# png(paste0(imgpref, 'ngwas_majority_group_reord.png'),res=450, units='in', width=8.5,height=4)
par(mar=c(5,3,1,2.5))
par(yaxs='i')
par(xaxs='i')
minor = 50
tcol = sapply(tdf$COLOR[pord], tsp.col)
barplot(tdf$uid[pord], space=0.1, width=1, col=tcol,
        ylim = c(0, (round(max(tdf$uid) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
par(new=TRUE)
barplot(ctsdf$uid[pord], space=0.1, width=1, col=tdf$COLOR[pord],
        ylim = c(0, (round(max(tdf$uid) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
abline(h=seq(minor,max(tdf$uid),minor), col='white', lwd=.5, lty='dotted')
# abline(h=seq(minor,max(tdf$uid),minor *2), col='white', lwd=2)
cs2 = cumsum(c(0.1, rep(1.1, (nrow(tdf)-2)),1.1))
txtoffset = diff(par()$usr[1:2]) / (2 * max(cs2))
txtloc = (cs2 / (max(cs2) + min(cs2)) * (diff(par()$usr[1:2]) - txtoffset * 1.5)) + txtoffset
axis(2, cex.axis=.75)
mtext("Total number of unique GWAS", side=2, cex=.75, line=2)
par(xpd=TRUE)
# text(txtloc, tdf$uid[pord] + minor / 10, paste0(tdf$uid, ' (', pct,'%)')[pord], cex=.5, srt=90,adj=0)
text(txtloc, tdf$uid[pord] + minor / 10, tdf$uid[pord], cex=.5, srt=90,adj=0)
# text(txtloc, tdf$uid + minor / 10, tdf$uid, cex=.5, srt=90,adj=0)
text(txtloc, ctsdf$uid[pord] + minor / 10, ctsdf$uid[pord], cex=.5, srt=90,adj=0)
text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     tdf$GROUP[pord], adj=1, srt=90,cex=.75)
par(xpd=FALSE)
par(new=TRUE)
plot(cs2 + 0.5, pct[pord], xaxt='n', yaxt='n',
     xlab='', ylab='', axes=F,
     col=tsp.col('indianred'), pch=19, type='p', cex=.5, lty=1,
     ylim=c(0,70),xlim=c(0.1, max(cs2) + 1), xpd=TRUE)
lines(cs2 + 0.5, pct[pord], col=tsp.col('indianred'))
text(cs2 + 0.5, 60, paste0(pct[pord], '%'), cex=.5, srt=90,adj=0, col='indianred')
atpar=seq(0,70, 10)
axis(4, at=atpar,labels=rep('', length(atpar)), cex.axis=.75, col='indianred')
text(x=par()$usr[2] + 0.025 * diff(par()$usr[1:2]), y=atpar,srt=90, xpd=TRUE,
     labels=paste0(atpar, '%'), col='indianred', cex=.75)
mtext("% Principal Tissue in GWAS", side=4, cex=.75, line=1.25, col='indianred')
dev.off()


# ---------------------------------------
# See #principal/# partner for all nodes:
# ---------------------------------------
tsdf = aggregate(log10p ~ uid + COLOR + GROUP + category + node, regdf, max)
# Groupwide max:
g.mtsdf = aggregate(log10p ~ uid + GROUP, tsdf, max)
names(g.mtsdf) = c('uid','GROUP.agg','max.lp')
g.mtsdf = merge(tsdf, g.mtsdf)
g.mtsdf = filter(g.mtsdf, GROUP != GROUP.agg)
# Max out-group p:
g.mtsdf = aggregate(max.lp ~ uid + COLOR + GROUP + node + log10p + category, g.mtsdf, max)
g.mtsdf = filter(g.mtsdf, log10p > max.lp)
mtsdf1 = merge(tsdf, g.mtsdf[,c('uid','node','log10p')])
# Don't let compete with others in group:
mtsdf = aggregate(log10p ~ uid, tsdf, max)
tsdf = merge(mtsdf, tsdf)
# All:
tsdf = unique(rbind(tsdf, mtsdf1[,colnames(tsdf)]))
tmpdf = tsdf[,c('uid','GROUP')]
names(tmpdf) = c('uid','max.GROUP')
tmprdf = merge(regdf, tmpdf)


# Exclude all where max is within the group:
tmprdf = tmprdf[tmprdf$GROUP != tmprdf$max.GROUP,]
tdf = aggregate(uid ~ COLOR + GROUP + category + node, tmprdf, function(x){length(unique(x))})
tdf = tdf[order(tdf$node),]
names(tdf)[5] = 'uid.any'
# Maximal:

# COUNT ALL PRIMARY OUTSIDE GROUP.
ctsdf = aggregate(uid ~ COLOR + GROUP + category + node, tsdf, function(x){length(unique(x))})
ctsdf = merge(ctsdf, tdf, all.y=T)
ctsdf$uid.any[is.na(ctsdf$uid.any)] = 0
ctsdf$uid[is.na(ctsdf$uid)] = 0
# Add back in the primary:
ctsdf$uid.any = ctsdf$uid + ctsdf$uid.any
ctsdf = ctsdf[order(ctsdf$node),]
ctsdf$pct = round(ctsdf$uid / ctsdf$uid.any * 100, 2)
ctsdf$GROUP = factor(ctsdf$GROUP, levels=gord)
# Order by orig order, then internally:
ctsdf = ctsdf[order(ctsdf$pct, decreasing=T),]
ctsdf = ctsdf[order(ctsdf$GROUP),]

# pdf(paste0(imgpref, 'ngwas_majority_group_reord_all.pdf'),width=20,height=4)
png(paste0(imgpref, 'ngwas_majority_group_reord_all.png'),res=450, units='in', width=20,height=4)
par(mar=c(5,3,1,2.5))
par(yaxs='i')
par(xaxs='i')
tcol = sapply(ctsdf$COLOR, tsp.col)
barplot(ctsdf$uid.any, space=0.01, width=1, col=tcol,
        ylim = c(0, (round(max(ctsdf$uid.any) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
par(new=TRUE)
minor = 5
barplot(ctsdf$uid, space=0.01, width=1, col=ctsdf$COLOR,
        ylim = c(0, (round(max(ctsdf$uid) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
abline(h=seq(minor,max(ctsdf$uid.any),minor), col='white', lwd=.5, lty='dotted')
# abline(h=seq(minor,max(tdf$uid),minor *2), col='white', lwd=2)
cs2 = cumsum(c(0.01, rep(1.01, (nrow(ctsdf)-2)),1.01))
txtoffset = diff(par()$usr[1:2]) / (2 * max(cs2))
txtloc = (cs2 / (max(cs2) + min(cs2)) * (diff(par()$usr[1:2]) - txtoffset * 1.5)) + txtoffset
axis(2, cex.axis=.75)
mtext("Total number of unique GWAS", side=2, cex=.75, line=2)
par(xpd=TRUE)
# text(txtloc, tdf$uid[pord] + minor / 10, paste0(tdf$uid, ' (', pct,'%)')[pord], cex=.5, srt=90,adj=0)
text(txtloc, ctsdf$uid.any + minor / 10, ctsdf$uid.any, cex=.2, srt=90,adj=0)
# text(txtloc, tdf$uid + minor / 10, tdf$uid, cex=.5, srt=90,adj=0)
text(txtloc, ctsdf$uid + minor / 10, ctsdf$uid, cex=.2, srt=90,adj=0)
text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     paste0(ctsdf$GROUP, ' - ', ctsdf$node), adj=1, srt=90,cex=.2)
par(xpd=FALSE)
par(new=TRUE)
plot(cs2 + 0.5, ctsdf$pct, xaxt='n', yaxt='n', xlab='', ylab='', axes=F,
     col=tsp.col('indianred'), pch=19, type='p', cex=.25, lty=1,
     ylim=c(0,100),xlim=c(0.01, max(cs2) + 1), xpd=TRUE)
lines(cs2 + 0.5, ctsdf$pct, col=tsp.col('indianred'))
text(cs2 + 0.5, 90, paste0(pct[pord], '%'), cex=.1, srt=90,adj=0, col='indianred')
atpar=seq(0,100, 25)
axis(4, at=atpar,labels=rep('', length(atpar)), cex.axis=.75, col='indianred')
text(x=par()$usr[2] + 0.01 * diff(par()$usr[1:2]), y=atpar,srt=90, xpd=TRUE,
     labels=paste0(atpar, '%'), col='indianred', cex=.75)
mtext("% Principal Tissue in GWAS", side=4, cex=.75, line=1.25, col='indianred')
dev.off()






# -----------------------------------------------------------
# - Number of hits for GWAS assigned to each tissue by max pval
# - Number of hits for GWAS assigned to each tissue if tissue involved.
# -----------------------------------------------------------
# Number of unique GWAS hit (with assignment to majority)
# Find max signal per tissue - or summed max?
tsdf = aggregate(log10p ~ uid + COLOR + GROUP + category, regdf, max)
mtsdf = aggregate(log10p ~ uid, tsdf, max)
tsdf = merge(mtsdf, tsdf)

nhdf = aggregate(GROUP ~ uid, regdf, function(x) length(unique(x)))
names(nhdf)[2] = 'nhit'
nhdf = merge(nhdf, tsdf) 

gp = ggplot(nhdf, aes(GROUP, nhit, fill=GROUP)) + geom_boxplot() +
    geom_jitter(cex=.1) +
    scale_fill_manual(values=odf$COLOR)
ggsave('~/test.png',gp, dpi=450, units='in', width=12, height=6)

nhdf2 = aggregate(GROUP ~ uid, regdf, function(x) length(unique(x)))
names(nhdf2)[2] = 'nhit'
nhdf2 = merge(nhdf2, regdf) 

gp = ggplot(nhdf2, aes(GROUP, nhit, fill=GROUP)) +
    geom_boxplot() +
    geom_jitter(cex=.1) +
    scale_fill_manual(values=odf$COLOR)
ggsave('~/test.png',gp, dpi=450, units='in', width=12, height=6)




mtsdf = aggregate(log10p ~ uid, tsdf, max)
tsdf = merge(mtsdf, tsdf)
dim(tsdf)
print(length(unique(tsdf$uid)))
ctsdf = aggregate(uid ~ COLOR + GROUP + category, tsdf, function(x){length(unique(x))})
ctsdf = merge(ctsdf, tdf[, c('GROUP','COLOR','category','space')], all.y=T)
ctsdf = ctsdf[order(ctsdf$GROUP),]
pct = round(ctsdf$uid / tdf$uid * 100, 1)
pct[is.na(pct)] = 0







# -----------------------------------
# Number of nodes hit per unique GWAS
# -----------------------------------
regdf$enid = 1:nrow(regdf) # For unique enrichments
tdf = aggregate(cbind(uid, enid) ~ COLOR + GROUP + category, regdf, function(x){length(unique(x))})
tdf$ratio = tdf$enid / tdf$uid
tdf = tdf[order(tdf$GROUP),]
tdf$space = c(0, diff(as.numeric(tdf$category))) * 1 + 0.1
tdf$space[1] = 0
cs = cumsum(tdf$space + 1)
cs = cs - cs[1]

png(paste0(imgpref, 'nnode_avg_group.png'),res=450,units='in',width=8,height=4)
par(mar=c(5,3,1,.5))
par(yaxs='i')
par(xaxs='i')
minor = 1
barplot(tdf$ratio, space=tdf$space, width=1, col=tdf$COLOR,
        ylim = c(0, (round(max(tdf$ratio) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
abline(h=seq(0,max(tdf$ratio),minor), col='white', lwd=.5, lty='dotted')
abline(h=seq(0,max(tdf$ratio),minor * 2), col='white', lwd=2)
txtoffset = diff(par()$usr[1:2]) / (2 * max(cs))
txtloc = (cs /max(cs) * (diff(par()$usr[1:2]) - txtoffset * 2)) + txtoffset
axis(2, cex.axis=.75)
mtext("Number of nodes per GWAS (if GWAS hits tissue)", side=2, cex=.75, line=2)
par(xpd=TRUE)
text(txtloc, tdf$ratio + minor / 10, round(tdf$ratio,2), cex=.5, srt=90, adj=0)
text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     tdf$GROUP, adj=1, srt=90,cex=.75)
dev.off()


# NORMALIZE BY # of nodes:
tdf = merge(tdf, nodetissue.stats)
tdf = tdf[order(tdf$GROUP),]
tdf$pct = tdf$ratio / tdf$node * 100

png(paste0(imgpref, 'pctnode_avg_group.png'),res=450,units='in',width=8,height=4)
par(mar=c(5,3,1,.5))
par(yaxs='i')
par(xaxs='i')
minor = 5
barplot(tdf$pct, space=tdf$space, width=1, col=tdf$COLOR,
        ylim = c(0, (round(max(tdf$pct) / minor) + 1) * minor),
        axes=F, border=NA, horiz=F)
abline(h=seq(0,max(tdf$pct), minor), col='white', lwd=.5, lty='dotted')
abline(h=seq(0,max(tdf$pct), minor *2), col='white', lwd=2)
txtoffset = diff(par()$usr[1:2]) / (2 * max(cs))
txtloc = (cs /max(cs) * (diff(par()$usr[1:2]) - txtoffset * 2)) + txtoffset
axis(2, cex.axis=.75)
mtext("Percent of nodes per GWAS (if GWAS hits tissue)", side=2, cex=.75, line=2)
par(xpd=TRUE)
text(txtloc, tdf$pct + minor / 10, paste0(round(tdf$pct,2),'%'), cex=.5, srt=90, adj=0)
text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
     tdf$GROUP, adj=1, srt=90,cex=.75)
dev.off()



plot.byleaf = FALSE
if (plot.byleaf){
    # By leaves:
    # Number of enrichments:
    tdf = aggregate(log10p ~ COLOR + GROUP + category + type, regdf, length)
    tdf = merge(tdf, merge(odf,
                           data.frame(type=c('Branch','Leaf'))), all.y=T)
    tdf$log10p[is.na(tdf$log10p)] = 0
    tdf = tdf[order(tdf$GROUP),]
    tdf$space = c(0, diff(as.numeric(tdf$category))) * 1 + 0.1
    tdf$space[1] = 0
    cs = cumsum(tdf$space + 1)
    cs = cs - cs[1]
    tdf$label = paste0(tdf$GROUP, ' - ', tdf$type)
    tdf$label[tdf$type == 'Leaf'] = 'Leaf'

    png(paste0(imgpref, 'nenrich_group_nodetype.png'),res=450,units='in',width=8,height=4)
    par(mar=c(4.5,4.5,1,1))
    par(yaxs='i')
    par(xaxs='i')
    barplot(tdf$log10p, space=tdf$space, width=1, col=tdf$COLOR, # names.arg=tdf$Project, 
            ylim = c(0, (round(max(tdf$log10p) / 100) + 1) * 100),
            ylab='Total number of enrichments', border=NA, horiz=F)
    abline(h=seq(0,max(tdf$log10p),250), col='white', lwd=.5, lty='dotted')
    abline(h=seq(0,max(tdf$log10p),500), col='white', lwd=2)
    txtoffset = diff(par()$usr[1:2]) / (2 * max(cs))
    txtloc = (cs /max(cs) * (diff(par()$usr[1:2]) - txtoffset * 2)) + txtoffset
    text(txtloc, tdf$log10p + 20, tdf$log10p, cex=.35)
    par(xpd=TRUE)
    text(txtloc, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
         tdf$label, adj=1, srt=90,cex=.35)
    dev.off()

    # Strength:
    png(paste0(imgpref, 'l10p_enrich_group_nodetype.png'),res=450,units='in',width=8,height=4)
    par(mar=c(5,3,1,1))
    par(yaxs='i')
    par(xaxs='i')
    atpar = cumsum(c(0, diff(as.numeric(tdf$category))) + 1)
    ylim = c(0,quantile(regdf$log10p,.999))
    boxplot(log10p ~ GROUP + type, data=regdf,
            at = atpar, width=rep(.01,nrow(tdf)), pch=19, cex=.25,
            ylim=ylim, axes=FALSE, col=tdf$COLOR, xaxs=FALSE)
    print(par()$usr)
    axis(2, cex.axis=.5)
    abline(h=seq(ylim[1],ylim[2],10), col='grey', lwd=.5, lty='dotted')
    par(xpd=TRUE)
    mtext("-log10 p-value of enrichment", side=2, cex=.75, line=2)
    text(atpar, par()$usr[3] - 0.01 * (diff(par()$usr[3:4])), 
         tdf$label, adj=1, srt=90,cex=.5)
    dev.off()
}


# ------------------------------
# Measures of tissue specificity
# ------------------------------
# Establish how specific certain signals are in GWAS
# How often a tissue has the top signal, etc.
# ------------------------------
# Transformation matrix - node to leaf
tmat = t(sapply(1:NN, function(x){labels(dend3) %in% declist$dec[[x]]}))
ntmat = 1 * t(sapply(1:NN, function(x){odf$GROUP %in% nodetissue$GROUP[nodetissue$node == x]}))
ltmat = 1 * t(sapply(1:NL, function(x){odf$GROUP %in% meta$GROUP[meta$id == labels(dend)[x]]}))
colnames(tmat) = labels(dend3)
colnames(ntmat) = odf$GROUP
colnames(ltmat) = odf$GROUP

# Make matrices (node, leaf) x GWAS:
cutoff = 3
kw = keptgwas[keptgwas %in% rownames(all.regmat)]
nodemat = 1 * (all.regmat[kw,] > cutoff)
nodemat = nodemat[guid,]
lnodemat = all.regmat[kw,]
lnodemat[lnodemat < cutoff] = 0
lnodemat = lnodemat[apply(lnodemat, 1, sum) > 0,]
# Remove all w/out signif:
# nodemat = 1 * (all.regmat >= 2)
leafmat = nodemat %*% tmat
# Turn leaf or node matrices to tissues:
nodetismat = nodemat %*% ntmat
lnodetismat = lnodemat %*% ntmat
leaftismat = leafmat %*% ltmat
nodetismat[nodetismat > 0] = 1
leaftismat[leaftismat > 0] = 1

# Calculate the number of intersections:
# NOTE: On node level not as pre-organized:
dnj = dist(t(nodemat), 'Jaccard')
# image(as.matrix(dnj), col=colspec)
dt = dnj
mat = as.matrix(dt)
ndt = rep("", NN)
ndt[nodetissue$node] = nodetissue$GROUP
ndt = make.unique(ndt)
rownames(mat) = ndt
colnames(mat) = ndt
ht <- hclust(dt, method='ward.D')
cocl <- order.optimal(dt, ht$merge)$order
mat = mat[cocl, cocl]
# NOTE: Order first by node type, then cocl
rn = sapply(rownames(mat), function(x){ sub("\\.[0-9]*$","", x) })
rn = factor(rn, levels=odf$GROUP)
ord = order(rn)
mat = mat[ord,ord]
mat[mat == 0] =NA
colramp=rev(colred)
colramp=colrb
colramp=colspec
diag(mat) = NA
# Plotting function:
zlim = range(mat, na.rm=T)
par(mar=c(.5,10,.5,.5))
par(yaxs='i')
par(xaxs='i')
minor = 1
image(mat, axes=F, col=colramp, useRaster=TRUE, zlim=zlim)

# Plot for sub cat:
for (cat in unique(odf$category)){
    cind = which(rn %in% odf$GROUP[odf$category == cat])
    cind = names(rn[cind])
    submat = mat[cind,cind]
    png(paste0(imgpref, 'node_', sub(" ", "_", cat), '_jaccard.png'),res=450,units='in',width=10,height=8)
    par(mar=c(.5,10,.5,.5))
    par(yaxs='i')
    par(xaxs='i')
    minor = 1
    image(submat, axes=F, col=colramp, useRaster=TRUE, zlim=zlim)
    yt = seq(0,1,length.out=ncol(submat))
    xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
    yaxlab=rownames(submat)
    text(y=yt, x=xt, labels=yaxlab, srt=0, adj=1, xpd=TRUE,cex=1 * 34 / ncol(submat))
    dev.off()
}


ht <- hclust(dnj, method='ward.D')
cocl <- order.optimal(dnj, ht$merge)$order
dnj = as.matrix(dnj)[cocl,cocl]
image(dnj, col=colspec)
# NOTE: - looks really weird for some??
# NOTE: Probably want to evaluate on the node level 
# or aggressively correct the leaf level.
# look just at one tissue group
# On an aggregate, then tissue basis?

plot.dt.sym = function(dt, nbreak=NULL){
    mat = as.matrix(dt)
    ht <- hclust(dt, method='ward.D')
    cocl <- order.optimal(dt, ht$merge)$order
    rn <- names(cocl)[cocl]
    mat = mat[rn,rn]
    colramp=rev(colred)
    colramp=colrb
    colramp=colspec
    diag(mat) = NA
    # Plotting function:
    par(mar=c(.5,10,.5,.5))
    par(yaxs='i')
    par(xaxs='i')
    minor = 1
    # image(mat, axes=F, zlim=c(0,1), col=colramp, useRaster=TRUE)
    image(mat, axes=F, col=colramp, useRaster=TRUE)
    yt = seq(0,1,length.out=ncol(mat))
    xt = par()$usr[1] - 0.01*diff(par()$usr[1:2])
    yaxlab=rownames(mat)
    text(y=yt, x=xt, labels=yaxlab, srt=0, adj=1, xpd=TRUE,cex=1)
    # Generate breaks:
    if (!(is.null(nbreak))){ 
        breaks = calc.breaks(ht, nbreak, cocl)
        abline(v=breaks,lty=2,lwd=.5)
        abline(h=breaks,lty=2,lwd=.5)
    }
}

# Figure without corrections, for leaves:
dlj = dist(t(leaftismat), 'Jaccard')
png(paste0(imgpref, 'leaf_tissue_jaccard.png'),res=450,units='in',width=10,height=8)
plot.dt.sym(dlj, nbreak=5)
dev.off()

# Show all nodes:
# nodemat 
dnj = dist(t(nodetismat), 'Jaccard')
png(paste0(imgpref, 'node_tissue_jaccard.png'),res=450,units='in',width=10,height=8)
plot.dt.sym(dnj, nbreak=5)
dev.off()

# Filter if one of top 3 
lnodetismat = lnodemat %*% ntmat
topn = 5
tnodetismat = apply(lnodetismat, 1, function(x, top=topn, cutoff=3){
                        topn = sort(x,decreasing=T)[top]
                        topn = max(topn, cutoff)
                        1 * (x >= topn) })
tnodetismat = t(tnodetismat)
sum(lnodetismat > 0)
sum(tnodetismat > 0)

# Show all nodes:
# nodemat 
dtj = dist(t(tnodetismat), 'Jaccard')
plot.dt.sym(dtj, nbreak=5)
png(paste0(imgpref, 'node_tissue_jaccard_top', topn,'.png'),res=450,units='in',width=10,height=8)
plot.dt.sym(dtj, nbreak=5)
dev.off()

# Looking at cosine dist instead?
lnodetismat = sweep(lnodetismat, 1, apply(lnodetismat, 1, sum), '/')
dne = dist(t(lnodetismat), 'cosine')
plot.dt.sym(dne, nbreak=5)

metric = 'cosine'
rnodetismat = nodemat %*% ntmat # [nid,]
dne = dist(t(rnodetismat), metric)
plot.dt.sym(dne, nbreak=5)

NT = nrow(dnj)

# Overly optimistic, doesn't take structure into account:
calc.indpt = FALSE
if (calc.indpt){
    library(jaccard)
    # Calc signif with jaccard test:
    pmat = matrix(NA, nrow=NT, ncol=NT)
    for (i in 1:(NT-1)){
        print(i)
        for (j in (i+1):NT){
            pmat[i,j] = jaccard.test.exact(nodetismat[,i], nodetismat[,j], verbose=FALSE)$pvalue
        }
    }
    rownames(pmat) = colnames(nodetismat)
    colnames(pmat) = colnames(nodetismat)
    # Adjust:
    pwide = data.frame(pmat)
    pwide$GROUP = rownames(pmat)
    plong = gather(pwide, GROUP2, pvalue, -GROUP)
    plong = plong[!is.na(plong$pvalue) ,]
    plong = plong[plong$GROUP != plong$GROUP2,]
    plong$padj = p.adjust(plong$pvalue)
    hist(plong$padj, 30)
    plong = plong[order(plong$padj), ]
    print(sum(plong$padj < 0.01))
    # This indicates ~260 interactions significant! Too much to consider.
    # On the other hand, our conservative permutations yield very few
    lpmat = -log10(pmat)
    # Tries to remove expectation, as above.
    dnj = dist(t(nodetismat), 'Jaccard')
    mg = apply(nodetismat, 2, mean)
    rsm = matrix(rep(mg, NT), nrow=NT,ncol=NT, byrow=T)
    ism = rsm * t(rsm)
    usm = rsm + t(rsm)
    jsm = ism / (usm - ism)
    diffnj = (as.matrix(dnj) + jsm)
    png(paste0(imgpref, 'node_tissue_jaccard_diff.png'),res=450,units='in',width=10,height=8)
    plot.dt.sym(as.dist(diffnj), nbreak=5)
    dev.off()
}

# NOTE: NODE-wise is much more accurate, other is confounded by tree structure. 
# ----------------------------------------------
# Perform permutations to evaluate significance:
# ----------------------------------------------
take.topn = function(x, top=topn, cutoff=3){
    topn = sort(x,decreasing=T)[top]
    topn = max(topn, cutoff)
    1 * (x >= topn) }

nodemat = 1 * (all.regmat[kw,] >= cutoff)
nodemat = nodemat[apply(nodemat, 1, sum) > 0,]
# Get the non-hit nodes:
nid = which(apply(nodemat, 2, sum) > 0)
gid = which(apply(nodemat, 1, sum) > 0)
nodemat = nodemat[gid, nid]

rnodetismat = nodemat %*% ntmat[nid,]
dne = dist(t(rnodetismat), 'cosine')
plot.dt.sym(dne, nbreak=5)

mg = apply(as.matrix(dne), 1, mean)
rsm = matrix(rep(mg, NT), nrow=NT,ncol=NT, byrow=T)
mdne = rsm * t(rsm)

diffne = as.dist(as.matrix(dne) - rsm)
# dne = dist(t(rnodetismat), 'cosine')
plot.dt.sym(diffne, nbreak=5)

nwide = data.frame(nodemat)
nwide$uid = rownames(nodemat)
nlong = gather(nwide, node, value, -uid)
nlong = nlong[nlong$value != 0,]
nr = nrow(nlong)
plong = nlong
cnode = unique(nlong$node)
ruid = unique(nlong$uid)

mr = apply(nodemat, 1, sum)
mc = apply(nodemat, 2, sum)
rnodetismat = nodemat %*% ntmat[nid,]
mtr = apply(rnodetismat, 1, sum)
mtc = apply(rnodetismat, 2, sum)

# If y is a vector:
jacc.dist.vecy = function(x, y=NULL, verbose=TRUE){
    if (is.null(y)){ 
        if (verbose){ print("Using t(x) as second matrix") }
        y = t(x) 
    }
    rsx = apply(x, 1, sum)
    rsy = apply(y, 2, sum)
    NTX = length(rsx)
    NTY = length(rsy)
    rsx = matrix(rep(rsx, NTY), nrow=NTX,ncol=NTY, byrow=T)
    rsy = matrix(rep(rsy, NTX), nrow=NTY,ncol=NTX, byrow=T)
    rsm = rsx + t(rsy)
    nint = x %*% y
    nunion = rsm - nint
    dist = 1.0 - nint / nunion 
    return(dist)
}

use.asym = TRUE
# metric = 'cosine'
# NPERM = 1000
metric = 'jaccard'
topn = 0
NPERM = 10000
permfile = paste0(eprefix, 'perm_node_tissue_', metric, '_t', topn, '_nperm', NPERM, suffix)
if (!file.exists(permfile)){
    if (metric == 'cosine'){
        rnodetismat = nodemat %*% ntmat[nid,]
        fixed.dt = as.matrix(dist(t(rnodetismat), metric))
    } else if (metric == 'jaccard'){
        if (topn != 0){
            # Filter if one of top topn
            lnodetismat = lnodemat %*% ntmat
            tnodetismat = apply(lnodetismat, 1, top=topn, take.topn)
            nodetismat = t(tnodetismat)
            fixed.dt = as.matrix(dist(t(nodetismat), metric))
            # mtr = apply(nodetismat, 1, sum)
            # mtc = apply(nodetismat, 2, sum)
        } else {
            nodetismat = nodemat %*% ntmat[nid,]
            nodetismat[nodetismat > 0] = 1
            fixed.dt = as.matrix(dist(t(nodetismat), metric))
        }
    }
    rnodetismat = nodemat %*% ntmat[nid,]
    mtr = apply(rnodetismat, 1, sum)
    mtc = apply(rnodetismat, 2, sum)
    perm.dt = fixed.dt
    perm.dt[] = 0
    for (t in 1:ncol(nodetismat)){
        tmtr = apply(rnodetismat[,-t], 1, sum)
        pb = txtProgressBar(min=0, max=NPERM, style = 3)
        for (i in 1:NPERM){
            setTxtProgressBar(pb, i)
            set.seed(i)
            # Permute all but the fixed tissue:
            # 1alt. Permute GWAS x tissue matrix instead:
            permtismat = r2dtable(1, tmtr, mtc[-t])[[1]]
            if (topn != 0){
                permtismat = t(apply(permtismat, 1,
                                     function(x, top=topn, cutoff=0.1){
                                         x = x + jitter(x, .001)
                                         topn = sort(x,decreasing=T)[top]
                                         topn = max(topn, cutoff)
                                         1 * (x >= topn) }))
            }
            # 3. Calculate the new distance matrix (asymmetrically)
            if (metric == 'cosine'){
                dt = cosine.dist(x=t(rnodetismat[,t]), y=permtismat)
            } else if (metric == 'jaccard'){
                permtismat[permtismat > 0] = 1
                dt = jacc.dist.vecy(x=t(permtismat), y=t(t(nodetismat[,t])))
            }
            # perm.dt[-t,t] = perm.dt[-t,t] + 1 * (dt < fixed.dt[-t,t])
            perm.dt[t,-t] = perm.dt[t,-t] + 1 * (t(dt) < fixed.dt[t,-t])
        }
    }
    close(pb)
    perm.dt = data.frame(perm.dt)
    perm.dt$G1 = rownames(perm.dt)
    permdf = gather(perm.dt, G2, below, -G1)
    permdf$pvalue = permdf$below / NPERM
    save(permdf, file=permfile)
} else {
    load(permfile)
}

# ------------------------------------------
# Compare permutations to the original data:
# NOTE: May need to correct for (33*34/2) tests.
# ------------------------------------------
g2map = data.frame(G2w=unique(permdf$G2), G2=odf$GROUP)
names(permdf)[2] = 'G2w'

# Alternatively: Evaluate number below for each:
if (metric == 'cosine'){
    rnodetismat = nodemat %*% ntmat[nid,]
    dt = as.matrix(dist(t(rnodetismat), metric))
} else if (metric == 'jaccard'){
    nodetismat = nodemat %*% ntmat[nid,]
    nodetismat[nodetismat > 0] = 1
    dt = as.matrix(dist(t(nodetismat), metric))
}
fdt = dt
dt = data.frame(dt)
colnames(dt) = rownames(dt)
dt$G1 = rownames(dt)
dt = gather(dt, G2, metric, -G1)
dt = dt[dt$G1 != dt$G2,]
dt = merge(dt, g2map)
dt = merge(dt, permdf)
# dt$pvalue =  qvalue(dt$pvalue)$qvalue

# Get ordering
dn = as.dist(fdt)
ht <- hclust(dn, method='ward.D')
cocl <- order.optimal(dn, ht$merge)$order
groupord <- names(cocl)[cocl]
# TEST:
dt$G1 = factor(dt$G1, levels=groupord)
dt$G2 = factor(dt$G2, levels=groupord)
dt$N1 = as.numeric(dt$G1)
dt$N2 = as.numeric(dt$G2)

dt$cex = 1
# dt$cex[dt$pvalue < 0.1 | dt$pvalue > 1 - 0.1] = 2
# dt$cex[dt$pvalue < 0.05 | dt$pvalue > 1 - 0.05] = 3
# dt$cex[dt$pvalue < 0.01 | dt$pvalue > 1 - 0.01] = 4
dt$cex[dt$pvalue < 0.1] = 2
dt$cex[dt$pvalue < 0.05] = 3
dt$cex[dt$pvalue < 0.01] = 4
dt$border = 'black'
dt$border[dt$pvalue >= 0.1] = 'white'
dt$label = ''
dt$jacc.label[dt$cex >= 4] = round(100 * (1 - dt$metric[dt$cex >= 4]),1)
# Label as number of intersected GWAS:
dt$label[dt$cex >= 4] = round(100 * (1 - dt$metric[dt$cex >= 4]),1)

nint = t(nodetismat) %*% nodetismat
nintdf = data.frame(nint)
nintdf$G1 = rownames(nintdf)
nintdf = gather(nintdf, G2, nint, -G1)
names(nintdf)[2] = 'G2w'
nintdf = merge(nintdf, g2map)
dt = merge(dt, nintdf)
dt$label[dt$cex >= 4] = dt$nint[dt$cex >= 4]


if (metric == 'cosine'){ 
    # palette = rev(viridis(100)) 
    # palette = rev(inferno(100)) 
    palette = rev(col1)
    txtcol = 'white'
} else { 
    # palette = colspec 
    palette = rev(col1)
    txtcol = 'black'
}
dt$bins <- cut(dt$metric, seq(min(dt$metric), 1, length.out=length(palette)), include.lowest=T) 
dt$col = palette[dt$bins]
linedf =data.frame(y0=1:NT, y1=1:NT, x0=1, x1=(1:NT)-1)
linedf = linedf[-1,]
# Get colors:
rownames(odf) = odf$GROUP
gcol = odf[groupord, 'COLOR']

# Update palette and legend for pvalues:
col_fun = function(x, pal=palette){
    palette = rev(pal)
    bin <- cut(x, seq(0, max(1-dt$metric) * 100, length.out=length(palette)), include.lowest=T) 
    palette[bin]
}

metric.legend = Legend(at = seq(0, max(1 - dt$metric) * 100, 10), 
                     col_fun=col_fun, title_position = "topleft", title = paste(capitalize(metric), "Similarity (%)"), direction = 'horizontal')
# Add the pvalue boxes as well
dt$cex[dt$pvalue < 0.1] = 2
dt$cex[dt$pvalue < 0.05] = 3
dt$cex[dt$pvalue < 0.01] = 4
plegend = packLegend(metric.legend)

if (use.asym) {
    subdt = dt
    subdt$X = subdt$N2
    subdt$Y = subdt$N1
    xlim=c(0, NT+1)
    ylim=c(0, NT+1)
} else {
    subdt = dt[dt$N1 < dt$N2,]
    subdt$X = subdt$N1
    subdt$Y = subdt$N2
    xlim=c(0,NT)
    ylim=c(1, NT+1)
}

# png(paste0(imgpref, 'node_tissue_', metric, '_t', topn, '_signif.png'),res=450,units='in',width=11,height=11.5)
pdf(paste0(imgpref, 'node_tissue_', metric, '_t', topn, '_signif.pdf'), width=11,height=11.5)
par(mar=c(2.5, 7, 7, 0.5))
par(yaxs='i')
par(xaxs='i')
minor = 1
plot(subdt$X, subdt$Y,
     pch=22, cex=subdt$cex, bg=subdt$col, axes=F,
     xlab='', ylab='', col=subdt$border, xlim=xlim, ylim=ylim)
if (use.asym){
    points(x=rep(0, NT), y=1:NT, pch=22, cex=4, col='white', bg=gcol)
    points(x=1:NT, y=rep(NT + 1, NT), pch=22, cex=4, col='white', bg=gcol)
} else {
    points(rep(0,NT-1), 2:NT, pch=22, cex=4, col='white', bg=gcol[-1])
    points(1:(NT-1), rep(NT + 1, NT-1), pch=22, cex=4, col='white', bg=gcol[-NT])
}
# Add the colors:
segments(x0=linedf$x0, x1=linedf$x1,
         y0=linedf$y0, y1=linedf$y1, lty='dashed', lwd=.5)
segments(y0=NT, y1=rev(NT - linedf$x1 + 1),
         x0=linedf$y0 - 1, x1=linedf$y1 - 1, lty='dashed', lwd=.5)
if (use.asym){
    # Flip all by subtracting lines from NT+1
    segments(x0=NT+1-linedf$x0, x1=NT+1-linedf$x1,
             y0=NT+1-linedf$y0, y1=NT+1-linedf$y1, lty='dashed', lwd=.5)
    segments(y0=NT+1-NT, y1=NT+1-rev(NT - linedf$x1 + 1),
             x0=NT+1-(linedf$y0 - 1), x1=NT+1-(linedf$y1 - 1), lty='dashed', lwd=.5)
}
par(new=TRUE)
plot(subdt$X, subdt$Y,
     pch=22, cex=subdt$cex, bg=subdt$col, axes=F,
     xlab='', ylab='', col=subdt$border, xlim=xlim, ylim=ylim)
# text(subdt$X, subdt$Y, labels=subdt$label, cex=.6, col=txtcol)
text(subdt$X, subdt$Y, labels=subdt$label, cex=.6, col=ifelse(subdt$metric < .7, 'white','black'))
loc = 1:NT
xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
yt = par()$usr[4] + 0.005*diff(par()$usr[3:4])
if (use.asym){
    text(y=loc, x=xt, labels=groupord, srt=0, adj=1, xpd=TRUE, cex=1)
    text(y=yt, x=loc, labels=groupord, srt=90, adj=0, xpd=TRUE, cex=1)
} else {
    text(y=loc[-1], x=xt, labels=groupord[-1], srt=0, adj=1, xpd=TRUE, cex=1)
    text(y=yt, x=loc[-NT], labels=groupord[-NT], srt=90, adj=0, xpd=TRUE, cex=1)
}
par(xpd=TRUE)
legend('bottom', legend=c('NS','p < 0.1','p < 0.05','p < 0.01'), title='p-value',
       col='black', fill='white', pch=22, pt.cex=c(1,2,3,4), bty='n', border=NA, cex=1, horiz=T, inset=c(0,-.04))
draw(plegend, x = unit(9.5,'in'), y=unit(0.6,'in'), just = "top")
dev.off()



# --------------------------
# Plot without significance.
# --------------------------
subdt = dt[dt$N1 < dt$N2,]
subdt$X = subdt$N1
subdt$Y = subdt$N2
xlim=c(0, NT)
ylim=c(1, NT+1)
subdt$cex = 4
subdt$border =NA
subdt$text.col = 'black'
subdt$text.col[subdt$metric < 0.8] = 'white'

png(paste0(imgpref, 'node_tissue_', metric, '_t', topn, '.png'),res=450,units='in',width=11,height=11.5)
par(mar=c(2.5, 7, 7, 0.5))
par(yaxs='i')
par(xaxs='i')
minor = 1
plot(subdt$X, subdt$Y,
     pch=22, cex=subdt$cex, bg=subdt$col, axes=F,
     xlab='', ylab='', col=subdt$border, xlim=xlim, ylim=ylim)
points(rep(0,NT-1), 2:NT, pch=22, cex=4, col='white', bg=gcol[-1])
points(1:(NT-1), rep(NT + 1, NT-1), pch=22, cex=4, col='white', bg=gcol[-NT])
par(new=TRUE)
plot(subdt$X, subdt$Y,
     pch=22, cex=subdt$cex, bg=subdt$col, axes=F,
     xlab='', ylab='', col=NA, xlim=xlim, ylim=ylim)
text(subdt$X, subdt$Y, labels=subdt$nint, cex=.8, col=subdt$text.col)
loc = 1:NT
xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
yt = par()$usr[4] + 0.005*diff(par()$usr[3:4])
text(y=loc[-1], x=xt, labels=groupord[-1], srt=0, adj=1, xpd=TRUE, cex=1)
text(y=yt, x=loc[-NT], labels=groupord[-NT], srt=90, adj=0, xpd=TRUE, cex=1)
par(xpd=TRUE)
points(1:NT, 1:NT, pch=22, cex=3, col='white', bg=gcol)
# legend('bottom', legend=c('NS','p < 0.1','p < 0.05','p < 0.01'), title='p-value',
#        col='black', fill='white', pch=22, pt.cex=c(1,2,3,4), bty='n', border=NA, cex=1, horiz=T, inset=c(0,-.04))
draw(plegend, x = unit(9.5,'in'), y=unit(0.6,'in'), just = "top")
dev.off()



# ----------------------------------
# Plot in the order of the barplots:
# ----------------------------------
groupord = rev(odf$GROUP)
dt$G1 = factor(dt$G1, levels=groupord)
dt$G2 = factor(dt$G2, levels=groupord)
dt$N1 = as.numeric(dt$G1)
dt$N2 = as.numeric(dt$G2)
# Get colors and locations:
rownames(odf) = odf$GROUP
gcol = odf[groupord, 'COLOR']
atpar = cumsum(c(0, diff(-rev(as.numeric(odf$category)))) + 1)
cuts = atpar[c(0, diff(-rev(as.numeric(odf$category)))) == 1]
starts = c(1, cuts)
ends = c(cuts, max(atpar))

# -----------------
# Line definitions:
# -----------------
linedf =data.frame(y0=atpar, y1=atpar, x0=1, x1=atpar)
linedf = merge(linedf, data.frame(start=starts, end=ends))
linedf = linedf[linedf$x1 > linedf$start,]
sel = (linedf$x1 < linedf$end) & (linedf$x1 > linedf$start)
linedf$end[sel] =  linedf$x1[sel] + 1
linedf$end[nrow(linedf)] = linedf$end[nrow(linedf)] + 1
# Vertical lines:
ends2 = c(cuts, max(atpar) + 2) - 2
starts2 = c(2, cuts)
vlinedf =data.frame(x0=rev(atpar), x1=rev(atpar), y0=max(atpar), y1 = max(atpar) - atpar + 1)
vlinedf = merge(vlinedf, data.frame(start=starts2, ends2))
vlinedf = vlinedf[vlinedf$x1 < vlinedf$end,]
sel = (vlinedf$x1 >= vlinedf$start) & (vlinedf$x1 <= vlinedf$end)
vlinedf$start[sel] =  vlinedf$x1[sel] + 1
vlinedf$start[nrow(vlinedf)] = vlinedf$start[nrow(vlinedf)] + 1
# For asymmetrical:
if (use.asym){
    # Line definitions:
    ends2 = c(cuts, max(atpar) + 2)
    linedf =data.frame(y0=atpar, y1=atpar, x0=1, x1=atpar)
    linedf = merge(linedf, data.frame(start=starts, end=ends2))
    starts2 = c(2, cuts)
    vlinedf =data.frame(x0=rev(atpar), x1=rev(atpar), y0=max(atpar), y1 = max(atpar) - atpar + 1)
    vlinedf = merge(vlinedf, data.frame(start=starts, ends2 - 2))
}


# Add the top GWAS For each of the top 5 or so intersections (by pval + intersection size)
# Sort by shared pvalue of the two over total: 
pnodemat = all.regmat >= 2
pnodetismat = pnodemat %*% ntmat
pnodetismat = pnodetismat[apply(pnodetismat,1, sum) > 0,]
pnodetismat = sweep(pnodetismat, 1, apply(pnodetismat, 1, sum), '/')
range(pnodetismat)

# png(paste0(imgpref, 'node_tissue_jaccard_signif_ord.png'),res=450,units='in',width=12,height=12)
if (use.asym) {
    subdt = dt
    subdt$X = subdt$N2
    subdt$Y = subdt$N1
    xlim=c(0, max(atpar) + 1)
    ylim=c(0, max(atpar) + 1)
} else {
    subdt = dt[dt$N1 < dt$N2,]
    subdt$X = subdt$N1
    subdt$Y = subdt$N2
    xlim=c(0, max(atpar))
    ylim=c(1, max(atpar)+1)
}

png(paste0(imgpref, 'node_tissue_', metric, '_t', topn, '_signif_ord.png'),res=450,units='in',width=12,height=12.5)
par(mar=c(2.5, 7, 7, 0.5))
par(yaxs='i')
par(xaxs='i')
minor = 1
plot(atpar[subdt$X], atpar[subdt$Y],
     pch=22, cex=subdt$cex, bg=subdt$col, axes=F,
     xlab='', ylab='', col=subdt$border, xlim=xlim, ylim=ylim)
par(xpd=FALSE)
if (use.asym){
    points(rep(0,NT), atpar, pch=22, cex=4, col='white', bg=gcol)
    points(atpar, rep(max(atpar) + 1, NT), pch=22, cex=4, col='white', bg=gcol)
} else {
    points(rep(0,NT-1), atpar[-1], pch=22, cex=4, col='white', bg=gcol[-1])
    points(atpar[-NT], rep(max(atpar) + 1, NT-1), pch=22, cex=4, col='white', bg=gcol[-NT])
}
# Add the colors:
segments(x0=vlinedf$x0, x1=vlinedf$x1,
         y0=vlinedf$start, y1=vlinedf$end, lty='dashed', lwd=.5)
segments(x0=linedf$start, x1=linedf$end - 2,
         y0=linedf$y0, y1=linedf$y1, lty='dashed', lwd=.5)
par(new=TRUE)
# plot(as.numeric(subdt$G1), as.numeric(subdt$G2), 
plot(atpar[subdt$X], atpar[subdt$Y],
     pch=22, cex=subdt$cex, bg=subdt$col, axes=F,
     xlab='', ylab='', col=subdt$border, xlim=xlim, ylim=ylim)
text(atpar[subdt$X], atpar[subdt$Y], labels=subdt$label, cex=.6, col=txtcol)
loc = atpar
xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
yt = par()$usr[4] + 0.005*diff(par()$usr[3:4])
if (use.asym){
    text(y=loc, x=xt, labels=groupord, srt=0, adj=1, xpd=TRUE, cex=1)
    text(y=yt, x=loc, labels=groupord, srt=90, adj=0, xpd=TRUE, cex=1)
} else {
    text(y=loc[-1], x=xt, labels=groupord[-1], srt=0, adj=1, xpd=TRUE, cex=1)
    text(y=yt, x=loc[-NT], labels=groupord[-NT], srt=90, adj=0, xpd=TRUE, cex=1)
}
points(atpar, atpar, pch=22, cex=3, col='white', bg=gcol)
# points(atpar, atpar, pch=22, cex=3, bg='white', col=gcol, lwd=4)
par(xpd=TRUE)
legend('bottom', legend=c('NS','p < 0.1','p < 0.05','p < 0.01'), title='p-value',
       col='black', fill='white', pch=22, pt.cex=c(1,2,3,4), bty='n', border=NA, cex=1, horiz=T, inset=c(0,-.04))
draw(plegend, x = unit(10,'in'), y=unit(0.6,'in'), just = "top")
dev.off()


# --------------------------
# Plot without significance.
# --------------------------
subdt = dt[dt$N1 < dt$N2,]
subdt$X = subdt$N1
subdt$Y = subdt$N2
xlim=c(0, max(atpar))
ylim=c(1, max(atpar)+1)
subdt$cex = 4
subdt$border ='white'
# subdt$border =NA

png(paste0(imgpref, 'node_tissue_', metric, '_t', topn, '_ord.png'),res=450,units='in',width=12,height=12.5)
par(mar=c(2.5, 7, 7, 0.5))
par(yaxs='i')
par(xaxs='i')
minor = 1
plot(atpar[subdt$X], atpar[subdt$Y],
     pch=22, cex=subdt$cex, bg=subdt$col, axes=F,
     xlab='', ylab='', col=subdt$border, xlim=xlim, ylim=ylim)
par(xpd=FALSE)
points(rep(0,NT-1), atpar[-1], pch=22, cex=4, col=NA, bg=gcol[-1])
points(atpar[-NT], rep(max(atpar) + 1, NT-1), pch=22, cex=4, col=NA, bg=gcol[-NT])
# Add the colors:
par(new=TRUE)
plot(atpar[subdt$X], atpar[subdt$Y],
     pch=22, cex=subdt$cex, bg=subdt$col, axes=F,
     xlab='', ylab='', col=subdt$border, xlim=xlim, ylim=ylim)
text(atpar[subdt$X], atpar[subdt$Y], labels=subdt$nint, cex=.6, col=txtcol)
loc = atpar
xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
yt = par()$usr[4] + 0.005*diff(par()$usr[3:4])
text(y=loc[-1], x=xt, labels=groupord[-1], srt=0, adj=1, xpd=TRUE, cex=1)
text(y=yt, x=loc[-NT], labels=groupord[-NT], srt=90, adj=0, xpd=TRUE, cex=1)
points(atpar, atpar, pch=22, cex=3, col='white', bg=gcol)
par(xpd=TRUE)
draw(plegend, x = unit(10,'in'), y=unit(0.6,'in'), just = "top")
dev.off()


# -----------------------------------------------
# Add some driving GWAS for pairs of enrichments:
# -----------------------------------------------
tlist = rbind(data.frame(tissue='Liver', against=c('Adipose','Digestive','HSC & B-cell', 'Blood & T-cell', 'Pancreas')), 
              # data.frame(tissue='Endocrine', against=c('Adipose','Heart','Epithelial')),
              # data.frame(tissue='Muscle', against=c('Bone','Myosat','Adipose','Heart', 'Lung')),
              data.frame(tissue='Heart', against=c('Bone','Adipose','Muscle', 'Endocrine')),
              data.frame(tissue='Adipose', against=c('Endothelial', 'Stromal','Heart','Muscle')))
              # data.frame(tissue='PNS', against=c('Sm. Muscle', 'Reproductive','Other','Bone', 'Adipose','Heart','Muscle')))
tlist$toptext = ''
NF = nrow(tlist)
tlist$tissue = as.character(tlist$tissue)
tlist$against = as.character(tlist$against)

for (toptype in c('byprincipal','any','either')){
    topn=5
    princip = colnames(pnodetismat)[apply(pnodetismat, 1, which.max)]
    for (i in 1:NF){
        main = as.character(tlist$tissue[i])
        second = as.character(tlist$against[i])
        if (toptype == 'byprincipal'){
            ind = which((princip == main) & pnodetismat[,second] > 0)
        } else if (toptype == 'either'){
            ind1 = which((princip == main) & pnodetismat[,second] > 0)
            ind2 = which((princip == second) & pnodetismat[,main] > 0)
            ind = unique(c(ind1, ind2))
        } else {
            ind = which(pnodetismat[,main] > 0 & pnodetismat[,second] > 0)
        }
        if (length(ind) > 1){
            smat = pnodetismat[ind, c(main,second)]
            # Get top 5 traits by percent of enrichments:
            smarg = apply(smat,1, sum)
            nord = names(sort(smarg, decreasing=T))
            nord = head(nord,topn)
            tord = sub("[0-9]* - ","", nord)
            pcts = paste0(round(smarg[nord] * 100,1), '%')
            ttxt = paste0(paste0(tord, ' (', pcts, ')'), collapse=", ")
            tlist$toptext[i] = split.text(ttxt, 110)
        }
    }
    # Grid: 1st, 2nd, top 5 enrichments, using textwidth:
    tlist$nlines = sapply(tlist$toptext, function(x){length(strsplit(x, "\n")[[1]])})
    tlist$COL1 = odf[tlist$tissue, 'COLOR']
    tlist$COL2 = odf[tlist$against, 'COLOR']

    # Number of facets to explore:
    # pdf('~/test.pdf', width=3.5, height=5.5)
    pdf(paste0(imgpref, 'node_', toptype, '_t', topn, '_top_expanded.pdf'),width=4,height=4)
    # layout(matrix(1:(NF * 3), NF, 3, byrow=TRUE), heights=tlist$nlines, widths=c(1,1, 5))
    par(xaxs='i')
    par(yaxs='i')
    fancylayout=FALSE
    sw=0.60
    sp=0.02
    if (fancylayout){
        midNF=round(NF/2)
        omat = as.numeric(sapply(1:(midNF*2), function(x){c(1,2,3,3) + (x-1) * 3}))
        omat = matrix(omat, ncol=2, byrow=TRUE)
        omat2 = cbind(omat[1:(2*midNF), ], rep(max(omat) + 1, 2*midNF), omat[(2*midNF+1):(midNF*4),])
        # Add primary secondary 
        layout(omat2, heights=rep(c(.05,.2), midNF), widths=c(sw,sw, .05, sw,sw), TRUE)
    } else {
        layout(matrix(1:(NF * 3), NF, 3, byrow=TRUE), heights=rep(.25, NF), widths=c(1,1, 5), TRUE)
    }
    for (i in 1:NF){
        # First color
        par(mar=rep(sp, 4))
        image(as.matrix(1), col=tlist$COL1[i], axes=F)
        text(mean(par()$usr[1:2]), mean(par()$usr[3:4]),
             tlist$tissue[i], cex=.5,col='white')
        par(mar=rep(sp, 4))
        image(as.matrix(1), col=tlist$COL2[i], axes=F)
        text(mean(par()$usr[1:2]), mean(par()$usr[3:4]),
             tlist$against[i], cex=.5, col='white')
        # Text:
        par(mar=c(sp*2, sp,sp,sp))
        plot(0, xlim=c(0,1), ylim=c(0,1), type='n', axes=F,
             ylab='', xlab='')
        # box(col='grey', lwd=0.25)
        if (tlist$nlines[i] == 2){ tcex=.35 } else { tcex = .35}
        text(0.01,.5, tlist$toptext[i], cex=tcex, adj=c(0,.5), xpd=TRUE)
    }
    dev.off()
}




# ---------------------------------
# Exploration of tissue-specificity
# ---------------------------------
# Marginal of the number of tissues per GWAS:
nmarg = apply(nodetismat, 1, sum)
nmarg = nmarg[nmarg > 0]

mdf = data.frame(uid=names(nmarg), marg=nmarg)
# By number of SNPs:
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
mdf = merge(nsnpdf, mdf)

par(mar=c(3.5, 7, 7, 0.5))
par(yaxs='i')
par(xaxs='i')
hist(mdf$marg, col='grey60',border='white', breaks=34, xlim=c(0,34))

png(paste0(imgpref, 'marginal_nsnp_cutoffs.png'),res=450,units='in',width=6,height=12)
par(mar=c(4, 4, 4, 0.5))
par(yaxs='i')
par(xaxs='i')
layout(matrix(c(1:4), 4,1), TRUE)
hist(mdf$marg[mdf$pValue <= 5], col='grey60',border='white', breaks=34, xlim=c(0,34))
hist(mdf$marg[mdf$pValue > 5], col='grey60',border='white', breaks=34, xlim=c(0,34))
hist(mdf$marg[mdf$pValue > 25], col='grey60',border='white', breaks=34, xlim=c(0,34))
hist(mdf$marg[mdf$pValue > 100], col='grey60',border='white', breaks=34, xlim=c(0,34))
dev.off()


# -----------------
# Height vs. value:
# -----------------

NIND = 50000
NIND = 40000 # 40k includes crohns, etc.
# NIND = 100000
keptgw = gwssdf[gwssdf$sampsize > NIND,]
keptgw = keptgw[order(keptgw$pubDate, decreasing=T),]
traits = unique(keptgw$trait)
keptuids = sapply(traits, function(x){
                      keptgw$uid[keptgw$trait == x][1]})
print(paste("Kept", length(keptuids), "GWAS (newest with sample > 40,000)"))

# gwssdf[gwssdf$trait == "Crohn's disease",]

keptuids = keptuids[keptuids %in% rownames(all.regmat)]
kuid = names(which(apply(all.regmat, 1, max) > 3))
kuid = kuid[kuid %in% keptuids]
print(paste("Kept", length(kuid), "GWAS (with signif/run)"))

# Image of the kept GWAS, clustered by node occurrence:
redmat = all.regmat[kuid,]
redmat[redmat < 3] = 0


# Numbers of GWAS
redmat = sweep(redmat, 1, apply(redmat, 1, sum), '/')
ntdf = data.frame(redmat)
ntdf$uid = rownames(redmat)
ntdf = gather(ntdf, node, value,-uid)
ntdf = ntdf[ntdf$value > 0,]
ntdf$node = as.numeric(sub("X", "", ntdf$node))
ntdf = merge(ntdf, data.frame(node=paste0(1:NN),
                              height=nodedf$X2))
ntdf = merge(ntdf, nodetissue)

plot(ntdf$value, ntdf$height, col=ntdf$COLOR)

# library(ggpubr)
gcol = odf$COLOR
names(gcol) = odf$GROUP

ggplot(ntdf, aes(value, height, color=GROUP)) + 
    facet_wrap(~GROUP) +
    # geom_point() + 
    geom_jitter() + 
    geom_smooth() + 
    theme_pubr() + 
    scale_color_manual(values=gcol)
ggsave('~/test.png',dpi=250,width=18,height=12, units='in')


#


