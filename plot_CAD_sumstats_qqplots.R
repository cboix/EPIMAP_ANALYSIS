#!/usr/bin/R
# -----------------------------------------------------
# Plot CAD against other traits -- qqplots
# A. qq-plots:
# 1. Extract a set of loci for each enriched sets
# 2. Intersect with sumstats for each trait - pull out the 
# - snps from CAD per tissue
# - snps from CAD tissue-specific
# - all enh (+/-2.5kb) in node
# - all snps in the other trait
# - abline
# 
# B. CAD vs. TRAIT:
# - Per tissue, plot CAD vs. trait. separate +/- beta effects.
# -----------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
library(qvalue)
library(gridExtra)
library(grid)
library(clusterProfiler)  # For GO
library(org.Hs.eg.db)
library(GenomicRanges)
# ------------------------
# Arguments for GWAS data:
# ------------------------
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
use.adj = TRUE
use.strict = TRUE
# use.strict = FALSE
MINP=3

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict) }
source(paste0(bindir, 'load_statistics_gwastree_enrichments.R'))

# ------------
# Directories:
# ------------
examp = 'CAD'
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
caddir = paste0(img, 'CAD_example/')
epref = paste0(usetree, '_e', tol, '_')
cadpref = paste0(caddir, examp, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, regdir, perdir, caddir)
system(cmd)


# ----------------
# Trait arguments:
# ----------------
suid = "29212778 - Coronary artery disease" # Have: CAD_UKBIOBANK.gz
slist = c("30595370 - Systolic blood pressure",  # Have: systolic_UKB2_sumstats.txt.gz
          "29507422 - High density lipoprotein cholesterol levels", # Can't find summ stats.
          "30061737 - Atrial fibrillation", # Here:
          "26426971 - Waist-to-hip ratio adjusted for BMI", # Split into m/f
          "30595370 - Waist-hip ratio"
) 

# -------------------------------------
# Load all gwas enrichment information:
# -------------------------------------
snpfile = paste0(regdir, epref, apref, '_logreg_all_wsnp', suffix)
load(snpfile)
cad.lp = rmat[suid,]
nodetissue = nodetissue[order(nodetissue$node),]
# Load module assignments:
flatset = 'modules'
enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
load(enhsetfile)
enhsetdf = ldply(1:300, function(i){data.frame(cls=paste0('c', i-1), row=enhsets[[i]])})

MINP=3
sub.qdf = qdf[qdf$uid == suid, ]
enhmarg = read.delim(gzfile(enhmargfile))

# For enhancers in significant nodes, 
# find composition of each module:
cad.nodes = which(cad.lp > 3)
NPLOT=min(20, length(cad.nodes))  # TOP N
nodegroup = nodetissue$GROUP[cad.nodes]
nodecol = nodetissue$COLOR[cad.nodes]
rawnodelp = cad.lp[cad.nodes]
nodelp = round(rawnodelp,1)
nodenames = paste0(nodegroup, ' - ', nodelp)
nord = tail(order(cad.lp[cad.nodes]), NPLOT)

# Get the names underneath the node:
leafrep = sapply(cad.nodes[nord], function(x){
                     # blacklist = c('GLAND', 'TISSUE', 'CELL')
                     x = sub("[0-9]*_","",declist$dec[[x]])
                     x = unique(x)
                     nx = length(x)
                     if (nx > 1){
                         term_words <- strsplit(x, "[ _,.]");
                         tab_all <- sort(table(unlist(term_words)));
                         tab_all = tab_all[tab_all > 1]
                         # tab_all = tab_all[!(names(tab_all) %in% blacklist)]
                         x = paste0(tolower(names(sort(tab_all))), collapse=', ')
                     } else {
                         x = tolower(x)
                     }
                     x = capitalize(x)
                     return(x)})
lind = which(leafrep == '')
leafrep[lind] = nodegroup[nord][lind]
# Alternatively label all with group:
# leafrep = nodegroup[nord]


# SNP matrices:
out = sapply(cad.nodes, function(x){cdll$diff[[x]]})
omat = sapply(out, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
out.cons = sapply(1:NN, function(x){cdll$cons[[x]]})
omat.cons = sapply(out.cons, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
colnames(omat) = nodenames
tform = make.tform(sub.qdf$queryHits)
smat = t(tform) %*% omat
smat.cons = t(tform) %*% omat.cons
smat[smat > 0] = 1
smat.cons[smat.cons > 0] = 1
snps = unique(sub.qdf$queryHits)
rownames(smat) = paste0('s', snps)
rownames(smat.cons) = paste0('s', snps)
colnames(smat.cons) = paste0(nodetissue$node, '-',nodetissue$GROUP)

# Get snp specificties and locations:
fmat = smat %*% make.tform(nodegroup)
print(apply(fmat > 0, 2, sum))
# Specifically for top nodes:
nid = rev(nord)
head(smat[, nid])

# Reduce-snps:
smarg = apply(smat, 1, sum)
keep.const = FALSE
if (keep.const){
    ksnp = names(which(smarg > 0))
} else { 
    ksnp = names(which(smarg <= 65 * .5 & smarg > 0))
}

# Get GR for each loci:
snpsets = lapply(nid, function(x){ 
                     x = names(which(smat[ksnp, x] > 0)); 
                     as.numeric(sub("s", "", x))})
snpgrs = lapply(snpsets, function(x, tol=2500){
                    GRanges(paste0('chr', gwdf$chrom[x]), 
                            IRanges(gwdf$chromStart[x] - tol, gwdf$chromEnd[x] + tol)) })

# Get GR for each node as well:
enhgrs = lapply(cad.nodes[nid], function(x){
                   es = cdll$diff[[x]]
                   es = enhmap[es]
                   return(dmgr[es]) })

# Load in data:
mgdir = 'gwas_sumstats/munged/'
preflist= c('SBP', 'HDL', 'AFB', 'WAB', 'WHR', 'CAD')

# FI

# Clean workspace:
rm(cdll, enhdf, all.regmat, rmat, out.cons, regwide)
gc()

# FIX WAB:
mainfile = paste0(mgdir, 'WAB.munged.tsv.gz')
if (!file.exists(mainfile)){
    fulldf = NULL
    wabprefs = c('WHR.WG5', 'WHR.WL5','WHR.MG5','WHR.ML5')
    for (pref in wabprefs){
        print(paste0("Reading in ", pref))
        ssfile = paste0(mgdir, pref, '.munged.tsv.gz')
        ssdf = read.delim(gzfile(ssfile), header=F, sep="\t")
        names(ssdf) = c('chr','pos',paste0(pref, '.beta'),paste0(pref, '.p'))
        ssdf = ssdf[!is.na(ssdf$pos),]
        ssdf = ssdf[as.character(ssdf$chr) != 'chrNA',]
        if (is.null(fulldf)){
            fulldf = ssdf[, c(1:4)]
        } else {
            fulldf = merge(fulldf, ssdf[, c(1:4)])
        }
    }
    # Make the full table from the minimal p-vals:
    pcols = paste0(wabprefs, '.p')
    bcols = paste0(wabprefs, '.beta')
    wm = apply(fulldf[,pcols], 1, which.min)
    dfid = cbind(row=1:nrow(fulldf), col=wm)
    minp = fulldf[,pcols][dfid]
    minbeta = fulldf[,bcols][dfid]
    # Write table to main file:
    ssdf = data.frame(chr=fulldf$chr, pos=fulldf$pos, beta=minbeta, p=minp)
    write.table(ssdf, gzfile(mainfile), col.names=F, quote=F, sep='\t', row.names=F)
}






NTOP=20
snp.pm = matrix(NA, nrow=5, ncol=NTOP)
node.pm = matrix(NA, nrow=5, ncol=NTOP)
nodevenh.pm = matrix(NA, nrow=5, ncol=NTOP)
snpvenh.pm = matrix(NA, nrow=5, ncol=NTOP)
enh.pm = matrix(NA, nrow=5, ncol=NTOP)
rownames(snp.pm) = preflist[1:5]
rownames(node.pm) = preflist[1:5]
rownames(nodevenh.pm) = preflist[1:5]
rownames(snpvenh.pm) = preflist[1:5]
rownames(enh.pm) = preflist[1:5]
png(paste0(cadpref, 'qqplots_all_vCAD_t', NTOP, '.png'), res=450, units='in', width=10, height=5.5)
layout(matrix(1:(NTOP*5),5,NTOP, byrow=TRUE), widths=c(1.2, rep(1, NTOP-1)), heights=c(1.2,rep(1,4)), TRUE)
par(yaxs='i')
par(xaxs='i')
sp=0.1
preflist = preflist[1:5]
# Read all prefixes:
for (pref in preflist){
    # Read sumstats:
    print(paste0("Reading in ", pref))
    sbfile = paste0(mgdir, pref, '.munged.binned.tsv.gz')
    ssfile = paste0(mgdir, pref, '.munged.tsv.gz')
    if (!file.exists(sbfile)){
        ssdf = read.delim(gzfile(ssfile), header=F, sep="\t")
        names(ssdf) = c('chr','pos','beta','p')
        ssdf = ssdf[!is.na(ssdf$pos),]
        ssdf = ssdf[as.character(ssdf$chr) != 'chrNA',]
        ssdf$bin = round(ssdf$pos / 1e4)
        sbdf = aggregate(p ~ chr + bin, ssdf, min)
        sbdf = sbdf[order(sbdf$chr),]
        write.table(sbdf, gzfile(sbfile), col.names=T, quote=F, sep='\t', row.names=F)
        rm(ssdf)
    } else {
        sbdf = read.delim(gzfile(sbfile), header=T, sep='\t')
    }

    # Put midpt +/- 2500 as the loc, then we only intersect once with 5kb regions:
    ssgr = GRanges(as.character(sbdf$chr), IRanges(sbdf$bin * 1e4 + 2.5e3, sbdf$bin * 1e4 + 2.5e3 ))
    lp = -log10(sbdf$p)
    mlp = max(lp[!is.infinite(lp)])
    lp[is.infinite(lp)] = mlp
    # All enhancers:
    odf = data.frame(findOverlaps(dmgr, ssgr))
    x = unique(odf$subjectHits)
    allenh.intp = sort(lp[x], decreasing=T)

    # Intersect against each GR:
    print(paste0("Intersecting ", pref))
    snp.intp = sapply(snpgrs, function(gr){
               odf = data.frame(findOverlaps(gr, ssgr))
               x = unique(odf$subjectHits)
               x = -log10(sbdf$p[x])
               return(sort(x, decreasing=T)) })
    enh.intp = sapply(enhgrs, function(gr){
               odf = data.frame(findOverlaps(gr, ssgr))
               x = unique(odf$subjectHits)
               x = -log10(sbdf$p[x])
               return(sort(x, decreasing=T)) })
    NP = length(lp)
    hx = -log10((1:NP) / (NP + 1))
    lp = sort(lp, decreasing=T)

    # Calculate ranges:
    lx = c(0, hx[1] * 1.1)
    ly = c(0, lp[1] * 1.1)
    height=1 + .25 * (pref=='SBP')

    # png(paste0(cadpref, 'qqplots_',pref, '_vCAD_t', NTOP, '_vtotal.png'), res=450, units='in', width=10.2, height=height)
    # layout(matrix(1:(NTOP),1,NTOP, byrow=TRUE), widths=c(1.2, rep(1, NTOP-1)), heights=c(height), TRUE)
    # par(yaxs='i')
    # par(xaxs='i')
    # sp=0.1
    # print(paste0("Plotting ", pref))
    for (i in 1:NTOP){
        par(mar=c(sp,sp + (1-sp) * (i==1), sp + (1-sp) * (pref=='SBP'), sp))
        icol = nodecol[rev(nord)][i]
        x = snp.intp[[i]]
        ix = round(seq(1, NP, length.out=length(x)))
        plot(lp[ix], x, col=icol, pch=19, cex=.1,
             axes=F, xlim=ly, ylim=ly)
        ut = wilcox.test(x, lp, alternative='greater')
        snp.pm[pref,i] = -log10(ut$p.value)
        text(0.8 * ly[2], 0.3 * ly[2], paste0("-log10p = ", round(-log10(ut$p.value),1)), cex=.4, col=icol)
        # add enh
        x = enh.intp[[i]]
        ix = round(seq(1, NP, length.out=length(x)))
        points(lp[ix], x, col=tsp.col(icol), pch='.',
             xlim=ly, ylim=ly)
        ut = wilcox.test(x, lp, alternative='greater')
        node.pm[pref,i] = -log10(ut$p.value)
        text(0.8 * ly[2], 0.2 * ly[2], paste0("-log10p = ", round(-log10(ut$p.value),1)), cex=.4, col=tsp.col(icol))
        # add all ENHANCERS:
        x = allenh.intp
        ix = round(seq(1, NP, length.out=length(x)))
        points(lp[ix], x, col='grey', pch='.',
             xlim=ly, ylim=ly)
        ut = wilcox.test(x, lp, alternative='greater')
        enh.pm[pref,i] = -log10(ut$p.value)
        text(0.8 * ly[2], 0.1 * ly[2], paste0("-log10p = ", round(-log10(ut$p.value),1)), cex=.4, col='grey')
        abline(0,1, lwd=.25)
        box(lwd=0.5)
        # Test snp vs. enhancers:
        ut = wilcox.test(snp.intp[[i]], allenh.intp, alternative='greater')
        snpvenh.pm[pref,i] = -log10(ut$p.value)
        # Test node vs. enhancers:
        ut = wilcox.test(enh.intp[[i]], allenh.intp, alternative='greater')
        nodevenh.pm[pref,i] = -log10(ut$p.value)
        if (pref == 'SBP'){ mtext(rev(leafrep)[i], side=3, cex=.5, col=icol) } 
        if (i == 1){ mtext(pref, side=2, cex=.5) }
    }
    # dev.off()
}
dev.off()


# MAKE P-VALUE FIGURE:
icols = unique(nodecol[rev(nord)][1:NTOP])
faclabels = factor(nodecol[rev(nord)][1:NTOP], levels=icols)
ext.preflist = c("Systolic\nblood pressure",
          "HDL\ncholesterol",
          "Atrial\nfibrillation",
          "Waist-to-hip ratio\n(adjusted for BMI)",
          "Waist-hip ratio")


# pltchoice = 'snps'
for (pltchoice in c('snps', 'nodes','nodevenh', 'snpvenh')){
    if (pltchoice == 'snps'){
        mat = snp.pm[rev(preflist),]
        PCUT = 10
    } else if (pltchoice == 'nodes'){
        mat = node.pm[rev(preflist),]
        PCUT = 100
    } else if (pltchoice == 'nodevenh'){
        mat = nodevenh.pm[rev(preflist),]
        PCUT = 50
    } else if (pltchoice == 'snpvenh'){
        mat = snpvenh.pm[rev(preflist),]
        PCUT = 10
    }
    scale = 2
    # png(paste0(cadpref, 'qqplots_all_pvals_', pltchoice,'_t', NTOP,'.png'), res=450, units='in', width=(NTOP) / scale, height=(2 + 5) / scale)
    pdf(paste0(cadpref, 'qqplots_all_pvals_', pltchoice,'_t', NTOP,'.pdf'), width=(NTOP) / scale, height=(2 + 5) / scale)
    rsp=8
    sp=.1
    layout(matrix(1:4, 2,2,byrow=TRUE), widths=c(5,NTOP), heights=c(5,5), TRUE)
    par(mar=c(sp, rsp, sp, sp))
    image(as.matrix(1), col='white', axes=F, ylab='', xlab='')
    text(x=.5, y=par()$usr[3] + .02 * diff(par()$usr[3:4]),
         labels='All SNPs in enhancers\nvs. all SNPs', col='black',
         cex=.85, xpd=TRUE, adj=c(0, 0), srt=90)
    # Labels:
    par(mar=c(sp, sp, 1.25 * rsp, sp))
    image(as.matrix(as.numeric(faclabels)), col=icols, axes=F)
    text(x=seq(0,1,length.out=NTOP), y=par()$usr[4] + .2 * diff(par()$usr[3:4]), 
         labels=rev(leafrep)[1:NTOP], col=nodecol[rev(nord)][1:NTOP],
         cex=.85, xpd=TRUE, adj=0, srt=90)
    # vs. enh.pm
    par(mar=c(sp, rsp, sp, sp))
    emat = enh.pm[,1]
    emat[emat > PCUT] = PCUT
    image(t(emat), col=col1, zlim=c(0,PCUT), axes=F)
    text(x=par()$usr[1] - .1 * diff(par()$usr[1:2]), 
         y=seq(0,1,length.out=length(preflist)), labels=rev(ext.preflist), 
         cex=1, xpd=TRUE, adj=1)
    text(x=mean(par()$usr[1:2]), y=seq(0,1,length.out=length(preflist)),
         labels=round(enh.pm, 1), cex=1, xpd=TRUE, col=ifelse(enh.pm > 0.25 * PCUT, 'grey85', 'grey15'))
    # vs. snp.pm
    par(mar=rep(sp, 4))
    mdf = data.frame(mat)
    mdf$pref = factor(rownames(mdf), levels=rownames(mdf))
    mlong = gather(mdf, col, value, -pref)
    mlong$col = factor(mlong$col, levels=unique(mlong$col))
    mat[mat > PCUT] = PCUT
    image(t(mat), col=col1, zlim=c(0,PCUT), axes=F, useRaster=TRUE)
    xat = seq(0,1,length.out=ncol(mat))
    yat = seq(0,1,length.out=nrow(mat))
    mlong = filter(mlong, value > 3)
    text(x=xat[as.numeric(mlong$col)], y=yat[as.numeric(mlong$pref)], 
         labels=round(mlong$value,1), col=ifelse(mlong$value >= .5 * PCUT, 'grey85', 'black'))
    dev.off()
}



# -----------------------------------------------------------
# PART B: Plot the CAD vs. trait SNPs - load the wCAD version
# -----------------------------------------------------------

NTOP=10
png(paste0(cadpref, 'beta_comp_all_vCAD_t', NTOP, '.png'), res=450, units='in', width=10, height=5.5)
layout(matrix(1:(NTOP*5),5,NTOP, byrow=TRUE), widths=c(1.2, rep(1, NTOP-1)), heights=c(1.2,rep(1,4)), TRUE)
par(yaxs='i')
par(xaxs='i')
sp=0.1
preflist = preflist[1:5]
# Read all prefixes:
for (pref in preflist){
    # Read sumstats:
    print(paste0("Reading in ", pref))
    ssdf = read.delim(gzfile(paste0(mgdir, pref, '.wCAD.tsv.gz')), header=F, )
    names(ssdf) = c('chr','pos','beta','p', 'cad.beta','cad.p')
    ssdf = ssdf[!is.na(ssdf$pos),]
    ssdf = ssdf[as.character(ssdf$chr) != 'chrNA',]
    ssgr = GRanges(as.character(ssdf$chr), IRanges(ssdf$pos, ssdf$pos))

    lp = -log10(ssdf$p)
    mlp = max(lp[!is.infinite(lp)])
    lp[is.infinite(lp)] = mlp
    cad.lp = -log10(ssdf$cad.p)
    cad.mlp = max(lp[!is.infinite(cad.lp)])
    cad.lp[is.infinite(cad.lp)] = cad.mlp

    col_fun = function(x, pal=colrb, range=c(0,100)){
        bin <- cut(x, seq(range[1], range[2], length.out=length(pal)), include.lowest=T) 
        pal[bin] }

    srx = c(-max(abs(ssdf$cad.beta)), max(abs(ssdf$cad.beta)))
    rx = c(-max(abs(ssdf$cad.beta)), max(abs(ssdf$cad.beta)))
    ry = c(-max(abs(ssdf$beta)), max(abs(ssdf$beta)))

    height=1 + .25 * (pref=='SBP')
    # png('~/test.png', res=450, units='in', width=10.2, height=height)
    # layout(matrix(1:(NTOP),1,NTOP, byrow=TRUE), widths=c(1.2, rep(1, NTOP-1)), heights=c(height), TRUE)
    # par(yaxs='i')
    # par(xaxs='i')
    # sp=0.1
    for (i in 1:NTOP){
        par(mar=c(sp,sp + (1.1-sp) * (i==1), sp + (1.1-sp) * (pref=='SBP'), sp))
        icol = nodecol[rev(nord)][i]
        x = enh.intx[[i]]
        # x = snp.intx[[i]]
        # Specific enhancer groups:
        x = x[which(-log10(ssdf$cad.p[x]) > 4 | -log10(ssdf$p[x]) > 4)]
        # ssdf$cad.beta[x]
        model = lm(beta ~ cad.beta, ssdf[x, ])
        ar2 = summary(model)$adj
        cc = coefficients(model)
        plot(ssdf$cad.beta[x], ssdf$beta[x], 
             col=tsp.col(icol), pch=19, cex=.1, axes=F, 
             xlim=rx, ylim=ry)
        text(.9 * rx[1], .9 * ry[1], paste0('adj.rsq = ', round(ar2, 4)), col=icol, cex=.5, adj=0)
        abline(cc[1], cc[2], lwd=.5, col='red')
        abline(h=0, lwd=.25)
        abline(v=0, lwd=.25)
        box(lwd=0.5)
        if (pref == 'SBP'){ mtext(rev(leafrep)[i], side=3, cex=.5, col=icol) } 
        if (i == 1){ mtext(pref, side=2, cex=.5) }
    }
    # dev.off()
}
dev.off()


    # Intersect against each GR:
    print(paste0("Intersecting ", pref))
    snp.intx = sapply(snpgrs, function(gr){
               odf = data.frame(findOverlaps(gr, ssgr))
               x = unique(odf$subjectHits)
               return(x) })

    enh.intx = sapply(enhgrs, function(gr){
               odf = data.frame(findOverlaps(gr, ssgr))
               x = unique(odf$subjectHits)
               return(x) })

    # Calculate ranges:
    lx = c(0, max(cad.lp) * 1.1)
    ly = c(0, max(lp) * 1.1)
    # range(ssdf$cad.beta)
    # range(ssdf$beta)

    png('~/test.png', res=450, units='in', width=10.2, height=height)
    height=1 + .25 * (pref=='SBP')
    NTOP = 10
    layout(matrix(1:(NTOP),1,NTOP, byrow=TRUE), widths=c(1.2, rep(1, NTOP-1)), heights=c(height), TRUE)
    par(yaxs='i')
    par(xaxs='i')
    sp=0.1
    for (i in 1:NTOP){
        par(mar=c(sp,sp + (1-sp) * (i==1), sp + (1-sp) * (pref=='SBP'), sp))
        icol = nodecol[rev(nord)][i]
        x = snp.intx[[i]]
        plot(-log10(ssdf$cad.p[x]), -log10(ssdf$p[x]),
             col=tsp.col(icol), pch=19, cex=.1, axes=F, 
             xlim=lx, ylim=ly)
        box(lwd=0.5)
        if (pref == 'SBP'){ mtext(rev(leafrep)[i], side=3, cex=.5, col=icol) } 
        if (i == 1){ mtext(pref, side=2, cex=.5) }
    }
    dev.off()


    png('~/test.png', res=450, units='in', width=10.2, height=height)
    height=1 + .25 * (pref=='SBP')
    NTOP = 10
    layout(matrix(1:(NTOP),1,NTOP, byrow=TRUE), widths=c(1.2, rep(1, NTOP-1)), heights=c(height), TRUE)
    par(yaxs='i')
    par(xaxs='i')
    sp=0.1
    for (i in 1:NTOP){
        par(mar=c(sp,sp + (1-sp) * (i==1), sp + (1-sp) * (pref=='SBP'), sp))
        icol = nodecol[rev(nord)][i]
        x = enh.intx[[i]]
        # Specific enhancer groups:
        plot(-log10(ssdf$cad.p[x]), -log10(ssdf$p[x]),
             col=tsp.col(icol), pch='.', cex=.1, axes=F, 
             xlim=lx, ylim=ly)
        box(lwd=0.5)
        if (pref == 'SBP'){ mtext(rev(leafrep)[i], side=3, cex=.5, col=icol) } 
        if (i == 1){ mtext(pref, side=2, cex=.5) }
    }
    dev.off()


}

