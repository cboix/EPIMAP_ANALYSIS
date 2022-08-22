#!/usr/bin/R
# ------------------------------------
# Compare DHSs to screen dELs and pELs
# in their overlap with GWAS hits:
# TODO: Also do an epigenome enrichment?
# ------------------------------------
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
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'auxiliary_gwastree_functions.R'))

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

# ------------------------------------------------------------------------
# Load in all of the GWAS analysis data for simple overlaps + enrichments:
# ------------------------------------------------------------------------
# Files for GWAS analysis:
ddir = 'DHS_Index_WM201902/'
dpref = 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19'
dmlfile = paste0(ddir, dpref, '.core.srt.txt')
dmlnamfile = paste0(ddir, dpref, '_r200_e0_names.core.srt.tsv')
dmlrdafile = paste0(ddir, dpref, '_enh.core.srt.Rda')
gwcatfile = 'gwascatalog_may03_2019_noquotes.txt'
gwrdafile = sub(".txt", ".Rda", gwcatfile)
gwrsidfile = sub(".txt", "_rsid.txt", gwcatfile)
gwintrdafile = sub(".txt", "_intersections.Rda", gwcatfile)
# Load  data:
load(dmlrdafile)
load(gwintrdafile)
load(gwrdafile)
tol = 2500
dmgr = GRanges(enhdf$chr, IRanges(enhdf$start - tol, enhdf$end + tol), name=enhdf$name)

# All DHS:
uqdf = read.delim(dmlfile, header=F)
names(uqdf) = c('chr','start','end', 'loc')
uqgr = GRanges(uqdf$chr, IRanges(uqdf$start - tol, uqdf$end + tol), name=uqdf$loc)

# Prefs:
today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "validation/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'screen_comparison_')

# --------------------------
# Questions on GWAS overlap:
# --------------------------
# - How much overall overlap with lead SNPs in each set
# - Both 2.5kb and right on top
# - Plot each distribution of nearest enh to each snp
# - How often a closer epimap enh. to the lead SNP

# -----------------------------------
# Read in the SCREEN comparison sets:
# -----------------------------------
denh = scan('epimap_screen_dELS.txt','c')
nenh = scan('epimap_screen_non_dELS.txt','c')
eenh = scan('epimap_non_screen_elements.txt','c')

# General overlap:
qdf$name = dmgr$name[qdf$subjectHits]
qdf$set = 'None'
qdf$set[qdf$name %in% denh] = 'dELS'
qdf$set[qdf$name %in% eenh] = 'Epimap-only'
qdf$set[qdf$name %in% nenh] = 'non-dELS'
qdf = qdf[qdf$set != 'None',]

ndf = data.frame(set = c('dELS','Epimap-only','non-dELS'),
                 num=c(length(denh), length(eenh), length(nenh)),
                 novl=c(length(unique(qdf$name[qdf$set == 'dELS'])),
                        length(unique(qdf$name[qdf$set == 'Epimap-only'])),
                        length(unique(qdf$name[qdf$set == 'non-dELS']))))


dhs.qdf = suppressWarnings(data.frame(findOverlaps(gwgr, uqgr)))
dhs.qdf$name = uqgr$name[dhs.qdf$subjectHits]

# Random DHSs from the full set:
set.seed(1)
NPERM=100
rnovl = c()
rnsnp = c()
for (i in 1:NPERM){
    print(i)
    rnum = sort(sample(1:length(uqgr), length(eenh), replace=FALSE))
    renh = uqgr$name[rnum]
    rnovl = c(rnovl, length(unique(dhs.qdf$name[dhs.qdf$name %in% renh])))
    rnsnp = c(rnsnp, length(unique(dhs.qdf$queryHits[dhs.qdf$name %in% renh])))
}

# Plot random samples, next to others:
rswe = rnsnp / length(gwgr)
rfws = rnovl / length(eenh)

# Random DHSs that are non-screen:
set.seed(1)
NPERM=100
nnovl = c()
nnsnp = c()
non.screen.dhs = uqgr$name[!(uqgr$name %in% c(nenh, denh))]
for (i in 1:NPERM){
    print(i)
    nsenh = sample(non.screen.dhs, length(eenh), replace=FALSE)
    nnovl = c(nnovl, length(unique(dhs.qdf$name[dhs.qdf$name %in% nsenh])))
    nnsnp = c(nnsnp, length(unique(dhs.qdf$queryHits[dhs.qdf$name %in% nsenh])))
}

# Plot random samples, next to others:
nsswe = nnsnp / length(gwgr)
nsfws = nnovl / length(eenh)

# All 
neswe = length(unique(dhs.qdf$queryHits[dhs.qdf$name %in% non.screen.dhs])) / length(gwgr)
nefws = length(unique(dhs.qdf$name[dhs.qdf$name %in% non.screen.dhs])) / length(non.screen.dhs)
print(c(nefws, neswe))

non.enh.dhs = uqgr$name[!(uqgr$name %in% c(nenh, denh, eenh))]
neswe = length(unique(dhs.qdf$queryHits[dhs.qdf$name %in% non.enh.dhs])) / length(gwgr)
nefws = length(unique(dhs.qdf$name[dhs.qdf$name %in% non.enh.dhs])) / length(non.enh.dhs)
print(c(nefws, neswe))

# Overall:
agg.rename(queryHits ~ set, qdf, length, 'nsnp')
# Unique:
tdf = merge(agg.rename(queryHits ~ set, qdf, function(x){length(unique(x))},'nsnp'), ndf)
tdf = rbind(tdf, data.frame(set='All',nsnp=length(unique(qdf$queryHits)), num=length(dmgr), novl=length(unique(qdf$name))))

# % SNPs covered by enhancers and % Enhancers with SNP:
tdf$snp.w.enh = tdf$nsnp / length(gwgr)
tdf$frac.with.snp = tdf$novl / tdf$num
#           set  nsnp     num   novl snp.w.enh frac.with.snp
# 1        dELS 64099  672226 113040 0.5639787     0.1681577  # dELS in epimap enhancers
# 2 Epimap-only 87297 1286740 206039 0.7680876     0.1601248
# 3    non-dELS 18100  110654  20054 0.1592539     0.1812316
# 4         All 93047 2069090 339133 0.8186793     0.1639044  # All epimap - 81% of SNPs covered

set.cols = c('All'='grey75', 'non-dELS' = '#ffa700ff',
            'dELS' = '#ffcd00ff', 'Epimap-only'=rgb(255,195,77, max=255))

library(RColorBrewer)
cs = colorRampPalette(c('white','orange'))(10)

set.cols = c('All'='orange', 'non-dELS' = cs[5],
            'dELS' = cs[5], 'Epimap-only'=cs[5])

# Plot these as barplots as well:
gp = ggplot(tdf, aes(set, snp.w.enh, fill=set)) + 
    geom_bar(stat='identity') + 
    labs(y='Percent of GWAS lead SNPs within 2.5kb of enh. centers', x='Set (out of 2.1M Epimap enhancers)') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0), lim=c(0,1)) + 
    scale_fill_manual(values=set.cols) + 
    coord_flip() + 
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'gwas_snp_recovery_pct.png'), gp, dpi=450, units='in', width=5, height=3)
ggsave(paste0(imgpref, 'gwas_snp_recovery_pct.pdf'), gp, width=5, height=3)

# Plot these as barplots as well:
gp = ggplot(tdf, aes(set, frac.with.snp, fill=set)) + 
    geom_bar(stat='identity') + 
    labs(y='Percent of enhancers with GWAS SNP within 2.5kb', x='Set (out of 2.1M Epimap enhancers)') + 
    scale_y_continuous(labels=scales::percent, expand=c(0,0), lim=c(0,.25)) + 
    scale_fill_manual(values=set.cols) + 
    coord_flip() + 
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'gwas_frac_w_snp_pct.png'), gp, dpi=450, units='in', width=5, height=3)
ggsave(paste0(imgpref, 'gwas_frac_w_snp_pct.pdf'), gp, width=5, height=3)


# -------------------------------
# Nearest enhancers to each SNP: 
# -------------------------------
dind = dmgr$name %in% denh
eind = dmgr$name %in% eenh
nind = dmgr$name %in% nenh

enhdf$mid = (enhdf$start + enhdf$end) / 2
dmgr2 = GRanges(enhdf$chr, IRanges(enhdf$mid, enhdf$mid), name=enhdf$name)
adist = as.data.frame(distanceToNearest(gwgr,dmgr2))$distance
ddist = as.data.frame(distanceToNearest(gwgr,dmgr2[dind]))$distance
edist = as.data.frame(distanceToNearest(gwgr,dmgr2[eind]))$distance
ndist = as.data.frame(distanceToNearest(gwgr,dmgr2[nind]))$distance
nrind = nearest(gwgr,dmgr2)
nearchunk = dmgr$name[nrind]
distdf = rbind(data.frame(set='All',dist=adist),
               data.frame(set='dELS',dist=ddist),
               data.frame(set='Epimap-only',dist=edist),
               data.frame(set='non-dELS',dist=ndist))

mdistdf = aggregate(dist ~ set, distdf, mean)
sdistdf = aggregate(dist ~ set, distdf, function(x){sum(x<=100) / length(ddist)})
s2distdf = aggregate(dist ~ set, distdf, function(x){sum(x<=2500) / length(ddist)})

gp = ggplot(distdf, aes(set, dist, fill=set)) + 
    geom_violin(alpha=1) + 
    labs(y='Distance (bp) of lead SNPs to nearest enh. center', x='Set (out of 2.1M Epimap enhancers)') + 
    geom_point(data=mdistdf, aes(set, dist), pch='+', cex=4) + 
    geom_text(data=mdistdf, aes(set, dist, label=round(dist,0)), pch='+', cex=4) + 
    geom_text(data=sdistdf, aes(set, 100, label=paste0(round(100 * dist,0),'%')), pch='+', cex=4) + 
    scale_y_log10(labels=scales::comma, expand=c(0,0)) + 
    coord_flip() + 
    scale_fill_manual(values=set.cols) + 
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'gwas_dist_nearest.png'), gp, dpi=450, units='in', width=5, height=3)
ggsave(paste0(imgpref, 'gwas_dist_nearest.pdf'), gp, width=5, height=3)

# Finally, how often is epimap, dELS, etc. the closest enhancer:
ncdf = data.frame(set = c('dELS','Epimap-only','non-dELS'),
                  nclose=c(sum(nearchunk %in% denh), sum(nearchunk %in% eenh), 
                           sum(nearchunk %in% nenh)))
ncdf = merge(tdf, ncdf)
ncdf$close.ratio = ncdf$nclose / ncdf$num

gp = ggplot(ncdf, aes(x=factor(1),y=nclose, fill=set)) + 
    geom_bar(position='stack',stat='identity') + 
    scale_y_continuous(labels=scales::comma, expand=c(0,0)) + 
    scale_fill_manual(values=c('grey55','grey70','grey85')) +
    labs(y='Number of GWAS lead SNPs with closest enh. in set', x='Set (out of 2.1M Epimap enhancers)') + 
    theme_pubr() + theme(legend.position='none')
ggsave(paste0(imgpref, 'gwas_top_enh_cts.png'), gp, dpi=450, units='in', width=2, height=3)
ggsave(paste0(imgpref, 'gwas_top_enh_cts.pdf'), gp, width=2, height=3)


# ------------------------------
# Questions on GWAS enrichments:
# ------------------------------
# - Percentage of enrichments lost when using only dELS matched to epimap
# - Correlation between results of dELS and epimap non-dELS

tol = 2500  # Plus/minus distance - window for enhancer overlaps
# Minimal hypergeometric test (for this type of data)
get_pval_hg_fast = function(suid, dflist, lenlist, 
                            qallnumdf, NF=NF, cutp=CUTP, olen=NENH){
    # Measure numbers of consensus SNPs:
    isubdf = dflist[dflist$uid == suid,]
    isnp = rep(0, NF)
    isnp[isubdf$node] = isubdf$nsnp
    # Measure numbers of enhancer snps:
    osnp = qallnumdf$queryHits[qallnumdf$uid == suid]
    if (is.null(olen)){ olen = length(enhind) }
    # Hyper-geometric testing:
    pvdf = cbind(q=isnp, draw=lenlist, m=osnp, N=olen)
    pout <- apply(pvdf, 1, run.hyper)
    lpout = -log10(pout)
    if (sum(is.infinite(lpout)) > 0){
        pvdf = as.data.frame(pvdf)
        pvdf$lp = lpout
        pvdf = pvdf[order(pvdf$lp, decreasing=T),]
        print(head(pvdf))
    }
    return(-log10(pout))
}

load(gwrdafile)

NIND = 10000
NSNP = 10
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
gwssdf = aggregate(sampsize ~ uid, gwssdf, max)
gwssdf = merge(gwssdf, nsnpdf)
keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]
kept.uids =  sort(unique(as.character(keptgw$uid)))
NUID = length(kept.uids)
print(paste("Number of kept uids:", NUID))
gwdf = gwdf[gwdf$uid %in% kept.uids,]
gwgr = GRanges(paste0('chr', gwdf$chrom), IRanges(gwdf$chromStart, gwdf$chromEnd))

TOTSNP = nrow(gwdf)
NCAT = 1000
print(paste0("Pruned to: ", TOTSNP, " snps in ", NUID, " GWAS."))

# ---------------------------------------------
# Load in the query df + process and intersect:
# TODO: Allow load in of index separately;
# ---------------------------------------------
fdir = 'flat_hg_data_112420/'
bedfile = paste0(fdir, 'epi_observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_full.bed.gz')
rdafile = sub("bed.gz", "Rda", bedfile)
df = read.delim(gzfile(bedfile), header=F, stringsAsFactors=F)
names(df) = c('loc','cls')
print(dim(df))
gc()

# Unique locations (should be cleaned, new enhancer list):
uqdf = read.delim(dmlfile, header=F)
names(uqdf) = c('chr','start','end', 'loc')
keep.uq = which(uqdf$loc %in% df$loc)
uqdf = uqdf[keep.uq,]
NENH = nrow(uqdf)
uqdf$id = 1:NENH

indmap = uqdf$id
names(indmap) = as.character(uqdf$loc)
df$id = indmap[df$loc]
dmgr = GRanges(uqdf$chr, IRanges(uqdf$start - tol, uqdf$end + tol))

# Each of the different sets:
dnum = uqdf$id[uqdf$loc %in% denh]
enum = uqdf$id[uqdf$loc %in% eenh]
nnum = uqdf$id[uqdf$loc %in% nenh]
# NOTE ONLY RANDOM ENH: SHOULD TEST ALL DHSs:
set.seed(2)
rnum = sort(sample(1:length(dmgr), length(enum), replace=FALSE))
renh = uqdf$loc[rnum]

rm(uqdf)
gc()


# Turn into list (enhsets) - takes a while for epi:
if (!file.exists(rdafile)){
    clslist = sort(unique(df$cls))
    enhsets = lapply(clslist, function(x){df$id[df$cls == x]})
    NF = length(enhsets)
    lenlist = sapply(enhsets, length)
    save(enhsets, lenlist, file=rdafile)
    gc()
} else { 
    load(rdafile) 
}

# Calculate intersections with GWAS:
qdf = suppressWarnings(data.frame(findOverlaps(gwgr, dmgr)))
qdf$uid = gwdf$uid[qdf$queryHits]
qdf$uid = as.character(qdf$uid)
qdf$uid = factor(qdf$uid, levels=kept.uids)
qdf$fuid = as.numeric(qdf$uid)

qdf$set = 'Epimap-only'
qdf$set[qdf$subjectHits %in% nnum] = 'non-dELS'
qdf$set[qdf$subjectHits %in% dnum] = 'dELS'
sqdf = qdf[qdf$subjectHits %in% rnum,]
sqdf$set = 'random'
# Aggregate by set:
set.qallnumdf = aggregate(queryHits ~ uid + set, qdf, function(x){length(x)})
set.qallnumdf = rbind(set.qallnumdf, aggregate(queryHits ~ uid + set, sqdf, function(x){length(x)}))
qallnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(x)})
qallnumdf$set = 'All'
set.qallnumdf = rbind(set.qallnumdf, qallnumdf)
sqwide = spread(set.qallnumdf, set, queryHits, fill=0)

# Evaluate the intersections (pre-hypergeom test)
tmpdf = df[,c('id','cls')]
rm(df) # Don't need full df anymore
gc() 
names(tmpdf) = c('subjectHits', 'node')
tmpdf$node = tmpdf$node + 1 # From 0-indexed sets --> 1-indexed
tmp.mat = NULL
tmp.mat = Matrix(0, nrow=max(tmpdf$subjectHits), ncol=max(tmpdf$node))
tmp.mat[as.matrix(tmpdf)] = 1
# Make the dflist dataframe:
if (is.null(tmp.mat)){
    enhtmp = merge(qdf, tmpdf)
    dflist = aggregate(queryHits ~ uid + node, enhtmp, length)
    names(dflist)[3] = 'nsnp'
    dflist = dflist[,c('uid','nsnp','node')]
} else {
    # Faster merge using the precomputed transformation matrix
    qmat = Matrix(0, nrow=length(kept.uids), ncol=max(tmpdf$subjectHits))
    qmat[as.matrix(qdf[,c('fuid','subjectHits')])] = 1
    dim(unique(qdf))
    gc()

    mmat = qmat %*% tmp.mat
    mmat = data.frame(as.matrix(mmat))
    gc()
    dmmat = qmat[,dnum] %*% tmp.mat[dnum,]
    dmmat = data.frame(as.matrix(dmmat))
    gc()
    emmat = qmat[,enum] %*% tmp.mat[enum,]
    emmat = data.frame(as.matrix(emmat))
    gc()
    nmmat = qmat[,nnum] %*% tmp.mat[nnum,]
    nmmat = data.frame(as.matrix(nmmat))
    gc()
    rmmat = qmat[,rnum] %*% tmp.mat[rnum,]
    rmmat = data.frame(as.matrix(rmmat))
    gc()
    mmat$uid = kept.uids
    dmmat$uid = kept.uids
    emmat$uid = kept.uids
    nmmat$uid = kept.uids
    rmmat$uid = kept.uids

    llist = list()
    llist[['All']] = colSums(tmp.mat)
    llist[['dELS']] = colSums(tmp.mat[dnum,])
    llist[['Epimap-only']] = colSums(tmp.mat[enum,])
    llist[['non-dELS']] = colSums(tmp.mat[nnum,])
    llist[['random']] = colSums(tmp.mat[rnum,])

    nlist = c('All' = nrow(tmp.mat),
              'dELS' = length(dnum),
              'Epimap-only' = length(enum),
              'non-dELS' = length(nnum),
              'random' = length(rnum))

    # mmat = qmat %*% tmp.mat 
    # mmat = data.frame(as.matrix(mmat))
    # mmat$uid = kept.uids
    # dflist = gather(mmat, node, nsnp, -uid)
    mdflist = gather(mmat, node, nsnp, -uid)
    mdflist$set = 'All'
    ddflist = gather(dmmat, node, nsnp, -uid)
    ddflist$set = 'dELS'
    edflist = gather(emmat, node, nsnp, -uid)
    edflist$set = 'Epimap-only'
    ndflist = gather(nmmat, node, nsnp, -uid)
    ndflist$set = 'non-dELS'
    rdflist = gather(rmmat, node, nsnp, -uid)
    rdflist$set = 'random'
    dflist = rbind(mdflist, ddflist, edflist, ndflist, rdflist)

    dflist = dflist[dflist$nsnp > 0,]
    dflist$node = as.numeric(sub("^X","", dflist$node))
}
gc()

# --------------------------------------
# 1. Calculate the hyper-geometric test:
# --------------------------------------
# Main:
sets = c('All','dELS','Epimap-only','non-dELS', 'random')
ralldf = c()
for (set in sets){
    pmat.file = paste0(fdir, '/hg_raw_enr.',set,'.Rda')
    if (!file.exists(pmat.file)){
        print(paste("[STATUS] Calculating enrichments for", set))
        NF = length(enhsets)
        raw.pmat = matrix(NA, nrow=NUID, ncol=NF, dimnames=list(kept.uids,NULL))
        sdflist = dflist[dflist$set == set,]
            sqallnumdf = set.qallnumdf[set.qallnumdf$set == set,]
        for (i in 1:NUID){
            suid = kept.uids[i]
            rawlp = try(get_pval_hg_fast(suid, dflist, llist[[set]], sqallnumdf, NF=NF, cutp=CUTP, olen=nlist[[set]]))
            raw.pmat[suid,] = rawlp
        }
        rdf = data.frame(raw.pmat)
        rdf$uid = rownames(rdf)
        rdf = gather(rdf, node, rawlp, -uid)
        rdf$set = set
        save(raw.pmat, rdf, file=pmat.file)
    } else {
        load(pmat.file)
    }
    ralldf = rbind(ralldf, rdf)
    print(paste0("Quantiles for enrichment values in ",set, ":"))
    print(quantile(rdf$rawlp, c(.98, .99, .999, .9999), na.rm=T))
}

raw.pmat[raw.pmat > 20] = 20
raw.pmat[is.na(raw.pmat)] = 0
rmat = reord(raw.pmat)
kuid = which(apply(rmat,2, max) > 12)
image(rmat[,kuid])

raw.pmat[raw.pmat > 20] = 20
raw.pmat[is.na(raw.pmat)] = 0
rmat = reord(raw.pmat)
kuid = which(apply(rmat,2, max) > 8)
image(rmat[,kuid])

# Margins:
mdf = aggregate(node ~ subjectHits, tmpdf, length)
mdf$set = 'Epimap-only'
mdf$set[mdf$subjectHits %in% nnum] = 'non-dELS'
mdf$set[mdf$subjectHits %in% dnum] = 'dELS'
gc()

ggplot(mdf, aes(set, node, fill=set)) + 
    geom_violin() + 
    scale_y_log10() + 
    theme_pubr()

aggregate(node ~ set, mdf[mdf$node > 200,], length)


rwide = spread(ralldf, set, rawlp, fill=0)
print(colnames(rwide))
colnames(rwide) = c('uid','node','All', 'dELS','EpiMap','non_dELS','random')
# rwide$status

gp = ggplot(rwide, aes_string('dELS', 'EpiMap')) + 
    geom_point(alpha=.5) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    geom_abline() + 
    theme_pubr()
ggsave(paste0(imgpref, 'gwas_enr_scatter.png'), gp, dpi=450, units='in', width=4, height=4)
ggsave(paste0(imgpref, 'gwas_enr_scatter.pdf'), gp, width=4, height=4)


gp = ggplot(rwide, aes_string('All', 'dELS')) + 
    geom_point(alpha=.5) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr()
ggsave(paste0(imgpref, 'gwas_enr_all_dELS_scatter.png'), gp, dpi=450, units='in', width=4, height=4)
ggsave(paste0(imgpref, 'gwas_enr_all_dELS_scatter.pdf'), gp, width=4, height=4)


gp = ggplot(rwide, aes_string('All', 'EpiMap')) + 
    geom_point(alpha=.5) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_pubr()


srdf = ralldf[ralldf$rawlp > 30 & ralldf$set == 'All',]

unique(srdf, 




gp = ggplot(rwide, aes_string('EpiMap','random')) + 
    geom_point(alpha=.5) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    geom_abline() + 
    theme_pubr()

gp = ggplot(rwide, aes_string('dELS','random')) + 
    geom_point(alpha=.5) + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    geom_abline() + 
    theme_pubr()






