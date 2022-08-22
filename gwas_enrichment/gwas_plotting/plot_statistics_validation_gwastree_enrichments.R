#!/usr/bin/R
# -----------------------------------------------------
# Plot epigenomes (or modules) with GWAS
# -----------------------------------------------------
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
library(ggpubr)
library(viridis)
library(tidyr)
library(dplyr)
library(UpSetR)
library(scales)

# -----------------------
# Defaults for gwas data:
# -----------------------
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE
use.adj = TRUE
# use.strict = TRUE
use.strict = FALSE

# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict) }
source(paste0(bindir, 'load_validation_gwastree_enrichments.R'))

rm(cdll)

# Plot prefix:
imgdir = paste0(img, 'validation/')
cmd = paste("mkdir -p", imgpref)
system(cmd)
imgpref = paste0(imgdir, 'onecut_')

# -------------------------------------
# Load all gwas enrichment information:
# -------------------------------------
suffixes = c('_adj1000_10.Rda','_adj1000_5.Rda', '_adj1000_1.Rda')
fdrs = c(1, 0.5, 0.1)
cutoff = 3
all.regmats = list()
regdf = c()
all.snpmats = list()
nldf = c()
for (suffix in suffixes){
    print(suffix)
    snpfiles = list()
    sets = c('enh','rdm','mod','epi')
    snpfiles[['enh']] = paste0(regdir, eprefix, apref, '_logreg_all_wsnp', suffix)
    snpfiles[['rdm']] = paste0(regdir, rprefix, apref, '_logreg_all_wsnp', suffix)
    snpfiles[['mod']] = paste0(regdir, eprefix, apref, '_modules_hg_all_wsnp', suffix)
    snpfiles[['epi']] = paste0(regdir, eprefix, apref, '_epigenomes_hg_all_wsnp', suffix)
    regmats = list()
    snpmats = list()
    for (dset in sets){
        load(snpfiles[[dset]])
        colnames(rmat) = paste0('c', c(1:ncol(rmat) - 1))
        rdf = data.frame(rmat)
        rdf$uid = rownames(rdf)
        rdf = gather(rdf, cls, logpval, -uid)
        rdf = rdf[rdf$logpval > cutoff,]
        rdf$suffix = suffix
        rdf$dset = dset
        # SNPs:
        sfilt = (rmat > 0) * smat
        sdf = data.frame(sfilt)
        sdf$uid = rownames(sdf)
        sdf = gather(sdf, cls, nsnp, -uid)
        sdf = sdf[sdf$nsnp > 0,]
        rdf = merge(rdf, sdf, all.x=TRUE)
        rdf$nsnp[is.na(rdf$nsnp)] = 0
        rdf = merge(rdf, data.frame(cls=paste0('c', c(1:ncol(rmat) - 1)),
                                    clen=lens[[dset]]))
        regmats[[dset]] = rmat
        snpmats[[dset]] = (rmat > 0) * smat
        regdf = rbind(regdf, rdf)
        # Numbers of unique snps and intersections:
        sndf = data.frame(nlist)
        sndf$uid = uids[as.numeric(rownames(sndf))]
        sndf = sndf[sndf$n.int > 0,]
        sndf$suffix = suffix
        sndf$dset = dset
        nldf = rbind(nldf, sndf)
    }
    all.regmats[[suffix]] = regmats
    all.snpmats[[suffix]] = snpmats
}


# Alternate load: 
all.regmats = list()
all.snpmats = list()
setprefs = c('_', '_modules_','_epigenomes_')
dsetmap = c('enh','mod','epi')
names(dsetmap) = setprefs
regdf = c()
nldf = c()
for (setpref in setprefs){
    aggpvfile = paste0(gtdir, 'agg_raw_pvals_snps', setpref, '20200218.Rda')
    dset = dsetmap[setpref]
    load(aggpvfile)
    for (lvl in c('1%', '0.1%')){
        rmat = rmatlist[[lvl]] 
        smat = smatlist[[lvl]]
        colnames(rmat) = paste0('c', c(1:ncol(rmat) - 1))
        rdf = data.frame(rmat)
        rdf$uid = rownames(rdf)
        rdf = gather(rdf, cls, logpval, -uid)
        rdf = rdf[rdf$logpval > cutoff,]
        rdf$suffix = lvl
        rdf$dset = dset
        # SNPs:
        sfilt = (rmat > 0) * smat
        sdf = data.frame(sfilt)
        sdf$uid = rownames(sdf)
        sdf = gather(sdf, cls, nsnp, -uid)
        sdf = sdf[sdf$nsnp > 0,]
        rdf = merge(rdf, sdf, all.x=TRUE)
        rdf$nsnp[is.na(rdf$nsnp)] = 0
        rdf = merge(rdf, data.frame(cls=paste0('c', c(1:ncol(rmat) - 1)),
                                    clen=lens[[dset]]))
        regmats[[dset]] = rmat
        snpmats[[dset]] = (rmat > 0) * smat
        regdf = rbind(regdf, rdf)
        # Numbers of unique snps and intersections:
        sndf = data.frame(nhitlist[[lvl]])
        sndf$uid = rownames(sndf)
        sndf = sndf[sndf$n.int > 0,]
        sndf$suffix = lvl
        sndf$dset = dset
        nldf = rbind(nldf, sndf)
    }
}


# -----------------------------------------------
# GWAS with at least 10k individuals and 10+ snps
# -----------------------------------------------
NIND = 10000 # Remove all the low GWAS
NSNP = 10
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
# NOTE: Some duplicates exist.
gwssdf2 = aggregate(sampsize ~ uid, gwssdf, max)
gwssdf = merge(gwssdf, gwssdf2)
gwssdf = unique(merge(gwssdf, nsnpdf))
# 803 GWAS: 
keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]
keptgwas = unique(sort(keptgw$uid)) # All GWAS

# Used to be 20k (name of col)
regdf$in.20k = 1 * (regdf$uid %in% keptgwas)
nldf$in.20k = 1 * (nldf$uid %in% keptgwas)
regdf = regdf[regdf$nsnp > 0,]
regdf = merge(regdf, data.frame(suffix=c('1%','0.1%'), 
                                fdr.level=c(1, 0.1)))
regdf$fdr.level = factor(regdf$fdr.level, c(0.1, 1))

# Count number of snps:
gwdf$in.20k = 1 * (gwdf$uid %in% keptgwas)
sum(gwdf$in.20k) # 66801 SNPs (71,379 for 10k + 10)


# setlabs = c('Tree\n(All)', 'Tree\n(Roadmap)', 'Flat\n(Modules)', 'Flat\n(Epigenomes)')
# setlabs = c('Tree', 'Flat\n(Modules)', 'Flat\n(Epigenomes)')
setlabs = c('Tree', 'Modules', 'Flat')
sets = c('enh','mod','epi')
setmap = data.frame(dset=sets, dlab=setlabs)
setmap$dlab = factor(setmap$dlab, levels=rev(setlabs))
regdf = merge(regdf, setmap)
regdf$dlab = factor(regdf$dlab, levels=rev(setlabs))

nldf = merge(nldf, data.frame(suffix=c('1%','0.1%'), 
                                fdr.level=c(1, 0.1)))
nldf$fdr.level = factor(nldf$fdr.level, c(0.1, 1))
nldf = merge(nldf, setmap)
nldf$dlab = factor(nldf$dlab, levels=rev(setlabs))

pcols = brewer.pal(n=12,name="Paired")
# svalues = pcols[c(2,4,6,8)]
svalues = pcols[c(2,6,8)]
names(svalues) = setlabs

urdf = unique(regdf[regdf$in.20k==1, c('uid', 'dset', 'dlab', 'fdr.level')])
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
urdf = merge(urdf, nsnpdf)
urdf$pValue[urdf$pValue > 100] = 100
urdf$dlab = factor(urdf$dlab, levels=rev(setlabs))

gt = ggplot(urdf, aes(pValue, alpha=fdr.level, fill=dlab)) + 
    facet_grid(fdr.level ~ dlab) + 
    geom_histogram() +
    # geom_text(data=tdf, aes(dlab, uid + 30, label=uid, alpha=fdr.level), 
    #           position = position_dodge(width = .9), size=3) +
    scale_fill_manual(values = svalues, name='Method (Dataset)') + 
    # scale_y_continuous(labels = scales::comma) + 
    scale_alpha_manual(values=c(1, 0.75,0.5), name='FDR level') + 
    labs(x="Number of lead SNPs per GWAS", y="# GWAS with signif. enrichments") +
    theme_pubr()

ggsave(paste0(imgpref, 'validation_nsnp.png'), gt, dpi=350, units='in', width=10, height=8)
ggsave(paste0(imgpref, 'validation_nsnp.pdf'), gt, width=10, height=8)

csrdf = aggregate(uid ~ pValue + fdr.level + dlab, urdf, length)
csrdf = csrdf[order(csrdf$pValue, decreasing=TRUE),]
csrdf$cs = 0

for(f in unique(csrdf$fdr.level)){
    for(e in unique(csrdf$dlab)){
        ind = which(csrdf$dlab== e & csrdf$fdr.level == f)
        csrdf$cs[ind] = cumsum(csrdf$uid[ind])
    }
}


gt2 = ggplot(csrdf, aes(pValue, cs, color=dlab)) + 
    # facet_grid(fdr.level ~ dlab) + 
    facet_grid(. ~ fdr.level) + 
    geom_line() +
    geom_hline(yintercept=0) + 
    scale_color_manual(values = svalues, name='Method (Dataset)') + 
    ylim(0,900) + 
    labs(x="Number of lead SNPs per GWAS", y="Cumulative # GWAS with signif. enrichments") +
    theme_pubr()

ggsave(paste0(imgpref, 'validation_nsnp_cs.png'), gt2, dpi=350, units='in', width=10, height=5)
ggsave(paste0(imgpref, 'validation_nsnp_cs.pdf'), gt2, width=10, height=5)


# ------------------------
# Plot overall statistics:
# ------------------------
# 1. Barplots of number of GWAS (of all + of kept)
# 2. Boxplots of numbers of enhancers per set (in general)
# 3. Boxplots of enrichment ratio (number of snps per 1k enhancers)
# 4. Boxplots of percent of sets hit per GWAS (specificity)
# ------------------------

# 1. Barplots of number of GWAS (of all + of kept):
tform = formula(uid ~ fdr.level + dset + dlab)
tdf = aggregate(tform, regdf, function(x){length(unique(x))})
tdf = aggregate(tform, regdf[regdf$in.20k == 1,], function(x){length(unique(x))})
tdf$dlab = factor(tdf$dlab, levels=rev(setlabs))
g1 = ggplot(tdf, aes(dlab, uid, alpha=fdr.level, fill=dlab)) + 
    geom_bar(stat='identity', position='dodge') +
    geom_text(data=tdf, aes(dlab, uid + 30, label=uid, alpha=fdr.level), 
              position = position_dodge(width = .9), size=3) +
    scale_fill_manual(values = svalues, name='Method (Dataset)') + 
    scale_y_continuous(labels = scales::comma) + 
    scale_alpha_manual(values=c(1, 0.75,0.5), name='FDR level') + 
    labs(x="", y="# GWAS with signif. enrichments") +
    theme_pubr()

# 2. Boxplots of numbers of enhancers in nodes hit:
g2 = ggplot(regdf, aes(dlab, clen / 1000, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    labs(x="", y="# Enhancers (in 1,000s)") +
    scale_fill_manual(values = svalues) + 
    scale_y_log10() + 
    theme_pubr()  +
    theme(legend.position='none')

# 3. Boxplots of enrichment ratio (number of snps per 1k enhancers)
regdf$ratio = regdf$nsnp / regdf$clen
g3 = ggplot(regdf, aes(dlab, ratio * 1000, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology",
         y="# SNPs per 1,000 enhancers") +
    theme_pubr() + 
    scale_y_log10() + 
    theme(legend.position='none')


# 4. Boxplots of percent of sets hit per GWAS (specificity)
tform = formula(cls ~ uid + fdr.level + dset + dlab)
tdf = aggregate(tform, regdf, function(x){length(unique(x))})
tdf = aggregate(tform, regdf[regdf$in.20k == 1,], function(x){length(unique(x))})
tdf = merge(tdf, data.frame(nf=nclist, dset=names(nclist)))
g4 = ggplot(tdf, aes(dlab, cls / nf, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology",
         y="Percent of sets signif. per GWAS") +
    scale_y_log10(labels = scales::percent) + 
    theme_pubr() + 
    theme(legend.position='none')

g4b = ggplot(tdf, aes(dlab, cls, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology",
         y="Number of sets signif. per GWAS") +
    scale_y_log10() + 
    theme_pubr() + 
    theme(legend.position='none')

# aggregate(cls ~ dlab + fdr.level, tdf, median)
# Number of SNP recovered per GWAS (potentially many overlaps):
tdf2 = aggregate(nsnp ~ uid + dlab + dset + fdr.level, regdf, sum)
g5 = ggplot(tdf2, aes(dlab, nsnp, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology",
         y="Number of interactions overall") +
    # scale_y_log10(labels = scales::percent) + 
    scale_y_log10() + 
    theme_pubr() + 
    theme(legend.position='none')

# Number of SNP recovered per GWAS (potentially many overlaps):
g5b = ggplot(regdf, aes(dlab, nsnp, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology",
         y="Number of interactions per set") +
    # scale_y_log10(labels = scales::percent) + 
    scale_y_log10() + 
    theme_pubr() + 
    theme(legend.position='none')



# Number of SNP recovered per GWAS (potentially many overlaps):
tdf2 = aggregate(nsnp ~ uid + dlab + dset + fdr.level, regdf, sum)

g6 = ggplot(nldf, aes(dlab, n.int / tot.int, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology",
         y="% SNP-Enh Intersections captured") +
    # scale_y_log10(labels = scales::percent) + 
    scale_y_continuous(labels = scales::percent) + 
    # scale_y_log10() + 
    theme_pubr() + 
    theme(legend.position='none')

g6b = ggplot(nldf, aes(dlab, n.uniq / tot.uniq, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology",
         y="% SNPs captured") +
    scale_y_continuous(labels = scales::percent) + 
    theme_pubr() + 
    theme(legend.position='none')


g7 = ggplot(nldf, aes(dlab, n.int, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology", y="# SNP-Enh Intersections captured") +
    ylim(0,90) + 
    # scale_y_log10() + 
    theme_pubr() + 
    theme(legend.position='none')

g7b = ggplot(nldf, aes(dlab, n.uniq, alpha=fdr.level, fill=dlab)) + 
    scale_alpha_manual(values=c(1, 0.75,0.5)) + 
    geom_boxplot() +
    scale_fill_manual(values = svalues) + 
    labs(x="Methodology", y="# SNPs captured") +
    ylim(0,90) + 
    # scale_y_log10() + 
    theme_pubr() + 
    theme(legend.position='none')

# TODO: TURN THESE into figures:
snp.tdf = aggregate(n.uniq ~ dlab + fdr.level, nldf[nldf$in.20k == 1,], sum)
snp.tdf$dlab = factor(snp.tdf$dlab, levels=rev(setlabs))
g8 = ggplot(snp.tdf, aes(dlab, n.uniq, alpha=fdr.level, fill=dlab)) + 
    geom_bar(stat='identity', position='dodge') +
    geom_text(data=snp.tdf, aes(dlab, n.uniq + 1000,
                                label=comma_format()(n.uniq), alpha=fdr.level), 
              position = position_dodge(width = .9), size=3) +
    scale_fill_manual(values = svalues, name='Method (Dataset)') + 
    scale_y_continuous(labels = scales::comma) + 
    scale_alpha_manual(values=c(1, 0.75,0.5), name='FDR level') + 
    labs(x="", y="# Lead SNPs in enriched annotations") +
    theme_pubr()

# aggregate(n.int ~ dlab + fdr.level, nldf, sum)


# Number of unique SNPs per GWAS in enriched sets 
# - show better distributes.

# Aggregate all images:
ggarrange(g1 + rremove("x.text"),
          g2 + rremove("x.text"), 
          g3, 
          g8,
          # g6,
          # labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, align='hv', 
          common.legend = TRUE, legend = "bottom")

ggsave(paste0(imgpref, 'validation_statistics.png'), dpi=350, units='in', width=10, height=8)
ggsave(paste0(imgpref, 'validation_statistics.pdf'), width=10, height=8)

ggsave(paste0(imgpref, 'validation_statistics_small.png'), dpi=300, units='in', width=8, height=6)
ggsave(paste0(imgpref, 'validation_statistics_small.pdf'), width=8, height=6)

# Aggregate all images:
ggarrange(g1, g8,
          ncol = 2, nrow = 1, align='hv', 
          common.legend = TRUE, legend = "bottom")

ggsave(paste0(imgpref, 'validation_statistics_short.png'), dpi=350, units='in', width=10, height=4)
ggsave(paste0(imgpref, 'validation_statistics_short.pdf'), width=10, height=4)


# Extended (TODO: add n-unique SNP + % SNP as well)
gp = ggarrange(g1 + rremove("x.text"),
               g2 + rremove("x.text"), 
               g3 + rremove("x.text"),
               g8 + rremove("x.text"),
               g4 + rremove("x.text"), 
               g4b + rremove("x.text"),
               g5 + rremove("x.text"), 
               g5b + rremove("x.text"), 
               g6, g6b, 
               g7, g7b, 
               # labels = c("A", "B", "C", "D", "E", "F"),
               ncol = 4, nrow = 3, align='hv', 
               common.legend = TRUE, legend = "bottom")

ggsave(paste0(imgpref, 'validation_statistics_ext.png'), gp, dpi=350, units='in', width=18, height=12)
ggsave(paste0(imgpref, 'validation_statistics_ext.pdf'), gp, width=18, height=12)



# Aggregate all images:
gp = ggarrange(g7, g7b, ncol = 2, nrow = 1, align='hv', 
          common.legend = TRUE, legend = "bottom")
ggsave(paste0(imgpref, 'validation_statistics_numintuq.png'), gp, dpi=350, units='in', width=8, height=4)



# ----------------------------------
# Concordance between idenfications:
# ----------------------------------
# 1. Look at which id by epigenomes, which by modules. Do they overlap more with leaves/nodes?
# 2. Agreement of tissue(s) of action - which is top tissue, which are next.
# ----------------------------------
# 1. Look at which id by epigenomes, which by modules. Do they overlap more with leaves/nodes?
# Answer - no, but modules do identify a lot that stand alone - not ID in other methods
ulist = list()
fl = 0.1
for (dset in sets){
    subdf = regdf[regdf$fdr.level == fl & regdf$in.20k == 1,]
    ulist[[dset]] = unique(subdf$uid[subdf$dset == dset])
}

png(paste0(imgpref, 'upset_methods_fdr', fl ,'.png'), res=250, units='in', width=6, height=4)
upset(fromList(ulist), order.by = "freq",empty.intersections = "on" )
dev.off()

# Separate by leaf/node and repeat upset analysis:
ntdf = nodetissue
ntdf = ntdf[order(ntdf$node),]
ntdf$cls = paste0('c', ntdf$node - 1)
ntdf$dset = 'enh'
ntdf$isleaf = declist$isleaf


r2df = merge(regdf, ntdf[,c('cls','dset','isleaf')], all.x=TRUE)
r2df$dset2 = r2df$dset
r2df$dset2[r2df$dset == 'enh' & r2df$isleaf == 1] = 'enh leaf'
r2df$dset2[r2df$dset == 'enh' & r2df$isleaf == 0] = 'enh node'

ulist2 = list()
fl = 0.1
for (dset in c('enh leaf', 'enh node', sets[2:4])){
    subdf = r2df[r2df$fdr.level == fl & r2df$in.20k == 1,]
    ulist2[[dset]] = unique(subdf$uid[subdf$dset2 == dset])
}

png(paste0(imgpref, 'upset_methods_leafnode_fdr', fl, '.png'), res=250, units='in', width=6, height=4)
upset(fromList(ulist2), order.by = "freq",empty.intersections = "on" )
dev.off()


# Look at strength of pvals for mod/epi identified x leaf/node
mdf = unique(regdf[regdf$dset %in% c('mod', 'epi'), c('uid','dset', 'fdr.level')])
names(mdf)[2] = 'pivotset'
mdf = merge(r2df, mdf)
mdf = mdf[mdf$dset == 'enh',]
mdf$isleaf = factor(mdf$isleaf)
# No difference - leaf vs. not:
ggplot(mdf, aes(isleaf,logpval, fill=pivotset)) +
    scale_y_log10() + 
    theme_pubr() + 
    geom_boxplot()


# Look at uids only in modules (not particularly crucial GWAS)
mlist = ulist[['mod']]
olist = unique(regdf$uid[regdf$dset %in% sets[c(1,2,4)]])
setdiff(mlist, olist)


# 2. Agreement of tissue(s) of action - which is top tissue, which are next.


