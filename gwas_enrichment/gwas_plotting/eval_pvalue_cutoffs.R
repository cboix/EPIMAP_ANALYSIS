#!/usr/bin/R
# -----------------------------------------------
# Look at pvalue cutoffs for full, modules, tree:
# -----------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
library(ggplot2)
library(ggpubr)

# Arguments for loading data:
usetree = 'enhancers'
# usetree = 'roadmap'  # For old epi only
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions

# Load in + process all of the relevant matrices/datasets:
# commandArgs <- function(trailingOnly=TRUE){
#     c(usetree, tol, singlematch, plotting.only) }
# source(paste0(bindir, 'load_gwastree_analysis.R'))
# rm(dflist)
# print(paste("Images go to: ", treeimgpref))

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
imgdir = paste0(img, "gwas_tree_analysis/examples/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
treeimgpref = paste0(imgdir, usetree, '_e', tol, '_')

# Under dbdir:
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, perdir)
system(cmd)

# Logistic regression arguments:
type = 'cons'
against = 'parent'
weighted = FALSE
apref = paste0(type, '_', against)
if (weighted){ 
    weights = sqrt(1 / matmarg[,2])
    apref = paste0(apref, '_weighted')
} else {
    weights = NULL
}
NCAT=1000


gwcatfile = 'gwascatalog_may03_2019_noquotes.txt'
gwrdafile = sub(".txt", ".Rda", gwcatfile)
load(gwrdafile)
NIND = 10000
NSNP = 10
nsnpdf = aggregate(pValue ~ uid, gwdf, length)
gwssdf = aggregate(sampsize ~ uid, gwssdf, max)
gwssdf = merge(gwssdf, nsnpdf)
keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]
kept.uids =  sort(unique(as.character(keptgw$uid)))
NUID = length(kept.uids)

# Counts of snps:
# nsnpdf = aggregate(pValue ~ uid, gwdf, length)
# npref = paste0(apref, sprintf("_permuted_%06d", nsnps))
# pmatfile = paste0(perpref, npref, '_', NCAT, '_lreg_pmat.Rda')
# setpref='_epigenomes'
setpref=''
for (setpref in c('', '_modules','_epigenomes')){
    print(setpref)
    # Get all files:
    pregex = paste0(epref, apref, setpref, '_permuted_.*', '_', NCAT, '_lreg_pmat.Rda')
    pmat.files = list.files(path=perdir, pattern=pregex)
    load(paste0(perdir, pmat.files[2]))
    NF = length(pmat.files)
    NN = ncol(pmat)
    all.pcut = NULL
    nsvec = NULL
    all.pmat = matrix(1, nrow=NF * nrow(pmat), ncol=ncol(pmat))
    for (i in 1:NF){
        print(i)
        pfile = pmat.files[i]
        ppref = sub(paste0(epref, apref, setpref, '_permuted_'),'', pfile)
        nsnp = sub(paste0('_', NCAT, '_lreg_pmat.Rda'), '', ppref)
        # Load and aggregate data:
        tc = try(load(paste0(perdir, pfile)))
        if (class(tc) == 'try-error'){
            pmat[] = 1
        }
        pidx = ((i-1)*nrow(pmat) + 1):((i)*nrow(pmat))
        all.pmat[pidx,] = pmat
        all.pcut = c(all.pcut, apply(pmat, 2, min))
        nsvec = c(nsvec, rep(as.numeric(nsnp), ncol(pmat)))
    }
    gc()

    # For FWER (minimum p-value in each rep, then 1-10% of those):
    min.cut = apply(all.pmat,1, min)
    mcmat = matrix(min.cut, nrow=243, ncol=1000, byrow=T)
    # At each quantile:
    qmat = t(apply(mcmat, 1, probs=c(0.01,0.02, 0.05, 0.1), quantile))
    l10qmat = -log10(qmat) # FWER cutoffs, by number of SNPs
    nsnplist = unique(nsvec)

    # Other metrics:
    nsdf = data.frame(p=-log10(all.pcut), n=nsvec)
    gplot = ggplot(nsdf, aes(factor(n), p)) + 
        geom_boxplot() + theme_pubr()
    
    pdf = data.frame(pmat)
    pdf$row = 1:nrow(pdf)
    pdf = gather(pdf, col, val, -row)
    pdf$val = -log10(pdf$val)
    mpdf = aggregate(val ~ col, pdf,median)
    mpdf = mpdf[order(mpdf$val), ]
    pdf$col = factor(pdf$col, levels=mpdf$col)

    gplot = ggplot(pdf, aes(factor(col), val)) + 
        geom_boxplot() + 
        theme_pubr()

    if (setpref != ''){
        flatset = sub("_", "", setpref)
        enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
        load(enhsetfile)
        ldf = data.frame(col=paste0('X', 1:length(lenlist)), len=lenlist)
        mpdf = merge(ldf, mpdf)
        gplot = ggplot(mpdf, aes(val, len)) +
            geom_point() + theme_pubr()
    }

    # Get the maximal p-val per row:
    # perm.maxp = apply(all.pmat, 1, min)
    perm.maxp = apply(all.pmat, 1, min)
    l10maxp = -log10(perm.maxp); gc()
    l10fdrmax2 = quantile(l10maxp, 0.98, na.rm=TRUE); gc()  # 50
    l10fdrmax = quantile(l10maxp, 0.99, na.rm=TRUE); gc()  # 58
    l10fdrmax1p = quantile(l10maxp, 0.999, na.rm=TRUE); gc()
    # Of all, only 42 pass l10fdrmax, 57 pass l10fdrmax2
    cat(l10fdrmax2,l10fdrmax,l10fdrmax1p)

    # Find FDR cutoff per node:
    fdrper = apply(all.pmat, 2, function(x, alpha=0.01){quantile(x, alpha)})
    fdrper1p = apply(all.pmat, 2, function(x, alpha=0.001){quantile(x, alpha)})
    l10pvec = as.numeric(-log10(all.pmat)); gc()
    l10fdrall = quantile(l10pvec, 0.99, na.rm=TRUE); gc() # 10.8116
    l10fdrall5p = quantile(l10pvec, 0.995, na.rm=TRUE); gc() # 
    l10fdrall1p = quantile(l10pvec, 0.999, na.rm=TRUE); gc() # 22.09788

    # Note that removing 1-10 yields about the same results (for tree):
    # l10pvec = as.numeric(-log10(all.pmat[10000:243000,])); gc()
    # l10fdrall = quantile(l10pvec, 0.99, na.rm=TRUE); gc() # 10.96692 
    # l10fdrall1p = quantile(l10pvec, 0.999, na.rm=TRUE); gc() # 22.34395
    l10fdrper = -log10(fdrper) #pernode 
    l10fdrper1p = -log10(fdrper1p) 
    l10fdr.each = -log10(all.pcut) #foreach
    xlim = c(0, max(l10fdr.each))
    breaks = seq(0, xlim[2], length.out=100)
    gc()

    l10fdrper[l10fdrper < 4] = 4
    l10fdrper1p[l10fdrper1p < 4] = 4

    # Write out the cutoffs:
    pcutdf = rbind(data.frame(node=1:NN, cut=l10fdrper, type='Per 1%'),
                   data.frame(node=1:NN, cut=l10fdrper1p, type='Per 0.1%'),
                   data.frame(node=0, cut=l10fdrall, type='All 1%'),
                   data.frame(node=0, cut=l10fdrall1p, type='All 0.1%'))
    write.table(pcutdf, paste0(gtdir, 'cutoffs_l10', setpref, '_all_pernode.tsv'),
                quote=F, row.names=F, col.names=T, sep="\t")


    # Plot the FDR cutoff distributions:
    png(paste0(img, 'pval_cutoffs', setpref, '_with_orig_catalogs_distr.png'), width=6, height=6, units='in', res=450)
    layout(matrix(1:2, nrow=2, ncol=1))
    par(mar=c(2,2,2,0))
    par(yaxs='i')
    hist(l10fdrper, breaks, xlim=xlim, col= 'darkgrey', border='white', main='Enhancer Tree: FDR cutoffs (1 per node)\nRed line = 1 for full catalog')
    abline(v=l10fdrall, lwd=2, lty='dashed', col='red')
    hist(l10fdr.each, breaks, xlim=xlim, col='darkgrey', border='white', main='Cutoffs - for each node + gwas')
    dev.off()

    # Plot changes on ggplot:
    lmat = matrix(l10fdr.each, nrow=NF, byrow=TRUE)
    perdf = data.frame(n=1:NN, p=l10fdrper)
    eachdf = data.frame(n=rep(1:NN, NF), p=l10fdr.each)

    perdf = perdf[order(perdf$p),]
    lvl.n = perdf$n
    perdf$n = factor(perdf$n, levels=lvl.n)
    eachdf$n = factor(eachdf$n, levels=lvl.n)


    # Plot the FDR cutoffs - raw values
    gplot = ggplot(eachdf, aes(n, p)) + 
        geom_boxplot(color='grey35', outlier.size=0.15, lwd=.15) + 
        geom_point(data=perdf, aes(n, p), color='goldenrod1', cex=.25) + 
        theme_pubr() + labs(x='Node', y='-log10 p-value for FDR < 1% cutoff') + 
        theme(axis.text.x=element_blank()) + 
        geom_hline(yintercept=l10fdrall, color='firebrick', lwd=1, lty='dashed')
    ggsave(paste0(img, 'pval_cutoffs', setpref, '_with_orig_catalogs_boxplots_pernode.png'), gplot, width=20, height=6, units='in', dpi=450)


    # ------------------------
    # Get the actual p-values:
    # ------------------------
    suffix = '.Rda' 
    gwcatfile = 'gwascatalog_may03_2019_noquotes.txt'
    gwrdafile = sub(".txt", ".Rda", gwcatfile)
    load(gwrdafile)
    uids = sort(as.character(unique(gwdf$uid)))
    if (setpref != ''){
        enh.files = list.files(path=regdir, pattern=paste0(epref, apref, setpref, ".*_hg", suffix))
        fnum = sub(paste0("^.*", apref, ".*_"), "", sub("_hg.*Rda", "", enh.files[length(enh.files)]))
    } else {
        enh.files = list.files(path=regdir, pattern=paste0(epref, apref, setpref, ".*_lreg", suffix))
        fnum = sub(paste0("^.*", apref, ".*_"), "", sub("_lreg.*Rda", "", enh.files[length(enh.files)]))
    }

    # Load in files:
    NUID = as.numeric(fnum)
    NN = ncol(pmat)
    aggpvfile = paste0(gtdir, 'agg_raw_pvals', setpref, '_20200218.Rda')
    if (!file.exists(aggpvfile)){
        rmat = matrix(0, nrow=NUID, ncol=NN)
        for (rfile in enh.files){
            if (setpref != ''){
                fnum = sub(paste0("^.*", apref, ".*_"), "", sub("_hg.*Rda", "", rfile))
            } else {
                fnum = sub(paste0("^.*", apref, ".*_"), "", sub("_lreg.*Rda", "", rfile))
            }
            fnum = as.numeric(fnum)
            load(paste0(regdir, rfile))
            if (class(ll) == 'list'){
                # rmat[fnum, ] = ll$rawlp  # Was broken.
                rmat[fnum, ] = -log10(ll$df$pout)
            }
        }
        rownames(rmat) = c(uids, paste0('X', 1:(nrow(rmat)-5454)))
        save(rmat, file=aggpvfile)
    } else {
        load(aggpvfile)
    }

    # mp = apply(rmat, 1, max)
    # print(sum(mp > l10fdrmax))

    # ---------------------
    # Thresholding by FWER:
    # ---------------------
    # Try FWER by number of pvalues:
    lrmat = rmat[kept.uids,]
    pvnum = keptgw[,'pValue']
    names(pvnum) = keptgw[,'uid']
    pvs = pvnum[rownames(lrmat)]
    mapv = sapply(pvs, function(x){which(nsnplist == x)}) # Map to FWER values:
    cutmat = l10qmat[mapv,]
    for (qt in colnames(cutmat)){
        cat(qt, '\t')
        cutqt = cutmat[,qt]
        zmat = sweep(lrmat, 1, cutqt, '/')
        zmat = 1 * (zmat >= 1)
        pass.uid = sum(apply(zmat, 1, sum) > 0)
        pass.loc = sum(apply(zmat, 2, sum) > 0)
        pass.names = rownames(zmat)[apply(zmat, 1, sum) > 0]
        cat(pass.uid, '\t', pass.loc,'\n')
        corr.mat = (zmat * lrmat)[pass.names,]
        # cmat = corr.mat[,apply(corr.mat, 2, sum) > 0]
        # Save the FWER corrected matrix:
        qtstr = sub('%','',qt)
        # write.table(corr.mat, paste0(gtdir, 'corrected_pernsnp', setpref, '_FWER_',qtstr,'.tsv'), sep="\t", row.names=T, col.names=T, quote=F)
    }

    # What happens if we keep all uids (except only 1 snp)?
    # FWER should do a good job of correcting bad GWAS out (esp. single-overlap)
    pvnum = nsnpdf[nsnpdf$pValue > 1,'pValue']
    names(pvnum) = nsnpdf[nsnpdf$pValue > 1,'uid']
    ann.uid = rownames(rmat)[rownames(rmat) %in% names(pvnum)]
    lrmat = rmat[ann.uid,]
    pvs = pvnum[rownames(lrmat)]
    mapv = sapply(pvs, function(x){which(nsnplist == x)}) # Map to FWER values:
    cutmat = l10qmat[mapv,]
    for (qt in colnames(cutmat)){
        cat(qt, '\t')
        cutqt = cutmat[,qt]
        zmat = sweep(lrmat, 1, cutqt, '/')
        zmat = 1 * (zmat >= 1)
        pass.uid = sum(apply(zmat, 1, sum) > 0)
        pass.loc = sum(apply(zmat, 2, sum) > 0)
        pass.names = rownames(zmat)[apply(zmat, 1, sum) > 0]
        cat(pass.uid, '\t', pass.loc,'\n')
        corr.mat = (zmat * lrmat)[pass.names,]
        # cmat = corr.mat[,apply(corr.mat, 2, sum) > 0]
        # Save the FWER corrected matrix:
        qtstr = sub('%','',qt)
        write.table(corr.mat, paste0(gtdir, 'corrected_pernsnp', setpref, '_FWER_',qtstr,'_alluids.tsv'), sep="\t", row.names=T, col.names=T, quote=F)
    }
    # See which are kept (many at 2-10 SNPs at 10% FDR):
    table(pvnum[pass.names])

    # -------------------------
    # Look at cutoffs - overall
    # -------------------------
    rmatall = rmat
    rmatall[rmatall < l10fdrall] = 0
    rmatall[rmatall >= l10fdrall] = 1
    sum(apply(rmatall, 1, max))
    # TOTAL: 984 - tree
    sum(apply(rmat >= l10fdrall, 1, max))
    sum(apply(rmat >= l10fdrall, 2, max))
    all1.mat = rmat >= l10fdrall
    # TOTAL: 199 - tree
    sum(apply(rmat >= l10fdrall1p, 1, max))
    sum(apply(rmat >= l10fdrall1p, 2, max))
    allp1.mat = rmat >= l10fdrall1p

    # Select only passing 10k + 10: 
    sum(apply(rmat[kept.uids,] >= l10fdrall, 1, max))
    sum(apply(rmat[kept.uids,] >= l10fdrall1p, 1, max))

    # ---------
    # Per node:
    # ---------
    reorder_mat = function(mat){
        dt <- dist(mat, 'eJaccard')
        ht <- hclust(dt, method='complete')
        cocl <- order.optimal(dt, ht$merge)$order
        dt2 <- dist(t(mat), 'eJaccard')
        ht <- hclust(dt2, method='complete')
        rocl <- order.optimal(dt2, ht$merge)$order
        mat = mat[cocl,rocl]
    }

    rmatper = rmat
    fper = l10fdrper
    fper[fper < 4] = 4
    rmatper = sweep(rmatper, 2, fper, '/')
    # TOTAL: 1453 - tree
    sum(apply(rmatper > 1, 1, max))
    sum(apply(rmatper > 1, 2, max))
    node1.mat = rmatper

    # Select only passing 10k + 10 (703):
    sum(apply(rmatper[kept.uids,] > 1, 1, max))

    idx = which(apply(rmatper > 1, 1, max) > 0)
    mat = 1 * (rmatper[idx,] > 1)
    mat = reorder_mat(mat)
    png(paste0(img, 'eval_keptgwas_pernode_fdr_1', setpref, '.png'), width=6, height=6, units='in', res=450)
    par(mar=c(1,1,1,0))
    image(mat, axes=F, ylab='', xlab='', col=c('white','grey25'))
    box(lwd=.5)
    mtext(paste0(nrow(mat), ' significant GWAS'), side=1, line=0)
    mtext(paste0(ncol(mat), ' tested sets'), side=2, line=0)
    mtext('FDR < 1%',side=3)
    dev.off()

    # Per node, 0.1%
    rmatper = rmat
    fper = l10fdrper1p
    fper[fper < 4] = 4
    rmatper = sweep(rmatper, 2, fper, '/')
    # TOTAL: 531 - tree 
    sum(apply(rmatper > 1, 1, max))
    sum(apply(rmatper > 1, 2, max))
    nodep1.mat = rmatper
    # Select only passing 10k + 10 (399)
    sum(apply(rmatper[kept.uids,] > 1, 1, max))

    idx = which(apply(rmatper > 1, 1, max) > 0)
    mat = 1 * (rmatper[idx,] > 1)
    mat = reorder_mat(mat)
    png(paste0(img, 'eval_keptgwas_pernode_fdr_p1', setpref, '.png'), width=6, height=6, units='in', res=450)
    par(mar=c(1,1,1,0))
    image(mat, axes=F, ylab='', xlab='', col=c('white','grey25'))
    box(lwd=.5)
    mtext(paste0(nrow(mat), ' significant GWAS'), side=1, line=0)
    mtext(paste0(ncol(mat), ' tested sets'), side=2, line=0)
    mtext('FDR < 0.1%',side=3)
    dev.off()


    # --------------------------------------------------------------
    # Plot the properties of the kept matrices for tree enrichments:
    # --------------------------------------------------------------
    # NIND = 10000
    # NSNP = 10
    # nsnpdf = aggregate(pValue ~ uid, gwdf, length)
    # gwssdf = aggregate(sampsize ~ uid, gwssdf, max)
    # gwssdf = merge(gwssdf, nsnpdf)
    # keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]
    # kept.uids =  sort(unique(as.character(keptgw$uid)))
    # NUID = length(kept.uids)
    # gwdf = gwdf[gwdf$uid %in% kept.uids,]
    # gwgr = GRanges(paste0('chr', gwdf$chrom), IRanges(gwdf$chromStart, gwdf$chromEnd))

    load(gwrdafile)
    nsnpdf = aggregate(pValue ~ uid, gwdf, length)
    # NOTE: Some duplicates exist.
    gwssdf = aggregate(sampsize ~ uid, gwssdf, max)
    gwssdf = unique(merge(gwssdf, nsnpdf))
    pass.f1 = names(which(apply(node1.mat > 1, 1, max) > 0))
    pass.fp1 = names(which(apply(nodep1.mat > 1, 1, max) > 0))
    pass.all.f1 = names(which(apply(all1.mat == 1, 1, max) > 0))
    pass.all.fp1 = names(which(apply(allp1.mat == 1, 1, max) > 0))
    # Make cutoffs:
    gw1 = gwssdf
    gwp1 = gwssdf
    gwa1 = gwssdf
    gwap1 = gwssdf
    gw1$pass = gw1$uid %in% pass.f1
    gwp1$pass = gwp1$uid %in% pass.fp1
    gwa1$pass = gw1$uid %in% pass.all.f1
    gwap1$pass = gwp1$uid %in% pass.all.fp1
    gw1$fdr.cutoff = '1%'
    gwp1$fdr.cutoff = '0.1%'
    gwa1$fdr.cutoff = 'All 1%'
    gwap1$fdr.cutoff = 'All 0.1%'
    gwtot = rbind(gw1, gwp1, gwa1, gwap1)

    # By sample size:
    gplot = ggplot(gwtot, aes(sampsize, fill=pass)) +
        facet_wrap(~fdr.cutoff, scales='free') + 
        theme_pubr() + 
        scale_x_log10(labels=scales::comma, expand=c(0,0)) + 
        scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
        labs(y='Percent of GWAS', x='Sample size') + 
        scale_fill_manual(values=c('grey90','grey50'), name='Significant at FDR cutoff') + 
        geom_histogram(position='fill') + 
        geom_vline(xintercept=c(1e4, 2e4), lwd=1, lty='dashed', col='red')
    ggsave(paste0(img, 'eval_keptgwas_pernode', setpref, '_bysamplesize.png'), gplot, width=12, height=4, units='in', dpi=450)

    # By sample size:
    gplot = ggplot(gwtot, aes(pValue, fill=pass)) +
        facet_wrap(~fdr.cutoff, scales='free') + 
        theme_pubr() + 
        scale_x_log10(labels=scales::comma, expand=c(0,0)) + 
        scale_y_continuous(labels=scales::percent, expand=c(0,0)) + 
        labs(y='Percent of GWAS', x='Sample size') + 
        scale_fill_manual(values=c('grey90','grey50'), name='Significant at FDR cutoff') + 
        geom_histogram(position='fill') + 
        # geom_histogram() + 
        geom_vline(xintercept=c(10, 50), lwd=1, lty='dashed', col='red')
    ggsave(paste0(img, 'eval_keptgwas_pernode', setpref, '_bynsnp.png'), gplot, width=12, height=4, units='in', dpi=450)

    subtot = gwtot[gwtot$pValue >= 10 & gwtot$sampsize >= 1e4,]
    aggregate(uid ~ pass + fdr.cutoff, subtot, length)

    subtot = gwtot[gwtot$sampsize >= 1e4,]
    aggregate(uid ~ pass + fdr.cutoff, subtot, length)

    subtot = gwtot[gwtot$pValue >= 10,]
    aggregate(uid ~ pass + fdr.cutoff, subtot, length)

    # Plot these:
    subtot = gwtot[gwtot$pValue >= 10 & gwtot$sampsize >= 1e4,]

    aa = subtot[subtot$pass & subtot$fdr.cutoff=='1%',]
    mat = 1 * (node1.mat[aa$uid,] > 1)
    mat = reorder_mat(mat)
    png(paste0(img, 'eval_keptgwas_pernode_fdr_1_10snp_10kind', setpref, '.png'), width=6, height=18, units='in', res=450)
    par(mar=c(1,10,1,1))
    image(t(mat), axes=F, ylab='', xlab='', col=c('white','grey25'))
    box(lwd=.5)
    mtext(paste0(nrow(mat), ' significant GWAS (10k+ individ. and 10+ lead SNPs)'), side=4, line=0)
    mtext(paste0(ncol(mat), ' tested sets'), side=1, line=0)
    mtext('Test FDR < 1%',side=3)
    text(x=par()$usr[1] - 0.002 * diff(par()$usr[1:2]),
         y=seq(0,1, length.out=nrow(mat)),
         labels=sapply(rownames(mat), function(x){sub('^[0-9]* - ','',x)}),
         adj=1, xpd=TRUE, cex=.15)
    dev.off()

    # Plot these:
    aa = subtot[subtot$pass & subtot$fdr.cutoff=='0.1%',]
    mat = 1 * (nodep1.mat[aa$uid,] > 1)
    mat = reorder_mat(mat)
    png(paste0(img, 'eval_keptgwas_pernode_fdr_p1_10snp_10kind', setpref, '.png'), width=6, height=11, units='in', res=450)
    par(mar=c(1,10,1,1))
    image(t(mat), axes=F, ylab='', xlab='', col=c('white','grey25'))
    box(lwd=.5)
    mtext(paste0(nrow(mat), ' significant GWAS (10k+ individ. and 10+ lead SNPs)'), side=4, line=0)
    mtext(paste0(ncol(mat), ' tested sets'), side=1, line=0)
    mtext('Test FDR < 0.1%',side=3)
    text(x=par()$usr[1] - 0.002 * diff(par()$usr[1:2]),
         y=seq(0,1, length.out=nrow(mat)),
         labels=sapply(rownames(mat), function(x){sub('^[0-9]* - ','',x)}),
         adj=1, xpd=TRUE, cex=.15)
    dev.off()

    # Plot these:
    aa = subtot[subtot$pass & subtot$fdr.cutoff=='All 1%',]
    mat = 1 * (all1.mat[aa$uid,] == 1)
    mat = reorder_mat(mat)
    png(paste0(img, 'eval_keptgwas_pernode_fdr_all1_10snp_10kind', setpref, '.png'), width=6, height=11, units='in', res=450)
    par(mar=c(1,10,1,1))
    image(t(mat), axes=F, ylab='', xlab='', col=c('white','grey25'))
    box(lwd=.5)
    mtext(paste0(nrow(mat), ' significant GWAS (10k+ individ. and 10+ lead SNPs)'), side=4, line=0)
    mtext(paste0(ncol(mat), ' tested sets'), side=1, line=0)
    mtext('Test FDR < 0.1%',side=3)
    text(x=par()$usr[1] - 0.002 * diff(par()$usr[1:2]),
         y=seq(0,1, length.out=nrow(mat)),
         labels=sapply(rownames(mat), function(x){sub('^[0-9]* - ','',x)}),
         adj=1, xpd=TRUE, cex=.15)
    dev.off()

    # Plot these:
    aa = subtot[subtot$pass & subtot$fdr.cutoff=='All 0.1%',]
    mat = 1 * (allp1.mat[aa$uid,] == 1)
    mat = reorder_mat(mat)
    png(paste0(img, 'eval_keptgwas_pernode_fdr_allp1_10snp_10kind', setpref, '.png'), width=6, height=11, units='in', res=450)
    par(mar=c(1,10,1,1))
    image(t(mat), axes=F, ylab='', xlab='', col=c('white','grey25'))
    box(lwd=.5)
    mtext(paste0(nrow(mat), ' significant GWAS (10k+ individ. and 10+ lead SNPs)'), side=4, line=0)
    mtext(paste0(ncol(mat), ' tested sets'), side=1, line=0)
    mtext('Test FDR < 0.1%',side=3)
    text(x=par()$usr[1] - 0.002 * diff(par()$usr[1:2]),
         y=seq(0,1, length.out=nrow(mat)),
         labels=sapply(rownames(mat), function(x){sub('^[0-9]* - ','',x)}),
         adj=1, xpd=TRUE, cex=.15)
    dev.off()

    write.table(node1.mat, paste0(gtdir, 'corrected_pernode', setpref, '_fdr_1.tsv'), sep="\t", row.names=T, col.names=T, quote=F)
    write.table(nodep1.mat, paste0(gtdir, 'corrected_pernode', setpref, '_fdr_p1.tsv'), sep="\t", row.names=T, col.names=T, quote=F)
    write.table(all1.mat, paste0(gtdir, 'corrected_all', setpref, '_fdr_1.tsv'), sep="\t", row.names=T, col.names=T, quote=F)
    write.table(allp1.mat, paste0(gtdir, 'corrected_all', setpref, '_fdr_p1.tsv'), sep="\t", row.names=T, col.names=T, quote=F)

    # Filter the full p-value matrix, use the all1.mat to adjust p-values as qvalues:
    rmatall = rmat
    amat = (rmat >= l10fdrall)
    write.table(amat, paste0(gtdir, 'corrected_rawp_all', setpref, '_fdr_1.tsv'), sep="\t", row.names=T, col.names=T, quote=F)

}















