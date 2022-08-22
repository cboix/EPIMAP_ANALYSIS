#!/usr/bin/R
# ------------------------------------------------------
# Get pvalue cutoffs for full, modules, tree using FWER:
# Updated: 04/13/21
# ------------------------------------------------------
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
tol = 2500

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

# ---------------------------
# Load catalog / filter uids:
# ---------------------------
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

# --------------------------------------------------
# Load permuted matrices, perform the FWER analysis:
# --------------------------------------------------
apref = 'cons_parent'
NCAT=1000
for (setpref in c('', '_modules','_epigenomes')){
    fullrda = paste0(perdir, epref, apref, setpref, '_permuted_', 'allaggregated', '_', NCAT, '_lreg_pmat.Rda')
    if (!file.exists(fullrda)){
        # Get all files:
        pregex = paste0(epref, apref, setpref, '_permuted_0.*', '_', NCAT, '_lreg_pmat.Rda')
        pmat.files = list.files(path=perdir, pattern=pregex)
        load(paste0(perdir, pmat.files[2]))
        NF = length(pmat.files)
        NN = ncol(pmat)
        all.pcut = NULL
        nsvec = NULL
        all.pmat = matrix(1, nrow=NF * nrow(pmat), ncol=ncol(pmat))
        # Load in one p-value matrix for 1000 permuted GWAS for each SNP-size
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
            nsvec = c(nsvec, rep(as.numeric(nsnp), nrow(pmat)))
        }
        gc()
        save(all.pmat, all.pcut, nsvec, NF, NN, file=fullrda)
    } else {
        load(fullrda)
    }

    # For FWER (minimum p-value in each rep, then 1-10% of those):
    min.cut = apply(all.pmat,1, min)
    mcmat = matrix(min.cut, nrow=243, ncol=1000, byrow=T)
    # At each quantile:
    qmat = t(apply(mcmat, 1, probs=c(0.01,0.02, 0.05, 0.1), quantile))
    l10qmat = -log10(qmat) # FWER cutoffs, by number of SNPs
    nsnplist = unique(nsvec)

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
        # Save the FWER corrected matrix:
        qtstr = sub('%','',qt)
        write.table(corr.mat, paste0(gtdir, 'corrected_pernsnp', setpref, '_FWER_',qtstr,'.tsv'), sep="\t", row.names=T, col.names=T, quote=F)
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
        # Save the FWER corrected matrix:
        qtstr = sub('%','',qt)
        write.table(corr.mat, paste0(gtdir, 'corrected_pernsnp', setpref, '_FWER_',qtstr,'_alluids.tsv'), sep="\t", row.names=T, col.names=T, quote=F)
    }
    # See which are kept (many at 2-10 SNPs at 10% FDR):
    table(pvnum[pass.names])

    # ----------------------------
    # Comparison to null catalogs:
    # ----------------------------
    pvs = as.numeric(pvnum[kept.uids]) # To match distr of #pvals for null catalog
    lrmat = rmat[kept.uids,]
    # FWER cutoffs:
    mapv = sapply(pvs, function(x){which(nsnplist == x)}) # Map to FWER values:
    cutmat = l10qmat[mapv,]
    # For sampling:
    imat = sapply(pvs, function(i){ idx = which(nsvec == i); })
    NNULL = 1000;
    nullstatdf = c()
    nint = list()
    ngwas = list()
    for (qt in colnames(cutmat)){
        cutqt = cutmat[,qt]
        zmat = sweep(lrmat, 1, cutqt, '/')
        zmat = 1 * (zmat >= 1)
        nint[[qt]] = sum(zmat)
        ngwas[[qt]] = sum(apply(zmat, 1, sum) > 0)
    }
    # Run null:
    for (j in 1:NNULL){
        cat(j,"\n")
        # Create catalog by matching null GWAS to each GWAS from nsvec:
        set.seed(j)
        ind = apply(imat, 2, function(x){ sample(x,1)})
        null.lrmat = -log10(all.pmat[ind,])
        for (qt in colnames(cutmat)){
            cutqt = cutmat[,qt]
            null.zmat = sweep(null.lrmat, 1, cutqt, '/')
            null.zmat = 1 * (null.zmat >= 1)
            # E[FDR] for this run:
            efdr.int = sum(null.zmat) / nint[[qt]]
            efdr.gwas = sum(apply(null.zmat, 1, sum) > 0) / ngwas[[qt]]
            nullstatdf = rbind(nullstatdf, data.frame(qt=qt, eint = efdr.int, egwas=efdr.gwas, j=j))
        }
    }
    nullstatdf$qt = factor(nullstatdf$qt, levels=colnames(cutmat))

    # Plot the distribution of EFDR:
    gplot = ggplot(nullstatdf, aes(qt, eint)) + 
        geom_boxplot() + 
        theme_pubr() + 
        geom_hline(yintercept=c(.01,.02,.05,.1), color='red',lty='dashed') + 
        scale_y_continuous(expand=c(0,0), labels=scales::percent) +
        labs(x='FWER cutoff', y='E[FDR] on Number of Significant GWAS x Samples')
    ggsave(paste0(img, 'efdr',setpref,'_null_catalogs_at_intersection_level.png'), gplot, dpi=400, units='in', width=5, height=8)

    # Plot the distribution of EFDR:
    gplot = ggplot(nullstatdf, aes(qt, egwas)) + 
        geom_boxplot() + 
        theme_pubr() + 
        geom_hline(yintercept=c(.01,.02,.05,.1), color='red',lty='dashed') + 
        scale_y_continuous(expand=c(0,0), labels=scales::percent) + 
        labs(x='FWER cutoff', y='E[FDR] on Number of Significant GWAS')
    ggsave(paste0(img, 'efdr',setpref,'_null_catalogs_at_gwas_level.png'), gplot, dpi=400, units='in', width=5, height=8)

    # Mean # significant tests + gwas per null catalog: 
    gdf = data.frame(nint=unlist(nint), ngwas=unlist(ngwas))
    gdf$qt = rownames(gdf)
    edf = aggregate(cbind(eint,egwas) ~ qt, nullstatdf, mean)
    edf = merge(edf, gdf)
    edf$null.nint = edf$eint * edf$nint
    edf$null.ngwas = edf$egwas * edf$ngwas
    print(edf)

    write.table(nullstatdf, paste0(gtdir, 'efdr',setpref,'_stats_on_null_catalogs.tsv'), sep="\t", row.names=F, col.names=T, quote=F)
    
}



