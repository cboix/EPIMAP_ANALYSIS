#!/usr/bin/R
# -----------------------------------------------
# Calculate flat enrichments for sets of regions:
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
source(paste0(bindir, 'auxiliary_gwastree_functions.R'))

library(GenomicRanges)
library(dplyr)
library(cba)
library(ggplot2)
library(argparser)

# argv = list()
# Non-overlapping DHS list from WM:
# argv$uqfile='/broad/compbio/cboix/EPIMAP_ANALYSIS/db/DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_nonovl_any_coords_hg19.core.srt.txt'
# argv$bedfile='/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/cls_merge2_wH3K27ac100_nonovl_raw/motif_enrichment/observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_nonovl_on_mixed_impobs_full.bed.gz'
# argv$resultdir='/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/cls_merge2_wH3K27ac100_nonovl_raw/'
# argv$tagline=''

# From downloaded:
# argv = list()
# argv$uqfile = 'gwas_hg_stats/rawdata/masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt'
# argv$bedfile = 'gwas_hg_stats/rawdata/observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_full.bed.gz'
# argv$resultdir = 'gwas_hg_stats/rawdata/'
# argv$tagline=''

# Create argument parser:
p <- arg_parser("Calculate GWAS enrichments for sets of regions")
p <- add_argument(p, "--bedfile", help="Region assignments (with chr start end cls OR identifer cls)")
p <- add_argument(p, "--resultdir", help="Output directory")
p <- add_argument(p, "--uqfile", help="List of unique regions + identifier", default="")
p <- add_argument(p, "--tagline", help="tag line", default="")
# Parse the command line arguments:
argv <- parse_args(p)

# Add tagline if empty
if (argv$tagline == ""){
    argv$tagline = sub(".*/", "", argv$bedfile)
    argv$tagline = paste("File:", sub("_seed.*", "", argv$tagline))
}


# Ensure directory exists:
cmd = paste("mkdir -p", argv$resultdir, paste0(argv$resultdir, 'permuted'))
system(cmd)

tol = 2500  # Plus/minus distance - window for enhancer overlaps
# tol = 0

# Minimal hypergeometric test (for this type of data)
get_pval_hg_fast = function(suid, dflist, lenlist, 
                            qdf, NF=NF, cutp=CUTP, olen=NENH){
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
    return(-log10(pout))
}

# --------------------------
# Load in the gwas datasets:
# -------------------------------------------------
# Only keep GWAS with 10+ lead SNPs (after pruning)
# and with at least 10k individuals.
# -------------------------------------------------
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
if (length(grep("gz$", argv$bedfile)) > 0){
    df = read.delim(gzfile(argv$bedfile), header=F, stringsAsFactors=F)
} else {
    df = read.delim(argv$bedfile, header=F, stringsAsFactors=F)
}

# Read uq file create indexing:
if (argv$uqfile == ""){
    names(df)[1:4] = c('chr','start','end', 'cls')
    print(dim(df))

    # Get unique enhancers:
    uqdf = unique(df[,c('chr','start','end')])
    NENH = nrow(uqdf)
    uqdf$id = 1:NENH
    dmgr = GRanges(uqdf$chr, IRanges(uqdf$start - tol, uqdf$end + tol))
    df = merge(df, uqdf)
} else {
    names(df) = c('loc','cls')
    print(dim(df))

    if (length(grep("gz$", argv$uqfile)) > 0){
        uqdf = read.delim(gzfile(argv$uqfile), header=F, stringsAsFactors=F)
    } else {
        uqdf = read.delim(argv$uqfile, header=F, stringsAsFactors=F)
    }

    # Keep only ones present in uqdf:
    names(uqdf) = c('chr','start','end', 'loc')
    keep.uq = which(uqdf$loc %in% df$loc)
    uqdf = uqdf[keep.uq,]
    NENH = nrow(uqdf)
    uqdf$id = 1:NENH

    # Get the id column for df:
    indmap = uqdf$id
    names(indmap) = as.character(uqdf$loc)
    df$id = indmap[df$loc]
}
dmgr = GRanges(uqdf$chr, IRanges(uqdf$start - tol, uqdf$end + tol))
rm(uqdf)
gc()

# Turn into list:
clslist = sort(unique(df$cls))
enhsets = lapply(clslist, function(x){df$id[df$cls == x]})
NF = length(enhsets)
lenlist = sapply(enhsets, length)
gc()

# Calculate intersections with GWAS:
qdf = suppressWarnings(data.frame(findOverlaps(gwgr, dmgr)))
qdf$uid = gwdf$uid[qdf$queryHits]
# From qdf: non-unique + unique counts:
# dim(qdf)
# dim(unique(qdf[,c(2,3)]))
# qdf = aggregate(queryHits ~ subjectHits + uid, )
qdf$uid = as.character(qdf$uid)
qdf$uid = factor(qdf$uid, levels=kept.uids)
qdf$fuid = as.numeric(qdf$uid)
# qnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(unique(x))})
qallnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(x)})
# Evaluate the intersections (pre-hypergeom test)
tmpdf = df[,c('id','cls')]
names(tmpdf) = c('subjectHits', 'node')
tmpdf$node = tmpdf$node + 1 # From 0-indexed --> 1-indexed
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
    mmat = qmat %*% tmp.mat 
    mmat = data.frame(as.matrix(mmat))
    mmat$uid = kept.uids
    dflist = gather(mmat, node, nsnp, -uid)
    dflist = dflist[dflist$nsnp > 0,]
    dflist$node = as.numeric(sub("^X","", dflist$node))
}
rm(df) # Don't need full df anymore
gc()

# --------------------------------------
# 1. Calculate the hyper-geometric test:
# --------------------------------------
if (tol == 2500){
    pmat.file = paste0(argv$resultdir, '/hg_raw_enr.Rda')
} else {
    pmat.file = paste0(argv$resultdir, '/hg_raw_enr_tol',tol,'.Rda')
}
if (!file.exists(pmat.file)){
    print("[STATUS] Calculating real enrichments")
    raw.pmat = matrix(NA, nrow=NUID, ncol=NF, dimnames=list(kept.uids,NULL))
    for (i in 1:NUID){
        suid = kept.uids[i]
        rawlp = try(get_pval_hg_fast(suid, dflist, lenlist, qdf, NF=NF, cutp=CUTP, olen=NENH))
        raw.pmat[suid,] = rawlp
    }
    save(raw.pmat, file=pmat.file)
} else {
    load(pmat.file)
}
print("Quantiles for enrichment values:")
print(quantile(raw.pmat, c(.98, .99, .999, .9999), na.rm=T))



# ----------------------------------------------
# 2. Run the same analysis on permuted catalogs:
# ----------------------------------------------
NPERM = 100
matlist = list()
for (j in 1:NPERM){
    if (tol == 2500){
        perm.file = paste0(argv$resultdir, '/permuted/hg_raw_perm_seed_', j, '.Rda')
    } else { 
        perm.file = paste0(argv$resultdir, '/permuted/hg_raw_perm_seed_', j, '_tol',tol,'.Rda')
    }
    if (!file.exists(perm.file)){
        t1 = proc.time()
        cat(j)
        # Shuffle UIDs:
        set.seed(j)
        shuff.uid = gwdf$uid[sample(1:TOTSNP, TOTSNP, replace=FALSE)]
        # Use same intersections object, change UID column: 
        qdf$uid = shuff.uid[qdf$queryHits]
        qdf$uid = factor(as.character(qdf$uid), levels=kept.uids)
        qdf$fuid = as.numeric(qdf$uid)
        # From qdf: non-unique + unique counts:
        qallnumdf = aggregate(queryHits ~ uid, qdf, function(x){length(x)})
        # Evaluate the intersections (pre-hypergeom test)
        cat(" - calculating intersections")
        # TODO: can we speed up this merge? Happens 100x and is very slow:
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
            mmat = qmat %*% tmp.mat 
            mmat = data.frame(as.matrix(mmat))
            mmat$uid = kept.uids
            dflist = gather(mmat, node, nsnp, -uid)
            dflist = dflist[dflist$nsnp > 0,]
            dflist$node = as.numeric(sub("^X","", dflist$node))
        }
        t2 = proc.time() 
        # print(t2- t1)
        # Calculate the hyper-geometric test:
        cat(" - calculating enrichments")
        perm.pmat = matrix(NA, nrow=NUID, ncol=NF, dimnames=list(kept.uids,NULL))
        for (i in 1:NUID){
            suid = kept.uids[i]
            rawlp = try(get_pval_hg_fast(suid, dflist, lenlist, qdf, NF=NF, cutp=CUTP, olen=NENH))
            perm.pmat[suid,] = rawlp
        }
        # TODO: Save file
        # print(proc.time() - t2)
        cat("- saving matrix\n")
        save(perm.pmat, file=perm.file)
        print("Quantiles for enrichment values:")
        print(quantile(perm.pmat, c(.98, .99, .999, .9999), na.rm=T))
    } else {
        load(perm.file)
    }
    matlist[[j]] = perm.pmat
}

# sapply(matlist, max)

# REMOVE in both observed + permutation - 0 overlap areas

# FWER by GWAS:
raw.pmat[is.na(raw.pmat)] = 0
minmat = sapply(matlist, function(x){apply(x,1, max)})
pcut = apply(minmat, 1, function(x){quantile(x, .9, na.rm=T)})
zmat = sweep(raw.pmat, 1, pcut, '/')
zmat[zmat < 1] = 0
zmat[zmat > 1] = 1
sum(apply(zmat, 1, sum) > 0)
sum(zmat)

rzmat = zmat * raw.pmat
rzmat[grep('Coronary artery',rownames(rzmat)),] -> tab

smat = rzmat[apply(rzmat, 1, sum) > 0,]
stop = apply(smat, 1, max)
head(sort(stop))

# As per Jason's method:
t = 150
t = 75
mnavg = mean(sapply(matlist, function(x){sum(x > t,na.rm=T)}))
nreal = sum(raw.pmat > t, na.rm=T)
mnavg / nreal
ap = apply(raw.pmat > t,1,sum)
# 13 enriched with tol = 0, t = 14.5 (~ 1.1% FDR)
# 18 enriched with tol = 2500, t = 150 (~ 1.1% FDR)
sum(ap > 0, na.rm=T) 
# ap[ap > 0]


# --------------------------------------
# 3. Adjust the raw enrichment p-values:
# --------------------------------------
all.pmat = do.call(rbind, matlist)
rm(matlist)
gc()

qtlist = c(.98, .99, .999, .9999)
cutall = quantile(all.pmat, qtlist, na.rm=T)
cutper = apply(all.pmat, 2, function(x){quantile(x, qtlist, na.rm=T)})

# One cutoff across all:
countdf = c()
calist = list()
for (i in 1:length(cutall)){
    ct = cutall[i]
    qt = names(cutall)[i]
    cat(paste0(qt, ": ", round(ct, 3), '\n')) 
    tmat = 1 * (raw.pmat > ct)
    ng = sum(apply(tmat, 1, sum) > 0)
    nn = sum(apply(tmat, 2, sum) > 0)
    nt = sum(tmat)
    cat(paste(ng, "gwas\n"))
    cat(paste(nn, 'nodes\n'))
    cat(paste(nt, 'total enrichments\n'))
    countdf = rbind(countdf, 
                    data.frame(qt=qt, ngwas=ng, nnodes=nn,
                               ntotal=nt, cutoff=ct, ct.type='Single cutoff'))
    calist[[qt]] = tmat * raw.pmat
}

# One cutoff per node:
cplist = list()
for (i in 1:nrow(cutper)){
    qt = names(cutall)[i]
    cat(paste0(qt, '\n')) 
    fper = cutper[i,]
    fper[fper < 4] = 4
    tmat = raw.pmat
    tmat = 1 * (sweep(tmat, 2, fper, '/') > 1)
    tmat[is.na(tmat)] = 0
    ng = sum(apply(tmat, 1, sum) > 0)
    nn = sum(apply(tmat, 2, sum) > 0)
    nt = sum(tmat)
    cat(paste(ng, "gwas"))
    cat(paste(nn, 'nodes'))
    cat(paste(nt, 'total enrichments'))
    countdf = rbind(countdf, 
                    data.frame(qt=qt, ngwas=ng, nnodes=nn,
                               ntotal=nt, cutoff=mean(fper), ct.type='Per node'))
    cplist[[qt]] = tmat * raw.pmat
}


# --------------------------
# Plot the basic statistics:
# --------------------------
clong = gather(countdf, metric, count, ngwas, nnodes, ntotal)
cl2 = aggregate(count ~ metric, clong, max)
names(cl2)[2] = 'max'
cl2$space = cl2$max * 0.06
clong = merge(clong, cl2)
cmap = data.frame(metric=c('ngwas','nnodes','ntotal'),
                  metric.nam= c('# GWAS', '# Sets', 'Total'))
clong = merge(clong, cmap)

qt.cols = colorRampPalette(brewer.pal(n=9,name="Blues"))(6)[3:6]
names(qt.cols) = as.character(unique(clong$qt))
totdf = data.frame(metric=c('nnodes','ngwas', 'ntotal'), 
                   tot=c(NF, NUID, 0))
totdf = merge(totdf, cmap)

gplot = ggplot(clong, aes(qt, count, fill=qt)) + 
    facet_grid(metric.nam ~ ct.type, scales='free') + 
    geom_bar(stat='identity', color=NA) + 
    geom_hline(data=totdf, aes(yintercept=tot), color='indianred',lty='dashed', lwd=.5) + 
    geom_text(aes(x=qt, y=count + space, label=count), cex=3, color='grey25') + 
    geom_blank(aes(qt, max + space * 2)) +  # For limits
    labs(x='FDR cutoff',y='# with sig. enrichment') + 
    scale_fill_manual(values=qt.cols) + 
    scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + 
    theme(legend.position='none')
if (argv$tagline != "") {gplot = gplot + labs(title=argv$tagline)}
ggsave(paste0(argv$resultdir, 'basic_enrich_stats.png'), gplot, 
       dpi=350, width=5, height=5, units='in')
print(paste("[STATUS] Basic stats in:", paste0(argv$resultdir, 'basic_enrich_stats.png')))

# -------------------------------
# Save the main matrices + stats:
# -------------------------------
write.table(clong, file=paste0(argv$resultdir, '/enr_stats_table.tsv'),
            quote=F, row.names=F, sep="\t")

data.file = paste0(argv$resultdir, '/hg_processed_enr.Rda')
save(cplist, calist, clong, file=data.file)

# ------------------------------------------
# Turn all results into one TSV + write out:
# ------------------------------------------
tsv.file = paste0(argv$resultdir, '/hg_processed_enr_long.tsv.gz')
enr.df = c()
for (fdr in names(cplist)){
    gmat = data.frame(cplist[[fdr]])
    colnames(gmat) = paste0('c', clslist)
    gmat$uid = rownames(gmat)
    glong = gather(gmat, cls, pvalue, -uid)
    glong = glong[glong$pvalue != 0,]
    glong$ct.type = 'Per'
    glong$qt = fdr
    enr.df = rbind(enr.df, glong)
}

for (fdr in names(calist)){
    gmat = data.frame(calist[[fdr]])
    colnames(gmat) = paste0('c', clslist)
    gmat$uid = rownames(gmat)
    glong = gather(gmat, cls, pvalue, -uid)
    glong = glong[glong$pvalue != 0,]
    glong$ct.type = 'All'
    glong$qt = fdr
    enr.df = rbind(enr.df, glong)
}

names(enr.df) = c('pmid','cluster','pvalue','ct.type','qt')
enr.df$trait = sapply(as.character(enr.df$pmid), function(x){ sub("^[0-9]* - ", "", x)})
write.table(enr.df, file=gzfile(tsv.file), quote=F, sep="\t", row.names=F)


print(paste("[STATUS] Finished calculating all uids for file", argv$bedfile))
