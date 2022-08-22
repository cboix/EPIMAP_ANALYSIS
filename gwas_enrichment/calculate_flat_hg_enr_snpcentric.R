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
options(width=170)

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
argv = list()
argv$uqfile = 'gwas_hg_stats/rawdata/masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt'
# argv$bedfile = 'gwas_hg_stats/rawdata/observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_full.bed.gz'
# argv$outpref = 'gwas_results_hg_snpcentric'
argv$bedfile = 'gwas_hg_stats/rawdata/observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_full_CLUSTERS.bed.gz'
argv$outpref = 'gwas_results_hg_snpcentric_clusters'
argv$resultdir = 'gwas_hg_stats/rawdata/'
argv$tagline=''

# # Create argument parser:
# p <- arg_parser("Calculate GWAS enrichments for sets of regions")
# p <- add_argument(p, "--bedfile", help="Region assignments (with chr start end cls OR identifer cls)")
# p <- add_argument(p, "--resultdir", help="Output directory")
# p <- add_argument(p, "--uqfile", help="List of unique regions + identifier", default="")
# p <- add_argument(p, "--tagline", help="tag line", default="")
# # Parse the command line arguments:
# argv <- parse_args(p)

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

# For pruning SNPs:
prune.snps = function(suid, df=gwdf, dist=5e3, quiet=TRUE){
    subdf = df[df$uid == suid,]
    keptdf = c()
    chrlist = c(as.character(1:22), 'X')
    keptlist = sapply(rep(0, 23), function(x){c()})
    names(keptlist) = chrlist
    for (i in 1:nrow(subdf)){
        chrom = as.character(subdf$chrom[i])
        loc = subdf$chromStart[i]
        # Check if any within range/add:
        if (is.null(keptlist[[chrom]])){
            keptlist[[chrom]] = loc
        } else {
            nclose = sum(abs(loc - keptlist[[chrom]]) < dist)
            if (nclose ==  0){
                keptlist[[chrom]] = c(keptlist[[chrom]], loc)
            }
        }
    }
    # Report out - chrom, loc, uid is sufficient:
    kdf = c()
    for (chrom in chrlist){
        if (!is.null(keptlist[[chrom]])){
            kdf = rbind(kdf, data.frame(chrom=chrom, 
                                        chromStart=keptlist[[chrom]],
                                        uid=suid))
        }
    }
    if (!quiet){ print(paste(suid, ': Kept', nrow(kdf), 'of',nrow(subdf))) }
    return(kdf)
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

run.pruning = FALSE
if (run.pruning){
    # Re-do pruning, this time to 1Mb:
    ut = as.character(unique(gwdf$uid))
    kept.snps = ldply(ut, df=gwdf, dist=1e6, quiet=FALSE, prune.snps)
    print(dim(gwdf))
    gwdf = merge(kept.snps, gwdf)
    print(dim(gwdf))
}

# Subset to high-quality GWAS:
subsetgwas = TRUE
if (subsetgwas){
    NIND = 10000
    NSNP = 10
    nsnpdf = aggregate(pValue ~ uid, gwdf, length)
    gwssdf = aggregate(sampsize ~ uid, gwssdf, max)
    gwssdf = merge(gwssdf, nsnpdf)
    keptgw = gwssdf[gwssdf$sampsize >= NIND & gwssdf$pValue >= NSNP,]
    kept.uids =  sort(unique(as.character(keptgw$uid)))
} else {
    kept.uids =  sort(unique(as.character(gwdf$uid)))
}
NUID = length(kept.uids)
print(paste("Number of kept uids:", NUID))
gwdf = gwdf[gwdf$uid %in% kept.uids,]
gwgr = GRanges(paste0('chr', gwdf$chrom), IRanges(gwdf$chromStart, gwdf$chromEnd))

TOTSNP = nrow(gwdf)
NCAT = 1000
print(paste0("Pruned to: ", TOTSNP, " snps in ", NUID, " GWAS."))

# nrow(unique(gwdf[,c('chrom', 'chromStart')]))

# ---------------------------------------------
# Load in the query df + process and intersect:
# ---------------------------------------------
t1 = proc.time()
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
    df$loc = NULL
    rm(indmap)
}
dmgr = GRanges(uqdf$chr, IRanges(uqdf$start - tol, uqdf$end + tol))
# rm(uqdf)
gc()
elt = proc.time() - t1
cat('Loading element sets:', paste0(round(elt[3], 2), 's elapsed\n'))

# -------------------------------------------
# Calculate enhancer intersections with GWAS:
# -------------------------------------------
t1 = proc.time()
qdf = suppressWarnings(data.frame(findOverlaps(gwgr, dmgr)))
qdf$uid = gwdf$uid[qdf$queryHits]
names(qdf)[2] = 'id'

# Reduce to number of SNPs per trait per biosample:
int.enh = sort(unique(qdf$id))
subdf = df[df$id %in% int.enh,]
# rm(df)
gcout = gc()

NSNP = nrow(gwdf)
NE = max(int.enh)
NC = max(df$cls) + 1

# Matrix: SNPs to Enhancers:
qmat = sparseMatrix(i = qdf[,1], j = qdf[,2], x = rep(1, nrow(qdf)), dims = c(NSNP,NE))
# Matrix: Sets to Enhancers.
smat = sparseMatrix(i = subdf[,2], j = subdf[,1] + 1, x = rep(1, nrow(subdf)), dims = c(NE, NC))
# Matrix: SNPs to UIDs.
NUID = length(kept.uids)
gwdf$uid = factor(gwdf$uid, levels=kept.uids)
umat = sparseMatrix(i = as.numeric(gwdf$uid), j = 1:nrow(gwdf), x = rep(1, nrow(gwdf)), dims = c(NUID, NSNP))

# Calculate overlaps of SNPs with CLS:
tmat = qmat %*% smat
tmat@x = rep(1, length(tmat@x)) # Set to unique SNPs captured only

# ------------------------------------------------
# Calculate the hypergeometric enrichment results:
# ------------------------------------------------
run.enhsnp.hyper = function(umat, tmat, seed=1, rand=FALSE){
    # Calculate numsnp in each GWAS x CLS intersection
    set.seed(seed)
    if (rand){
        # Shuffle SNPs, leaving the same # of SNPs per GWAS:
        sind = sample(1:ncol(umat), ncol(umat), replace=FALSE)
        ns.mat = umat[,sind] %*% tmat
    } else { ns.mat = umat %*% tmat }
    # Number of captured snps:
    totuid = rowSums(umat)
    totcls = colSums(tmat)
    # Data frame of number of SNPs for each calc:
    ctdf = as.data.frame(summary(ns.mat))
    ctdf$j = ctdf$j - 1
    names(ctdf) = c('i','cls','nsnp')
    # Alternatively, if we want to keep + test the 0 intersections:
    # matdf = data.frame(as.matrix(ns.mat))
    # matdf$i = 1:nrow(matdf)
    # ctdf = gather(matdf, cls, nsnp, -i)
    # ctdf$cls = as.numeric(sub("^X","",ctdf$cls)) - 1
    ctdf$uid = kept.uids[ctdf$i]
    ctdf$totcls = totcls[ctdf$cls + 1]
    ctdf$totuid = totuid[ctdf$i]
    ctdf$totsnp = length(gwgr)
    # Run SNP-centric hypergeometric test:
    pvdf = ctdf[,c('nsnp','totuid','totcls','totsnp')]
    pout <- apply(pvdf, 1, run.hyper)
    ctdf$p = pout
    ctdf$lp = -log10(pout)
    ctdf = ctdf[order(ctdf$p),]
    ctdf$padj = p.adjust(ctdf$p,'BH') # Correct by BH
    return(ctdf)
} 


# Actual results:
t1 = proc.time()
resdf = run.enhsnp.hyper(umat, tmat, seed=1, rand=FALSE)
elt = proc.time() - t1
cat(paste0(round(elt[3], 2), 's elapsed\n'))

cat("At p.adj < 0.05:", length(unique(resdf$uid[resdf$padj < 0.05])),'sig. GWAS\n') # 325
cat("At p.adj < 0.02:", length(unique(resdf$uid[resdf$padj < 0.02])),'sig. GWAS\n') # 256
cat("At p.adj < 0.01:",length(unique(resdf$uid[resdf$padj < 0.01])),'sig. GWAS\n') # 226
cat("At p.adj < 0.05:", length(resdf$uid[resdf$padj < 0.05]),'sig. GWAS\n') # 22898
cat("At p.adj < 0.02:", length(resdf$uid[resdf$padj < 0.02]),'sig. GWAS\n') # 15866
cat("At p.adj < 0.01:",length(resdf$uid[resdf$padj < 0.01]),'sig. GWAS\n') # 12494
mp5cut = max(resdf$p[resdf$padj < 0.05])
mp2cut = max(resdf$p[resdf$padj < 0.02])
mp1cut = max(resdf$p[resdf$padj < 0.01])
cat("Uncorrected p-value for p.adj < 0.05:", mp5cut, '\n')
cat("Uncorrected p-value for p.adj < 0.02:", mp2cut, '\n')
cat("Uncorrected p-value for p.adj < 0.01:", mp1cut, '\n')


# Number of loci:
siguid = unique(resdf$uid[resdf$padj < 0.02])
sum(gwdf$uid %in% siguid)
sum(gwdf$uid %in% siguid) / nrow(gwdf)



# Add the tissue to results + save results:
# -----------------------------------------
source(paste0(bindir, 'load_metadata.R'))
bssids = scan('Enhancer_matrix_names.txt', 'c', quiet=T)
if (max(resdf$cls) > 300){
    resdf$id = bssids[resdf$cls + 1]
    resdf = merge(resdf, meta[,c('id','ct','infoline')])
} else {
    resdf$id = resdf$cls + 1
}

# Also sanity check:
resdf = resdf[order(resdf$p),]
outdf = resdf[resdf$padj < 0.02,]
# outdf = outdf[, c('uid','infoline','nsnp', 'lp')]

# Sanity check on adding tissue:
head(outdf[grep('Coronary', outdf$uid),], 5)
head(outdf[grep('Alzheimer', outdf$uid),], 5)
head(outdf[grep('Cholesterol', outdf$uid),], 5)

# Save these results:
keep.cols = c('uid','id','nsnp','totcls','totuid','totsnp','p','lp','padj', 'ct','infoline')
keep.cols = keep.cols[keep.cols %in% colnames(resdf)]
outdf = resdf[,keep.cols]
write.table(outdf, file=gzfile(paste0(argv$outpref, '.tsv.gz')), quote=F, row.names=F, sep="\t")
saveRDS(outdf, file=paste0(argv$outpref, '.Rds'))



# Run permutations to evaluate FDR control:
# -----------------------------------------
NPERM = 100
mplist = c()
allplist = c()
uidplist = c()
trackuid = TRUE # If we run this, we need more than 100 permutations.
pmat = matrix(NA, nrow=NUID, ncol=NC)
for (i in 1:NPERM){
    t1 = proc.time()
    randdf = run.enhsnp.hyper(umat, tmat, seed=i, rand=TRUE)
    mplist = c(mplist, min(randdf$p)) # Collect the smallest p-value each time:
    allplist = c(allplist, randdf$p) # Collect all of the p-values:
    # For evaluating FWER per GWAS:
    if (trackuid){
        pmat = sparseMatrix(i = randdf$i, j = randdf$cls + 1, x = -log10(randdf$p), dims = c(NUID,NC))
        uid.marg = apply(pmat, 1, max)
        cls.marg = apply(pmat, 2, max)
        uidplist = rbind(uidplist, uid.marg)
    }
    nsig = length(unique(randdf$uid[randdf$padj < 0.05]))
    elt = proc.time() - t1
    cat(i, paste0(round(elt[3], 2), 's elapsed.\n'))
    pct5 = round(sum(randdf$p < mp5cut) / sum(resdf$p < mp5cut) * 100,2)
    pct1 = round(sum(randdf$p < mp1cut) / sum(resdf$p < mp1cut) * 100,2)
    cat("At the p.adj < 0.05 cutoff:", sum(randdf$p < mp5cut) ,'sig. intersections', pct5, '%\n')
    cat("At the p.adj < 0.01 cutoff:", sum(randdf$p < mp1cut) ,'sig. intersections', pct1, '%\n')
}

# Should be minimum by FDR by Jason's suggestion:
p1cut = quantile(mplist, .01, na.rm=T)
p5cut = quantile(mplist, .05, na.rm=T)

cat("int @5%:", length(resdf$uid[resdf$p < p5cut]),'\n') # 3516 intersections
cat("GWAS @5%:", length(unique(resdf$uid[resdf$p < p5cut])),'\n') # 101 sig. GWAS
cat("int @1%:", length(resdf$uid[resdf$p < p1cut]),'\n') # 3302 intersections
cat("GWAS @1%:", length(unique(resdf$uid[resdf$p < p1cut])),'\n') # 100 sig. GWAS

# This protocol is probably too strict?
(sum(allplist < p1cut) / NPERM) / sum(resdf$p < p1cut) * 100 # 0.0018%
(sum(allplist < p5cut) / NPERM) / sum(resdf$p < p5cut) * 100 # 0.0028%

# Using the p-adjust threshold:
cat("int @5%:", length(resdf$uid[resdf$padj < 0.05]),'\n') # 22898 intersections for epigenomes, 348 modules
cat("GWAS @5%:", length(unique(resdf$uid[resdf$padj < 0.05])),'\n') # 325 sig. GWAS epi, 146 modules
cat("int @2%:", length(resdf$uid[resdf$padj < 0.02]),'\n') # 15866 intersections epi, 234 modules
cat("GWAS @2%:", length(unique(resdf$uid[resdf$padj < 0.02])),'\n') # 256 sig. GWAS epi, 104 modules
cat("int @1%:", length(resdf$uid[resdf$padj < 0.01]),'\n') # 12494 intersections epi, 183 modules
cat("GWAS @1%:", length(unique(resdf$uid[resdf$padj < 0.01])),'\n') # 226 sig. GWAS epi, 93 modules

# The p-adjust threshold by BH seems like it gives a reasonable FDR,
# controlling the estimated percentage of false discoveries within the total number of significant tests:
(sum(allplist < mp1cut) / NPERM) / sum(resdf$p < mp1cut) * 100 # 0.6% for epigenomes, 1.02% for modules
(sum(allplist < mp2cut) / NPERM) / sum(resdf$p < mp2cut) * 100 # 1.27% for epigenomes, 2.15% for modules
(sum(allplist < mp5cut) / NPERM) / sum(resdf$p < mp5cut) * 100 # 3.3% for epigenomes, 5.64% for modules
# e.g. at the adjusted p-value threshold of 5%, an estimated 3.3% of the significant tests are false
# and at the adjusted p-value threshold of 1%, an estimated 0.6% of the significant tests are false

if (trackuid){
    # Finally, apply different FWER threshold for each GWAS: 
    up1cut = apply(uidplist, 2, function(x){quantile(10^-x, .01, na.rm=T)})
    up5cut = apply(uidplist, 2, function(x){quantile(10^-x, .05, na.rm=T)})
    up1cut = apply(uidplist, 2, max)
    pmat = sparseMatrix(i = resdf$i, j = resdf$cls + 1, x = -log10(resdf$p), dims = c(NUID,NC))
    z1mat = 1 * (sweep(pmat, 1, up1cut, '/') > 1)
    sum(z1mat)
    sum(rowSums(z1mat) > 0)

    # Compare to the last random catalog:
    pmat = sparseMatrix(i = randdf$i, j = randdf$cls + 1, x = -log10(randdf$p), dims = c(NUID,NC))
    z1mat = 1 * (sweep(pmat, 1, up1cut, '/') > 1)
    sum(z1mat)
    sum(rowSums(z1mat) > 0)

    # Compare to a new random catalog (not stable unless we calculate cutoffs based on ~1000+ null catalogs)
    outdf = run.enhsnp.hyper(umat, tmat, seed=-1, rand=TRUE)
    pmat = sparseMatrix(i = outdf$i, j = outdf$cls + 1, x = -log10(outdf$p), dims = c(NUID,NC))
    z1mat = 1 * (sweep(pmat, 1, up1cut, '/') > 1)
    sum(z1mat)
    sum(rowSums(z1mat) > 0)
}


# From Jason:
# q: 'white balls drawn without replacement from an urn which contains both black and white balls.' - # of trait SNPs overlapping an enhancer DHS in tissue (minus 1)
# m: 'the number of white balls in the urn.' - # of pruned GWAS catalog SNPs overlapping enhancer DHS in tissue
# n: 'the number of black balls in the urn' - # of pruned GWAS catalog SNPs not overlapping enhancer DHS in tissue
# k: 'the number of balls drawn from the urn' - # of trait SNPs
# One could still extend enhancers by 2.5kb, if doing things from a SNP perspective and thus counting an overlap at most once, though I note I did not see the 2.5kb extension documented in the paper.

