#!/usr/bin/R
# ------------------------------------------------
# Load, process, normalize, and write out RNA-seq:
# Make PCA plot
# Write log2FPKM matrix
# Write QN log2FPKM matrix
# - for matrices, write: PC genes, and all genes.
# Make RNA-seq modules with KCCA and write out
# ------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
library(dplyr)
library(preprocessCore)  # For QN
# library(ggplot2)
# library(viridis)
# library(ggrepel)

# Locations:
imgdir = paste0(img, "RNAseq/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'RNAseq_qc_')

# -------------------
# Read in data files:
# -------------------
rnadir = 'RNA-seq/files/RNA-seq/tsv/'
colfile = paste0(rnadir, 'gene_columns.tsv')
lenfile = paste0(rnadir, 'gene_lengths.tsv')
exprfile = paste0(rnadir, 'all_expression.tsv.gz')
genes = scan(colfile, 'c')
exprdf = read.delim(gzfile(exprfile), sep="\t", header=F)
names(exprdf) = c('fpkm', 'dataset')

# Lengths:
lendf = read.delim(lenfile, sep="\t", header=F, stringsAsFactors=F)
names(lendf) = c('gene','length')
lendf = lendf[lendf$gene %in% genes, ]
lendf$length = as.numeric(lendf$length)

# Add genes column:
exprdf$gene = rep(genes, ND)
datasets = unique(exprdf$dataset)
ND = length(datasets)

# Metadata:
bssid = sapply(colnames(mat), function(x){sub("_.*", "",x)})
bcol = meta[bssid,'COLOR']

# ------------------
# QC, normalization:
# ------------------
# Pivot table and make matrix:
mat = pivot.tomatrix(exprdf, 'dataset', 'fpkm')

# NOTE: the huge ones are raw counts, not FPKM. Convert to FPKM
# (Fragments per kilobase of transcript per million reads)
# If no fractional counts, is feature counts:
resid = mat %% 1
mmarg = apply(resid, 2, sum)
cid = which(mmarg == 0)  # 63 datasets
mreads = apply(mat, 2, sum)
# Turn into fpkm:
m2 = sweep(mat[,cid], 2, mreads[cid] / 1e6, '/')
m2 = sweep(m2, 1, lendf$length / 1e3, '/')
mat[,cid] = m2

# Run PCA on the log2 FPKM matrix:
lmat = log2(mat + 1)
pc = prcomp(t(lmat), rank=10)
dim(pc$x)

print(summary(pc))

# Plot PCA diagnostics - looks good - almost no batch effect or otherwise.
pdf(paste0(imgpref,'PCA.pdf'),width=14.5,height=12, onefile=T)
ve = pc$sdev^2 / sum(pc$sdev^2)
plot(1:10, ve[1:10] * 100, type='b', xlab='PC', ylab='% variance explained', ylim=c(0,100))
plot(1:50, cumsum(ve[1:50] * 100), type='b', xlab='PC', ylab='Cummulative % variance explained', ylim=c(0,100))
plot(pc$x[, 1], pc$x[, 2], col = bcol, main = "PCA", xlab = "PC1", ylab = "PC2", pch=19)
plot(pc$x[, 1], pc$x[, 3], col = bcol, main = "PCA", xlab = "PC1", ylab = "PC3", pch=19)
plot(pc$x[, 2], pc$x[, 3], col = bcol, main = "PCA", xlab = "PC2", ylab = "PC3", pch=19)
dev.off()

# ---------------------------------------------
# Average replicates and write log2fpkm matrix:
# NOTE: averaging log values to reduce strong outliers
# ---------------------------------------------
lwide = data.frame(lmat)
lwide$gene = rownames(lmat)
ldf = gather(lwide, dataset, log2fpkm, -gene)
idmap = data.frame(id=bssid, dataset=colnames(mat))
ldf = merge(ldf, idmap)
ldf = aggregate(log2fpkm ~ id + gene, ldf, mean)

# Write this matrix out:
lmatfile = paste0(rnadir, 'merged_log2fpkm.mtx.gz')
write.table(ldf, file=gzfile(lmatfile), quote=F, row.names=F, sep="\t")

# Reduce to just the protein coding genes + write:
pcgenes = scan('Annotation/pc_genelist.tsv', 'c')
pcldf = ldf[ldf$gene %in% pcgenes,]
lmatpcfile = paste0(rnadir, 'merged_log2fpkm.pc.mtx.gz')
write.table(pcldf, file=gzfile(lmatpcfile), quote=F, row.names=F, sep="\t")

# ---------------------------------------------
# Quantile normalize, average, and write out:
# ---------------------------------------------
qmat = normalize.quantiles(mat)
rownames(qmat) = rownames(mat)
colnames(qmat) = colnames(mat)
qlmat = log2(qmat + 1)
pc = prcomp(t(qlmat), rank=10)
dim(pc$x)

print(summary(pc))

# Plot PCA diagnostics - looks good - almost no batch effect or otherwise.
pdf(paste0(imgpref,'qn_PCA.pdf'),width=14.5,height=12, onefile=T)
ve = pc$sdev^2 / sum(pc$sdev^2)
plot(1:10, ve[1:10] * 100, type='b', xlab='PC', ylab='% variance explained', ylim=c(0,100))
plot(1:50, cumsum(ve[1:50] * 100), type='b', xlab='PC', ylab='Cummulative % variance explained', ylim=c(0,100))
plot(pc$x[, 1], pc$x[, 2], col = bcol, main = "PCA", xlab = "PC1", ylab = "PC2", pch=19)
plot(pc$x[, 1], pc$x[, 3], col = bcol, main = "PCA", xlab = "PC1", ylab = "PC3", pch=19)
plot(pc$x[, 2], pc$x[, 3], col = bcol, main = "PCA", xlab = "PC2", ylab = "PC3", pch=19)
dev.off()

# Average and write this matrix out:
qlwide = data.frame(qlmat)
qlwide$gene = rownames(qlmat)
qldf = gather(qlwide, dataset, log2fpkm, -gene)
qldf = merge(qldf, idmap)
qldf = aggregate(log2fpkm ~ id + gene, qldf, mean)

# Write this matrix out:
qlmatfile = paste0(rnadir, 'merged_qn_log2fpkm.mtx.gz')
write.table(qldf, file=gzfile(qlmatfile), quote=F, row.names=F, sep="\t")

# Reduce to just the protein coding genes + write:
pcqldf = qldf[qldf$gene %in% pcgenes,]
qlmatpcfile = paste0(rnadir, 'merged_qn_log2fpkm.pc.mtx.gz')
write.table(pcqldf, file=gzfile(qlmatpcfile), quote=F, row.names=F, sep="\t")

# Remove all previous matrices::
rm(mat, lmat, qmat, ldf, exprdf, qldf, pcqldf)

# --------------------------
# Compute gene modules
# - On protein coding genes
# - Using log2fpkm, no qnorm
# --------------------------
lmat = pivot.tomatrix(pcldf, 'id', 'log2fpkm')
 
library(flexclust)

NCENTERS=300
fam = 'kmedians'
cl1 <- kcca(lmat, k=NCENTERS, family=kccaFamily(fam))


# Get kcca cluster attributes:
acl <- attributes(cl1)
dcl <- dist(acl$centers,'ejaccard')
hcl <- hclust(dcl)

names(acl)

# Order similar clusters together:
cocl <- order.optimal(dcl, hcl$merge)$order
cls <- ordered(acl$cluster,levels=cocl)
ord <- order(cls) # will order according to factor

# Number of regions (not bp size):
r <- range(acl$cluster)
c.uniq <- r[1]:r[2]
sizes <- sapply(c.uniq,function(x){sum(cls == x)})

# Write out centers:
write.table(cbind(acl$centers,sizes),file=paste0('BySTATE/Centers_KCCA_',fam,'_N_',NCENTERS,'_state_',state,MOTIFS,'.tsv'),sep='\t',row.names=FALSE,quote=F)




