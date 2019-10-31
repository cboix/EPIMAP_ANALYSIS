#!/usr/bin/R
# ---------------------------------------------------
# Preliminary analysis for tree GWAS - 
# Counts on a fixed matrix from epigenomic similarity
# ---------------------------------------------------
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
library(viridis)
library(ggrepel)
library(dendextend)
library(circlize)
library(ComplexHeatmap)
library(gridBase)
library(GenomicRanges)
library(ape)
library(dplyr)
options(scipen=45) # So we dont get issues writing integers into bedfiles
# options(repos='http://cran.rstudio.com/')

# Load specific distance matrices:
source(paste0(bindir, 'load_metadata.R'))
setprefix = 'gwas_'

today <- format(Sys.time(), "%m%d%y")
imgdir = paste0(img, "gwas_tree_analysis/")
cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, setprefix)


# ----------------------------------
# Load gwas matrix for tree creation
# ----------------------------------
gwasfile = 'observed_aux_18_on_mixed_impobs_QCUT_ENH_bin_on_mixed_impobs_5000_enrich.tsv'
filepref = 'cls_merge2_wH3K27ac100_raw'
gwlindf = read.delim(gwasfile, header=F)
names(gwlindf) = c('pvalue','cluster','pmid','trait',
              'counthit','countall','fold')
namesfile = paste0(filepref, '_names.tsv')
epinames = scan(namesfile,'c')
gwlindf$pmt = paste0(gwlindf$pmid, '_', gwlindf$trait)
gwlindf$cls = paste0('c', gwlindf$cluster)
gwlindf$logpval = -log10(gwlindf$pvalue)
gwlong = aggregate(logpval ~ cls + pmt, gwlindf, max)
wide = spread(gwlong, pmt, logpval, fill=0)
gwmat = as.matrix(wide[,-1])
rownames(gwmat) = wide$cls
gwmat[gwmat < 1] = 0

# Threshold for plotting:
zmax=12
zmin=2
gwmat[gwmat > zmax] <- zmax
gwmat[gwmat < zmin] <- 0

clsn = paste0('c',1:length(epinames) - 1)
names(epinames) = clsn
epimat = gwmat
rownames(epimat) = epinames[rownames(gwmat)]


# -----------------------
# Plotting the gwas tree:
# -----------------------
# Keep only gwas with at least X hit:
uids = colnames(epimat)
# Reduce to maximal: 
cutoff = 3
keptu = uids[which(apply(epimat, 2, max) > cutoff)]
repimat = epimat[,keptu]
# Aggregate to trait level:
rdf = data.frame(repimat)
colnames(rdf) = colnames(repimat)
rdf$id = rownames(repimat)
rlong = gather(rdf, uid, value, -id)
rlong$trait = sub("^[0-9]*_","", rlong$uid)
rrlong = aggregate(value ~ id+ trait, rlong, mean)
# Make matrix again:
rrwide = spread(rrlong, trait, value)
remat = as.matrix(rrwide[,-1])
rownames(remat) = rrwide$id

# Cluster matrix:
metric='euclidean'
# metric='jaccard'
dt = dist(t(remat), method=metric)
if (metric == 'jaccard') { dt = dist(t(remat) > 0, method=metric) }
# method = 'complete'
method = 'ward.D'
ht <- hclust(dt, method=method)
ht$order <- order.optimal(dt, ht$merge)$order

# Make tree from clustering:
dend = as.dendrogram(ht)
lab = labels(dend)
head(lab)
NLAB = length(lab)
NCLUST=20
colpair = colorRampPalette(brewer.pal(n=12,name="Paired"))(NCLUST)
dend <- color_branches(dend, k=NCLUST, col=colpair)
dend = set(dend, "labels_cex", .18)

pdf(paste0(imgpref,"tree_basic.pdf"), width=14.5, height=12, onefile=T)
# plot.new()
# Plot figure:
circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
circos.initialize(factors = "single", xlim = c(0, NLAB)) 
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
                 circos.text(1:NLAB-0.5, rep(0, NLAB), labels(dend), col = labels_colors(dend), cex=.25,
                             facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) }, bg.border = NA, track.height = 0.1)
max_height = max(attr(dend, "height"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                 circos.dendrogram(dend, max_height = max_height)
                 }, track.height = 0.5, bg.border = NA)
circos.clear()
circos.clear()
dev.off()


# ----------------------------------------------------------------------------
# Do hypergeometric (or permutation-based?) enrichment of gwas for annotation:
# ----------------------------------------------------------------------------
# Run hypergeometric for each comb, make table:
# - Aka which clusters are more adult vs. which top (1-2) CT
# For each cluster center, evaluate against each metadata column
# ---------------------------------------------------
run.hyper <- function(y){phyper(q=y[1]-1, m=y[3], n=y[4] - y[3], k=y[2], lower.tail=FALSE)} 
run.fisher <- function(y){
    table <-  matrix(as.numeric(c(y[1], y[2], y[3], y[4])), ncol = 2, byrow = TRUE)
    if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
    return(p)
}

# Set up the kept elements
threshold = zmin
hdf = filter(rrlong, value > threshold)
covariates = c('lifestage','sex','type','Project','GROUP')
hdf = merge(hdf, meta[,c('id',covariates)]) 

# For each cluster (relabel) - test
keptid = unique(hdf$id)
submeta = meta[keptid,]

hgdf = NULL
enrbreaks = c(0)
for (cov in covariates){
    lvls = unique(hdf[[cov]])
    # lvls = levels(hdf[[cov]])
    cat(paste0('COVARIATE: ', cov, 
               '\nLEVELS: ', paste0(lvls, collapse=", "), "\n"))
    enrbreaks = c(enrbreaks, enrbreaks[length(enrbreaks)] + length(lvls))
    # Do multiple fishers or a hypergeom?
    # TEST: For single cluster, which covariate levels are more enriched - 1 vs. 0 
    # For each level, test:
    for (lv in lvls){
        draws = table(hdf[, 'trait'])
        subhits = table(hdf[hdf[[cov]] == lv, 'trait'])
        hits = draws * 0
        hits[names(subhits)] = subhits
        df = cbind(q=hits, draw=draws,
                   m=nrow(submeta[submeta[[cov]] == lv,]),
                   N=nrow(submeta))
        pout <- apply(df, 1, run.hyper)
        df = data.frame(df, p.value=pout, trait=rownames(df), covariate=cov, level=lv)
        rownames(df) = NULL
        if (is.null(hgdf)){hgdf = df} else {hgdf = rbind(hgdf, df)}
    }
}

# Turn into wide matrix and plot:
hgdf$log10p = -log10(hgdf$p.value)
hgdf$log10p[hgdf$p.value < 10^-12] <- 12  # Cap inf
hgdf$cl = paste0(hgdf$covariate, '\n', hgdf$level)
hgwide = spread(hgdf[, c('trait','covariate','cl','log10p')], trait, log10p)
hgmat = as.matrix(hgwide[,-c(1,2)])
CUTOFF=5
hgmat[hgmat > CUTOFF] <- CUTOFF
hgmat[hgmat < 2] <- 0
rownames(hgmat) = hgwide$cl

# Reduce the plot to top 2-3 cells + top of others.
kept.traits = sort(unique(rrlong$trait))
NT = length(kept.traits)
enrichmat = matrix("unknown/mixed", nrow=NT, ncol=ncol(metamat),
                   dimnames = list(kept.traits, colnames(metamat)))

for (cov in covariates){
    # Min p.value per cluster:
    sdf = filter(hgdf, covariate == cov)
    mdf = aggregate(p.value ~ trait, sdf, min)
    mdf = merge(mdf, sdf)
    # Remove all log10p lower than 2
    mdf = filter(mdf, log10p > 2)
    # Reduce to one per (tie breaking by lvl)
    mdf = aggregate(level ~ trait,  mdf, function(x){head(x,1)})
    covname = tolower(cov)
    enrichmat[as.character(mdf$trait), covname] = as.character(mdf$level)
}



# ----------------------------------------
# Plot the gwas tree with the enrichments:
# ----------------------------------------
mll = meta.image(enrichmat[lab,5:1], colvals=colvals, cex=.6, horiz=TRUE, return.mat=TRUE)
nummat = mll[[1]]
colsmeta = mll[[2]]
colmat = matrix(colsmeta[nummat], nrow=nrow(nummat))
labels_dend <- labels(dend)
if (as.logical(anyDuplicated(labels_dend))) {
    labels(dend) <- paste0(seq_along(labels_dend), "_", labels_dend)
    labels_dend <- labels(dend) }

dend = set(dend, "labels_cex", .75)

# pdf(paste0(imgpref,"tree_annotated.pdf"), width=14.5, height=12, onefile=T)
pdf(paste0(imgpref,"tree_annotated.pdf"), width=14.25, height=11, onefile=T)
# Plot figure:
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
    just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
# FIGURE:
circos.par(cell.padding = c(0, 0, 0, 0), track.margin=rep(0.001,2))
circos.initialize(factors = "single", xlim = c(0, NLAB)) 
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
                 circos.text(1:NLAB-0.5, rep(0, NLAB), labels(dend), col = labels_colors(dend), cex=.5,
                             facing = "clockwise", niceFacing=TRUE, adj = c(0, 0.5)) }, bg.border = NA, track.height = 0.1)
circos.track(ylim = c(0, 5), bg.border = NA, panel.fun = function(x, y) {
                 m = t(colmat[,5:1])
                 nr = nrow(m)
                 nc = ncol(m)
                 for(i in 1:nr) {
                     circos.rect(1:nc - 1, rep(nr - i, nc), 
                                 1:nc, rep(nr - i + 1, nc), 
                                 border = NA, col = m[i, ])
                 } }, track.height=0.125) 
max_height = max(attr(dend, "height"))
circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
                 circos.dendrogram(dend, max_height = max_height)
                 }, track.height = 0.5, bg.border = NA)
circos.clear()
circos.clear()
upViewport()
draw(pd.legend, x = circle_size, just = "left")
dev.off()


