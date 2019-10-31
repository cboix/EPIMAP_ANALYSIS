#!/usr/bin/R
# -----------------------------------------------------
# Explore and annotate main examples of complex traits:
# - As part of main figure for paper
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
library(flexclust)
library(gridExtra)
library(kableExtra)  # For latex table.
library(grid)
library(clusterProfiler)  # For GO
library(org.Hs.eg.db)
library(GenomicRanges)

# ----------------------
# Arguments for modules:
# ----------------------
filepref = 'cls_merge2_wH3K27ac100_300'
tagline = 'ChromHMM Enhancers'
imgdir = paste0(img, "clusters/") 
# Load in and process data (saves to matrices):
commandArgs <- function(trailingOnly=TRUE){
    c(filepref, tagline, imgdir) }
source(paste0(bindir, 'load_modules_data.R'))


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

# ------------------
# Load RNA-seq data:
# ------------------
# Read in the RNA-seq protein-coding genes:
rnadir = 'RNA-seq/files/RNA-seq/tsv/'
lmatpcfile = paste0(rnadir, 'merged_log2fpkm.pc.mtx.gz')
availfile = paste0(rnadir, 'avail_rna.txt')
if (!file.exists(availfile)){
    pcdf = read.delim(gzfile(lmatpcfile), header=T, sep="\t", stringsAsFactors=F)
    pcid = unique(pcdf$id)
    write.table(pcid, availfile, quote=F, row.names=F, sep="\t")
} else { 
    pcid = scan(availfile, 'c', sep="\n")
}
rna.avail = 1 * (labels(dend) %in% pcid)

# (very rough) average per tissue:
avgfile = paste0(rnadir, 'avg_rna.Rda')
if (!file.exists(avgfile)){
    pcdf = merge(pcdf, meta[,c('id','GROUP')])
    red.pcdf = aggregate(log2fpkm ~ gene + GROUP, pcdf, function(x){x = 2^x; log2(mean(x))})
    rm(pcdf)
    save(red.pcdf, file=avgfile)
} else {
    load(avgfile)
}

# Load in TSS annotation:
tssfile = 'Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz'
tssdf = read.delim(gzfile(tssfile), header=F, sep="\t", stringsAsFactors=F)
names(tssdf) = c('chr','tss','t2','gene')

gmfile = 'Annotation/gencode.gene.name.mapping.tsv'
gmdf = read.delim(gmfile, header=F, sep="\t", stringsAsFactors=F)
names(gmdf) = c('gene', 'symbol')

# ----------------
# Trait arguments:
# ----------------
# Main + 5 comparison traits:
if (examp == 'CAD'){
    suid = "29212778 - Coronary artery disease" # Have: CAD_UKBIOBANK.gz
    slist = c("30595370 - Systolic blood pressure",  # Have: systolic_UKB2_sumstats.txt.gz
              "29507422 - High density lipoprotein cholesterol levels", # Can't find summ stats.
              # "24097068 - HDL cholesterol", # http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_HDL.txt.gz
              "30061737 - Atrial fibrillation", # Here:
              # http://csg.sph.umich.edu/willer/public/afib2018/nielsen-thorolfsdottir-willer-NG2018-AFib-gwas-summary-statistics.tbl.gz
              # "24952745 - QT interval",
              # "30578418 - Diastolic blood pressure", ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs001672/analyses/phs001672.pha004730.txt
              # "30578418 - Pulse pressure", # ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs001672/analyses/phs001672.pha004731.txt
              # "30578418 - Systolic blood pressure", # ftp://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs001672/analyses/phs001672.pha004732.txt
              "26426971 - Waist-to-hip ratio adjusted for BMI", # Split into m/f
              # http://portals.broadinstitute.org/collaboration/giant/images/4/48/WHRADJ.WOMEN.LE50.publicrelease.txt.gz
              # http://portals.broadinstitute.org/collaboration/giant/images/6/62/WHRADJ.WOMEN.GT50.publicrelease.txt.gz
              # http://portals.broadinstitute.org/collaboration/giant/images/6/60/WHRADJ.MEN.LE50.publicrelease.txt.gz
              # http://portals.broadinstitute.org/collaboration/giant/images/6/6f/WHRADJ.MEN.GT50.publicrelease.txt.gz
              # Alternatively, this is UKB2 meta-analysis data: 
              # https://zenodo.org/record/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz
              "30595370 - Waist-hip ratio"
              # UKB2: # https://zenodo.org/record/1251813/files/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz
    ) 
} else if (examp == 'AD'){
    # Different AD GWAS
    # Show T2D?
    # Show brain traits?
} else if (examp == 'CD'){
    suid = "28067908 - Crohn's disease"
    # CD, UC, IBD - agreement.
    slist = c("28067908 - Crohn's disease",
              "29507422 - High density lipoprotein cholesterol levels")
} else if (examp == 'MAC'){
    suid = "30535121 - Macular thickness"
    slist = c("30595370 - Systolic blood pressure",  # UKB2 data 
              "29507422 - High density lipoprotein cholesterol levels", # Can't find summ stats.
              "30061737 - Atrial fibrillation",
              "26426971 - Waist-to-hip ratio adjusted for BMI", # Split into m/f
              "30595370 - Waist-hip ratio") # UKB2 data
}

pdf(paste0(img, 'leadsnp_pvalues.pdf'), width=10, height=4)
par(xaxs='i')
par(yaxs='i')
par(mar=c(2,4,1,0))
hist(-log10(gwdf$pValue), 200, col='grey25', border=NA, main='')
dev.off()

plist = -log10(gwdf$pValue)
plist[plist > 50] = 50
pdf(paste0(img, 'leadsnp_pvalues_cut.pdf'), width=10, height=4)
par(xaxs='i')
par(yaxs='i')
par(mar=c(4,4,1,1))
hist(plist, 200, col='grey25', border=NA, main='', xlim=c(0, 50), xlab='-log10 p-value')
axis(1, at=c(seq(0,20,2), seq(50, 300,50)))
abline(v=-log10(5e-8), lwd=2, col='red', lty='dashed')
dev.off()

sum(plist < -log10(5e-8)) / length(plist)
sum(plist >= -log10(5e-8)) / length(plist)


# -- may need to move up consensus tree

# -----------------------------
# Functional evaluation figure:
# -----------------------------
# Panels:
# - top enr + pval (color + #)
# - enrichments of enhancers (start with the enrichments of the modules?)
# (enrichmat)
# - set size
# - enhancer set size (enhmargfile)
# - GO enrichments for nearest gene
# - GO enrichments for nearest expressed gene
# - Clustered SNP matrix (by all cons)
# - GO enrichments for clustered matrix

# -------------------------------------
# Load all gwas enrichment information:
# -------------------------------------
snpfile = paste0(regdir, eprefix, apref, '_logreg_all_wsnp', suffix)
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

cad.dlen = cdlenlist$diff[cad.nodes]
cad.clen = cdlenlist$cons[cad.nodes]
# Boxplot of specificity:
cad.setdf = ldply(cad.nodes, function(x){data.frame(node=x, row=cdll$diff[[x]])})
cad.setdf = merge(cad.setdf, enhsetdf)
cad.setdf = merge(cad.setdf, enhmarg)

# Set attributes:
cad.setcls = aggregate(row ~ node + cls, cad.setdf, length)
enrdf = data.frame(cbind(enrichmat, cls=rownames(enrichmat)))
enrdf = gather(enrdf, facet, set, -cls)
cad.setcls = merge(cad.setcls, enrdf)
cad.setcls$col = sapply(1:nrow(cad.setcls), function(i){colvals[[cad.setcls$facet[i]]][cad.setcls$set[i]]})

cad.setcls = cad.setcls[order(cad.setcls$row, decreasing=T),]
cad.setcls = cad.setcls[order(cad.setcls$facet),]
cad.setcls = cad.setcls[order(cad.setcls$node),]

# Get colors for the enhsets:
dist.order = plot.centers(centlist, tagline, rngroup, counts=counts, 
                          cls=acutgroup.nam, title=FALSE, ablwd=.5, calc.only=TRUE)
maxclsdf = data.frame(COLOR=c('grey', colset)[dist.order[[3]] + 1], cls=dist.order[[1]])
maxclsdf = merge(maxclsdf, odf, all.x=TRUE)
maxclsdf$GROUP = as.character(maxclsdf$GROUP)
maxclsdf$GROUP[is.na(maxclsdf$GROUP)] = 'Non-specific'
maxclsdf$GROUP = factor(maxclsdf$GROUP, levels=c(odf$GROUP, 'Non-specific'))


subdf = cad.setcls[cad.setcls$facet == 'group',]
subdf = merge(subdf, maxclsdf)
subdf = merge(subdf, data.frame(node=cad.nodes, name=nodenames))
gcol = c(colvals[['group']], 'Non-specific'='grey')
gp = ggplot(subdf, aes(name, row, fill=GROUP)) + 
    geom_bar(stat='identity', position='stack') + 
    coord_flip() + 
    scale_fill_manual(values = gcol) + 
    theme_minimal()
ggsave(paste0('~/', examp,'_module_groupprop.png'), gp, dpi=450, units='in', width=10, height=9)


ls.subdf = cad.setcls[cad.setcls$facet == 'lifestage',]
ls.gcol = c(colvals[['lifestage']], 'unknown/mixed'='grey')
gp = ggplot(ls.subdf, aes(factor(node), row, fill=set)) + 
    geom_bar(stat='identity', position='fill') + 
    scale_fill_manual(values = ls.gcol) + 
    theme_minimal()
ggsave(paste0('~/', examp,'_module_lsprop.png'), gp, dpi=450, units='in', width=10, height=9)

# For boxplots of specificity:
sub.cad.setdf = cad.setdf[cad.setdf$node %in% cad.nodes[nord],]
sub.cad.setdf$node = factor(sub.cad.setdf$node, levels = cad.nodes[nord])

# Find the nearest genes (using genomic ranges)
cad.enhloc = cad.setdf[,c('row','node')]
names(cad.enhloc) = c('cls','node')
cad.enhloc = merge(cad.enhloc, enhdf)
enhgr = GRanges(cad.enhloc$chr, IRanges(cad.enhloc$start, cad.enhloc$end))
tssgr = GRanges(tssdf$chr, IRanges(tssdf$tss, tssdf$tss))
ndf = c(nearest(enhgr, tssgr))
cad.enhloc$gene = tssdf$gene[ndf]

write.table(cad.enhloc[,c('chr','start','end','node')], 
            file=paste0(examp, '_top', NPLOT,'_nodes_location.bed'),
            row.names=F, quote=F, sep="\t", col.names=F)
cad.enhloc = cad.enhloc[,c('cls','node','gene')]
cad.enhloc = merge(cad.enhloc, gmdf)

# Perform over-representation test enrichment:
cad.gofile = paste0(regpref, examp,'_go_nodediff_enh.Rda')
if (!file.exists(cad.gofile)){
    cad.godf = c()
    for (ni in cad.nodes){
        x = unique(cad.enhloc$symbol[cad.enhloc$node == ni])
        gdf = go.enr(x, allx=gmdf$symbol)
        gdf$node = ni
        cad.godf = rbind(cad.godf, gdf)
    }
    save(cad.godf, file=cad.gofile)
} else { 
    load(cad.gofile)
}

# Prune the obvious terms?
sub.godf = cad.godf[, c('node','Description', 'qvalue')]
sub.godf = sub.godf[sub.godf$qvalue < 1e-8,]
ndesc = aggregate(node ~ Description, sub.godf, length)
ndesc = ndesc[order(ndesc$node, decreasing=T), ] 
# go.cutoff = NPLOT / 2
go.cutoff = 5
go.desc = ndesc$Description[ndesc$node <= go.cutoff]
# Keep 1-3 terms per term, make dataframe with all incidences:
sub.godf = sub.godf[sub.godf$Description %in% go.desc,]
sub.godf = sub.godf[order(sub.godf$qvalue),]
go.keep = c()
for(ni in cad.nodes){ go.keep = c(go.keep, head(sub.godf$Description[sub.godf$node == ni], 2)) }
go.keep = unique(go.keep)

# Get order for GO terms:
ssdf = sub.godf[sub.godf$Description %in% go.keep,]
ssdf$log10q = -log10(ssdf$qvalue)
sswide = spread(ssdf[,c('Description','node','log10q')], node, log10q, fill=0)
mat = as.matrix(sswide[,-1])
rownames(mat) = sswide$Description
mat = mat[, as.character(rev(cad.nodes[nord]))]
rd = diag.mat2(t(mat) > 0)
go.ord = rd[[2]]
# Make nicer plot:
ssdf$node = factor(ssdf$node, levels=cad.nodes[nord])
ssdf$Description = factor(ssdf$Description, levels=rev(go.ord))
lindf = rbind(data.frame(x0=1:length(go.ord), x1=1:length(go.ord), y0=1, y1=NPLOT), 
              data.frame(x0=1, x1=length(go.ord), y0=1:NPLOT, y1=1:NPLOT))
ylindf = data.frame(x0=1:length(go.ord), x1=1:length(go.ord), y0=1, y1=NPLOT)

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

# -------------------------------
# Display shared top enrichments:
# -------------------------------
shmat = rmat[keptgwas, cad.nodes[nord]]
shmat = shmat[apply(shmat, 1, max) > 3,]
shmat2 = shmat[apply(shmat > 0, 1, sum) > 1,]
shmat2 = t(diag.mat2(t(shmat2))[[1]])

# Fix row:
repl = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
              "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
rownames(shmat2)[rownames(shmat2) == repl] = "26974007 - Chronic inflammatory diseases (pleiotropy)"

# Add legend
enr.col_fun = function(x, pal=hcl.colors(12, "YlOrRd", rev=TRUE)){
    palette = pal
    bin <- cut(x, seq(3, 50, length.out=length(palette)+1), include.lowest=T, right=TRUE) 
    palette[bin] }

enr.legend = Legend(at=round(seq(3, 50, length.out=7), 2),
                     labels_gp = gpar(fontsize=3.5),
                     title_gp = gpar(fontsize=3.5, fontface='bold'),
                     col_fun=enr.col_fun, title_position = "topleft", 
                     title='-log10p', direction = 'vertical')
plegend = packLegend(enr.legend)

pdf(paste0(cadpref, 'all_compare_uid.pdf'), width=3.75, height=4)
sp=.15
bsp=13
par(mar=c(sp, bsp, 4, sp))
SR = nrow(shmat2)
SC = ncol(shmat2)
shrn = rev(1:SR)
shcn = rev(1:SC)
shmat2[shmat2 > 50] = 50
image(t(shmat2[shrn,shcn]), zlim=c(3,50), axes=F, useRaster=T)
xat = seq(0,1,length.out=NPLOT)
text(y=par()$usr[4] + .005 * diff(par()$usr[3:4]), x=xat, 
     rev(leafrep), col=rev(nodecol[nord]), srt=90, adj=0, xpd=TRUE, cex=.3)
box(lwd=.5)
xat = seq(0, 1, length.out=SC)
yat = seq(0, 1, length.out=SR)
text(x=par()$usr[1]-.01 * diff(par()$usr[1:2]), y=yat,
     rev(rownames(shmat2)), col='black', 
     # rev(sub("^[0-9]* - ", "", rownames(shmat))), col='black', 
     srt=0, adj=1, xpd=TRUE, cex=.35)
draw(plegend, x = unit(0.25,'in'), y=unit(2.05,'in'), just = "top")
dev.off()

shdf = data.frame(shmat[slist,])
cord = colnames(shdf)
rord = rownames(shdf)
shdf$uid = rownames(shdf)
shdf = gather(shdf, node, value, -uid)
shdf = shdf[shdf$value > 0,]
shdf$col = 'black'
shdf$col[shdf$value > 40] = 'white'
shdf$uid = factor(shdf$uid, levels=rord)
shdf$node = factor(shdf$node, levels=cord)
shdf$node2 = as.numeric(sub("X", "", as.character(shdf$node)))
shdf$node2 = cad.nodes[nord][shdf$node2 ]
shmat[shmat > max(cad.lp)] = max(cad.lp)
write.table(shdf[,c('uid','node2')], 
            file=paste0(examp, '_top', NPLOT,'_shared_nodes.tsv'),
            row.names=F, quote=F, sep="\t", col.names=F)


# Make the SNP matrices:
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

set.seed(1)
smat = smat[apply(smat, 1, sum) > 0,]
dnj.cons = dist(smat, 'jaccard')
# dnj.cons = dist(smat.cons, 'jaccard')
hnj.cons <- hclust(dnj.cons, method='ward.D')
cocl = order.optimal(as.dist(dnj.cons), hnj.cons$merge)$order
snp.cons.rn <- names(cocl)[cocl]
snp.cons.rn = snp.cons.rn[snp.cons.rn %in% rownames(smat)]
snp.breaks = calc.breaks(hnj.cons, 10, cocl)
snp.cto = cutree(hnj.cons, 10)[cocl]


# ----------------------
# GO terms for clusters:
# ----------------------
snpdf = data.frame(queryHits=as.numeric(sub("s","",names(snp.cto))), 
                   snp=names(snp.cto), cls=snp.cto)
snpdf$chr = paste0('chr',gwdf$chrom[snpdf$queryHits])
snpdf$start = gwdf$chromStart[snpdf$queryHits]
snpgr = GRanges(snpdf$chr, IRanges(snpdf$start, snpdf$start))
snpdf$gene = tssdf$gene[nearest(snpgr, tssgr)]
stssdf = merge(snpdf, gmdf)
# Aggregate GO test sets per cluster:
tree.setsdf = aggregate(symbol ~ cls, stssdf,
                        function(x){paste(unique(x), collapse=", ")})
mapsnpdf = aggregate(symbol ~ snp, stssdf, function(x)x[1])
mapsnp = mapsnpdf$symbol
names(mapsnp) = mapsnpdf$snp


plot_snp_partial = function(mat, numsnp=NULL){
    bsp=3
    tsp=6
    # Side 1
    xat = seq(0,1,length.out=NPLOT)
    par(mar=c(sp, bsp, tsp, sp))
    image(t(mat), col=c('white','grey25'), axes=F, useRaster=T)
    text(y=par()$usr[4] + .005 * diff(par()$usr[3:4]), x=xat, 
         rev(leafrep), col=rev(nodecol[nord]), 
         srt=90, adj=0, xpd=TRUE, cex=.5)
    box(lwd=.5)
    msnp = mapsnp[rownames(mat)]
    if (!is.null(numsnp)){ msnp = paste0(rev(numsnp) , '. ', msnp) }
    yat = seq(0, 1, length.out=length(msnp))
    text(x=par()$usr[1]-.01 * diff(par()$usr[1:2]), y=yat,
         msnp, col='black', srt=0, adj=1, xpd=TRUE, cex=.35)
}

# Plot the SNPs + closest genes as is:
pdf(paste0(cadpref, 'all_snps_top20.pdf'), width=6, height=6)
layout(matrix(c(1:3),1,3))
par(xaxs='i')
par(yaxs='i')
sp=.15
smat.reord = smat[snp.cons.rn, rev(nodenames[nord])]
smat.reord = smat.reord[apply(smat.reord,1, sum) > 0,]
midSR = round(nrow(smat.reord) / 3)
plot_snp_partial(smat.reord[1:midSR,], numsnp=1:midSR)
plot_snp_partial(smat.reord[(midSR+1):(2*midSR),], numsnp=(midSR+1):(2*midSR))
plot_snp_partial(smat.reord[(2*midSR+1):nrow(smat.reord),], numsnp=(2*midSR+1):nrow(smat.reord))
dev.off()


# -----------------------------------
# Perform GO enrichment for clusters:
# -----------------------------------
x = tree.setsdf$symbol[1]
x = strsplit(x, ", ")[[1]]
print("[STATUS] Running GO on cls")
ts.snpfile = paste0(regpref, examp, '_go_snpcls_', suffix)
if (!file.exists(ts.snpfile)){
    fullts.resdf = c()
    ts.resdf = c()
    for (i in 1:nrow(tree.setsdf)){
        x = tree.setsdf$symbol[i]
        cls = tree.setsdf$cls[i]
        x = strsplit(x, ", ")[[1]]
        df = go.enr(x)
        if (!is.null(df)){
            df$cls = cls
            fullts.resdf = rbind(fullts.resdf, df)
            ts.resdf = rbind(ts.resdf, df[1:min(3,nrow(df)),])
        }
    }
    save(fullts.resdf, ts.resdf, file=ts.snpfile)
} else {
    load(ts.snpfile)
}

tree.mid = apply(cbind(c(0,snp.breaks), c(snp.breaks, 1)),1,mean)
descdf = aggregate(Description ~ cls, ts.resdf, function(x){paste(x, collapse=', ')})
descdf$Description = sapply(descdf$Description, width=40, split.text)
txtdesc = descdf$Description
tdloc = tree.mid[descdf$cls]


# ----------------------
# Statistics on snp mat:
# ----------------------
sum(apply(smat,1, sum) > ncol(smat) / 2) # 25 of 339
# hist(apply(smat,1, sum))
smat.group = smat %*% make.tform(nodegroup)
apply(smat.group > 0, 2, sum)


# --------------------------------
# Aggregate GO test sets per node:
# --------------------------------
print("[STATUS] Getting nearest expressed gene")
smatdf = data.frame(t(smat[,nodenames[nord]]))
smatdf$node = rownames(smatdf)
smatdf = gather(smatdf, snp, value, -node)
smatdf = smatdf[smatdf$value > 0,]
smatdf$pval = as.numeric(sub(".* - ", "",smatdf$node))
smatdf$GROUP = sub(" - .*", "",smatdf$node)
smarg = apply(smat, 1, sum)
smargdf = data.frame(snp=names(smarg), marg=smarg)
smatdf = merge(smatdf, smargdf)
smatdf$queryHits = as.numeric(sub("s","",smatdf$snp))
smatdf$chr = paste0('chr',gwdf$chrom[smatdf$queryHits])
smatdf$start = gwdf$chromStart[smatdf$queryHits]
# Add gene expr, tss:
stssdf = merge(tssdf, smatdf)
stssdf$dist = abs(stssdf$tss - stssdf$start)
# hist(red.pcdf$log2fpkm, 50)
red.pcdf = red.pcdf[red.pcdf$log2fpkm > 2,]
stssdf = merge(stssdf, red.pcdf)
# Add nearest (expressed gene)
stssdf = merge(stssdf, aggregate(dist ~ snp, stssdf, min))
stssdf = merge(stssdf, gmdf)
# If any snp missing, add the nearest gene:
smatdf = merge(smatdf, mapsnpdf)
names(smatdf)[ncol(smatdf)] = 'nearest.symbol'
stssdf = merge(smatdf, stssdf, all.x=TRUE)
sind = is.na(stssdf$symbol)
stssdf$symbol[sind] = stssdf$nearest.symbol[sind]
red.stssdf = stssdf[stssdf$marg < 0.5 * ncol(smat),]
# Node centered sets:
tree.nodesdf = aggregate(symbol ~ node + pval, stssdf, function(x){paste(unique(x), collapse=", ")})
tree.nodesdf = tree.nodesdf[order(tree.nodesdf$pval, decreasing=T),]

print("[STATUS] Running GO on nodes")
tn.nodefile = paste0(regpref, examp, '_go_nodes_', suffix)
tn.resdf = c()
fulltn.resdf = c()
if (!file.exists(tn.nodefile)){
    for (i in 1:nrow(tree.nodesdf)){
        x = tree.nodesdf$symbol[i]
        x = strsplit(x, ", ")[[1]]
        df = go.enr(x, minsize=10, maxsize=150)
        if (!is.null(df)){
            node = tree.nodesdf$node[i]
            df$node = node
            df$GO = rownames(df)
            fulltn.resdf = rbind(fulltn.resdf, df)
            tn.resdf = rbind(tn.resdf, df[1:min(3,nrow(df)),])
        }
    }
    tn.resdf$geneID = gsub("/",", ", tn.resdf$geneID)
    tn.resdf$geneID = sapply(tn.resdf$geneID, width=50, split.text)
    tn.resdf$Description = sapply(tn.resdf$Description, width=90, split.text)
    tn.resdf$pval = as.numeric(sub("^.* - ", "", tn.resdf$node))
    tn.resdf = tn.resdf[order(tn.resdf$pval, decreasing=T), ] 
    save(fulltn.resdf, tn.resdf, file=tn.nodefile)
} else {
    load(tn.nodefile)
}

tn.resdf$node[tn.resdf$node == 'None - 38.3'] = 'Multiple - 38.3'
tn.resdf$node[tn.resdf$node == 'None - 26.8'] = 'Multiple - 26.8'
fulltn.resdf$node[fulltn.resdf$node == 'None - 38.3'] = 'Multiple - 38.3'
fulltn.resdf$node[fulltn.resdf$node == 'None - 26.8'] = 'Multiple - 26.8'


# Prune terms by qvalue first:
sub.godf = fulltn.resdf[, c('node','Description', 'qvalue', 'GO')]
if (examp == 'CAD'){
    sub.godf = sub.godf[sub.godf$qvalue < 1e-2,]
} else {
    sub.godf = sub.godf[sub.godf$qvalue < 1e-1,]
}
sub.godf = sub.godf[nchar(sub.godf$Description) < 60,]
ndesc = aggregate(node ~ GO + Description, sub.godf, length)
ndesc = ndesc[order(ndesc$node, decreasing=T), ] 
go.cutoff = 5
go.desc = ndesc$Description[ndesc$node <= go.cutoff]
go.id = ndesc$GO[ndesc$node <= go.cutoff]

# Get similarity matrix of GO terms for pruning:
go.simple = FALSE
if (go.simple){
    library(GOSemSim)
    d <- godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
    NGO = length(go.id)
    comp.go = matrix(NA, nrow=NGO, ncol=NGO, dimnames=list(go.id, go.id))
    for (i in 1:NGO){
        for(j in 1:NGO){
            comp.go[i,j] = goSim(go.id[i], go.id[j], semData=d, measure="Wang")
        }
    }
    image(comp.go, useRaster=T)
}

# Keep 1-3 terms per term, make dataframe with all incidences:
sub.godf = sub.godf[order(sub.godf$qvalue),]
sub.godf.old = sub.godf
sub.godf = sub.godf[sub.godf$Description %in% go.desc,]
# go.keep = go.desc
go.keep = c()
if (examp == 'CAD'){ NKEEP = 2 } else {NKEEP = 4 }
for(ni in nodenames[nord]){ 
    go.keep = c(go.keep, head(sub.godf$Description[sub.godf$node == ni], NKEEP)) }
go.keep = unique(go.keep)

process_godf = function(ssdf){
    ssdf$log10q = -log10(ssdf$qvalue)
    sswide = spread(ssdf[,c('Description','node','log10q')], node, log10q, fill=0)
    mat = as.matrix(sswide[,-1])
    rownames(mat) = sswide$Description
    rn = as.character(rev(nodenames[nord]))
    rn = rn[rn %in% colnames(mat)]
    mat = mat[, rn]
    rd = diag.mat2(t(mat) > 0)
    go.ord = rd[[2]]
    # Make nicer plot:
    ssdf$node = factor(ssdf$node, levels=nodenames[nord])
    ssdf$Description = factor(ssdf$Description, levels=rev(go.ord))
    ssdf$cex = 1.5
    ssdf$cex[ssdf$qvalue < 1e-4] = 1
    ssdf$cex[ssdf$qvalue < 1e-3] = .75
    lindf = rbind(data.frame(x0=1:length(go.ord), x1=1:length(go.ord), y0=1, y1=NPLOT), 
                  data.frame(x0=1, x1=length(go.ord), y0=1:NPLOT, y1=1:NPLOT))
    ylindf = data.frame(x0=1:length(go.ord), x1=1:length(go.ord), y0=1, y1=NPLOT)
    return(list(ssdf, lindf, ylindf, go.ord))
}

# Get order for GO terms:
sub.godf = sub.godf[sub.godf$Description %in% go.keep,]
ll = process_godf(sub.godf)
ssdf = ll[[1]]
lindf = ll[[2]]
ylindf = ll[[3]]
go.ord = ll[[4]]

ll.old = process_godf(sub.godf.old) # For plotting all
ssdf.old = ll.old[[1]]
lindf.old = ll.old[[2]]
ylindf.old = ll.old[[3]]
go.ord.old = ll.old[[4]]


# -----------------------------
# Put the full figure together:
# -----------------------------
NCAT=1000
lvl=1
plot.tree =FALSE
# png(paste0(cadpref, 'test.png'), res=450, units='in', width=7, height=10 * plot.tree + 2 * (NPLOT / 20) + 1)
pdf(paste0(cadpref, 'test.pdf'), width=7, height=10 * plot.tree + 2 * (NPLOT / 20) + 1)
# Set up:
if (plot.tree){
    layout(matrix(c(1:12),2,6), heights=(10, 3), widths=c(1.5,1,2,.75, 4,3))
} else {
    layout(matrix(c(1:6),1,6), widths=c(1.5,1,2,.75,4,3))
}
par(xaxs='i')
par(yaxs='i')
# 1. Top enriched nodes + pval (color + #)
bsp = 10
sp = 0.2
if (plot.tree){
    # Dendrogram:
    plotind = which(uids == suid)
    ipref = paste0(apref, sprintf("_%05d", plotind))
    cregfile = paste0(regpref, ipref, '_lreg_adj', NCAT, '_', lvl, '.Rda')
    # cregfile = paste0(regpref, ipref, '_hg_adj', NCAT, '_', lvl, '.Rda')
    load(cregfile)
    dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x=0, y=0.5, width=circle_size, 
                          height=circle_size, just=c("left", "center")))
    par(omi=gridOMI(), new=TRUE)
    circleplot(dend=dd$dend, lab=lab, udf=dd$udf,
               fractional=FALSE, hit.track=dd$hitdf, leaf.pval=rmat[,suid])
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x=circle_size, just="left")
    title(suid)
}
par(mar=c(bsp, 5, sp,sp))
image(t(as.matrix(nodelp[nord])), axes=F, zlim=c(3,max(nodelp)), useRaster=T)
box(lwd=.5)
text(x=par()$usr[1] - .25 * diff(par()$usr[1:2]), y=seq(0,1,length.out=NPLOT), 
     nodegroup[nord], col=nodecol[nord], srt=0, adj=1, xpd=TRUE, cex=.75)
text(x=mean(par()$usr[1:2]), y=seq(0,1,length.out=NPLOT), 
     round(rawnodelp[nord],1), col='black', srt=0, xpd=TRUE, cex=.5)
# 2. Enrichments of enhancers (start with the enrichments of the modules?)
# - scrapped - see the earlier figure. could do surprisal?
# (enrichmat)
# - set size
par(mar=c(bsp, sp, sp,sp))
barplot(cad.clen[nord], horiz=T, col=NA, border=NA, axes=F)
minor=50000
abline(v=seq(minor, max(cad.clen), minor), col='grey', lwd=.5, lty='dotted')
barplot(cad.clen[nord], horiz=T, col=sapply(nodecol[nord], tsp.col), border=NA, axes=F, add=TRUE)
barplot(cad.dlen[nord], horiz=T, col=nodecol[nord], border=NA, axes=F, add=TRUE)
axis(1, at=seq(0,max(cad.dlen[nord]), 1e5), las=2, cex.axis=.75, lwd=.5)
# 3. Enhancer specificity (enhmargfile)
par(mar=c(bsp, sp, sp,sp))
boxplot(col ~ node, sub.cad.setdf, cex=.15, pch=19, ylab='', xlab='', axes=F,
        horizontal=T, width=rep(.1,NPLOT), col=nodecol[nord], lwd=.5) 
atpar = seq(0,800,200)
abline(v=atpar, col='grey', lwd=.5, lty='dotted')
boxplot(col ~ node, sub.cad.setdf, cex=.15, pch=19, ylab='', xlab='', axes=F,
        horizontal=T, width=rep(.1,NPLOT), col=nodecol[nord], lwd=.5, add=TRUE) 
axis(1, at=atpar,las=2, xpd=TRUE, cex.axis=.75, lwd=.5)
# - Shared enrichments (tissue) with other traits:
image(shmat, zlim=c(3,max(shmat)), axes=F, useRaster=T)
xat = seq(0, 1, length.out=nrow(shmat))
yat = seq(0, 1, length.out=ncol(shmat))
text(x=xat[as.numeric(shdf$uid)], y=yat[as.numeric(shdf$node)], 
     labels=round(shdf$value,0), col=shdf$col, cex=.45)
text(y=par()$usr[3] - .005 * diff(par()$usr[3:4]), x=xat,
     sub("^[0-9]* - ", "", rownames(shmat)), col='black', srt=90, adj=1, xpd=TRUE, cex=.5)
box(lwd=.5)
# - GO enrichments for nearest gene
# - (OR) GO enrichments for nearest expressed gene
par(mar=c(bsp, sp, sp,sp))
plot(1, type='n', axes=F, ylab='', xlab='',
     ylim=c(.5,NPLOT+.5), xlim=c(.5, length(go.ord) + .5))
segments(x0=lindf$x0, x1=lindf$x1, 
         y0=lindf$y0, y1=lindf$y1, 
         lwd=.25, lty='dotted')
# abline(v=1:length(go.ord), lwd=.25) # , lty='dotted')
# abline(h=1:NPLOT, lwd=.25)
par(xpd=TRUE)
points(as.numeric(ssdf$Description), as.numeric(ssdf$node), 
       ylim=c(.5,NPLOT+.5), xlim=c(.5, length(go.ord) + .5),
       ylab='', xlab='', pch=22, cex=1.75, lwd=.5, col=NA, bg='firebrick')
par(xpd=FALSE)
# bg=ssdf$log10q, axes=F)
text(x=1:length(go.ord), y=par()$usr[3] - 0.01 * diff(par()$usr[3:4]),
     # sapply(rev(go.ord), width=40,split.text), srt=90, adj=1, xpd=TRUE, cex=0.5)
     rev(go.ord), srt=90, adj=1, xpd=TRUE, cex=0.5)
box(lwd=.5)
# - Clustered SNP matrix (by all cons)
par(mar=c(bsp, sp, sp,sp))
image(smat[snp.cons.rn, nodenames[nord]], col=c('white','darkblue'), axes=F, useRaster=T)
abline(v=snp.breaks, lwd=.25)
abline(h=seq(0,1,length.out=NPLOT), lwd=.25, lty='dotted')
# plot(1, type='n', axes=F, ylab='', xlab='')
box(lwd=.5)
yt = par()$usr[3] - 0.005*diff(par()$usr[3:4])
text(tdloc, y=yt, labels=txtdesc, srt=90, adj=1, xpd=TRUE,cex=.3)
dev.off()


# --------------------------------
# Put the reduced figure together:
# --------------------------------
NCAT=1000
lvl=1
# png(paste0(cadpref, 'reduced_panel.png'), res=450, units='in', width=4, height=3)
pdf(paste0(cadpref, 'reduced_panel.pdf'), width=4, height=3)
layout(matrix(c(1:5),1,5), widths=c(1.5,1,.75,.75,2.5))
par(xaxs='i')
par(yaxs='i')
# 1. Top enriched nodes + pval (color + #)
bsp = 12
sp = 0.1
par(mar=c(bsp, 6, sp,sp))
image(t(as.matrix(nodelp[nord])), axes=F, zlim=c(3,max(nodelp)), useRaster=T)
box(lwd=.5)
yat = seq(0,1,length.out=NPLOT)
text(x=par()$usr[1] - .25 * diff(par()$usr[1:2]), y=yat, 
     leafrep, col=nodecol[nord], srt=0, adj=1, xpd=TRUE, cex=.5)
text(x=mean(par()$usr[1:2]), y=seq(0,1,length.out=NPLOT), 
     round(rawnodelp[nord]), col=ifelse(rawnodelp[nord] > 40, 'white','black'), srt=0, xpd=TRUE, cex=.5)
# 2. Enrichments of enhancers (start with the enrichments of the modules?)
# - scrapped - see the earlier figure. could do surprisal?
# (enrichmat)
# - set size
par(mar=c(bsp, sp, sp,sp))
barplot(cad.clen[nord], horiz=T, col=NA, border=NA, axes=F)
minor=50000
abline(v=seq(minor, max(cad.clen), minor), col='grey', lwd=.5, lty='dotted')
barplot(cad.clen[nord], horiz=T, col=sapply(nodecol[nord], tsp.col), border=NA, axes=F, add=TRUE)
barplot(cad.dlen[nord], horiz=T, col=nodecol[nord], border=NA, axes=F, add=TRUE)
atpar = seq(0,max(cad.dlen[nord]), 1e5)
axis(1, at=atpar, labels = rep('', length(atpar)), las=2, cex.axis=.75, lwd=.5)
text(x=atpar, y=par()$usr[3] - 0.065 * diff(par()$usr[3:4]), 
     labels=comma_format()(atpar), cex=.65,srt=90, adj=1, xpd=NA)
box(lwd=.5)
# - Shared enrichments (tissue) with other traits:
image(shmat, zlim=c(3,max(shmat)), axes=F, useRaster=T)
xat = seq(0, 1, length.out=nrow(shmat))
yat = seq(0, 1, length.out=ncol(shmat))
text(x=xat[as.numeric(shdf$uid)], y=yat[as.numeric(shdf$node)], 
     labels=round(shdf$value,0), col=shdf$col, cex=.45)
text(y=par()$usr[3]-.01 * diff(par()$usr[3:4]), x=xat,
     sub("^[0-9]* - ", "", rownames(shmat)), col='black', 
     srt=90, adj=1, xpd=TRUE, cex=.65)
box(lwd=.5)
plot(1,1, axes=F, type='n', ylab='', xlab='')
text(1,1, 'Pleiotropy', cex=.5)
box(lwd=.5)
# - GO enrichments for nearest gene
# - (OR) GO enrichments for nearest expressed gene
par(mar=c(bsp, sp, sp, sp))
plot(1, type='n', axes=F, ylab='', xlab='',
     ylim=c(.5,NPLOT+.5), xlim=c(.5, length(go.ord) + .5))
segments(x0=lindf$x0, x1=lindf$x1, 
         y0=lindf$y0, y1=lindf$y1, 
         lwd=.1, lty='solid')
par(xpd=TRUE)
points(as.numeric(ssdf$Description), as.numeric(ssdf$node), 
       ylim=c(.5,NPLOT+.5), xlim=c(.5, length(go.ord) + .5),
       ylab='', xlab='', pch=22, cex=ssdf$cex, lwd=.25,
       col='black', bg=nodecol[nord][as.numeric(ssdf$node)])
par(xpd=FALSE)
text(x=1:length(go.ord), y=par()$usr[3] - 0.01 * diff(par()$usr[3:4]),
     rev(go.ord), srt=90, adj=1, xpd=TRUE, cex=0.5)
box(lwd=.5)
dev.off()



# --------------------------------
# Flip the reduced figure:
# --------------------------------
NCAT=1000
lvl=1
# png(paste0(cadpref ,'reduced_panel_flipped.png'), res=450, units='in', width=3, height=3.5)
pdf(paste0(cadpref, 'reduced_panel_flipped.pdf'), width=3, height=2.75)
layout(matrix(c(1:3),3,1), heights=c(1.7,1,2.5))
par(xaxs='i')
par(yaxs='i')
# 1. Top enriched nodes + pval (color + #)
bsp = 12
sp = 0.1
par(mar=c(sp,bsp, 6,sp))
image(as.matrix(nodelp[rev(nord)]), axes=F, zlim=c(3,max(nodelp)), useRaster=T)
box(lwd=.5)
xat = seq(0,1,length.out=NPLOT)
text(y=par()$usr[4] + .25 * diff(par()$usr[3:4]), x=xat, 
     rev(leafrep), col=rev(nodecol[nord]), srt=90, adj=0, xpd=TRUE, cex=.5)
text(y=mean(par()$usr[3:4]), x=xat, 
     rev(round(rawnodelp[nord])), srt=0, xpd=TRUE, cex=.45,
     col=rev(ifelse(rawnodelp[nord] > 40, 'white','black')))
# 2. Enrichments of enhancers (start with the enrichments of the modules?)
# - scrapped - see the earlier figure. could do surprisal?
# (enrichmat)
# - set size
par(mar=c(sp, bsp, sp,sp))
barplot(rev(cad.clen[nord]), horiz=F, col=NA, border=NA, axes=F)
minor=50000
abline(h=seq(minor, max(cad.clen), minor), col='grey', lwd=.5, lty='dotted')
barplot(rev(cad.clen[nord]), horiz=F, 
        col=rev(sapply(nodecol[nord], tsp.col)), border=NA, axes=F, add=TRUE)
barplot(rev(cad.dlen[nord]), horiz=F, col=rev(nodecol[nord]), border=NA, axes=F, add=TRUE)
atpar = seq(0,max(cad.dlen[nord]), 1e5)
axis(2, at=atpar, labels = rep('', length(atpar)), las=1, cex.axis=.75, lwd=.5)
text(y=atpar, x=par()$usr[1] - 0.065 * diff(par()$usr[1:2]), 
     labels=comma_format()(atpar), cex=.65,srt=0, adj=1, xpd=NA)
box(lwd=.5)
# - Shared enrichments (tissue) with other traits:
# SR = nrow(shmat)
# SC = ncol(shmat)
# shrn = rev(1:SR)
# shcn = rev(1:SC)
# image(t(shmat[shrn,shcn]), zlim=c(3,max(shmat)), axes=F, useRaster=T)
# xat = seq(0, 1, length.out=SC)
# yat = seq(0, 1, length.out=SR)
# text(y=1 - yat[as.numeric(shdf$uid)], x=1 - xat[as.numeric(shdf$node)], 
#      labels=round(shdf$value,0), col=shdf$col, cex=.45)
# text(x=par()$usr[1]-.01 * diff(par()$usr[1:2]), y=yat,
#      rev(sub("^[0-9]* - ", "", rownames(shmat))), col='black', 
#      srt=0, adj=1, xpd=TRUE, cex=.65)
# box(lwd=.5)
# - GO enrichments for nearest gene
# - (OR) GO enrichments for nearest expressed gene
par(mar=c(sp, bsp, sp, sp))
plot(1, type='n', axes=F, ylab='', xlab='',
     xlim=c(.5,NPLOT+.5), ylim=c(.5, length(go.ord) + .5))
segments(y0=lindf$x0, y1=lindf$x1, 
         x0=lindf$y0, x1=lindf$y1, 
         lwd=.1, lty='solid')
par(xpd=TRUE)
points(y=length(go.ord) - as.numeric(ssdf$Description) + 1, 
       x=NPLOT - as.numeric(ssdf$node) + 1, 
       xlim=c(.5,NPLOT+.5), ylim=c(.5, length(go.ord) + .5),
       ylab='', xlab='', pch=22, cex=ssdf$cex, lwd=.25,
       col='black', bg=nodecol[nord][as.numeric(ssdf$node)])
par(xpd=FALSE)
text(y=1:length(go.ord), x=par()$usr[3] - 0.01 * diff(par()$usr[3:4]),
     go.ord, srt=0, adj=1, xpd=TRUE, cex=0.55)
box(lwd=.5)
dev.off()


# ---------------------------------------
# Plot the GO terms as a separate figure:
# ---------------------------------------
plot_go_partial = function(rnrow){
    bsp=12
    tsp=5
    ssdf.old$row = length(go.ord.old) - as.numeric(ssdf.old$Description) + 1 
    sub.ssdf = ssdf.old[ssdf.old$row %in% rnrow,]
    go.ord.sub = go.ord.old[rnrow]
    lindf.old = rbind(data.frame(x0=rnrow, x1=rnrow, y0=1, y1=NPLOT), 
                      data.frame(x0=min(rnrow), x1=max(rnrow), y0=1:NPLOT, y1=1:NPLOT))
    par(mar=c(sp, bsp, tsp, sp))
    ylim=c(min(rnrow) - .5, max(rnrow) + .5)
    print(ylim)
    plot(1, type='n', axes=F, ylab='', xlab='', 
         xlim=c(.5,NPLOT+.5), ylim=ylim)
    segments(y0=lindf.old$x0, y1=lindf.old$x1, 
             x0=lindf.old$y0, x1=lindf.old$y1, 
             lwd=.1, lty='solid')
    print(head(go.ord.sub))
    text(y=rnrow, 
         x=par()$usr[1] - 0.01 * diff(par()$usr[1:2]),
         go.ord.sub, srt=0, adj=1, xpd=TRUE, cex=0.45)
    par(xpd=TRUE)
    points(y=sub.ssdf$row, 
           x=NPLOT - as.numeric(sub.ssdf$node) + 1, 
           xlim=c(.5,NPLOT+.5), ylim=ylim,
           ylab='', xlab='', pch=22, cex=sub.ssdf$cex/1.5, lwd=.25,
           col='black', bg=nodecol[nord][as.numeric(sub.ssdf$node)])
    par(xpd=FALSE)
    # Node names:
    xat = 1:NPLOT
    text(y=par()$usr[4] + .005 * diff(par()$usr[3:4]), x=xat, 
         rev(leafrep), col=rev(nodecol[nord]), 
         srt=90, adj=0, xpd=TRUE, cex=.5)
    box(lwd=.5)
}


# png(paste0(cadpref, 'all_go_term.png'), res=300, units='in', width=8, height=6)
pdf(paste0(cadpref, 'all_go_term.pdf'), width=8, height=6)
layout(matrix(c(1:2),1,2, byrow=TRUE))
par(xaxs='i')
par(yaxs='i')
# 1. Top enriched nodes + pval (color + #)
midSR = round(length(go.ord.old) / 2)
# - GO enrichments for nearest gene
# - (OR) GO enrichments for nearest expressed gene
plot_go_partial((midSR +1):length(go.ord.old))
plot_go_partial(1:midSR)
dev.off()






# ------------------------------
# Plot the CAD gwas figure (6A):
# ------------------------------
# All Znam in all GWAS past certain 
# for (suid in suidlist) {
# Load regression for suid:
suid = "29212778 - Coronary artery disease" # Have: CAD_UKBIOBANK.gz
plotind = which(uids == suid)
strait = unique(gwdf$trait[gwdf$uid == suid])
spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
sgwssdf = gwssdf[gwssdf$uid == suid,]
ssamp = sgwssdf$sampsize[order(sgwssdf$sampsize, decreasing=T)][1]
sinit = sgwssdf$initSample[order(sgwssdf$sampsize, decreasing=T)][1]
s2init = split.text(sinit, width=90)
ipref = paste0(apref, sprintf("_%05d", plotind))
print(paste(suid, '----- output with prefix', ipref))
regfile = paste0(regpref, ipref, '_lreg', suffix)
load(regfile)

# epigeAt 0.1%
snpfile = paste0(regdir, eprefix, apref, '_epigenomes_hg_all_wsnp', '_adj1000_1.Rda')
load(snpfile)
epi.lp = rmat[suid,]

NTOP=20
dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=3, altline=3, altlwd=c(.5, .25), ntop=NTOP)
# dd$dend = try(set(dd$dend, "branches_lwd", .5))
# Label the top 5 nodes with tissue: 
ldf = ll$df
ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
ldf = ldf[1:NTOP,] # TOP N only
ldf = ldf[ldf$rawlog10p > MINP,]
ntdf = merge(ldf, nodetissue)
ntdf = ntdf[order(ntdf$rawlog10p, decreasing=T),]
ntdf$GROUP = paste0(1:NTOP,'. ', ntdf$GROUP, "\n(", rev(leafrep)[1:NTOP], ")")
ntdf = ntdf[,c('node','GROUP','COLOR')]
names(ntdf) = c('node','symbol','color')
ntdf = merge(ntdf, nodedf)
ntdf$cex=.4

# Setup dendrogram:
# png(paste0(cadpref, 'main_example.png'), res=450, units='in', width=5, height=5)
pdf(paste0(cadpref, 'main_example.pdf'), width=5, height=5)
# Plot small dendrogram
par(mar=c(0,0,1,0))
set.seed(2) # For legend locations - find a seed that looks ok.
circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
           plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=FALSE, 
           nodedf=ntdf, scale=1.5, leaf.pval=epi.lp)
circos.clear()
mtext(strait, side=3, line=0, cex=.7)
mtext(spmid, side=3, line=-.75, cex=.5)
dev.off()





# Get the leaf reduced information (bag of words, max 3 terms):
leafrep = sapply(1:NN, function(x, max.terms=3){
                     blacklist = c('GLAND', 'TISSUE', 'CELL')
                     x = sub("[0-9]*_","",declist$dec[[x]])
                     # x = unique(x)
                     nx = length(x)
                     if (nx > 1){
                         term_words <- strsplit(x, "[ _,.]");
                         tab_all <- sort(table(unlist(term_words)));
                         tab_all = tab_all[tab_all > 1]
                         tab_all = tab_all[!(names(tab_all) %in% blacklist)]
                         x = paste0(tolower(head(names(sort(tab_all, decreasing=T)), max.terms)), collapse=', ')
                     } else { x = tolower(x) }
                     x = capitalize(x)
                     return(x)})
lind = which(leafrep == '')
leafrep[lind] = nodetissue$GROUP[lind]


# Potential small example figures:
suidlist = c("29777097 - Maternal history of Alzheimer's disease",
             "30617256 - Alzheimer's disease or family history of Alzheimer's disease",
             "29777097 - Paternal history of Alzheimer's disease",
             # Immune:
             "30595370 - Hypothyroidism",
             "26192919 - Ulcerative colitis",
             "28067908 - Crohn's disease",
             "28067908 - Inflammatory bowel disease",
             # Eye:
             "30535121 - Macular thickness",
             '23326517 - Age-related macular degeneration',
             "30054594 - Glaucoma", 
             # Other
             "26343387 - Myocardial infarction",
             "30679814 - QT interval",
             # RELATED TO CAD:
             "30595370 - Systolic blood pressure", 
             "29507422 - High density lipoprotein cholesterol levels",
             "30061737 - Atrial fibrillation", 
             "24952745 - QT interval",
             "30578418 - Diastolic blood pressure", 
             "30578418 - Pulse pressure", 
             "30578418 - Systolic blood pressure", 
             "26426971 - Waist-to-hip ratio adjusted for BMI", 
             "30595370 - Waist-hip ratio") 


pdf(paste0(treeimgpref,apref, '_examples_for_main_figures.pdf'), width=2.75,height=3, onefile=T)
for (suid in suidlist) {
    # Load regression for suid:
    plotind = which(uids == suid)
    strait = unique(gwdf$trait[gwdf$uid == suid])
    spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
    sgwssdf = gwssdf[gwssdf$uid == suid,]
    ssamp = sgwssdf$sampsize[order(sgwssdf$sampsize, decreasing=T)][1]
    ipref = paste0(apref, sprintf("_%05d", plotind))
    print(paste(suid, '----- output with prefix', ipref))
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    load(regfile)
    if (sum(ll$log10p) > 0){
        # Setup dendrogram:
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)
        dd$dend = try(set(dd$dend, "branches_lwd", .5))
        # Label the top 5 nodes with tissue: 
        NTOP=5
        ldf = ll$df
        ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
        ldf = ldf[1:NTOP,] # TOP N only
        ldf = ldf[ldf$rawlog10p > MINP,]
        ntdf = merge(ldf, nodetissue)
        ntdf = ntdf[order(ntdf$rawlog10p, decreasing=T),]
        ntdf$GROUP = paste0(1:min(NTOP, nrow(ntdf)),'. ', ntdf$GROUP)
        ntdf = ntdf[,c('node','GROUP','COLOR')]
        names(ntdf) = c('node','symbol','color')
        ntdf = merge(ntdf, nodedf)
        ntdf$cex=.5
        # Plot small dendrogram
        par(mar=c(0,0,1,0))
        circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
                   plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
        circos.clear()
        mtext(strait, side=3, line=0, cex=.7)
        mtext(spmid, side=3, line=-.75, cex=.5)
    }
}
dev.off()


fsuidlist = c("24952745 - QT interval",
             "30535121 - Macular thickness",
             "29777097 - Maternal history of Alzheimer's disease",
             "26426971 - Waist-to-hip ratio adjusted for BMI")

pdf(paste0(treeimgpref,apref, '_examples_for_main_figures_final.pdf'), width=2.75,height=3, onefile=T)
for (suid in fsuidlist) {
    # Load regression for suid:
    plotind = which(uids == suid)
    strait = unique(gwdf$trait[gwdf$uid == suid])
    spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
    sgwssdf = gwssdf[gwssdf$uid == suid,]
    ssamp = sgwssdf$sampsize[order(sgwssdf$sampsize, decreasing=T)][1]
    ipref = paste0(apref, sprintf("_%05d", plotind))
    print(paste(suid, '----- output with prefix', ipref))
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    load(regfile)
    if (sum(ll$log10p) > 0){
        # Setup dendrogram:
        dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)
        dd$dend = try(set(dd$dend, "branches_lwd", .5))
        # Label the top 5 nodes with tissue: 
        NTOP=5
        ldf = ll$df
        ldf = ldf[order(ldf$log10p <= 0, ldf$pout),]
        ldf = ldf[1:NTOP,] # TOP N only
        ldf = ldf[ldf$rawlog10p > MINP,]
        ldf$rank = 1:nrow(ldf)
        ntdf = merge(ldf, nodetissue)
        ntdf$GROUP = paste0(ntdf$rank,'. ', ntdf$GROUP, "\n(", leafrep[ntdf$node], ")")
        ntdf = ntdf[,c('node','GROUP','COLOR')]
        names(ntdf) = c('node','symbol','color')
        ntdf = merge(ntdf, nodedf)
        ntdf$cex=.4
        # Plot small dendrogram
        par(mar=c(0,0,1,0))
        circleplot(dend=dd$dend, lab=lab, fractional=FALSE, hit.track=dd$hitdf, 
                   plot.labels=FALSE, add.metadata=TRUE, reduced.metadata=TRUE, nodedf=ntdf, scale=2)
        circos.clear()
        mtext(strait, side=3, line=0, cex=.7)
        mtext(spmid, side=3, line=-.75, cex=.5)
    }
}
dev.off()








