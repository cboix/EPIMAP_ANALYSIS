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

# Arguments for loading data:
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
gtdir = "gwas_tree_analysis/"
regdir = paste0(gtdir, "regressions/")
perdir = paste0(gtdir, "permuted_catalogs/")
epref = paste0(usetree, '_e', tol, '_')
regpref = paste0(regdir, epref)
perpref = paste0(perdir, epref)
cmd = paste('mkdir -p', gtdir, regdir, perdir)
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

# Make RNA-seq consensus at a variety of log2fpkm cutoffs:
rnallfile = paste0(rnadir, 'rna_consensus_object_', usetree, '_070819.Rdata')
print("[STATUS] Loading rnall from file")
load(rnallfile)

# -------------------------------
# Plotting specific (large) GWAS:
# -------------------------------
# Choose one gwas: 
suid = "29212778 - Coronary artery disease"
slist = c("29212778 - Coronary artery disease",
          "26343387 - Myocardial infarction",
          '30367059 - Thyroid stimulating hormone',
          '30275531 - Triglycerides',
          '22581228 - Fasting blood glucose',
          '20081858 - Fasting blood insulin',
          "30617256 - Alzheimer's disease or family history of Alzheimer's disease")

# -------------------------------------
# Load all gwas enrichment information:
# -------------------------------------
snpfile = paste0(regdir, eprefix, apref, '_epigenomes_hg_all_wsnp', suffix)
load(snpfile)
colnames(rmat) = paste0('c', c(1:ncol(rmat) - 1))
# Order the top studies:
filepref = 'cls_merge2_wH3K27ac100_raw'
namesfile = paste0(filepref, '_names.tsv')
epinames = scan(namesfile,'c')
clsn = paste0('c', 1:length(epinames) - 1)
names(epinames) = clsn
rmat = t(rmat)
epimat = rmat
rownames(epimat) = epinames[rownames(rmat)]


pdf(paste0(treeimgpref, 'exploratory_logreg_example_large.pdf'),width=14.5,height=12, onefile=T)
for (suid in slist){
    MINP=3
    plotind = which(uids == suid)
    ipref = paste0(apref, sprintf("_%05d", plotind))
    print(paste(suid, '----- output with prefix', ipref))
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    load(regfile)
    # Setup dendrogram:
    dd  = setup_dendrogram(dend3, ll, udf, declist=declist, bcutoff=MINP, altline=0)
    plot.new()
    circle_size = unit(1, "snpc") # snpc unit gives you a square region
    pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                          just = c("left", "center")))
    par(omi = gridOMI(), new = TRUE)
    circleplot(dend=dd$dend, lab=lab, udf=dd$udf, fractional=FALSE, hit.track=dd$hitdf, leaf.pval=rmat[,suid])
    upViewport()
    circos.clear()
    draw(pd.legend.ext, x = circle_size, just = "left")
    title(paste(suid, '-', ll$title))
}
dev.off()




# ---------------------------------
# Explore diff traits specifically:
# ---------------------------------
slist = c("29212778 - Coronary artery disease",
          "26343387 - Myocardial infarction",
          "29507422 - High density lipoprotein cholesterol levels",
          "29507422 - Low density lipoprotein cholesterol levels",
          "29507422 - Total cholesterol levels",
          "29507422 - Triglycerides",
          #
          "30804560 - Lung function (FEV1/FVC)",
          "30804558 - Autism and major depressive disorder (MTAG)",
          "29403010 - Mean arterial pressure",
          #
          "30595370 - Systolic blood pressure",
          "27841878 - Systolic blood pressure",
          "30595370 - Cardiovascular disease",
          "29892015 - Atrial fibrillation",
          "30061737 - Atrial fibrillation",
          "24952745 - QT interval", # All heart?
          #
          "29531354 - Ischemic stroke",
          "29531354 - Stroke",
          "27588450 - Glomerular filtration rate",
          "25628336 - Motion sickness", 
          #
          # "29187730 - Mood instability",
          "30595370 - Mean corpuscular hemoglobin",
          "30038396 - Cognitive performance (MTAG)",
          "30038396 - Highest math class taken (MTAG)",
          "30038396 - Educational attainment (years of education)",
          "30038396 - Educational attainment (MTAG)",
          "30289880 - Lipid traits (pleiotropy) (HIPO component 1)",
          "27182965 - Allergy",
          "30535121 - Macular thickness",
          #
          "30054594 - Intraocular pressure",
          "30617256 - Alzheimer's disease or family history of Alzheimer's disease")



# Preparation:
tree.enhsets = cdll$diff
tree.enhlens = cdlenlist$diff
tree.par.enhsets = cdll$cons
tree.par.enhlens = cdlenlist$cons
# 
nodetissue = nodetissue[order(nodetissue$node),]
flatset = 'epigenomes'
lvl = 10
NCAT = 1000
enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
load(enhsetfile)
MINP=3

# "25628336 - Motion sickness"
# "28107422 - Fibrinogen levels"
# "29777097 - Family history of Alzheimer's disease"
# suid = "30595370 - Heel bone mineral density"
# suid = "30595370 - Waist−hip ratio"


for (suid in slist){
    plotind = which(uids == suid)
    ipref = paste0(apref, sprintf("_%05d", plotind))

    pdf(paste0(treeimgpref, ipref, 'exploratory_analysis.pdf'),width=14.5,height=12, onefile=T)

    # --------------------
    # Load the regression:
    # --------------------
    print(paste(suid, '----- output with prefix', ipref))
    regfile = paste0(regpref, ipref, '_lreg', suffix)
    nsnpdf = aggregate(pValue ~ uid, gwdf, length)
    nsnps = nsnpdf$pValue[nsnpdf$uid == suid]
    nhits = qallnumdf$queryHits[qallnumdf$uid == suid]
    ksnps = qnumdf$queryHits[qallnumdf$uid == suid]
    sub.qdf = qdf[qdf$uid == suid, ]
    load(regfile)
    # Counts of snps:
    nepi = sum(rmat[,suid] > 3)
    ntree = sum(ll$rawlp > 3)
    # Descriptions:
    snp.str = paste(ksnps, 'of', nsnps, 'snps in', nhits, 'enhancers')
    num.str = paste0(nepi," of 833 epi. (", round(nepi / 833 * 100, 1), "%) vs. ",
                     ntree, " of 1665 nodes (", round(ntree/1665 * 100, 1), "%)")
    print(snp.str)
    print(num.str)

    # ----------------------
    # Figure 0: Dendrogram: 
    # ----------------------
    # Setup dendrogram:
    print("[STATUS] Plotting dendrogram")
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

    # -------------------------------
    # Get the snp sets for each node:
    # -------------------------------
    rind = which(ll$rawlp > 3)
    nodes = ll$df$node[rind]
    lp = ll$rawlp[rind]
    parents = declist$parent[nodes]
    siblings = sapply(nodes, function(x){
                          p = declist$parent[x]
                          d = which(declist$parent == p)
                          d[d != x] })
    lp.par = - log10(ll$df$pout)[parents]
    pn = parents[parents %in% nodes]

    # ----------------------------------------
    # Get node metadata, make the snps matrix:
    # ----------------------------------------
    nodegroup = nodetissue$GROUP[nodes]
    nodecol = nodetissue$COLOR[nodes]
    nodelp = round(ll$rawlp[rind],1)
    nodenames = paste0(nodegroup, ' - ', nodelp)

    # Get enhancers in nodes/sets:
    out = sapply(nodes, function(x){tree.enhsets[[x]]})
    omat = sapply(out, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
    out.all = sapply(1:NN, function(x){tree.enhsets[[x]]})
    omat.all = sapply(out.all, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
    out.cons = sapply(1:NN, function(x){tree.par.enhsets[[x]]})
    omat.cons = sapply(out.cons, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 

    colnames(omat) = nodenames
    names(nodecol) = nodenames
    tform = make.tform(sub.qdf$queryHits)
    smat = t(tform) %*% omat
    smat.all = t(tform) %*% omat.all
    smat.cons = t(tform) %*% omat.cons
    smat[smat > 0] = 1
    smat.all[smat.all > 0] = 1
    smat.cons[smat.cons > 0] = 1
    snps = unique(sub.qdf$queryHits)
    rownames(smat) = paste0('s', snps)
    rownames(smat.all) = paste0('s', snps)
    rownames(smat.cons) = paste0('s', snps)
    colnames(smat.all) = paste0(nodetissue$node, '-',nodetissue$GROUP)
    colnames(smat.cons) = paste0(nodetissue$node, '-',nodetissue$GROUP)
    smat.cons.sig = smat.cons[,rind]
    colnames(smat.cons.sig) = nodenames
    # sind = which(apply(smat,1,sum) > 0)
    # smat = smat[sind,]

    out2 = sapply(parents, function(x){tree.enhsets[[x]]})
    omat2 = sapply(out2, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
    out3 = sapply(parents, function(x){tree.par.enhsets[[x]]})
    omat3 = sapply(out3, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 
    out4 = sapply(siblings, function(x){tree.enhsets[[x]]})
    omat4 = sapply(out4, function(x){ x = enhmap[x]; return(1 * (sub.qdf$subjectHits %in% x)) }) 

    r1 = apply(omat,2,sum) / tree.enhlens[nodes]
    r2 =  apply(omat2, 2, sum) / tree.enhlens[parents]
    r3 =  apply(omat3, 2, sum) / tree.par.enhlens[parents]
    r4 =  apply(omat4, 2, sum) / tree.enhlens[siblings]

    plot(r1, r2, ylim =c(0,max(r1,r2)), xlim=c(0,max(r1,r2)), pch=19)
    text(r1, r2, nodenames, col=nodecol)
    abline(0,1)

    plot(r1, r4, ylim =c(0,max(r1,r4)), xlim=c(0,max(r1,r4)), pch=19)
    text(r1, r4, nodenames, col=nodecol)
    abline(0,1)

    png('~/test.png', res=450, units='in', width=9, height=9)
    plot(r1, r3, ylim =c(0,max(r1,r3)), xlim=c(0,max(r1,r3)), 
         pch=21, cex = lp/10, bg=sapply(nodecol, tsp.col), lwd = lp.par > 10,
         xlab='# hits in node diff / # enhancers',
         ylab='# hits in parent consensus/ # enhancers', 
         main='Signif. nodes for Coronary Artery Disease\n(size=node log10p, text=parent log10p)')
    # text(r1, r3, nodenames, col=nodecol, cex=lp/50)
    text(r1, r3, round(lp.par,1), col=nodecol, cex=lp/50)
    abline(0,1)
    dev.off()

    # -------------------------------------
    strait = unique(gwdf$trait[gwdf$uid == suid])
    spmid = unique(gwdf$pubMedID[gwdf$uid == suid])
    fig.title = paste0(strait,'\n', spmid)
    sgwssdf = gwssdf[gwssdf$uid == suid,]
    sinit = sgwssdf$initSample[order(sgwssdf$sampsize, decreasing=T)][1]
    s2init = split.text(sinit, width=30)
    # Figure 1: similarity matrix of nodes:
    print("[STATUS] Plotting similarity matrix")
    dnj = jacc.dist(t(smat))
    plot.dt.sym(as.dist(dnj), nbreak=ifelse(nrow(dnj) > 10, 3, NULL), txtcol=nodecol)
    text(.5,.75, fig.title, cex=4, col=rgb(0,0,0,.5))
    text(.5,.5, 'Node Similarity\n(SNPs shared)', cex=4, col=rgb(0,0,0,.5))
    text(.5,.25, s2init, cex=2, col=rgb(0,0,0,.5))
    # -------------------------------------


    hnj <- hclust(as.dist(dnj), method='ward.D')
    cocl = order.optimal(as.dist(dnj), hnj$merge)$order
    tree.rn <- names(cocl)[cocl]
    tree.breaks = calc.breaks(hnj, ifelse(nrow(dnj) > 20, 10, 2), cocl)
    # sord = order(apply(smat,1, sum))
    # smat = smat[sord,]
    #--------------------------
    # dsj = dist(smat, 'jaccard')
    # dsj = dist(smat.all, 'jaccard')
    dsj = dist(smat.cons, 'jaccard')
    hsj <- hclust(dsj, method='ward.D')
    snp.cocl = order.optimal(as.dist(dsj), hsj$merge)$order
    snp.rn = names(snp.cocl)[snp.cocl]
    # Cut snps into groups:
    NCENTERS = max(round(nrow(smat) / 30) + 1, 2)
    snp.cto = cutree(hsj, NCENTERS)[snp.cocl]
    tree.sbreaks = calc.breaks.acut(snp.cto)
    tree.smat = smat[snp.rn,tree.rn]
    plot.dt.sym(as.dist(dsj), nbreak=ifelse(nrow(dsj) > 10, 3, NULL))

    # ----------------------------
    # Calculate the nearest genes:
    # ----------------------------
    snpdf = data.frame(queryHits=as.numeric(sub("s","",names(snp.cto))), 
                       snp=names(snp.cto), cls=snp.cto)
    snpdf$chr = paste0('chr',gwdf$chrom[snpdf$queryHits])
    snpdf$start = gwdf$chromStart[snpdf$queryHits]
    stssdf = merge(tssdf, snpdf)
    stssdf$dist = abs(stssdf$tss - stssdf$start)
    stssdf = merge(stssdf, aggregate(dist ~ snp, stssdf, min))
    stssdf = merge(stssdf, gmdf)
    # Mapping for plot:
    mapsnpdf = aggregate(symbol ~ snp, stssdf, function(x)x[1])
    mapsnp = mapsnpdf$symbol
    names(mapsnp) = mapsnpdf$snp


    png('~/test.png', res=450, units='in', width=14, height=9)
    # -----------------------
    # Figure 2: Snp clusters (add nearest gene:
    print("[STATUS] Plotting SNP matrix")
    layout(matrix(c(1:3), 3,1), heights=c(5,5,5), widths=14, TRUE)
    layout(matrix(c(1:1), 1,1), heights=5, widths=14, TRUE)
    sp = 0.1
    par(mar=c(2, 5, sp, 5))
    image(tree.smat, col=c('white','navy'), axes=F)
    box()
    abline(v=tree.sbreaks)
    abline(h=tree.breaks,lty=2,lwd=.5)
    yt = seq(0,1,length.out=ncol(smat))
    xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
    xt2 = par()$usr[2] + 0.005*diff(par()$usr[1:2])
    text(y=yt, x=xt, labels=tree.rn, srt=0, adj=1, 
         xpd=TRUE,cex=.5, col=nodecol[tree.rn])
    text(y=yt, x=xt2, labels=tree.rn, srt=0, adj=0, 
         xpd=TRUE,cex=.5, col=nodecol[tree.rn])
    xt = seq(0,1,length.out=nrow(smat))
    yt = par()$usr[3] - 0.005*diff(par()$usr[3:4])
    xaxlab=mapsnp[rownames(smat[snp.rn,])]
    text(y=yt, x=xt, labels=xaxlab, srt=90, adj=1, 
         xpd=TRUE,cex=.3, col=ifelse(snp.cto %% 2 == 1, 'black', 'darkgrey') )
    # -----------------------
    dev.off()

    dnj.all = dist(t(smat.all), 'jaccard')
    hnj.all <- hclust(dnj.all, method='ward.D')
    cocl = order.optimal(as.dist(dnj.all), hnj.all$merge)$order
    tree.all.rn <- names(cocl)[cocl]
    ncol.all = nodetissue$COLOR[cocl]

    dnj.cons = dist(t(smat.cons), 'jaccard')
    hnj.cons <- hclust(dnj.cons, method='ward.D')
    cocl = order.optimal(as.dist(dnj.cons), hnj.cons$merge)$order
    tree.cons.rn <- names(cocl)[cocl]
    ncol.cons = nodetissue$COLOR[cocl]
    tree.cons.breaks = calc.breaks(hnj.cons, ifelse(nrow(dnj.cons) > 50, 30, 5), cocl)


    png('~/CAD_snp_node_heatmaps.png', res=450, units='in', width=7, height=10)
    # -----------------------
    # Figure 2: Snp clusters:
    print("[STATUS] Plotting SNP matrix")
    layout(matrix(c(1:4), 4,1), heights=c(5,5,5,5), widths=14, TRUE)
    sp = 0.1
    par(mar=c(sp, 5, sp, 5))
    image(tree.smat, col=c('white','navy'), axes=F)
    text(.5,.5, 'Significant nodes\nDifference in consensus at node', cex=2, col=rgb(0,0,0,.5))
    box()
    abline(v=tree.sbreaks)
    abline(h=tree.breaks,lty=2,lwd=.5)
    yt = seq(0,1,length.out=ncol(smat))
    xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
    xt2 = par()$usr[2] + 0.005*diff(par()$usr[1:2])
    text(y=yt, x=xt, labels=tree.rn, srt=0, adj=1, 
         xpd=TRUE,cex=.5, col=nodecol[tree.rn])
    text(y=yt, x=xt2, labels=tree.rn, srt=0, adj=0, 
         xpd=TRUE,cex=.5, col=nodecol[tree.rn])
    # Cons snps:
    par(mar=c(sp, 5, sp, 5))
    image(smat.cons.sig[snp.rn,tree.rn], col=c('white','navy'), axes=F)
    text(.5,.5, 'Significant nodes\nconsensus at node', cex=2, col=rgb(0,0,0,.5))
    abline(v=tree.sbreaks)
    abline(h=tree.breaks,lty=2,lwd=.5)
    text(y=yt, x=xt, labels=tree.rn, srt=0, adj=1, 
         xpd=TRUE,cex=.5, col=nodecol[tree.rn])
    text(y=yt, x=xt2, labels=tree.rn, srt=0, adj=0, 
         xpd=TRUE,cex=.5, col=nodecol[tree.rn])
    box()
    # All diff:
    image(smat.all[snp.rn,tree.all.rn], col=c('white','navy'), axes=F)
    yt = seq(0,1,length.out=ncol(smat.all))
    xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
    text(y=yt, x=xt, labels=tree.all.rn, srt=0, adj=1, 
         xpd=TRUE,cex=.25, col=ncol.all)
    text(.5,.5, 'All nodes\nDifference in consensus at node', cex=2, col=rgb(0,0,0,.5))
    box()
    abline(v=tree.sbreaks)
    # All cons:
    par(mar=c(sp, 5, sp, 5))
    image(smat.cons[snp.rn,tree.cons.rn], col=c('white','navy'), axes=F)
    yt = seq(0,1,length.out=ncol(smat.cons))
    xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
    text(y=yt, x=xt, labels=tree.cons.rn, srt=0, adj=1, 
         xpd=TRUE,cex=.25, col=ncol.cons)
    text(.5,.5, 'All nodes\nConsensus at the node', cex=2, col=rgb(0,0,0,.5))
    box()
    abline(v=tree.sbreaks)
    dev.off()


    # -----------------------
    # Aggregate GO test sets per cluster:
    tree.setsdf = aggregate(symbol ~ cls, stssdf,
                            function(x){paste(unique(x), collapse=", ")})
    # Perform GO enrichment:
    x = tree.setsdf$symbol[1]
    x = strsplit(x, ", ")[[1]]
    print("[STATUS] Running GO on cls")
    ts.snpfile = paste0(regpref, ipref, '_go_snpcls_', suffix)
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
    # ---------------------------------
    # Figure 3: GO enrichments on SNPs:
    png('~/CAD_cluster_GOtable.png', res=450, units='in', width=15, height=10)
    plot.new()
    grid.table(ts.resdf[,c('cls','Description', 'geneID')])
    dev.off()
    # ---------------------------------

    # -----------------------
    # Figure 4: Snp clusters + meta (ts.res)
    png('~/CAD_snp_node_GO_mainheatmaps.png', res=450, units='in', width=10, height=8)
    print("[STATUS] Plotting SNP matrix (with GO)")
    tree.mid = apply(cbind(c(0,tree.sbreaks), c(tree.sbreaks, 1)),1,mean)
    descdf = aggregate(Description ~ cls, ts.resdf, function(x){paste(x, collapse=', ')})
    descdf$Description = sapply(descdf$Description, width=40, split.text)
    txtdesc = descdf$Description
    tdloc = tree.mid[descdf$cls]
    # Figure:
    layout(matrix(c(1:2), 2,1), heights=c(5,5), widths=12, TRUE)
    sp = 0.1
    par(mar=c(sp, 5, 5, 5))
    image(tree.smat, col=c('white','navy'), axes=F)
    box()
    abline(v=tree.sbreaks, lwd=.5)
    abline(h=tree.breaks,lty=2,lwd=.5)
    text(.5,.5, 'Significant nodes\nDifference in consensus at node', cex=2, col=rgb(0,0,0,.5))
    yt = seq(0,1,length.out=ncol(smat))
    xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
    xt2 = par()$usr[2] + 0.005*diff(par()$usr[1:2])
    text(y=yt, x=xt, labels=tree.rn, srt=0, adj=1, 
         xpd=TRUE,cex=.5, col=nodecol[tree.rn])
    text(y=yt, x=xt2, labels=tree.rn, srt=0, adj=0, 
         xpd=TRUE,cex=.5, col=nodecol[tree.rn])
    yt2 = par()$usr[4] + 0.005*diff(par()$usr[3:4])
    text(tdloc, y=yt2, labels=txtdesc, srt=90, adj=0, xpd=TRUE,cex=.3)
    # All cons:
    par(mar=c(sp, 5, sp, 5))
    image(smat.cons[snp.rn,tree.cons.rn], col=c('white','navy'), axes=F)
    yt = seq(0,1,length.out=ncol(smat.cons))
    xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
    text(y=yt, x=xt, labels=tree.cons.rn, srt=0, adj=1, 
         xpd=TRUE,cex=.15, col=ncol.cons)
    text(.5,.5, 'All nodes\nConsensus at the node', cex=2, col=rgb(0,0,0,.5))
    box()
    abline(v=tree.sbreaks, lwd=.5)
    par(xpd=NA)
    abline(h=tree.cons.breaks, lwd=.5, lty='dotted')
    par(xpd=TRUE)
    dev.off()


    # -----------------------

    # --------------------------------
    # Aggregate GO test sets per node:
    # --------------------------------
    print("[STATUS] Getting nearest expressed gene")
    smatdf = data.frame(t(smat))
    smatdf$node = rownames(t(smat))
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
    red.pcdf = red.pcdf[red.pcdf$log2fpkm > 1,]
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
    tree.red.nodesdf = aggregate(symbol ~ node + pval, red.stssdf, function(x){paste(unique(x), collapse=", ")})
    tree.red.nodesdf = tree.red.nodesdf[order(tree.red.nodesdf$pval, decreasing=T),]
    rownames(tree.nodesdf) = NULL
    if (ncol(smat) < 25){ tree.red.nodesdf = tree.nodesdf } # Otherwise everything excluded.

    print("[STATUS] Running GO on nodes")
    tn.nodefile = paste0(regpref, ipref, '_go_nodes_', suffix)
    tn.resdf = c()
    fulltn.resdf = c()
    if (!file.exists(tn.nodefile)){
        for (i in 1:min(10, nrow(tree.red.nodesdf))){
            x = tree.red.nodesdf$symbol[i]
            x = strsplit(x, ", ")[[1]]
            df = go.enr(x)
            if (!is.null(df)){
                node = tree.red.nodesdf$node[i]
                df$node = node
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
    # ---------------------------------
    # Figure 3 (GO enrichments - table)
    plot.new()
    grid.table(tn.resdf[,c('node','Description', 'geneID')])
    # ---------------------------------

    dev.off()
}



# -------------------------------
# Show this but with the epidata:
# -------------------------------

ipref = paste0(apref, '_', flatset, sprintf("_%05d", plotind))
cregfile = paste0(regpref, ipref, '_hg_adj', NCAT, '_', lvl, '.Rda')
rind = which(rmat[,suid] > 3)
nodes = rind

# Get enhancers in nodes/sets:
out = sapply(nodes, function(x){enhsets[[x]]})
sub.qdf = qdf[qdf$uid == suid, ]

# Get enhhits in nodes/sets:
out = sapply(nodes, function(x){enhsets[[x]]})
omat = sapply(out, function(x){
                  x = enhmap[x]
                  return(1 * (sub.qdf$subjectHits %in% x)) }) 


# Add names as node+pval
nodetissue = nodetissue[order(nodetissue$node),]
nodegroup = nodetissue$GROUP[nodes]
nodecol = nodetissue$COLOR[nodes]
nodelp = round(rmat[rind,suid],1)
nodenames = paste0(nodegroup, ' - ', nodelp)
colnames(omat) = nodenames

dnj = jacc.dist(t(omat))
plot.dt.sym(as.dist(dnj))

tform = make.tform(sub.qdf$queryHits)
smat = t(tform) %*% omat
smat[smat > 0] = 1
smat = smat[apply(smat,1,sum) > 0,]

dnj = jacc.dist(t(smat))
plot.dt.sym(as.dist(dnj))

hnj <- hclust(as.dist(dnj), method='ward.D')
cocl = order.optimal(as.dist(dnj), hnj$merge)$order
rn <- names(cocl)[cocl]
breaks = calc.breaks(hnj, 10, cocl)

sord = order(apply(smat,1, sum))
smat = smat[sord,]

# -------------------------
# Cluster snps into groups:
# -------------------------
set.seed(1)
fam <- 'ejaccard'
NCENTERS = 20
cl1 <- kcca(smat[,rn], k=NCENTERS, family=kccaFamily(fam))
summary(cl1)

# Get kcca cluster attributes:
acl <- attributes(cl1)
dcl <- dist(acl$centers,'ejaccard')
hcl <- hclust(dcl)
# Order by max:
cocl = order(apply(acl$centers, 1, which.max), decreasing=TRUE)
cls <- ordered(acl$cluster,levels=cocl)
ord <- order(cls) # will order according to factor
cto = as.numeric(cls[ord])
sbreaks = calc.breaks.acut(cto)
names(nodecol) = nodenames


# Diagonalize the clusters of SNPs:
png('~/CAD_snps_rough_cluster_epigenomes.png', res=450, units='in', width=15, height=4)
layout(matrix(c(1:2), 1,2), widths=c(11,4), TRUE)
sp= 0.1
par(mar=c(sp, 5, sp, 5))
image(smat[ord,rn], col=c('white','navy'), axes=F)
box()
abline(v=sbreaks)
abline(h=breaks,lty=2,lwd=.5)
yt = seq(0,1,length.out=ncol(smat))
xt = par()$usr[1] - 0.005*diff(par()$usr[1:2])
xt2 = par()$usr[2] + 0.005*diff(par()$usr[1:2])
yaxlab=rn
txtcex=0.25
text(y=yt, x=xt, labels=yaxlab, srt=0, adj=1, 
     xpd=TRUE,cex=txtcex, col=nodecol[rn])
text(y=yt, x=xt2, labels=yaxlab, srt=0, adj=0, 
     xpd=TRUE,cex=txtcex, col=nodecol[rn])
# Add heatmap:
par(mar=c(sp, sp, sp, sp))
image(dnj[rn,rn], col=colspec, axes=F)
abline(v=breaks,lty=2,lwd=.5)
abline(h=breaks,lty=2,lwd=.5)
box()
dev.off()



