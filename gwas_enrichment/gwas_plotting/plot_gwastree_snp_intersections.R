#!/usr/bin/R
# -----------------------------------------------------------
# Script to plot/explore large + small (precomputed) examples
# Q: Plot only at 0.1% - or highlight 0.1%, show 1%(?)
# NOTE: Updated 07/09/20 to plot FDR corrections by a single p-value cutoff
# -----------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
library(qvalue)
library(stringr)
library(scales)

# Arguments for loading data:
usetree = 'enhancers'
tol = 2500
singlematch = FALSE
plotting.only = FALSE  # No, need to run the regressions
use.adj = TRUE
use.strict = FALSE
use.onecutoff = TRUE
MINP=3


parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}


# Load in + process all of the relevant matrices/datasets:
commandArgs <- function(trailingOnly=TRUE){
    c(usetree, tol, singlematch, plotting.only, use.adj, use.strict, use.onecutoff) }
source(paste0(bindir, 'load_statistics_gwastree_enrichments.R'))

rm(regwide)
gc()

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

# 800 * 1665

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
nodetissue = nodetissue[order(nodetissue$node),]
leafrep[lind] = nodetissue$GROUP[lind]


# TODO: for each GWAS, pull down the SNPs, get the enh intersections for all the SNPs
# TODO: Calc distance from enhancer in overlap:
# TODO: use the mapping to the node (any underneath??)

qdf = qdf[qdf$uid %in% guid,]
qdf$snploc = gwgr@ranges@start[qdf$queryHits]
dhsmid = (dmgr@ranges@start + dmgr@ranges@width / 2)
qdf$dhsloc = dhsmid[qdf$subjectHits]
qdf$dist = abs(qdf$snploc - qdf$dhsloc)

# ------------------------------------------------
# Get prioritized genes in tissue-specific manner:
# ------------------------------------------------
gmapfile = 'Annotation/gencode.gene.name.mapping.tsv'
tssfile = 'Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz'
genemap = read.delim(gmapfile, header=F)
names(genemap) = c('gene','symbol')
tssdf = read.delim(gzfile(tssfile), header=F, sep="\t", stringsAsFactors=F)
names(tssdf) = c('chr','tss','tss2','gene')
tssdf = tssdf[!(tssdf$chr %in%  c('chrY', 'chrM')),c('chr','tss','gene')]
tssdf = merge(tssdf, genemap)

# Get nearest gene per enh:
tssgr = GRanges(seqnames=tssdf$chr, IRanges(start=tssdf$tss, end=tssdf$tss), name=tssdf$gene)
eout = nearest(dmgr, tssgr)
enhdf$nearest = tssdf$symbol[eout]

# Get nearest gene per snp:
eout = nearest(gwgr, tssgr)
gwdf$nearest = tssdf$symbol[eout]
gwdf$locstr = paste0(gwdf$chrom, '_', gwdf$chromStart)

# For main figures:
GWASids = c('CAD' = '29212778 - Coronary artery disease',
            'BRCA' = '29059683 - Breast cancer',
            'SCZ' = '26198764 - Schizophrenia')

load(paste0(extpref, 'kept_allgwas_ordered', cutsuff, suffix))
# suid = guid[103]
# Want to map to closest active, tissue-specific enhancer
# for (j in 1:length(Znam)) {
#     suid = rev(Znam)[j]

for (j in 1:3){
    suid = GWASids[j]

    print(paste(j, suid))
    trait = sub(".* - ", "", suid)
    pmid = sub(" - .*", "", suid)
    traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
    traitstrnoparen = gsub("\\)" ,"", gsub("\\(","",traitstr))
    if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
        fullid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                       "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
    } else { fullid = suid }
    lp = all.regmat[fullid,]
    lind = which(lp > 0)
    lind = lind[order(lp[lind])] # Sort
    lind = tail(lind, 25)
    sqdf = qdf[qdf$uid %in% fullid,c('queryHits','subjectHits','dist')]


    # Assign closest active enhancer for the node:
    snpdf = ldply(lind, function(i){ 
                      x = cdll$cons[[i]] 
                      x = enhmap[x]; 
                      df = sqdf[sqdf$subjectHits %in% x,]
                      df2 = aggregate(dist ~ queryHits, df, min)
                      df = merge(df, df2)
                      df$node = i
                      return(df) }) 


    snpdf$loc = paste0(gwdf$chrom[snpdf$queryHits],"_", 
                       gwdf$chromStart[snpdf$queryHits])
    snpdf$snp.nearest = gwdf$nearest[snpdf$queryHits]
    snpdf$p = gwdf$pValue[snpdf$queryHits]
    swide = spread(unique(snpdf[,c('node','loc','dist')]),node, dist, fill=5000)
    smat = as.matrix(swide[,-1, drop=F])
    # Set the NA color to white, but not if none are NA:
    if (sum(smat == 5000) > 0){
        matcols = c(rev(col1), 'white')
    } else { matcols = rev(col1[20:100]) }
    if (ncol(smat) > 1){ smat = smat[,as.character(lind), drop=F] }
    rownames(smat) = swide[,1]
    enr.mat = as.matrix(round(lp[lind],1))
    # pdf = unique(snpdf[,c('loc','p')])
    pdf = unique(snpdf[,c('loc','p', 'snp.nearest')])
    pdf$chrom = sub("_.*", "", pdf$loc)
    pdf$pos = comma_format()(as.numeric(sub(".*_", "", pdf$loc)))
    pdf = pdf[order(pdf$p),]
    pdf$fancyloc = paste0('chr',pdf$chrom, ':', pdf$pos, " (", pdf$snp.nearest, ")")
    rownames(pdf) = pdf$loc
    snp.mat = as.matrix(pdf$p)
    rownames(snp.mat) = pdf$loc
    smat = smat[rev(pdf$loc),, drop=F]

    nodetissue = nodetissue[order(nodetissue$node),]
    strait = sub(".* - ", "", suid)
    spmid = sub(" - .*", "", suid)

    rows = (seq(0, ncol(smat), 5) - .5) / (ncol(smat) - 1)
    cols = (seq(0, nrow(smat), 10) - .5) / (nrow(smat) - 1)

    # TODO: Add legends
    mxchar = max(nchar(paste0(20, '. ', leafrep)))
    leaflabs = paste0(1:length(lind), '. ', rev(leafrep[lind]))
    topheight = max(nchar(leaflabs)) / mxchar * 20
    # heights = c(nrow(smat), 2, topheight)
    heights = c(3, topheight, 2, nrow(smat))
    widths= c(6,3.5, ncol(smat) + 2)
    w = sum(widths) / 12
    h = sum(heights) / 12
    strait2 = paste0(split.text(strait, width=35), '\n(PubMedID ', spmid, ')')

    # pdf(paste0(imgpref, 'snp_intersections_', traitstrnoparen, "_", pmid, '.pdf'), width=w, h=h)
    png(paste0(imgpref, 'snp_intersections_', sprintf('page_%d.png', j)), units='in', res=450, width=w, h=h)
    layout(matrix(c(1,1,1, 2:10), nrow=4, ncol=3, byrow=TRUE), widths=widths, heights=heights, TRUE)
    sp = 0.1
    rsp = .25
    par(mar=c(sp,sp,sp,sp))
    # Add title:
    plot(1,1,type='n', axes=F, ylab='', xlab='', xlim=c(0,1), ylim=c(0,1))
    txtcex = 1
    while (strwidth(strait2, cex=txtcex) > 1){ txtcex = 0.95 * txtcex }
    while (strheight(strait2, cex=txtcex) > 1){ txtcex = 0.95 * txtcex }
    text(.5, .5, strait2, cex=txtcex)
    plot(1,type='n', axes=F, ylab='', xlab='')
    plot(1,type='n', axes=F, ylab='', xlab='')
    par(mar=c(sp,sp,sp,rsp))
    image(t(smat), col='white', useRaster=T, axes=F)
    xat = seq(0,1,length.out=ncol(smat))
    yat = seq(0,1,length.out=nrow(smat))
    text(x=xat, y=parpos(2, 0), labels=leaflabs,   
         col = nodetissue$COLOR[rev(lind)], xpd=TRUE, srt=90, cex=.7, adj=0)
    par(mar=c(sp,sp,sp,sp))
    plot(1,type='n', axes=F, ylab='', xlab='')
    plot(1,type='n', axes=F, ylab='', xlab='')
    par(mar=c(sp,sp,sp,rsp))
    image(enr.mat[rev(1:length(lind)),,drop=F], col=col2, useRaster=T, axes=F)
    abline(v=rows, lwd=.25, lty='dashed')
    enr.vals = round(lp[rev(lind)],1)
    text(x=xat, y=parpos(2,-.5), labels=enr.vals,  xpd=TRUE, srt=90,cex=.5, adj=.5,
         col=ifelse(enr.vals > (max(enr.vals) - diff(range(enr.vals))[1] * .5), 'white','black'))
    box(lwd=.5)
    # Labels:
    par(mar=c(sp,sp,sp,sp))
    image(t(smat), col='white', useRaster=T, axes=F)
    text(x=parpos(1, -.95), y=yat, labels=pdf[rownames(smat),'fancyloc'], xpd=TRUE, cex=.5, adj=1)
    # SNPs:
    par(mar=c(sp,sp,sp,sp))
    pltsnp.mat = -log10(snp.mat[rownames(smat),,drop=F])
    image(t(pltsnp.mat), col=col3, useRaster=T, axes=F)
    abline(h=cols, lwd=.25, lty='dashed')
    text(x=parpos(1,-.5), y=yat, labels=sprintf('%0.1e', snp.mat[rownames(smat),]),
         xpd=TRUE, cex=.5, adj=.5, 
         col=ifelse(pltsnp.mat > (max(pltsnp.mat) - diff(range(pltsnp.mat)) * .5), 'white','black'))
    box(lwd=.5)
    # Image
    par(mar=c(sp,sp,sp,rsp))
    if (ncol(smat) == 1){
        plt.mat = smat
    } else {
        plt.mat = smat[,rev(as.character(lind)), drop=F]
    }
    image(t(plt.mat), col=matcols, useRaster=T, axes=F)
    abline(h=cols, lwd=.25, lty='dashed')
    abline(v=rows, lwd=.25, lty='dashed')
    box(lwd=.5)
    # Text on matrix:
    ind = which(plt.mat < 5000, arr.ind=T)
    vals = plt.mat[ind]
    vstr = as.character(round(vals,-1))
    vstr[vals > 1000] = paste0(round(vals[vals >1000]/1000,1),'k')
    text(x=xat[ind[,2]], y=yat[ind[,1]],
         labels=vstr, xpd=TRUE, cex=.3, srt=0, adj=.5, col='white')
    dev.off()

    # -------------------------------------------------------------------
    # Restrict to top N, change legends to have the SNP id + nearest gene
    # -------------------------------------------------------------------
    TOPN = 20
    TOPC = 5
    if (suid == '29212778 - Coronary artery disease'){
        TOPN = 30
        TOPC = 10
    } 
    keep.cols = (ncol(smat)-TOPC + 1):ncol(smat)
    red.smat = tail(smat[,keep.cols], TOPN)
    rows = (seq(0, ncol(red.smat), 5) - .5) / (ncol(red.smat) - 1)
    cols = (seq(0, nrow(red.smat), 10) - .5) / (nrow(red.smat) - 1)

    # TODO: Add legends
    mxchar = max(nchar(paste0(20, '. ', leafrep)))
    red.lind = as.numeric(colnames(red.smat))
    leaflabs = paste0(1:ncol(red.smat), '. ', rev(leafrep[red.lind]))
    topheight = max(nchar(leaflabs)) / mxchar * 20
    # heights = c(nrow(smat), 2, topheight)
    heights = c(3, topheight, 2, nrow(red.smat))
    widths= c(6,3.5, ncol(red.smat) + 2)
    w = sum(widths) / 12
    h = sum(heights) / 12
    strait2 = paste0(split.text(strait, width=35), '\n(PubMedID ', spmid, ')')

    pdf(paste0(imgpref, 'snp_intersections_topn_', TOPN, '_', traitstrnoparen, "_", pmid, '.pdf'), width=w, h=h)
    # png(paste0(imgpref, 'snp_intersections_', sprintf('page_%d.png', j)), units='in', res=450, width=w, h=h)
    layout(matrix(c(1,1,1, 2:10), nrow=4, ncol=3, byrow=TRUE), widths=widths, heights=heights, TRUE)
    sp = 0.1
    rsp = .25
    par(mar=c(sp,sp,sp,sp))
    # Add title:
    plot(1,1,type='n', axes=F, ylab='', xlab='', xlim=c(0,1), ylim=c(0,1))
    txtcex = 1
    while (strwidth(strait2, cex=txtcex) > 1){ txtcex = 0.95 * txtcex }
    while (strheight(strait2, cex=txtcex) > 1){ txtcex = 0.95 * txtcex }
    text(.5, .5, strait2, cex=txtcex)
    plot(1,type='n', axes=F, ylab='', xlab='')
    plot(1,type='n', axes=F, ylab='', xlab='')
    par(mar=c(sp,sp,sp,rsp))
    image(t(red.smat), col='white', useRaster=T, axes=F)
    xat = seq(0,1,length.out=ncol(red.smat))
    yat = seq(0,1,length.out=nrow(red.smat))
    text(x=xat, y=parpos(2, 0), labels=leaflabs,   
         col = nodetissue$COLOR[rev(red.lind)], xpd=TRUE, srt=90, cex=.7, adj=0)
    par(mar=c(sp,sp,sp,sp))
    plot(1,type='n', axes=F, ylab='', xlab='')
    plot(1,type='n', axes=F, ylab='', xlab='')
    par(mar=c(sp,sp,sp,rsp))
    image(enr.mat[rev(1:length(red.lind)),,drop=F], col=col2, useRaster=T, axes=F)
    abline(v=rows, lwd=.25, lty='dashed')
    enr.vals = round(lp[rev(red.lind)],1)
    text(x=xat, y=parpos(2,-.5), labels=enr.vals,  xpd=TRUE, srt=90,cex=.5, adj=.5,
         col=ifelse(enr.vals > (max(enr.vals) - diff(range(enr.vals))[1] * .5), 'white','black'))
    box(lwd=.5)
    # Labels:
    par(mar=c(sp,sp,sp,sp))
    image(t(red.smat), col='white', useRaster=T, axes=F)
    text(x=parpos(1, -.95), y=yat, labels=pdf[rownames(red.smat),'fancyloc'], xpd=TRUE, cex=.5, adj=1)
    # SNPs:
    par(mar=c(sp,sp,sp,sp))
    pltsnp.mat = -log10(snp.mat[rownames(red.smat),,drop=F])
    image(t(pltsnp.mat), col=col3, useRaster=T, axes=F)
    abline(h=cols, lwd=.25, lty='dashed')
    text(x=parpos(1,-.5), y=yat, labels=sprintf('%0.1e', snp.mat[rownames(red.smat),]),
         xpd=TRUE, cex=.5, adj=.5, 
         col=ifelse(pltsnp.mat > (max(pltsnp.mat) - diff(range(pltsnp.mat)) * .5), 'white','black'))
    box(lwd=.5)
    # Image
    par(mar=c(sp,sp,sp,rsp))
    if (ncol(red.smat) == 1){
        plt.mat = red.smat
    } else {
        plt.mat = red.smat[,rev(as.character(red.lind)), drop=F]
    }
    image(t(plt.mat), col=matcols, useRaster=T, axes=F)
    abline(h=cols, lwd=.25, lty='dashed')
    abline(v=rows, lwd=.25, lty='dashed')
    box(lwd=.5)
    # Text on matrix:
    ind = which(plt.mat < 5000, arr.ind=T)
    vals = plt.mat[ind]
    vstr = as.character(round(vals,-1))
    vstr[vals > 1000] = paste0(round(vals[vals >1000]/1000,1),'k')
    text(x=xat[ind[,2]], y=yat[ind[,1]],
         labels=vstr, xpd=TRUE, cex=.3, srt=0, adj=.5, col='white')
    dev.off()

}

# Legend for enhancers:
col_fun = function(x, pal=col1){
    palette = rev(pal)
    bin <- cut(x, seq(0, 5000, length.out=length(palette)), include.lowest=T) 
    palette[bin] }

plegend = Legend(at = seq(0, 5000, 500), 
                    labels_gp = gpar(fontsize=5),
                    title_gp = gpar(fontsize=5, fontface='bold'),
                    col_fun=col_fun, title_position = "topleft", 
                    title='Distance to nearest enhancer (bp)', direction = 'vertical')

pdf(paste0(imgpref,'dist_plegend.pdf'), width=7.75, height=4.25)
draw(plegend, x = unit(4.25,'in'), y=unit(4,'in'), just = "top")
dev.off()




# ----------------------------------------------------------
# Pre-calculate the necessary SNPs to intersect for linking:
# ----------------------------------------------------------
# TODO: SAVE as idset:
prelinkrda = 'linking_data/precalc_gwas_SNP_links.Rda'
if (!file.exists(prelinkrda)){ 
    reqdf = c()
    for (j in 1:length(Znam)) {
        # Want to map to closest active, tissue-specific enhancer
        suid = rev(Znam)[j]
        print(paste(j, suid))
        if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
            fullid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                            "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
        } else { fullid = suid }
        lp = all.regmat[fullid,]
        lind = which(lp > 0)
        lind = lind[order(lp[lind])] # Sort
        lind = tail(lind, 25)
        sqdf = qdf[qdf$uid %in% fullid,c('queryHits','subjectHits','dist')]
        # Assign closest active enhancer for the node:
        snpdf = ldply(lind, function(i){ 
                          x = cdll$cons[[i]] 
                          x = enhmap[x]; 
                          df = sqdf[sqdf$subjectHits %in% x,]
                          df2 = aggregate(dist ~ queryHits, df, min)
                          df = merge(df, df2)
                          df$node = i
                          return(df) }) 
        snpdf$name = enhdf$name[snpdf$subjectHits]
        # Add name x node intersections:
        reqdf = rbind(reqdf, unique(snpdf[,c('node','name')]))
    }

    # Map the indices with the links:
    allids = c()
    nodes = sort(unique(reqdf$node))
    idmap = c()
    for (k in nodes){
        x = declist$dec[[k]] # For breast cancer
        ids = leafmeta[leafmeta$label %in% x,'id']
        idmap = rbind(idmap, data.frame(node=k, id=ids))
    }

    # Load the relevant links:
    linkdf = merge(idmap, reqdf)
    allids = sort(unique(linkdf$id))
    all.links = c()
    t1 = proc.time()
    for (i in 1:length(allids)){
        cat(i,"\t")
        id = allids[i]
        cat(id,"\t")
        tab = read.delim(paste0('linking_data/predictions/collated/', id, '_collated_pred.tsv.gz'), header=F)
        names(tab) = c('chr','start','end','gene','score','enh')
        cat(nrow(tab),"\t")
        # Specific loc:
        enam = linkdf$name[linkdf$id == id]
        sub.enhdf = enhdf[enhdf$name %in% enam,]
        tab = merge(tab, sub.enhdf)
        if (nrow(tab) > 0){
            tab$id = id
            all.links = rbind(all.links, tab[,c('gene','score','nearest', 'enh','name', 'id')])
        }
        t2 = proc.time() - t1
        cat(round(t2[3],1),"s\n")
    }
    save(linkdf, all.links, file=prelinkrda)
} else {
    load(prelinkrda)
}


# -----------------------------------------
# Pre-compute or load group-specific links:
# To be used in plotting around GWAS loci.
# -----------------------------------------
all.group.rda = 'linking_data/all_links_by_group.Rda'
if (!file.exists(all.group.rda)){
    # Enhancer mapping:
    enhdf$loc = with(enhdf, paste(chr, start, end, sep="_"))
    nmapdf = enhdf[,c('loc','name')]
    namemap = enhdf$name
    names(namemap) = enhdf$loc
    # Reduce metadata:
    smeta = meta[meta$id %in% cellorder,]
    groups = as.character(odf$GROUP)
    for (group in groups){
        groupstr = tolower(gsub("__", "", str_replace_all(group, "[^[:alnum:]]", "_")))
        sub.group.rda = paste0('linking_data/links_by_group.', groupstr, '.Rda')
        if (!file.exists(sub.group.rda)){
            relevant.ids = smeta$id[smeta$GROUP == group]
            cat(group, "\t", length(relevant.ids), "\n")
            # TODO: Merge in as we load data:
            t1 = proc.time()
            gdf = c()
            for (i in 1:length(relevant.ids)){
                id = relevant.ids[i]
                cat(paste0(i, '. ', id,"\t"))
                tab = read.delim(paste0('linking_data/predictions/collated/', id, '_collated_pred.tsv.gz'), header=F)
                names(tab) = c('chr','start','end','gene','score','enh')
                tab$loc = with(tab, paste(chr, start, end, sep="_"))
                cat(nrow(tab),"\t")
                tab$name = namemap[tab$loc]
                if (nrow(tab) > 0){
                    tab$id = id
                    gdf = rbind(gdf, tab[,c('gene','score','name','id')])
                }
                t2 = proc.time() - t1
                cat(round(t2[3],1),"s\n")
            }
            gdf = aggregate(score ~ gene + name + id, gdf, mean) 
            gdf = aggregate(score ~ gene + name, gdf, sum) 
            gdf$score = gdf$score / length(relevant.ids)
            gdf$group = group
            save(gdf, file=sub.group.rda)
        } else { 
            load(sub.group.rda) 
        }
    }

    # TODO: Merge all:
    save(group.linksdf, file=all.group.rda)
} else {
    load(all.group.rda)
}



# Test the linking figures:
# suid = guid[103]
all.gwlinked.rda = 'linking_data/all_gwas_SNP_links.Rda'
if (!file.exists(all.gwlinked.rda)){
    alldf = c()
    rownames(enhdf) = enhdf$name
    for (j in 1:length(Znam)) {
        # Want to map to closest active, tissue-specific enhancer
        suid = rev(Znam)[j]
        print(paste(j, suid))
        if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
            fullid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                            "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
        } else { fullid = suid }
        lp = all.regmat[fullid,]
        lind = which(lp > 0)
        lind = lind[order(lp[lind])] # Sort
        lind = tail(lind, 25)
        sqdf = qdf[qdf$uid %in% fullid,c('queryHits','subjectHits','dist')]

        # Assign closest active enhancer for the node:
        snpdf = ldply(lind, function(i){ 
                          x = cdll$cons[[i]] 
                          x = enhmap[x]; 
                          df = sqdf[sqdf$subjectHits %in% x,]
                          df2 = aggregate(dist ~ queryHits, df, min)
                          df = merge(df, df2)
                          df$node = i
                          return(df) }) 
        snpdf$loc = paste0(gwdf$chrom[snpdf$queryHits],"_", 
                           gwdf$chromStart[snpdf$queryHits])
        snpdf$name = enhdf$name[snpdf$subjectHits]

        # Expanded linkdf:
        sldf = unique(merge(linkdf, snpdf))
        nndf = aggregate(id ~ node, unique(sldf[,c('id','node')]), length)
        sldf = unique(merge(all.links, sldf))
        if (nrow(sldf) > 0){
            sldf = aggregate(score ~ id + node + name + gene + nearest + dist + subjectHits + loc, sldf, mean)
            sldf = aggregate(score ~ node + name + gene + nearest + dist + subjectHits + loc, sldf, sum)
            sldf = merge(sldf, nndf)
            sldf$score = sldf$score / sldf$id
            snpdf = merge(snpdf, sldf, all.x=TRUE)
        } else {
            snpdf = merge(snpdf, sldf[,c('node','score','name','gene','nearest','dist','subjectHits','loc')], all.x=TRUE)
        }
        # Merge back into snpdf, add attributes:
        snpdf = merge(snpdf, enhdf[enhdf$name %in% snpdf$name,], all.x=TRUE)
        snpdf$p = gwdf$pValue[snpdf$queryHits]
        snpdf = merge(snpdf, tssdf, all.x=TRUE)
        snpdf$linkdist = with(snpdf, (start + end) / 2 - tss)
        snpdf$nearest = enhdf[snpdf$name, 'nearest']
        snpdf$enr.p = round(lp[snpdf$node], 2)
        snpdf = merge(snpdf, data.frame(node=lind, node.name=leafrep[lind], node.rank=rev(1:length(lind))))
        # Reorder table:
        snpdf = snpdf[order(snpdf$enr.p, decreasing=T),] 
        snpdf = snpdf[order(snpdf$p),] 

        # Rename and add to full table:
        scols = c('node.rank', 'node.name', 'enr.p', 'loc','p', 'dist', 'nearest','symbol','score','linkdist')
        subdf = snpdf[,scols]
        subdf$uid = suid
        alldf = rbind(alldf, subdf)
    }
    save(alldf, file=all.gwlinked.rda)
} else {
    load(all.gwlinked.rda)
}
print(dim(alldf))


# ------------------------------------------------------------------------
# Load in the gene expression matrix, get expression for the linked genes:
# ------------------------------------------------------------------------





# -----------------------------------------------
# Make per-SNP linking plots for a specific GWAS:
# -----------------------------------------------
source(paste0(bindir, 'load_linking_objects.R'))

enhind = enhdf$cls
enhmap = rep(0, max(enhind))
enhmap[enhind] = 1:length(enhind)
j = which(rev(Znam) == '29212778 - Coronary artery disease')
j = which(rev(Znam) == '29892016 - Prostate cancer')
# j = which(rev(Znam) == '26198764 - Schizophrenia')
for (j in 3:length(Znam)) {

    # Want to map to closest active, tissue-specific enhancer
    suid = rev(Znam)[j]
    trait = sub(".* - ", "", suid)
    pmid = sub(" - .*", "", suid)
    traitstr = gsub("'","_", gsub(" ", "_", tolower(trait)))
    traitstrnoparen = gsub("\\)" ,"", gsub("\\(","",traitstr))
    print(paste(j, suid))
    if (suid == "26974007 - Chronic inflammatory diseases (pleiotropy)"){
        fullid = paste0("26974007 - Chronic inflammatory diseases (ankylosing spondylitis, ",
                        "Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)")
    } else { fullid = suid }
    lp = all.regmat[fullid,]
    lind = which(lp > 0)
    lind = lind[order(lp[lind])] # Sort
    lind = tail(lind, 25)
    sqdf = qdf[qdf$uid %in% fullid,c('queryHits','subjectHits','dist')]

    # Assign closest active enhancer for the node:
    snpdf = ldply(lind, function(i){ 
                      x = cdll$cons[[i]] 
                      x = enhmap[x]; 
                      df = sqdf[sqdf$subjectHits %in% x,]
                      df2 = aggregate(dist ~ queryHits, df, min)
                      df = merge(df, df2)
                      df$node = i
                      return(df) }) 

    snpdf$chr = gwdf$chrom[snpdf$queryHits]
    snpdf$pos = gwdf$chromStart[snpdf$queryHits]
    snpdf$loc = paste0(gwdf$chrom[snpdf$queryHits],"_", 
                       gwdf$chromStart[snpdf$queryHits])
    snpdf$name = enhdf$name[snpdf$subjectHits]
    snpdf$enr.p = round(lp[snpdf$node], 2)
    snpdf = merge(snpdf, nodetissue, all.x=TRUE)

    # Get the nearest genes to the GWAS loci:
    sub.gwdf = unique(gwdf[gwdf$uid == suid, c('chrom','chromStart','pValue')])
    sub.gwdf = sub.gwdf[order(sub.gwdf$chrom),]
    sub.gwdf$chr = paste0('chr', sub.gwdf$chrom)
    locgr = with(sub.gwdf, GRanges(paste0('chr', chrom), IRanges(chromStart, chromStart)))
    sub.gwdf$gene = as.character(tssgr$name[nearest(locgr, tssgr)])
    sub.gwdf = sub.gwdf[order(sub.gwdf$chr),]
    center.genes = unique(sub.gwdf$gene)

    # Also determine top tissues for this GWAS:
    subdf = alldf[alldf$uid %in% suid,]
    ndf = unique(snpdf[,c('node','GROUP','enr.p')])
    ndf = aggregate(enr.p ~ GROUP, ndf, sum)
    ndf = ndf[order(ndf$enr.p, decreasing=T),]
    ndf = ndf[!(ndf$GROUP %in% c('Multiple', 'Cancer', 'Other')),]
    add.groups = head(ndf$GROUP, 10)

    ldf = ndf[ndf$enr.p >= max(ndf$enr.p) /2 ,]
    link.groups = head(ldf$GROUP,4)

    if (length(add.groups) < 10){
        add.groups = unique(c(add.groups, 'HSC & B-cell','Muscle','Liver','Brain','Kidney',
                              'Pancreas','Heart', 'Stromal', 'ESC', 'Blood & T-cell'))[1:10]
    }

    # Make UID specific directory for images:
    traitdir = paste0(gwasdir, traitstrnoparen, '_', pmid, '/')
    cmd = paste('mkdir -p', traitdir)
    system(cmd)

    sldf = unique(linkdf[linkdf$node %in% snpdf$node,c('node','id')])

    # Names that we should pull down:
    window = 2e6
    genetss <- with(tssdf[tssdf$gene %in% center.genes,], GRanges(seqnames=chr, ranges=IRanges(start=tss - window/2, end=tss + window/2)))
    locgr = with(locdf, GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), name=name))
    uidovldf = as.data.frame(findOverlaps(genetss, locgr))
    enam = locgr$name[unique(uidovldf$subjectHits)]
    tssovldf = as.data.frame(findOverlaps(genetss, tssgr))
    egene = tssgr$name[unique(tssovldf$subjectHits)]

    # TODO: HOW AND WHAT TO PLOT HERE.

    # TODO: Need to pull down + aggregate ALL the predicted links 
    # in the neighborhood of the enhancers/SNPs we care about
    sub.enhdf = enhdf[enhdf$name %in% enam,]

    print("[STATUS] Getting group-relevant links:")
    brdf = c()
    for (group in link.groups){
        print(group)
        groupstr = tolower(gsub("__", "", str_replace_all(group, "[^[:alnum:]]", "_")))
        sub.group.rda = paste0('linking_data/links_by_group.', groupstr, '.Rda')
        load(sub.group.rda)
        gdf = merge(gdf, sub.enhdf)
        brdf = rbind(brdf, gdf)
    }
    brdf$GROUP = brdf$group
    
    # --------------------------------
    # For each SNP, plot neighborhood:
    # --------------------------------
    # TODO: Make sure all center.genes work...
    # TODO: Make sure the browser ACTUALLY CONTAINS THE SNP:
    # gene: ADGRL2
    center.genes = center.genes[center.genes %in% anno$gene]
    cpairsdf = NULL
    for (k in 1:length(center.genes)){
        # Use locgr:
        # k = which(sapply(center.genes, function(x){as.character(anno$symbol[anno$gene == x])}) == 'ADGRG1')
        # k = which(sapply(center.genes, function(x){as.character(anno$symbol[anno$gene == x])}) == 'IGSF9B')
        # k = which(sapply(center.genes, function(x){as.character(anno$symbol[anno$gene == x])}) == 'RPGRIP1L')

        ensg = center.genes[k]
        gene = as.character(anno$symbol[anno$gene == ensg])
        print(paste(k, gene))
        infoline = anno[anno$symbol == gene,]
        chrom = infoline$chr

        # Get overlapping genes:
        genetss <- with(infoline, GRanges(seqnames=chr, ranges=IRanges(start=tss - window/2, end=tss + window/2)))

        tssgr = GRanges(seqnames=tssdf$chr, IRanges(tssdf$tss, tssdf$tss), name=tssdf$gene)
        anno.ovl = data.frame(findOverlaps(genetss, tssgr))
        kept.gene = as.character(tssgr$name[anno.ovl$subjectHits])
        gloc = tssdf$tss[anno.ovl$subjectHits]
        gnames = as.character(sapply(kept.gene, function(x){genemap$symbol[genemap$gene == x]}))
        ord = order(gloc) 
        gnames = gnames[ord]
        kept.gene = kept.gene[ord]
        gloc = gloc[ord]
        names(gloc) = kept.gene

        # Load tested pairs:
        if (is.null(cpairsdf) || cpairsdf$chrom[1] != chrom){
            print(paste("Loading links for", chrom))
            cfile = paste0(lddir, chrom, '_pairs_ENH_ovl_df.tsv.gz')
            # genes = sub.gwdf$gene[sub.gwdf$chr == chrom]
            cpairsdf = read.delim(cfile, header=T, sep="\t", stringsAsFactors=F)
            cpairsdf$mind = cpairsdf$mind + 1  # Was 0-indexed
            cpairsdf$chrom = chrom
        }
        sub.ovldf = cpairsdf[cpairsdf$gene %in% kept.gene,]
        sub.ovldf$chr = chrom
        sub.ovldf$corr = h1[sub.ovldf$ind]
        sub.ovldf$corr_k4 = h2[sub.ovldf$ind]

        # Get all links:
        sub.ovldf$name = locdf$name[sub.ovldf$mind]
        sub.ovldf$mid = locdf$mid[sub.ovldf$mind]
        sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]
        sub.ovldf$dist = sub.ovldf$mid - infoline$tss 
        sub.ovldf = sub.ovldf[abs(sub.ovldf$dist) < (window / 2),]
        sub.ovldf = sub.ovldf[order(sub.ovldf$mid),]
        # Kept locations:
        kept.locdf = unique(sub.ovldf[,c('mid','mind', 'name')])
        kept.dhs = kept.locdf$mind
        kept.mid = kept.locdf$mid
        rownames(kept.locdf) = kept.locdf$name
        kept.locdf$id = 1:nrow(kept.locdf)

        # Slice locations from hdf5:
        h5f = H5Fopen(H3K27ac.file)
        h5d = h5f&"matrix"
        sub.mmat = h5d[,matind[kept.dhs]]
        H5Dclose(h5d)
        H5Fclose(h5f)
        sub.mmat = t(sub.mmat)
        colnames(sub.mmat) = mnames

        # Make average dhs tracks:
        tform = make.tform(meta[mnames,'GROUP'], norm=TRUE, u=odf$GROUP)
        sub.avg.mmat = sub.mmat %*% tform 

        # Average expression:
        tform = make.tform(meta[colnames(gmat),'GROUP'], norm=TRUE, u=odf$GROUP)
        sub.gmat = gmat[kept.gene,]
        sub.avg.expr = sub.gmat %*% tform
        sub.avg.expr[is.na(sub.avg.expr)] = 0

        # Color each enhancer (will be to link) by only H3K27ac
        group.max = apply(sub.avg.mmat, 1, which.max)
        group.col = odf$COLOR[group.max]
        group.nam = odf$GROUP[group.max]

        # TODO: Experiment with H3K27ac or H3K4me1 or combined scores
        THRESHOLD = 0.4
        # plot(sub.ovldf$corr, sub.ovldf$corr_k4, pch=19,
        #      col=ifelse(sub.ovldf$corr < THRESHOLD, 'grey50','indianred'))
        # abline(h=0); abline(v=0)
        sublinkdf = sub.ovldf[sub.ovldf$corr > THRESHOLD, c('gene','name', 'mid', 'corr')]
        sublinkdf$gene = as.character(sublinkdf$gene)
        sublinkdf$tss = gloc[sublinkdf$gene]
        # sublinkdf$mid = kept.locdf[sublinkdf$name, 'mid']
        sublinkdf$id = kept.locdf[sublinkdf$name, 'id']
        sublinkdf$color = group.col[sublinkdf$id]
        sublinkdf$nam = group.nam[sublinkdf$id]
        sublinkdf$symbol = as.character(sapply(sublinkdf$gene, function(x){genemap$symbol[genemap$gene == x]}))
        # write.table(sublinkdf, paste0('linking_data/corr_links/gene_links_rawcorr_H3K27ac_', gene, '_w', window,'.tsv'), sep="\t", quote=F, row.names=F)

        # TODO: FIRST PASS PLOT THE KEPT LINKS with intersections:
        # if (suid == '26198764 - Schizophrenia' || suid == '29212778 - Coronary artery disease'){
        use.grouponly = TRUE
        # } else { use.grouponly = FALSE }
        # use.grouponly = FALSE

        # TODO: Return plots by two different score thresholds:
        if (use.grouponly){
            preddf = brdf[brdf$gene %in% kept.gene,]
            preddf = preddf[preddf$score > 0.01,]
            # preddf$GROUP = 'Brain'
            preddf = merge(preddf, odf)
        } else {
            preddf = all.links[all.links$gene %in% kept.gene,]
            preddf = merge(preddf, meta[,c('id','COLOR','GROUP')])
        }
        # TODO: Reduce properly:
        if (nrow(preddf) > 0){
            npdf = aggregate(score ~ GROUP + COLOR + gene + name, preddf, length)
            npdf = merge(npdf, locdf)
            npdf$mid = (npdf$start + npdf$end) / 2
            npdf = merge(npdf, tssdf)
            npdf$color = as.character(npdf$COLOR)
            npdf$nam = as.character(npdf$GROUP)
            # TODO: TESTING PLOT LIKE THIS:
            sublinkdf = npdf
        } else {
            sublinkdf = preddf[,c('score','GROUP','COLOR','gene','name')] 
        }

        plot.red = FALSE
        if (plot.red){
            # sublinkdf = sublinkdf[sublinkdf$nam %in% add.groups,]
            if (j == 386){
                sublinkdf = sublinkdf[sublinkdf$nam %in% c('Cancer','Reproductive'),]
            } else if (suid == '26198764 - Schizophrenia'){
                sublinkdf = sublinkdf[sublinkdf$nam %in% c('Brain'),]
            } else {
                sublinkdf = sublinkdf[sublinkdf$nam %in% c('Liver','Heart','Endothelial'),]
            }
        }

        if (nrow(sublinkdf) > 0){
            # Calc link parameters:
            link.center = with(sublinkdf, (tss + mid) / 2)
            link.radius = abs(link.center - sublinkdf$tss)
        }
        # GWAS SNPs:
        gwloc = sub.gwdf$chromStart[sub.gwdf$chr == chrom]
        # TODO: Find linked SNPs - for each gene (using gwdf)

        # TODO: FINE MAP each SNP to the closest RELEVANT ENHANCER 
        # Set up the data limits:
        xlim = with(infoline, c(tss - window /2 , tss + window / 2))
        # ylim.atac = c(0, max(sub.avg.mmat, na.rm=T))
        ylim.atac = c(0, 20)

        # To convolve signal slightly
        # yconv = c(0.01, .05, .15, .30, .4, .30, .15, .05, 0.01) / .4
        # xconv = c(-200, -150, -100, -50, 0, 50, 100, 150, 200)
        # For axis:
        step = ((window * 2) / 5)
        mx = (floor(xlim[1] / step) + 1) * step
        atpar = seq(mx, xlim[2], by=step)
        lbls = comma_format()(atpar)

        add.genes = TRUE
        if (add.genes){
            genesuf = '_with_tx'
            sub.txdf = aggregate(txid ~ ., txdf[txdf$symbol %in% gnames,c('chr','start','end','strand','txid', 'symbol'),], function(x){x[1]})
            # Pad for name or other
            sub.txdf$pad = 40000
            sub.txdf = within(sub.txdf, pad.start <- start - pad * (strand == '+'))
            sub.txdf = within(sub.txdf, pad.end <- end + pad * (strand == '-'))
            sub.txdf = sub.txdf[order(sub.txdf$pad.start),]
            # Figure out how many lines we need:
            stgr = with(sub.txdf, GRanges(seqnames=chr, IRanges(pad.start, pad.end)))
            cvg = coverage(stgr)
            nlines = max(cvg)
            genes.height = .5 * nlines + .1
            # Order of genes:
            tx.ovl = data.frame(findOverlaps(stgr, stgr))
            tx.ovl = tx.ovl[tx.ovl$queryHits != tx.ovl$subjectHits,]
            # tx.ovl$qn = sub.txdf$symbol[tx.ovl$queryHits]
            NGENES = nrow(sub.txdf)
            sub.txdf$line = rep(0, NGENES)
            poss.lines = rev(1:nlines)
            for (i in 1:NGENES){
                if (i %in% tx.ovl$queryHits){
                    pdf = tx.ovl[tx.ovl$queryHits == i,]
                    comp = sort(pdf[,2])
                    notline = c()
                    for (j in comp){
                        if (i > j){
                            notline = c(notline, sub.txdf$line[j])
                        }
                    }
                    if (length(notline) > 0){
                        l = poss.lines[!(poss.lines %in% notline)][1]
                    } else { l = nlines }
                    sub.txdf$line[i] = l
                } else {
                    sub.txdf$line[i] = 1
                }
            }
            sub.exdf = merge(exdf, sub.txdf[,c('txid','symbol', 'line')])
            sub.uxdf = merge(uxdf, sub.txdf[,c('txid','symbol', 'line')])
        } else { genesuf = ''; genes.height = 0}

        if (plot.red){
            # add.groups = c('Heart','Sm. Muscle','Endothelial','HSC & B-cell')
            genesuf = paste0(genesuf, '_reduced')
        }
        if (use.grouponly){
            genesuf = paste0(genesuf, '_grouponly')
        }

        NCT = length(add.groups)

        top.arc = TRUE
        plot.width=8
        plot.height=1 + 3 * (NCT / 5) + 1 * add.genes * genes.height / 1.5
        plot.height=2 / 1.2
        png(paste0(traitdir, 'gene_link_vis_', gene, '_w', window, genesuf, '_testing.png'), units='in', res=450, width=plot.width, height=plot.height)
        # arc.sect=3
        arc.sect=2.5
        if (top.arc){
            lmat = matrix(c(NCT+2 + 1 * add.genes, 1:(NCT + 1 + add.genes * 1)), nrow=NCT + 2 + 1 * add.genes, ncol=1)
            heights = c(arc.sect, rep(1, NCT + 1))
        } else {
            lmat = matrix(1:(NCT + 1), nrow=NCT + 1 + 1 * add.genes, ncol=1)
            heights = c(rep(1, NCT), arc.sect)
        }
        if (add.genes) { heights = c(heights, genes.height) }
        layout(lmat, heights=heights)
        par(xaxs='i')
        sp=0
        par(mar=rep(sp, 4))
        for (i in 1:NCT){ 
            sgroup = add.groups[i]
            x = kept.locdf$mid
            y = c(sub.avg.mmat[,sgroup])
            ylim = ylim.atac
            plot(x, y, type='n', axes=F, ylab='', xlab='', ylim=ylim, xlim=xlim)
            # Add highlights:
            if (i == NCT){ 
                axis(1, at=atpar, labels=rep('', length(atpar)), lwd=.25, tck=-0.15)
                text(x=atpar, y=parpos(2, .55 - .2 * plot.red), labels=paste0(chrom, ":", lbls), xpd=NA, cex=.6)
            }
            # Add the highlighted linked locations:
            tol = 500
            if (nrow(sublinkdf) > 0){
                link.atac = sublinkdf$mid
                rect(xleft=link.atac - tol, xright=link.atac + tol, 
                     ybottom=ylim.atac[1], ytop=ylim.atac[2], col='grey90', border=NA)
            }
            # abline(v=gwloc, col='slateblue', lwd=.5, lty='dashed')
            diffpar = seq(-10, 10, by=2.5) * 1e5 
            abline(v=infoline$tss + diffpar, lwd=.5, col='grey50', lty='dashed')
            abline(v=gloc, lwd=.5, col='red', lty='dashed')
            # Plot actual track:
            rect(xleft=x-tol, ybottom=0, xright=x+tol, ytop=y, col=colvals$group[sgroup], border=NA)
            groupcex = .8
            text(x=xlim[1] + 0.028 * diff(xlim), y=mean(ylim.atac), sgroup, 
                 font=2, cex=groupcex, col=colvals$group[sgroup], adj=0)
            box(lwd=.25)
            arrowcol='grey25'
            if (i == 1){
                if (!(add.genes)) {
                    text(gloc + -1 * diff(xlim) * 1e-3, parpos(2, -.4),
                         labels=gnames, adj=1, col=arrowcol, font=2, cex=.9)
                }
                for (dd in diffpar){
                    tstr = paste(dd / 1e6,'Mb')
                    adj = 1.2 * (dd > 0) - .1
                    if (dd == 0){ tstr = 'TSS'; adj=-.05 }
                    text(x=infoline$tss + dd, y = ylim[1] + .8 * (diff(ylim.atac)),
                         labels=tstr, cex=.5, adj=adj, col='grey30')
                }
            }
        }
        # TADs and SNPs on the bottom:
        par(mar=rep(sp, 4))
        plot(0, 1, xlim = xlim, ylim=c(0,1), type='n', axes=F)
        # Plot snps:
        gwind = which(gwloc >= infoline$tss - window/2 & gwloc <= infoline$tss + window/2)
        if (length(gwind) >  0){
            kept.gwloc = gwloc[gwind]
        } else {
            kept.gwloc = gwloc[gwloc >= infoline$tss - window]
        }
        rect(xleft=kept.gwloc-tol * 2, xright=kept.gwloc + tol * 2, 
             ybottom=0.8, ytop=1, border=NA, col='slateblue', lwd=1)
        tcex = 0.6
        tlab = paste(trait, 'GWAS lead SNPs')
        twidth = strwidth(tlab, cex=tcex)
        tloc = xlim[1] + .01 * diff(xlim)
        text(tloc, y=.15,labels=tlab,
             font=2, cex=tcex, col='slateblue', adj=0)
        min.snp = min(kept.gwloc) # First possible snp.
        max.snp = max(kept.gwloc) # First possible snp.
        # Bottom:
        yloc = c(0.1, 0.75)
        segments(x0=tloc + twidth + .008 *  diff(xlim),
                 x1=max.snp - 0.005 *  diff(xlim),
                 y0=yloc[1],y1=yloc[1], col='slateblue', lwd=.5, xpd=TRUE)
        segments(x0=kept.gwloc - 0.005 *  diff(xlim), x1=kept.gwloc,
                 y0=yloc[1],y1=yloc[2], col='slateblue', lwd=.5, xpd=TRUE)
        segments(x0=kept.gwloc - 0.001 * diff(xlim),
                 x1=kept.gwloc + 0.001 * diff(xlim),
                 y0=yloc[2],y1=yloc[2], col='slateblue', lwd=.5, xpd=TRUE)
        if (add.genes){
            par(mar=rep(sp, 4))
            plot(0, 1, xlim = xlim, ylim=c(.5,nlines + .5), type='n', axes=F)
            # box(lwd=.5)
            genecol = 'midnightblue'
            tpad = 1000
            for (i in 1:nrow(sub.txdf)){
                genecex = 0.45
                with(sub.txdf[i,], text(x=(start - tpad) * (strand == '+') + (end + tpad) * (strand == '-'),
                                        y=line, labels=symbol, cex=genecex, xpd=TRUE, adj=(strand == '+') * 1))
                # segments(x0=sub.txdf$start[i], x1=sub.txdf$end[i],
                #          y0=sub.txdf$line[i], y1=sub.txdf$line[i], 
                #          lwd=.25, col='slateblue')
                start = with(sub.txdf[i,], start * (strand == '+') + end * (strand == '-'))
                end = with(sub.txdf[i,], start * (strand == '-') + end * (strand == '+'))
                multarrows(x0=start, x1=end,
                           y0=sub.txdf$line[i], y1=sub.txdf$line[i], 
                           n_arr = round((sub.txdf$end[i] - sub.txdf$start[i]) / 5000),
                           length=.02,
                           lwd=.25, col='midnightblue')
            }
            lpad =0.2
            # Add exons:
            rect(xleft=sub.exdf$start, xright=sub.exdf$end,
                 ybottom=sub.exdf$line - lpad, ytop=sub.exdf$line + lpad, 
                 lwd=.25, col='midnightblue', border=NA)
            rect(xleft=sub.uxdf$start, xright=sub.uxdf$end,
                 ybottom=sub.uxdf$line - lpad /2, ytop=sub.uxdf$line + lpad /2, 
                 lwd=.25, col='midnightblue', border=NA)
        }
        # TODO: Add rsid
        if (nrow(sublinkdf) > 0){
            if (top.arc){ 
                par(mar=rep(sp, 4))
                plot(0, 1, xlim = xlim, ylim=c(0,1), type='n', axes=F)
            }
            asp.ratio = plot.width / (plot.height * (arc.sect / (NCT + add.genes * genes.height + top.arc + arc.sect)))
            mod.link.radius = abs(link.radius) / cos(pi/4)
            yvar = mod.link.radius / diff(xlim) * asp.ratio / 1.03
            # yvar = mod.link.radius / diff(xlim) * asp.ratio / (1.03 + 0.04 * add.genes * (genes.height / 1.5)**2)
            # if (top.arc){ ang1 = pi/4 }
            draw.arc(x=link.center, y= 1 * (!top.arc) + (2 * (!top.arc) - 1) * (yvar) * (pi/4), xpd=TRUE, 
                     lwd=(sublinkdf$gene == ensg) * .5 + .5,
                     radius=mod.link.radius, angle1=pi * (!top.arc) + pi/4,angle2=pi * (!top.arc)+ pi - pi/4, col=sublinkdf$col)
        }
        dev.off()

    }

}





