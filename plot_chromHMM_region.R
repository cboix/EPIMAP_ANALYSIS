#!/usr/bin/R
# -----------------------------------------------
# Given a matrix of ChromHMM states (loc by file)
# Plot the region:
# -----------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))

options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")
plot.extra = FALSE # Only if we want to plot extra, not-fixed order views

# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need calls dataframe filename.")
} else {        
    filename = args[1]
    tagline = args[2]
    if(length(args) > 2){
        rdfilename = args[3]
        plot.rdmap = TRUE
    } else { plot.rdmap = FALSE }
}

# filename = '~/extract_region_25364/imputed_region_calls.bed'
# tagline = 'imputed_observed_aux_18_on_mixed_impobs_QCUT'
# rdfilename = '~/extract_region_25364/roadmap_region_calls.bed'
# plot.rdmap = TRUE

# Directories:
imgpref = paste0(img, "ChromHMM/")
chmmdir = 'ChromHMM/'
system(paste('mkdir -p', imgpref))

# Read in calls:
df = read.delim(filename, header=F)
names(df) <- c('chr','start','end','BSSID','stateid')
range = c(min(df$start), max(df$end))
df$state = as.numeric(sub("E","", df$stateid))
nstates = max(df$state)
regionline = paste(df$chr[1], range[1], range[2], sep="_")
infoline = paste0(regionline,'_',tagline,"_n",nstates)
# From load_metadata:
df = df[df$BSSID %in% cellorder,]

if (nstates == 15){
    model = 'observed_15' 
} else if (nstates == 18){
    model = 'observed_aux_18'
} else {
    model = 'imputed12_25'
}

# Get state colors:
dcol <- read.table(paste0('CHMM_', model, '_colors.tsv'),header=F)
names(dcol) <- c('state','name','color')
dcol$col <- sapply(dcol$color,function(x){
                       s <- unlist(strsplit(as.character(x),split=','))
                       s <- as.numeric(s)
                       rgb(s[1],s[2],s[3],max=255) })

# Set up matrix for plotting:
cells = as.character(unique(df$BSSID))
mat = matrix(0, nrow=diff(range)/200, ncol=length(cells))
colnames(mat) <- cells

# Fill in matrix from dataframe
df$start = (df$start - range[1]) / 200 + 1
df$end = (df$end - range[1]) / 200
df$BSSID = factor(as.character(df$BSSID), levels=cells)
df$col = as.numeric(df$BSSID)
print("Populating matrix")
for (i in 1:nrow(df)){
    mat[df$start[i]:df$end[i], df$col[i]] = df$state[i]
}

# Check that everything was filled in:
print("Checking no values are left as 0:")
print(which(mat == 0))

# Matrix element of plot:
mat.element = function(mat, df, add.names=FALSE, axlab=TRUE, add.axis=TRUE,
                       margins=NULL, noleft=FALSE, noright=FALSE, title=TRUE){
    if (is.null(margins)){
        par(mar=c(3,(.5 + 2.5 * axlab + 6* add.names) * (!noleft), 2, .5 * (!noright)))
    } else { par(mar=margins) }
    image(mat, axes=FALSE, col=dcol$col, useRaster=TRUE)
    if (add.names){
        text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]), labels=colnames(mat), srt=0, adj=1, xpd=TRUE,cex=.2)
    }
    if (add.axis){
        # Add location:
        sq = seq(range[1], range[2], length.out=6)
        loc = (sq - range[1]) / (200 * nrow(mat))
        axis(side=1, at=loc, labels=FALSE)
        axislabels = as.character(sapply(sq / 10^6, function(x){sprintf("%0.2f",x)}))
        text(x=loc, y=par()$usr[3] - 0.025 * diff(par()$usr[3:4]), labels=axislabels, xpd=NA, cex=1)
    }
    ns = max(df$state)
    if(title){
        mtext(paste0(df$chr[1],":", range[1],"-", range[2],' (',ns,' states)'),side=3, line=.25, cex=1.5)
        mtext(paste0(df$chr[1], ' Location (Mbp)'),side=1, line=2, cex=.8)
    }
    if (axlab){
        mtext(paste0('Cell Types/States'), side=2, line= 1 + 6 * add.names, cex=1.5)
    }
}

# ------------------------------
# Get matrix in the fixed order:
# ------------------------------
print("Ordering cells by fixed order")
ord = read.delim('Annotation/bssid_order_frz20190326.tsv', header=F, stringsAsFactors=F)
cord = ord[,1]
cord = cord[cord %in% cellorder]
cord = cord[cord %in% colnames(mat)]
reord = mat[, cord]

labels = meta[as.character(cord), 'GROUP']
faclabels = as.matrix(as.numeric(labels))
colset = as.character(rdcol$COLOR)
# Label runs (that are not NONE):
lablist = label.runs(faclabels, labels, rdcol)

# Plot matrix:
print("Plotting reordered matrix + colors.")
png(paste0(imgpref, 'chromHMM_region_', infoline, '_col_fixord.png'), res=300, width=8, height=11, units='in')
layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
par(yaxs="i")
par(xaxs="i")
par(mar=c(3,5,2,.25))
image(t(faclabels), axes=F, col=colset, useRaster=T)
abline(h=c(0,1),lty=1,lw=0.5)
abline(v=c(-1,1),lty=1,lw=0.5)
text(y=lablist[[1]],
     x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
     labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
text(x=0, labels='Type/Tissue\n(As Roadmap)',
     y=par()$usr[1]-0.01 * (par()$usr[2]-par()$usr[1]), 
     srt=0, adj=1, xpd=TRUE,cex=.65, col='black')
mat.element(reord, df, axlab=FALSE)
dev.off()


# Plot roadmap:
if (plot.rdmap){
    print("Getting roadmap matrix.")
    # Read in calls:
    rdf = read.delim(rdfilename, header=F)
    names(rdf) <- c('chr','start','end','eid','stateid')
    rdf$state = as.numeric(sub("E","", rdf$stateid))
    if (nstates != max(df$state)){
        print("[ERROR] Roadmap is not is same model")
    }
    eidlist = as.character(unique(rdf$eid))
    print(paste("There are", length(unique(rdf$eid)), "Roadmap epigenomes"))
    eidlist[!(eidlist %in% rd.matching$eid)]
    rdf = merge(rdf, rd.matching)

    # Fill in matrix from dataframe
    rdf$start = (rdf$start - range[1]) / 200 + 1
    rdf$end = (rdf$end - range[1]) / 200
    # Set up matrix for plotting:
    rdcells = as.character(unique(rdf$id))
    print(paste("There are", length(rdcells), "matched epigenomes"))
    rmat = matrix(nstates, nrow=diff(range)/200, ncol=length(cells))
    colnames(rmat) <- cord
    # map:
    rdf$id = factor(as.character(rdf$id), levels=colnames(rmat))
    rdf$col = as.numeric(rdf$id)
    print("Populating matrix")
    for (i in 1:nrow(rdf)){
        rmat[rdf$start[i]:rdf$end[i], rdf$col[i]] = rdf$state[i]
    }

    # Plot matrix:
    print("Plotting reordered matrix with ROADMAP - staggered")
    png(paste0(imgpref, 'chromHMM_region_', infoline, '_withrdmap_col_fixord.png'), res=300, width=10, height=11, units='in')
    outer=1.75
    inner=10
    layout(matrix(c(1:3),1,3), widths=c(inner, inner, outer), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    mat.element(rmat, rdf, axlab=FALSE)
    mat.element(reord, df, axlab=FALSE)
    # Plot labels on right
    par(mar=c(3,0,2,5))
    image(t(faclabels), axes=F, col=colset,useRaster=TRUE)
    abline(h=c(0,1),lty=1,lw=0.5)
    abline(v=c(-1,1),lty=1,lw=0.5)
    text(y=lablist[[1]],
         x=par()$usr[2]+0.1*(par()$usr[2]-par()$usr[1]), 
         labels=lablist[[2]], srt=0, adj=0, xpd=TRUE,cex=.8, col=lablist[[3]])
    dev.off()

    # ----------------------------------
    # Make matrix for plotting together:
    # ----------------------------------
    # Center:
    rdcord = cord[cord %in% rdcells]
    notrdcell = cord[!(cord %in% rdcells)]
    NNOT = length(notrdcell)
    # Make collapsed rdmap:
    crmat = matrix(nstates, nrow=diff(range)/200, ncol=length(cells))
    colnames(crmat) <- c(notrdcell[1:(round(NNOT/2))], rdcord, notrdcell[(round(NNOT/2) + 1):NNOT])
    # Map:
    rdf$id = factor(as.character(rdf$id), levels=colnames(crmat))
    rdf$ccol = as.numeric(rdf$id)
    # Fill in matrix from dataframe
    print("Populating matrix")
    for (i in 1:nrow(rdf)){
        crmat[rdf$start[i]:rdf$end[i], rdf$ccol[i]] = rdf$state[i]
    }

    # Lines are:
    linedf = unique(rdf[,c('col','ccol')])

    # Plot matrix:
    print("Plotting reordered matrix with ROADMAP - with lines")
    png(paste0(imgpref, 'chromHMM_region_', infoline, '_collrdmap_col_fixord.png'), res=300, width=11, height=11, units='in')
    outer=1.75
    inner=10
    layout(matrix(c(1:4),1,4), widths=c(inner, outer * 1.5, inner, outer), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    mat.element(crmat, rdf, axlab=FALSE, noright=FALSE)
    # Make segments
    par(mar=c(3,0,2,0))
    plot(1,1, type='n', xlim=c(0,1), ylim=c(0,length(cells)),axes=F)
    segments(rep(0, nrow(linedf)), linedf$ccol, rep(1,nrow(linedf)), linedf$col, col='gray')
    # Plot the full imputed:
    mat.element(reord, df, axlab=FALSE, noleft=F)
    # Plot labels on right
    par(mar=c(3,0,2,5))
    image(t(faclabels), axes=F, col=colset)
    abline(h=c(0,1),lty=1,lw=0.5)
    abline(v=c(-1,1),lty=1,lw=0.5)
    text(y=lablist[[1]],
         x=par()$usr[2]+0.1*(par()$usr[2]-par()$usr[1]), 
         labels=lablist[[2]], srt=0, adj=0, xpd=TRUE,cex=.8, col=lablist[[3]])
    dev.off()

    # ------------- PLOT ACCORDING TO COLORS: --------------------
    odf = rdcol  # Copy don't change rdcol
    odf = rbind(odf, c('None','grey80','Other'))
    clvs = c('Tissue', 'Primary Cell', 'Immune', 'Neural', 'Stem-like', 'Other')
    odf$category = factor(odf$category, levels=clvs)
    # Sort alpha then group:
    odf = odf[order(odf$GROUP),]
    odf = odf[order(odf$category),]

    # 1. Establish new order:
    # First by rdmp colors then by fixord:
    labels = factor(labels, levels=odf$GROUP)
    ord2 = order(labels, decreasing=TRUE)
    cord2 = cord[ord2]
    reord2 = mat[,cord2]
    labels2 = meta[as.character(cord2), 'GROUP']
    colcut = as.numeric(labels2)
    col.breaks = calc.breaks.acut(colcut)
    faclabels2 = as.matrix(colcut)
    colcut.nam = colcut
    names(colcut.nam) = cord2
    # Loc labels. Adjust vals if very close to each other.
    ldf = data.frame(GROUP=labels2, quant=1:length(labels2) * 1.0/ length(labels2) )
    labeldf = aggregate(quant ~ GROUP, ldf, mean)
    labeldf = merge(labeldf, odf)
    lord = order(labeldf$quant, decreasing=FALSE)
    labeldf = labeldf[lord,]
    rmat2 = rmat[,cord2]

    # Plot matrix:
    print("Plotting reordered matrix with ROADMAP - staggered")
    png(paste0(imgpref, 'chromHMM_region_', infoline, '_withrdmap_byrdcol_fixord.png'), res=300, width=10, height=11, units='in')
    outer=1.75
    inner=10
    layout(matrix(c(1:3),1,3), widths=c(inner, inner, outer), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    mat.element(rmat2, rdf, axlab=FALSE)
    mat.element(reord2, df, axlab=FALSE)
    # Plot labels on right
    par(mar=c(3,0,2,5))
    image(t(faclabels2), axes=F, col=colset)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    text(y=labeldf$quant,
         x=par()$usr[2]+0.1*(par()$usr[2]-par()$usr[1]), 
         labels=labeldf$label, srt=0, adj=0, xpd=TRUE,cex=.7, col=colset)
    dev.off()

    # ----------------------------------
    # Make matrix for plotting together:
    # ----------------------------------
    # Center:
    rdcord2 = cord2[cord2 %in% rdcells]
    notrdcell = cord2[!(cord2 %in% rdcells)]
    NNOT = length(notrdcell)
    # Make collapsed rdmap:
    crmat2 = matrix(nstates, nrow=diff(range)/200, ncol=length(cells))
    # colnames(crmat2) <- c(notrdcell[1:(round(NNOT/2))], rdcord2, notrdcell[(round(NNOT/2) + 1):NNOT])
    colnames(crmat2) <- c(notrdcell, rdcord2)
    # Map:
    rdf$id = factor(as.character(rdf$id), levels=colnames(crmat2))
    rdf$ccol2 = as.numeric(rdf$id)
    # Fill in matrix from dataframe
    print("Populating matrix")
    for (i in 1:nrow(rdf)){
        crmat2[rdf$start[i]:rdf$end[i], rdf$ccol2[i]] = rdf$state[i]
    }

    # Corresponds to rdcord2
    linedf = data.frame(rdord = sapply(rdcord2, function(x){which(colnames(crmat2) == x)}),
                        neword = sapply(rdcord2, function(x){which(cord2 == x)}))
    linedf$COLOR = COLOR=meta[rdcord2, 'COLOR']
    use.bmart=FALSE
    if(use.bmart){
        library(Gviz)
        library(GenomicRanges)
        library(biomaRt)
        bm <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
        itrack <- IdeogramTrack(genome='hg19', chromosome='chr19', name='chr19', cex=.7)
    }

    # Plot matrix:
    print("Plotting reordered matrix with ROADMAP - with lines")
    png(paste0(imgpref, 'chromHMM_region_', infoline, '_collrdmap_byrdcol_fixord.png'), res=300, width=7, height=7, units='in')
    # pdf(paste0(imgpref, 'chromHMM_region_', infoline, '_collrdmap_byrdcol_fixord.pdf'), width=7, height=7)
    sp=0.15
    bsp=1.5
    outer=4
    inner=10
    heights=c(6,1)
    widths=c(inner, outer, inner, outer * 1.2)
    layout(matrix(c(1:8),2,4, byrow=TRUE), heights=heights, widths=widths)
    par(yaxs="i")
    par(xaxs="i")
    mat.element(crmat2, rdf, axlab=FALSE, noright=TRUE, noleft=TRUE, 
                add.axis=FALSE, title=FALSE, margins=c(3,0,sp,0))
    # Make segments
    par(mar=c(bsp,0,sp,0))
    plot(1,1, type='n', xlim=c(0,1), ylim=c(0,length(cells)),axes=F)
    rx = c(0.01, 0.06, 0.94, 0.99)
    xx = par()$usr[1] + rx * diff(par()$usr[1:2])
    NL = nrow(linedf)
    segments(x0=rep(xx[1], NL), y0=linedf$rdord, lwd=0.5,
             x1=rep(xx[2], NL), y1=linedf$rdord, col=linedf$COLOR)
    segments(x0=rep(xx[2], NL), y0=linedf$rdord, lwd=0.5,
             x1=rep(xx[3], NL), y1=linedf$neword, col=linedf$COLOR)
    segments(x0=rep(xx[3], NL), y0=linedf$neword, lwd=0.5,
             x1=rep(xx[4], NL), y1=linedf$neword, col=linedf$COLOR)
    # Plot the full imputed:
    mat.element(reord2, df, axlab=FALSE, noleft=TRUE, noright=TRUE, 
                title=FALSE, margins=c(bsp,0,sp,0))
    box(lwd=.25)
    # Plot labels on right
    par(mar=c(bsp,sp,sp,7.7))
    image(t(faclabels2), axes=F, col=colset, useRaster=TRUE)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    box.pad = 0.018
    yy = labeldf$quant
    xx = space.1d(yy, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    xx = space.1d(xx, box.pad=box.pad, lim=c(0 + 1e-3, 1 - 1e-3))
    rx = c(0.1, 0.3, 0.7, 0.9, 1)
    x = par()$usr[2] + rx * diff(par()$usr[1:2])
    text(y=xx, x=x[5], labels=labeldf$GROUP, srt=0, adj=0,
         xpd=TRUE, cex=1, col=labeldf$COLOR)
    par(xpd=TRUE)
    segments(x0=x[1], y0=yy, x1=x[2], y1=yy, col=labeldf$COLOR, lwd=.5)
    segments(x0=x[2], y0=yy, x1=x[3], y1=xx, col=labeldf$COLOR, lwd=.5)
    segments(x0=x[3], y0=xx, x1=x[4], y1=xx, col=labeldf$COLOR, lwd=.5)
    par(xpd=FALSE)
    text(x=x[1], y=par()$usr[3] - 0.025 * diff(par()$usr[3:4]), 
         'Mbp (chr19)', cex=1.3, xpd=TRUE, adj=0)
    pushViewport(viewport(layout=grid.layout(3, 4, widths=widths, heights=c(heights[1], .8 * heights[2], .2 * heights[2]))))
    pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
    biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome='chr19', start=range[1], end=range[2], 
                                        stacking="squish",
                                        # filter=list(biotype='protein_coding'), name="ENSEMBL", biomart=bm)
                                        filter=list(with_refseq_mrna=TRUE, biotype='protein_coding'), name="ENSEMBL", biomart=bm)
    par(xpd=NA)
    plotTracks(biomTrack, # col.line=NULL, # col=NULL, 
               lwd=.5, stackHeight=0.7, panel.only=TRUE)
    popViewport(1)
    pushViewport(viewport(layout.pos.col=3:4, layout.pos.row=3))
    plotTracks(itrack, # col.line=NULL, # col=NULL, 
               chromosome='chr19',
               from=range[1], to=range[2], margin=1,
               lwd=.5, stackHeight=0.7, panel.only=TRUE)
    popViewport(1)
    dev.off()

}

# -----------------------------
# Plot extra views/reorderings:
# -----------------------------
if (plot.extra){
    # Cluster dataset: 
    print("Ordering cells by mahattan distance.")
    dt <- dist(t(mat), 'manhattan')
    ht = hclust(dt)
    colord = ht$order
    breaks = calc.breaks(ht, 20, colord)
    reord = mat[, colord]
    ct.list = as.character(cells)
    labels = meta[as.character(cells[colord]), 'GROUP']
    faclabels = as.matrix(as.numeric(labels))
    colset = as.character(rdcol$COLOR)
    # Label runs (that are not NONE):
    lablist = label.runs(faclabels, labels, rdcol)

    # Plot matrix:
    print("Plotting reordered matrix + colors.")
    png(paste0(imgpref, 'chromHMM_region_', infoline, '_col_reord.png'), res=300, width=8, height=11, units='in')
    layout(matrix(c(1,2),1,2), widths=c(1.75,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(mar=c(3,5,2,.25))
    image(t(faclabels), axes=F, col=colset)
    abline(h=c(0,1),lty=1,lw=0.5)
    abline(v=c(-1,1),lty=1,lw=0.5)
    text(y=lablist[[1]],
         x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
         labels=lablist[[2]], srt=0, adj=1, xpd=TRUE,cex=.8, col=lablist[[3]])
    text(x=0, labels='Type/Tissue\n(As Roadmap)',
         y=par()$usr[1]-0.01 * (par()$usr[2]-par()$usr[1]), 
         srt=0, adj=1, xpd=TRUE,cex=.65, col='black')
    mat.element(reord, df, axlab=FALSE)
    abline(h=breaks,lty=1,lw=0.5)
    dev.off()

    # Sort first by rdmp colors then by clusters
    ord2 = order(labels, decreasing=TRUE)
    colord2 = colord[ord2]
    reord2 = mat[,colord2]
    labels2 = meta[as.character(cells[colord2]), 'GROUP']
    colcut = as.numeric(labels2)
    col.breaks = calc.breaks.acut(colcut)
    faclabels2 = as.matrix(colcut)
    colcut.nam = colcut
    names(colcut.nam) = cells[colord2]
    # Loc labels. Adjust vals if very close to each other.
    ldf = data.frame(label=labels2, quant=1:length(labels2) * 1.0/ length(labels2) )
    labeldf = aggregate(quant ~ label, ldf, mean)

    # Plot matrix:
    print("Plotting BYrdcol matrix + colors.")
    png(paste0(imgpref, 'chromHMM_region_', infoline, '_col_BYrdcol.png'), res=300, width=8, height=11, units='in')
    layout(matrix(c(1,2),1,2), widths=c(1.5,10), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(mar=c(3,4,2,.25))
    image(t(faclabels2), axes=F, col=colset, useRaster=T)
    abline(v=par()$usr[1:2],lty=1,lw=0.5)
    abline(h=par()$usr[3:4],lty=1,lw=0.5)
    text(y=labeldf$quant,
         x=par()$usr[1]-0.1*(par()$usr[2]-par()$usr[1]), 
         labels=labeldf$label, srt=0, adj=1, xpd=TRUE,cex=.7, col=colset)
    text(x=0, labels='Type/Tissue\n(As Roadmap)',
         y=par()$usr[1]-0.01 * (par()$usr[2]-par()$usr[1]), 
         srt=0, adj=1, xpd=TRUE,cex=.65, col='black')
    mat.element(reord2, df, axlab=FALSE)
    dev.off()
}

print("Plotting finished.")
