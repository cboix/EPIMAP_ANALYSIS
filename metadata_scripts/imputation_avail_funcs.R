#!/usr/bin/R
# Functions for plotting imputation results:
library(RColorBrewer)
library(scales)
library(jsonlite)
library(reshape2)
library(httr)
library(cba)
library(ape)
library(infotheo)
library(dendextend)

# Colors:
spec <- colorRampPalette(brewer.pal(n=9,name="Spectral"))(100)
rdgy <- colorRampPalette(brewer.pal(n=11,name="RdGy"))(100)

# df to table:
to.table <- function(df, epord='', ctord='', attribute=''){ 
    formula='accession ~ assay_target + biosample_term_name'
    if (attribute != ''){
        formula = paste(formula, '+', attribute)
    }
    udf <- aggregate(as.formula(formula), length, data=df)
    # Reorder if necessary:
    if (epord == ''){ 
        epdf <- aggregate(accession ~ assay_target, length, data=udf)
        epord <- with(epdf, assay_target[order(accession, decreasing=T)]) 
    }
    if (ctord == ''){ 
        ctdf <- aggregate(accession ~ biosample_term_name, length, data=udf)
        ctord <- with(ctdf, biosample_term_name[order(accession, decreasing=T)]) 
    }
    udf$biosample_term_name <- factor(as.character(udf$biosample_term_name), levels=ctord) 
    udf$assay_target <- factor(as.character(udf$assay_target), levels=epord) 
    udf <- udf[order(udf$biosample_term_name, udf$assay_target),]
    # Reshape to wide:
    if (attribute %in% c('award', 'status', 'bw_file', 'prio')){
        wide <- dcast(udf, assay_target ~ biosample_term_name, value.var=attribute, 
                      fun.aggregate=function(x){ sort(x, decreasing=T)[1]})  # Prioritizes any roadmap expt
    } else if (attribute == '') {
        wide <- dcast(udf, assay_target ~ biosample_term_name, value.var='accession')
    } else {
        wide <- dcast(udf, assay_target ~ biosample_term_name, value.var=attribute)
    }
    mat <- as.matrix(wide[,-1])
    rownames(mat) = wide[,1]
    return(list(mat, epord, ctord))
}

# Plot wide table from above:
plot.avail <- function(mat, colramp = heat.colors(12), title='', with.rownames=FALSE,
                       add.cutoffs=TRUE, horiz=TRUE, legend='', cex.lab='auto', highlightct=c()){
    if (horiz){ 
        mat = t(mat) 
        par(mar=c(4 + with.rownames * 7, 7, 3, 2)) 
    } else {
        par(mar=c(7, 4 + with.rownames * 7, 3, 2)) 
    }
    image(mat, axes=FALSE, main='', col=colramp)
    if (title != '') { mtext(title, side=3, cex=2.25, line=.75)} 
    mtext(paste0('Samples (', ncol(mat) * (!horiz) + nrow(mat) * horiz, ')'), side=2 - horiz * 1, cex=2, line=1.5 + with.rownames * 8)
    mtext(paste0('Genomic Assays (', nrow(mat) * (!horiz) + ncol(mat) * horiz, ')'), side=1 + horiz * 1, cex=2, line=5)
    if (cex.lab == 'auto'){
        NS = nrow(mat)
        if (NS > 80){ cex.lab = 0.35 } else { cex.lab = 1 - 0.65 * NS / 80.0 }
    }
    if (horiz){
        if (NS > 80){ coeff = 0.01 } else { coeff = -0.05 + 0.06 * NS / 80.0 } 
        text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[3] + coeff*(par()$usr[4]-par()$usr[3]), labels=colnames(mat), srt=0, adj=1, xpd=TRUE,cex=.75)
    } else { 
        text(x=seq(0,1,length.out=nrow(mat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(mat), srt=90, adj=1, xpd=TRUE,cex=.75)
    }
    if (with.rownames * horiz){
        xt = seq(0,1,length.out=nrow(mat))
        yt = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3])
        ids = which(rownames(mat) %in% highlightct)
        if (length(ids) > 0){
            if (length(rownames(mat)[ids]) > 0){
                text(x=xt[ids], y=yt, labels=rownames(mat)[ids], col='red', srt=90, adj=1, xpd=TRUE,cex=cex.lab) }
            if (length(rownames(mat)[-ids]) > 0){
                text(x=xt[-ids], y=yt, labels=rownames(mat)[-ids], srt=90, adj=1, xpd=TRUE,cex=cex.lab) }
        } else { text(x=xt, y=yt, labels=rownames(mat), srt=90, adj=1, xpd=TRUE,cex=cex.lab) }
    } else if (with.rownames){
        text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]), labels=colnames(mat), srt=0, adj=1, xpd=TRUE,cex=cex.lab)
    }
    if (add.cutoffs){
        if (horiz){ mat <- t(mat) } # revert }
        marks <- rownames(mat)
        N.EPI <- length(marks)
        pct <- rep(0, N.EPI)
        num <- rep(0, N.EPI)
        cuts <- 1:N.EPI
        for (cutoff in cuts){
            mat2 <- mat[marks[1:cutoff],]
            num[cutoff] = sum(!(is.na(mat2)))
            if (cutoff == 1){
                pct[cutoff] = num[cutoff] / length(mat2)
            } else {
                pct[cutoff] = num[cutoff] / prod(dim(mat2))
            }
        }
        breaks = c(7,12,16,19,nrow(mat))
        bcols <- (breaks - .5)/ (N.EPI-1)
        bpct = pct[breaks]
        bnum = num[breaks]
        ll <- diff(c(0,breaks))
        mids <- ll / 2 + c(0, breaks)[-(length(breaks) +1)]
        textpos <- (mids -.5) / (N.EPI- 1)
        ypos = 0.55 + horiz * 0.2
        projtext = paste0(breaks, ' assays (', round(bpct*100, 1), '%: ', bnum, ')')
        if (horiz) {
            abline(h=bcols,col='darkred',lty='dashed',lwd=2)
            text(y=textpos, x=ypos, labels=projtext, srt=0, adj=0, xpd=TRUE, col=alpha('darkred',.85), cex=2)
        } else { 
            abline(v=bcols,col='darkred',lty='dashed',lwd=2)
            text(x=textpos, y=ypos, labels=projtext, srt=90, adj=0, xpd=TRUE, col=alpha('darkred',.85), cex=2)
        }
    }
    if (legend != ''){
        if (horiz) {
            legend(.1, 1.03, legend=legend, col=colramp, lty=1, lwd=2, cex=1.5, box.lty=0)
        } else {
            legend(.1, 1.03, legend=legend, col=colramp, lty=1, lwd=2, cex=1.5, box.lty=0)
        }
    }
}

# Plot symmetric heatmap
plot.sym <- function(mat, labels=c(), colramp=heat.colors(12), breaks=c(),
                     yaxlab=c(), quant=0.05, cex.lab=.3, highlight=c()){ 
    if (length(labels) > 0){
        set <- colorRampPalette(brewer.pal(n=length(unique(labels)),
                                           name="Set1"))(length(unique(labels)))
        layout(matrix(c(1,2),1,2), widths=c(10, .5), TRUE)
    } 
    par(yaxs="i")
    par(xaxs="i")
    if (length(yaxlab) > 0){
        par(mar=c(2,10,2,.25))
    } else { 
        par(mar=c(4,1,2,.25))
    }
    # Cut mat:
    q1 = quantile(mat, quant)
    q2 = quantile(mat, 1 - quant)
    mat[mat < q1] <- q1
    mat[mat > q2] <- q2
    image(mat, axes=F, zlim=c(q1,q2), col=colramp, useRaster=TRUE)  # Raster more efficient?
    if (length(yaxlab) > 0){
        # grid(nx=nrow(mat), ny=ncol(mat),col='grey',lty='solid',lwd=.25)
        yt = seq(0,1,length.out=ncol(mat))
        xt = par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3])
        if (length(highlight) > 0){
            ids = which(yaxlab %in% highlight)
            print(here)
            text(y=yt[-ids], x=xt, col='grey', labels=yaxlab[-ids], srt=0, adj=1, xpd=TRUE, cex=cex.lab)
            text(y=yt[ids],  x=xt, col='red',  labels=yaxlab[ids],  srt=0, adj=1, xpd=TRUE, cex=cex.lab)
        } else { text(y=yt, x=xt, labels=yaxlab, srt=0, adj=1, xpd=TRUE,cex=cex.lab) }
    }
    if (length(breaks) > 0){
        abline(v=breaks,lty=2,lwd=.5)
        abline(h=breaks,lty=2,lwd=.5)
    }
    if (length(labels) > 0){
        par(mar=c(4,0,2,0))
        image(t(labels), axes=F, col=set)
    }
}

get.breaks <- function(ht, order, nclust){
    acut <- cutree(ht, nclust)[order]
    cuts <- cumsum(rle(acut)$lengths)
    step <- head(diff(seq(0,1,length.out=length(order))),1)
    cuts <-  cuts[-length(cuts)] -1
    return(cuts * step + step/2)
}

pickToN <- function(mat, set, N){
    rn = 1:nrow(mat)
    print(set)
    if (is.null(set)){
        idx = c()
        set = c()
    } else {
        idx = which(colnames(mat) %in% set)
        set = colnames(mat)[idx]
    }
    # TODO break if idx is depleted....
    while(length(idx) < N){
        if (length(idx) == 0){
            top = which.max(apply(mat, 2, min))
            set = names(top)
            idx = rn[top]
        } else {
            if (length(idx) == 1){
                top = which.max(mat[idx, -idx])
            } else { 
                top = which.max(apply(mat[idx, -idx], 2, min))
            }
            set = c(set, names(top))
            idx = c(idx, rn[-idx][top])
        }
    }
    print(set)
}

# Function to plot data:
# TODO SPLIT:
plot.set <- function(set, tag, desc){
    # Plot data: core + extra set
    coredf = df[df$biosample_term_name %in% set,]
    tt <- to.table(coredf)
    core.mat <- tt[[1]]

    # PLOT 1:
    plot.avail(core.mat> 0, col=rdgy[95], with.rownames=TRUE, title=paste(desc, 'Datasets - Availability'), cex.lab=.6, highlightct=newct)

    # PLOT 2: BW availability
    bw_file.mat = to.table(coredf, epord=tt[[2]], ctord=tt[[3]], attribute='bw_file')[[1]]
    bw_file.mat[bw_file.mat == 'signal'] <- 1
    bw_file.mat[bw_file.mat == 'other'] <- 2
    bw_file.mat[bw_file.mat == 'none'] <- 3
    class(bw_file.mat) <- 'numeric'
    plot.avail(bw_file.mat, col=c('darkblue','firebrick', 'orange'), with.rownames=TRUE, add.cutoffs=TRUE, title=paste(desc, 'Datasets - Processed BigWig'),
               legend=c("Signal p-value", "Other (such as read-depth normalized signal)", "No BigWig"), cex.lab=.6, highlightct=newct)

    # PLOT 3: Distances context:
    plot.sym(ordmat, colramp=spec, breaks=breaks, yaxlab=colnames(ordmat), quant=0.02, highlight=set)
    title(main=paste(desc, 'datasets within distance matrix'))
}

