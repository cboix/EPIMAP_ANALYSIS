#!/usr/bin/R
# Pick a diverse subset from the imputation data
# TODO diff domains
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
library(jsonlite)
library(httr)
library(infotheo)
library(ape)
library(dendextend)
prefix = paste0(img, 'Imputation/DataSummary_')

# Functions:
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
                       add.cutoffs=TRUE, horiz=TRUE, legend='', cex.lab=0.35, highlightct=c()){
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
    # grid(nx=nrow(mat), ny=ncol(mat),col='grey',lty='solid',lwd=.25)
    if (horiz){
        if (cex.lab > 0.4){
            text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=colnames(mat), srt=0, adj=1, xpd=TRUE,cex=.75)
        } else {
            text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[3]+0.01*(par()$usr[4]-par()$usr[3]), labels=colnames(mat), srt=0, adj=1, xpd=TRUE,cex=.75)
        }
    } else { 
        text(x=seq(0,1,length.out=nrow(mat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(mat), srt=90, adj=1, xpd=TRUE,cex=.75)
    }
    if (with.rownames * horiz){
        xt = seq(0,1,length.out=nrow(mat))
        yt = par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3])
        if (length(highlightct) > 0){
            ids = which(rownames(mat) %in% highlightct)
            text(x=xt[ids], y=yt, labels=rownames(mat)[ids], col='red', srt=90, adj=1, xpd=TRUE,cex=cex.lab)
            text(x=xt[-ids], y=yt, labels=rownames(mat)[-ids], srt=90, adj=1, xpd=TRUE,cex=cex.lab)
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

# =================
# Load annotations:
# =================
df = read.delim(paste0(prefix, 'ordered_table.tsv'), header=F, stringsAsFactors=F)
names(df) <- c('biosample_term_name','accession','assay_title','assay_target','award', 'bw_file','status','biosample_expt_count','priority','prio')
core_assays <- c('ChIP-seq','ATAC-seq', 'DNase-seq', 'polyA RNA-seq', 'total RNA-seq')
df <- df[df$assay_title %in% core_assays,]
map <- read.delim('Annotation/all_submitted_released_biosample_mapping.tsv',header=F)
rownames(map) <- map[,2]
newbio <- unlist(read.delim('Annotation/new_biosample_names.tsv', header=F, stringsAsFactors=F))
newunique <- unique(newbio)

# Keep direct matches between two tissues:
newanno <- read.delim('Annotation/new_samples_mapping.tsv',header=T)
newanno <- newanno[as.character(newanno$Snyder.samples) != '',]
newanno <- newanno[as.character(newanno$ENTEx) != '',]
# mapentex = lapply(newanno$ENTEx, function(x){grep(x, map$V1)})

# Roughly; map back biosamples to BSS IDs
ct = as.character(unique(df$biosample_term_name))
mapback = lapply(ct, function(x){grep(x, map$V1)})
newback = lapply(newunique, function(x){grep(x, map$V1)})
newct = lapply(newunique, function(x){grep(x, ct)})
bss = sort(unique(unlist(mapback)))

# Matrix to merge data:
m2 = sapply(mapback, function(x){
                y = rep(0, max(bss))
                if (length(x) > 0){ y[x] = 1 }
                return(y) })
tform = matrix(unlist(m2), nrow=length(ct), ncol=max(bss), byrow=T)
# Normalize:
tform = sweep(tform, 1, apply(tform, 1, sum), '/')
tform[is.na(tform)] <- 0

# ==========================
# Get all distance matrices:
# ==========================
fnames <- list.files(path='ChromImpute',pattern='distance_all_*')
marks <- c(sub('.tsv','',sub('distance_all_','',fnames)),'Full')
NMARKS <- length(marks) - 1
ll <- list()

ct.list = as.character(map$V2)
N = length(ct.list)
full <- matrix(0.0, nrow=N, ncol=N, dimnames=list(ct.list, ct.list))
full_med <- matrix(0.0, nrow=N, ncol=N, dimnames=list(ct.list, ct.list))
full_norm <- matrix(0.0, nrow=N, ncol=N, dimnames=list(ct.list, ct.list))
used <- matrix(0.0, nrow=N, ncol=N, dimnames=list(ct.list, ct.list))
tmp <- full
for (i in 1:NMARKS){
    mark <- marks[i]
    print(paste0('Evaluating ',mark))
    mat <- as.matrix(read.delim(paste0('ChromImpute/',fnames[i]),sep="\t",header=T))
    rn = colnames(mat)
    rownames(mat) <- rn
    idx <- mat == 0
    mat[idx] <- t(mat)[idx]
    # Fill:
    md <- median(as.numeric(unlist(mat)),na.rm=T) # TODO median or mean?
    tmp[] <- NA
    tmp[rn,rn] <- mat[rn,rn]
    tmp[is.na(tmp)] <- md
    full_med <- full_med + tmp
    full[rn,rn] <- full[rn, rn] + mat[rn, rn]
    used[rn, rn] = used[rn, rn] + 1
    full_norm[rn,rn] <- full_norm[rn, rn] + (mat[rn, rn] - md) / diff(range(mat[rn,rn])) # TODO div by sd works 
}
mat <- full / used
mat[is.na(mat)] <- 0
mat <- (mat - diag(diag(mat)))^2
plot.sym(mat, col=spec)

matN <- full_norm / used
matN[is.na(matN)] <- 0
#matN <- (matN - diag(diag(matN)))^2
plot.sym(matN, col=spec)

mat2 = full_med / NMARKS
mat2 <- (mat2 - diag(diag(mat2)))^2
plot.sym(mat2, col=spec)

tmat = tform %*% mat %*% t(tform)
tmat2 = tform %*% mat2 %*% t(tform)
tmatN = tform %*% matN %*% t(tform)

# Save figure
mat = tmatN
kept = which(rowSums(mat) != 0)
matK= mat[kept, kept]
kept_ct = ct[kept] 
rownames(matK) = kept_ct
colnames(matK) = kept_ct
mat = matK
mat[mat == 0] <- median(mat[mat > 0])
dt = dist(1 - mat, 'cosine') 
ht <- hclust(dt, 'ward.D2')
cocl <- order.optimal(dt, ht$merge)$order 
NBREAKS=30
breaks <- get.breaks(ht, cocl, NBREAKS)
FULLMAT = mat # Save

pdf(paste0(prefix, 'distances_', NBREAKS, '.pdf'), width=18, height=15)
plot.sym(mat[cocl,cocl], breaks=breaks, colramp=spec, yaxlab=colnames(mat)[cocl], quant=0.02)
dev.off()

write.table(mat[cocl,cocl], paste0(prefix, 'ordered_distances.tsv'), row.names=T, col.names=T, quote=F, sep="\t")

# Use same breaks for coloring tree:
pdf(paste0(prefix, 'tree.pdf'), width=18, height=8)
par(mar=c(2,2,2,0))
plot(ht, cex=.4, lwd='.5',lty='dashed', yaxt='n')
dev.off()

# using piping to get the dend
dend = as.dendrogram(ht)

# Highlight priority targets:
prioct = ct[unique(unlist(newct))]
idx = which(labels(dend) %in% prioct)
labcol = rep("black", length(labels(dend)))
labcol[idx] = "red"
dend = set(dend, "labels_col", labcol)
dend = set(dend, "labels_cex", .3)

pdf(paste0(prefix, 'cut_tree_', NBREAKS, '.pdf'), width=8, height=18)
par(mar=c(2,2,2,10))
plot(color_branches(dend, k=NBREAKS), horiz=T, cex=.3, yaxt='n')
rect.dendrogram(dend, k=NBREAKS, horiz=T)
abline(v = heights_per_k.dendrogram(dend)[str(NBREAKS)] + .6, lwd = .5, lty = 2, col = "black")
dev.off()

# ===============
# Make choices:
# ===============
seed = c('IMR-90', 'H1-hESC','H9', 'trophoblast cell', 'mesenchymal stem cell', 'neural stem progenitor cell', 'K562', 'adrenal gland', 'heart left ventricle')
newsamp = lapply(newanno$Snyder, function(x){which(ct %in% x)})
matches = unlist(newsamp)
newct = ct[matches]
core = c(seed, newct)
indexes = sapply(core, function(x){which(ct %in% x)})

# =========================
# Plot data reduced to core
# =========================
coredf = df[df$biosample_term_name %in% core,]
tt <- to.table(coredf)
core.mat <- tt[[1]]
desc = 'Core + New Samples'
plot.avail(core.mat> 0, col=rdgy[95], with.rownames=TRUE, title=paste(desc, 'Datasets - Availability'), cex.lab=.7, highlightct=newct)

bw_file.mat = to.table(coredf, epord=tt[[2]], ctord=tt[[3]], attribute='bw_file')[[1]]
bw_file.mat[bw_file.mat == 'signal'] <- 1
bw_file.mat[bw_file.mat == 'other'] <- 2
bw_file.mat[bw_file.mat == 'none'] <- 3
class(bw_file.mat) <- 'numeric'
plot.avail(bw_file.mat, col=c('darkblue','firebrick', 'orange'), with.rownames=TRUE, add.cutoffs=TRUE, title=paste(desc, 'Datasets - Processed BigWig'),
           legend=c("Signal p-value", "Other (such as read-depth normalized signal)", "No BigWig"), cex.lab=.7, highlightct=newct)

# Save figure
ordmat = FULLMAT[cocl,cocl]
subidx = which(colnames(ordmat) %in% core)
submat = ordmat[subidx,subidx]

# FUNCTION TO PLOT DATA:
plot.set <- function(set, tag, desc){
    # Plot data: core + extra set
    coredf = df[df$biosample_term_name %in% set,]
    tt <- to.table(coredf)
    core.mat <- tt[[1]]
    pdf(paste0(prefix, tag, '_subset_availability.pdf'), width=15, height=8)
    plot.avail(core.mat> 0, col=rdgy[95], with.rownames=TRUE, title=paste(desc, 'Datasets - Availability'), cex.lab=.6, highlightct=newct)
    dev.off()
    # BW availability
    pdf(paste0(prefix, tag, '_subset_hasBW.pdf'), width=15, height=8)
    bw_file.mat = to.table(coredf, epord=tt[[2]], ctord=tt[[3]], attribute='bw_file')[[1]]
    bw_file.mat[bw_file.mat == 'signal'] <- 1
    bw_file.mat[bw_file.mat == 'other'] <- 2
    bw_file.mat[bw_file.mat == 'none'] <- 3
    class(bw_file.mat) <- 'numeric'
    plot.avail(bw_file.mat, col=c('darkblue','firebrick', 'orange'), with.rownames=TRUE, add.cutoffs=TRUE, title=paste(desc, 'Datasets - Processed BigWig'),
               legend=c("Signal p-value", "Other (such as read-depth normalized signal)", "No BigWig"), cex.lab=.6, highlightct=newct)
    dev.off()
    # Distances context:
    pdf(paste0(prefix, tag, '_subset_distances.pdf'), width=18, height=15)
    plot.sym(ordmat, colramp=spec, breaks=breaks, yaxlab=colnames(ordmat), quant=0.02, highlight=set)
    title(main=paste(desc, 'datasets within distance matrix'))
    dev.off()
    # Write set:
    write.table(coredf, paste0(prefix, tag, '_subset_experiments.tsv'), row.names=F, col.names=T, quote=F, sep="\t")
}

for (NMORE in c(40, 50, 60)){
    # Choose N more according to the MAXMIN criterion:
    # Pick the one whose min dist with the data is maximal.
    subidx = which(colnames(ordmat) %in% core)
    rn = 1:nrow(ordmat)
    extra = c()
    while(length(subidx) < NMORE){
        top = which.max(apply(ordmat[subidx, -subidx], 2, min))
        extra = c(extra, names(top))
        subidx = c(subidx, rn[-subidx][top])
    }
    coreset = c(core, extra)
    tag = paste0('core_', NMORE)
    plot.set(coreset, tag, 'Core/New + Diverse')

    # Choose most diverse subset by measure, starting from H1:
    subidx = which(colnames(ordmat) == 'H1-hESC')
    rn = 1:nrow(ordmat)
    extra = c()
    for (i in 2:NMORE){
        if (length(subidx) > 1){
            top = which.max(apply(ordmat[subidx, -subidx], 2, min))
        } else { 
            top = which.max(ordmat[subidx, -subidx])
        }
        extra = c(extra, names(top))
        subidx = c(subidx, rn[-subidx][top])
    }
    diverse = c('H1-hESC', extra)
    tag = paste0('diverse_', NMORE)
    plot.set(diverse, tag, 'Most Diverse (MaxMin criterion)')
}

