#!/usr/bin/R
library(infotheo)
library(ape)
library(ggplot2)
source('~/data/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')

# Functions:
perc.rank <- function(x) trunc(rank(x))/length(x)

plot.imputed = FALSE
plot.imputed = TRUE
if (plot.imputed){
    fnames <- list.files(path='ChromImpute',pattern='imp_distance_all_*')
    marks <- c(sub('.tsv','',sub('imp_distance_all_','',fnames)),'Full')
    imgdir = paste0(img, "imp_distance/")
} else {
    fnames <- list.files(path='ChromImpute',pattern='distance_all_*')
    marks <- c(sub('.tsv','',sub('distance_all_','',fnames)),'Full')
    imgdir = paste0(img, "distance/")
}

cmd = paste0('mkdir -p ', imgdir)
system(cmd)
imgpref = paste0(imgdir, 'distances_')

map <- read.delim('Annotation/all_submitted_released_biosample_mapping_20180924.tsv',header=F)
rownames(map) <- map[,2]

# PARAMETERS:
NMARKS <- length(marks) - 1
NCLUST=10

# ==========================
# Get all distance matrices:
# ==========================
ll <- list()
ct.list <- c()
for (i in 1:NMARKS){
    mat <- read.delim(paste0('ChromImpute/',fnames[i]),sep="\t",header=T)
    mark <- marks[i]
    if (dim(mat)[1] > 2){
        idx = grep("^BSS", colnames(mat))
        mat <- mat[idx, idx]
    }
    if (dim(mat)[1] > 2){
        rownames(mat) <- colnames(mat)
        print(paste0('Evaluating ',mark))
        # Fill in backwards:
        idx <- mat == 0
        mat[idx] <- t(mat)[idx]
        mat[is.na(mat)] <- 0
        mat[mat == 0] <- NA
        # Fill in missing with median guess: 
        medians1 = apply(mat, 1, mean, na.rm=T)
        medians2 = apply(mat, 2, mean, na.rm=T)
        repl = outer(medians1, medians2, '+') / 2
        mat[is.na(mat)] =  repl[is.na(mat)]
        diag(mat) = NA
        # Replace names:
        colnames(mat) <- map[colnames(mat),1]
        rownames(mat) <- colnames(mat)
        ct.list <- sort(unique(c(ct.list,colnames(mat))))
        mat <- as.matrix(mat)
        mat.orig <- mat
        mat <- reord(mat)
        mat <- mat[,rownames(mat)]
    }
    ll[[mark]] <- mat
}

# ============================
# Rough merge of all datasets:
# ============================
N <- length(ct.list)
full  <- matrix(0, nrow=N, ncol=N, dimnames=list(ct.list,ct.list))
tmp <- full
ll2 <- list()
for (i in 1:NMARKS){
    mat2 <- matrix(NA, nrow=N, ncol=N, dimnames=list(ct.list,ct.list))
    mark <- marks[i]
    mat <- ll[[mark]]
    if (dim(mat)[1] > 2){
        nn<- rownames(mat)
        mat[is.na(mat)] <- 0
        md <- median(as.numeric(unlist(mat)),na.rm=T)
        tmp[] <- NA
        tmp[nn,nn] <- mat[nn,nn]
        tmp[is.na(tmp)] <- md
        full <- full + tmp
        mat2[nn,nn] <- mat[nn,nn]
    }
    ll2[[mark]] <- mat2  # Re-sized Matrices for comparison
}

mat <- full / NMARKS
mat <- (mat - diag(diag(mat)))^2
mat <- reord(mat)
mat <- mat[,rownames(mat)]
ORD = rownames(mat)
ll[['Full']] <- mat

# ----------------------------------------
# Get availability from sample/mark table:
# ----------------------------------------
mat = read.delim('ChromImpute/sample_mark_table.tsv', header=F)
names(mat) <- c('CellType','Epitope','file')
mat$file=1

# Reordering by availability: 
cells <- aggregate(file ~ CellType, mat, length)
epitopes <- aggregate(file ~ Epitope, mat, length)
ordCT <- as.character(with(cells,CellType[order(file, decreasing=TRUE)]))
ordE <- as.character(with(epitopes,Epitope[order(file, decreasing=TRUE)]))

wm <- spread(mat, CellType, file, fill=0)
rownames(wm) <- wm$Epitope 
wm <- wm[ordE,ordCT]

# Get mark availability
matcol=heat.colors(12)
matcol=c('white','indianred')
colnames(wm) <- map[colnames(wm),1]

# Main marks: 
main.marks = marks[unlist(lapply(ll, nrow)) > 20]
main.marks = main.marks[-length(main.marks)]  # Remove full.
NMAIN = length(main.marks)
crossnames = outer(ct.list, ct.list, FUN=paste, sep="\n")
comparison = as.character(crossnames)

# Percentiles across all
pc = list()
dc = list()
for (i in 1:NMAIN){
    m1 = main.marks[i]
    x1 = ll2[[m1]]
    x1[x1 == 0] = NA
    perc = x1
    idx = !is.na(x1)
    perc[idx] = perc.rank(x1[idx])
    near = apply(perc, 1, min, na.rm=T)
    near[near == Inf] = NA
    near2 = apply(x1, 1, min, na.rm=T)
    near2[near2 == Inf] = NA
    pc[[m1]] = near
    dc[[m1]] = near2
}

pcdf = as.data.frame(pc)
# Ordered by min:
mins = apply(pcdf, 1, min , na.rm=T)
pcdf = pcdf[order(mins, decreasing=F),]
pcmat = t(pcdf)

png(paste0(imgpref,'percentiles_all.png'), res=450, units='in', width=6, height=11)
par(mar=c(2,10,2,2))
par(yaxs='i', xaxs='i', cex.lab=1.2, cex.axis=.30) 
boxplot(pcmat, horizontal=TRUE, main='All Marks - Nearest Percentiles', col='darkgrey', ylim=c(0,1),
        names=colnames(pcmat), las=2, cex=.5)
abline(h=seq(1,ncol(pcmat),5),lty=2,lwd=.5)
dev.off()

dcdf = as.data.frame(dc)
mins = apply(dcdf, 1, min , na.rm=T)
dcdf = dcdf[order(mins, decreasing=F),]
dcmat = t(dcdf)
ranked = 1:length(names(mins))
names(ranked) = names(sort(mins))
top10 =  rev(tail(names(sort(mins)), 10))

png(paste0(imgpref,'nearest_all.png'), res=450, units='in', width=6, height=11)
par(mar=c(2,10,2,2))
par(yaxs='i', xaxs='i', cex.lab=1.2, cex.axis=.30) 
boxplot(dcmat, horizontal=TRUE, main='All Marks - Nearest Distances', col='darkgrey', ylim=c(0,1),
        names=colnames(dcmat), las=2, cex=.5)
abline(h=seq(1,ncol(dcmat),5),lty=2,lwd=.5)
dev.off()


# Plotting:
for (i in 1:(NMARKS +1)){
    cx <- c(.45,rep(.5, NMARKS - 1),.25)[i]
    mark <- marks[i]
    print(paste0('Plotting ',mark))
    mat <- ll[[mark]]
    tag = paste0(mark,' (', dim(mat)[1], ')')
    if (dim(mat)[1] > 10){
        mat[is.na(mat)] <-0 
        mat.save <- mat
        mat[mat == 0] <- NA
        nn <- rownames(mat)
        ord <- ORD[ORD %in% nn]
        mat <- mat[ord,ord]

        # TODO distance measure needs to be clear.
        # ==============
        # 1. PLOT MATRIX
        # ==============
        png(paste0(imgpref,'distance_',mark,'.png'),res=450,units='in',width=12.5,height=11)
        par(mar=c(2,14,2,2))
        image(mat, axes=FALSE,main=paste(tag,"Dataset Distance"),col=colryb)
        grid(nx=nrow(mat), ny=ncol(mat),col='grey',lty='solid',lwd=.25)
        text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), labels=colnames(mat), srt=0, adj=1, xpd=TRUE,cex=cx)
        dev.off()

        # 1b. DETERMINE CUTOFFS:
        # 2. CLUSTER AND PLOT
        # ========================================
        # 3. PLOT + LABEL DENDROGRAM FROM DISTANCE
        # ========================================
        tr <- nj(mat.save)
        png(paste0(imgpref,'dendro_',mark,'.png'),res=450,units='in',width=7,height=11)
        par(mar=c(.5,.5,4,0))
        plot(tr,'phylogram', main=paste(tag,"Dataset Distance"),
             cex=cx+.2, edge.color = "darkblue",
             edge.width = 1, edge.lty = 2, tip.color = "darkblue")
        dev.off()

        # -----------------------------------
        # 3c. Alternatively plot hclust tree:
        # -----------------------------------
        dt <- dist(mat,'eJaccard')
        ht <- hclust(dt)
        cocl <- order.optimal(dt, ht$merge)$order
        ht$order <- cocl
        pt <- as.phylo(ht)

        png(paste0(imgpref,'hclusttree_',mark,'.png'),res=450,units='in',width=5,height=11)
        par(mar = c(0,0,2,0))
        plot(pt, type = "cladogram", cex=cx+.2,
             main=paste(tag,"Dataset Distance"),
             edge.color = "darkblue", edge.width = 2, edge.lty = 2,
             tip.color = "darkblue")
        dev.off()

        # CUTOFF MATRIX:
        if (dim(mat)[1] > 100) {
            NCLUST = 20 
        } else if (dim(mat)[1] < 20){
            NCLUST=4 
        } else { NCLUST = 10}
        acut <- cutree(ht,NCLUST)[cocl]
        cuts <- cumsum(rle(acut)$lengths)
        step <- head(diff(seq(0,1,length.out=nrow(mat))),1)
        cuts <-  cuts[-length(cuts)] -1
        breaks <- cuts * step + step/2
        mat.reclust <- mat[cocl,cocl]

        # Make availability:
        avail = matrix(NA, nrow=dim(mat)[1], ncol= NMARKS)
        rownames(avail) = rownames(mat.reclust)  # Ordered as cluster)
        colnames(avail) = marks[-length(marks)]

        rn = rownames(avail)[rownames(avail) %in% colnames(wm)]
        avail[rn, colnames(avail)] = t(wm[colnames(avail), rn])
        marg = apply(avail,2, sum, na.rm=T)
        avail = avail[, marg > 5]
        avail = t(avail)
        # If in mark, should be marked as available.
        if (mark != 'Full'){ avail[mark, is.na(avail[mark,])] = 0}

        # --------------------------------
        # 3b. PLOT MATRIX WITH DENDROGRAM:
        # --------------------------------
        png(paste0(imgpref,'matrixtree_',NCLUST,'_',mark,'.png'),res=300,units='in',width=12,height=11)
        layout(matrix(c(1,2,3,4),2,2), widths=c(10,1), heights=c(1,8), TRUE)
        par(yaxs="i")
        par(xaxs="i")
        par(mar=c(0,14,2,.25))
        plot(ht,main=paste(tag,"Dataset Similarity"), col='black',lwd=1,lty=1,hang=-1,labels=FALSE,axes=F,ylab='')
        par(mar=c(3,14,0,.25))
        image(mat.reclust, axes=FALSE,col=colryb)
        grid(nx=nrow(mat.reclust), ny=ncol(mat.reclust),col='grey',lty='solid',lwd=.25)
        text(y=seq(0,1,length.out=ncol(mat.reclust)), x=par()$usr[3]-0.001*(par()$usr[4]-par()$usr[3]), labels=colnames(mat.reclust), srt=0, adj=1, xpd=TRUE,cex=cx)
        abline(v=breaks,lty=2,lwd=.5)
        abline(h=breaks,lty=2,lwd=.5)
        # Column 2
        par(mar=c(0,0,2,2))
        plot.new()
        par(mar=c(3,0,0,2))
        image(avail, axes=FALSE,col=matcol) # , zlim=c(0, max(avail)))
        grid(nx=nrow(avail), ny=ncol(avail),col='grey',lty='solid',lwd=.25)
        abline(h=breaks,lty=2,lwd=.5)
        text(x=seq(0,1,length.out=nrow(avail)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(avail), srt=90, adj=1, xpd=TRUE,cex=.5)
        dev.off()

        # --------------------------------
        # 3b. PLOT MATRIX WITH DENDROGRAM AND RANKING.
        # --------------------------------
        png(paste0(imgpref,'matrixtree_v3_',NCLUST,'_',mark,'.png'),res=300,units='in',width=12,height=11)
        layout(matrix(c(1,2,3,4),1,4), widths=c(1.5,10,3,.25), TRUE)
        par(yaxs="i")
        par(xaxs="i")
        par(mar=c(4,1,14,.25))
        plot(pt, type = "phylogram", cex=.0001,
             edge.color = "black", edge.width = 1, edge.lty = 1)
        par(mar=c(4,0,14,0))
        image(mat.reclust, axes=FALSE,col=colryb)
        grid(nx=nrow(mat.reclust), ny=ncol(mat.reclust),col='grey',lty='solid',lwd=.25)
        abline(v=breaks,lty=2,lwd=.5)
        abline(h=breaks,lty=2,lwd=.5)
        title(paste(tag,"Dataset Similarity"), outer = TRUE, line = -2)
        text(x=seq(0,1,length.out=ncol(mat.reclust)), y=par()$usr[4]+0.001*(par()$usr[4]-par()$usr[3]), labels=colnames(mat.reclust), srt=-90, adj=1, xpd=TRUE,cex=cx)
        # Column 3
        par(mar=c(4,10,14,.25)) 
        image(avail, axes=FALSE,col=matcol) # , zlim=c(0, max(avail)))
        grid(nx=nrow(avail), ny=ncol(avail),col='grey',lty='solid',lwd=.25)
        abline(h=breaks,lty=2,lwd=.5)
        text(x=seq(0,1,length.out=nrow(avail)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(avail), srt=90, adj=1, xpd=TRUE,cex=.75)
        rankings = as.matrix(ranked)^2
        rankings = as.matrix(rankings[colnames(mat.reclust),])
        cxlist = 1 * rankings / max(rankings)
        cxlist[cxlist < 0.6] = 0.01
        cxlist[cxlist > 0.6] = 0.4
        collist = rep('black', length(cxlist))
        collist[cxlist > 0.1] = 'red'
        text(y=seq(0,1,length.out=ncol(mat.reclust)), x=par()$usr[3]-0.07*(par()$usr[4]-par()$usr[3]), labels=colnames(mat.reclust), srt=0, adj=1, xpd=TRUE,cex=cxlist)
        par(mar=c(4,0,14,1))
        image(t(rankings), axes=FALSE,col=colramp) 
        dev.off()

        # --------------------------------
        # RANKING AND AVAILABILITY
        # --------------------------------
        png(paste0(imgpref,'matrixtree_v4_',NCLUST,'_',mark,'.png'),res=300,units='in',width=4,height=13)
        layout(matrix(c(1,2),1,2), widths=c(3,.75), TRUE)
        par(yaxs="i")
        par(xaxs="i")
        if (mark == 'DNase-seq'){ 
            par(mar=c(6,11,1,.25)) 
        } else {
            par(mar=c(6,9,1,.25)) 
        }
        image(avail, axes=FALSE,col=matcol) # , zlim=c(0, max(avail)))
        grid(nx=nrow(avail), ny=ncol(avail),col='grey',lty='solid',lwd=.25)
        abline(h=breaks,lty=2,lwd=.5)
        text(x=seq(0,1,length.out=nrow(avail)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(avail), srt=90, adj=1, xpd=TRUE,cex=.75)
        rankings = as.matrix(ranked)^2
        rankings = as.matrix(rankings[colnames(mat.reclust),])
        cxlist = 1 * rankings / max(rankings)
        cxlist[cxlist < 0.6] = 0.01
        cxlist[cxlist > 0.6] = 0.4
        collist = rep('black', length(cxlist))
        collist[cxlist > 0.1] = 'red'
        text(y=seq(0,1,length.out=ncol(mat.reclust)), x=par()$usr[3]-0.07*(par()$usr[4]-par()$usr[3]), labels=colnames(mat.reclust), srt=0, adj=1, xpd=TRUE,cex=.30, col=collist)
        par(mar=c(6,0,1,1))
        avail.sum = colSums(exp(avail), na.rm=T)
        avail.sum[avail.sum == 0] <- 1
        priomat = cbind(sqrt(rankings), avail.sum, sqrt(rankings) / avail.sum) 
        image(t(priomat), axes=FALSE,col=colramp) 
        text(x=seq(0,1,length.out=ncol(priomat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), 
             labels=c('Ranking', 'Prior Assays', 'Ranking/Assays'), srt=90, adj=1, xpd=TRUE,cex=.75)
        dev.off()

        # For reccomendations:
        top20 =  rev(tail(names(sort(priomat[,1])), 20))
        top20d =  rev(tail(names(sort(priomat[,3])), 20))

        # ============
        # 5. PLOT MDS:
        # ============
        c <- cmdscale(dist(t(mat)))
        c <- data.frame(c,rep=rownames(c),clust=factor(acut[rownames(c)]))
        gp <- ggplot(c,aes(X1,X2,label=rep,color=clust)) + geom_point(alpha=.25,cex=2) +
            geom_text(vjust='inward',hjust='inward',nudge_y=0.0001,check_overlap=T,cex=2) +
            labs(x="", y="", title=paste0("MDS of ", mark, " Dataset Distance")) +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5),legend.position='none')
        ggsave(paste0(imgpref,'MDS_',NCLUST,'_',mark,'.png'), gp,dpi=300,units='in',width=12,height=11)

        # ====================
        # 6. Plot distribution
        # ====================
        png(paste0(imgpref,'distr_',mark,'.png'),res=450,units='in',width=7,height=5)
        par(mar = c(2,2,2,0))
        plot(hist(mat), main = paste0(tag,": Distances between Datasets"), xlim=c(0,1), col='darkgrey', border='white', ylab = '', xlab='')
        dev.off()

    }
}

# TODO label as ROADMAP VS NOT ROADMAP.

# Comparison of distances on multiple axes:

# Calc rsq
rsqmat = matrix(NA, NMAIN, NMAIN)
for (i in 1:(NMAIN - 1)){
    for (j in (i+1):NMAIN) {
        m1 = main.marks[i]
        m2 = main.marks[j]
        x1 = as.numeric(ll2[[m1]])
        x2 = as.numeric(ll2[[m2]])
        x1[x1 == 0] <- NA
        x2[x2 == 0] <- NA
        df <- data.frame(x1, x2, comparison = comparison)
        fit = lm(x1 ~ x2, df)
        rsqmat[i,j] = summary(fit)$adj.r.squared
        rsqmat[j,i] = rsqmat[i,j] 
    }
}

# Calculate mutual information: - unclear how much this helps.
MImat = matrix(NA, NMAIN, NMAIN)
for (i in 1:(NMAIN - 1)){
    for (j in (i+1):NMAIN) {
        m1 = main.marks[i]
        m2 = main.marks[j]
        x1 = as.numeric(ll2[[m1]])
        x2 = as.numeric(ll2[[m2]])
        x1[x1 == 0] <- NA
        x2[x2 == 0] <- NA
        df <- data.frame(x1, x2, comparison = comparison)
        df$x12 <- with(df,x1 + x2)
        df <- df[!is.na(df$x12),]
        d.df <- discretize(df[,c('x1','x2')], 'equalwidth', nbins=100)
        mi = mutinformation(d.df)
        e1 = mi[1,1]
        e2 = mi[2,2]
        ce = condentropy(d.df$x1, d.df$x2)
        MImat[i,j] = mi[1,2] / ((e1 + e2) / 2)
        MImat[j,i] = MImat[i,j] 
    }
}

# Reorder by mean
rsqmeans = apply(rsqmat,1, mean, na.rm=T)
rsq.reordered = order(rsqmeans, decreasing=T)

# TODO MDS of cell types across marks simultaneously
# Plot all against all
png(paste0(imgpref,'comparison_all.png'),res=450,units='in',width=25,height=16)
layout(matrix(1:((NMAIN + 1)^2), NMAIN + 1, NMAIN +1), TRUE)
par(yaxs="i")
par(xaxs="i")
breaks = seq(0,1,.05)
for (i in c(rsq.reordered, 0)){
    for (j in c(0, rsq.reordered)) {
        par(mar=c(0,0,0,0))
        if (i == 0 || j == 0){
            if (i == 0){
                if (j == 0){
                    plot.new()
                } else{ 
                    m2 = main.marks[j]
                    m = ll2[[m2]]
                    m[m == 0] <- NA
                    m = m[!is.na(m)] / max(m,na.rm=T)
                    xhist = hist(m, plot=FALSE,breaks=breaks)
                    barplot(xhist$density, axes=FALSE, horiz=TRUE, col='darkgrey', border='white')
                }
            } else if (j == 0){
                m1 = main.marks[i]
                m = ll2[[m1]]
                m[m == 0] <- NA
                m = m[!is.na(m)] / max(m,na.rm=T)
                yhist <- hist(m, plot=FALSE, breaks=breaks) 
                barplot(yhist$density, axes=FALSE, horiz=FALSE, col='darkgrey', border='white')
            }
        } else {
            # Plot normal:
            m1 = main.marks[i]
            m2 = main.marks[j]
            if (i == j){
                plot.new()
                rsq.mean = mean(rsqmat[i,], na.rm=T)
                tag = paste0(m1, '\n', round(rsq.mean,3))
                text(.5,.5, tag, cex=3.5)
            } else {
                x1 = as.numeric(ll2[[m1]])
                x2 = as.numeric(ll2[[m2]])
                x1[x1 == 0] <- NA
                x2[x2 == 0] <- NA
                plot(x1, x2, ylab='', xlab='', pch=19, col=alpha('black', 0.05),
                     xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n')
                if (rsqmat[i, j] < 0.1){
                    text(.5,.5, round(rsqmat[i,j], 3) , cex=3.5, col=alpha('red', 0.75))
                } else {
                    text(.5,.5, round(rsqmat[i,j], 3) , cex=3.5, col=alpha('darkgreen', 0.75))
                }
                abline(0,1, 'grey')
            }
        }
    }
}
dev.off()

# Plot individual comparisons: 
for (i in 1:(NMAIN -1)){
    for (j in (i+1):NMAIN) {
        m1 = main.marks[i]
        m2 = main.marks[j]

        x1 = as.numeric(ll2[[m1]])
        x2 = as.numeric(ll2[[m2]])
        x1[x1 == 0] <- NA
        x2[x2 == 0] <- NA
        df <- data.frame(x1, x2, comparison = comparison)
        if (rsqmat[i, j] < 0.1){ lblcolor='red' } else { lblcolor='darkgreen' }

        gp <- ggplot(df, aes(x1, x2, label=comparison)) + geom_point(alpha=.2,cex=2) +
            geom_abline(intercept=0, slope=1, color='grey') + 
            geom_smooth(method = "lm", se = FALSE) +
            geom_text(vjust='outward',hjust='outward',nudge_y=0.0001,check_overlap=T,cex=3) +
            labs(x=paste(m1, "Distance"), y=paste(m2,"Distance"), title=paste0(m1, ' vs. ', m2)) +
            geom_text(x=.5, y=0.1, label=round(rsqmat[i,j],3), cex=20, color=alpha(lblcolor, .5)) +
            lims(y=c(0,1), x=c(0,1)) + theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5), legend.position='none')
        ggsave(plot=gp, filename=paste0(imgpref,'comparison_',m1,'_',m2,'.png'), dpi=300,
               units='in', width=10,height=10)
    }
}



# TODO plot for each of the marks.
