#!/usr/bin/R
# Plot the processing status report for ENCODE3/ROADMAP DATA
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')

# All processing files:
files <- list.files(path='Annotation',pattern='processing_status_matrix_.*.tsv')
datstr <- sub(".tsv","",sub("processing_status_matrix_","",files))
dates <- as.POSIXct(datstr,format="%Y%m%d-%H%M")

# ==============================================
# Get most recent processing status information:
# ==============================================
id <- which.max(dates)
date <- dates[id]
dstr <- datstr[id]
day <- as.character(date,'%m/%d/%Y')
mat <- read.delim(paste0('Annotation/',files[id]),header=T,sep=" ")
marks <- c('WCE',scan('Annotation/histone_mark_list','c'),'DNase-seq','CTCF')
# TODO also plot core histone mark list?

# Read map file:
mapping <- read.delim('Annotation/all_submitted_released_biosample_mapping.tsv')
names(mapping) <- c('CellTypeName','CellType')
mat <- merge(mat, mapping)

# Reordering by availability: 
cells <- aggregate(CellRepl ~ CellType,mat,sum)
epitopes <- aggregate(CellRepl ~ Epitope,mat,sum)
ordCT <- as.character(with(cells,CellType[order(CellRepl,decreasing=TRUE)]))
ordE <- as.character(with(epitopes,Epitope[order(CellRepl,decreasing=TRUE)]))

wide <- list()
for (nam in colnames(mat)[-c(1:3)]){ 
    mat[,nam] <- as.numeric(as.character(mat[,nam]))
    tmp <- dcast(mat,Epitope ~ CellType,value.var=nam)
    rownames(tmp) <- tmp$Epitope 
    wide[[nam]] <- tmp[ordE,ordCT]
}

# ===================
# Availability plots:
# ===================
abv <- c('CellRepl','RepsProcessed','PoolSCCA','PeakCalls','SignalTracks','Depth','SubsampleDepth','Fraglen','FraglenSubsample')
descriptions <- c('Replicates','Replicates Processed','Pooled and SCCA complete','Peaks Called','Signal Tracks Generated','Total Read Depth','Subsampled Read Depth','Fragment Length (full data)','Fragment Length (subsampled data)')

# Make processing report directory:
procdir <- paste0('processing_status_',dstr)
dir.create(file.path(img, procdir), showWarnings = FALSE)
for (i in 1:length(abv)){ 
    nam <- abv[i]
    desc <- descriptions[i]
    print(desc)
    wm <- as.matrix(wide[[nam]])
    # Modify depending on plot:
    matcol=c('white','darkblue')
    if (nam %in% c('CellRepl','RepsProcessed')){ 
        wm <- log(wm)
        matcol=heat.colors(12)
    } else if (nam %in% c('Depth','SubsampleDepth')){ 
        wm[wm == 0] <- NA
        wm[wm >= 30] <- 30
        wm[wm < 30] <- 1
        matcol=c('darkblue','red')
    } else if (nam %in% c('Fraglen','FraglenSubsample')) { 
        wm[wm <= 0] <- NA
        rn <- range(wm,na.rm=T)
        desc <- paste0(desc,' - [',rn[1],',',rn[2],']')
        matcol=heat.colors(12)
        # TODO check if any below expected - if so change params to SCCA.
    }

    # Mark off every 50:
    cols <- (seq(0,nrow(wm),50) - .5)/ (nrow(wm) -1)
    rows <- (seq(0,ncol(wm),50) - .5)/ (ncol(wm) -1)

    # ADD TEXT: DIMENSIONS(AXES), SPARSITY, and % COMPLETE.
    png(paste0(img,'/',procdir,'/',i,'-',nam,'.png'),res=450,units='in',width=17,height=11)
    par(mar=c(5,4,3,2))
    image(wm, axes=FALSE,main='',col=matcol)
    mtext(paste0(desc,' - ',day), side=3, cex=2.25, line=.75)
    mtext(paste0('Biosamples (', ncol(wm), ')'), side=2, cex=2, line=1)
    mtext(paste0('ChIP-seq Epitopes and DNase-seq (', nrow(wm), ')'), side=1, cex=2, line=3.5)
    # image(wm, axes=FALSE,main=paste0(desc,' - ',day))
    grid(nx=nrow(wm), ny=ncol(wm),col='grey',lty='solid',lwd=.25)
    text(x=seq(0,1,length.out=nrow(wm)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(wm), srt=90, adj=1, xpd=TRUE,cex=.35)
    # text(y=seq(0,1,length.out=ncol(wm)), x=par()$usr[3]+0.001*(par()$usr[4]-par()$usr[3]), labels=colnames(wm), srt=0, adj=1, xpd=TRUE,cex=.35)
    abline(h=rows,col='black',lty='solid',lwd=.5)
    abline(v=cols,col='black',lty='dashed',lwd=.5)
    dev.off()

    # Only histone marks:
    wm.rows <- marks[marks %in% rownames(wm)]
    wm.marks <- wm[wm.rows,]
    cols <- (seq(0,nrow(wm.marks),5) - .5)/ (nrow(wm.marks) -1)
    rows <- (seq(0,ncol(wm.marks),50) - .5)/ (ncol(wm.marks) -1)

    png(paste0(img,'/',procdir,'/',i,'-',nam,'_marks.png'),res=450,units='in',width=11,height=17)
    par(mar=c(8,4,3,1))
    image(wm.marks, axes=FALSE,main='',col=matcol)
    mtext(paste0(desc,' - ',day), side=3, cex=2.25, line=.75)
    mtext(paste0('Biosamples (', ncol(wm.marks), ')'), side=2, cex=2, line=1)
    mtext(paste0('Histone Marks, DNase-seq, and CTCF (', nrow(wm.marks), ')'), side=1, cex=2, line=6)
    # image(wm.marks, axes=FALSE,main=paste0(desc,' - ',day))
    grid(nx=nrow(wm.marks), ny=ncol(wm.marks),col='grey',lty='solid', lwd=.2)
    text(x=seq(0,1,length.out=nrow(wm.marks)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(wm.marks), srt=90, adj=1, xpd=TRUE,cex=.85)
    # text(y=seq(0,1,length.out=ncol(wm.marks)), x=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=colnames(wm.marks), srt=0, adj=1, xpd=TRUE,cex=.35)
    abline(h=rows,col='black',lty='solid',lwd=.5)
    abline(v=cols,col='black',lty='dashed',lwd=.5)
    dev.off()
}

# ====================================
# Finally cut off at different depths:
# what happens if we reject at different sample depths?
# ====================================
depths <- c(0,30,45,60,75)
for (depth in depths){
    i = 6
    nam <- abv[i]
    desc <- descriptions[i]
    wm <- as.matrix(wide[[nam]])
    # Modify depending on plot:
    matcol=c('darkblue','red')

    # Only histone marks:
    wm.rows <- marks[marks %in% rownames(wm)]
    wm.marks <- wm[wm.rows,]
    cols <- (seq(0,nrow(wm.marks),5) - .5)/ (nrow(wm.marks) - 1)
    rows <- (seq(0,ncol(wm.marks),50) - .5)/ (ncol(wm.marks) - 1)

    # Depth cutoff: 
    wm.marks[wm.marks < depth] <- 0
    wm.marks[wm.marks >= depth] <- 1
    idx <- which(colSums(wm.marks,na.rm=T) > 0)
    wm2 <- cbind(wm.marks[,idx],wm.marks[,-idx])
    total.num <- length(idx)
    cutcell <- (total.num - .5)/ (ncol(wm.marks) -1)

    # TODO ADD DESCRIPTOR OF PERCENT DONE.
    titlestr <- paste0(total.num,' cell types at depth ',depth,' - ',day)
    png(paste0(img,'/',procdir,'/Depth_',depth,'_marks.png'),res=450,units='in',width=11,height=17)
    par(mar=c(8,4,3,1))
    image(wm2, axes=FALSE,main='',col=matcol)
    mtext(titlestr, side=3, cex=2.25, line=.75)
    mtext(paste0('Biosamples (', ncol(wm2), ')'), side=2, cex=2, line=1)
    mtext(paste0('Histone Marks, DNase-seq, and CTCF (', nrow(wm.marks), ')'), side=1, cex=2, line=6)
    grid(nx=nrow(wm2), ny=ncol(wm2),col='grey',lty='solid',lwd=.25)
    text(x=seq(0,1,length.out=nrow(wm2)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(wm2), srt=90, adj=1, xpd=TRUE,cex=.85)
    # text(y=seq(0,1,length.out=ncol(wm2)), x=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels=colnames(wm2), srt=0, adj=1, xpd=TRUE,cex=.35)
    abline(h=rows,col='black',lty='solid',lwd=.5)
    abline(h=cutcell,col='red',lty='dotted',lwd=3)
    abline(v=cols,col='black',lty='dashed',lwd=.5)
    dev.off()
}


# ================
# Adding in IHEC: 
# ================
add <- read.delim('Annotation/IHEC_bw_histonemods.tsv',header=F)
ihec <- add[,c(2,6,7)]
ihec2 <- t(apply(ihec,1,function(x){
          if (length(grep(x[3],x[2])) > 0) {
              x[2] <- sub(paste0('_',x[3]),'',x[2])
          }
          return(x) } 
))

ihec <- data.frame(ihec2)
names(ihec) <- c('Project','CellType','Epitope')
ihec$Epitope <- as.character(ihec$Epitope)
ihec$Epitope[ihec$Epitope == 'Input'] <- 'WCE'
ihec$CellRepl <- 1
mat$Project <- 'ENCODE/Roadmap'
rnam <- c('CellType','Epitope','CellRepl','Project')
full <- rbind(mat[,rnam],ihec[,rnam])
nameproj = unique(full[,c('CellType','Project')])
rownames(nameproj) = nameproj$CellType

# Reordering by availability: 
cells <- aggregate(CellRepl ~ CellType,full,sum)
epitopes <- aggregate(CellRepl ~ Epitope,full,sum)
ordCT <- as.character(with(cells,CellType[order(CellRepl,decreasing=TRUE)]))
ordE <- as.character(with(epitopes,Epitope[order(CellRepl,decreasing=TRUE)]))

tmp <- dcast(full,Epitope ~ CellType,value.var='CellRepl', fun.aggregate=sum)
rownames(tmp) <- tmp$Epitope 
wihec <- tmp[ordE,ordCT]

# ===============
# Plot with IHEC:
# ===============
wm <- as.matrix(wihec)
wm.rows <- marks[marks %in% rownames(wm)]
wm <- wm[wm.rows,]
wm[wm > 0] <- 1

nam <- 'Availability'
matcol=c('white','darkblue')

# Only histone marks:
cols <- (seq(0,nrow(wm),5) - .5)/ (nrow(wm) -1)
rows <- (seq(0,ncol(wm),50) - .5)/ (ncol(wm) -1)

# Order and cut by project:
project <- nameproj[colnames(wm),2]
reord <- order(project, decreasing=TRUE)
uniq.proj <- unique(project[reord])
wm <- wm[,reord]

# Positioning lines / text:
counts <- rle(nameproj[colnames(wm),2])
cs = cumsum(counts$lengths) 
mids <- counts$lengths / 2 + c(0, cs)[-(length(cs) +1)]
cutcell <- (cs - .5) / (ncol(wm) - 1)
textpos <- (mids -.5) / (ncol(wm) - 1)
xpos = rep(c(0.05,.5),length(cs))[1:length(cs)]
projtext = paste0(uniq.proj, ' (', counts$lengths, ')')

png(paste0(img,'/',procdir,'/',nam,'_full_marks.png'),res=450,units='in',width=11,height=17)
par(mar=c(8,4,3,1))
image(wm, axes=FALSE,main='',col=matcol)
mtext(paste0(ncol(wm),' biosamples from ENCODE3 and IHEC - ',day), side=3, cex=2.25, line=.75)
mtext(paste0('Biosamples (', ncol(wm2), ')'), side=2, cex=2, line=1)
mtext(paste0('Histone Marks, DNase-seq, and CTCF (', nrow(wm.marks), ')'), side=1, cex=2, line=6)
text(x=seq(0,1,length.out=nrow(wm)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(wm), srt=90, adj=1, xpd=TRUE,cex=.85)
text(x=xpos, y=textpos, labels=projtext, srt=0, adj=0, xpd=TRUE, col=alpha('darkgrey',.85), cex=4)
abline(h=cutcell,col='red',lty='dotdash',lwd=1.4)
abline(h=rows,col='black',lty='solid',lwd=.5)
abline(v=cols,col='black',lty='dashed',lwd=.5)
dev.off()
