#!/usr/bin/R
# Plot reads distribution (for thresholding) for ENCODE data: 
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')

# Read in depths of different file: 
marks <- c('WCE',scan('Annotation/histone_mark_list','c'),'CTCF')
dist <- list()
for (mark in marks){ 
    fname = paste0('ChIP-seq/files/',mark,'/qc/',mark,'_read_totals.tsv')
    if (file.exists(fname)){
        depths <- as.numeric(scan(fname,'c')) / 10^6
        dist[[mark]] <- depths[depths != 0] 
    }
} 
all.depths <- unlist(dist)
cutoffs <- c(30,45,60,75) 

rn <- c(0,cutoffs,max(all.depths) + 1)
loc <- rn[-length(rn)] + mean(diff(cutoffs)) / 2
wl <- hist(all.depths,breaks=rn, plot=FALSE)
wlabs <- wl$counts

# Plot total: 
w200 <- all.depths[all.depths <= 200]
png(paste0(img,'total_reads_distr.png'),res=300,width=9,height=6,units='in')
h <- hist(w200,50,plot=FALSE)
plot(h,col='grey',border='black',xlim=c(0,200),main='All Pooled Samples',ylab='Count',xlab='Read Depth (millions)')
abline(v=cutoffs,col='blue',lwd=2,lty='dashed')
text(x=loc,y=.85 * max(h$counts), labels=wlabs,col='darkred',cex=1.3)
dev.off()

# =============
# Add in DNase:
# =============
mark <- 'DNase-seq'
fname = paste0(mark,'/files/qc/',mark,'_read_totals.tsv') # by experiment!
if (file.exists(fname)){
    depths <- as.numeric(scan(fname,'c'))
    dist[[mark]] <- depths[depths != 0] 
}

rn <- c(0,cutoffs,max(unlist(dist)) + 1) # Max will have changed with DHS
for (mark in names(dist)){
    depths = dist[[mark]]
    wl <- hist(depths,breaks=rn, plot=FALSE)
    wlabs <- wl$counts

    w200 <- depths[depths <= 200]
    png(paste0(img,mark,'_reads_distr.png'),res=300,width=9,height=6,units='in')
    h <- hist(w200,50,plot=FALSE)
    plot(h,col='grey',border='black',xlim=c(0,200),main=paste0(mark,' Pooled Samples'),ylab='Count',xlab='Read Depth (millions)')
    abline(v=cutoffs,col='blue',lwd=2,lty='dashed')
    text(x=loc,y=.85 * max(h$counts), labels=wlabs,col='darkred',cex=1.3)
    dev.off()
}

