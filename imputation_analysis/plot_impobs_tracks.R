#!/usr/bin/R
# -------------------------------------------------
# Plot observed vs. imputed tracks given: 
# 1) Bedfile listing a set of observed tracks 
# 2) Chr 
# 3) Start (default=0)
# 4) End (default=100MB)
# -------------------------------------------------
# CONVERTED_FILE=${INPUTFILE}.wig.gz
# IMPUTED_FILE=impute_${SAMPLE}_${MARK}.wig.gz
# OUTPUT_FILE=${EVAL_DIR}/${SAMPLE}_${MARK}_eval.txt

domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    source('~/data/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
} else {
    source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
}
source('~/data/EPIMAP_ANALYSIS/bin/load_metadata.R')
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")

# Defaults
chr = 'chr10'
start = 0
end = 10 * 10^6
end = 100 * 10^6
label = TRUE
cap = NULL
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need table filename.")
} else {        
    tsvfile = args[1]
    if (length(args) > 1){ chr = args[2] }
    if (length(args) > 2){ start = as.numeric(args[3]) }
    if (length(args) > 3){ end = as.numeric(args[4]) }
    if (length(args) > 4){ label = as.logical(args[5]) }
    if (length(args) > 5){ cap = as.numeric(args[6]) }
}
tag = sub("\\.tsv","", gsub(".*/","", tsvfile))

# Functions:
read.track = function(pref, slice){
    as.numeric(scan(gzfile(paste0(pref, '.wig.gz')), 'c', skip=2))[slice[1]:slice[2]]
}

plot.track = function(y, color='royalblue', cap=NULL, 
                      axlab=NULL, sublabel=FALSE){
    # Plot polygon
    if (is.null(cap)){
        yaxis = c(0, max(y))
    } else {
        yaxis = c(0, cap)
    }
    if (length(y) < 100000){
        n = length(y)
        x = c(1, 1:n, n)
        y = c(0, y, 0)
        plot(range(x), yaxis, type='n', axes=F, ylab='', xlab='')
        polygon(x, y, col=color, border=NA)
    } else {
        plot(y, col=color, type='l', axes=F, cex.axis=.5, ylim=yaxis, ylab='')
    }
    if (!is.null(axlab)){
        if (sublabel){
            text(labels=axlab, y=mean(yaxis),
                 x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]),
                 srt=0, adj=1, xpd=TRUE, cex=.48, col='darkgrey')
        } else {
            text(labels=axlab, y=mean(yaxis),
                 x=par()$usr[1]-0.005*(par()$usr[2]-par()$usr[1]),
                 srt=0, adj=1, xpd=TRUE,cex=.65)
        }
    }
    abline(h=0, col='black', lw=0.5)
}

# Directories:
imgpref = paste0(img, "impobs_tracks/")
system(paste('mkdir -p', imgpref))
convpref = paste0('ChromImpute/converted/', chr, '_')
imppref = paste0('ChromImpute/imputed/', chr, '_')

# Annotations:
map <- read.delim('Annotation/all_submitted_released_biosample_mapping_20180924.tsv',header=F)
rownames(map) <- map[,2]
names(map) <- c('celltype','BSSID')

# Read in information:
df = read.delim(tsvfile, header=F, stringsAsFactors=F)
plot.diff=FALSE
if (ncol(df) == 3){
    names(df) <- c('BSSID','mark','file')
    df$impfile = paste0("impute_",df$BSSID,"_",df$mark)
} else if (ncol(df) == 4){
    names(df) <- c('BSSID','mark','file', 'impfile')
} else if (ncol(df) == 6){
    names(df) <- c('BSSID','mark','file', 'impfile', 'difffile', 'matchfile')
    plot.diff=TRUE
    df$matchlab = sub('_', '  ', sub('impute_', '', df$matchfile))
}
slice = c(round(start / 25) + 1, round(end / 25))
df$lab = paste0(df$BSSID, '  ', df$mark)

# Plot image (normal or difference image):
if (!plot.diff){
    N = nrow(df)
    mat = matrix(0, ncol=N * 2, nrow=slice[2] - slice[1] + 1)
    for (i in 1:N){
        print(i)
        obs = as.numeric(scan(gzfile(paste0(convpref, df$file[i], '.wig.gz')), 'c', skip=2))[slice[1]:slice[2]]
        imp = as.numeric(scan(gzfile(paste0(imppref, df$impfile[i], '.wig.gz')), 'c', skip=2))[slice[1]:slice[2]]
        mat[,2 * (i-1) + 1] = obs
        mat[,2 * (i-1) + 2] = imp
    }

    prefix = paste0(imgpref, chr, "_", start, "_", end, "_", tag)
    # png(paste0(prefix, "_alternating.png"), units='in', res=300, width=(1+3.25), height=N / 4)
    pdf(paste0(prefix, "_alternating.pdf"), width=(1+3.25), height=N / 4)
    layout(matrix(c(1:(4 * N)),N * 2, 2, byrow=TRUE), widths=c(1.75,7), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    sp = 0.01
    rsp = 0.01
    for (i in 1:N){
        obs = mat[,2 * (i-1) + 1]
        imp = mat[,2 * (i-1) + 2]
        gcol = meta[df$BSSID[i], 'COLOR']
        tcex=.85
        par(mar=rep(sp, 4))
        if (is.null(cap)){quant = max(min(quantile(obs,.998), 20),5)} else {quant=cap}
        # if (label) {lab = df$mark[i]} else {lab = NULL}
        image(as.matrix(1), axes=F, ylab='', xlab='', col=gcol)
        text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), 
             df$mark[i], cex=tcex, col='grey85')
        par(mar=rep(sp, 4))
        plot.track(obs, 'royalblue', cap=quant)
        # Second:
        par(mar=rep(sp, 4))
        if (is.null(cap)){quant = max(min(quantile(imp,.998), 20),5)} else {quant=cap}
        # if (label) {lab = df$BSSID[i]} else {lab = NULL}
        image(as.matrix(1), axes=F, ylab='', xlab='', col=gcol)
        text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), 
             meta[df$BSSID[i],'infoline'], cex=tcex * .6, col='grey85')
        par(mar=rep(sp, 4))
        plot.track(imp, 'indianred', cap=quant)
    }
    dev.off()

    plot.grouped=FALSE
    if (plot.grouped){
        png(paste0(prefix, "_grouped.png"), res=300, units='in',width=5, height=N / 4)
        # png(paste0(prefix, "_grouped.png"), res=300, units='in',width=16,height=2 * N / 4)
        layout(matrix(c(1:(2 * N)),N * 2, 1), TRUE)
        par(yaxs="i")
        par(xaxs="i")
        quant = min(quantile(mat,.998), 20)
        for (i in 1:N){
            obs = mat[,2 * (i-1) + 1]
            par(mar=c(.1,6.5,.1,1))
            if (is.null(cap)){quant = min(quantile(obs,.998), 20)} else {quant = cap}
            if (label) {lab = df$lab[i]} else {lab = NULL}
            plot.track(obs, 'royalblue', cap=quant, axlab=lab)
        }
        for (i in 1:N){
            imp = mat[,2 * (i-1) + 2]
            par(mar=c(.1,6.5,.1,1))
            if (is.null(cap)){quant = min(quantile(imp,.998), 20)} else {quant = cap}
            if (label) {lab = df$lab[i]} else {lab = NULL}
            plot.track(imp, 'indianred', cap=quant, axlab=lab)
        }
        dev.off()
    }

} else {
    # Plot the difference:
    N = nrow(df)
    mat = matrix(0, ncol=N * 4, nrow=slice[2] - slice[1] + 1)
    for (i in 1:N){
        print(i)
        obs = read.track(paste0(convpref, df$file[i]), slice)
        imp = read.track(paste0(imppref, df$impfile[i]), slice)
        diff = read.track(paste0(imppref, df$difffile[i]), slice)
        match = read.track(paste0(imppref, df$matchfile[i]), slice)

        mat[,4 * (i-1) + 1] = obs
        mat[,4 * (i-1) + 2] = imp
        mat[,4 * (i-1) + 3] = diff
        mat[,4 * (i-1) + 4] = match
    }

    prefix = paste0(imgpref, 'diff_', chr, "_", start, "_", end, "_", tag)
    # TODO: make different size??
    scale = 1.5
    png(paste0(prefix, "_alternating.png"), res=300, units='in',width=16,height=2 * N / 4 * scale)
    layout(matrix(c(1:(2*5 * N)),N * 5, 2, byrow=TRUE), widths=c(1,15), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    for (i in 1:N){
        obs = mat[,4 * (i-1) + 1]
        imp = mat[,4 * (i-1) + 2]
        diff = mat[,4 * (i-1) + 3]
        match = mat[,4 * (i-1) + 4]
        sp=0.05
        rsp=sp
        tcex=.7
        # t1
        par(mar=c(sp,rsp,sp,sp))
        quant = min(quantile(obs,.998), 20)
        image(as.matrix(1), axes=F, ylab='', xlab='', 
              col=ifelse(i%%2==1, 'grey25', 'grey85'))
        text(mean(par()$usr[1:2]), mean(par()$usr[3:4]),
             df$lab[i], cex=tcex, col=ifelse(i%%2==1, 'white', 'black'))
        par(mar=c(sp,rsp,sp,sp))
        plot.track(obs, 'royalblue', cap=quant)
        # t2
        par(mar=c(sp,rsp,sp,sp))
        quant = min(quantile(imp,.998), 20)
        image(as.matrix(1), axes=F, ylab='', xlab='', 
              col=ifelse(i%%2==1, 'grey25', 'grey85'))
        text(mean(par()$usr[1:2]), mean(par()$usr[3:4]),
             df$lab[i], cex=tcex, col=ifelse(i%%2==1, 'white', 'black'))
        par(mar=c(sp,rsp,sp,sp))
        plot.track(imp, 'indianred', cap=quant)
        # t3
        par(mar=c(sp,rsp,sp,sp))
        quant = min(quantile(diff,.998), 20)
        image(as.matrix(1), axes=F, ylab='', xlab='', 
              col=ifelse(i%%2==1, 'grey25', 'grey85'))
        text(mean(par()$usr[1:2]), mean(par()$usr[3:4]),
             df$lab[i], cex=tcex, col=ifelse(i%%2==1, 'white', 'black'))
        par(mar=c(sp,rsp,sp,sp))
        plot.track(diff, 'darkorchid4', cap=quant)
        #t4
        par(mar=c(sp,rsp,sp,sp))
        quant = min(quantile(match,.998), 20)
        image(as.matrix(1), axes=F, ylab='', xlab='', 
              col=ifelse(i%%2==1, 'grey25', 'grey85'))
        text(mean(par()$usr[1:2]), mean(par()$usr[3:4]),
             df$matchlab[i], cex=tcex, col=ifelse(i%%2==1, 'white', 'black'))
        par(mar=c(sp,rsp,sp,sp))
        plot.track(match, 'forestgreen', cap=quant)
        # t5 (BLANK)
        par(mar=rep(sp, 4))
        plot(1, type='n', ylab='', xlab='', axes=F)
        par(mar=rep(sp, 4))
        plot(1, type='n', ylab='', xlab='', axes=F)
    }
    dev.off()
}

print("Plotting finished.")
