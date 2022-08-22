#!/usr/bin/R
# ------------------------------
# Plot the chromHMM frequencies 
# for 15 and 18 state models for
# different cells.
# ------------------------------
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

# Functions:
perc.rank <- function(x) trunc(rank(x))/length(x)

plot.freq <- function(frac, avail, epitopes=NULL, include.names=TRUE, maxy=1, include.title=TRUE){
    # Set up area:
    sp =0.1
    if (maxy == 1){
        scaleup = 1.5
        layout(matrix(c(1,2),2,1), heights=c(3,.5 + 1 * include.names), TRUE)
    } else {
        scaleup = 1.5 * maxy
        layout(matrix(c(1,2),2,1), heights=c(3 * maxy,.5 + 1 * include.names), TRUE)
    }
    par(yaxs="i")
    par(xaxs="i")
    cols <- (seq(0, nrow(avail), 100) - .5)/ (nrow(avail) - 1)
    pr = barplot(frac, col=dcol$col, border='NA', las=1.5, cex.axis=.5, plot=FALSE, ylim=c(0,maxy))
    scale = (max(pr) + min(pr)) / nrow(avail)
    bcols = seq(0, nrow(avail), 100) * scale
    # Plot frequencies:
    par(mar=c(sp,4,sp + 2.5 * include.title,sp))
    colnames(frac) = NULL
    barplot(frac, col=dcol$col, border='NA', las=1.5, cex.axis=.5, ylim=c(0,maxy))
    abline(v=bcols, col='grey',lty='dashed',lwd=.5)
    abline(h=seq(.1, maxy, .1),col='grey',lty='dashed',lwd=.5)
    if (include.title){ mtext(paste0('State Frequencies (model: ',model, ')'), side=3, cex=1.2 * scaleup, line=.5)}
    if (!is.null(epitopes)){
        keep.cols = colnames(avail) %in% epitopes
        avail = avail[,keep.cols]
    }
    rows <- (seq(0, ncol(avail), 5) - .5)/ (ncol(avail) - 1)
    # Plot sample sheet:
    par(mar=c(2 * scaleup + 6 * include.names, 4, sp, sp))
    image(avail, axes=FALSE,col=c('white','indianred'), useRaster=TRUE)
    abline(v=cols,col='grey',lty='dashed',lwd=.5)
    abline(h=rows,col='grey',lty='dashed',lwd=.5)
    lbcex = 1 * scaleup - .25 * scaleup * (!is.null(epitopes))
    mtext(paste0('Assays (', ncol(avail), ')'), side=2, cex=lbcex, line=1 + 1.5 * scaleup)
    mtext(paste0('Epigenomes (', nrow(avail), ')'), side=1, cex=lbcex, line=.5 * scaleup + 4 * scaleup * include.names)
    text(y=seq(0,1,length.out=ncol(avail)), x=par()$usr[1]-0.002*(par()$usr[2]-par()$usr[1]), labels=colnames(avail), srt=0, adj=1, xpd=TRUE,cex=.35)
    if (include.names){
        text(x=seq(0,1,length.out=nrow(avail)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(avail), srt=90, adj=1, xpd=TRUE,cex=.25)
    }
}

# Directories:
imgpref = paste0(img, "model_diagnostics/")
chmmdir = 'ChromHMM/'
freqdir = paste0(chmmdir, "freqs/")
system(paste('mkdir -p', imgpref))

# Annotations:
map <- read.delim('Annotation/all_submitted_released_biosample_mapping_20180924.tsv',header=F)
rownames(map) <- map[,2]
names(map) <- c('celltype','BSSID')
samplesheet <- read.delim('ChromImpute/sample_mark_table.tsv',header=F)
names(samplesheet) <- c('BSSID','epitope','file')
samplesheet <- merge(samplesheet, map)
samplesheet$file = 1
widesample = spread(samplesheet, epitope, file)
widesample[is.na(widesample)] <- 0
rownames(widesample) = widesample$BSSID

# Models that have frequency files:
models = list.files(path=freqdir)
for (model in models){
    modeldir = paste0(freqdir, model,'/')
    tabnames = list.files(path=modeldir,pattern='table*')

    # Read in frequency files:
    tablist = list()
    tabdf = NULL
    for (file in tabnames){
        cell = sub(".txt","",sub("table_","",file))
        filename = paste0(modeldir, file)
        if (file.size(filename) > 0){
            tablist[[cell]] = read.delim(filename, sep=" ", header=F)
            names(tablist[[cell]]) <- c('state', cell)
            rownames(tablist[[cell]]) = tablist[[cell]]$state
            if (is.null(tabdf)){ 
                tabdf = data.frame(state=1:max(tablist[[cell]]$state))
                states = data.frame(state=1:max(tablist[[cell]]$state))
            }
            tmpdf = merge(states, tablist[[cell]], all.x=TRUE)
            tabdf[[cell]] = tmpdf[[cell]]
        }
    }
    tabdf[is.na(tabdf)] <- 0
    write.table(tabdf, paste0(modeldir, 'state_counts_', today, '.tsv'), sep="\t", row.names=F, col.names=T)


    modelpref = sub('_on_.*$','',model)
    modelpref = sub('_QCUT$','',modelpref)
    # Plot state counts as bars:
    dcol <- read.table(paste0('CHMM_', modelpref, '_colors.tsv'),header=F)
    names(dcol) <- c('state','name','color')
    dcol$col <- sapply(dcol$color,function(x){
                           s <- unlist(strsplit(as.character(x),split=','))
                           s <- as.numeric(s)
                           rgb(s[1],s[2],s[3],max=255) })
    mat = as.matrix(tabdf[,-1])
    ord = order(mat[nrow(mat),], decreasing=T)
    mat = mat[,ord]
    frac = sweep(mat,2,apply(mat,2,sum), '/')
    bssids = colnames(frac)
    bssids = bssids[bssids %in% cellorder]
    frac = frac[,bssids]
    # Average coverage:
    # boxplot(apply(frac[-18,], 2, sum))
    # mean(apply(frac[-18,], 2, sum))

    # ----------------------------------------
    # Boxplots per state with jittered counts:
    # ----------------------------------------
    fdf = data.frame(t(frac))
    fdf$id = rownames(fdf)
    fdf = gather(fdf, state, fraction, -id)
    fdf$state = sub("X", "", fdf$state)
    fdf = merge(fdf, dcol)
    fdf = merge(fdf, meta[,c('id','GROUP','COLOR')])
    # fdf = fdf[fdf$state != 18,]

    # png(paste0(imgpref, 'state_fraction_boxplots_', model, '.png'), res=450, width=8, height=3.5, units='in')
    pdf(paste0(imgpref, 'state_fraction_boxplots_', model, '.pdf'), width=8, height=3.5)
    par(mar=c(3,5,.5,.5))
    par(yaxs='i')
    par(xaxs='i')
    par(xpd=FALSE)
    ylim = c(-0.005,max(fdf$fraction) * 1.1)
    NS = max(as.numeric(fdf$state))
    atpar = 1:NS
    fdf$state = factor(fdf$state, levels = as.character(1:NS))
    boxplot(fraction ~ state, data=fdf,
            at = rev(atpar), width=rep(.01,NS), pch=19, cex=.25,
            ylim=ylim, axes=FALSE, col=NA, xaxs=FALSE, 
            xlab='',ylab='', lwd=.5, pars=list(outpch=''),
            horizontal=TRUE)
    ataxis = c(seq(0,0.20,.05), .50, .75, 1.00)
    axis(1, at=ataxis, labels=paste0(round(ataxis * 100,0), '%'))
    abline(v=ataxis, lty='dashed')
    box()
    points(fdf$fraction, NS + 1 - jitter(as.numeric(fdf$state), 2), 
           # pch=19, col=sapply(fdf$COLOR, tsp.col), cex=.1)
           pch=19, col=fdf$COLOR, cex=.1)
    boxplot(fraction ~ state, data=fdf,
            at = rev(atpar), width=rep(.01,NS), pch=19, cex=.25,
            ylim=ylim, axes=FALSE, col=sapply(dcol$col, tsp.col), 
            xaxs=FALSE, xlab='',ylab='', 
            lwd=.5, pars=list(outpch=''),
            horizontal=TRUE, add=TRUE)
    mtext('Percent of Epigenome', side=1, line=2)
    text(x=par()$usr[1] - 0.01 * diff(par()$usr[1:2]), y=rev(atpar),
         label=dcol$name, col=c(dcol$col[-18], 'black'), xpd=TRUE, adj=1)
    dev.off()

    # For stats:
    aggregate(fraction ~ state, fdf, function(x){round(mean(x) * 100,2)})

    # colnames(frac) = map[bssids, 1]
    avail = as.matrix(widesample[bssids,-c(1,2)])
    # rownames(avail) = map[rownames(avail), 1]

    # Reord avail:
    dt <- dist(t(avail),'manhattan')
    ht <- hclust(dt)
    cocl = ht$order
    colord <- colnames(avail)[cocl]
    avail = avail[,colord]

    # Subset of most important epitopes
    eplist = c('DNase-seq','H3K4me3','H3K36me3',
               'H3K4me1', 'H3K9me3', 'H3K27ac',
               'H3K27me3', 'H3K9ac', 'ATAC-seq')

    png(paste0(imgpref, 'state_counts_', model, '.png'), res=300, width=17, height=11, units='in')
    plot.freq(frac, avail)
    dev.off()

    png(paste0(imgpref, 'state_counts_', model, '_clean.png'), res=300, width=17, height=11, units='in')
    plot.freq(frac, avail, include.names=FALSE, epitopes=eplist)
    dev.off()

    fc = frac[nrow(frac),]
    fc = fc[fc > .1]
    top = 1 - min(fc)
    png(paste0(imgpref, 'state_counts_', model, '_clean_cut.png'), res=300, width=12, height=3, units='in')
    plot.freq(frac, avail, include.names=FALSE, epitopes=eplist, maxy = max(0.4, top))
    dev.off()

    pdf(paste0(imgpref, 'pdf/state_counts_', model, '.pdf'), width=17, height=11)
    plot.freq(frac, avail)
    dev.off()

    pdf(paste0(imgpref, 'pdf/state_counts_', model, '_clean.pdf'), width=9, height=4)
    plot.freq(frac, avail, include.names=FALSE, epitopes=eplist, include.title=FALSE, maxy=max(0.4,top))
    dev.off()

    # Make the emissions figure:
    if (modelpref == 'observed_aux_18'){
        # Get emissions:
        emfile = 'emissions_18_core_K27ac.txt'
        emdf = read.delim(emfile)
        emat = as.matrix(emdf[,-1])
        rownames(emat) = emdf[,1]
        # Get mnemonics:
        mfile = 'CHMM_observed_aux_18_mnemonic.tsv'
        mndf = read.delim(mfile, header=F)
        names(mndf) = c('state','name','string','colstr','color')
        dcol = merge(dcol, mndf)
        dcol = dcol[order(dcol$state),]
        rownames(dcol) = NULL
        pltord = rev(1:18)

        # Labels in figure:

        mot.col_fun = function(x, palette=col1){
            bin <- cut(x, seq(0, 100, length.out=length(palette)), include.lowest=T) 
            palette[bin] }
        mot.legend = Legend(at=seq(0,100, 25), 
                            labels_gp = gpar(fontsize=5),
                            title_gp = gpar(fontsize=5, fontface='bold'),
                            col_fun=mot.col_fun, title_position = "topleft", 
                            title='Emission %', direction = 'horizontal')
        plegend = packLegend(mot.legend)

        pdf(paste0(imgpref, 'pdf/', modelpref, '_emissions.pdf'), width=3.25, height=2.9)
        sp=0
        bsp=4.5
        layout(matrix(1:3, 1,3), widths=c(1,1,2), heights=3.5, TRUE)
        # Text Mnemonic:
        par(mar=c(bsp,sp,sp,sp))
        image(t(as.matrix(pltord)), col=dcol$col, axes=F, useRaster=TRUE)
        yat = seq(0, 1, length.out=18)
        text(-.95, yat, paste(pltord, as.character(dcol$name)[pltord]), adj=0)
        # Emissions:
        par(mar=c(bsp,sp,sp,sp))
        image(t(emat[pltord,]), col=col1, axes=F, useRaster=TRUE)
        text(x=seq(0,1, length.out=ncol(emat)),
             y=par()$usr[3]-0.01 * diff(par()$usr[3:4]),
             colnames(emat), adj=1, srt=90, xpd=TRUE)
        box(lwd=.5)
        # String name
        par(mar=c(bsp,sp,sp,sp))
        image(t(as.matrix(pltord)), col=dcol$col, axes=F, useRaster=TRUE)
        text(0, yat, as.character(dcol$string)[pltord], adj=.5)
        draw(plegend, x = unit(2.5,'in'), y=unit(.5,'in'), just = "top")
        dev.off()
    }

}



