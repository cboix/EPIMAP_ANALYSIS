#!/usr/bin/R
# -------------------------------------------------
# Plot obs/imp tracks according to a specific order/
# + dendrogram
# Plot observed vs. imputed tracks given: 
# 1) Group # to plot
# 2) Chr 
# 3) Start (default=0)
# 4) End (default=100MB)
# 5) Mark to plot
# -------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    source('~/data/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
} else {
    source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
}
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")

# Defaults
chr = 'chr10'
start = 0
end = 10 * 10^6
end = 100 * 10^6
markslist = NULL
# tsvfile='ChromImpute/sample_bed.bed'
# Arguments:
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need group number")
} else {        
    taskid = as.numeric(args[1])
    if (length(args) > 1){ chr = args[2] }
    if (length(args) > 2){ start = as.numeric(args[3]) }
    if (length(args) > 3){ end = as.numeric(args[4]) }
    if (length(args) > 4){ markslist = args[4] }
}

# Functions:
read.track = function(pref, slice){
    as.numeric(scan(gzfile(paste0(pref, '.wig.gz')), 'c', skip=2))[slice[1]:slice[2]]
}

plot.track = function(y, color='royalblue', cap=NULL, axlab=NULL){
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
        text(y=mean(yaxis),
             x=par()$usr[1]-0.0067*(par()$usr[2]-par()$usr[1]),
             labels=axlab, srt=0, adj=1, xpd=TRUE,cex=.65)
    }
    abline(h=0, col='black', lw=0.5)
}

# Tracks metadata:
imptab = read.delim('ChromImpute/impobs_table_fordist.tsv', header=F, stringsAsFactors=F)
names(imptab) = c('id','mark','file')
imptab = imptab[grep("FINAL_", imptab$file, invert=TRUE),]
obstab = read.delim('ChromImpute/sample_mark_table.tsv', header=F, stringsAsFactors=F)
names(obstab) = c('id','mark','file')
treeord = read.delim('tree_ordered_metadata.tsv', header=T, stringsAsFactors=F)
grouplist = sort(unique(treeord$group))

if (taskid == 0){
    group = 'ALL'
    ids = treeord$id
} else {
    group = grouplist[taskid]
    ids = treeord$id[treeord$group == group]
}
if (is.null(markslist)){ markslist = as.character(sort(unique(imptab[,2]))) }

tagline = paste0(group, ' (', length(ids),' samples)')
print(paste0('[STATUS] Plotting tracks for ', tagline))

# Submit plotting tracks from ALL + by metadata groups:
# For each mark (6-12) x different resolutions/areas (3)
# NOTE: PLOT IMPUTED EVEN WHEN THERE IS NO OBSERVED.
# NOTE: PLOT AVERAGE BY GROUP (OR BY DENDROGRAM CLUSTER?)

# Directories:
imgpref = paste0(img, "ordered_tree_tracks/")
system(paste('mkdir -p', imgpref))
convpref = paste0('ChromImpute/converted/', chr, '_')
imppref = paste0('ChromImpute/imputed/', chr, '_')

# Make tracks of only imputed, 
# and of both together.

# General Metadata, colors, mapping:
map <- read.delim('Annotation/all_submitted_released_biosample_mapping_20180924.tsv',header=F)
rownames(map) <- map[,2]
names(map) <- c('sample name','id')
rdcol = read.delim('Annotation/updated_tissue_colors.tsv', header=T, stringsAsFactors=F)
names(rdcol) = c('GROUP','COLOR')
rdcol = rdcol[order(rdcol$COLOR, decreasing=T),]

# Metadata for imputation datasets:
rd.map = read.delim('Annotation/updated_imputation_metadata.tsv',header=T, stringsAsFactors=F)
meta = rd.map[,c('BSSID','sample.name','newgroup','Extended.Info','origin','perturb','embryonic','age','sex','type')]
names(meta)[1:4] = c('id','ct','GROUP','infoline')
rdcol = rdcol[rdcol$GROUP %in% unique(meta$GROUP),]
meta = merge(meta, rdcol)
meta[is.na(meta)] <- ''
rdcol = rbind(rdcol, data.frame(GROUP='Modified',COLOR='lightgrey'))

# In case some aren't labeled:
ct.list = sort(unique(imptab[,1]))
meta = merge(meta, map[ct.list,], all.y=TRUE, all.x=TRUE)
meta$GROUP[is.na(meta$GROUP)] = 'Modified'
meta$COLOR[meta$GROUP == 'Modified'] = 'white'
meta$GROUP = factor(meta$GROUP, levels = as.character(rdcol$GROUP))
rownames(meta) = meta$id
meta$type[meta$type == 'NA'] = 'neither'
meta$type[is.na(meta$type)] = 'neither'
meta$sex[meta$sex == 'NA'] = 'neither'
meta$sex[meta$sex == ''] = 'neither'
meta$sex[is.na(meta$sex)] = 'neither'
meta$embryonic[meta$embryonic == 'NA'] = 'neither'
meta$embryonic[meta$embryonic == ''] = 'neither'
meta$embryonic[is.na(meta$embryonic)] = 'neither'

# Colors:
colvals = list()
colvals[['sex']] = c('male' = 'royalblue', 'female' = 'forestgreen', 'female and male' = 'lightgrey', 'neither' = 'lightgrey')
colvals[['type']] = c('TISSUE' = 'darkgoldenrod1', 'LINE' = 'royalblue', 'neither' = 'lightgrey')
colvals[['embryo']] = c('embryo' = 'forestgreen', 'adult' = 'royalblue', 'newborn' = 'indianred', 'child' = 'mediumpurple3', 'neither' = 'lightgrey')
colvals[['group']] = rdcol$COLOR
names(colvals[['group']]) = rdcol$GROUP
if (group == "ALL"){ gcol='black' } else { gcol=colvals[['group']][group] }

# Function to plot tracks from df.
plot.track.df = function(trackdf, chr, start, end, title=NULL, color='black'){
    N = nrow(trackdf)
    if (is.null(title)){
        guessmark = unique(trackdf$mark)
        guessid = unique(trackdf$id)
        if(length(guessmark) <= length(guessid)){
            title = paste0(guessmark, collapse=", ")
        } else {
            title = paste0(guessid, collapse=", ")
        }
    }
    slice = c(round(start / 25) + 1, round(end / 25))
    mat = matrix(0, ncol=N, nrow=slice[2] - slice[1] + 1)
    convpref = paste0('ChromImpute/converted/', chr, '_')
    imppref = paste0('ChromImpute/imputed/', chr, '_')
    # Load each track in:
    print("[STATUS] Loading tracks in...")
    for (i in 1:N){
        if (trackdf$type[i] == 'imputed'){pref = imppref} else {pref = convpref}
        track = as.numeric(scan(gzfile(paste0(pref, trackdf$file[i], '.wig.gz')),
                                'c', skip=2))[slice[1]:slice[2]]
        mat[,i] = track
    }
    # Plot the tracks:
    print("[STATUS] Plotting tracks...")
    layout(matrix(c(1:(N+2)),(N+2),1), heights=c(2,rep(1,N),3), TRUE)
    par(yaxs="i")
    par(xaxs="i")
    par(mar=c(0.05, 5, .05, 2))
    plot(c(0,1),c(0,1), col='white',axes=F, ylab='', xlab='')
    # Title:
    text(x=mean(par()$usr[1:2]), y=mean(par()$usr[3:4]), labels=title, cex=2.5, col=color)
    for (i in 1:N){
        track = mat[,i]
        par(mar=c(0.05, 5, .05, 2))
        quant = min(quantile(track,.998), 20)
        plot.track(track, color=as.character(trackdf$color[i]), cap=quant)
        # Label each axis:
        text(y=mean(par()$usr[3:4]),
             x=par()$usr[1]-0.0067*(par()$usr[2]-par()$usr[1]),
             labels=paste0(trackdf$id[i],'\n',trackdf$infoline[i]), srt=0, adj=1, xpd=TRUE,cex=.5)
        if(i == N) {
            sq = seq(start, end, length.out=6)
            sq2 = sq
            sq2[1] = start + 25
            loc = (sq2 - start) / (25)
            axis(side=1, at=loc, labels=FALSE)
        }
    }
    # Genomic locations:
    par(mar=c(2, 5, 0, 2))
    plot(c(0,1),c(0,1), col='white',axes=F, ylab='', xlab='')
    loc = (sq - start) / (25 * nrow(mat))
    text(x=loc, y=mean(par()$usr[3:4]), labels=sq, xpd=TRUE, cex=1)
    mtext(paste0(chr, ' Location (bp)'),side=1, line=0, cex=1)
}

# -----------------------
# Plot each mark in list:
# -----------------------
info = meta[ids, 'infoline']
for (mark in markslist){
    print(mark)
    # Make table, imputed
    trackdf = NULL
    for (id in ids){
        orow = merge(obstab, data.frame(id=id, mark=mark))
        irow = merge(imptab, data.frame(id=id, mark=mark))
        if(nrow(orow) > 0){ trackdf = rbind(trackdf, data.frame(orow, type='observed', color='royalblue')) }
        if(nrow(irow) > 0){ trackdf = rbind(trackdf, data.frame(irow, type='imputed', color='indianred')) }
    }
    trackdf = merge(trackdf, meta[,c('id','infoline')])
    trackdf$id = factor(trackdf$id, levels= ids)
    trackdf = trackdf[order(trackdf$id),]

    SPLIT = 150
    print(paste0("Plotting:", nrow(trackdf), " tracks."))
    if (nrow(trackdf) < SPLIT){
        # Plot all tracks:
        prefix = paste0(imgpref, chr, "_", start, "_", end, "_", mark, "_group", taskid)
        png(paste0(prefix, "_alltracks.png"), res=250, units='in',width=16,height=nrow(trackdf)/4)
        plot.track.df(trackdf, chr, start, end, title=paste0(tagline, "\t", mark), color= gcol)
        dev.off()

        # Plot only imputed tracks:
        imptrackdf = trackdf[trackdf$type == 'imputed',]
        prefix = paste0(imgpref, chr, "_", start, "_", end, "_", mark, "_group", taskid)
        png(paste0(prefix, "_imptracks.png"), res=250, units='in',width=16,height=nrow(imptrackdf)/4)
        plot.track.df(imptrackdf, chr, start, end, title=paste0(tagline,"\t", mark), color=gcol)
        dev.off()
    } else {
        print(paste0("Breaking image up by ", SPLIT," sized segments"))
        NR = nrow(trackdf)
        for (i in 1:(floor(NR/SPLIT) + 1)){
            print(i)
            selrows = ((i-1) * SPLIT + 1):(min(i * SPLIT, NR))
            subtrackdf = trackdf[selrows,]
            # Plot all tracks:
            prefix = paste0(imgpref, chr, "_", start, "_", end, "_", mark, "_group", taskid)
            png(paste0(prefix, "_alltracks_part", i, ".png"), res=250, units='in',width=16,height=nrow(subtrackdf)/4)
            plot.track.df(subtrackdf, chr, start, end, title=paste0(tagline, "\t", mark), color= gcol)
            dev.off()
        }
        # Plot only imputed tracks:
        imptrackdf = trackdf[trackdf$type == 'imputed',]
        NR = nrow(imptrackdf)
        for (i in 1:(floor(NR/SPLIT) + 1)){
            print(i)
            selrows = ((i-1) * SPLIT + 1):(min(i * SPLIT, NR))
            subtrackdf = imptrackdf[selrows,]
            prefix = paste0(imgpref, chr, "_", start, "_", end, "_", mark, "_group", taskid)
            png(paste0(prefix, "_imptracks_part",i,".png"), res=250, units='in',width=16,height=nrow(subtrackdf)/4)
            plot.track.df(subtrackdf, chr, start, end, title=paste0(tagline,"\t", mark), color=gcol)
            dev.off()
        }
    }
}

# NOTE: To see all marks together, collate as PDF.

print("Plotting finished.")
