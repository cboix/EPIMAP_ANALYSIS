#!/usr/bin/R
# ------------------
# Plot motifs logos:
# ------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'load_metadata.R'))
today <- format(Sys.time(), "%m%d%y")
# library(seqLogo) # Awful colors
# library(motifStack) # Awful interface
library(ggseqlogo) # Compromise...
library(ggplot2)
library(cowplot)

imgdir = paste0('motifs/img/')
cmd = paste0('mkdir -p ', imgdir)
system(cmd)

# -----------------------
# Read in the motif pfms:
# -----------------------
mtfile = 'motifs/collated_motifs.txt'
mtlines = readLines(mtfile)
motlist = list()
curr = ""
for (i in 1:length(mtlines)){
    line = mtlines[i]
    if (length(grep("^>", line)) > 0){
        curr = line
        motlist[[curr]] = c()
    } else {
        line = as.numeric(strsplit(line, " ")[[1]][2:5])
        motlist[[curr]] = cbind(motlist[[curr]], line)
    }
}

# ---------------------------------------------------------
# Plot each of the ones required for the motifs as png/pdf:
# ---------------------------------------------------------
kept.motifs = scan('motifs/kept_160repr_082720_ordered.txt', 'c')

headers = names(motlist)
mnames = sub(">", "", sub(" .*", "", headers))
names(motlist) = mnames
motifs = list()
for (i in 1:length(kept.motifs)){
    print(i)
    mnam = kept.motifs[i]
    mat = motlist[[mnam]]
    rownames(mat) = c('A','C','G','T')
    h = .5
    w = dim(mat)[2] / 4
    gplot = ggseqlogo(mat) + 
        theme_nothing() + 
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0))
    ggsave(paste0('motifs/img/logo_', mnam, '.png'), gplot, dpi=300, width=w, height=h, units='in')
}

# Arrange them all onto a pdf: (for now)

library(png)
library(stringr)

pdf('motifs/img/logo_160_ordered.pdf', width=5, height=20)
par(mar=rep(0,4))
plot(1,1, type='n', axes=F, ylab='', xlab='', xlim=c(0,1), ylim=c(0,1))
for (i in 1:length(kept.motifs)){
    print(i)
    mnam = kept.motifs[i]
    mat = motlist[[mnam]]
    rownames(mat) = c('A','C','G','T')
    h = .5
    w = dim(mat)[2] / 4
    pngfile = paste0('motifs/img/logo_', mnam, '.png')
    img <- try(readPNG(pngfile))
    xj = 0.25
    ysize = .008 * h * 2
    xsize = .5 * w * 2 * .1
    # Image is 600 x 800
    if (class(img) != 'try-error'){
        yat = seq(0,1, length.out=160)
        yj = yat[i]
        rasterImage(img, xleft=xj, xright=xj + xsize,
                    ybottom=yj-ysize, ytop=yj + ysize)
    }
}
dev.off()



# -----------------------------------
# Plot all motifs, for tables online:
# -----------------------------------
headers = names(motlist)
mnames = sub(">", "", sub(" .*", "", headers))
names(motlist) = mnames
motifs = list()
for (i in 1:length(mnames)){
    print(i)
    mnam = mnames[i]
    mat = motlist[[mnam]]
    rownames(mat) = c('A','C','G','T')
    h = .5
    w = dim(mat)[2] / 4
    gplot = ggseqlogo(mat) + 
        theme_nothing() + 
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0))
    ggsave(paste0('motifs/img/logo_', mnam, '.png'), gplot, dpi=300, width=w, height=h, units='in')
}


