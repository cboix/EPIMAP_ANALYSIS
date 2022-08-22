#!usr/bin/R
# ---------------------------------------------------------
# Calculate the positions of motifs for all 2.1M enhancers:
# ---------------------------------------------------------
library(epimapAUX)
set_proj('EPIMAP_ANALYSIS')
library(chromVARmotifs)
library(Matrix)
library(motifmatchr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))


motdir = paste0(dbdir, 'motif_matrices/')
system(paste('mkdir -p', motdir))

# 1. Load motif data:
data("human_pwms_v1")
data("encode_pwms")
lsos()

# 2. Load the enhancer locations as genomic ranges:
ddir = 'DHS_Index_WM201902/'
dpref = 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19'
dmlfile = paste0(ddir, dpref, '.core.srt.txt')
dmlrdafile = paste0(ddir, dpref, '_enh.core.srt.Rda')
load(dmlrdafile)
# To gr object:
enhgr = GRanges(enhdf$chr, IRanges(enhdf$start, enhdf$end))

# Preliminary test set for metabolism + phase separation:
tfset = c('YY1', 'RXRA', 'RARB', 'RARA',
          'PPARG', 'PPARA', 'PGC1A', 'NCOA1',
          'MED1', 'HNF4A', 'ESRRA', 'ESR1',
          'EP300', 'CREBBP', 'CREB1')

# Just interest:
matfile = paste0(motdir, 'metab_motifs.mtx')
matnamfile = paste0(motdir, 'metab_motifs_names.txt')
if (!file.exists(matfile)){
    tfnam = sapply(tfset, function(x){names(encode_pwms)[grep(x, names(encode_pwms))]})
    tfnam = unlist(tfnam)
    motif_ix <- matchMotifs(encode_pwms[tfnam], enhgr, genome = "hg19") 
    mat = motifMatches(motif_ix)
    writeMM(mat, matfile, quote=F)
    write.table(tfnam, matnamfile, quote=F, row.names=F, col.names=F)
} else { 
    mat = readMM(matfile)
    tfnam = scan(matnamfile, 'c')
    colnames(mat) = tfnam
}

dist.to.df = function(dt){
    jdf = data.frame(as.matrix(dt))
    jdf$M1 = rownames(jdf)
    jdf = gather(jdf, M2, value, -M1)
}


raw.jt = jacc.dist(x=t(mat))
jdf = dist.to.df(raw.jt)


# Raw:
epimapAUX:::plot.dt.sym(as.dist(raw.jt))

# Load in the enhancer sets:
flatset = 'epigenomes'
enhsetfile = paste0(flatset, '_enhancer_sets.Rda')
load(enhsetfile)
NF = length(enhsets)

# enh-mapping:
enhind = as.numeric(scan('Enhancer_indices.txt', 'c')) + 1
enhmap = rep(0, max(enhind))
enhmap[enhind] = 1:length(enhind)

# vs. other enhancers
fulldf = c()
for (i in 1:NF){
    em = enhmap[enhsets[[i]]]
    sub.jt = jacc.dist(x=t(mat[em,]))
    # epimapAUX:::plot.dt.sym(as.dist(sub.jt))
    subdf = dist.to.df(sub.jt)
    subdf$set = i
    fulldf = rbind(fulldf, subdf)
}

jtenhfile = paste0(motdir, 'metab_motifs_byenh.tsv')
write.table(fulldf, jtenhfile, quote=F, row.names=F, sep="\t")

mfdf = aggregate(value ~ M1 + M2, fulldf, mean)
sfdf = aggregate(value ~ M1 + M2, fulldf, sd)

max(sfdf$value)

# Figure out which interactions change - and in what cells
# mean interaction, sd of int. zscore + cell of action



# ------------------------------------------
# Get motif matches for motifs in enhancers:
# ------------------------------------------
# TODO: ALL matches:
for (i in 1:length(encode_pwms)){
    enam = names(encode_pwms)[i]
    motif_ix <- matchMotifs(encode_pwms[[i]], enhgr, genome = "hg19") 
    # TODO: write as we go:
}


