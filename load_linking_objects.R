#!/usr/bin/R
# ---------------------------------------
# Load objects for linking
# ------------------------------------------------
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))

library(rhdf5)
library(GenomicRanges)
library(Matrix)
library(scales)
library(plotrix)

#' Calculate position relative to par()$usr 
#'
#' @param axis 1 is x; 2 is y;
#' @param shift percent of axis to shift over
#' @return position shifted from start of x or y axis
#' @export
parpos = function(axis, shift){
    # NOTE: par()$usr is x1, x2, y1, y2
    if (axis == 1) { # X axis 
        par()$usr[1] - shift * diff(par()$usr[1:2])
    } else { # Y axis
        par()$usr[3] - shift * diff(par()$usr[3:4])
    }
}

# For plotting multiple arrows on a segment:
multarrows <- function(x0, y0, x1, y1, n_arr, ...) {
    x <- seq(x0, x1, length=n_arr + 1)
    y <- seq(y0, y1, length=n_arr + 1)
    arrows(x[-length(x)], y[-length(y)], x[-1], y[-1], ...) }

# ----------------------------------------
# Make sub-directories for data and plots:
# ----------------------------------------
gwasdir = paste0(img, "linking/gwas/")
cmd = paste('mkdir -p', gwasdir)
system(cmd)

# -------------------------------------------
# Read and process the GTF for gene plotting:
# -------------------------------------------
parse.gtf.info = function(x){
    x = strsplit(x, "; ")[[1]]
    x = unlist(strsplit(x, " "))
    out = x[seq(2,length(x), by=2)]
    names(out) = x[seq(1,length(x), by=2)] 
    return(out) }


gtf.file = 'Annotation/gencode.v30lift37.basic.annotation.gtf.gz'
gtf.rdafile= 'Annotation/gencode.v30lift37.basic.annotation.Rda'
if (!file.exists(gtf.rdafile)){
    tab = read.delim(gtf.file, sep="\t", skip=5, header=F, stringsAsFactors=F)
    # Find longest tx matching tss in each gene:
    txdf = tab[tab[,3] == 'transcript',c(1,4,5,7,9)]
    names(txdf) = c('chr','start','end','strand','info')
    txdf$tss = ifelse(txdf$strand == '+', txdf$start, txdf$end)
    tlist = sapply(as.character(txdf$info), parse.gtf.info)
    txdf$txid = sapply(tlist, function(x){x[['transcript_id']]})
    txdf$symbol = sapply(tlist, function(x){x[['gene_name']]})
    txdf$txtype = sapply(tlist, function(x){x[['transcript_type']]})
    txdf$length = txdf$end - txdf$start
    # NOT possible:
    txdf = txdf[txdf$txtype == 'protein_coding',]
    # Reduce to longest transcript only:
    txdf = merge(txdf, aggregate(length ~ symbol, txdf, max))
    txdf$tss = with(txdf,start * (strand == '+') + end * (strand == '-'))
    # Match exons to these transcripts:
    exdf = tab[tab[,3] == 'exon',c(1,4,5,7,9)]
    names(exdf) = c('chr','start','end','strand','info')
    exdf = exdf[grep('protein_coding', exdf$info),]
    elist = sapply(as.character(exdf$info), parse.gtf.info)
    exdf$txid = sapply(elist, function(x){x[['transcript_id']]})
    exdf$exnum = sapply(elist, function(x){x[['exon_number']]})
    exdf$exnum = as.numeric(exdf$exnum)
    exdf = exdf[exdf$txid %in% txdf$txid,]
    exdf = merge(exdf, txdf[,c('txid','symbol')])
    # Match UTRs to these transcripts:
    uxdf = tab[tab[,3] == 'UTR',c(1,4,5,7,9)]
    names(uxdf) = c('chr','start','end','strand','info')
    uxdf = uxdf[grep('protein_coding', uxdf$info),]
    ulist = sapply(as.character(uxdf$info), parse.gtf.info)
    uxdf$txid = sapply(ulist, function(x){x[['transcript_id']]})
    uxdf = uxdf[uxdf$txid %in% txdf$txid,]
    uxdf = merge(uxdf, txdf[,c('txid','symbol')])
    save(exdf, txdf, uxdf, file=gtf.rdafile)
} else {
    load(gtf.rdafile)
}


# ---------------------
# Load expression data:
# ---------------------
lddir = 'linking_data/'
# TODO: LOAD THE v28lif37???
tssfile = 'Annotation/gencode.gene.tss.v30lift37.basic.annotation.gtf.gz'
matfile = paste0(lddir, 'merged_log2fpkm.pc.mtx.gz')
generdafile = paste0(lddir, 'gene_mat_tss.Rda')
gmapfile = 'Annotation/gencode.gene.name.mapping.tsv'

if (!file.exists(generdafile)){
    tssdf = read.delim(gzfile(tssfile), header=F, sep="\t", stringsAsFactors=F)
    names(tssdf) = c('chr','tss','tss2','gene')
    tssdf = tssdf[!(tssdf$chr %in%  c('chrY', 'chrM')),]
    # Matrix to wide:
    matlong = read.delim(gzfile(matfile), header=T)
    matwide = spread(matlong, id, log2fpkm)
    gmat = as.matrix(matwide[,-1])
    rownames(gmat) = matwide[,1]
    rm(matwide, matlong)
    # Subset to shared genes:
    kept.genes = sort(tssdf$gene[tssdf$gene %in% rownames(gmat)])
    tssdf = tssdf[tssdf$gene %in% kept.genes,]
    tssdf$gene = factor(tssdf$gene, levels=kept.genes)
    tssdf = tssdf[order(tssdf$gene), ]
    gmat = gmat[kept.genes,]
    gc()
    save(gmat, tssdf, file=generdafile)
} else {
    load(generdafile)
}

genemap = read.delim(gmapfile, header=F)
names(genemap) = c('gene','symbol')
anno = merge(tssdf, genemap) # Note: not same as tssdf file.


# --------------------------------------
# Read in the dhs locations and indices:
# --------------------------------------
lddir = 'linking_data/'
mapfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r200_e0_names.core.srt.tsv'
locfile = 'DHS_Index_WM201902/masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt'
indfile = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_masterlist_indices.tsv'
marknamesfile = paste0(lddir, 'mark_matrix_names.txt')
dhsrdafile = paste0(lddir, 'dhs_loc.Rda')

# Load data:
if (!file.exists(dhsrdafile)){
    # Locations:
    locmap = read.delim(mapfile, header=T, stringsAsFactors=F)
    locdf = read.delim(locfile, header=F, stringsAsFactors=F)
    colnames(locdf) = c('chr','start','end','name')
    head(locmap)
    head(locdf)
    # Merge mapping:
    locdf = merge(locmap, locdf)
    locdf = locdf[order(locdf$cls),]
    locdf$mid = (locdf$start + locdf$end) / 2
    rm(locmap)
    gc()
    # Matrix names + indices:
    mnames = scan(marknamesfile,'c')
    matind = as.numeric(scan(indfile, 'c')) + 1
    locdf = locdf[matind,]
    # Save for loading:
    save(mnames, locdf, matind, file=dhsrdafile)
} else {
    load(dhsrdafile)
}
gc()

# -------------------------------------
# Read in the precomputed correlations:
# -------------------------------------
corrfile = paste0(lddir, 'mark_precomputed_corr.Rda')  # For fast reload
if (!file.exists(corrfile)){
    h1 = read.delim(gzfile(paste0(lddir, 'H3K27ac_precomputed_corr.tsv.gz')), header=F)
    h1 = as.numeric(h1$V1)
    h2 = read.delim(gzfile(paste0(lddir, 'H3K4me1_precomputed_corr.tsv.gz')), header=F)
    h2 = as.numeric(h2$V1)
    save(h1, h2, file=corrfile)
} else {
    load(corrfile)
}

# TODO: Get the TSS from here.
# Data files:
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K4me1_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'




