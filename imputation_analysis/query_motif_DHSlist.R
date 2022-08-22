#!/usr/bin/R
# -------------------------------------------------
# Quickly query locations of a motif in a DHS list:
# -------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (length(domain) == 0) domain = ''
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))

library(motifmatchr)  # Loads slowly?
library(GenomicRanges)
library(rhdf5)

# Files for GenomicRanges objects:
dhsdir = 'DHS_Index_WM201902/'
dhspref = paste0(dhsdir, 'masterlist_DHSs_733samples_WM20180608_all_coords_hg19')
mapfile = paste0(dhspref, '_r200_e0_names.core.srt.tsv')
locfile = paste0(dhspref, '.core.srt.txt')
dgrfile = paste0(dhspref, '_GRanges.Rda')
indfile = 'masterlist_matindices/matindices_on_mixed_impobs_MODEL_observed_aux_18_on_mixed_impobs_QCUT_intersect_H3K27ac_ENH_masterlist_indices.tsv'
marknamesfile = 'linking_data/mark_matrix_names.txt'

if (!file.exists(dgrfile)){
    # Locations:
    locmap = read.delim(mapfile, header=T, stringsAsFactors=F)
    locdf = read.delim(locfile, header=F, stringsAsFactors=F)
    colnames(locdf) = c('chr','start','end','name')
    head(locmap)
    head(locdf)
    # Merge mapping:
    locdf = merge(locmap, locdf)
    locdf = locdf[order(locdf$cls),]
    rownames(locdf) = NULL
    rm(locmap)
    gc()
    # Matrix names + indices:
    mat.names = scan(marknamesfile,'c')
    enh.ind = as.numeric(scan(indfile, 'c')) + 1
    # Make a full and an enhancer GRanges objects:
    dhs.gr = with(locdf, GRanges(seqnames=chr, IRanges(start, end), name=name))
    enh.gr = with(locdf[enh.ind,], GRanges(seqnames=chr, IRanges(start, end), name=name))
    # Save for loading:
    save(dhs.gr, enh.gr, mat.names, enh.ind, file=dgrfile)
    rm(locdf)
    gc()
} else {
    load(dgrfile)
}


# load some example motifs
# data(example_motifs, package = "motifmatchr") 

# TODO: Generalize to any motif.
trps.pwm = read.delim('motif_matrices/TRPS1_pwm_homer.txt', skip=1, sep="\t", header=F)
colnames(trps.pwm) = c('A','C','G','T')
seqLogo::seqLogo(t(trps.pwm))

trps.pfm = TFBSTools::PFMatrix(ID="TRPS1_Homer", name="Trps1", 
                               matrixClass="GATA-like", strand="+",
                               bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                               tags=list(family="ZNF", tax_group="vertebrates"),
                               profileMatrix=t(trps.pwm) * 1000)

icm = TFBSTools::toICM(trps.pfm)
# Match in defined peaks:
motif_ix <- matchMotifs(trps.pfm, enh.gr, genome = "hg19") 
mat = motifMatches(motif_ix)
motif_loc = which(mat[,1])
# Alternatively, get motif positions:
# motif_pos <- matchMotifs(trps.pfm, enh.gr, 
#                          genome = "hg19", out = "positions") 

# Score the motifs by tissue:
H3K27ac.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'
H3K4me1.file = 'linking_data/H3K27ac_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'

# H3K27ac
h5f = H5Fopen(H3K27ac.file)
h5ls(h5f)
h5d = h5f&"matrix"
h27mat = h5d[,enh.ind[motif_loc]]
H5Dclose(h5d)
H5Fclose(h5f)
rownames(h27mat) = mat.names
hmat = t(h27mat)
rm(h27mat)
gc()


cmarg = colSums(hmat > 2)
cdf = data.frame(id=names(cmarg), marg=cmarg)
cdf = merge(cdf, meta[,c('id','infoline','GROUP','lifestage')])
cdf = cdf[order(cdf$marg, decreasing=T),]

head(cdf, 20)





