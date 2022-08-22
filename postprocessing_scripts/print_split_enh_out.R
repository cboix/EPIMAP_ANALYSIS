#!/usr/bin/R
# --------------------------------------
# Write BED files of enhancer locations:
# Updated 05/11/21
# --------------------------------------
library(rhdf5)

# Read in metadata:
metadata = read.delim('../public_metadata_released/main_metadata_table.tsv', header=T)
mmap = read.delim('mnemonic_mapping.tsv', header=F)
names(mmap) = c('id','mn')
rownames(mmap) = mmap$id

# Load DHS locations:
dhsdf = read.delim('masterlist_DHSs_733samples_WM20180608_all_coords_hg19.core.srt.txt', header=F)
eind = as.numeric(scan('ENH_masterlist_indices_0indexed.tsv','c')) + 1
pind = as.numeric(scan('PROM_masterlist_indices_0indexed.tsv','c')) + 1
names(dhsdf) = c('chr','start','end','name')

# ----------------------------------------
# Make all of the calls for the enhancers:
# ----------------------------------------
# Load the enhancer by epigenomes matrix: 
enames = scan('Enhancer_H3K27ac_intersect_matrix.names.tsv','c')
hd.file = 'Enhancer_H3K27ac_intersect_matrix.hdf5'
h5ls(hd.file)
system('mkdir -p enhancers_bysample/')

chunksize = 25
nchunk = (length(enames) %/% chunksize) + 1
for (i in 1:nchunk){
    cat(paste('Chunk',i,'of',nchunk, '\n'))
    ind = (chunksize * (i-1) + 1):min(c(i *chunksize, length(enames)))
    # Slice the hdf5 matrix:
    h5f = H5Fopen(hd.file)
    h5d = h5f&"matrix"
    mat = h5d[ind,]
    H5Dclose(h5d)
    H5Fclose(h5f)
    rownames(mat) = enames[ind]
    mat = t(mat)
    print(colSums(mat))
    # Write each out:
    for (j in 1:ncol(mat)){
        cat('.')
        id = colnames(mat)[j]
        mn = as.character(mmap[id,'mn'])
        # Get loc:
        hits = which(mat[,id] != 0)
        # Intersect with enhancer ind:
        hits = hits[hits %in% eind]
        df = dhsdf[hits,]
        # Write out:
        write.table(df, gzfile(paste0('enhancers_bysample/', id, '_', mn, '_hg19_enhancer_list.bed.gz')), quote=F, sep="\t", row.names=F)
    }
    cat ('\n')
}

# ----------------------------------------
# Make all of the calls for the promoters:
# ----------------------------------------
# Load the promoter by epigenomes matrix: 
pnames = scan('Promoter_H3K27ac_intersect_matrix.names.tsv','c')
prom.hd.file = 'Promoter_H3K27ac_intersect_matrix.hdf5'
h5ls(prom.hd.file)
system('mkdir -p promoters_bysample/')

chunksize = 25
nchunk = (length(pnames) %/% chunksize) + 1
for (i in 1:nchunk){
    cat(paste('Chunk',i,'of',nchunk, '\n'))
    ind = (chunksize * (i-1) + 1):min(c(i *chunksize, length(pnames)))
    # Slice the hdf5 file:
    h5f = H5Fopen(prom.hd.file)
    h5d = h5f&"matrix"
    mat = h5d[ind,]
    H5Dclose(h5d)
    H5Fclose(h5f)
    rownames(mat) = pnames[ind]
    mat = t(mat)
    print(colSums(mat))
    # Write each out:
    for (j in 1:ncol(mat)){
        cat('.')
        id = colnames(mat)[j]
        mn = as.character(mmap[id,'mn'])
        # Get loc:
        hits = which(mat[,id] != 0)
        # Intersect with promoter ind:
        hits = hits[hits %in% pind]
        df = dhsdf[hits,]
        # Write out:
        write.table(df, gzfile(paste0('promoters_bysample/', id, '_', mn, '_hg19_promoter_list.bed.gz')), quote=F, sep="\t", row.names=F)
    }
    cat ('\n')
}




