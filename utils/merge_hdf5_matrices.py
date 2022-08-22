# !/usr/bin/python
# -----------------------------------------------------------------
# Merge all datasets, so that we can perform PCA, etc. on all data:
# -----------------------------------------------------------------
import h5py
import numpy as np

# Relative to EPIMAP_ANALYSIS/db ($DBDIR)
marks = ['H3K27ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac','DNase-seq']
lddir = 'linking_data/'
datasuf = '_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5'

merge_hdf5 = lddir + 'merged6mark' + datasuf
with h5py.File(merge_hdf5, mode='w') as h5fw:
    col1 = 0
    # for h5name in glob.glob('file*.h5'):
    for i, mark in enumerate(marks):
        print(i,mark)
        data_hdf5 = lddir + mark + datasuf
        h5fr = h5py.File(data_hdf5,'r')
        arr_data = h5fr['matrix'][:]
        nrows = arr_data.shape[0]
        ncols = arr_data.shape[1]
        if col1 == 0:
            h5fw.create_dataset('alldata', dtype="f",
                                shape=(nrows,ncols*len(marks)),
                                maxshape=(nrows, ncols*len(marks)))
        if col1+ncols <= h5fw['alldata'].shape[1] :
            h5fw['alldata'][:, col1:col1+ncols] = arr_data[:]
        else :
            # Resize array to fit:
            h5fw['alldata'].resize((nrows, col1+ncols) )
            h5fw['alldata'][:, col1:col1+ncols] = arr_data[:]
        print(h5fw['alldata'].shape)
        h5fr.close()
        col1 += ncols


