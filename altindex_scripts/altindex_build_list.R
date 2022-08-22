#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Method for generating a master list / Index of DNaseI hypersensitivity sites.
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.
# Adapted for use with ENCODE Imputation by Carles Boix
# ------------------------------------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'altindex_auxfunctions.R'))
library(caTools)
options(scipen=20)
options(warn=2)

args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied.")
} else {
    workdir = args[[1]]
}
chunknum = as.numeric(Sys.getenv("SGE_TASK_ID"))

# Go to directory:
setwd(workdir)

dir.create("DHSs_all", showWarnings=FALSE, recursive=TRUE)
dir.create("peaks_all", showWarnings=FALSE, recursive=TRUE)

chunk <- paste("chunk", sprintf("%05d", chunknum), ".bed", sep="");
print(paste("Loading in data from chunk file:", chunk))

# Load in data per chunk, each separated by at least 10kb
peaks <- read.delim(paste(workdir,chunk,sep='/'), header=FALSE, as.is=T)
colnames(peaks) <- c("seqname", "start", "end", "sampleID", "max", "wavelet_summit")
peaks$score = 10^(-peaks$max)
peaks <- peaks[order(peaks$wavelet_summit),] # Order peaks by wavelet_summit first

DHSs <- NULL; # Final list of delineated DHSs

# Overwrite peak IDs, chunked based on a >25bp separation of wavelet_summits
# NOTE: Changed this value due to 25bp resolution after imputation.
loc_ID_seps <- c(0, which(diff(peaks$wavelet_summit) > 25), nrow(peaks))
loc_IDs <- rep(1:(length(loc_ID_seps)-1), times=diff(loc_ID_seps))
peaks$ID <- loc_IDs;

# Iterate over peak clumps
for (loc_ID in unique(peaks$ID)) {
    # Localize data, in particular wavelet_summit coordinates
    sel <- which(peaks$ID == loc_ID);
    loc_peaks <- peaks[sel,]
    xlim <- range(loc_peaks$wavelet_summit)
    loc_peaks$wavelet_summit_loc <- loc_peaks$wavelet_summit - xlim[1] + 1
    n <- diff(xlim)+1

    # Obtain cut points, to be processed independently
    cuts <- get_cut_points(loc_peaks$wavelet_summit_loc);

    # Select localized data for FWHM determination 
    loc_DHSs <- NULL
    for (i in 1:(length(cuts)-1)) {
        idx <- which(loc_peaks$wavelet_summit_loc > cuts[i] & loc_peaks$wavelet_summit_loc <= cuts[i+1])
        if (length(idx) == 0) next;
        new_DHS <- get_FWHM(loc_peaks[idx,])
        new_DHS$ID <- paste(loc_ID, i, sep="_")
        loc_peaks$ID[idx] <- paste(loc_ID, i, sep="_")
        loc_DHSs <- rbind(loc_DHSs, new_DHS) # Save resultant (candidate) DHS
    }

    if (nrow(loc_DHSs) > 1) {
        # Resolve summit overlaps (summit of one element in FWHM of another)
        num <- 1;
        while(num > 0) {
            ovl_rm <- merge_overlap(loc_DHSs, loc_peaks, type="summit")
            loc_DHSs <- ovl_rm$DHSs
            loc_peaks <- ovl_rm$peaks
            num <- ovl_rm$num
        }
    }

    DHSs <- rbind(DHSs, loc_DHSs)
    peaks$ID[sel] <- loc_peaks$ID
}

# Final pass across all data in this chunk
if (nrow(DHSs) > 1) {
    # Resolve summit overlaps across originally defined "peak clumps"
    num <- 1;
    while(num > 0) {
        ovl_rm <- merge_overlap(DHSs, peaks, type="summit")
        DHSs <- ovl_rm$DHSs
        peaks <- ovl_rm$peaks
        num <- ovl_rm$num
        print(num)
    }
}
DHSs <- DHSs[order(DHSs$seqname, DHSs$start, DHSs$end),]

peaks$ID <- paste(gsub(".bed", "", chunk), peaks$ID, sep="_")
DHSs$ID <- paste(gsub(".bed", "", chunk), DHSs$ID, sep="_")

write.table(DHSs, file=paste("DHSs_all", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(peaks, file=paste("peaks_all", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

cat("[STATUS] Finished building chunk")
