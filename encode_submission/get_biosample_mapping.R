#!/usr/bin/R
# ----------------------------------
# Get ENCODE file metadata:
# Run after get_experiments_ENCODE.R
# ----------------------------------
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
library(jsonlite)
library(httr)
today = format(Sys.time(), "%Y%m%d")
today = '20180924' # FROZEN metadata get

args=(commandArgs(TRUE))
if (length(args)==0) {
    genus = 'Homo'
    species = 'sapiens'
    print("No arguments supplied: Using Homo sapiens.")
} else {        
    print(args)
    genus = args[1]
    species = args[2]
}

# ----------------
# Naming prefixes:
# ----------------
nam = 'all_submitted_released' 
if (genus == 'Homo'){ 
    prefix = ''
} else { 
    prefix = tolower(paste0(substr(genus,0,1),
                            substr(species,0,1)))
    nam = paste(genus,species,nam,sep='_')
}


if (prefix == 'mm') { assemblies = c('mm10','mm9','mm10-minimal') } 
if (prefix == '') { assemblies = c('hg19', 'GRCh38') } 
outputs = c('alignments', 'redacted alignments', 
            'unfiltered alignments', 'redacted unfiltered alignments')


# =============
# Get metadata:
# =============
outfile <- paste0('Annotation/',nam,'_metadata.tsv')
metadata <- read.delim(outfile,sep="\t",header=T)
metadata$UID <- paste(metadata$Experiment.accession,
                      metadata$Biological.replicate.s.,
                      metadata$Technical.replicate,
                      sep=".")






# Get the biosamples mapped to the original experiments used

# Map one biosample to each BSSID




