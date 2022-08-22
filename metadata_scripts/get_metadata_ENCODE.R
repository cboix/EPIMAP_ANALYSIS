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

# ================
# Naming prefixes:
# ================
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

################################################################################################################
consumer_key = "R7GOMSZL";
consumer_secret = "roctnxvjxljlfrwg";
secret <- RCurl::base64(paste(consumer_key, consumer_secret, sep=":"));
################################################################################################################

# =============
# Get metadata:
# =============
inner = paste0('type=Experiment&status=released&status=submitted&replicates.library.biosample.donor.organism.scientific_name=',genus,'+',species)
url <- paste("https://www.encodeproject.org/metadata/", inner, "/metadata.tsv", sep="");
# NOTE: Is correct, may need to wait to download - very large table
outfile <- paste0('Annotation/',nam,'_metadata.tsv')
if (!file.exists(outfile)){ 
    auth <- paste(consumer_key, consumer_secret, sep=":")
    cmd <- paste0("curl -L -u '",auth,"' -o ",outfile," ",url)
    print(cmd)
    system(cmd)
}

metadata <- read.delim(outfile,sep="\t",header=T)
metadata$UID <- paste(metadata$Experiment.accession,
                      metadata$Biological.replicate.s.,
                      metadata$Technical.replicate,
                      sep=".")

# THE FOLLOWING TABLE COMES FROM:
# https://www.encodeproject.org/report/?type=Experiment&field=%40id&field=accession&field=assay_title&field=target.label&field=biosample_summary&field=biosample_term_name&field=award.project&field=status&field=replicates.library.biosample.accession&field=replicates.library.biosample.life_stage&field=replicates.library.biosample.age&field=replicates.library.biosample.age_units&field=replicates.library.biosample.treatments.treatment_term_name&field=replicates.library.biosample.treatments.concentration&field=replicates.library.biosample.treatments.concentration_units&field=replicates.library.biosample.treatments.duration&field=replicates.library.biosample.treatments.duration_units&field=replicates.library.biosample.synchronization&field=biosample_type&field=date_released&field=lab&field=month_released&field=replication_type&field=replicates.library.biosample.donor.accession&field=replicates.library.biosample.sex&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&sort=-date_released
extdf = read.delim('Annotation/experiment_metadata_20180128.tsv',header=T, sep="\t", skip=1)
iddf = metadata[,c('File.accession','File.download.URL','Experiment.accession')]
names(iddf) <- c('File','URL','Accession')
# GET RNA:
rnadf = extdf[extdf$Assay.Nickname %in% c('polyA RNA-seq', 'RNA microarray', 'total RNA-seq') ,]
rnadf = merge(iddf, rnadf)
# Get the file - experiment + metadata mapping table:
extdf = extdf[extdf$Assay.Nickname %in% c('ChIP-seq','DNase-seq','ATAC-seq'),]
iddf = merge(iddf, extdf)
write.table(iddf, 'Annotation/file_extended_metadata_20180128.tsv', row.names=F, col.names=T, sep="\t", quote=F)

# Write the award/project into separate file:
awarddf = metadata[,c('File.accession','File.download.URL','Project')]
write.table(awarddf, 'file_project_origin.tsv',sep="\t", row.names=F, col.names=T, quote=F)

get.unique.perUID <- function(df){ 
    uids <- unique(df$UID)
    eids <- c()
    for (uid in uids){ 
        ids <- which(df$UID == uid)
        if (length(ids) == 1){
            eids <- c(eids, ids)
        } else {
            # Take only the first eid:
            eids <- c(eids, ids[1])
        }
    }
    return(df[eids,])
}

# Load experiments list (from get_experiments_ENCODE.R)
load(file=paste("Annotation/experiments_", nam, "_", today, ".RData", sep=""))
mapping = acc.df[, c('ID','uid', 'Biosample.summary')]
print(paste("Number of experiments:", dim(mapping)[1]))

# PREPROCESSING OF DATASETS:
# 1. Remove all genetically modified:
modif.idx = grep("modif", mapping$Biosample.summary)
mapping.nomod = mapping[-modif.idx,]

# ================================
# Get assay-specific data:
# TODO also add RNA-seq/polyRNA/etc?
# ================================
assays <- c('ChIP-seq','DNase-seq','ATAC-seq')
for (i in 1:length(assays)){ 
    assay = assays[i]
    print(paste("Getting metadata and links for",assay))

    # Output directory:
    links_dir = paste(assay,"file_links", nam, sep="/");
    dir.create(links_dir, showWarnings=FALSE, recursive=TRUE);

    # Filter for relevant datasets:
    df <- metadata[metadata$Assay == assay,]
    df <- df[df$File.format %in% c('bam','tagAlign'),]
    df <- df[df$File.Status %in% c('released','submitted'),]
    data <- df[,c('File.accession','Output.type','UID', 'Biosample.term.name',
                  'Biosample.treatments','Experiment.target','Lab','Mapped.read.length',
                  'File.download.URL','Assembly','Audit.NOT_COMPLIANT', 'Experiment.accession', 'Experiment.date.released')]
    data$ID = paste0("/experiments/", data$Experiment.accession, "/")
    #data = merge(data, mapping)
    data = merge(data, mapping.nomod) # Remove all data from modified cells
    print(paste("Start # BSSIDs:", length(sort(unique(data$uid)))))

    # Prioritize processed and newer assemblies: 
    print(paste0("Before filtering for outputs and assemblies: ", paste(dim(data), collapse=", ")))
    data = data[data$Output.type %in% outputs,]
    data = data[data$Assembly %in% assemblies,]
    print(paste0("After filtering for outputs and assemblies: ", paste(dim(data), collapse=", ")))
    data$Output.type <- factor(data$Output.type, levels=outputs)
    data$Assembly <- factor(data$Assembly,levels=assemblies)
    # NOTE: THERE ARE 
    data <- data[order(data$Assembly,data$Output.type,data$Audit.NOT_COMPLIANT),]
    if (prefix == 'mm'){ 
        data$epitope = sub('-mouse','',data$Experiment.target)
    } else {
        data$epitope = data$Experiment.target
    }
    # Pick out the top one:
    # TODO: Ensure correct with tagAligns:
    data.uniq <- get.unique.perUID(data)
    print(paste("End # BSSIDs:", length(sort(unique(data.uniq$uid)))))

    # Use uid as cell type and clean text:
    data.uniq$cell_type = data.uniq$uid
    data.uniq$biosamp <- gsub("\'", "_", gsub("\\/", "_", 
                                              gsub("%", "pct", gsub("\ ", "_",data.uniq$Biosample.summary))))
    data.uniq$biosamp <- stringi::stri_trans_general(data.uniq$biosamp, 'greek-latin')
    if (assay == 'ChIP-seq'){ 
        # Rename Control to WCE:
        data.uniq$epitope = gsub("-human","", data.uniq$epitope)
        data.uniq$epitope <- as.character(data.uniq$epitope)
        cid <- which(data.uniq$epitope == 'Control')
        if (length(cid) > 0) data.uniq$epitope[cid] <- 'WCE'
        data.uniq$epitope <- factor(data.uniq$epitope)
        # NOTE: Treatments for ChIP-seq are handled in get experiments script.
        data.uniq$submitted_file_name= 'none'
        out.df <- data.uniq[,c('epitope','File.download.URL','UID','Assembly','Output.type','uid', 'submitted_file_name', 'biosamp', 'Experiment.date.released')]
        epitopes <- sort(unique(out.df$epitope))
        colnames(out.df) <- c('epitope','file','lbl','assembly','output_type','cell_type','submitted_file_name', 'biosamp', 'date')
        # Write each of epitopes out:
        pb = txtProgressBar(min = 0, max = length(epitopes), style = 3)
        for (j in 1:length(epitopes)){
            epitope = epitopes[j]
            edf <- out.df[out.df$epitope == epitope,-1]
            aggdf = aggregate(file ~ cell_type + biosamp, edf, length)
            aggdf = aggdf[order(aggdf$file, decreasing=TRUE),]
            write.table(edf, file=paste0(links_dir, "/", epitope,"_bam.csv"), quote=FALSE, row.names=FALSE, sep="\t");
            setTxtProgressBar(pb, j)
        }
        close(pb)
        # ChIP-seq epitopes metadata:
        data.uniq$count <- 1:nrow(data.uniq)
        mdf <- aggregate(count ~ epitope,data.uniq,length)
        mdf <- mdf[order(mdf$count,decreasing=TRUE),]
        write.table(mdf,file=paste0('Annotation/',nam,'_',assay,'_files_abundance.tsv'),row.names=F,col.names=F,quote=F,sep="\t")
    } else { 
        # Write assay links out for DNase/ATAC/other:
        data.uniq$submitted_file_name= 'none'
        out.df <- data.uniq[,c('File.download.URL','UID','Assembly','Output.type','cell_type', 'submitted_file_name', 'biosamp')]
        colnames(out.df) <- c('file','lbl','assembly','output_type','cell_type','submitted_file_name', 'biosamp')
        write.table(out.df, file=paste0(links_dir, "/", assay, "_bam.csv"), quote=FALSE, row.names=FALSE, sep="\t");
    }
}


# --------------------------------------
# Also get matching RNA-seq/polyRNA/etc:
# --------------------------------------
# Filter for relevant datasets:
assay = 'RNA-seq'
df <- metadata[metadata$Assay == assay,]
print(paste("Getting metadata and links for",assay))

# Output directory:
links_dir = paste(assay,"file_links", nam, sep="/");
dir.create(links_dir, showWarnings=FALSE, recursive=TRUE);

# rnatypes = c('polyA RNA-seq', 'RNA microarray', 'total RNA-seq')
df <- df[df$File.Status %in% c('released','submitted'),]
df <- df[df$File.format %in% c('tsv','gtf','bam'),]


data <- df[,c('File.accession','Output.type','UID', 'Biosample.term.name',
              'Biosample.treatments','Experiment.target','Lab','Mapped.read.length',
              'File.download.URL','Assembly','Audit.NOT_COMPLIANT', 'Experiment.accession', 'Experiment.date.released')]
data$ID = paste0("/experiments/", data$Experiment.accession, "/")
data = merge(data, rnadf[, c('ID','Biosample.summary')])
data = unique(data)

# Merge:
uniqmap = unique(mapping.nomod[,c('uid','Biosample.summary')])
data = merge(data, uniqmap) # Remove all data from modified cells
print(paste("Start # BSSIDs:", length(sort(unique(data$uid)))))

# Pick preferentially processed?
unique(data$Output.type)
rnaouts= c('gene quantifications', 'transcript quantifications')
# 'transcriptome alignments', 'alignments')  # NOTE: Alignments don't add any BSSID

# Prioritize processed and newer assemblies: 
print(paste0("Before filtering for outputs and assemblies: ", paste(dim(data), collapse=", ")))
data = data[data$Output.type %in% rnaouts,]
data = data[data$Assembly %in% assemblies,]
print(paste0("After filtering for outputs and assemblies: ", paste(dim(data), collapse=", ")))
data$Output.type <- factor(data$Output.type, levels=rnaouts)
# NOTE: Reverse assembly choice here:
data$Assembly <- factor(data$Assembly,levels=rev(assemblies))
data <- data[order(data$Output.type,data$Assembly,data$Audit.NOT_COMPLIANT),]

# Pick out the top one:
data.uniq <- get.unique.perUID(data)
print(paste("End # BSSIDs:", length(sort(unique(data.uniq$uid)))))

data.uniq$Output.type

# Use uid as cell type and clean text:
data.uniq$cell_type = data.uniq$uid
data.uniq$biosamp <- gsub("\'", "_", gsub("\\/", "_", 
                                          gsub("%", "pct", gsub("\ ", "_",data.uniq$Biosample.summary))))
data.uniq$biosamp <- stringi::stri_trans_general(data.uniq$biosamp, 'greek-latin')
# Write assay links out for DNase/ATAC/other:
data.uniq$submitted_file_name= 'none'
out.df <- data.uniq[,c('File.download.URL','UID','Assembly','Output.type','cell_type', 'submitted_file_name', 'biosamp')]
colnames(out.df) <- c('file','lbl','assembly','output_type','cell_type','submitted_file_name', 'biosamp')
write.table(out.df, file=paste0(links_dir, "/", assay, "_tsv.csv"), quote=FALSE, row.names=FALSE, sep="\t");

