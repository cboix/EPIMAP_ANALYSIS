#!/usr/bin/R
# ----------------------------------------
# Make the bigwig metadata table
# 1) Annotation objects
# 2) File objects
# Same size, matching rows so we can slice
# ----------------------------------------
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

# -----------------------------
# All bigwig tracks to release:
# -----------------------------
df = read.delim('all_released_tracks.tsv', header=F)
names(df) = c('id','mark','prefix')
df$set = 'observed'
df$set[grep('impute', df$prefix)] = 'imputed'
# Metadata + aliases:
df = merge(df, meta, all.x=TRUE)
df$aliases = paste0('manolis-kellis:chromImpute_EpiMap_', df$id, '_', df$mark, '_', df$set)
df$award = 'U24HG009446'
df$lab = 'manolis-kellis'
df$assay_term_name = 'ChIP-seq'
df$assay_term_name[df$mark == 'DNase-seq'] = 'DNase-seq'
df$assay_term_name[df$mark == 'ATAC-seq'] = 'ATAC-seq'
df$targets = paste0(df$mark,'-human')
df$targets[df$mark == 'DNase-seq'] = NA
df$targets[df$mark == 'ATAC-seq'] = NA
df$annotation_type = 'imputation'
df$software_used = 'manolis-kellis:chrom_impute_1_03'
df$organism = 'human'

# Read biosample derivations:
tab = read.delim('all_encode_file_imputation_identifiers_mappings.tsv', header=F, sep="\t")
names(tab) = c('mark','url','accid','genome','align','id','other','ct')
tab$Accession = sub("\\..*","",tab$accid)
bioaccmap = read.delim('Annotation/all_submitted_released_biosample_accession_mapping_20180924.tsv', header=T, sep="\t")
tab = merge(tab, bioaccmap, all.x=TRUE)
tab$file = sub("/.*","", sub(".*/files/","",tab$url))

# File derivations:
collate.unique = function(x){
    x = c(unlist(sapply(x, function(x){strsplit(x,",")[[1]]})))
    x = sort(unique(x))
    paste(x, collapse=",")
}

# Ontology:
ibdf = read.delim('Annotation/biosample_ontology_bssid_mapping.tsv', header=T)

# At the file level:
repfiledf = aggregate(file ~ id + mark, tab, collate.unique)
names(repfiledf) = c('id','mark', 'derived_from')

# At sample level
impfiledf = aggregate(file ~ id, tab, collate.unique)
names(impfiledf) = c('id','derived_from')

# Merge observed at file level and imputed at sample level:
idf = df[df$set == 'imputed',]
odf = df[df$set == 'observed',]
idf = merge(idf, impfiledf, all.x=TRUE)
odf = merge(odf, repfiledf, all.x=TRUE)
df = rbind(idf, odf)
df = merge(df, ibdf, all.x=TRUE)

# Description
df$description = ''
df$description[df$set == 'observed'] = paste0('observed EpiMap hg19 track of ', df$mark,' in ', df$id, ': ',df$name, ' from donor(s) ',df$Donor)[df$set == 'observed']
df$description[df$set == 'imputed'] = paste0('chromImpute EpiMap hg19 imputation of ', df$mark,' in ', df$id, ': ',df$name,' from donor(s) ',df$Donor)[df$set == 'imputed']
df$encyclopedia_version = 'ENCODE v3'

# For files:
df$file_format = 'bigWig'
df$output_type = 'signal p-value'
df$submitted_file_name = paste0(df$prefix, '.bigWig')
df$assembly = 'hg19'
df$filealiases = paste0(df$aliases, '_bigwig')

# -------------------------------------------------
# Make the table for submitting Annotation objects:
# -------------------------------------------------
acols = c('aliases',
          'award','lab',
          'annotation_type',
          'assay_term_name',
          'biosample_ontology',
          'encyclopedia_version',
          'description', 
          'software_used',
          'organism',
          'targets')
anndf = df[df$set == 'imputed',acols]

# -------------------------------------------
# Make the table for submitting File objects:
# -------------------------------------------
fcols = c('aliases',
          'filealiases',
          'award','lab',
          'file_format',
          'output_type',
          'submitted_file_name',
          'derived_from',
          'assembly')
filedf = df[df$set == 'imputed',fcols]
names(filedf)[1:2] = c('dataset','aliases')

# Note: need to be in same directory as file to submit.
write.table(anndf, 'Annotation/encode_bigwig_annotation_submission_metatable.tsv', quote=F, sep="\t", row.names=F)
write.table(filedf, 'Annotation/encode_bigwig_file_submission_metatable.tsv', quote=F, sep="\t", row.names=F)

