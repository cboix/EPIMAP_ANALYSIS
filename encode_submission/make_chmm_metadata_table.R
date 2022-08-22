#!/usr/bin/R
# ----------------------------------------
# Make the ChromHMM metadata table
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

# -------------------------
# Reduce to 833 epigenomes:
# -------------------------
rownames(meta) = meta$id
df = meta[sort(cellorder),]
df$aliases = paste0('manolis-kellis:chromHMM_EpiMap_18-state_', df$id, '_', df$ct)
df$award = 'U24HG009446'
df$lab = 'manolis-kellis'
df$assay_term_name = 'ChIP-seq'
df$description = paste0('ChromHMM 18-state model of ', df$id, 
                        ': ',df$name,' from donor(s) ',df$Donor)
df$encyclopedia_version = 'ENCODE v3'
df$software_used = "encode:ChromHMM 1.12"
df$annotation_type = 'chromatin state'
df$organism = 'human'
df$file_format_type = 'bed9'

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

# Merge at sample level
impdf = aggregate(Biosample.accession ~ id, tab, collate.unique)
names(impdf) = c('id','biosample_accession')
df = merge(df, impdf, all.x=TRUE)

# all 
# all.bioacc = c(unlist(sapply(df$biosample_accession, function(x){strsplit(x,",")[[1]]})))
all.bioacc = c(unlist(sapply(impdf$biosample_accession, function(x){strsplit(x,",")[[1]]})))
# names(all.bioacc) = NULL
# all.bioacc = unique(all.bioacc)

# Get ontology:
library(jsonlite)
obj <- fromJSON('Annotation/biosample_type_ontologylist.json')
bdf = obj$'@graph' # Ontology table (954 obj)
bdf$biosample_ontology = sub('/','',sub('/biosample-types/','',bdf$'@id'))
# Map to biosample??

obj = fromJSON('Annotation/biosample_all_121720.json')
bmapdf = obj$'@graph' # Ontology table (954 obj)
rm(obj); gc()
bmapdf$term_name = bmapdf$biosample_ontology$term_name
bmapdf$classification = bmapdf$biosample_ontology$classification
# bmapdf = bmapdf[bmapdf$accession %in% all.bioacc,]
bmdf = bmapdf[bmapdf$accession %in% all.bioacc,c('accession','term_name','classification')]
bmdf = merge(bmdf, data.frame(accession=all.bioacc), all.y=TRUE)
bmdf = merge(bmdf, bdf[,c('biosample_ontology','term_name','classification')], all.x=TRUE)

# Duplicates:
# cdf = aggregate(term_name ~ accession, bmdf, length)
# dup.bioacc = cdf$accession[cdf$term_name > 1]
# mis.bioacc = all.bioacc[!(all.bioacc %in% bmapdf$accession)]
bmdf$biosample_ontology[bmdf$biosample_ontology == 'primary_cell_NTR_0000498'] = 'primary_cell_CL_0000897'
bmdf$biosample_ontology[bmdf$biosample_ontology == 'primary_cell_NTR_0000503'] = 'primary_cell_CL_0000909'

# Populate missed:
bmdf$biosample_ontology[bmdf$accession %in% c('ENCBS921TFG','ENCBS955MUX')] = 'primary_cell_CL_1001608'
bmdf$biosample_ontology[bmdf$accession %in% c('ENCBS324GKA', 'ENCBS343AKO', 'ENCBS445SNK', 'ENCBS806KLM')] = 'in_vitro_differentiated_cells_CL_0000746'
bmdf = unique(bmdf)
rownames(bmdf) = NULL

# Map the biosample_ontology to samples
ibdf = ldply(1:nrow(impdf), function(i){
                                   samp = impdf$id[i]
                                   x = impdf$biosample_accession[i]
                                   x = strsplit(x,",")[[1]]
                                   acc = bmdf$biosample_ontology[bmdf$accession %in% x]
                                   return(data.frame(id=samp, biosample_ontology=unique(acc)))
                        })

# Resolve bad maps, fixing multiples:
df$biosample_ontology = NULL
cmdiff = ibdf$id[ibdf$biosample_ontology == 'in_vitro_differentiated_cells_CL_0002664']
ibdf$biosample_ontology[ibdf$id == cmdiff] = 'in_vitro_differentiated_cells_CL_0002664'
cmdiff = ibdf$id[ibdf$biosample_ontology == 'primary_cell_CL_0011022']
ibdf$biosample_ontology[ibdf$id == cmdiff] = 'primary_cell_CL_0011022'
ibdf$biosample_ontology[ibdf$id == 'BSS00084'] = 'tissue_UBERON_8300001'
ibdf$biosample_ontology[ibdf$id == 'BSS01154'] = 'tissue_UBERON_8300003'
ibdf$biosample_ontology[ibdf$id == 'BSS01251'] = 'cell_line_EFO_0005696'
ibdf = unique(ibdf)
rownames(ibdf) = NULL
df = merge(df, ibdf, all.x=TRUE)
# df[df$id == 'BSS00084',] # NOTE was relabeled
# df[df$id == 'BSS01154',] # ALSO relabeled

# Write out the biosample mapping
write.table(ibdf, 'Annotation/biosample_ontology_bssid_mapping.tsv', quote=F, sep="\t", row.names=F)

# Get derivation:
tab = read.delim('ChromImpute/mixobs_table_fordist.tsv', header=F)
names(tab) = c('id','mark','set')
tab = tab[tab$id %in% df$id,]
tab = tab[tab$mark %in% c('H3K27ac','H3K4me1','H3K4me3','H3K9me3','H3K36me3', 'H3K27me3'),]
tab$set = 'observed'
tab$set[grep('impute', tab$prefix)] = 'imputed'
tab$derived_from = paste0('manolis-kellis:chromImpute_EpiMap_', tab$id, '_', tab$mark, '_', tab$set, '_bigwig')
tab = aggregate(derived_from ~ id, tab, function(x){paste(x, collapse=',')})
df = merge(df, tab, all.x=TRUE)

# Fix aliases with commas:
df$aliases = gsub("\\,","_", df$aliases)
# df$aliases[grep(",", df$aliases)]

# -------------------------------------------------
# Make the table for submitting Annotation objects:
# -------------------------------------------------
acols = c('aliases',
          'award','lab',
          'assay_term_name',
          'annotation_type',
          'organism',
          'biosample_ontology',
          'encyclopedia_version',
          'software_used',
          'description')
anndf = df[,acols]


# -------------------------------------------
# Make the table for submitting File objects:
# -------------------------------------------
fcols = c('aliases',
          'filealiases',
          # 'final_derived',
          'award','lab',
          'file_format',
          'output_type',
          'submitted_file_name',
          'file_format_type',
          'assembly')

# Files:
filedf = c()
# df$final_derived = df$derived_from
derivdf = c()
for (ext in c('bed','bigBed')){
    df$file_format = ext
    if (ext == 'bed'){
        df$submitted_file_name = paste0(df$id, '_18_CALLS_segments.bed.gz')
    } else {
        df$submitted_file_name = paste0(df$id, '_18_CALLS_segments.bb')
    }
    for (genome in c('hg19','hg38')){
        # chmmdir = paste0('https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_', genome, '/CALLS/bigBed/')
        # For files:
        df$output_type = 'semi-automated genome annotation'
        if (genome == 'hg38'){
            df$assembly = 'GRCh38'
        } else {
            df$assembly = genome
        }
        df$filealiases = paste0(df$aliases, '_', genome, '_', ext)
        # if (ext == 'bigBed'){
        #     df$final_derived = paste0(df$derived_from, ',', df$aliases, '_', genome, '_bed')
        # } else {
        #     if (genome == 'hg38'){
        #         df$final_derived = paste0(df$derived_from, ',', df$aliases, '_hg19_bed')
        #     }
        # }
        derivdf = rbind(derivdf, df[,c(fcols, 'derived_from')])
        filedf = rbind(filedf, df[,fcols])
    }
}
names(filedf)[1:2] = c('dataset','aliases') # , 'derived_from')
# NOTE: NOT USING DERIVED FROM BECAUSE WE ARE NOT UPLOADING OBSERVED DATA.

# Note: need to be in same directory as file to submit.
write.table(anndf, 'Annotation/encode_chmm_annotation_submission_metatable.tsv', quote=F, sep="\t", row.names=F)
write.table(filedf, 'Annotation/encode_chmm_file_submission_metatable.tsv', quote=F, sep="\t", row.names=F)

write.table(derivdf, 'Annotation/encode_chmm_file_submission_metatable_with_deriv.tsv', quote=F, sep="\t", row.names=F)
write.table(df[,c(acols, 'age','age.units')], 'Annotation/encode_chmm_annotation_submission_metatable_with_age.tsv', quote=F, sep="\t", row.names=F)



