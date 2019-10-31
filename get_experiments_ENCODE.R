#!/usr/bin/R
# --------------------------------------------------
# Obtain new biosample mapping and file annotations:
# Follow this by running get_metadata_ENCDODE.R
# Which will get the files corresponding to each project.
# --------------------------------------------------
source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
library(jsonlite)
library(httr)
today = format(Sys.time(), "%Y%m%d")

################################################################################################################
consumer_key = "R7GOMSZL";
consumer_secret = "7sccfx4rtdb6sqse";
secret <- RCurl::base64(paste(consumer_key, consumer_secret, sep=":"));
################################################################################################################
# IF DOESNT WORK, GET EXPERIMENTS AS:
# curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o ~/EPIMAP_ANALYSIS/db/Annotation/all_submitted_released.json https://www.encodeproject.org/search/\?type\=Experiment\&status\=released\&status\=submitted\&replicates.library.biosample.donor.organism.scientific_name\=Homo+sapiens\&format\=json\&limit\=all
################################################################################################################

# ### This is where we obtain most of the metadata
nam = 'all_submitted_released' 
inner = 'type=Experiment&status=released&status=submitted&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens'
url <- paste("https://www.encodeproject.org/search/?", inner, "&format=json&limit=all", sep="");
# cmd = paste0("curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o Annotation/", nam, ".json ", url)
# source(cmd)

annfile = paste0('Annotation/all_experiments_', today, '.tsv')
if (!file.exists(annfile)){
    # # Get and save JSONS for tagAlign and Bam:
    for (filetype in c('bam', 'tagAlign')){
        url_files <- paste0("https://www.encodeproject.org/search/?type=file&status=released&status=submitted&file_format=", filetype, "&frame=object&format=json&limit=all", sep="")
        cmd = paste0("curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o Annotation/all_files_", filetype, ".json ", url_files)
        print(cmd)
        # system(cmd)
    }
    # TSV report (works best):
    url = 'https://www.encodeproject.org/report.tsv?type=Experiment'
    cmd = paste0("curl -L -u 'R7GOMSZL:7sccfx4rtdb6sqse' -o Annotation/all_experiments_", today, ".tsv ", url)
    print(cmd)
    system(cmd)
}

exp.df = read.delim(annfile, skip=1, header=T, sep="\t")
exp.df = exp.df[exp.df$Species == 'Homo sapiens',]
exp.df = exp.df[exp.df$Status %in% c('released','submitted'),]
print(dim(exp.df))

cols = c('Accession','Lab','Project','Biosample.accession','Life.stage','Age','Age.Units','Replicates', 'Treatment')
otherdf = exp.df[,cols]

# ============================================================
# Make the info files: prioritized marks + list all cell types
# Load in original mapping:
# ============================================================
mapfile = paste0('Annotation/',nam,'_biosample_mapping.tsv')
newmapfile = paste0('Annotation/',nam,'_biosample_mapping_', today,'.tsv')
maptable <- read.delim(mapfile, header=F, sep="\t")
names(maptable) <- c('Biosample.summary', 'uid')

# --------------------------
# Fix the "treated" naming: 
# --------------------------
treat.id <- grep('treat', exp.df$Description)
fix.id = grep('treat', exp.df$Biosample.summary[treat.id], invert=TRUE)
# Remove other treatments:
fix.id = fix.id[grep(' 0 ', exp.df$Description[treat.id[fix.id]], invert=TRUE)]
fix.id = fix.id[grep('control', exp.df$Description[treat.id[fix.id]], invert=TRUE)]
fix.id = fix.id[grep('phase', exp.df$Description[treat.id[fix.id]], invert=TRUE)]

fix.desc = gsub(".* on [human ]*", "", exp.df$Description[treat.id[fix.id]])
fix.desc = gsub("cell line ", "", fix.desc)
exp.df$Biosample.summary = as.character(exp.df$Biosample.summary)
exp.df$Biosample.summary[treat.id[fix.id]] = gsub("mins\\.", "minutes", fix.desc)

print("Fixed Biosample summary names:")
print(exp.df$Biosample.summary[treat.id])

# Drop anything not submitted/released:
sums = sort(unique(exp.df$Biosample.summary))
map <- data.frame(Biosample.summary=sums)  # Helps remove bad characters

# Contains old and new replacements?
newmap <- merge(map, maptable, all.x=TRUE)
newnames <- newmap[is.na(newmap$uid), 'Biosample.summary']
oldsums <- maptable$Biosample.summary
oldnames <- oldsums[!(oldsums %in% newmap$Biosample.summary)]

print(which(is.na(newmap$uid)))
# NOTE: Will just update new map as of 09/18/18.

UIDS = sprintf("BSS%05d", 1:length(sums))
map <- data.frame(Biosample.summary=sums, uid=UIDS)  # Helps remove bad characters
exp.df <- merge(exp.df, map)

acc.id = which(exp.df$Assay.Type %in% c('ChIP-seq', 'ATAC-seq', 'DNase-seq'))
acc.df = exp.df[acc.id,]
control.id = which(acc.df$Target.label == "Control")
hist.id = (1:nrow(acc.df) %in% grep("^H[0-4]", acc.df$Target.label))
acc.df$is.hist = hist.id
chip.id = which(acc.df$Assay.Type == 'ChIP-seq')
open.id = which(acc.df$Assay.Type %in% c('ATAC-seq', 'DNase-seq'))
all.id = sort(unique(c(chip.id, control.id, open.id)))
# all.id = sort(unique(c(hist.id, control.id, open.id)))
acc.df = acc.df[all.id,]

# Count biosamples and cells:
cts = sort(unique(acc.df$Biosample))
bss = sort(unique(acc.df$Biosample.summary))
print(paste("CT #:", length(cts)))
print(paste("BSS #:", length(bss)))

tab = acc.df[grep('treat', acc.df$Description),
             c('Biosample.summary','Description')] 

# Marks; for prioritization:
sel <- which(acc.df$Assay.Type == 'ChIP-seq')
df <- data.frame(count=sel, marks=acc.df$Target.label[sel])
df <- aggregate(count ~ marks,df,length)
df <- df[order(df$count,decreasing=TRUE),]

write.table(map, file=newmapfile,row.names=F,col.names=F,quote=F,sep="\t")
write.table(df, file=paste0('Annotation/',nam,'_ChIP_abundance_', today, '.tsv'),row.names=F,col.names=F,quote=F,sep="\t")

ct <- sort(unique(acc.df$Biosample.summary))
write.table(ct, file=paste0('Annotation/',nam,'_cell_types_', today, '.tsv'),row.names=F,col.names=F,quote=F,sep="\t")

save(acc.df, file=paste("Annotation/experiments_", nam, "_", today ,".RData", sep=""))


# NOTE: OLD CODE. REMOVE WHEN SURE WORKS ABOVE:
# if (length(grep(paste("^Annotation/experiments_", nam, ".RData$", sep=""), dir())) == 0) {

#     call <- GET(url, config(c("Authorization" = paste("Basic", secret))))
#     obj <- fromJSON(rawToChar(call$content))
#     status <- obj$facets[6,"terms"][[1]]

#     if (status[status$key=='submitted','doc_count'] == 0){ 
#         print('Couldn\'t access submitted data, loading previously stored JSON file') 
#         obj <- fromJSON(paste0('Annotation/',nam,'.json'))
#     }

#     # Parse some of the data and save.
#     terms <- obj$facets[[2]]
#     assays <- obj$facets[3,"terms"][[1]]
#     dates <- obj$facets[22,"terms"][[1]]
#     filetypes <- obj$facets[15,"terms"][[1]];
#     experiments <- obj[["@graph"]]

#     # ============================================================
#     # Make the info files: prioritized marks + list all cell types
#     # ============================================================
#     # LOAD IN ORIGINAL MAPPING:
#     mapfile = paste0('Annotation/',nam,'_biosample_mapping.tsv')
#     maptable <- read.delim(mapfile, header=F, sep="\t")
#     names(maptable) <- c('biosample_summary', 'uid')
#     # Drop anything not submitted/released:
#     experiments <- experiments[experiments$status %in% c('submitted','released'),]
#     sums = sort(unique(experiments$biosample_summary))
#     map <- data.frame(biosample_summary=sums)  # Helps remove bad characters
#     # Contains old and new replacements?
#     newmap <- join(map, maptable)
#     newnames <- newmap[is.na(newmap$uid), 'biosample_summary']
#     oldsums <- maptable$biosample_summary
#     oldnames <- oldsums[!(oldsums %in% newmap$biosample_summary)]
#     print(which(is.na(newmap$uid)))

#     UIDS = sprintf("BSS%05d", 1:length(sums))
#     map <- data.frame(biosample_summary=sums, uid=UIDS)  # Helps remove bad characters
#     experiments <- join(experiments, map)

#     # Disambiguate treated/normal!
#     desc <- experiments$description
#     treat <- experiments[grep('treat',desc),] 
#     experiments <- experiments[-grep('treat',desc),] 

#     # Marks; for prioritization:
#     sel <- which(experiments$assay_term_name == 'ChIP-seq')
#     df <- data.frame(count=sel, marks=experiments$target[sel,1])
#     df <- aggregate(count ~ marks,df,length)
#     df <- df[order(df$count,decreasing=TRUE),]

#     write.table(map,file=mapfile,row.names=F,col.names=F,quote=F,sep="\t")
#     write.table(df,file=paste0('Annotation/',nam,'_ChIP_abundance.tsv'),row.names=F,col.names=F,quote=F,sep="\t")

#     ct <- sort(unique(experiments$biosample_term_name))
#     write.table(ct,file=paste0('Annotation/',nam,'_cell_types.tsv'),row.names=F,col.names=F,quote=F,sep="\t")

#     save(assays, dates, filetypes, experiments, terms, treat, file=paste("Annotation/experiments_", nam, ".RData", sep=""))
# }

