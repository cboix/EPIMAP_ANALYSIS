#!/usr/bin/R
# -----------------------------------------------------
# Make the bigwigs trackhub for the average tracks:
# Have observed, imputed, and by sample slicing
# Color tracks by either sample type or by obs/imputed?
# -----------------------------------------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
source(paste0(bindir, 'general_EPIMAP_ANALYSIS.R'))
source(paste0(bindir, 'auxiliary_chromImpute_functions.R'))
source(paste0(bindir, 'load_metadata.R'))
library(tidyr)
library(dplyr)
options(scipen=45) # so we dont get issues writing integers into bedfiles
today <- format(Sys.time(), "%m%d%y")


# TODO: Make an average trackhub (core marks, average tracks)
# TODO: Make a per-group trackhub??

# -----------------------
# Read in average tracks:
# -----------------------
avgdir = 'https://personal.broadinstitute.org/cboix/epimap/averagetracks_bygroup/'
avdf = read.delim('average_trackslist.txt', header=F)
names(avdf) = 'file'
avdf$prefix = sub(".bigWig", "", sub("average_","",avdf$file))
smat = t(sapply(avdf$prefix, function(x){ x = sub("Sm_M","SmM",x); x = strsplit(x, split="_")[[1]] }))
avdf$mark = smat[,1]
avdf$dataset = sapply(smat[,2], capitalize)
avdf$groupstr = smat[,3]
avdf$url = paste0(avgdir, avdf$file)

# Merge with color:
odf$groupstr = sub(" & ","and", odf$GROUP)
odf$groupstr = sub("\\. ","", odf$groupstr)
avdf = merge(avdf, odf, all.x=TRUE)
avdf$groupsp = gsub(" ","_", avdf$GROUP)

# Number of supporting tracks:
tdf = read.delim('data_to_release/all_released_tracks_withgroup.tsv', header=F)
names(tdf) = c('id','mark','track','GROUP')
tdf$dataset = 'Observed'
tdf$dataset[grep("impute", tdf$track)] = "Imputed"
totdf = aggregate(id ~ mark + GROUP + dataset, tdf, length)
# TODO: FIX ALL single tracks - failed.
# t2df = merge(totdf, avdf, all.x=TRUE)
# t2df$id[is.na(t2df$groupstr)]
avdf = merge(avdf, totdf, all.x=TRUE)

# TODO: Do DNase-seq and Histone marks have diff. properties
avdf$track = paste0("track BW_",  avdf$groupstr, "_", avdf$mark, "_", substr(avdf$dataset,1,3))
avdf$shortLabel = paste("shortLabel", substr(avdf$dataset,1,3), avdf$GROUP, avdf$mark)
avdf$longLabel = paste0("longLabel Average of ", avdf$id, " ", ifelse(avdf$id == 1, "track", "tracks"),
                        " for ", avdf$GROUP, " ", avdf$mark, ' (', avdf$dataset, ')')
avdf$bigDataUrl = paste0('bigDataUrl ', avdf$url)
avdf$attr = paste(c("type bigWig 0 30", 
                    "viewLimits 0.0:30.0",
                    "noInherit on",
                    "autoScale on",
                    "yLineOnOff on",
                    "maxHeightPixels 64:32:16",
                    "priority 1"), collapse='\n')
avdf$color = sapply(avdf$COLOR, function(x){x = col2rgb(x); paste0('color ',x[1], ',', x[2],',', x[3])})
avdf$metadata = paste('metadata', paste0('"Group"="<span style="color:', avdf$COLOR,'">', avdf$GROUP, '</span>"'), sep=" ")
avdf$subGroups = paste0('subGroups view=Sig sampleType=', avdf$groupstr, ' assayType=', avdf$mark, ' dataType=', avdf$dataset) 
avdf$parent = paste0('parent EpiMapBWSigView off')
track.meta = paste(avdf$track, avdf$shortLabel, avdf$longLabel, avdf$bigDataUrl,
                   avdf$attr, avdf$color, avdf$metadata, avdf$subGroups, avdf$parent, sep="\n")

# Assays and donors:
atypes = sort(unique(avdf$mark))

# Overall metadata:
ovl.meta = c('track EpiMapBWSig',
             'compositeTrack on',
             'shortLabel Average Tracks',
             'longLabel All observed and imputed signal tracks from EpiMap',
             'group epimap',
             'subGroup1 view Views Sig=Signal',
             paste(c('subGroup2 sampleType Sample_Type', paste0(unique(avdf$groupstr), '=', unique(avdf$groupsp))), collapse=' '),
             paste(c('subGroup3 assayType Assay_Type ', paste0(atypes, '=', atypes)), collapse=' '),
             'subGroup4 dataType Data_Type Observed=Observed Imputed=Imputed',
             'dimensions dimensionX=assayType dimensionY=sampleType dimA=dataType',
             'filterComposite on',
             'sortOrder sampleType=+ assayType=+ dataType=+',
             'dragAndDrop on',
             'visibility hide',
             'priority 21',
             'type bed 3',
             'noInherit on')

# View track:
view.meta = c('track EpiMapBWSigView',
              'shortLabel Signal',
              'view Sig',
              'maxHeightPixels 64:32:16',
              'parent EpiMapBWSig',
              'type bigWig',
              'visibility dense')


genome = 'hg19'
hub.meta = c('hub EpiMapAverage',
             'shortLabel EpiMap release v0.9 (Group Averages)',
             'longLabel Average tracks per mark and sample group tracks for EpiMap v0.9',
             'useOneFile on',
             'email cboix@mit.edu',
             '\n', 
             paste('genome', genome))


# --------------------------------
# Write out the chromHMM trackHub:
# --------------------------------
thfile = 'trackHub_groupaverages.txt'

hub.line = paste(hub.meta, collapse='\n')
ovl.line = paste(ovl.meta, collapse='\n')
view.line = paste(view.meta, collapse='\n')
track.line = paste(track.meta, collapse='\n\n')

writeLines(c(hub.line, ovl.line, view.line, track.line), con=thfile, sep = "\n\n")


# ------------------------
# Make the WUSTL trackhub:
# ------------------------
avdf$shortLabel = paste("shortLabel", substr(avdf$dataset,1,3), avdf$GROUP, avdf$mark)
avdf$longLabel = paste0("longLabel Average of ", avdf$id, " ", ifelse(avdf$id == 1, "track", "tracks"),
                        " for ", avdf$GROUP, " ", avdf$mark, ' (', avdf$dataset, ')')

# TODO: Make the metadata definitions:
sdf = avdf
sdf$colorstr = sapply(sdf$COLOR, function(x){x = col2rgb(x); paste0('"pr": ',x[1], ', "pg": ', x[2],', "pb": ', x[3])})
# Base string attributes:
base.info = '{"type": "native_track", "list": [{"name": "refGene" ,"mode":"full"}]}, {"type":"coordinate_override","coord":"chr9,36329955,chr9,37537411"}'
# Make the per-mark strings:
sdf$trackline = paste0("{", paste0('"name": "', paste(avdf$GROUP, avdf$mark, paste0('(', substr(avdf$dataset,1,3), ')')), '", '),
                       # sdf$trackline = paste0("{", paste0('"name": "', substr(sdf$dataset,1,3), " ", sdf$id, " ", sdf$mark,'", '),
                       paste0('"url": "', sdf$url, '", '),
                       '"type": "bigwig", "height": 30, "mode":1, "qtc":{"anglescale":1, "height":30, "summeth":2, "thtype":1, "thmin":0, "thmax":20, "smooth":3,',
                       sdf$colorstr, '}', 
                       ', "metadata": {', 
                       paste0('"Group": "', sdf$GROUP,'",'), 
                       paste0('"Assay": "', sdf$mark,'"'), '}, ',
                       paste0('"options": {"color": "', sdf$COLOR, '"}'),
                       "}")
json.line = paste("[", base.info, ",\n ", paste(sdf$trackline, collapse=",\n "), "]")

thfile = 'wustlhub.groupaverages.json'
writeLines(json.line, con=thfile)

