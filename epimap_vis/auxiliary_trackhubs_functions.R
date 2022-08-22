#!/usr/bin/R
# ---------------------------------
# Functions for creating trackhubs:
# ---------------------------------

make.ucsc = function(sdf){
    # Make subtrack metadata:
    sdf$track = paste0("track BW_",  sdf$id, "_", sdf$mark,"_", substr(sdf$dataset,1,3))
    sdf$shortLabel = paste("shortLabel", substr(sdf$dataset,1,3), sdf$infoline)
    sdf$longLabel = paste0("longLabel ", sdf$id, " ", sdf$mark, " ", sdf$name, ' (', sdf$dataset, ')')
    sdf$bigDataUrl = paste0('bigDataUrl ', sdf$dir, sdf$prefix, '.bigWig')
    sdf$attr = paste(c("type bigWig 0 30", 
                       "viewLimits 0.0:30.0",
                       "noInherit on",
                       "autoScale on",
                       "yLineOnOff on",
                       "maxHeightPixels 64:32:16",
                       "priority 1"), collapse='\n')
    sdf$color = sapply(sdf$COLOR, function(x){x = col2rgb(x); paste0('color ',x[1], ',', x[2],',', x[3])})
    sdf$metadata = paste('metadata', 
                           paste0('"Epigenome ShortName"="', sdf$infoline,'"'),
                           paste0('"Epigenome FullName"="', sdf$name,'"'),
                           paste0('"Group"="<span style="color:', sdf$COLOR,'">', sdf$GROUP, '</span>"'),
                           paste0('"Type"="', sdf$type,'"'),
                           paste0('"Sex"="', sdf$sex,'"'),
                           # paste0('"Age"="', sdf$fullage,'"'),
                           paste0('"Donor"="', sdf$Donor,'"'),
                           sep=" ")
    sdf$subGroups = paste0('subGroups view=Sig sampleType=', sdf$id, ' assayType=', sdf$mark, ' dataType=', sdf$dataset, ' donor=', sdf$Donor) 
    sdf$parent = paste0('parent EpiMapBWSigView off')
    track.meta = paste(sdf$track, sdf$shortLabel, sdf$longLabel, sdf$bigDataUrl,
                       sdf$attr, sdf$color, sdf$metadata, sdf$subGroups, sdf$parent, sep="\n")

    # Assays and donors:
    adf = aggregate(id ~ mark, sdf, length)
    adf = adf[order(adf$id, decreasing=T),]
    atypes = adf$mark

    ddf = aggregate(id ~ Donor, sdf, length)
    ddf = ddf[order(ddf$id, decreasing=T),]
    dtypes = ddf$Donor

    idf = unique(sdf[,c('id','ct')])
    idf$ctspace = gsub(" ","_", idf$ct)

    # Overall metadata:
    ovl.meta = c('track EpiMapBWSig',
                 'compositeTrack on',
                 'shortLabel All Signal Tracks',
                 'longLabel All observed and imputed signal tracks from EpiMap',
                 'group epimap',
                 'subGroup1 view Views Sig=Signal',
                 paste(c('subGroup2 sampleType Sample_Type', paste0(idf$id, '=', idf$ctspace)), collapse=' '),
                 paste(c('subGroup3 assayType Assay_Type ', paste0(atypes, '=', atypes)), collapse=' '),
                 'subGroup4 dataType Data_Type Observed=Observed Imputed=Imputed',
                 paste(c('subGroup5 donor Donor ', paste0(dtypes, '=', dtypes)), collapse=' '),
                 'dimensions dimensionX=assayType dimensionY=sampleType dimA=dataType dimB=donor',
                 'filterComposite on',
                 'sortOrder cellType=+ assayType=+ donor=+',
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

    # Hub description; use one file for the full hub:
    hub.meta = c('hub EpiMapCustom',
                 'shortLabel EpiMap release v0.9 (Custom TrackHub)',
                 'longLabel Custom trackHub per mark and sample group tracks for EpiMap v0.9',
                 'useOneFile on',
                 'email cboix@mit.edu',
                 '\n', 
                 'genome hg19')

    hub.line = paste(hub.meta, collapse='\n')
    ovl.line = paste(ovl.meta, collapse='\n')
    view.line = paste(view.meta, collapse='\n')
    track.line = paste(track.meta, collapse='\n\n')
    ucschub = c(hub.line, ovl.line, view.line, track.line)

    ucschub = paste(ucschub, collapse='\n\n')
    return(ucschub)
}


# TODO: Can we automatically make data go into a dense visualization if more than X tracks 
# TODO: Allow choosing track color by datatype, group, (lifestage, etc must be careful with grey/unknown.)

make.legacy.wustl = function(sdf){
    sdf$colorstr = sapply(sdf$COLOR, function(x){x = col2rgb(x); paste0('"pr": ',x[1], ', "pg": ', x[2],', "pb": ', x[3])})
    # Base string attributes:
    base.info = '{"type": "native_track", "list": [{"name": "refGene" ,"mode":"full"}]}, {"type":"coordinate_override","coord":"chr9,36329955,chr9,37537411"}'
    # TODO: make the metadata definitions:

    # Make the per-mark strings:
    sdf$trackline = paste0("{", paste0('"name": "', substr(sdf$dataset,1,3), " ", sdf$id, " ", sdf$mark,' (', sdf$ct, ')", '),
                           # sdf$trackline = paste0("{", paste0('"name": "', substr(sdf$dataset,1,3), " ", sdf$id, " ", sdf$mark,'", '),
                           paste0('"url": "', sdf$dir, sdf$prefix, '.bigWig", '),
                           '"type": "bigwig", "height": 30, "mode":1, "qtc":{"anglescale":1, "height":30, "summeth":2, "thtype":1, "thmin":0, "thmax":20, "smooth":3,',
                           sdf$colorstr, '}', 
                           ', "metadata": {',
                           paste0('"Assay": "', sdf$mark,'", '),
                           paste0('"Sample": "', sdf$ct,'", '),
                           paste0('"Group": "', sdf$GROUP,'", '),
                           paste0('"Lifestage": "', sdf$lifestage,'", '),
                           paste0('"Sex": "', sdf$sex,'", '),
                           paste0('"Type": "', sdf$type,'"'),"},",
                           paste0('"options": {"color": "', sdf$COLOR, '"}'),
                           "}")
    json.line = paste("[", base.info, ",\n ", paste(sdf$trackline, collapse=",\n "), "]")
    return(json.line)
}
