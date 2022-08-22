#!/usr/bin/R
# Check for full data availability (all + submitted)
# Up-to-date as of 06-2017
source('~/data/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
# source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
library(jsonlite)
library(httr)
type='all_submitted_released'
args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied: Need assay and datatype.")
} else {        
    print(args)
    assay = args[1]
    datatype = args[2] # NOTE: Currently doesnt work - submitter keys should be passed through curl -u
}
# TODO pull tagAlign links as well

links_dir <- paste(assay,"file_links", type, sep="/");
dir.create(links_dir, showWarnings=FALSE, recursive=TRUE);

################################################################################################################
consumer_key = "R7GOMSZL";
consumer_secret = "roctnxvjxljlfrwg";
secret <- RCurl::base64(paste(consumer_key, consumer_secret, sep=":"));
# TODO figure out how to get links. (outside of brute force)
ckey2 <- "GHSMED6N"
csecret2 <- "o3bxdwvcdxpuspoc"
secret <- RCurl::base64(paste(ckey2, csecret2, sep=":"));

################################################################################################################
load(file=paste("Annotation/experiments_", type, ".RData", sep=""))

# For if we want to merge treated/not treated:
# experiments <- rbind(experiments,treat)

# For now load full:  (process + get info about file)
f.file='Annotation/all_links_bam_tagAlign_ENCODE.rda'
if (length(grep(f.file, dir())) == 1) {
    obj_files <- fromJSON('Annotation/all_files_bam.json')  # TODO currently doesnt exist.
    files <- obj_files[["@graph"]] 
    files <- files[files$status %in% c('released','in progress', 'submitted'),]
    # Remove anything which isn't just one replicate! (most of those are compiled 1,2 or 2,3 or 1,3)
    rid <- which(sapply(files$biological_replicates,length) != 1)
    files <- files[-rid,]

    obj_files <- fromJSON('Annotation/all_files_tagAlign.json')
    tafiles <- obj_files[["@graph"]] 
    tafiles <- tafiles[tafiles$status %in% c('released','in progress', 'submitted'),]
    save(files,tafiles,file=paste0(f.file))
} else { 
    load(f.file)
}

##################################################################################################################
# Data for the specific datatype + assay
rerun=TRUE
hfile = paste0(links_dir, "/hrefs_",datatype,"_",assay,".RData")
if (assay == 'ChIP-seq'){ 
    if (length(grep(hfile, dir())) == 1 || ! rerun) {
        load(hfile)
    } else { 
        # ==========
        # Get hrefs:
        # ==========
        # Columns of interest:
        cols <- c('href','file_size','dataset','assembly','output_type','submitted_file_name','biological_replicates')

        print(paste('Getting HREFS for bamfiles of assay',assay))
        sel <- which(experiments$assay_term_name == assay);
        hrefs <- data.frame();
        pb = txtProgressBar(min = 0, max = length(sel), style = 3)
        for (i in 1:length(sel)){ 
            id <- sel[i]
            for (dataset in experiments[["@id"]][id]) {
                if (is.null(files)){ 
                    url_files <- paste("https://www.encodeproject.org/search/?type=file&dataset=", dataset, "&file_format=", filetype, "&frame=object&format=json&limit=all", sep="")
                    call <- GET(url_files, config(c("Authorization" = paste("Basic", secret))))
                    obj_files <- fromJSON(rawToChar(call$content))
                    dfiles <- obj_files[["@graph"]]
                } else { 
                    # If dataset not in files
                    dfiles <- files[files$dataset == dataset,]
                    tfiles <- tafiles[tafiles$dataset == dataset,]
                }

                # add to hrefs!
                if ((nrow(tfiles) + nrow(dfiles)) > 0){
                    hrefs <- rbind(hrefs, cbind(rbind(tfiles[, cols], dfiles[, cols]),
                                                epitope=experiments$target$label[id],
                                                cell_type=experiments$uid[id], biosamp=experiments$biosample_summary[id], 
                                                desc=experiments$description[id]))
                }
            }
            setTxtProgressBar(pb, i)
        }
        hrefs$biosamp <- gsub("\'", "_", gsub("\\/", "_", gsub("\ ", "_",hrefs$biosamp)))
        # Rename Control to WCE:
        cid <- which(hrefs$epitope == 'Control')
        hrefs$epitope <- as.character(hrefs$epitope)
        if (length(cid) > 0){
            hrefs$epitope[cid] <- 'WCE'
            egfp <- grep('eGFP', hrefs$desc[cid])
            if (length(egfp) > 0) hrefs$epitope[cid[egfp]] <- 'WCE-eGFP'
        }
        hrefs$epitope <- factor(hrefs$epitope)
        save(hrefs, file=hfile);
        close(pb)
    }

    # ==================================
    # Prioritize unprocessed and GRCh38: 
    # ==================================
    hrefs$output_type <- factor(hrefs$output_type,levels=c('alignments','unfiltered alignments'))
    hrefs$assembly <- factor(hrefs$assembly,levels=c('hg19','GRCh38'))
    hrefs <- hrefs[order(hrefs$assembly,hrefs$output_type),]
    hrefs$lbl <- paste0(hrefs$dataset,hrefs$biological_replicates)

    # ===========================
    # Write each of epitopes out:
    # ===========================
    epitopes <- sort(unique(hrefs$epitope))
    printcols <- c('file','lbl','assembly','output_type','cell_type','submitted_file_name', 'biosamp')
    pb = txtProgressBar(min = 0, max = length(epitopes), style = 3)
    for (i in 1:length(epitopes)){
        epitope = epitopes[i]
        edf <- hrefs[hrefs$epitope == epitope,]
        edf$file <- paste0("https://www.encodeproject.org", edf$href)

        # Get one file for each Experiment + Replicate
        lbls <- unique(edf$lbl) 
        df <- ldply(lbls,file=edf,function(x,file){head(file[file$lbl == x,],1) })

        # Write:
        write.table(df[,printcols], file=paste0(links_dir, "/", epitope,"_",datatype,".csv"), quote=FALSE, row.names=FALSE, sep="\t");
        setTxtProgressBar(pb, i)
    }
    close(pb)

} else { 
    # ==========================
    # For DNase-seq, WGBS, RRBS:
    # ==========================
    # NOTE most methylation studies redirect you to GEO. TODO pipe to get that output?? 
    epitope <- assay
    if (epitope == 'WGBS') assay <- 'whole-genome shotgun bisulfite sequencing'

    if (length(grep(hfile, dir())) == 1 || ! rerun) {
        load(hfile)
    } else { 
        print(paste('Getting HREFS for bamfiles of assay',assay))
        sel <- which(experiments$assay_term_name == assay);
        cols <- c('href','file_size','dataset','assembly','output_type','submitted_file_name','biological_replicates')

        hrefs <- data.frame();
        pb = txtProgressBar(min = 0, max = length(sel), style = 3)
        for (i in 1:length(sel)){ 
            id <- sel[i]
            for (dataset in experiments[["@id"]][id]) {
                if (is.null(files)){ 
                    url_files <- paste("https://www.encodeproject.org/search/?type=file&dataset=", dataset, "&file_format=", filetype, "&frame=object&format=json&limit=all", sep="")
                    call <- GET(url_files, config(c("Authorization" = paste("Basic", secret))))
                    obj_files <- fromJSON(rawToChar(call$content))
                    dfiles <- obj_files[["@graph"]]
                } else { 
                    dfiles <- files[files$dataset == dataset,]
                    tfiles <- tafiles[tafiles$dataset == dataset,]
                }

                # add to hrefs!
                if ((nrow(tfiles) + nrow(dfiles)) > 0){
                    #hrefs <- rbind(hrefs, cbind(rbind(tfiles[,cols],dfiles[,cols]),celltype=experiments$uid[id]))
                    hrefs <- rbind(hrefs, cbind(rbind(tfiles[,cols],dfiles[,cols]),cell_type=experiments$uid[id], biosamp=experiments$biosample_summary[id]))
                }
            }
            setTxtProgressBar(pb, i)
        }
        hrefs$biosamp <- gsub("\'", "_", gsub("\\/", "_", gsub("\ ", "_",hrefs$biosamp)))
        save(hrefs, file=hfile);
        close(pb)
    }

    # ==============================
    # Prioritize processed and hg19: 
    # ==============================
    hrefs$output_type <- factor(hrefs$output_type,levels=c('alignments', 'unfiltered alignments'))
    hrefs$assembly <- factor(hrefs$assembly,levels=c('hg19', 'GRCh38'))
    hrefs <- hrefs[order(hrefs$assembly,hrefs$output_type),]
    hrefs$lbl <- paste0(hrefs$dataset,hrefs$biological_replicates)

    lbls <- unique(hrefs$lbl) 

    # ===========================
    # Write each of epitopes out:
    # ===========================
    printcols <- c('file','lbl','assembly','output_type','cell_type','submitted_file_name', 'biosamp')
    hrefs$file <- paste0("https://www.encodeproject.org", hrefs$href)

    # Get one file for each Experiment + Replicate
    df <- ldply(lbls,file=hrefs,function(x,file){head(file[file$lbl == x,],1) })
    write.table(df[,printcols], file=paste0(links_dir, "/", epitope, "_",datatype,".tsv"), quote=FALSE, row.names=FALSE, sep="\t");
}
