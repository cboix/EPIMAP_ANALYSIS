source('~/EPIMAP_ANALYSIS/bin/general_EPIMAP_ANALYSIS.R')
library(jsonlite)
library(httr)

################################################################################################################
consumer_key = "R7GOMSZL";
consumer_secret = "7sccfx4rtdb6sqse";
ua <- paste(consumer_key, consumer_secret, sep=":")
secret <- RCurl::base64(paste(consumer_key, consumer_secret, sep=":"));
################################################################################################################

### This is where we obtain most of the metadata
nam = 'all_submitted_released' 
inner = 'type=Experiment&status=released&status=submitted&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens'
url <- paste("https://www.encodeproject.org/search/?", inner, "&format=json&limit=all", sep="");

if (length(grep(paste("^Annotation/experiments_TX_", nam, ".RData$", sep=""), dir())) == 0) {

    call <- GET(url, config(c("Authorization" = paste("Basic", secret))))
    obj <- fromJSON(rawToChar(call$content))
    status <- obj$facets[6,"terms"][[1]]

    if (status[status$key=='submitted','doc_count'] == 0){ 
        print('Couldn\'t access submitted data, loading previously stored JSON file') 
        obj <- fromJSON(paste0('Annotation/',nam,'.json'))
    }

    # Parse some of the data and save.
    experiments <- obj[["@graph"]]

    # =======================
    # Get all TX experiments:
    # =======================
    # Drop anything not submitted/released 
    experiments <- experiments[experiments$status %in% c('submitted','released'),]
    assays <- c('polyA mRNA RNA-seq','RNA-seq','RNA microarray')
    experiments <- experiments[experiments$assay_title %in% assays,]

    # Disambiguate treated/normal!
    desc <- experiments$description
    treat <- experiments[grep('treat',desc),] 
    experiments <- experiments[-grep('treat',desc),] 

    ct <- sort(unique(experiments$biosample_term_name))
    write.table(ct,file=paste0('Annotation/',nam,'_TX_cell_types.tsv'),row.names=F,col.names=F,quote=F,sep="\t")

    save(experiments, treat, file=paste("Annotation/experiments_TX_", nam, ".RData", sep=""))


}

# Get the links: 
links_dir <- paste(assay,"file_links", type, sep="/");
dir.create(links_dir, showWarnings=FALSE, recursive=TRUE);


