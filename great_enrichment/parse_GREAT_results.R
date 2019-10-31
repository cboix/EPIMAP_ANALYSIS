#!/usr/bin/env Rscript
# -----------------------
# Parse raw GREAT outputs
# Adapted from WM
# -----------------------
domain = system("hostname -d", intern=TRUE)
if (domain == 'broadinstitute.org'){
    bindir='~/data/EPIMAP_ANALYSIS/bin/'
} else {
    bindir='~/EPIMAP_ANALYSIS/bin/'
}
library(caTools)
library(gtools)
options(scipen=20)
options(warn=2)

# cpref='cls_merge2_wH3K27ac100_raw'
cpref='cls_merge2_wH3K27ac100_300'
calldir = '/broad/compbio/cboix/EPIMAP_ANALYSIS/db/ChromHMM/calls'
mid = '/observed_aux_18_on_mixed_impobs_QCUT/ENH/clust/'
infofile = paste0(calldir, mid, cpref, '/enrichment_infofile.tsv')
resultsdir = paste0(calldir, mid , cpref, '/go_enrichment')

args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied.")
} else {
    infofile = args[[1]]
    resultsdir = args[[2]]
}

# Directories:
outputdir <- paste(resultsdir, "parsed_results", sep="/");
dir.create(outputdir, showWarnings=FALSE, recursive=TRUE);

info = read.delim(infofile, header=F, stringsAsFactors=F) 
names(info) = c('file', 'name')

ID_select <- c();
ID_select_per_onto <- list();
nams <- c();
for (i in 1:nrow(info)){
    filename = paste0(sub(".*/","", sub(".gz$", "", info$file[i])),"_results_noregions.gz")
    nam = as.character(info$name[i])
    nams = c(nams, nam)
    # Read in file:
    fn = paste0(resultsdir, '/', filename)
    if (file.exists(fn)){
        res <- try(read.delim(gzfile(fn), sep="\t", header=F, skip=1, comment.char="#", as.is=T), silent=TRUE)
        if (class(res) != 'try-error'){
            names(res) = c("Ontology", "ID", "Desc", "BinomRank", "BinomP", "BinomBonfP", "BinomFdrQ", "RegionFoldEnrich", "ExpRegions",
                           "ObsRegions", "GenomeFrac", "SetCov", "HyperRank", "HyperP", "HyperBonfP", "HyperFdrQ", "GeneFoldEnrich", 
                           "ExpGenes", "ObsGenes", "TotalGenes", "GeneSetCov", "TermCov")  # , "Regions", "Genes")
            # Filter only FC > 2 and both FDR < 0.01:
            sel <- which(res$RegionFoldEnrich > 2 & res$BinomFdrQ < 0.01 & res$HyperFdrQ < 0.01);
            print(paste(nam, length(sel), sep=" - "))
            signif_res <- res[sel,];
            ID_select <- unique(c(ID_select, res$ID[sel])); # save for later.
            signif_res_per_onto <- list();
            res <- signif_res;
            if (!is.null(res)) {
                ontos <- unique(res$Ontology);
                bogus <- sapply(ontos, function(x) {
                                    if (is.null(signif_res_per_onto[[x]])) { 
                                        signif_res_per_onto[[x]] <<- vector("list", 1);
                                    }
                                    if (is.null(ID_select_per_onto[[x]])) { 
                                        ID_select_per_onto[[x]] <<- c();
                                    }
                                    signif_res_per_onto[[x]] <<- res[res$Ontology==x,];
                                    ID_select_per_onto[[x]] <<- unique(c(ID_select_per_onto[[x]], res$Desc[res$Ontology==x]));
                           });
            }
            save(signif_res_per_onto, file=paste(outputdir, "/signif_res_per_onto_", nam, ".RData", sep=""));
        } else { print(filename) }
    } else { print(filename) }
}  


res_all <- list();
for (onto in names(ID_select_per_onto)) {
    print(onto);
    ID_select <- ID_select_per_onto[[onto]];
    res_all[[onto]] <- matrix(0, nrow=length(ID_select), ncol=length(nams), dimnames=list(ID_select,nams));
    for (i in 1:nrow(info)){
        filename = paste0(sub(".*/","", sub(".gz$", "", info$file[i])),"_results_noregions.gz")
        nam = as.character(info$name[i])
        # Obtain and store enrichment results.
        fn = paste0(resultsdir, '/', filename)
        if (file.exists(fn)){
            res <- try(read.delim(gzfile(fn), sep="\t", header=F, skip=1, comment.char="#", as.is=T), silent=TRUE)
            if (class(res) != 'try-error'){
                names(res) = c("Ontology", "ID", "Desc", "BinomRank", "BinomP", "BinomBonfP", "BinomFdrQ", "RegionFoldEnrich", "ExpRegions",
                               "ObsRegions", "GenomeFrac", "SetCov", "HyperRank", "HyperP", "HyperBonfP", "HyperFdrQ", "GeneFoldEnrich", 
                               "ExpGenes", "ObsGenes", "TotalGenes", "GeneSetCov", "TermCov")  # , "Regions", "Genes")
                pvals_max <- pmax(res$BinomFdrQ, res$HyperFdrQ);
                # WM20160301: changed to the highest common p-value
                res_all[[onto]][,nam] <- -log10(pmin(1, pvals_max[match(ID_select, res$Desc)]+1e-9)); 
                print(nam)
            } else { print(filename) }
        } else { print(filename) }
    }
    res_all[[onto]][is.na(res_all[[onto]])] <- 0;
}
save(res_all, file=paste(outputdir, "/res_all_per_onto.RData", sep=""))
