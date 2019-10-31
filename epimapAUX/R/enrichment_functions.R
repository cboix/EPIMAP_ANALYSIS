#!/usr/bin/R
# -------------------------------
# Auxiliary enrichment functions:
# -------------------------------

#' Run hypergeometric test:
#' 
#' Give vector with number of hits, number of draws, number of white balls(hittable), number of total balls:
#' @param y Test vector of 4 values 
#' @return enrichment p-value based on hypergeometric
#' @export
run.hyper <- function(y){phyper(q=y[1]-1, m=y[3], n=y[4] - y[3], k=y[2], lower.tail=FALSE)} 


#' Run fisher's exact test:
#' 
#' @param y Vector of contingency table values
#' @return Enrichment p-value based on fisher's exact test 
#' @export
run.fisher <- function(y){
    table <-  matrix(as.numeric(c(y[1], y[2], y[3], y[4])), ncol = 2, byrow = TRUE)
    if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
    return(p)
}
