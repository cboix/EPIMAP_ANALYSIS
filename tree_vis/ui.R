#!/usr/bin/R
# -------------------------------------
# UI side for GWAS - tree visualization
# -------------------------------------
library(shiny)
options(scipen=45)
domain <- Sys.getenv('DOMAINNAME')
if (domain == 'broadinstitute.org'){
    maindir = '/broad/compbio/cboix/EPIMAP_ANALYSIS/'
} else {
    maindir = '~/EPIMAP_ANALYSIS/'
}
dbdir  <- paste0(maindir, 'db/')
img <- paste0(maindir, 'img/')
bin <- paste0(maindir, 'bin/')

# Read in the traits:
autolist = scan(paste0(dbdir,'expanded_traitlist.txt'), 'c', sep="\n")
uidlist = scan(paste0(dbdir,'full_uidlist.txt'), 'c', sep="\n")
snpuidlist = scan(paste0(dbdir,'full_snpuidlist.txt'), 'c', sep="\n")
# Reduce uid list based on files that are processed
pfiles = list.files(pattern='*lrtree.pdf', path='~/EPIMAP_ANALYSIS/img/gwas_tree_analysis/examples/')
pind = sapply(pfiles, function(x){
                  x = sub("_lrtree.pdf", "", x);
                  x = sub(".*parent_", "", x);
                  return(as.numeric(x)) })
uidlist = uidlist[pind]
snpuidlist = snpuidlist[pind]

shinyUI(fluidPage(
                  titlePanel("Imputation GWAS Analysis"),
                  sidebarLayout(sidebarPanel( 
                                             h3("Display pre-computed trait + PubMedID combinations:"),
                                             selectizeInput('uid', 'Select Trait/Disease',
                                                            choices = snpuidlist,
                                                            multiple = FALSE, options = list(create = TRUE)),
                                             checkboxInput("minp", "Lower threshold to -log10p > 2", FALSE),
                                             checkboxInput("genes", "Add linked genes (may not be computed)", FALSE),
                                             checkboxInput("allgenes", "Report all SNPs (sub-threshold)", FALSE),
                                             h3("Dynamically plot traits:"),
                                             selectizeInput('trait', 'Select Trait/Disease',
                                                            choices = autolist,
                                                            multiple = FALSE, options = list(create = TRUE)),
                                             h3("Plotting options:"),
                                             checkboxInput("addtree", "Plot tree (faster without)?", TRUE),
                                             checkboxInput("meta", "Add metadata", TRUE),
                                             checkboxInput("leafpval", "Add leaf hg pvals", TRUE),
                                             checkboxInput("legend", "Add legend", TRUE),
                                             checkboxInput("fulltrack", "Add track count", TRUE),
                                             checkboxInput("branch", "Hide non-significant branches", TRUE),
                                             # textInput('title', 'Title:'),
                                             h3("Download figures:"),
                                             downloadButton('downTree', label="Download Tree")
                                             ),
                                mainPanel(tabsetPanel(type = "tabs",
                                                      tabPanel("Tree (pdf)", uiOutput("Pdf")),
                                                      tabPanel("Tree (dynamic)", plotOutput("Tree")))
                                )
                  )
)
)
