#!/usr/bin/R
# UI side for imputation visualization:
library(shiny)
options(scipen=45)
domain <- Sys.getenv('DOMAINNAME')
if (domain == 'broadinstitute.org'){
    #setwd('/broad/compbio/cboix/EPIMAP_ANALYSIS/db/')
    img <- '/broad/compbio/cboix/EPIMAP_ANALYSIS/img/'
    bin <- '/broad/compbio/cboix/EPIMAP_ANALYSIS/bin/'
} else {
    #setwd('~/EPIMAP_ANALYSIS/db/')
    img <- '~/EPIMAP_ANALYSIS/img/'
    bin <- '~/EPIMAP_ANALYSIS/bin/'
}

# Read in autocomplete list:
autolist = as.character(unlist(read.delim(paste0('data/celltypes.tsv'), header=F)))


shinyUI( fluidPage(titlePanel("Imputation Dataset Selection"),
                   sidebarLayout(
                                 sidebarPanel(
                                              h3("Add cell type sets:"),
                                              checkboxInput("gt2", "1. At least 2 assays", TRUE),
                                              checkboxInput("miss", "2. Missing core marks", TRUE),
                                              checkboxInput("sign", "3. Have signal tracks", TRUE),
                                              selectizeInput('txt', '4. Manually add cell types',
                                                             choices = autolist,
                                                             multiple = TRUE, options = list(create = TRUE)),
                                              checkboxInput("core", "5. Add core set", TRUE),
                                              checkboxInput("sugg", "6. Add validation set", TRUE),
                                              sliderInput("diverse", "7. Total cell types (fill rest with diverse cells):",
                                                          min = 0, max = length(autolist), step=1, value = 50),
                                              h3("Plotting:"),
                                              textInput('title', 'Title:'),
                                              selectInput("attr", "Attribute:", 
                                                          list("Availability" = "accession",
                                                               "Status" = "status",
                                                               "Project" = "award",
                                                               "BigWig" = "bw_file",
                                                               "Priority" = "prio"), "accession"),
                                              sliderInput("nbreaks", "Breaks on heatmap:",
                                                          min = 0, max = 100, step=1, value = 30),
                                              h3("Download figures/table:"),
                                              downloadButton('downAvail', label="Download Availability"),
                                              downloadButton('downHeat', label="Download Heatmap"),
                                              downloadButton('downTable', label="Download Table")
                                              ),

                                 mainPanel(
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("Data Matrix", plotOutput("Avail")),
                                                       tabPanel("Dataset Heatmap", plotOutput("Heatmap")), 
                                                       tabPanel("Data Table", dataTableOutput('Table')))
                                           )
                                 )
                   )
)
