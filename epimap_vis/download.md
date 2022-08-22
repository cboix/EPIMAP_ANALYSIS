#### About:

* Main website: [http://compbio.mit.edu/epimap/](http://compbio.mit.edu/epimap/)
* Paper: Regulatory genomic circuitry of human disease loci by integrative epigenomics [(Boix _et al._, 2021)](https://doi.org/10.1038/s41586-020-03145-z)
* Original manuscript on bioRxiv: https://www.biorxiv.org/content/10.1101/810291v2
* Contact: cboix@mit.edu

#### Tracks and Metadata

* Sample metadata table ([xlsx](https://personal.broadinstitute.org/cboix/epimap/metadata/Imputation_Metadata.xlsx), [tsv](https://personal.broadinstitute.org/cboix/epimap/metadata/main_metadata_table.tsv)).
* Track download ([bigWig files hosted by WUSTL](https://epigenome.wustl.edu/epimap/data/)).
* Average per-group and per-mark tracks ([bigWig files hosted by the Broad](https://personal.broadinstitute.org/cboix/epimap/averagetracks_bygroup/)).
* Track hub files ([WUSTL - json; and UCSC - txt](https://personal.broadinstitute.org/cboix/epimap/trackhubs/)).
* Chromatin state calls using ChromHMM ([hg19](https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg19/CALLS/), [hg38](https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/)), calculated using the [18-state Roadmap model](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#exp_18state). 

####  Enhancers, Modules, and Enrichments:

* Enhancer locations: [directory](https://personal.broadinstitute.org/cboix/epimap/mark_matrices/) for enhancer/promoter/dyadic element lists; [website](http://www.meuleman.org/project/dhsindex/) for DHS indices from [Meuleman _et al._, (2020)](https://www.nature.com/articles/s41586-020-2559-3)
* Enhancer modules: [assignments](https://personal.broadinstitute.org/cboix/epimap/extended_data/clusters_data/cls_merge2_wH3K27ac100_300_assignments.loc.bed), [centers](https://personal.broadinstitute.org/cboix/epimap/extended_data/clusters_data/cls_merge2_wH3K27ac100_300_centers.tsv.gz), and [metadata directory](https://personal.broadinstitute.org/cboix/epimap/extended_data/clusters_data)
* Enrichments: [archetype motifs](https://personal.broadinstitute.org/cboix/epimap/extended_data/clusters_data/motifs_bkgdhs_cls_enrich.reduced.tsv.gz), [full motifs](https://personal.broadinstitute.org/cboix/epimap/extended_data/clusters_data/motifs_bkgdhs_cls.enrich.full.tsv.gz), [GO terms](https://personal.broadinstitute.org/cboix/epimap/extended_data/clusters_data/modules_GO_enrichments_matrix_ordered.tsv.gz)
* Motif enrichments are on JASPAR 2018, HOCOMOCOv11, and  HT-SELEX (Jolma, 2013). 
* Motifs are split into archetypes according to [Vierstra (2020)](https://www.vierstra.org/resources/motif_clustering)

#### GWAS enrichments:

* GWAS Enrichments: \***Updated** snp-centric enr. ([tsv](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/gwas_results_hg_snpcentric.tsv.gz) and [Rds](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/gwas_results_hg_snpcentric.Rds)), old enr. [tsv (1%)](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/all_enrichment_summaries_alone.tsv.gz), [tsv (0.1%)](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/all_enrichment_summaries_alone_p1fdr.tsv.gz), and [Rda](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/all_enrichment_summaries_alone.Rda)
* GWAS Figures: [tree enrichments and SNP heatmaps](https://personal.broadinstitute.org/cboix/epimap/gwas_smallfigures/) and [GWAS loci](https://personal.broadinstitute.org/cboix/epimap/gwas_smallfigures/gwas_locus_figures/).
    * All plotted links in GWAS loci are given per-folder as either tsv.gz or Rda files.
* Enhancers near lead SNPs for enriched epigenomes: [tsv](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/all_snpxenhancer_intersections.tsv.gz) or [Rda](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/all_snpxenhancer_intersections.Rda)
* Directly linked genes: [tsv](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/all_gwas_SNP_links.tsv.gz) or [Rda](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/all_gwas_SNP_links.Rda)

#### Linking:
* For each sample: [gene-enhancer links](https://personal.broadinstitute.org/cboix/epimap/links/)
* Aggregated for each sample group (Adipose, Brain, etc.): [gene-enhancer links](https://personal.broadinstitute.org/cboix/epimap/links/pergroup/)


\***2022-07-02:** We have added snp-centric GWAS enrichments to the table. Other methodologies may show inflated enrichment statistics due to enhancer overlap, which is corrected to the nominal p-value level by our correction but is not completely FDR corrected.

We have re-calculated our enrichments (here in [tsv](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/gwas_results_hg_snpcentric.tsv.gz) and [Rds](https://personal.broadinstitute.org/cboix/epimap/gwas_resources/gwas_results_hg_snpcentric.Rds) formats) with a snp-centric hypergeometric test as in [Roadmap Epigenomics (2015)](https://www.nature.com/articles/nature14248) to better control FDR. These are enriched for 226 GWAS at FDR of 1% (256 at 2%, 325 at 5%). Please contact cboix at mit dot edu with any questions.
