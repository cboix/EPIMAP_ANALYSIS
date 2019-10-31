1. Ran Carles' script 'create_func_bed.sh' on the original cluster/chromstate pair BED files, to end up with 'functional clusters', as now stored in the directory 'BED_files'
2. Run submit_GREAT_jobs.sh
3. Run submit_parse_GREAT_results_jobs.sh
4. Run code_parse_GREAT_results_for_chromstates.R; to obtain matched result matrices with (majority) chromatin states, for each enriched gene-set.
5. Run submit_plot_GREAT_results.sh


