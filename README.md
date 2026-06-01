# CS185-Skitty-Project: Code for C&S Bio 185 Project

CONTENT Pipeline Code that covers grabbing the data from the OneK1K cohort and GENCODE gene annotations, creating pseudobulked matrices of the single-cell RNA sequencing data from OneK1K, simulating genotype data following CONTENT paper parameters, running the CONTENT framework, and outputting various summary statistics files.

Run order: get_data.sh (download Onek1k and gene annotations, then create gene annotation .tsv file) -> make_pseudobulk.py (makes inverted pseudobulk matrices) -> fixing_pseudobulk.py (transposing pseudobulk matrices to correct CONTENT input orientation) -> run_pipeline_FINAL.sh (runs simulating_real_genotypes_FINAL.R to simulate genotype data, then run_content_FINAL.R to run the three CONTENT models and Context-by-Context framework 4 total times, one for the full-context run and once for each of 3 ablation runs, then training_summaries.R, testing_performance.R, testing_summaries.R which corrects errors in testing_performance.R, then genetic_sums.R to summarize and compute statistics for all four models).


## AI Use/Disclaimer: Use logs are not available, but the following summarizes AI usage across created code:

simulating_real_genotypes_FINAL: Gemini was used to understand how to apply reference paper simulation parameters to our simulation when drawing from

run_content_FINAL.R: Claude was used to make corrections to the script to follow reference paper parameters and procedures more exactly.

Generally: AI tools were often used to make script-wide naming convention/formatting edits when doing so manually would have taken exceedingly long, as well as to clarify the purpose of various parameters in the reference paper. Since our simulation drew from real donor IDs, information, and cell-type data while the original paper's simulation was completely disconnected from reality, AI was used to clarify the best procedures to integrate our data with the original simulation parameters. Additionally, AI provided advice on functions, packages, etc. to use to increase computational efficiency (ie. fread for large amounts of data and mclapply to run parallel workers on machine cores). 
