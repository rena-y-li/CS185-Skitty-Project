
# This file runs the whole pipeline beginning from simulating genotypes and ending at summary statistics

#!/bin/bash
set -e

MAINDIR=$(pwd)


echo "=== Step 1: Running Simulation ==="

# Check if the final ground truth file for the simulation already exists
if [ -f "sim_ground_truth.rds" ]; then
    echo "Simulation files already exist. Skipping simulation to save time..."
else
# If it doesn't exist, run genotype simulation
    echo "Running simulation..."
    Rscript simulating_real_genotypes_FINAL.R
fi


echo "=== Step 2: Running CONTENT models ==="

# Run the CONTENT model training file 4 times total, first on all data 

Rscript run_content_FINAL.R \
  --expr-dir    sim_expr \
  --cov-dir     sim_covariates \
  --geno-prefix sim_geno/sim_content \
  --annot       sim_gene_annotation.tsv \
  --output-dir  sim_results_full_F

# Run again using the optional contexts parameter but remove NK_Proliferating
Rscript run_content_FINAL.R \
  --expr-dir    sim_expr \
  --cov-dir     sim_covariates \
  --geno-prefix sim_geno/sim_content \
  --annot       sim_gene_annotation.tsv \
  --output-dir  sim_results_no_nkp_F \
  --contexts    "ASDC,B_intermediate,B_memory,B_naive,CD14_Mono,CD16_Mono,CD4_CTL,CD4_Naive,CD4_Proliferating,CD4_TCM,CD4_TEM,CD8_Naive,CD8_Proliferating,CD8_TCM,CD8_TEM,Doublet,Eryth,HSPC,ILC,MAIT,NK,NK_CD56bright,Plasmablast,Platelet,Treg,cDC1,cDC2,dnT,gdT,pDC"

# Remove NK_Proliferating and NK_CD56Bright
Rscript run_content_FINAL.R \
  --expr-dir    sim_expr \
  --cov-dir     sim_covariates \
  --geno-prefix sim_geno/sim_content \
  --annot       sim_gene_annotation.tsv \
  --output-dir  sim_results_no_nkp_nkcd56_F \
  --contexts    "ASDC,B_intermediate,B_memory,B_naive,CD14_Mono,CD16_Mono,CD4_CTL,CD4_Naive,CD4_Proliferating,CD4_TCM,CD4_TEM,CD8_Naive,CD8_Proliferating,CD8_TCM,CD8_TEM,Doublet,Eryth,HSPC,ILC,MAIT,NK,Plasmablast,Platelet,Treg,cDC1,cDC2,dnT,gdT,pDC"

# Remove NK_Proliferating, NK_CD56Bright, and NK (all NK types)
Rscript run_content_FINAL.R \
  --expr-dir    sim_expr \
  --cov-dir     sim_covariates \
  --geno-prefix sim_geno/sim_content \
  --annot       sim_gene_annotation.tsv \
  --output-dir  sim_results_no_all_nk_F \
  --contexts    "ASDC,B_intermediate,B_memory,B_naive,CD14_Mono,CD16_Mono,CD4_CTL,CD4_Naive,CD4_Proliferating,CD4_TCM,CD4_TEM,CD8_Naive,CD8_Proliferating,CD8_TCM,CD8_TEM,Doublet,Eryth,HSPC,ILC,MAIT,Plasmablast,Platelet,Treg,cDC1,cDC2,dnT,gdT,pDC"


echo "=== Step 3: Copying and renaming TRAINING performance CSVs ==="

# Training performance csvs are all called "all_genes_performance.csv" and are within the output-dir subfolders
# Need to rename them based on their subfolder name and move them into the main working directory

for suffix in full no_nkp no_nkp_nkcd56 no_all_nk; do
  src="sim_results_${suffix}_F/all_genes_performance.csv"
  dst="${MAINDIR}/all_genes_performance_${suffix}_F.csv"
  
  if [ -f "$src" ]; then
    cp "$src" "$dst"
    echo "  Copied $src -> $dst"
  else
    echo "  WARNING: $src not found"
  fi
done


echo "=== Step 4: Running Training Summaries ==="

# Run script to summarize training statistics

Rscript training_summaries.R

echo "=== Step 5: Running Test Performance Evaluation ==="

# Run script to evaluate testing performance
Rscript testing_performance.R

echo "=== Step 6: Running TEST Summaries ==="

# Run script to summarize testing performance
  
Rscript testing_summaries.R

echo "=== Step 7: Running Pure Signal Summaries ==="

# Run script to summarize prediction accuracy on pure genetic signal from ground truth
  
Rscript genetic_sums.R

echo "=== Pipeline complete ==="
