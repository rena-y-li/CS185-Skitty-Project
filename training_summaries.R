# This script summarizes training performance of all four models

#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(readr)


# Defining runs
runs <- list(
  list(name="31_contexts", file="all_genes_performance_full_F.csv"),
  list(name="30_contexts", file="all_genes_performance_no_nkp_F.csv"),
  list(name="29_contexts", file="all_genes_performance_no_nkp_nkcd56_F.csv"),
  list(name="28_contexts", file="all_genes_performance_no_all_nk_F.csv")
)

# Manually order output
method_order <- c("CONTENT_Full", "CONTENT_Shared", "CONTENT_Specific", "Context_by_Context")

# Loading all summary .csvs

all_summaries <- list()

cat("Processing training performance CSVs...\n")

for (run in runs) {
  if (!file.exists(run$file)) {
    cat(sprintf("  Warning: File '%s' not found. Skipping %s.\n", run$file, run$name))
    next
  }
  
  cat(sprintf("  Reading and summarizing %s...\n", run$file))
  
  df <- read_csv(run$file, show_col_types = FALSE)

  # Corrective R^2 calculation: instead of dropping NAs (failures), set their values to 0

  df_clean <- df %>%
    mutate(
      r2 = ifelse(is.na(r2) | n_snps_nz == 0, 0, r2),
      p_value = ifelse(is.na(p_value) | n_snps_nz == 0, 1, p_value),
      n_snps_nz = ifelse(is.na(n_snps_nz), 0, n_snps_nz)
    )
  
  # Calculate metrics
  
  run_summary <- df_clean %>%
    group_by(method) %>%
    summarise(
      n_gene_contexts = n(),
      mean_donors     = mean(n_samples, na.rm = TRUE), 
      mean_r2         = mean(r2),
      median_r2       = median(r2),
      mean_p          = mean(p_value),
      pct_r2_gt0      = mean(r2 > 0) * 100,
      mean_n_snps_nz  = mean(n_snps_nz),
      .groups = "drop"
    ) %>%
    mutate(run = run$name)
  
  all_summaries[[run$name]] <- run_summary
}

# Compute and put all summaries into one csv

if (length(all_summaries) > 0) {
  final_training_table <- bind_rows(all_summaries) %>%
    select(run, method, n_gene_contexts, mean_donors, mean_r2, median_r2, mean_p, pct_r2_gt0, mean_n_snps_nz) %>%
    arrange(run, match(method, method_order))
  
  cat("\n=== UNBIASED TRAINING PERFORMANCE SUMMARY TIBBLE ===\n")
  print(as_tibble(final_training_table), width = Inf)
  
  # Save to file for easy comparison plots later
  write_csv(final_training_table, "training_r2_unbiased_summary_F.csv")
  cat("\nSummary table saved to 'training_r2_unbiased_summary_F.csv'\n")
  
} else {
  stop("Error: None of the specific performance CSV files could be found. Please double-check your working directory.")
}
