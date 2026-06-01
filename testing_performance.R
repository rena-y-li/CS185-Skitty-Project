# This script runs performance metrics using the ground truth file created during simulation, but the method of R^2 calculation is slightly biased

library(data.table)
library(dplyr)
library(parallel)
library(readr)
library(tidyr)

# Extract ground truth and variables

cat("Loading ground truth...\n")
gt <- readRDS("sim_ground_truth.rds")

G_test        <- gt$G_test
test_donors   <- gt$test_donors
all_expr_test <- gt$all_expr_test
contexts      <- gt$contexts
gene_ids      <- gt$gene_ids
snp_ids       <- gt$snp_ids

cat("Test donors:", length(test_donors), "\n")
cat("Genes:", length(gene_ids), "\n")

# Define runs

runs <- list(
  list(name="31_contexts", dir="sim_results_full_F",
       contexts=contexts, outfile="test_r2_existing_models_long.csv"),
  list(name="30_contexts", dir="sim_results_no_nkp_F",
       contexts=setdiff(contexts, "NK_Proliferating"), outfile="test_r2_no_nkp_long.csv"),
  list(name="29_contexts", dir="sim_results_no_nkp_nkcd56_F",
       contexts=setdiff(contexts, c("NK_Proliferating", "NK_CD56bright")), outfile="test_r2_no_nkp_nkcd56_long.csv"),
  list(name="28_contexts", dir="sim_results_no_all_nk_F",
       contexts=setdiff(contexts, c("NK_Proliferating", "NK_CD56bright", "NK")), outfile="test_r2_no_all_nk_long.csv")
)


# Load SNP weights from each directory

load_all_weights <- function(run_dir, gene_ids, contexts) {
  cat("  Preloading SNP weights from", run_dir, "...\n")
  
  weights <- mclapply(gene_ids, function(gene_id) {
    gene_dir <- file.path(run_dir, gene_id)
    if (!dir.exists(gene_dir)) return(NULL)
    
    ctx_weights <- lapply(contexts, function(ctx) {
      wgt_file <- file.path(gene_dir, paste0("snp_weights__", ctx, ".csv"))
      if (!file.exists(wgt_file)) return(NULL)
      wgt <- tryCatch(fread(wgt_file, data.table=FALSE), error=function(e) NULL)
      if (is.null(wgt) || nrow(wgt) == 0) return(NULL)
      wgt
    })
    names(ctx_weights) <- contexts
    ctx_weights
  }, mc.cores=11)
  
  names(weights) <- gene_ids
  weights
}

# Compute R^2 but this is biased! Drops NAs instead of setting it to 0 if the model failed. This is corrected in testing_summaries.R which runs on the outputs of this script

compute_r2 <- function(w, G_test, snp_ids, y_true) {
  if (all(w == 0)) return(NA_real_)
  common <- intersect(names(w), snp_ids)
  if (length(common) == 0) return(NA_real_)
  G_sub <- G_test[, match(common, snp_ids), drop=FALSE]
  pred  <- as.numeric(G_sub %*% w[common])
  ct    <- tryCatch(cor.test(pred, y_true), error=function(e) NULL)
  if (is.null(ct)) return(NA_real_)
  ct$estimate^2
}



cat("Evaluating test R² across completed runs...\n")
all_results <- list()

for (run in runs) {
  cat("\nRun:", run$name, "\n")
  run_dir      <- run$dir
  run_contexts <- run$contexts
  
  if (!dir.exists(run_dir)) {
    cat(sprintf("  Skipping %s: directory %s does not exist yet.\n", run$name, run_dir))
    next
  }

  # Load SNP weights
  all_weights <- load_all_weights(run_dir, gene_ids, run_contexts)


  run_results <- mclapply(gene_ids, function(gene_id) {
    true_expr   <- all_expr_test[[gene_id]]
    gene_weights <- all_weights[[gene_id]]
    if (is.null(true_expr) || is.null(gene_weights)) return(NULL)

    gene_rows <- lapply(run_contexts, function(ctx) {
      if (!ctx %in% names(true_expr)) return(NULL)
      wgt <- gene_weights[[ctx]]
      if (is.null(wgt)) return(NULL)

      y_true <- true_expr[[ctx]][test_donors]
      if (is.na(var(y_true)) || var(y_true) == 0) return(NULL)

      # CBC weights
      cbc_w        <- setNames(wgt$cbc_weight, wgt$snp_id)
      r2_cbc       <- compute_r2(cbc_w, G_test, snp_ids, y_true)

      # CONTENT Full weights
      full_w       <- setNames(wgt$content_full_weight, wgt$snp_id)
      r2_full      <- compute_r2(full_w, G_test, snp_ids, y_true)

      # CONTENT Shared weights 
      shared_w     <- setNames(wgt$content_shared_weight, wgt$snp_id)
      r2_shared    <- compute_r2(shared_w, G_test, snp_ids, y_true)

      # CONTENT Specific weights 
      specific_w   <- setNames(wgt$content_specific_weight, wgt$snp_id)
      r2_specific  <- compute_r2(specific_w, G_test, snp_ids, y_true)

      data.frame(
        gene    = gene_id,
        context = ctx,
        cbc     = r2_cbc,
        full    = r2_full,
        shared  = r2_shared,
        specific = r2_specific,
        stringsAsFactors = FALSE
      )
    })
    bind_rows(Filter(Negate(is.null), gene_rows))
  }, mc.cores=11)

  run_df <- bind_rows(Filter(Negate(is.null), run_results))
  
  if (nrow(run_df) == 0) {
    cat("  No data parsed for this run.\n")
    next
  }

  cat("  Gene-context pairs evaluated:", nrow(run_df), "\n")
  cat("  Mean test R² CBC:            ", round(mean(run_df$cbc,  na.rm=TRUE), 4), "\n")
  cat("  Mean test R² CONTENT Full:   ", round(mean(run_df$full, na.rm=TRUE), 4), "\n")

  # Add run annotation back for the comprehensive final file
  all_results[[run$name]] <- run_df %>% mutate(run = run$name)


  # Pivots methods out to rows to cleanly integrate with testing_summaries.R
  run_long <- run_df %>%
    pivot_longer(cols = c(cbc, full, shared, specific), 
                 names_to = "method", 
                 values_to = "r2") %>%
    mutate(method = case_when(
      method == "cbc"      ~ "Context_by_Context",
      method == "full"     ~ "Full_Model",
      method == "shared"   ~ "Shared_Only",
      method == "specific" ~ "Specific_Only",
      TRUE ~ method
    ))
  
  write_csv(run_long, run$outfile)
  cat(sprintf("  Saved dedicated long-form results to %s\n", run$outfile))
  
  # Free memory before next iteration
  rm(all_weights, run_df, run_long)
  gc()
}


if (length(all_results) > 0) {
  cat("\nCombining comprehensive master summaries...\n")
  combined <- bind_rows(all_results)
  
  existing_runs <- intersect(c("31_contexts","30_contexts","29_contexts","28_contexts"), unique(combined$run))
  combined$run <- factor(combined$run, levels=existing_runs)
  
  summary_test <- combined %>%
    group_by(run) %>%
    summarise(
      n_gene_ctx         = n(),
      mean_r2_cbc        = mean(cbc,  na.rm=TRUE),
      mean_r2_full       = mean(full, na.rm=TRUE),
      median_r2_cbc      = median(cbc,  na.rm=TRUE),
      median_r2_full     = median(full, na.rm=TRUE),
      .groups = "drop"
    )
  
  cat("\n=== RUNNING TEST SET R² SUMMARY ===\n")
  print(summary_test)
  
  # Summarize for non-NK contexts only
  non_nk <- subset(combined, !context %in% c("NK", "NK_CD56bright", "NK_Proliferating"))
  if (nrow(non_nk) > 0) {
    summary_non_nk <- non_nk %>%
      group_by(run) %>%
      summarise(
        mean_r2_cbc  = mean(cbc,  na.rm=TRUE),
        mean_r2_full = mean(full, na.rm=TRUE),
        .groups = "drop"
      )
    cat("\n=== NON-NK CONTEXTS ONLY ===\n")
    print(summary_non_nk)
    write.csv(summary_non_nk, "test_r2_non_nk_summary_snps.csv", row.names=FALSE)
  }
  
  write.csv(combined,     "test_r2_results.csv",       row.names=FALSE)
  write.csv(summary_test, "test_r2_summary.csv",       row.names=FALSE)
  cat("\nMaster files saved successfully.\n")
} else {
  cat("\nNo processing directories were available. No outputs written.\n")
}
