
# This file is an Rscript run from the terminal that creates summary statistics for model performance when predicting the genetic component only (no noise)

library(data.table)
library(dplyr)
library(parallel)

# Load ground truth file and create variables that will be used later

gt <- readRDS("sim_ground_truth.rds")

# Extracting columns as their own variables
G_test        <- gt$G_test
test_donors   <- gt$test_donors
all_effects   <- gt$all_effects
contexts      <- gt$contexts
gene_ids      <- gt$gene_ids
snp_ids       <- gt$snp_ids

# Coding NK contexts
nk_contexts <- c("NK", "NK_CD56bright", "NK_Proliferating")

# Defining the runs and where they are
runs <- list(
  list(name  =  "31_contexts", dir  =  "sim_results_full_F",
       contexts  =  contexts),
  list(name  =  "30_contexts", dir  =  "sim_results_no_nkp_F",
       contexts  =  setdiff(contexts, "NK_Proliferating")),
  list(name = "29_contexts", dir  =  "sim_results_no_nkp_nkcd56_F",
       contexts  =  setdiff(contexts, c("NK_Proliferating", "NK_CD56bright"))),
  list(name  =  "28_contexts", dir  =  "sim_results_no_all_nk_F",
       contexts  =  setdiff(contexts, c("NK_Proliferating", "NK_CD56bright", "NK")))
)

# Calculate what the true genetic component is to all the genes

true_genetic <- mclapply(gene_ids, function(g) {
  eff <- all_effects[[g]]
  if (is.null(eff)) return(NULL)

  # Only keep SNPs that are actually in our genotype matrix (common snps)
  
  cis_snps    <- eff$cis_snps
  common_snps <- intersect(cis_snps, snp_ids)
  if (length(common_snps)   ==  0) return(NULL)

  # Extracts cis-SNPs for each gene from the genotype matrix and gets indices
  G_cis   <- G_test[, match(common_snps, snp_ids), drop = FALSE]
  cis_pos <- match(common_snps, cis_snps)

  # Calculating the shared signal, which is the same for all the contexts
  # Matrix multiplication between cis-SNPs for the donor and the shared effect size
  
  shared_signal <- as.numeric(G_cis %*% eff$shared[cis_pos])

  # Add our specific and truly specific signals to the context-specific signals in the same way
  
  ctx_signals <- lapply(contexts, function(ctx) {
    specific_signal       <- as.numeric(G_cis %*% eff$specific[cis_pos, ctx])
    truly_specific_signal <- as.numeric(G_cis %*% eff$truly_specific[cis_pos, ctx])
    # Total signal = shared + both specific components
    signal <- shared_signal + specific_signal + truly_specific_signal
    # Mapping to test donors
    names(signal) <- test_donors
    signal
  })
  names(ctx_signals) <- contexts
  ctx_signals
}, mc.cores = 11) # Use all cores 

names(true_genetic) <- gene_ids
# Drops cis-SNPs not in the genotype matrix (failed)
true_genetic <- Filter(Negate(is.null), true_genetic)
cat("Genes with true genetic components:", length(true_genetic), "\n")


# Get the SNP weights by loading everything at once to save compute time

load_all_weights <- function(run_dir, gene_ids, contexts) {
  cat("  Loading SNP weights from", run_dir, "...\n")
  weights <- mclapply(gene_ids, function(gene_id) {
    gene_dir <- file.path(run_dir, gene_id)
    if (!dir.exists(gene_dir)) return(NULL)
    ctx_weights <- lapply(contexts, function(ctx) {
      wgt_file <- file.path(gene_dir, paste0("snp_weights__", ctx, ".csv"))
      if (!file.exists(wgt_file)) return(NULL)
      tryCatch(fread(wgt_file, data.table = FALSE), error = function(e) NULL)
    })
    # Makes nested list of all genes and all SNP weights per context within each gene
    names(ctx_weights) <- contexts
    ctx_weights
  }, mc.cores = 11)
  names(weights) <- gene_ids
  weights
}


# Helper function to get the correlation coefficient between predicted and observed   
compute_correlation_r2 <- function(pred, y) {
  if (is.null(pred) || is.null(y)) return(0)
  if (length(pred)  !=  length(y)) return(0)
  if (var(pred, na.rm = TRUE)  ==  0 || var(y, na.rm = TRUE)  ==  0) return(0)
  
  ct <- tryCatch(cor.test(pred, y), error = function(e) NULL)
  if (is.null(ct)) return(0)
  return(as.numeric(ct$estimate^2))
}


all_results <- list()
                 
# Compute the values using created functions
for (run in runs) {
  cat("\nRun:", run$name, "\n")
  run_dir      <- run$dir
  run_contexts <- run$contexts

  all_weights <- load_all_weights(run_dir, gene_ids, run_contexts)

  # First gets the simulated true signal + trained weights for each gene
  run_results <- mclapply(gene_ids, function(gene_id) {
    true_gen     <- true_genetic[[gene_id]]
    gene_weights <- all_weights[[gene_id]]

    if (is.null(true_gen)) return(NULL)

    gene_rows <- lapply(run_contexts, function(ctx) {
      # Expression vector for each donor and cell type weights
      y_gen <- true_gen[[ctx]][test_donors]
      wgt   <- if(!is.null(gene_weights)) gene_weights[[ctx]] else NULL

      # NA for genes that the mask filtered out
      has_gen <- !is.null(y_gen) && length(y_gen) > 0 && var(y_gen, na.rm = TRUE) > 0
      has_wgt <- !is.null(wgt) && nrow(wgt) > 0

      # Set all prediction variables to 0
      pred_cbc <- pred_full <- pred_shared <- pred_specific <- NULL
      if (has_wgt) {
        common_snps <- intersect(wgt$snp_id, snp_ids)
        if (length(common_snps) > 0) {
          # Matches model SNPs and matrix SNPs
          G_sub <- G_test[, match(common_snps, snp_ids), drop = FALSE]

          # Prediction for each model, using the weights produced by each
          
          cbc_w   <- wgt$cbc_weight[match(common_snps, wgt$snp_id)]
          pred_cbc <- as.numeric(G_sub %*% cbc_w)

          full_w   <- wgt$content_full_weight[match(common_snps, wgt$snp_id)]
          pred_full <- as.numeric(G_sub %*% full_w)

          shared_w   <- wgt$content_shared_weight[match(common_snps, wgt$snp_id)]
          pred_shared <- as.numeric(G_sub %*% shared_w)

          specific_w   <- wgt$content_specific_weight[match(common_snps, wgt$snp_id)]
          pred_specific <- as.numeric(G_sub %*% specific_w)
        }
      }

      # Compute correlation for each gene x context x run it calculates R^2 for each model
      
      data.frame(
        run                  =  run$name,
        gene                 =  gene_id,
        context              =  ctx,
        r2_cbc_genetic       =  if(has_gen) compute_correlation_r2(pred_cbc, y_gen) else 0,
        r2_full_genetic      =  if(has_gen) compute_correlation_r2(pred_full, y_gen) else 0,
        r2_shared_genetic    =  if(has_gen) compute_correlation_r2(pred_shared, y_gen) else 0,
        r2_specific_genetic  =  if(has_gen) compute_correlation_r2(pred_specific, y_gen) else 0,
        stringsAsFactors     =  FALSE
      )
    })
    bind_rows(gene_rows)
  }, mc.cores = 11) 

  run_df <- bind_rows(Filter(Negate(is.null), run_results))
  all_results[[run$name]] <- run_df

  # Reports mean for each model while running
  cat("  Mean R² CBC (genetic tracking):              ", round(mean(run_df$r2_cbc_genetic), 4), "\n")
  cat("  Mean R² Full (genetic tracking):             ", round(mean(run_df$r2_full_genetic), 4), "\n")
  cat("  Mean R² Shared (genetic tracking):           ", round(mean(run_df$r2_shared_genetic), 4), "\n")
  cat("  Mean R² Specific (genetic tracking):         ", round(mean(run_df$r2_specific_genetic), 4), "\n")

  rm(all_weights)
  gc()
}

# Combine all the results to final csv files 
cat("\nCombining results...\n")
combined <- bind_rows(all_results)

combined$run <- factor(combined$run,
                       levels = c("31_contexts","30_contexts",
                                "29_contexts","28_contexts"))

summary_test <- combined %>%
  group_by(run) %>%
  summarise(
    n_gene_ctx                =  n(),
    mean_r2_cbc_genetic       =  mean(r2_cbc_genetic),
    mean_r2_full_genetic      =  mean(r2_full_genetic),
    mean_r2_shared_genetic    =  mean(r2_shared_genetic),
    mean_r2_specific_genetic  =  mean(r2_specific_genetic),
    .groups  =  "drop"
  )

cat("\n  ==  =  TRUE GENETIC SIGNAL PURSUIT SUMMARY   ==  = \n")
print(summary_test)

write.csv(combined,     "test_genetic_results_pure.csv", row.names = FALSE)
write.csv(summary_test, "test_genetic_summary_pure.csv", row.names = FALSE)
cat("\nSaved test_genetic_results_pure.csv and test_genetic_summary_pure.csv\n")


