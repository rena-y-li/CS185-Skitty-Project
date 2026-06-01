# This file runs the CONTENT pipeline following the original paper's GitHub

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(glmnet)
  library(BEDMatrix)
  library(optparse)
  library(parallel)
})

# Creating arguments 

opt_list <- list(
  make_option("--expr-dir",       type="character", default=NULL),
  make_option("--cov-dir",        type="character", default=NULL),
  make_option("--geno-prefix",    type="character", default=NULL),
  make_option("--annot",          type="character", default=NULL),
  make_option("--output-dir",     type="character", default="results/content_models"),
  make_option("--gene",           type="character", default=NULL),
  make_option("--contexts",       type="character", default=NULL),
  make_option("--cis-window-kb",  type="integer",   default=1000),
  make_option("--alpha",          type="double",    default=0.5),
  make_option("--nfolds",         type="integer",   default=5),
  make_option("--min-snps",       type="integer",   default=2),
  make_option("--write-fusion",   action="store_true", default=FALSE),
  make_option("--seed",           type="integer",   default=185)
)

opt <- parse_args(OptionParser(option_list=opt_list))
set.seed(opt$seed)


log_msg <- function(fmt, ...) {
  cat(sprintf(paste0("[%s] ", fmt, "\n"), format(Sys.time(), "%H:%M:%S"), ...))
}

# Loading inputs (contexts from pseudobulk files)

load_contexts <- function(expr_dir, requested = NULL) {
  
  # Scan the folder to find cell types based on file names
  
  files <- list.files(expr_dir, pattern = "^pseudobulk_int__.*\\.tsv\\.gz$")
  contexts <- sub("^pseudobulk_int__(.+)\\.tsv\\.gz$", "\\1", files)
  if (!is.null(requested)) {
    contexts <- intersect(contexts, requested)
    if (length(contexts) == 0)
      stop("None of the requested contexts found in expr-dir.")
  }
  contexts
}

# Loading expression data (from pseudobulk expression)

load_all_expression <- function(expr_dir, contexts) {
  log_msg("Loading all expression data into memory...")
  expr_data <- list()
  
  # Loop through contexts and read files into memory to speed up later lookup
  
  for (ctx in contexts) {
    f  <- file.path(expr_dir, paste0("pseudobulk_int__", ctx, ".tsv.gz"))
    dt <- fread(f, data.table=FALSE)
    rownames(dt) <- dt[[1]]
    dt <- dt[, -1, drop=FALSE]
    expr_data[[ctx]] <- dt
    log_msg("  Loaded %s: %d genes x %d donors", ctx, nrow(dt), ncol(dt))
  }
  expr_data
}

# Function to get gene expression data

get_gene_expression <- function(gene_id, ctx_data) {
  # Pull out a single gene vector and skip if there are not enough donors
  if (!gene_id %in% rownames(ctx_data)) return(NULL)
  vals <- as.numeric(ctx_data[gene_id, ])
  names(vals) <- colnames(ctx_data)
  vals <- vals[!is.na(vals)]
  if (length(vals) < 10) return(NULL)
  vals
}

# Loading covariates per context

load_all_covariates <- function(cov_dir, contexts) {
  log_msg("Loading all covariate data into memory...")
  cov_data <- list()
  
  # Load variables and drop any that have no variance
  
  for (ctx in contexts) {
    f <- file.path(cov_dir, paste0("covariates__", ctx, ".tsv"))
    if (!file.exists(f)) next
    df <- fread(f, data.table=FALSE)
    rownames(df) <- as.character(df[[1]])
    df <- df[, -1, drop=FALSE]
    df <- df[, apply(df, 2, function(x) var(x, na.rm=TRUE) > 0), drop=FALSE]
    cov_data[[ctx]] <- df
  }
  cov_data
}

# Residualizing expression

residualize <- function(expr_vec, cov_df, min_donors = 20) {
  
  # Use linear regression to remove non genetic background noise
  
  common <- intersect(names(expr_vec), rownames(cov_df))
  if (length(common) < min_donors) return(NULL)
  y <- expr_vec[common]
  X <- as.matrix(cov_df[common, , drop=FALSE])
  fit <- lm(y ~ X)
  r <- residuals(fit)
  names(r) <- common
  r
}

decompose_expression <- function(resid_by_context) {
  ctx_names  <- names(resid_by_context)
  all_donors <- unique(unlist(lapply(resid_by_context, names)))

  mat <- matrix(NA_real_, nrow=length(all_donors), ncol=length(ctx_names),
                dimnames=list(all_donors, ctx_names))
  for (ctx in ctx_names) {
    v <- resid_by_context[[ctx]]
    mat[names(v), ctx] <- v
  }

  # Isolate donors that appear in at least two contexts
  
  n_ctx_per_donor <- rowSums(!is.na(mat))
  shared_donors   <- all_donors[n_ctx_per_donor >= 2]
  if (length(shared_donors) < 20) {
    warning("Fewer than 20 donors in >= 2 contexts; shared model skipped.")
    return(list(shared=NULL,
                specific=setNames(vector("list", length(ctx_names)), ctx_names)))
  }

  # Calculate the cross context shared baseline signal for each donor
  
  shared <- rowMeans(mat[shared_donors, , drop=FALSE], na.rm=TRUE)

  # Calculate the context specific deviation by subtracting the shared baseline signal
  
  specific <- lapply(ctx_names, function(ctx) {
    donors_ok <- shared_donors[!is.na(mat[shared_donors, ctx])]
    if (length(donors_ok) < 20) return(NULL)
    v <- mat[donors_ok, ctx] - shared[donors_ok]
    names(v) <- donors_ok
    v
  })
  names(specific) <- ctx_names

  list(shared=shared, specific=specific)
}

extract_cis_genotypes <- function(gene_chrom, gene_tss, bed, bim,
                                  cis_window_kb, min_snps) {
  
  # Filter for SNPs within the local search window around the gene
  
  window_bp <- cis_window_kb * 1000L
  idx <- which(as.character(bim$chrom) == as.character(gene_chrom) &
               bim$pos >= gene_tss - window_bp &
               bim$pos <= gene_tss + window_bp)
  if (length(idx) < min_snps) return(NULL)

  geno <- as.matrix(bed[, idx, drop=FALSE])
  snps <- bim[idx, ]

  # Replace missing genotypes with the mean value
  
  for (j in seq_len(ncol(geno))) {
    na_j <- is.na(geno[, j])
    if (any(na_j)) geno[na_j, j] <- mean(geno[!na_j, j])
  }

  # Drop SNPs that are completely flat across all individuals
  
  keep <- apply(geno, 2, var) > 0
  if (sum(keep) < min_snps) return(NULL)

  list(geno=geno[, keep, drop=FALSE], snps=snps[keep, ])
}

fit_enet <- function(geno, expr_vec, alpha, nfolds, min_n = 30) {
  common <- intersect(rownames(geno), names(expr_vec))
  if (length(common) < min_n) return(NULL)
  X <- geno[common, , drop=FALSE]
  y <- expr_vec[common]
  if (var(y) == 0) return(NULL)

  foldid <- sample(rep(1:nfolds, length.out = length(y)))

  # Fit elastic net model using cross validation to find predictive SNPs
  
  fit <- tryCatch(
    cv.glmnet(X, y, alpha=alpha, nfolds=nfolds, foldid=foldid, standardize=TRUE, keep=TRUE),
    error = function(e) { warning("cv.glmnet failed: ", e$message); NULL }
  )
  if (is.null(fit)) return(NULL)

  # Collect prediction statistics and optimal model coefficients
  
  lm_idx   <- which(fit$lambda == fit$lambda.min)
  cv_r2    = fit$glmnet.fit$dev.ratio[lm_idx]
  
  coef_vec <- as.numeric(coef(fit, s="lambda.min"))
  names(coef_vec) <- c("(Intercept)", colnames(X))
  
  cv_pred <- as.numeric(fit$fit.preval[, lm_idx])
  names(cv_pred) <- common

  ct <- tryCatch(cor.test(cv_pred, y), error=function(e) NULL)
  list(
    coefs      = coef_vec,
    cv_pred    = cv_pred,
    cv_r2      = cv_r2,
    pearson_r2 = if (!is.null(ct)) ct$estimate^2 else NA_real_,
    p_value    = if (!is.null(ct)) ct$p.value    else NA_real_,
    n_snps_nz  = sum(coef_vec[-1] != 0),
    n_samples  = length(common),
    fit_obj    = fit
  )
}

fit_content_full_cv <- function(shared_m, specific_m, observed_expr, snp_ids, nfolds=5, min_n = 20) {
  if (is.null(shared_m) || is.null(specific_m)) return(NULL)
  
  common <- Reduce(intersect, list(names(observed_expr),
                                   names(shared_m$cv_pred),
                                   names(specific_m$cv_pred)))
  if (length(common) < min_n) return(NULL)

  y   <- observed_expr[common]
  p_s <- shared_m$cv_pred[common]
  p_k <- specific_m$cv_pred[common]

  folds <- sample(rep(1:nfolds, length.out = length(y)))
  oof_predictions <- rep(NA_real_, length(y))
  names(oof_predictions) <- common

  # Blend shared and specific cross validation predictions using linear regression
  
  for(f in 1:nfolds) {
    test_idx <- which(folds == f)
    train_df <- data.frame(y = y[-test_idx], p_s = p_s[-test_idx], p_k = p_k[-test_idx])
    test_df  <- data.frame(p_s = p_s[test_idx], p_k = p_k[test_idx])
    
    meta_lm <- lm(y ~ p_s + p_k, data = train_df)
    oof_predictions[test_idx] <- predict(meta_lm, newdata = test_df)
  }

  ct <- tryCatch(cor.test(oof_predictions, y), error=function(e) NULL)

  # Run final refit to get the final mixing weights
                 
  full_lm_fit <- lm(y ~ p_s + p_k)
  beta_s  <- coef(full_lm_fit)["p_s"]; if (is.na(beta_s)) beta_s <- 0
  beta_k  <- coef(full_lm_fit)["p_k"]; if (is.na(beta_k)) beta_k <- 0

  # Combine individual SNP effects into a single composite weight value
                 
  w_s <- shared_m$coefs[-1][snp_ids]; w_s[is.na(w_s)] <- 0
  w_k <- specific_m$coefs[-1][snp_ids]; w_k[is.na(w_k)] <- 0
  combined_weights <- (beta_s * w_s) + (beta_k * w_k)
  names(combined_weights) <- snp_ids

  list(
    combined_weights = combined_weights,
    beta_shared      = beta_s,
    beta_specific    = beta_k,
    pearson_r2       = if (!is.null(ct)) ct$estimate^2 else NA_real_,
    p_value          = if (!is.null(ct)) ct$p.value    else NA_real_,
    n_samples        = length(common)
  )
}

write_fusion_weights <- function(gene_id, context, snps_df, cbc_m, full_m,
                                 shared_m, specific_m, gene_out) {
  
  # Export weights in the specific format required by FUSION software (for TWAS analysis further down, we did not include this in the project)
  
  snp_ids <- snps_df$snp_id
  n_snps  <- length(snp_ids)

  wgt_cols <- list()
  if (!is.null(cbc_m))
    wgt_cols[["enet_cbc"]] <- cbc_m$coefs[-1][snp_ids]
  if (!is.null(full_m))
    wgt_cols[["enet_content"]] <- full_m$combined_weights[snp_ids]
  if (length(wgt_cols) == 0) return(invisible(NULL))

  wgt.matrix <- do.call(cbind, lapply(wgt_cols, function(w) {
    out <- rep(0, n_snps); names(out) <- snp_ids
    w2  <- w[!is.na(w)]; out[names(w2)] <- w2; out
  }))
  rownames(wgt.matrix) <- snp_ids

  snps <- data.frame(SNP=snps_df$snp_id, CHR=snps_df$chrom,
                     BP=snps_df$pos, A1=snps_df$a1, A2=snps_df$a2,
                     stringsAsFactors=FALSE)

  perf_rows <- list()
  if (!is.null(cbc_m))
    perf_rows[["enet_cbc"]]     <- c(rsq=cbc_m$pearson_r2,  pval=cbc_m$p_value)
  if (!is.null(full_m))
    perf_rows[["enet_content"]] <- c(rsq=full_m$pearson_r2, pval=full_m$p_value)
  cv.performance <- do.call(rbind, perf_rows)

  save(wgt.matrix, snps, cv.performance,
       file=file.path(gene_out, paste0(gene_id, "__", context, ".wgt.RDat")))
}

run_gene <- function(gene_id, gene_chrom, gene_tss,
                     contexts, expr_data, cov_data,
                     bed, bim,
                     cis_window_kb, alpha, nfolds, min_snps,
                     write_fusion, out_dir) {

  gene_out <- file.path(out_dir, gene_id)
  dir.create(gene_out, recursive=TRUE, showWarnings=FALSE)

  # Check expression across contexts for this gene
  
  expr_by_ctx <- list()
  for (ctx in contexts) {
    if (!ctx %in% names(expr_data)) next
    v <- get_gene_expression(gene_id, expr_data[[ctx]])
    if (!is.null(v)) expr_by_ctx[[ctx]] <- v
  }
  valid_ctx <- names(expr_by_ctx)
  if (length(valid_ctx) < 2) {
    log_msg("  SKIP %s: expression in < 2 contexts", gene_id)
    return(NULL)
  }

  # Run residualization to clean raw data
  
  resid_by_ctx <- list()
  for (ctx in valid_ctx) {
    if (!ctx %in% names(cov_data)) next
    r <- residualize(expr_by_ctx[[ctx]], cov_data[[ctx]])
    if (!is.null(r)) resid_by_ctx[[ctx]] <- r
  }
  if (length(resid_by_ctx) < 2) {
    log_msg("  SKIP %s: residualization left < 2 contexts", gene_id)
    return(NULL)
  }

  # Split cleaned signals into shared and context components
  
  decomp <- decompose_expression(resid_by_ctx)

  # Pull nearby genetic data
  
  geno_result <- extract_cis_genotypes(gene_chrom, gene_tss, bed, bim,
                                       cis_window_kb, min_snps)
  if (is.null(geno_result)) {
    log_msg("  SKIP %s: < %d cis-SNPs within %d kb", gene_id, min_snps, cis_window_kb)
    return(NULL)
  }
  geno_mat <- geno_result$geno
  snps_df  <- geno_result$snps
  snp_ids  <- snps_df$snp_id

  # Fit the baseline shared model
  shared_model <- if (!is.null(decomp$shared)) {
    fit_enet(geno_mat, decomp$shared, alpha, nfolds)
  } else NULL

  perf_rows   <- list()

  if (!is.null(shared_model)) {
    perf_rows[["SHARED"]] <- data.frame(
      gene=gene_id, context="SHARED", method="CONTENT_Shared",
      r2=shared_model$pearson_r2, p_value=shared_model$p_value,
      cv_dev_r2=shared_model$cv_r2, n_snps_nz=shared_model$n_snps_nz,
      n_samples=shared_model$n_samples, stringsAsFactors=FALSE)
  }

  # Loop contexts to train Context by Context, Specific and Full models
  
  for (ctx in names(resid_by_ctx)) {
    cbc_m  <- fit_enet(geno_mat, resid_by_ctx[[ctx]], alpha, nfolds)
    spec_m <- if (!is.null(decomp$specific[[ctx]])) {
      fit_enet(geno_mat, decomp$specific[[ctx]], alpha, nfolds)
    } else NULL
    
    full_m <- fit_content_full_cv(shared_model, spec_m,
                                  resid_by_ctx[[ctx]], snp_ids, nfolds=nfolds)

    if (!is.null(cbc_m))
      perf_rows[[paste0(ctx, "_CBC")]] <- data.frame(
        gene=gene_id, context=ctx, method="Context_by_Context",
        r2=cbc_m$pearson_r2, p_value=cbc_m$p_value,
        cv_dev_r2=cbc_m$cv_r2, n_snps_nz=cbc_m$n_snps_nz,
        n_samples=cbc_m$n_samples, stringsAsFactors=FALSE)

    if (!is.null(spec_m))
      perf_rows[[paste0(ctx, "_SPEC")]] <- data.frame(
        gene=gene_id, context=ctx, method="CONTENT_Specific",
        r2=spec_m$pearson_r2, p_value=spec_m$p_value,
        cv_dev_r2=spec_m$cv_r2, n_snps_nz=spec_m$n_snps_nz,
        n_samples=spec_m$n_samples, stringsAsFactors=FALSE)

    if (!is.null(full_m))
      perf_rows[[paste0(ctx, "_FULL")]] <- data.frame(
        gene=gene_id, context=ctx, method="CONTENT_Full",
        r2=full_m$pearson_r2, p_value=full_m$p_value,
        cv_dev_r2=NA_real_,
        n_snps_nz=sum(abs(full_m$combined_weights) > 0, na.rm=TRUE),
        n_samples=full_m$n_samples, stringsAsFactors=FALSE)


    wgt_df <- data.frame(
      snp_id=snp_ids, chrom=snps_df$chrom, pos=snps_df$pos,
      a1=snps_df$a1, a2=snps_df$a2,
      cbc_weight            = if (!is.null(cbc_m))  cbc_m$coefs[-1][snp_ids]          else NA_real_,
      content_full_weight   = if (!is.null(full_m)) full_m$combined_weights[snp_ids]  else NA_real_,
      content_shared_weight = if (!is.null(shared_model)) shared_model$coefs[-1][snp_ids] else NA_real_,
      content_specific_weight = if (!is.null(spec_m)) spec_m$coefs[-1][snp_ids]       else NA_real_,
      stringsAsFactors=FALSE
    )
    
    # Secure structural zeroes for all model types uniformly
      
    wgt_df$cbc_weight[is.na(wgt_df$cbc_weight)]                   <- 0
    wgt_df$content_full_weight[is.na(wgt_df$content_full_weight)] <- 0
    wgt_df$content_shared_weight[is.na(wgt_df$content_shared_weight)]     <- 0
    wgt_df$content_specific_weight[is.na(wgt_df$content_specific_weight)] <- 0
    
    write.csv(wgt_df,
              file.path(gene_out, paste0("snp_weights__", ctx, ".csv")),
              row.names=FALSE)

    if (write_fusion)
      write_fusion_weights(gene_id, ctx, snps_df, cbc_m, full_m,
                           shared_model, spec_m, gene_out)
  }

  # Export model stats for this gene
      
  perf_df <- bind_rows(perf_rows)
  write.csv(perf_df, file.path(gene_out, "model_performance.csv"),
            row.names=FALSE)

  log_msg("  DONE %s — %d contexts, %d cis-SNPs, %d model rows",
          gene_id, length(resid_by_ctx), ncol(geno_mat), nrow(perf_df))
  perf_df
}

main <- function() {

  # Verify user command line arguments
  
  required <- c("expr-dir", "cov-dir", "geno-prefix", "annot")
  missing_args <- required[sapply(required, function(a) is.null(opt[[a]]))]
  if (length(missing_args) > 0)
    stop("Missing required arguments: --", paste(missing_args, collapse=", --"))

  expr_dir    <- opt$`expr-dir`
  cov_dir     <- opt$`cov-dir`
  geno_prefix <- opt$`geno-prefix`
  gene_annot  <- opt$annot
  out_dir     <- opt$`output-dir`

  # Make sure input files physically exist
                                  
  for (p in c(expr_dir, cov_dir, gene_annot))
    if (!file.exists(p)) stop("Path does not exist: ", p)
  for (ext in c(".bed", ".bim", ".fam")) {
    f <- paste0(geno_prefix, ext)
    if (!file.exists(f)) stop("Missing PLINK file: ", f)
  }

  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  # Figure out which cell types or contexts to process
                                  
  requested_ctx <- if (!is.null(opt$contexts)) {
    trimws(strsplit(opt$contexts, ",")[[1]])
  } else NULL
  contexts <- load_contexts(expr_dir, requested_ctx)
  log_msg("Contexts (%d): %s", length(contexts), paste(contexts, collapse=", "))

  # Load data blocks upfront
                                  
  expr_data <- load_all_expression(expr_dir, contexts)
  cov_data  <- load_all_covariates(cov_dir, contexts)

  # Check gene annotation layout mapping
                                  
  annot <- fread(gene_annot, data.table=FALSE)
  req_cols <- c("gene_id", "chrom", "tss")
  miss_cols <- setdiff(req_cols, names(annot))
  if (length(miss_cols) > 0)
    stop("--annot missing columns: ", paste(miss_cols, collapse=", "))

  # Set script to process one gene or all genes in the table
                                  
  genes_to_run <- if (!is.null(opt$gene)) {
    sub <- annot[annot$gene_id == opt$gene, ]
    if (nrow(sub) == 0) stop("Gene not found: ", opt$gene)
    sub
  } else annot

  # Connect to PLINK files
  log_msg("Opening PLINK genotypes: %s", geno_prefix)
  bed <- BEDMatrix(geno_prefix, simple_names=TRUE)
  bim <- fread(paste0(geno_prefix, ".bim"), data.table=FALSE,
               col.names=c("chrom","snp_id","cm","pos","a1","a2"))

  log_msg("Running CONTENT: %d genes x %d contexts",
          nrow(genes_to_run), length(contexts))

  # Use all CPU cores
                                  
  all_perf <- mclapply(seq_len(nrow(genes_to_run)), function(i) {
    row  <- genes_to_run[i, ]
    gid  <- row$gene_id
    gchr <- row$chrom
    gtss <- as.integer(row$tss)
    log_msg("[%d/%d] %s  chr%s:%d", i, nrow(genes_to_run), gid, gchr, gtss)

    tryCatch(
      run_gene(gid, gchr, gtss, contexts, expr_data, cov_data,
               bed, bim,
               opt$`cis-window-kb`, opt$alpha, opt$nfolds, opt$`min-snps`,
               opt$`write-fusion`, out_dir),
      error = function(e) {
        log_msg("  ERROR in %s: %s", gid, conditionMessage(e))
        NULL
      }
    )
  }, mc.cores = 11)

  names(all_perf) <- genes_to_run$gene_id
  all_perf <- Filter(Negate(is.null), all_perf)

  # Merge final results and export summary
  if (length(all_perf) > 0) {
    combined <- bind_rows(all_perf)
    out_csv  <- file.path(out_dir, "all_genes_performance.csv")
    write.csv(combined, out_csv, row.names=FALSE)
    log_msg("Summary written: %s", out_csv)
  }

  log_msg("Finished. %d / %d genes completed.", length(all_perf), nrow(genes_to_run))
}

main()
