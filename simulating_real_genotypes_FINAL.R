
# This file goes through the entire genotype simulation

#!/usr/bin/env Rscript

library(data.table)
library(MASS)

params <- list(
  n_donors        = 981,
  n_genes         = 1000,
  snps_per_gene   = 500,     # 500 SNPs each
  h2_shared       = 0.15,    # 15% of expression variance driven by shared genetics
  h2_specific     = 0.05,    # 5% of variance driven by background context variation
  h2_nk           = 0.10,    # 10% variance (heritability) specifically for NK cells
  pi_causal       = 0.05,    # 5% of 500 SNPs = 25 true causal SNPs per gene
  prop_specific   = 0.8,     # 80% of shared causal SNPs also have a context-specific effect
  lambda          = 1,       
  rho             = 0.3,     # 30% environmental correlation between contexts
  train_frac      = 0.8,     # Train test split
  causal_contexts = c("NK", "NK_CD56bright", "NK_Proliferating"),
  seed            = 185
)

set.seed(params$seed)


# Load real donor information from pseudobulk matrices (cell types as contexts)

cat("Loading contexts and donors...\n")
pb_dir <- "fixed_pseudobulk"
expr_files <- list.files(pb_dir, pattern="pseudobulk_int__.*\\.tsv\\.gz", full.names=TRUE)
contexts   <- gsub(".*pseudobulk_int__(.+)\\.tsv\\.gz", "\\1", expr_files)
n_contexts <- length(contexts)
cat("Contexts:", n_contexts, "\n")

# Get unique donor IDs across all expression files
cat("Extracting all unique donor IDs...\n")
all_donors_list <- lapply(expr_files, function(f) {
  dt <- fread(f, data.table=FALSE, nrows=1)
  colnames(dt)[-1] # Exclude the gene_id column
})
donors <- unique(unlist(all_donors_list))
n_donors <- length(donors)

# Troubleshooting
if(n_donors != params$n_donors) {
  cat("Note: Found", n_donors, "unique donors, which differs from params$n_donors (", params$n_donors, ")\n")
} else {
  cat("Donors successfully matched:", n_donors, "\n")
}

# Get real donor counts per context
cat("Getting real donor counts per context...\n")
real_donor_counts <- sapply(expr_files, function(f) {
  dt <- fread(f, data.table=FALSE, nrows=1)
  ncol(dt) - 1
})
names(real_donor_counts) <- contexts

# Make sure the causal contexts exist and can be matched
missing_ctx <- setdiff(params$causal_contexts, contexts)
if (length(missing_ctx) > 0)
  stop("Causal contexts not found: ", paste(missing_ctx, collapse=", "))
cat("Causal NK contexts:", paste(params$causal_contexts, collapse=", "), "\n")

# Load the gene annotation file 

cat("Loading simulation annotation...\n")
annot_sim <- read.table("sim_gene_annotation.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
gene_ids  <- annot_sim$gene_id
n_genes   <- nrow(annot_sim)
cat("Genes:", n_genes, "\n")

# Create the empirical gene mask to simulate real biological sparsity

mask_file <- "empirical_gene_mask.rds"
if (!file.exists(mask_file)) {
  cat("empirical_gene_mask.rds not found. Generating dynamically from pseudobulk...\n")
  
  # Load reference file to get the list of real gene Ensembl IDs
  ref_df <- fread(expr_files[1], data.table=FALSE)
  all_genome_genes <- ref_df$gene_id[!is.na(ref_df$gene_id) & ref_df$gene_id != ""]
  
  if (length(all_genome_genes) < n_genes) {
    stop("Error")
  }
  
  # Select a random subset of genes across chromosomes (1000 for this project)
  real_ensembl_ids <- sample(all_genome_genes, n_genes)
  
  # Empirical mask structure, all 0's for now (Rows = generic simulated IDs, Cols = contexts)
  mask <- matrix(0, nrow=n_genes, ncol=n_contexts, dimnames=list(gene_ids, contexts))
  
  for (ctx in contexts) {
    ctx_file <- file.path(pb_dir, paste0("pseudobulk_int__", ctx, ".tsv.gz"))
    if (!file.exists(ctx_file)) next

    # Rows are the sampled gene IDs
    df <- fread(ctx_file, data.table=FALSE)
    rownames(df) <- df$gene_id

    # Calculates variance of gene expression across people, if it's 0 it is considered inactive in the cell type
    common_genes <- intersect(real_ensembl_ids, rownames(df))
    gene_vars <- apply(df[common_genes, -1, drop=FALSE], 1, var, na.rm=TRUE)
    active_real_genes <- common_genes[!is.na(gene_vars) & gene_vars > 0]

    # These are the locations of the ensembl genes that are active
    active_indices <- which(real_ensembl_ids %in% active_real_genes)
    active_sim_genes <- gene_ids[active_indices]
    
    mask[active_sim_genes, ctx] <- 1
    cat("  Context:", ctx, "| Active genes:", length(active_sim_genes), "/", n_genes, "\n")
  }
  saveRDS(mask, mask_file)
  cat("Saved 'empirical_gene_mask.rds'\n")
} else {
  cat("Loading existing biological expression mask...\n")
  mask <- readRDS(mask_file)
}

# Simulate genotypes with base pair window
cat("Simulating targeted cis-window genotypes...\n")
snps_per_gene <- params$snps_per_gene
window_bp <- 1000000L

# Variable setup
snp_chrom <- character()
snp_pos   <- integer()
snp_ids   <- character()
G_list    <- list()

for (i in seq_len(n_genes)) {
  gene_chr <- as.character(annot_sim$chrom[i])
  gene_tss <- annot_sim$tss[i]
  
  pos <- sample((gene_tss - window_bp):(gene_tss + window_bp), snps_per_gene, replace=FALSE)
  pos <- sort(pos)
  
  snp_chrom <- c(snp_chrom, rep(gene_chr, snps_per_gene))
  snp_pos   <- c(snp_pos, pos)
  
  current_ids <- paste0("rs_sim_", gene_chr, "_", pos)
  snp_ids <- c(snp_ids, current_ids)
  
  maf <- runif(snps_per_gene, 0.05, 0.50)
  G_gene <- sapply(maf, function(p) rbinom(n_donors, 2, p))
  G_list[[i]] <- G_gene
}

G <- do.call(cbind, G_list)
G <- scale(G)
rownames(G) <- donors
colnames(G) <- snp_ids
cat("Total targeted SNPs simulated:", ncol(G), "\n")

# Picking train and test donors randomly

cat("Splitting train/test...\n")
n_train      <- floor(n_donors * params$train_frac)
train_idx    <- sample(seq_len(n_donors), n_train)
test_idx     <- setdiff(seq_len(n_donors), train_idx)
train_donors <- donors[train_idx]
test_donors  <- donors[test_idx]
G_train      <- G[train_idx, ]
G_test       <- G[test_idx,  ]


# Assigning donors per context based on real donor counts

cat("Assigning donors per context...\n")
context_donors <- lapply(contexts, function(ctx) {
  n_ctx <- round(real_donor_counts[ctx] * params$train_frac)
  n_ctx <- min(n_ctx, length(train_donors))
  sample(train_donors, n_ctx, replace=FALSE)
})
names(context_donors) <- contexts


# Simulating effect size and expression per gene using empirical mask

simulate_gene <- function(gene_id, G_train, G_test, annot_sim,
                          contexts, context_donors, params,
                          snp_chrom, snp_pos, snp_ids, mask, train_donors, test_donors) {

  gene_row  <- annot_sim[annot_sim$gene_id == gene_id, ]
  gene_chr  <- as.character(gene_row$chrom)
  gene_tss  <- gene_row$tss
  window_bp <- 1000000L

  cis_idx <- which(as.character(snp_chrom) == gene_chr &
                   snp_pos >= gene_tss - window_bp &
                   snp_pos <= gene_tss + window_bp)
  if (length(cis_idx) < 2) return(NULL)

  G_cis_train <- G_train[, cis_idx, drop=FALSE]
  G_cis_test  <- G_test[,  cis_idx, drop=FALSE]
  M           <- ncol(G_cis_train)
  n_ctx       <- length(contexts)

  # Causal identifiers, sampling from binomial
  Im         <- rbinom(M, 1, params$pi_causal)
  causal_idx <- which(Im == 1)
  if (length(causal_idx) == 0) causal_idx <- sample(M, 1)
  non_causal_idx <- setdiff(seq_len(M), causal_idx)

  # Shared regulatory background, Normal
  sd_shared   <- sqrt(params$h2_shared / (M * params$pi_causal))
  beta_shared <- rep(0, M)
  beta_shared[causal_idx] <- rnorm(length(causal_idx), 0, sd_shared)

  # Context-specific modifications on shared SNPs
  beta_specific <- matrix(0, nrow=M, ncol=n_ctx)
  colnames(beta_specific) <- contexts

  for (ctx in contexts) {
    h2_specific_ctx <- if (ctx %in% params$causal_contexts) params$h2_specific + params$h2_nk else params$h2_specific
    sd_specific_ctx <- sqrt(h2_specific_ctx / (params$lambda * M * params$pi_causal))

    n_sp <- rbinom(1, length(causal_idx), params$prop_specific)
    if (n_sp > 0) {
      sp_idx <- sample(causal_idx, min(n_sp, length(causal_idx)))
      beta_specific[sp_idx, ctx] <- rnorm(length(sp_idx), mean=0, sd=sd_specific_ctx)
    }
  }

  # Truly context-specific alternative architectures
  beta_truly_specific <- matrix(0, nrow=M, ncol=n_ctx)
  colnames(beta_truly_specific) <- contexts

  if (length(non_causal_idx) > 0) {
    for (ctx in contexts) {
      n_true_sp <- if (ctx %in% params$causal_contexts) rpois(1, 2) else rpois(1, 1)
      if (n_true_sp > 0) {
        n_true_sp <- min(n_true_sp, length(non_causal_idx))
        true_sp_idx <- sample(non_causal_idx, n_true_sp)
        
        h2_specific_ctx <- params$h2_specific
        sd_specific_ctx <- sqrt(h2_specific_ctx / (params$lambda * M * params$pi_causal))
        beta_truly_specific[true_sp_idx, ctx] <- rnorm(n_true_sp, mean=0, sd=sd_specific_ctx)
      }
    }
  }

  # Noise covariance setup
  sigma2_mat <- matrix(0, nrow=n_ctx, ncol=n_ctx)
  for (ci in seq_along(contexts)) {
    ctx <- contexts[ci]
    h2_total <- params$h2_shared + params$h2_specific + (if (ctx %in% params$causal_contexts) params$h2_nk else 0)
    sigma2_mat[ci, ci] <- max(0.01, 1 - h2_total)
  }

  Sigma <- matrix(0, nrow=n_ctx, ncol=n_ctx)
  for (ci in seq_along(contexts)) {
    for (cj in seq_along(contexts)) {
      if (ci == cj) {
        Sigma[ci, cj] <- sigma2_mat[ci, ci]
      } else {
        Sigma[ci, cj] <- params$rho * sqrt(sigma2_mat[ci, ci] * sigma2_mat[cj, cj])
      }
    }
  }

  # Sub-generator enforcing empirical sparsity matrix masks
  make_expr <- function(G_cis, target_donors, context_donors_local, is_test=FALSE) {
    shared_signal         <- G_cis %*% beta_shared
    specific_signal       <- G_cis %*% beta_specific
    truly_specific_signal <- G_cis %*% beta_truly_specific
    epsilon               <- mvrnorm(n=nrow(G_cis), mu=rep(0, n_ctx), Sigma=Sigma)

    E <- as.numeric(shared_signal) + specific_signal + truly_specific_signal + epsilon
    rownames(E) <- rownames(G_cis)
    colnames(E) <- contexts

    expr_by_ctx <- lapply(contexts, function(ctx) {
      # If marked inactive by mask, force expression vector to absolute zero variance
      if (mask[gene_id, ctx] == 0) {
        out_zeros <- rep(0, if(is_test) length(target_donors) else length(intersect(context_donors_local[[ctx]], rownames(E))))
        names(out_zeros) <- if(is_test) target_donors else intersect(context_donors_local[[ctx]], rownames(E))
        return(out_zeros)
      }
      
      if (is_test) {
        return(E[target_donors, ctx])
      } else {
        ctx_donors <- intersect(context_donors_local[[ctx]], rownames(E))
        return(E[ctx_donors, ctx])
      }
    })
    names(expr_by_ctx) <- contexts
    expr_by_ctx
  }

  list(
    expr_train          = make_expr(G_cis_train, train_donors, context_donors, is_test=FALSE),
    expr_test           = make_expr(G_cis_test,  test_donors, 
                                    setNames(lapply(contexts, function(ctx) test_donors), contexts), is_test=TRUE),
    beta_shared        = beta_shared,
    beta_specific      = beta_specific,
    beta_truly_specific = beta_truly_specific,
    cis_snp_ids        = snp_ids[cis_idx],
    causal_idx         = causal_idx
  )
}

cat("Simulating expression for", n_genes, "genes...\n")
all_effects    <- list()
all_expr_train <- list()
all_expr_test  <- list()

for (i in seq_len(n_genes)) {
  if (i %% 100 == 0) cat("  Gene", i, "/", n_genes, "\n")
  result <- simulate_gene(gene_ids[i], G_train, G_test, annot_sim,
                          contexts, context_donors, params,
                          snp_chrom, snp_pos, snp_ids, mask, train_donors, test_donors)
  if (is.null(result)) {
    cat("  SKIP", gene_ids[i], "— no cis-SNPs\n")
    next
  }
  all_effects[[gene_ids[i]]] <- list(
    shared           = result$beta_shared,
    specific         = result$beta_specific,
    truly_specific   = result$beta_truly_specific,
    cis_snps         = result$cis_snp_ids,
    causal_idx       = result$causal_idx
  )
  all_expr_train[[gene_ids[i]]] <- result$expr_train
  all_expr_test[[gene_ids[i]]]  <- result$expr_test
}

valid_genes <- names(all_expr_train)
cat("Genes with cis-SNPs:", length(valid_genes), "\n")

# Creating expression files for training donors

cat("Writing expression files...\n")
dir.create("sim_expr", showWarnings=FALSE)

for (ctx in contexts) {
  ctx_donors <- context_donors[[ctx]]

  expr_mat <- sapply(valid_genes, function(g) {
    v <- all_expr_train[[g]][[ctx]]
    out <- rep(NA_real_, length(ctx_donors))
    names(out) <- ctx_donors
    out[names(v)] <- v
    out
  })
  expr_df <- as.data.frame(t(expr_mat))
  expr_df <- cbind(gene_id=valid_genes, expr_df)
  colnames(expr_df) <- c("gene_id", ctx_donors)

  gz <- gzfile(file.path("sim_expr", paste0("pseudobulk_int__", ctx, ".tsv.gz")), "w")
  write.table(expr_df, gz, sep="\t", row.names=FALSE, quote=FALSE)
  close(gz)
}

# Subsetting covariates to only have donors who actually appear for each context

cat("Writing covariate files...\n")
dir.create("sim_covariates", showWarnings=FALSE)

for (ctx in contexts) {
  cov_path <- file.path("covariates_per_context", paste0("covariates__", ctx, ".tsv"))
  if(!file.exists(cov_path)) next
  cov <- read.table(cov_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  ctx_donors <- context_donors[[ctx]]
  cov_sub    <- cov[cov[[1]] %in% ctx_donors, ]
  write.table(cov_sub,
              file.path("sim_covariates", paste0("covariates__", ctx, ".tsv")),
              sep="\t", row.names=FALSE, quote=FALSE)
}

# Creating PLINK files for CONTENT

cat("Writing PLINK files...\n")
dir.create("sim_geno", showWarnings=FALSE)

write_plink_binary_v2 <- function(geno_mat, snp_chrom, snp_pos, out_prefix) {
  n_samp  <- nrow(geno_mat)
  n_snps  <- ncol(geno_mat)
  donors  <- rownames(geno_mat)
  snp_ids <- colnames(geno_mat)

  fam <- data.frame(FID=donors, IID=donors, PAT=0L, MAT=0L, SEX=0L, PHENO=-9L)
  write.table(fam, paste0(out_prefix, ".fam"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")

  bim <- data.frame(CHR=snp_chrom, SNP=snp_ids, CM=0L, POS=snp_pos, A1="A", A2="G")
  write.table(bim, paste0(out_prefix, ".bim"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

  geno_raw <- round(sweep(sweep(geno_mat, 2, attr(geno_mat, "scaled:scale"), "*"), 2, attr(geno_mat, "scaled:center"), "+"))
  geno_raw <- pmax(0L, pmin(2L, as.integer(geno_raw)))
  dim(geno_raw) <- dim(geno_mat)

  encode <- function(x) ifelse(is.na(x), 1L, ifelse(x == 0L, 0L, ifelse(x == 1L, 2L, 3L)))
  bytes_per_snp <- ceiling(n_samp / 4L)
  pad_n         <- bytes_per_snp * 4L - n_samp

  con <- file(paste0(out_prefix, ".bed"), "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (j in seq_len(n_snps)) {
    bits <- encode(geno_raw[, j])
    if (pad_n > 0L) bits <- c(bits, rep(0L, pad_n))
    m         <- matrix(bits, nrow=4L)
    byte_vals <- m[1,] + m[2,] * 4L + m[3,] * 16L + m[4,] * 64L
    writeBin(as.raw(byte_vals), con)
  }
  close(con)
  cat("PLINK files written to", out_prefix, "\n")
}

write_plink_binary_v2(G, snp_chrom, snp_pos, "sim_geno/sim_content")

# Saving everything into the ground truth file

cat("Saving ground truth...\n")
saveRDS(list(
  params              = params,
  all_effects         = all_effects,
  all_expr_train      = all_expr_train,
  all_expr_test       = all_expr_test,
  train_donors        = train_donors,
  test_donors         = test_donors,
  context_donors      = context_donors,
  G_train             = G_train,
  G_test              = G_test,
  gene_ids            = valid_genes,
  contexts            = contexts,
  snp_ids             = snp_ids,
  snp_chrom           = snp_chrom,
  snp_pos             = snp_pos,
  real_donor_counts   = real_donor_counts
), "sim_ground_truth.rds")

cat("\n=== SIMULATION COMPLETE ===\n")


