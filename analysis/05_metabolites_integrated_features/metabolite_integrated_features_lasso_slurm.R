#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(glmnet)
  library(rsq)
})

source("utils/milk_metabolomics_env.R")

# -----------------------------------------------------------------------------
# Run metabolite ~ integrated features lasso for a single metabolite.
#
# Usage:
#   Rscript metabolite_integrated_features_lasso_slurm.R <met_idx> <outdir> <nreps>
#
# If <met_idx> is not provided, the script falls back to SLURM_ARRAY_TASK_ID.
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

met_idx <- NA_integer_
outdir <- "analysis/05_metabolites_integrated_features/output/by_metabolite"
nreps <- 10L

if (length(args) >= 1 && nzchar(args[1])) {
  met_idx <- suppressWarnings(as.integer(args[1]))
}
if (length(args) >= 2 && nzchar(args[2])) {
  outdir <- args[2]
}
if (length(args) >= 3 && nzchar(args[3])) {
  nreps <- suppressWarnings(as.integer(args[3]))
}

if (is.na(met_idx) || met_idx < 1L) {
  met_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_)))
}

if (is.na(met_idx) || met_idx < 1L) {
  stop("No metabolite index provided. Use the first command-line argument or SLURM_ARRAY_TASK_ID.")
}
if (is.na(nreps) || nreps < 1L) {
  stop("nreps must be a positive integer.")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Prepare data
# -----------------------------------------------------------------------------
met.resid <- readRDS("data/lasso_input_combined_metabolites_genoPCresid.rds")
met.resid.mat <- as.matrix(met.resid %>% select(-randID))
rownames(met.resid.mat) <- met.resid$randID

# combined DHQ & maternal features
dhq.pt.mat <- readRDS("data/lasso_input_combined_dhq_metadat.rds")
# genotype data
mqtl.gt.dat <- readRDS("data/lasso_input_combined_mqtl_genotypes.rds")
# list of mqtl lead variants by metabolite
mqtl.met.pairs <- readRDS("data/lasso_input_combined_mqtl_variant_list.rds")

# create feature matrix for this metabolite - depending on if it has mqtl or not, include mqtl lead var GTs
mqtl.match = mqtl.met.pairs %>% filter(metabolite==m)
if (nrow(mqtl.match)==0){
  # feature matrix w/o genotypes
  feat.tbl = dhq.pt.mat %>% na.omit() %>% mutate_at(vars(b_matage,b_mat_bmi,b_parity),RankNorm)
} else if (nrow(mqtl.match)>0){
  # feature matrix w/ genotypes
  feat.tbl = dhq.pt.mat %>% left_join(mqtl.gt.dat %>% select(randID, all_of(mqtl.match$rs))) %>%
    na.omit() %>% mutate_at(vars(b_matage,b_mat_bmi,b_parity),RankNorm)
}

# create matrix for input features
feat.mat = as.matrix(feat.tbl %>% select(-randID))
rownames(feat.mat) = feat.tbl$randID

# filter met.resid.mat to same inds w/ complete data
met.resid.mat = met.resid.mat[rownames(feat.mat),]

samp.size = nrow(met.resid.mat)

if (met_idx > ncol(met.resid.mat)) {
  stop(sprintf("met_idx=%d exceeds number of metabolites=%d", met_idx, ncol(met.resid.mat)))
}

m <- colnames(met.resid.mat)[met_idx]
message("Running metabolite: ", m, " (index ", met_idx, ")")

y <- met.resid.mat[, m]

# -----------------------------------------------------------------------------
# Run 10 replicates of 10-fold CV lasso
# -----------------------------------------------------------------------------
r.corr <- rep(NA_real_, nreps)
r.nsigfeat <- rep(NA_integer_, nreps)
r.feat <- vector("list", nreps)

lambda_to_col <- function(cvfit) {
  which.min(abs(cvfit$lambda - cvfit$lambda.min))
}

for (r in seq_len(nreps)) {
  message("  replicate ", r, "/", nreps)
  set.seed(r)

  cv <- glmnet::cv.glmnet(
    x = feat.mat,
    y = y,
    alpha = 1,
    nfolds = 10,
    type.measure = "mse",
    keep = TRUE
  )

  coef_min <- coef(cv, s = "lambda.min")
  nfeat <- nrow(coef_min)
  selected_idx <- which(as.vector(coef_min[2:nfeat, , drop = TRUE]) != 0)
  nsigfeat <- length(selected_idx)
  r.nsigfeat[r] <- nsigfeat

  if (nsigfeat > 0L) {
    pred_col <- lambda_to_col(cv)
    pred_oof <- cv$fit.preval[, pred_col]
    r.corr[r] <- summary(lm(y ~ pred_oof))$r.squared
    r.feat[[r]] <- colnames(feat.mat)[selected_idx]
  }
}

# -----------------------------------------------------------------------------
# Skip metabolites where no stable model exists
# -----------------------------------------------------------------------------
if (sum(r.nsigfeat, na.rm = TRUE) == 0L || isTRUE(all.equal(stats::median(r.corr, na.rm = TRUE), 0))) {
  message(m, " no model")
  out <- list(
    metabolite = m,
    met_idx = met_idx,
    status = "no_model",
    nreps = nreps,
    r.nsigfeat = r.nsigfeat,
    r.corr = r.corr
  )

  safe_name <- paste0(sprintf("%04d_", met_idx), gsub("[^A-Za-z0-9_.-]", "_", m), ".rds")
  saveRDS(out, file = file.path(outdir, safe_name))
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# Stable features and final linear model
# -----------------------------------------------------------------------------
r.feat.stab <- r.feat[r.nsigfeat > 0]
stab.feat <- names(which(table(unlist(r.feat.stab)) > (length(r.feat.stab) / 2)))

if (length(stab.feat) == 0L) {
  message(m, " no stable features after replicate filtering")
  out <- list(
    metabolite = m,
    met_idx = met_idx,
    status = "no_stable_features",
    nreps = nreps,
    r.nsigfeat = r.nsigfeat,
    r.corr = r.corr
  )

  safe_name <- paste0(sprintf("%04d_", met_idx), gsub("[^A-Za-z0-9_.-]", "_", m), ".rds")
  saveRDS(out, file = file.path(outdir, safe_name))
  quit(save = "no", status = 0)
}

tmp <- data.frame(feat.mat[, stab.feat, drop = FALSE])
colnames(tmp) <- stab.feat
all.fit <- lm(y ~ ., data = tmp)
fitsumm <- summary(all.fit)

save.fit <- data.frame(
  n_features = length(stab.feat),
  n_inds = samp.size,
  rsq = fitsumm$r.squared,
  adj.rsq = fitsumm$adj.r.squared,
  F_pval = pf(
    fitsumm$fstatistic[1],
    fitsumm$fstatistic[2],
    fitsumm$fstatistic[3],
    lower.tail = FALSE
  ),
  metabolite = m,
  stringsAsFactors = FALSE
)
save.fit$q <- p.adjust(save.fit$F_pval, method = "BH")

met.mod.coef <- summary(all.fit)$coefficients
met.rsq.partial <- rsq.partial(all.fit)

feat.ve <- as_tibble(cbind(
  feature = met.rsq.partial$variable,
  feat_r2 = met.rsq.partial$partial.rsq,
  feat_est = met.mod.coef[2:nrow(met.mod.coef), 1],
  feat_p = met.mod.coef[2:nrow(met.mod.coef), 4]
)) %>%
  mutate(metabolite = m)

oof.df <- tibble(
  n_features = r.nsigfeat,
  oof_rsq = r.corr,
  rep = seq_len(nreps),
  metabolite = m
)

out <- list(
  metabolite = m,
  met_idx = met_idx,
  status = "ok",
  stable_features = stab.feat,
  fit = save.fit,
  feature_ve = feat.ve,
  oof = oof.df
)

safe_name <- paste0(sprintf("%04d_", met_idx), gsub("[^A-Za-z0-9_.-]", "_", m), ".rds")
saveRDS(out, file = file.path(outdir, safe_name))
message("Saved: ", file.path(outdir, safe_name))
