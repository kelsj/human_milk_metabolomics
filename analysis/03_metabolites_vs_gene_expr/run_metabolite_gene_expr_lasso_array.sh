#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${PROJECT_DIR}/metabolite_gene_expr_lasso_slurm.R"
OUTDIR="${PROJECT_DIR}/analysis/03_metabolites_vs_gene_expr/output/by_metabolite"
LOGDIR="${PROJECT_DIR}/analysis/03_metabolites_vs_gene_expr/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

N_METABOLITES="$(Rscript -e 'x <- readRDS("data/milk_metabolomics_covarResiduals.rds"); cat(ncol(dplyr::select(x, -randID)))')"

sbatch \
  --job-name=met_gene_lasso \
  --array=1-"${N_METABOLITES}" \
  --cpus-per-task=1 \
  --mem=8G \
  --time=08:00:00 \
  --output="${LOGDIR}/met_gene_lasso_%A_%a.out" \
  --error="${LOGDIR}/met_gene_lasso_%A_%a.err" \
  --export=ALL \
  --wrap="Rscript '${R_SCRIPT}' '${SLURM_ARRAY_TASK_ID}' '${OUTDIR}' 10"
