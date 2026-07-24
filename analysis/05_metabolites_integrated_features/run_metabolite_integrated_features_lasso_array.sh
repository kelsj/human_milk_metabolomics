#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/project2/kelseyj3_1714/kelseyj3/milk_metabolomics/human_milk_metabolomics"
cd "${PROJECT_DIR}"

R_SCRIPT="${PROJECT_DIR}/analysis/05_metabolites_integrated_features/metabolite_integrated_features_lasso_slurm.R"
OUTDIR="${PROJECT_DIR}/analysis/05_metabolites_integrated_features/output/by_metabolite"
LOGDIR="${PROJECT_DIR}/analysis/05_metabolites_integrated_features/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

N_METABOLITES="$(Rscript -e 'x <- readRDS("data/lasso_input_combined_metabolites_genoPCresid.rds"); cat(ncol(dplyr::select(x, -randID)))')"

sbatch \
  --job-name=met_gene_lasso \
  --array=1-"${N_METABOLITES}" \
  --cpus-per-task=1 \
  --mem=8G \
  --time=08:00:00 \
  --chdir="${PROJECT_DIR}" \
  --output="${LOGDIR}/met_gene_lasso_%A_%a.out" \
  --error="${LOGDIR}/met_gene_lasso_%A_%a.err" \
  --export=ALL \
  --wrap="Rscript '${R_SCRIPT}' \"\$SLURM_ARRAY_TASK_ID\" '${OUTDIR}' 10"
