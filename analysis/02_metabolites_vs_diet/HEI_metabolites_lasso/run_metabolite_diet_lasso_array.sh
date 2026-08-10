#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/project2/kelseyj3_1714/kelseyj3/milk_metabolomics/human_milk_metabolomics"
cd "${PROJECT_DIR}"

R_SCRIPT="${PROJECT_DIR}/analysis/02_metabolites_vs_diet/HEI_metabolites_lasso/metabolite_diet_lasso_slurm.R"
OUTDIR="${PROJECT_DIR}/analysis/02_metabolites_vs_diet/HEI_metabolites_lasso/output/by_metabolite"
LOGDIR="${PROJECT_DIR}/analysis/02_metabolites_vs_diet/HEI_metabolites_lasso/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

N_HEI=14

sbatch \
  --job-name=met_gene_lasso \
  --array=1-"${N_HEI}" \
  --cpus-per-task=1 \
  --mem=8G \
  --time=08:00:00 \
  --chdir="${PROJECT_DIR}" \
  --output="${LOGDIR}/hei_met_lasso_%A_%a.out" \
  --error="${LOGDIR}/hei_met_lasso_%A_%a.err" \
  --export=ALL \
  --wrap="Rscript '${R_SCRIPT}' \"\$SLURM_ARRAY_TASK_ID\" '${OUTDIR}' 10"
