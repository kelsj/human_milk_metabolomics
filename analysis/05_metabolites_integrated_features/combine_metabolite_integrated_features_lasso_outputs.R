suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

indir <- "analysis/05_metabolites_integrated_features/output/by_metabolite"

files <- list.files(indir, pattern = "\\.rds$", full.names = TRUE)

res <- lapply(files, readRDS)

fit_list <- lapply(res, function(x) if (!is.null(x$fit)) x$fit else NULL)
fit_list <- fit_list[!vapply(fit_list, is.null, logical(1))]
met.fit.ss <- bind_rows(fit_list) %>% as_tibble()

ve_list <- lapply(res, function(x) if (!is.null(x$feature_ve)) x$feature_ve else NULL)
ve_list <- ve_list[!vapply(ve_list, is.null, logical(1))]
met.ve.ss <- bind_rows(ve_list) %>% as_tibble() %>%
  mutate_at(vars(2:4),as.numeric)

oof_list <- lapply(res, function(x) if (!is.null(x$oof)) x$oof else NULL)
oof_list <- oof_list[!vapply(oof_list, is.null, logical(1))]
met.oof.fit <- bind_rows(oof_list) %>% as_tibble()

saveRDS(met.fit.ss,
        file = "analysis/05_metabolites_integrated_features/output/metabolite_integrated_linRegModelFit.rds")

saveRDS(met.ve.ss,
        file = "analysis/05_metabolites_integrated_features/output/metabolite_integrated_stabSelFeat_featureVE.rds")

saveRDS(met.oof.fit,
        file = "analysis/05_metabolites_integrated_features/output/metabolite_integrated_oof_rsq.rds")
