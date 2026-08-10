### KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome

source("utils/milk_metabolomics_env.R")

# get HEI values from metadata
metadat = readRDS("data/metabolite_wgcna_metadata.rds")
hei = metadat %>% select(randID, all_of(grep("^HEI",colnames(metadat),value=T)))

# load covar residualized metabolomics data
met.resid = readRDS("data/milk_metabolomics_covarResiduals.rds")

# rank norm HEI values
# combine w/ metabolomics data
# corr b/t HEI2015C3_TOTALFRUIT and HEI2015C4_WHOLEFRUIT = 0.8253173
# -> remove HEI2015C3_TOTALFRUIT, keep only HEI2015C4_WHOLEFRUIT
hei.met = hei %>% select(-HEI2015C3_TOTALFRUIT) %>% na.omit() %>%
  mutate_at(vars(2:14),RankNorm) %>% left_join(met.resid, by="randID")

# pairwise correlations
hei.met.corr = rcorr(as.matrix(hei.met %>% select(-randID)), type="pearson")
hei.met.corr.tbl = as_tibble(hei.met.corr$r)
hei.met.corr.tbl$t1 = rownames(hei.met.corr$r)
hei.met.corr.tbl = hei.met.corr.tbl %>% gather(key=t2, value=r, 1:ncol(hei.met.corr$r))
hei.met.corr.tblp = as_tibble(hei.met.corr$P)
hei.met.corr.tblp$t1 = rownames(hei.met.corr$r)
hei.met.corr.tblp = hei.met.corr.tblp %>% gather(key=t2, value=p, 1:ncol(hei.met.corr$r))

hei.met.corr.tbl = hei.met.corr.tbl %>% left_join(hei.met.corr.tblp) %>%
  filter(grepl("^HEI",t1)==T, grepl("^HEI",t2)==F) %>% arrange(p)

hei.met.corr.tbl$q = qvalue(hei.met.corr.tbl$p)$qvalue
hei.met.corr.tbl %>% filter(q<0.05)

# get hierarchical clustering order for HEI categories & milk metabolites
metab.hc.order = hc_order(as.matrix(met.resid %>% select(-randID)),method="pearson")
hei.hc.order = hc_order(as.matrix(hei %>% select(-randID,-HEI2015C3_TOTALFRUIT)),
                        method="pearson")
hei.met.corr.tbl$t1 = factor(hei.met.corr.tbl$t1, levels = hei.hc.order)
hei.met.corr.tbl$t2 = factor(hei.met.corr.tbl$t2, levels = metab.hc.order)

saveRDS(hei.met.corr.tbl %>% rename(HEI_scale = t1, metabolite=t2),
        file="analysis/02_metabolites_vs_diet/output/hei_metab_pwise_corrs.rds")

