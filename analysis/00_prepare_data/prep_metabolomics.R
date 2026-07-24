source("utils/milk_metabolomics_env.R")

## normalize metabolite values and regress out covariates

# metabs loaded in R env

# load metadata
metadat = read.csv("data/sample_metadata.csv") %>% as_tibble()

# log normalize metabolomics data
metabs.lognorm = metabs %>% mutate_at(vars(2:ncol(metabs)), log) %>%
  mutate_at(vars(2:ncol(metabs)), scale) %>% mutate_at(vars(2:ncol(metabs)), as.numeric)

# n=349 inds with all needed metadata
# study site, m1_infantage, storage time, m1_milkvol
met.mat.resid = metabs %>% select(randID) %>% 
  left_join(metadat %>% select(randID, enrollment_site, m1_infantage, storage_years, m1_milkvol)) %>% 
  na.omit() %>% select(randID)

for (m in colnames(metabs %>% select(-randID))){
  # 342 have all metadata cols
  mydat = metabs.lognorm %>% select(randID, all_of(m)) %>% rename(mymet=!!m) %>%
    left_join(metadat %>% select(randID, enrollment_site, m1_infantage, storage_years, 
                                       m1_milkvol)) %>% na.omit() %>%
    mutate_at(vars(m1_infantage, storage_years, m1_milkvol), scale) %>%
    mutate_at(vars(m1_infantage, storage_years, m1_milkvol), as.numeric)
  mydat$resid = lm(mydat$mymet ~ ., data = mydat %>% select(-randID, -mymet))$residuals
  colnames(mydat)[ncol(mydat)] = m
  met.mat.resid = met.mat.resid %>% left_join(mydat[,c(1,ncol(mydat))])
}

saveRDS(met.mat.resid, "data/milk_metabolomics_covarResiduals.rds")

