source("utils/milk_metabolomics_env.R")

# load unproccessed DHQ data
dhq.obs = readRDS("data/DHQ_diet_data_unprocessed.rds")

# check number of zeroes per feature/tmpt
missing.zeros = dhq.obs %>% group_by(var,redcap_event_name,dhq_type_component) %>% 
  summarise(nObs = n(), fracZero = sum(val==0,na.rm=T)/nObs) %>%
  ungroup()
missing.zeros %>% arrange(-fracZero)
table(missing.zeros$nObs)

# remove features that have >90% zeroes in 1+ timepoint/study version
# don't remove fped_a_drinks_usda [expected to be zero for most for pregnancy timepoints!! ]
rm.vars.zeros = missing.zeros %>% filter(fracZero>0.9, var!="fped_a_drinks_usda") %>% 
  pull(var) %>% unique()
# removes 12 features

# calculate log normalized values for each variable within each tmpt
# replace zeros with (minimum value / 10) to calculate log norm
dhq.obs.sc = dhq.obs %>% filter(!(var %in% rm.vars.zeros)) %>%
  na.omit() %>% 
  group_by(var) %>% 
  mutate(grp.min = min(val[val > 0])) %>% 
  ungroup() %>% 
  mutate(val.imp = ifelse(val==0, grp.min/10, val)) %>%
  group_by(redcap_event_name,var) %>% mutate(val.sc = as.numeric(scale(log(val.imp)))) %>% 
  ungroup()

# get each individual's values averaged across all timepoints
# record # timepoints going into average (nTmpt)
# then normalize w/in each variable
dhq.ind.means = dhq.obs.sc %>% select(randID, dhq_type_component, redcap_event_name, var, val.sc) %>% 
  group_by(randID,dhq_type_component,var) %>% summarise(mean.val = mean(val.sc), nTmpt = n()) %>% 
  ungroup() %>% 
  group_by(var) %>%
  mutate(mean.val = as.numeric(scale(mean.val))) %>% ungroup()

table(dhq.ind.means$nTmpt) # people with 1, 2, or 3 survey timepoints are all present
table(dhq.ind.means$dhq_type_component) # both study versions present

# std normalize within # of timepoints used to calculate DHQ mean 
# AND by DHQ survey version
dhq.ind.means = dhq.ind.means %>% group_by(var, nTmpt, dhq_type_component) %>%
  mutate(mean.val = as.numeric(scale(mean.val))) %>% ungroup()

saveRDS(dhq.ind.means, "data/DHQ_indMeans.rds")


### prepare input data for metabolite ~ diet feature lasso 

## prepare metabolomics data ##

# load covar residualized metabolomics data
met.resid = readRDS("data/milk_metabolomics_covarResiduals.rds")

# make matrix of metabolite residuals
met.resid.matrix = as.matrix(met.resid %>% select(-randID))
rownames(met.resid.matrix) = met.resid$randID

## prepare diet data ## 

# load DHQ data
dhq = readRDS("data/DHQ_indMeans.rds")
# make into matrix
dhq.obs = dhq %>% select(-nTmpt,-dhq_type_component) %>%
  spread(key=var, value=mean.val)
# n=334 inds
dhq.met.inds = dhq.obs %>% filter(randID %in% met.resid$randID)
# make matrix - exclude missing data. n=325
dhq.met.inds = dhq.met.inds %>% na.omit() 
dhq.pt.mat = as.matrix(dhq.met.inds %>% select(-randID))
rownames(dhq.pt.mat) = dhq.met.inds$randID

## filter mets for inds w/ all DHQ data (n=325)
met.resid.matrix = met.resid.matrix[rownames(dhq.pt.mat),]

## get rid of colinear predictors from dhq.pt.mat
dhq.pt.corrs = rcorr(dhq.pt.mat)
dhq.pt.corr.tbl = as_tibble(dhq.pt.corrs$r) 
dhq.pt.corr.tbl$t1 = rownames(dhq.pt.corrs$r)
dhq.pt.corr.tbl = dhq.pt.corr.tbl %>% gather(key=t2, value=r, 1:ncol(dhq.pt.corrs$r)) %>%
  filter(t1<t2)

tight.corr = dhq.pt.corr.tbl %>% filter(r>0.99)

rm.traits = c()
while (nrow(tight.corr)>0){
  rm.t = tight.corr %>% arrange(-r) %>% pull(t1) %>% head(n=1)
  print(rm.t)
  rm.traits = c(rm.traits, rm.t)
  tight.corr = dhq.pt.corr.tbl %>% filter(r>0.99, !(t1 %in% rm.traits), !(t2 %in% rm.traits))
}
length(rm.traits) # removes 25 traits
rm.traits.index = which(colnames(dhq.pt.mat) %in% rm.traits)
dhq.pt.mat = dhq.pt.mat[,-rm.traits.index]
dim(dhq.pt.mat) #[1] 325 174

## rank normalize all diet features & metabolite values
dhq.pt.mat = apply(dhq.pt.mat, 2, RankNorm)
met.resid.matrix = apply(met.resid.matrix, 2, RankNorm)

# save input datasets
saveRDS(dhq.pt.mat, file="data/lasso_input_dhq_features.rds")
saveRDS(met.resid.matrix, file="data/lasso_input_metabolites.rds")


