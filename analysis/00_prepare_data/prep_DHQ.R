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

