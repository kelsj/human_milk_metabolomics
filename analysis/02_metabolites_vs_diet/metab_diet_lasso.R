### KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome

source("utils/milk_metabolomics_env.R")

### metabolite ~ diet feature lasso 

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
saveRDS(dhq.pt.mat, file="analysis/02_metabolites_vs_diet/lasso_input_dhq_features.rds")
saveRDS(met.resid.matrix, file="analysis/02_metabolites_vs_diet/lasso_input_metabolites.rds")


## lasso feature selection for each metabolite ##

# loop through each metabolite
met.fit = data.frame()
met.ve = data.frame()

# for (m in 1:ncol(met.resid.matrix)){
for (m in 1:10){
    # m = 282 # KYNURENINE
  mymet = colnames(met.resid.matrix)[m]
  print(m)
  print(mymet)

  ## lasso variable selection for diet features
  # run 10 replicates of 10-fold cross validation
  nreps=10
  r.corr = vector(length=nreps)
  r.nsigfeat = vector(length=nreps)
  r.feat = list()
  for (r in 1:nreps){
    # print(r)
    set.seed(r)
    cv=glmnet::cv.glmnet(dhq.pt.mat, met.resid.matrix[,m], alpha=1, nfolds=10, 
                         type.measure='mse', keep=T)
    nfeat = length(coef(cv, s='lambda.min'))
    nsigfeat = length(which(as.vector(coef(cv, s='lambda.min')[2:nfeat])!=0))
    r.nsigfeat[r] = nsigfeat
    if (nsigfeat>0){
      # calculate R2 with real data
      r.corr[r] = cor(met.resid.matrix[,m], cv$fit.pre[,which(cv[["lambda"]]==cv$lambda.min)])
      r.feat[[r]] = colnames(dhq.pt.mat)[which(as.vector(coef(cv, s='lambda.min')[2:nfeat])!=0)]
    }
  }
  # if all replicates have 0 features selected OR no positive correlation between predicted & real metabolites
  # skip this metabolite
  if (sum(r.nsigfeat)==0 | sum(r.corr>0)==0){
    print(paste(mymet,"no model"))
    next
  }
  # save features from best replicate (by correlation)
  best.rep = which(r.corr==max(r.corr))[1] # just take first in case of tie
  best.feat = r.feat[[best.rep]]

  # save best.feat model fit
  tmp=data.frame(dhq.pt.mat[,best.feat])
  colnames(tmp) = best.feat
  all.fit=lm(met.resid.matrix[,m]~.,tmp)
  fitsumm = summary(all.fit)
  save.fit=c(nsigfeat,fitsumm$r.squared,fitsumm$adj.r.squared,
             pf(fitsumm$fstatistic[1],fitsumm$fstatistic[2],fitsumm$fstatistic[3], 
                lower.tail = F),mymet)
  met.fit = rbind(met.fit, save.fit)
  ## double check below
  
  # feature partial rsq
  met.mod.coef = summary(all.fit)$coefficients
  met.rsq.partial = rsq.partial(all.fit)
  feat.ve = as_tibble(cbind(met.rsq.partial$variable, met.rsq.partial$partial.rsq,
                            met.mod.coef[2:nrow(met.mod.coef),1], met.mod.coef[2:nrow(met.mod.coef),4])) %>%
    mutate(metabolite=mymet)
  met.ve = rbind(met.ve, feat.ve)
}

colnames(met.fit) = c("n_features","rsq","adj.rsq","F_pval","metabolite")
met.fit = met.fit %>% as_tibble() %>% mutate_at(vars(1:4),as.numeric)
met.fit$q = p.adjust(met.fit$F_pval, method='BH')
saveRDS(met.fit, file="analysis/02_metabolites_vs_diet/metabolite_diet_lassoFeat_linRegModelFit.rds")

colnames(met.ve) = c("feature","feat_r2","feat_est","feat_p","metabolite")
met.ve = as_tibble(met.ve) %>% mutate_at(vars(2:4),as.numeric)
saveRDS(met.ve, file="analysis/02_metabolites_vs_diet/metabolite_diet_lassoFeat_featureVE.rds")


