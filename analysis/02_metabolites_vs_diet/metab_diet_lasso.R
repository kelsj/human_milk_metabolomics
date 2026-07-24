### KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome

source("utils/milk_metabolomics_env.R")


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


