### KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome

source("utils/milk_metabolomics_env.R")

### metabolite ~ gene expression lasso 

## prepare data ##

# load covar residualized metabolomics data
met.resid = readRDS("data/milk_metabolomics_covarResiduals.rds")

# make matrix of metabolite residuals
met.resid.matrix = as.matrix(met.resid %>% select(-randID))
rownames(met.resid.matrix) = met.resid$randID

# load gene expression data
tpm.1mo.resid = readRDS("data/gene_expression_lasso_inputDat.rds")

# same inds in each matrix
met.resid.mat = as.matrix(met.resid %>% select(-randID))
rownames(met.resid.mat) = met.resid$randID

met.resid.mat = met.resid.mat[colnames(tpm.1mo.resid),]

tpm.1mo.resid = t(tpm.1mo.resid)

samp.size = nrow(tpm.1mo.resid) # 299

## run lasso for each metabolite
met.fit.ss = data.frame()
met.ve.ss = data.frame()
met.oof.fit = data.frame()

# run 10 replicates of 10-fold cross validation for each metabolite
nreps=10

for (i in 1:ncol(met.resid.mat)){
  # i = 282 # KYNURENINE
  m = colnames(met.resid.mat)[i]
  print(m)

  ## lasso variable selection
  r.corr = vector(length=nreps) # vector of out-of-fold prediction r2 values for each rep
  r.nsigfeat = vector(length=nreps) # vector of number selected features for each rep
  r.feat = list() # list of selected features for each rep
  for (r in 1:nreps){
    print(r)
    set.seed(r)
    cv=glmnet::cv.glmnet(tpm.1mo.resid, met.resid.mat[,m], alpha=1, nfolds=10, type.measure='mse', keep=T)
    nfeat = length(coef(cv, s='lambda.min'))
    nsigfeat = length(which(as.vector(coef(cv, s='lambda.min')[2:nfeat])!=0))
    r.nsigfeat[r] = nsigfeat
    if (nsigfeat>0){
      # calculate R2 b/t out-of-fold predictions and real data
      r.corr[r] = summary(lm(met.resid.mat[,m] ~ cv$fit.pre[,which(cv[["lambda"]]==cv$lambda.min)]))$r.squared
      # store selected features
      r.feat[[r]] = colnames(tpm.1mo.resid)[which(as.vector(coef(cv, s='lambda.min')[2:nfeat])!=0)]
    }
  }
  # if all replicates have 0 features selected
  # or median oof rsq == 0
  # skip this metabolite
  if (sum(r.nsigfeat)==0 | median(r.corr)==0){
    print(paste(m,"no model"))
    next
  }
  # save oof prediction r2
  met.oof.fit = rbind(met.oof.fit, cbind(r.nsigfeat, r.corr, 1:nreps, rep(m,nreps)))

  # save features selected in >half of replicates
  # only include reps with >0 features selected
  r.feat.stab = r.feat[r.nsigfeat>0]
  stab.feat = names(which(table(unlist(r.feat.stab))>(length(r.feat.stab)/2)))

  # save selected feature model fit
  tmp=data.frame(tpm.1mo.resid[,stab.feat])
  colnames(tmp) = stab.feat
  all.fit=lm(met.resid.mat[,m]~.,tmp)
  fitsumm = summary(all.fit)
  save.fit=c(length(stab.feat),fitsumm$r.squared,fitsumm$adj.r.squared,
             pf(fitsumm$fstatistic[1],fitsumm$fstatistic[2],fitsumm$fstatistic[3], lower.tail = F),m)
  met.fit.ss = rbind(met.fit.ss, save.fit)
  # get partial r2 for selected features
  met.mod.coef = summary(all.fit)$coefficients
  met.rsq.partial = rsq.partial(all.fit)
  feat.ve = as_tibble(cbind(met.rsq.partial$variable, met.rsq.partial$partial.rsq, 
                            met.mod.coef[2:nrow(met.mod.coef),1], met.mod.coef[2:nrow(met.mod.coef),4])) %>%
    mutate(metabolite=m)
  met.ve.ss = rbind(met.ve.ss, feat.ve)
  
}

colnames(met.fit.ss) = c("n_features","rsq","adj.rsq","F_pval","metabolite")
met.fit.ss = met.fit.ss %>% as_tibble() %>% mutate_at(vars(1:4),as.numeric)
met.fit.ss$q = p.adjust(met.fit.ss$F_pval, method='BH')
saveRDS(met.fit.ss, file="analysis/03_metabolites_vs_gene_expr/output/metabolite_geneExpr_stabSelFeat_linRegModelFit.rds")

colnames(met.ve.ss) = c("feature","feat_r2","feat_est","feat_p","metabolite")
met.ve.ss = as_tibble(met.ve.ss) %>% mutate_at(vars(2:4),as.numeric)
saveRDS(met.ve.ss, file="analysis/03_metabolites_vs_gene_expr/output/metabolite_geneExpr_stabSelFeat_featureVE.rds")

colnames(met.oof.fit) = c("n_features","oof_rsq","rep","metabolite")
met.oof.fit = as_tibble(met.oof.fit) %>% mutate_at(vars(1:2),as.numeric)
saveRDS(met.oof.fit, file="analysis/03_metabolites_vs_gene_expr/output/metabolite_geneExpr_oof_rsq.rds")


