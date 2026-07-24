## KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome
## Figure 2

source("utils/milk_metabolomics_env.R")

## 2A: histogram of variance exlained by diet features selected by lasso reg

# load model fit 
# met.fit = readRDS("~/Documents/postdoc/milk_genomics/milk_metabolome_paper_2024/diet_metabolites_March2025/lasso_var_expl_dietOnly_v2/metabolite_diet_lassoFeat_linRegModelFit.rds")
met.fit = readRDS("analysis/02_metabolites_vs_diet/output/metabolite_diet_lassoFeat_linRegModelFit.rds")

met.fit$metabolite = factor(met.fit$metabolite, levels = met.fit %>% arrange(-rsq) %>% pull(metabolite))
median.rsq = met.fit %>% pull(rsq) %>% median()

f2a = ggplot(met.fit %>% filter(q<=0.05)) +
  geom_bar(aes(x=metabolite,y=rsq,fill=n_features),stat="identity",width=1) + 
  geom_hline(yintercept=median.rsq, linetype="dotted")+ 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank()) + 
  xlab("Metabolite") + 
  ylab("Variance explained") + labs(fill="# features") +
  scale_fill_gradient(low="#BDE8F5",high="#005461")
ggsave(f2a, filename = "figures/Fig2A.png",width=7,height=2)


## 2B-2D: examples of metabolites well explained by diet

# load feature variance explained
# met.ve = readRDS("~/Documents/postdoc/milk_genomics/milk_metabolome_paper_2024/diet_metabolites_March2025/lasso_var_expl_dietOnly_v2/diet_lasso_feature_model_stats.rds")
met.ve = readRDS("analysis/02_metabolites_vs_diet/output/metabolite_diet_lassoFeat_featureVE.rds")

# load DHQ3 data
# dhq3.db = read_excel("~/Documents/postdoc/milk_genomics/MILK_data/DHQ/dhq3.vers0516.nutdb.xlsx")
dhq3.db = read_excel("data/dhq3.vers0516.nutdb.xlsx")

# load lasso input data
dhq.pt.mat = readRDS("analysis/02_metabolites_vs_diet/output/lasso_input_dhq_features.rds")
met.resid.matrix = readRDS("analysis/02_metabolites_vs_diet/output/lasso_input_metabolites.rds")

# metabolites to plot
plot.mets = list("B"="QUINATE","C"="PROLINE BETAINE","D"="S1P")

for (p in names(plot.mets)){
  mymet = plot.mets[[p]]
  
  # first plot: var. expl. by top 10 diet features
  met.feat = met.ve %>% filter(metabolite==mymet) %>% arrange(-feat_r2)
  # add direction of effect in mv regression model
  met.feat = met.feat %>% mutate(Direction=ifelse(feat_est<0,"neg.","pos."))
  # if only one feature vs. multiple
  # plot
  met.feat$feature = factor(met.feat$feature, levels = rev(met.feat$feature))
  p2aa = ggplot(met.feat %>% head(10)) + 
    geom_bar(aes(y=feature,x=feat_r2,fill=Direction),stat="identity",width=1) +
    xlab("Partial R-squared") + 
    ylab("Diet feature") +
    scale_fill_manual(values=c("#161E54","#DA498D"))
  
  # second plot: scatter of metabolite vs. top diet feature
  myfeat = met.feat %>% arrange(-feat_r2) %>% head(1) %>% pull(feature) %>% as.character()
  plot.dat = as_tibble(cbind(dhq.pt.mat[,myfeat],met.resid.matrix[,mymet]))
  p2ab = ggplot(plot.dat, aes(x=V1, y=V2)) +
    geom_smooth(method="lm",se=T,color="black",alpha=0.2,linetype="dotted",linewidth=0.5) +
    geom_point(color="#3B9AB2",size=1) +
    xlab(myfeat) + ylab(mymet)
  
  # third plot: top foods contributing to this feature
  feat.dhq.name = grep(myfeat,colnames(dhq3.db),ignore.case=T,value=T)
  if (length(feat.dhq.name)==0){
    print(paste0("no match ",myfeat))
    next
  }
  feat.foods = dhq3.db %>% filter(Sex==1, Size==1) %>%
    select(Name, `Age Group`, Sex, Size, all_of(feat.dhq.name)) %>% 
    rename(feat=!!feat.dhq.name) %>% mutate(feat=as.numeric(feat)) %>%
    filter(feat>0) %>% arrange(-feat) %>% head(n=10)
  feat.foods$Name = factor(feat.foods$Name, levels = rev(feat.foods$Name))
  p2ac = ggplot(feat.foods %>% head(n=10)) +
    geom_bar(aes(y=Name, x=feat), stat="identity", fill="#B0C5A4") +
    xlab(myfeat) + ylab("Food")

  # put all 3 together
  plot_grid(plotlist=list(p2aa,p2ab,p2ac), nrow=1, rel_widths=c(4,2,3), align="h")
  ggsave(paste0("figures/Fig2",p,".png"), width=10, height=2)
}


## 2E: pairwise correlations b/t HEI & metabolites
hei.metab.corrs = readRDS("analysis/02_metabolites_vs_diet/output/hei_metab_pwise_corrs.rds")

ggplot(hei.metab.corrs) + 
  geom_tile(aes(x=metabolite, y=HEI_scale, fill=r)) +
  xlab("Milk metabolite") + ylab("HEI component") +
  scale_fill_gradient2(low="#060771",mid="#f7f7f7",high="#E67E22",midpoint = 0) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid=element_blank())
ggsave("figures/Fig2E.png",width=7,height=2)



