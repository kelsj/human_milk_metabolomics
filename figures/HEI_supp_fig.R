## KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome
## supp fig: HEI score var expl. by HEI ~ metabolites lasso

source("utils/milk_metabolomics_env.R")

# load out-of-fold r2
hei.fit.oof = readRDS("analysis/02_metabolites_vs_diet/HEI_metabolites_lasso/output/HEI_metabolite_oof_rsq.rds")
oof.summary = hei.fit.oof %>% group_by(metabolite) %>% 
  slice_max(oof_rsq, n=1, with_ties=F) %>% ungroup() %>% arrange(-oof_rsq)

oof.summary$metabolite = factor(oof.summary$metabolite, levels = oof.summary$metabolite)
median.rsq = oof.summary %>% pull(oof_rsq) %>% median()

ggplot(oof.summary) +
  geom_bar(aes(x=metabolite,y=oof_rsq,fill=n_features),stat="identity",width=1) + 
  geom_hline(yintercept=median.rsq, linetype="dotted")+ 
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1)) + 
  xlab("HEI component") + 
  ylab("Variance explained") + labs(fill="# metabolites") +
  scale_fill_gradient(low="#BDE8F5",high="#005461")
ggsave(filename = "figures/HEI_supp_fig.png",width=6,height=4, dpi=600)


