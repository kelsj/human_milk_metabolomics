## KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome
## Figure 3

source("utils/milk_metabolomics_env.R")

## 3A: histogram of variance exlained by genes expression selected by lasso reg

# load var expl per metabolite output
# met.fit = readRDS("~/Documents/postdoc/milk_genomics/milk_metabolome_paper_2024/metab_vs_geneExpr/metabolite_vs_geneExpr_lasso/stabSel_v2/metabolite_geneExpr_stabSelFeat_linRegModelFit.rds")
met.fit = readRDS("analysis/03_metabolites_vs_gene_expr/output/metabolite_geneExpr_stabSelFeat_linRegModelFit.rds")

met.fit$metabolite = factor(met.fit$metabolite, levels = met.fit %>% arrange(-rsq) %>% 
                              pull(metabolite))
median.rsq = met.fit %>% pull(rsq) %>% median()

ggplot(met.fit %>% filter(q<=0.05)) +
  geom_bar(aes(x=metabolite,y=rsq,fill=n_features),stat="identity",width=1) + 
  geom_hline(yintercept=median.rsq, linetype="dotted")+ 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank()) + 
  xlab("Metabolite") + 
  ylab("Variance explained") + labs(fill="# features") +
  scale_fill_gradient(low="#BDE8F5",high="#4477AA")
ggsave("figures/Fig3A.png",width=7,height=2)


## 3B-C: load gene expression & metabolite data
# set up data
tpm.1mo.resid = readRDS("data/gene_expression_lasso_inputDat.rds")

met.resid = readRDS("data/milk_metabolomics_covarResiduals.rds")
met.resid.matrix = as.matrix(met.resid %>% select(-randID))
rownames(met.resid.matrix) = met.resid$randID
met.resid.mat = met.resid.mat[colnames(tpm.1mo.resid),]

tpm.1mo.resid = t(tpm.1mo.resid)

## 3B: QPRT (ENSG00000103485) expression vs quinolate in milk
plot.dat = as_tibble(cbind(tpm.1mo.resid[,"ENSG00000103485"],met.resid.mat[,"QUINOLINATE"]))
ggplot(plot.dat, aes(x=V1, y=V2)) +
  geom_smooth(method="lm",se=T,color="black",alpha=0.2,linetype="dotted",linewidth=0.5) +
  geom_point(color="#4477AA",size=1) +
  xlab("QPRT expression") + ylab("Quinolinate")
ggsave("figures/Fig3B.png",width=3,height=3)

## 3C: LALBA (ENSG00000167531) expression vs palmitoleic acid in milk
plot.dat = as_tibble(cbind(tpm.1mo.resid[,"ENSG00000167531"],met.resid.mat[,"PALMITOLINOLEATE"]))
ggplot(plot.dat, aes(x=V1, y=V2)) +
  geom_smooth(method="lm",se=T,color="black",alpha=0.2,linetype="dotted",linewidth=0.5) +
  geom_point(color="#4477AA",size=1) +
  xlab("LALBA expression") + ylab("Palmitoleic acid")
ggsave("figures/Fig3C.png",width=3,height=3)

## 3D: DPEP1 (ENSG00000015413) eQTL -> cysteine-glycine



