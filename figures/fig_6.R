## KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome
## Figure 6

source("utils/milk_metabolomics_env.R")
mypal = c("#3B9AB2", "#EBCC2A", "#F21A00")

## 6A: histogram of variance explained by integrated features selected by lasso reg

# load out-of-fold r2
met.fit.oof = readRDS("analysis/05_metabolites_integrated_features/output/metabolite_integrated_oof_rsq.rds")
oof.summary = met.fit.oof %>% group_by(metabolite) %>% 
  slice_max(oof_rsq, n=1, with_ties=F) %>% ungroup() %>% arrange(-oof_rsq)

# filter for metabolites with r2>1% (179 metabolites)
keep.mets = oof.summary %>% filter(oof_rsq>0.02) %>% pull(metabolite) %>% as.character()

# load feature-level r2
feat.ve = readRDS("analysis/05_metabolites_integrated_features/output/metabolite_integrated_stabSelFeat_featureVE.rds")
# load feature types
dhq.feat = colnames(readRDS("data/lasso_input_dhq_features.rds"))
mqtl.feat = colnames(readRDS("data/lasso_input_combined_mqtl_genotypes.rds") %>% select(-randID))
# classify feature type & scale partial r2 by oof r2
feat.ve = feat.ve %>% mutate(type=ifelse(feature %in% dhq.feat,"Maternal diet",
                               ifelse(feature %in% mqtl.feat, "mQTL","Other")))
# scale by oof r2 and filter metabolites
feat.ve.type = feat.ve %>% group_by(metabolite,type) %>%
  reframe(type_r2 = sum(feat_r2)) %>% left_join(oof.summary %>% select(metabolite,oof_rsq)) %>%
  group_by(metabolite) %>% mutate(type_r2_prop = type_r2/sum(type_r2), type_r2_sc = type_r2_prop*oof_rsq) %>%
  ungroup() %>% filter(metabolite %in% keep.mets)

# arrange metabolites by majority type, then by total rsq
pct.by.type.top.sort = feat.ve.type %>% group_by(metabolite) %>%
  slice_max(type_r2_sc, n=1, with_ties=F) %>% ungroup() %>% 
  arrange(type, -oof_rsq)
feat.ve.type$metabolite = factor(feat.ve.type$metabolite, levels = pct.by.type.top.sort$metabolite)

f6a = ggplot(feat.ve.type) +
  geom_bar(aes(x=metabolite,y=type_r2_sc,fill=type),position="stack",stat="identity",width=1) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank()) + 
  xlab("Metabolite") + 
  ylab("Variance explained") + labs(fill="Feature type") +
  scale_fill_manual(values=mypal)
ggsave(f6a, file="figures/Fig6A.png",width=7,height=2,dpi=600)


## 6B: var expl for metabolites with all 3 feature types selected

# metabolites with all 3 feature types selected
all3plot = feat.ve.type %>% group_by(metabolite) %>% filter(n()==3) %>% ungroup() %>%
  left_join(met.annot %>% select(metabolite,other_name)) %>%
  mutate(other_name = ifelse(metabolite=="2-NONENOYLCARNITINE","2-Nonenoylcarnitine",other_name))

f6b = ggplot(all3plot) +
  geom_bar(aes(x=other_name,y=type_r2_sc,fill=type),stat="identity",position="dodge",width=0.5) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.ticks.x = element_blank(),panel.grid = element_blank(),
        legend.position="none") + 
  xlab("Metabolite") + 
  ylab("Variance explained") + labs(fill="Feature type") +
  scale_fill_manual(values=mypal) + coord_flip()

plot_grid(plotlist=list(f6a + theme(legend.position="none"), f6b), nrow=1, rel_widths = c(2,1.5))
ggsave(file="figures/Fig6AB.png",width=10,height=2.5,dpi=600)



## Fig 6C: 2-nonenoylcarnitine vs. its top features of each type

# scale partial r2 by oof_rsq
# get top feature from each feature type
feat.2nonenoyl = feat.ve %>% filter(metabolite=="2-NONENOYLCARNITINE") %>% 
  left_join(oof.summary %>% select(metabolite,oof_rsq)) %>%
  mutate(feat_r2_prop = feat_r2/sum(feat_r2), feat_r2_sc = feat_r2_prop*oof_rsq) %>%
  group_by(type) %>%
  slice_max(feat_r2,n=1) %>% ungroup()

# get input feature & metabolomics data
dhq.pt.dat = readRDS("data/lasso_input_combined_dhq_metadat.rds")
mqtl.gt.dat = readRDS("data/lasso_input_combined_mqtl_genotypes.rds")
metab.dat = readRDS("data/lasso_input_combined_metabolites_genoPCresid.rds")

dat.2nonenoyl = metab.dat %>% 
  left_join(dhq.pt.dat) %>% left_join(mqtl.gt.dat) %>% 
  select(randID, `2-NONENOYLCARNITINE`, all_of(feat.2nonenoyl$feature))

f6ca = ggplot(dat.2nonenoyl, aes(x=sfa_22_0_behenic_acid_g_ndsr, y=`2-NONENOYLCARNITINE`)) +
  geom_smooth(color="gray",linetype="dotted",alpha=0.2,method="lm") +
  geom_point(color=mypal[1],size=0.5) +
  xlab("sfa_22_0_behenic_acid_g_ndsr") + ylab("2-Nonenoylcarnitine")
    
f6cb = ggplot(dat.2nonenoyl %>% filter(!is.na(cesarean)), aes(x=as.factor(cesarean), y=`2-NONENOYLCARNITINE`)) +
  geom_beeswarm(color=mypal[3],size=0.5,cex=2) +
  geom_boxplot(color=mypal[3],fill=NA,outlier.shape=NA) +
  xlab("cesarean") + ylab("2-Nonenoylcarnitine")

f6cc = ggplot(dat.2nonenoyl %>% filter(!is.na(rs10804178)), 
              aes(x=as.factor(rs10804178), y=`2-NONENOYLCARNITINE`)) +
  geom_beeswarm(color=mypal[2],size=0.5,cex=2) +
  geom_boxplot(color=mypal[2],fill=NA,outlier.shape=NA) +
  xlab("rs10804178") + ylab("2-Nonenoylcarnitine")

# contributing foods to behenic acid
dhq3.db = read_excel("data/dhq3.vers0516.nutdb.xlsx")

f = "sfa_22_0_behenic_acid_g_ndsr"
feat.dhq.name = grep(f,colnames(dhq3.db),ignore.case=T,value=T)
feat.foods = dhq3.db %>% filter(Sex==1, Size==1) %>%
  select(Name, `Age Group`, Sex, Size, all_of(feat.dhq.name)) %>% 
  rename(feat=!!feat.dhq.name) %>% mutate(feat=as.numeric(feat)) %>%
  filter(feat>0) %>% arrange(-feat) %>% head(n=10)
feat.foods$Name = factor(feat.foods$Name, levels = rev(feat.foods$Name))
f6cd = ggplot(feat.foods %>% head(n=10)) +
  geom_bar(aes(y=Name, x=feat), stat="identity", fill="#B0C5A4") +
  xlab(f) + ylab("Food") +
  theme(panel.grid=element_blank()) +
  scale_x_continuous(expand=c(0,0))

f6c_all = plot_grid(plotlist = list(f6ca,f6cb,f6cc,f6cd), nrow=1, rel_widths = c(1,1,1,1.5))
ggsave(f6c_all,file="figures/fig6C.png",width=10,height=2.5,dpi=600)
