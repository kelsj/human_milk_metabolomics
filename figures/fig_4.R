## KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome
## Figure 4

source("utils/milk_metabolomics_env.R")

## load mQTL summary statistics (filtered for p-value < 1e-5)
res = readRDS("data/gwas_metabs_gemma_maf10_p1e-5.rds")

## bonferroni threshold: 5e-8 / 458 = 1.091703e-10
bonf.thresh = 5e-8/458

## Fig 4A: manhattan plot of all metabolites. filter for MAF>10%
met.res = res
met.res.cum <- met.res %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(ps)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)

met.res <- met.res %>% 
  inner_join(met.res.cum, by = "chr") %>% 
  mutate(bp_cum = ps + bp_add)

axis_set <- met.res %>% 
  group_by(chr) %>% 
  summarise(center = mean(bp_cum))

# plot top loci
ggplot(met.res) + 
  geom_hline(yintercept = -log10(bonf.thresh),linetype="dashed",color="black") +
  geom_hline(yintercept = -log10(5e-8),linetype="dashed",color="darkgray") +
  geom_point(aes(x=bp_cum, y=-log10(p_wald), color=factor(chr)), size=0.5) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  theme(legend.position="None",panel.grid = element_blank(),
        axis.text.x=element_text(size=rel(0.75))) + 
  xlab("Chromosome") +
  ylab("-log10(P-value)")
ggsave("figures/fig4A.png",width=6.5,height=3,dpi=600)


## 4C: plasma vs. milk mQTL effect estimates
milk_v_plasma_qtl = readRDS("data/milk_vs_chen_plasma_mQTL_effects.rds")

ggplot(milk_v_plasma_qtl) + # %>% filter(power>0.95, p_wald>0.05)) + 
  geom_hline(yintercept=0,linetype="dashed",color="gray") +
  geom_vline(xintercept=0,linetype="dashed",color="gray") +
  geom_point(aes(x=plasma_b, y=beta, color=power),size=0.5) +
  geom_text_repel(aes(x=plasma_b, y=beta, color=power,label=milk_met),size=3,max.overlaps=15) +
  geom_errorbar(aes(x=plasma_b, ymin=beta-1.96*se, ymax=beta+1.96*se, color=power),alpha=0.2,width=0) +
  geom_errorbar(aes(y=beta, xmin=plasma_b-1.96*plasma_se, xmax=plasma_b+1.96*plasma_se, color=power),alpha=0.2,width=0) +
  scale_color_gradient(low="darkblue",high="magenta") +
  xlab("Plasma effect estimate (95% CI)") +
  ylab("Milk effect estimate (95% CI)") +
  labs(color="Replication\npower") +
  theme(panel.grid=element_blank())

ggsave("figures/fig4C.png",width=5,height=4,dpi=600)


## 4D: plasma vs. milk metabolite PGS R^2
pgs_comp = readRDS("data/plasmaPGS_replication_milk_powerEst.rds")

ggplot() + 
  scale_fill_gradient(low="darkblue",high="magenta") +
  geom_point(data=pgs_comp, aes(x=Internal_R2,y=R2, fill=power),shape=21,size=3, alpha=0.4) +
  geom_point(data=pgs_comp %>% filter(p.value<=0.05), 
             aes(x=Internal_R2,y=R2, fill=power),shape=21,size=3, alpha=1) +
  geom_text_repel(data=pgs_comp, aes(x=Internal_R2,y=R2, label=Query),
                  size=rel(2),point.padding = 0, max.overlaps=8) +
  xlab("PGS plasma R-squared") +
  ylab("PGS milk R-squared") + 
  labs(fill="Replication\npower") + 
  scale_x_log10() + scale_y_log10() +
  theme(panel.grid=element_blank())
ggsave("figures/fig4D.png",width=5,height=4,dpi=600)



