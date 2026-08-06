## KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome
## Figure 1

source("utils/milk_metabolomics_env.R")

## 1C: metabolite mean abundance vs. coefficient of variation

# calculate mean & c.o.v.
met.summary = metabs %>% gather(key=metabolite, value=abund, 2:ncol(metabs)) %>%
  group_by(metabolite) %>% reframe(metmean=mean(abund),metsd=sd(abund),cv=metsd/metmean) %>% 
  left_join(met.annot %>% select(metabolite, other_name, Sub_Class))

# if metab sub class not in color code, rename as "other"
named.levels = as.character(color.code$level[1:(nrow(color.code)-1)])
met.summary = met.summary %>% mutate(Sub_Class=ifelse(Sub_Class %in% named.levels, Sub_Class, "Other"))
# order sub classes to match color palette
met.summary$Sub_Class = factor(met.summary$Sub_Class, levels = color.code$level)

# plot mean abund vs. c.o.v.
f1c = ggplot(met.summary) + 
  geom_point(aes(x=metmean,y=cv,fill=Sub_Class),color="black",shape=21,size=2.5) + 
  geom_text_repel(aes(x=metmean,y=cv,label=other_name),max.overlaps = 12, color="black",size=3) + 
  scale_x_log10() + scale_y_log10() +
  scale_fill_manual(values=color.code$colorcode) +
  xlab("Metabolite mean abundance") +
  ylab("Metabolite coef. of var.") + labs(fill="") +
  theme(legend.position="none", panel.grid=element_blank())
ggsave(f1c, filename="figures/Fig1C.png",width=4,height=4,dpi=600)



## Fig 1D. plot metabolite modules from wgcna

# load module assignments
net = readRDS("analysis/01_metabolite_wgcna/output/metabolite_modules_network.rds")

# count of metabolites per subclass/module
# remove module 0 (unassigned module) + add sub class assignments from met.annot
metab.modules = as_tibble(cbind(names(net$colors), net$colors)) %>% 
  rename(metabolite=V1, module = V2) %>%
  filter(module!="0") %>% 
  left_join(met.annot %>% select(metabolite, other_name, Sub_Class))

# if metab sub class not in color code, rename as "Other"
named.levels = as.character(color.code$level[1:(nrow(color.code)-1)])
metab.modules = metab.modules %>% 
  mutate(Sub_Class=ifelse(Sub_Class %in% named.levels, Sub_Class, "Other"))

# if not in color code, rename as "Other"
module.cat.counts = metab.modules  %>% 
  group_by(module,Sub_Class) %>% reframe(nMets=n())

# order levels by color code
module.cat.counts$Sub_Class = factor(module.cat.counts$Sub_Class, levels = color.code$level)
keep.levels = as.character(unique(module.cat.counts$Sub_Class))
color.code = color.code[which(color.code$level %in% keep.levels),]

f1d = ggplot(module.cat.counts) + 
  geom_bar(aes(x=factor(module), y=nMets, fill=Sub_Class), stat="identity", position="stack") +
  scale_fill_manual(values=color.code$colorcode) + 
  xlab("Metabolite module") + ylab("Number of metabolites") +
  labs(fill="") +
  theme(panel.grid=element_blank())
ggsave(f1d, filename="figures/Fig1D.png",width=6,height=4,dpi=600)


## Fig 1E. metabolite subclass overrepresentation in metabolite modules
fisher.res = readRDS("analysis/01_metabolite_wgcna/output/metabolite_module_subclass_overrep_res.rds")

f1e = ggplot(fisher.res) + 
  geom_point(aes(x=module, y=Sub_Class, fill=Sub_Class, size = -log10(q), color=sig), shape=22) +
  scale_fill_manual(values=rev(color.code$colorcode)) + 
  scale_color_manual(values=c("white","black")) +
  xlab("Metabolite module") + ylab("Metabolite subclass") +
  scale_x_discrete(labels=1:8) +
  guides(fill="none", color="none") +
  scale_size(range=c(1,6)) +
  theme(axis.ticks.y=element_blank())
ggsave(f1e, filename="figures/Fig1E.png",width=5,height=4,dpi=600)


## Fig 1F. correlations between metab modules & metadata
# load pairwise corrs
all.trait.corrs.tbl = readRDS("analysis/01_metabolite_wgcna/output/wgcna_module_trait_corrs.rds")

# arrange traits by type
mat.traits = c("Maternal BMI","Maternal age","Parity","GDM","PPWL (1mo)","OGCT result", 
               "EBF (3mo)","EBF (6mo)")
diet.traits = unique(grep("^HEI",all.trait.corrs.tbl$trait, value=T))
all.trait.corrs.tbl = all.trait.corrs.tbl %>% 
  mutate(type = ifelse(trait %in% mat.traits,"Other",
                       ifelse(trait %in% diet.traits,"Diet (HEI)","Milk components")))
all.trait.corrs.tbl$type = factor(all.trait.corrs.tbl$type, 
                                  levels = c("Diet (HEI)","Other","Milk components"))

# only include factors that have at least one sig assoc
sig.traits = all.trait.corrs.tbl %>% filter(q<0.05) %>% pull(trait) %>% unique()

f1f = ggplot(all.trait.corrs.tbl %>% filter(trait %in% sig.traits)) + 
  geom_point(aes(y=trait, x=module, fill=r, color=sig, size=-log10(p.edit)), shape=21) +
  scale_fill_gradient2(low="#060771",mid="#f7f7f7",high="#E67E22",midpoint = 0) + 
  xlab("Metabolite module") + ylab("Traits") +
  labs(size="-log10(p-value)") +
  theme(axis.text.y = element_text(hjust=1,vjust=0.5),
        axis.title.y = element_blank(),
        legend.box = "vertical", panel.grid = element_blank()) +
  scale_x_discrete(labels=1:7) +
  scale_color_manual(values=c("white","black")) +
  scale_size_continuous(range=c(1,6),breaks=seq(3,14,3)) +
  guides(fill = guide_legend(order = 1), # Assign an order to the legends
         size = guide_legend(order = 2), color="none") +
  facet_grid(rows = vars(type), scale="free_y", space="free_y")
ggsave(f1f, filename="figures/fig1F.png",width=5,height=4,dpi=600)

## metabolite subclass legend
subclass.legend = get_legend(f1d + theme(legend.key.size = unit(0.4, "cm"),
                                         legend.text = element_text(size=rel(0.8))))
plot_grid(subclass.legend)
ggsave(filename="figures/fig1_legend.png",width=4,height=4,dpi=600)






