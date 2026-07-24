## KE Johnson et al. 2026. Maternal diet and genetics shape the human milk metabolome

## WGCNA analysis of milk metabolites

source("utils/milk_metabolomics_env.R")

library(WGCNA)
library(flashClust)

options(stringsAsFactors = FALSE)

## load metabolomics data

# metabolomics data with tech covars regressed out: data/milk_metabolomics_covarResiduals.rds
met = readRDS("data/milk_metabolomics_covarResiduals.rds")

#make into matrix with metabolites as cols
log.met.resid = as.matrix(met %>% select(-randID))
rownames(log.met.resid) = met$randID

## format data for WGCNA

#check for missing values
gsg = goodSamplesGenes(log.met.resid, verbose = 3)
gsg$allOK # TRUE

#cluster samples & check for outliers
sampleTree = hclust(dist(log.met.resid), method = "average")

#Determine clusters below the line - removes 2 samples
clust = cutreeStatic(sampleTree, cutHeight = 50, minSize = 10)
table(clust)

#remove sample outliers
keepSamples = (clust==1)
rmSamples = rownames(log.met.resid)[clust==0] # [1] "23828" "60040"

log.met.resid = log.met.resid[keepSamples, ]

nGenes = ncol(log.met.resid) # 458
nSamples = nrow(log.met.resid) # 347

#choose soft thresholding power
power.choice = 4

# construct network
set.seed(54321)
# net = blockwiseModules(log.met.resid, power = power.choice,
#                        TOMType = "unsigned", minModuleSize = 10, maxBlockSize = 20000,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "milkMetabNet",
#                        verbose = 3)
# 
# # save network
# saveRDS(net, file="analysis/01_metabolite_wgcna/output/metabolite_modules_network.rds")

net = readRDS("analysis/01_metabolite_wgcna/output/metabolite_modules_network.rds")

## test overrepresentation of metabolite subclasses in metabolite modules
# fisher's exact test

# remove module 0 (unassigned module) + add sub class assignments from met.annot
metab.modules = as_tibble(cbind(names(net$colors), net$colors)) %>% 
  rename(metabolite=V1, module = V2) %>%
  filter(module!="0") %>% 
  left_join(met.annot %>% select(metabolite, other_name, Sub_Class))

# if metab sub class not in color code, rename as "other"
named.levels = as.character(color.code$level[1:(nrow(color.code)-1)])
metab.modules = metab.modules %>% 
  mutate(Sub_Class=ifelse(Sub_Class %in% named.levels, Sub_Class, "Other"))

# get counts of metabolites by subclass
# if sub class not in metabolite color code, rename as "other"
module.cat.counts = metab.modules  %>% 
  group_by(module,Sub_Class) %>% reframe(nMets=n())

fisher.res = data.frame()
for (m in unique(module.cat.counts$module)){
  for (sc in unique(module.cat.counts$Sub_Class)){
    sc.in = sum(module.cat.counts %>% filter(module==m, Sub_Class==sc) %>% pull(nMets))
    sc.out = sum(module.cat.counts %>% filter(module!=m, Sub_Class==sc) %>% pull(nMets))
    other.in = sum(module.cat.counts %>% filter(module==m, Sub_Class!=sc) %>% pull(nMets))
    other.out = sum(module.cat.counts %>% filter(module!=m, Sub_Class!=sc) %>% pull(nMets))
    mod.sc.mat = matrix(c(sc.in,other.in,sc.out,other.out),nrow=2, 
                        dimnames = list(Sub_Class = c("T","F"), Module = c("T","F")))
    test.res = fisher.test(mod.sc.mat, alternative="greater")
    fisher.res = rbind(fisher.res, c(test.res$estimate, test.res$p.value, m, sc))
  }
}
colnames(fisher.res) = c("OR","p","module","Sub_Class")
fisher.res = as_tibble(fisher.res) %>% mutate_at(vars(1:2),as.numeric) %>% arrange(p)

fisher.res$q = p.adjust(fisher.res$p, method="BH")
fisher.res$Sub_Class = factor(fisher.res$Sub_Class, levels = rev(color.code$level))

# replace inf OR w/ max value
max.OR = max(fisher.res %>% filter(!is.infinite(OR)) %>% pull(OR))
fisher.res = fisher.res %>% mutate(OR = ifelse(is.infinite(OR), max.OR, OR),
                                   sig = ifelse(q<0.05,T,F))

# save results
saveRDS(fisher.res, file="analysis/01_metabolite_wgcna/output/metabolite_module_subclass_overrep_res.rds")


## correlations between modules eigengenes & cohort metadata

# load metadata
metadat = readRDS("data/metabolite_wgcna_metadata.rds")

# metabolite module eigenvalues
module.vals = net$MEs %>% as_tibble()
module.vals$randID = rownames(net$MEs)

# combine module vals with metadata
comb.all = module.vals %>% left_join(metadat) %>% select(randID, everything())

# pairwise correlations
trait.corrs = rcorr(as.matrix(comb.all %>% select(-randID)))

module.names = colnames(net$MEs)
trait.names = colnames(comb.all %>% select(-randID, -all_of(module.names)))

trait.corrs.tbl = as_tibble(trait.corrs$r[module.names,trait.names])
trait.corrs.tbl$module = module.names
trait.corrs.tbl = trait.corrs.tbl %>% gather(key=trait, value=r, 1:(ncol(trait.corrs.tbl)-1))
trait.corrs.tblP = as_tibble(trait.corrs$P[module.names,trait.names])
trait.corrs.tblP$module = module.names
trait.corrs.tblP = trait.corrs.tblP %>% gather(key=trait, value=p, 1:(ncol(trait.corrs.tblP)-1))

## remove ME0 (the 'leftover' module)
trait.corrs.tbl = trait.corrs.tbl %>% left_join(trait.corrs.tblP) %>% 
  filter(module!="ME0") %>% arrange(p)

## need to edit p/q-values of 0
min.nonz.p = trait.corrs.tbl %>% filter(p>0) %>% pull(p) %>% min()

trait.corrs.tbl = trait.corrs.tbl %>% mutate(p.edit = ifelse(p==0, min.nonz.p/2, p))
trait.corrs.tbl$q = qvalue(trait.corrs.tbl$p.edit)$qvalue

trait.corrs.tbl$trait = factor(trait.corrs.tbl$trait, levels = trait.names)
trait.corrs.tbl$module = factor(trait.corrs.tbl$module, levels = paste0("ME",1:7))

trait.corrs.tbl = trait.corrs.tbl %>%
  mutate(sig = ifelse(q<0.05, T, F))

saveRDS(trait.corrs.tbl, file="analysis/01_metabolite_wgcna/output/wgcna_module_trait_corrs.rds")





