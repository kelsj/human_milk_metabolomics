### R env for milk metabolomics paper ###

# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(stringr)
library(Hmisc)
library(ggcorrplot)
library(ggrepel)
library(readxl)
library(glmnet)
library(RNOmni)
library(qvalue)
library(cowplot)
library(glmnet)
library(rsq)
theme_set(theme_bw())

## hierarchical clustering function of matrix
hc_order = function(X, method){
  # 1. Compute correlation matrix
  cor_mat <- cor(X, method = method, use = "complete.obs")
  # 2. Convert correlation to distance (1 - correlation)
  dist_mat <- as.dist(1 - cor_mat)
  # 3. Hierarchical clustering
  hc <- hclust(dist_mat, method = "average")
  return(colnames(X)[hc$order])
}

## load metabolite annotations
met.annot = readRDS("data/metabolite_annotations.rds")

## load unprocessed metabolomics data
metabs = readRDS("data/milk_metabolomics.rds")

## load metabolites color scheme - organized by metabolite sub class
grp.assign = readRDS("utils/metabolite_subclass_color_assignments.rds")
color.code = readRDS("utils/metabolite_subclass_color_palette.rds")



