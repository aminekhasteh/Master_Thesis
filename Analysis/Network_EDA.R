# Importing libraries we need:
library(data.table)
library(openxlsx)
library(rms)
library(ggplot2)
library(factoextra)
library(diptest)
library(corrplot)
library(dendextend)
library(colorspace) # get nice colors
library(harrietr) # dist to long format
library(igraph)
library(ape) # For the dendogram

library("ggplot2")
library("ggdendro")
library("reshape2")

# Reading datasets:
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated.csv")

geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_Resid_PCA_Clust/",full.names = T)

# read in PRS files
results1 <- list()
for (prsfile in prs_filenames[1:10]) {
                results1[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
}

#############################################################################################
###################################   |                  |   ################################                                    
###################################   | Network Analysis |   ################################
###################################   |                  |   ################################
#############################################################################################

#graph_df <- graph_from_data_frame(dist_df_passed, directed = FALSE)
#kc <- cluster_fast_greedy(graph_df)
#length(kc)
#sizes(kc)

summary(dist_df$dist)
sum(dist_df$dist<1)
sum(duplicated(signif(dist_df$dist,digits=11)))

check <- dist_df[which(duplicated(signif(dist_df$dist,digits=8))),]
dist_df1 <- dist_df[which(!duplicated(signif(dist_df$dist,digits=8))),]



graph_df <- simplify(graph_df)
g_degree <- degree(graph_df)
which(g_degree>1172)
hist(g_degree)
meta_pheno[which(meta_pheno$phenocode_annotate_lst == "CA_Summary_Operations_h0.20346_n812"),]$full_description


