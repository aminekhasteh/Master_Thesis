# Importing libraries we need:
library(data.table)
library(openxlsx)
library(rms)
library(ggplot2)
library(ggrepel) # for overlapping labels in ggplot
library(corrplot)
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(colorspace) # get nice colors
library(cluster)    # clustering algorithms
library(Hmisc)
library(dplyr) # data manipulation
library(tidyverse)  # data manipulation
library(summarytools)

# Reading datasets:
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated_new.csv")
geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

# -----------------------------# Here we can alternate between the PRS with no MHC and with MHC # ----------------------------------------#                                                                                       #|
                                                                                                                                          #|                                                                                                     
#prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_Resid_PCA_Clust/"        #|
#prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_No_MHC_Resid_PCA_Clust/"  #|
Meta_PRS_No_MHC <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/MetaPRS/Top_6/Meta-PRS_No_MHC.csv")                                                                                                                                         #|  #|
                                                                                                                                          #|
#-----------------------------------------------------------------------------------------------------------------------------------------#

# read in PRS files
# results.S <- list()
# prs_filenames <- list.files(prs_path,full.names = T)
# for (prsfile in prs_filenames) {
#                 results.S[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
# }

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#############################################################################################
################################   |                   |   ##################################                                   
################################   | PCA on Phenotypes |   ##################################
################################   |                   |   ##################################
#############################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This loop will generate a table to better underestand the PCs


# pca_results <- matrix(ncol = 7)
# for (pc in names(results.S)){
#                 respca <- results.S[[pc]]$pca
#                 pc1_cont <- sum(summary(respca)$importance[2,1])*100
#                 pc2_cont <- sum(summary(respca)$importance[2,2])*100
#                 top_5_cont <- sum(summary(respca)$importance[2,1:5])*100
#                 top_10_cont <- sum(summary(respca)$importance[2,1:10])*100
#                 for(r in 5:dim(summary(respca)$importance)[2]){
#                                 if(sum(summary(respca)$importance[2,1:r])>0.80){
#                                                 top_r_80_explained <- r
#                                                 break
#                                 }
#                 }
#                 row <- c(pc,pc1_cont,pc2_cont,top_5_cont,top_10_cont,r,dim(summary(respca)$importance)[2])
#                 pca_results <- rbind(pca_results,row)
# }
# 
# colnames(pca_results) <- c("P-value_thresh","PC1_var_explained(%)","PC2_var_explained(%)","PC1-PC5_var_explained(%)",
#                            "PC1-PC10_var_explained(%)","top_r_comps_explain_80%","total_components")
# pca_results <- pca_results[2:11,]
# pca_results <- as.data.frame(pca_results)
# write.csv(pca_results,'/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Manuscript/pca_results_no_mhc.csv')

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This loop will find the p-values of each PC for each selected phenotype in ROSMAP, and generate plots for them
# results.G <- list()
# for (pc in names(results.S)){
#
#
# 
#                 results.G[[pc]] <- list(resultsols=assocres,
#                                         plot=g)
#                 
# }

respca <- prcomp(Meta_PRS_No_MHC[,1:1100])
res.pcs <- as.data.frame(respca$x)
res.pcs$IID <- Meta_PRS_No_MHC$IID
pheno_pcs <- merge(res.pcs,ROSmaster,by="IID")
pheno_pcs2 <- merge(pheno_pcs,geno_pcs,by="IID")
## identify relationships of PCs with phenotypes of interest
pheno.list <- c("amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt",
                "neuroticism_12","anxiety_20items","cesdsum_lv","educ","mf3123","it3123","vm3123","pput3123")

PClist <- paste0("PC",seq(1:20))
index <- 1
pvalues <- NULL
bvalues <- NULL
nvalues <- NULL
phenovalues <- NULL
pcvalues <- NULL
for (pheno in pheno.list) {
                for (PC in PClist) {
                                form <- formula(paste(pheno,"~",PC))
                                mod <- lm(data=pheno_pcs2,form)
                                pvalues[index] <- anova(mod)[1,5]
                                bvalues[index] <- coef(mod)[2]
                                nvalues[index] <- dim(mod$model)[1]
                                pcvalues[index] <- PC
                                phenovalues[index] <- pheno
                                index <- index + 1
                }
}
assocres <- data.frame(pheno=phenovalues,
                       pc=pcvalues,
                       b=bvalues,
                       p=pvalues,
                       n=nvalues,
                       fdr=p.adjust(pvalues))

g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                geom_bar(stat="identity",position="dodge")+
                geom_hline(yintercept = -log10(0.05),col="red",lty=2)+  #/nrow(assocres)
                geom_text(data=subset(assocres,p<0.05),aes(label=pc))+  #/nrow(assocres)
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                ggtitle('Meta-PRSs')

# if((pc == 'PCAresults_p_val_1.rds')|
#    (pc == "PCAresults_p_val_1e-05.rds")|
#    (pc == "PCAresults_p_val_1e-06.rds")|
#    (pc == "PCAresults_p_val_1e-07.rds")|
#    (pc == "PCAresults_p_val_5e-08.rds")){
#                 
#                 
# } else {
#                 g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
#                                 geom_bar(stat="identity",position="dodge")+
#                                 geom_hline(yintercept = -log10(0.05/nrow(assocres)),col="red",lty=2)+ 
#                                 geom_text(data=subset(assocres,p<0.05/nrow(assocres)),aes(label=pc))+
#                                 theme_minimal()+
#                                 theme(axis.text.x=element_text(angle = -45, hjust = 0))+
#                                 ggtitle(paste0('PRSs at P-value threshhold of',"\n",
#                                                paste0(unlist(strsplit(pc, "[_.]"))[4], '.',
#                                                       unlist(strsplit(pc, "[_.]"))[5])))
# }

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

res.var <- get_pca_var(respca)

pcnum <- 2
sum(summary(respca)$importance[2,pcnum])
### get contributing PRS for given PC
### get contributing PRS for given PC10
pc_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc_top20$phenotypes <- rownames(pc_top20)
names(pc_top20) <- c("Cos2","Contribution","phenotypes")
# represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc_top20, aes(x=reorder(phenotypes, Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum))+
                coord_flip()

# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc_top20, aes(x=reorder(phenotypes, Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum))+
                coord_flip()




### get cos2 for top 5 components

corrplot(results.S$PCAresults_p_val_0_0001.rds$res.var$cos2[1:10,1:5],is.corr = F,tl.cex = 0.3)




pheno_pcs3 <- merge(pheno_pcs2,prs_list$p_val_1e_05.txt,by="IID")
summary(lm(data=pheno_pcs3, cogn_global_random_slope ~ PC9))
summary(lm(data=pheno_pcs3, cogn_global_random_slope ~ icd10_C71_C71_Malignant_neoplasm_of_brain))
