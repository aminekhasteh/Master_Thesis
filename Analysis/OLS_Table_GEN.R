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
library(ggrepel) # for overlapping labels in ggplot

# Reading ROS/MAP phenotype dataset
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/ROSMAP_Phenotype/ROSmaster.rds")
# Reading Filtered PNUKBB manifest
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated.csv")
# Reading PCA of the ROS/MAP Genotype dataset
geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function will generate a table for OLS

assoc_table_GEN <- function(MHC,
                            APOE,
                            cum=TRUE,
                            PCnum=12,
                            label_subset = levels(as.factor(meta_pheno$category_manual)),
                            trait_type = levels(as.factor(meta_pheno$trait_type))){
                if(isTRUE(cum)){
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_MHC/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/Associativity/Cum/No_MHC/"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/Associativity/Cum/No_APOE/"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_MHC_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/Associativity/Cum/No_MHC_APOE/"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/With_MHC_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/Associativity/Cum/With_MHC_APOE/"
                                }
                } else {
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_MHC/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/Associativity/Non_cum/No_MHC/"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/Associativity/Non_cum/No_APOE/"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_MHC_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/Associativity/Non_cum/No_MHC_APOE/"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/With_MHC_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/Associativity/Non_cum/With_MHC_APOE/"
                                }
                }
                
                results.G <- readRDS(paste0(path,"PCA_results_all_p-vals",".rds"))
                for (p in names(results.G)){
                                print(p)
                                respca <- results.G[[p]]$pc
                                res.var <- get_pca_var(respca)
                                
                                if (dim(summary(respca)$importance)[2] <= PCnum){
                                                N <- dim(summary(respca)$importance)[2]
                                } else {
                                                N <- PCnum
                                }
                                
                                for (pcnum in 1:N){
                                                var <-sum(summary(respca)$importance[2,pcnum])*100 
                                                ### get contributing PRS for given PC
                                                pc_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
                                                pc_top20$phenotypes <- rownames(pc_top20)
                                                names(pc_top20) <- c("Cos2","Contribution","phenotypes")
                                                # represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
                                                
                                                # contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
                                                g <- ggplot(pc_top20, aes(x=reorder(phenotypes, Contribution), y=Contribution)) +
                                                                geom_bar(stat="identity",position="dodge")+
                                                                theme_minimal()+
                                                                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                                                                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum, 
                                                                               ' (',var,'% variation explained',')'))+
                                                                coord_flip()
                                                ggsave(paste0(path_to_save,"top_",PCnum,"_contribution_at_", p ,"_threshold_PC_",pcnum,".jpg"), width = 30, height = 37, units = "cm")
                                }
                }
}