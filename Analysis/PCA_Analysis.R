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
# This function generates the info we need for PCA results:
pca_results_GEN <- function(MHC,
                            APOE,
                            cum=TRUE,
                            label_subset = levels(as.factor(meta_pheno$category_manual)),
                            trait_type = levels(as.factor(meta_pheno$trait_type))){
                
                if(isTRUE(cum)){
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Cum/No_MHC/"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_APOE_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Cum/No_APOE/"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC_APOE_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Cum/No_MHC_APOE/"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Cum/With_MHC_APOE/"
                                }
                } else {
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Non_Cum/No_MHC/"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_APOE_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Non_Cum/No_APOE/"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC_APOE_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Non_Cum/No_MHC_APOE/"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Non_Cum/With_MHC_APOE/"
                                }
                }
                
                # read in residuals of PRSs files
                lmatrix <- list()
                prs_filenames <- list.files(path,full.names = T)
                for (prsfile in prs_filenames[which(grepl(".rds", prs_filenames))]) {
                                lmatrix[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
                }
                
                pca_results <- matrix(ncol = 7)
                for (prs_matrix in names(lmatrix)){
                                print(prs_matrix)
                                matrix <- lmatrix[[prs_matrix]]$residuals[,-dim(lmatrix[[prs_matrix]]$residuals)[2]]
                                # including manual categories-----------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
                                # including trait type-----------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
                                #--------------------------------------------------------------------------------------------------------------------------
                                
                                respca <- prcomp(matrix)
                                pc1_cont <- sum(summary(respca)$importance[2,1])*100
                                pc2_cont <- sum(summary(respca)$importance[2,2])*100
                                top_5_cont <- sum(summary(respca)$importance[2,1:5])*100
                                top_10_cont <- sum(summary(respca)$importance[2,1:10])*100
                                for(r in 5:dim(summary(respca)$importance)[2]){
                                                if(sum(summary(respca)$importance[2,1:r])>0.80){
                                                                top_r_80_explained <- r
                                                                break
                                                }
                                }
                                row <- c(prs_matrix,pc1_cont,pc2_cont,top_5_cont,top_10_cont,r,dim(summary(respca)$importance)[2])
                                pca_results <- rbind(pca_results,row)
                }
                
                colnames(pca_results) <- c("P-value_thresh","PC1_var_explained(%)","PC2_var_explained(%)","PC1-PC5_var_explained(%)",
                                           "PC1-PC10_var_explained(%)","top_r_comps_explain_80%","total_components")
                pca_results <- pca_results[2:dim(pca_results)[1],]
                pca_results <- as.data.frame(pca_results)
                write.csv(pca_results,paste0(path_to_save,'pca_summary.csv'))
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function will find the p-values of each PC for each selected phenotype in ROSMAP, and generate plots for them
pca_association <- function(MHC,
                            APOE,
                            cum=TRUE,
                            PCnum=12,
                            label_subset = levels(as.factor(meta_pheno$category_manual)),
                            trait_type = levels(as.factor(meta_pheno$trait_type))){
                if(isTRUE(cum)){
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_MHC/"
                                                path_to_save_plot = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Cum/No_MHC/"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_APOE_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_APOE/"
                                                path_to_save_plot = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Cum/No_APOE/"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC_APOE_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_MHC_APOE/"
                                                path_to_save_plot = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Cum/No_MHC_APOE/"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/With_MHC_APOE/"
                                                path_to_save_plot = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Cum/With_MHC_APOE/"
                                }
                } else {
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_MHC/"
                                                path_to_save_plot = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Non_cum/No_MHC/"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_APOE_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_APOE/"
                                                path_to_save_plot = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Non_cum/No_APOE/"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC_APOE_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_MHC_APOE/"
                                                path_to_save_plot = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Non_cum/No_MHC_APOE/"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_Resid/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/With_MHC_APOE/"
                                                path_to_save_plot = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Non_cum/With_MHC_APOE/"
                                }
                }
                
                # read in residuals of PRSs files
                lmatrix <- list()
                prs_filenames <- list.files(path,full.names = T)
                for (prsfile in prs_filenames[which(grepl(".rds", prs_filenames))]) {
                                lmatrix[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
                }
                
                results.G <- list()
                for (pc in names(lmatrix)){
                                print(pc)
                                matrix <- lmatrix[[pc]]$residuals[,-dim(lmatrix[[pc]]$residuals)[2]]
                                # including manual categories-----------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
                                # including trait type-----------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
                                #--------------------------------------------------------------------------------------------------------------------------
                                
                                respca <- prcomp(matrix)
                                res.pcs <- as.data.frame(respca$x)
                                res.pcs$IID <- lmatrix[[pc]]$residuals$IID
                                pheno_pcs <- merge(res.pcs,ROSmaster,by="IID")
                                pheno_pcs2 <- merge(pheno_pcs,geno_pcs,by="IID")
                                ## identify relationships of PCs with phenotypes of interest
                                pheno.list <- c("amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt",
                                                "neuroticism_12","anxiety_20items","cesdsum_lv","educ","mf3123","it3123","vm3123","pput3123",
                                                "stroke_ever","thyroid_ever","headinjrloc_ever","diabetes_sr_rx_ever","cancer_ever","hypertension_ever")
                                
                                PClist <- paste0("PC",seq(1:PCnum))
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
                                
                                if((pc == 'RESIDresults_p_val_1.rds')|
                                   (pc == "RESIDresults_p_val_1e-05.rds")|
                                   (pc == "RESIDresults_p_val_1e-06.rds")|
                                   (pc == "RESIDresults_p_val_1e-07.rds")|
                                   (pc == "RESIDresults_p_val_5e-08.rds")|
                                   (pc == "RESIDresults_p_val_1e-08.rds")|
                                   (pc == "RESIDresults_p_val_1e-09.rds")|
                                   (pc == "RESIDresults_p_val_5e-05.rds")|
                                   (pc == "RESIDresults_p_val_5e-06.rds")|
                                   (pc == "RESIDresults_p_val_5e-07.rds")|
                                   (pc == "RESIDresults_p_val_5e-09.rds")){
                                                p_val <- paste0(unlist(strsplit(pc, "[_.]"))[4])
                                                assocres$colour <- ifelse(assocres$b < 0, "Negative effect","Positive effect")
                                                g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                                                                geom_bar(stat="identity",position="dodge",aes(fill = colour))+
                                                                geom_hline(aes(yintercept = -log10(0.05/nrow(assocres)),color="P-value < 0.00022"),lty=2)+ 
                                                                geom_hline(aes(yintercept = -log10(0.05),color="P-value < 0.05"),lty=2)+ 
                                                                scale_linetype_manual(name = "limit", values = c(2, 2), 
                                                                                      guide = guide_legend(override.aes = list(color = c("green", "orange"))))+
                                                                geom_text_repel(data=subset(assocres,p<0.05),aes(label=pc))+
                                                                theme_minimal()+
                                                                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                                                                ggtitle(paste0('PRSs at P-value threshhold of',"\n",
                                                                               p_val))
                                                ggsave(paste0(path_to_save_plot,"Associations_top_",PCnum,"_PCs_at_", p_val,".jpg"), 
                                                       width = 54, height = 27, units = "cm")
                                } else {
                                                p_val <- paste0(unlist(strsplit(pc, "[_.]"))[4], '.',
                                                               unlist(strsplit(pc, "[_.]"))[5])
                                                assocres$colour <- ifelse(assocres$b < 0, "Negative effect","Positive effect")
                                                g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                                                                geom_bar(stat="identity",position="dodge",aes(fill = colour))+
                                                                geom_hline(aes(yintercept = -log10(0.05/nrow(assocres)),color="P-value < 0.00022"),lty=2)+ 
                                                                geom_hline(aes(yintercept = -log10(0.05),color="P-value < 0.05"),lty=2)+ 
                                                                scale_linetype_manual(name = "limit", values = c(2, 2), 
                                                                                      guide = guide_legend(override.aes = list(color = c("green", "orange"))))+
                                                                geom_text_repel(data=subset(assocres,p<0.05),aes(label=pc))+
                                                                theme_minimal()+
                                                                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                                                                ggtitle(paste0('PRSs at P-value threshhold of',"\n",
                                                                               p_val))
                                                ggsave(paste0(path_to_save_plot,"Associations_top_",PCnum,"_PCs_at_", p_val,".jpg"), 
                                                       width = 54, height = 27, units = "cm")
                                                
                                }
                                
                                results.G[[p_val]] <- list(pc= respca,
                                                                resultsols=assocres,
                                                                plot=g)
                }
                saveRDS(results.G,file=paste0(path_to_save,"PCA_results_all_p-vals",".rds"))
                #return(results.G)
                
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function will plot contributions of top n PCs for each selected PRS matrix
pca_contribution <- function(MHC,
                             APOE,
                             cum=TRUE,
                             PCnum=12,
                             label_subset = levels(as.factor(meta_pheno$category_manual)),
                             trait_type = levels(as.factor(meta_pheno$trait_type))){
                if(isTRUE(cum)){
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_MHC/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Cum/No_MHC/"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Cum/No_APOE/"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/No_MHC_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Cum/No_MHC_APOE/"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PCA/With_MHC_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Cum/With_MHC_APOE/"
                                }
                } else {
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_MHC/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Non_cum/No_MHC/"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Non_cum/No_APOE/"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/No_MHC_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Non_cum/No_MHC_APOE/"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PCA/With_MHC_APOE/"
                                                path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Plots/PCA/Non_cum/With_MHC_APOE/"
                                }
                }
                
                results.G <- readRDS(paste0(path,"PCA_results_all_p-vals",".rds"))
                for (p in names(results.G)){
                                print(p)
                                respca <- results.G[[p]]$pc
                                res.var <- get_pca_var(respca)
                                for (pcnum in 1:PCnum){
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
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################################################################################################################################################
pca_results_GEN(MHC = TRUE, APOE = TRUE)
pca_association(MHC = TRUE, APOE = TRUE, PCnum = 20)
pca_contribution(MHC = TRUE, APOE = TRUE, PCnum = 20)
