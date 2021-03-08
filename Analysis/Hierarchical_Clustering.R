# if (!requireNamespace("BiocManager", quietly = TRUE))
#                 install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")
# Importing libraries we need:
library(data.table)
library(openxlsx)
library(rms)
library(ggplot2)
library(corrplot)
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(colorspace) # get nice colors
library(cluster)    # clustering algorithms
library(Hmisc)
library(dplyr) # data manipulation
library(tidyverse)  # data manipulation
library(summarytools)
library(RColorBrewer)


# Reading datasets:
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated_new.csv")

# -----------------------------# Here we can alternate between the PRS with no MHC and with MHC # ----------------------------------------#                                                                                       #|
                                                                                                                                          #|                                                                                                     
#prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_Resid_PCA_Clust/"        #|
prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_No_MHC_Resid_PCA_Clust/"  #|
                                                                                                                                          #|
#prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PLINK_Resid_PCA_Clust/"         #|
#prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PLINK_No_MHC_Resid_PCA_Clust/"  #|
                                                                                                                                          #|
#-----------------------------------------------------------------------------------------------------------------------------------------#

# read in PRS files
results.S <- list()
prs_filenames <- list.files(prs_path,full.names = T)
for (prsfile in prs_filenames) {
                results.S[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
}
                

# Let's categorize the phenotypes based on their heritability scale:
names <- strsplit(meta_pheno$phenocode_annotate_lst,'_')  
label_h <- NULL
for(i in 1:length(names)){
                label_h <- append(label_h,as.numeric(strsplit(names[[i]][(length(names[[i]])-1)],'h')[[1]][2]))
}
label <- as.factor(ntile(label_h, 4))
# Add new label to meta_pheno:

meta_pheno$herit_label <- label

#############################################################################################
################################   |                        |   #############################                                    
################################   | Clusters on Phenotypes |   #############################
################################   |                        |   #############################
#############################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function generates all of the hc-plots for us:
plot_hc <- function(lmatrix,
                    scaling_method = "Scaled",
                    hc_method = c("Agglomerative","Divisive"),
                    dis_mat_method=c("Manhattan","Euclidean","Correlation"),
                    label_subset = levels(as.factor(meta_pheno$category_manual)),
                    agg_method = "ward.D2",
                    trait_type = levels(as.factor(meta_pheno$trait_type)),
                    width=1800,
                    height=900){ # lmatrix : list of our PRS matrices
                for(prs_matrix in names(lmatrix)){
                                matrix <- lmatrix[[prs_matrix]]$residuals[,-dim(lmatrix[[prs_matrix]]$residuals)[2]]
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
                                # including trait type-----------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
                                #--------------------------------------------------------------------------------------------------------------------------
                                names <- strsplit(colnames(matrix),'_')
                                # label_trait_type <- NULL
                                label_manual <- NULL
                                # label_h <- NULL
                                # for(i in 1:length(names)){
                                #                 label_trait_type <- append(label_trait_type,names[[i]][1])
                                #                 label_h <- append(label_h,as.numeric(strsplit(names[[i]][(length(names[[i]])-1)],'h')[[1]][2]))
                                # }
                                for(i in names(matrix)){
                                                label_manual <- append(label_manual,meta_pheno$category_manual[which(meta_pheno$new_pheno_annot==i)])
                                }
                                
                                # label <- rev(levels(as.factor(label_trait_type))) # This label didn't work well
                                # label <- rev(levels(as.factor(ntile(label_h, 4)))) # This label didn't work as well!
                                label <- rev(levels(as.factor(label_manual)))
                                
                                if(hc_method=="Agglomerative"){
                                                if(dis_mat_method == "Manhattan"){
                                                                clusters <- hclust(dist(t(matrix),method="manhattan"),method = agg_method)
                                                                dend <- as.dendrogram(clusters)
                                                                # order it the closest we can to the order of the observations:
                                                                dend <- rotate(dend, 1:length(clusters$labels))
                                                } else if (dis_mat_method == "Euclidean"){
                                                                clusters <- hclust(dist(t(matrix),method="euclidean"),method = agg_method)
                                                                dend <- as.dendrogram(clusters)
                                                                # order it the closest we can to the order of the observations:
                                                                dend <- rotate(dend, 1:length(clusters$labels))
                                                } else if (dis_mat_method == "Correlation"){
                                                                cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                                                                clusters <- hclust(as.dist(1-cols.cor),method = agg_method)
                                                                dend <- as.dendrogram(clusters)
                                                                # order it the closest we can to the order of the observations:
                                                                dend <- rotate(dend, 1:length(clusters$labels))
                                                }
                                }
                                
                                if(hc_method=="Divisive"){
                                                if(dis_mat_method == "Manhattan"){
                                                                clusters <- diana(dist(t(matrix),method="manhattan"))
                                                                dend <- as.dendrogram(clusters)
                                                                # order it the closest we can to the order of the observations:
                                                                dend <- rotate(dend, 1:(length(clusters$height)+1))
                                                } else if (dis_mat_method == "Euclidean"){
                                                                clusters <- diana(dist(t(matrix),method="euclidean"))
                                                                dend <- as.dendrogram(clusters)
                                                                # order it the closest we can to the order of the observations:
                                                                dend <- rotate(dend, 1:(length(clusters$height)+1))
                                                } else if (dis_mat_method == "Correlation"){
                                                                cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                                                                clusters <- diana(as.dist(1-cols.cor))
                                                                dend <- as.dendrogram(clusters)
                                                                # order it the closest we can to the order of the observations:
                                                                dend <- rotate(dend, 1:(length(clusters$height)+1))
                                                }
                                }
                                
                                # Color the branches based on the clusters:
                                dend <- color_branches(dend,k= length(label),groupLabels = label) #, k= length(levels(label))
                                # Manually match the labels, as much as possible, to the real classification of the flowers:
                                labels_colors(dend) <-
                                                rainbow_hcl(length(label))[sort_levels_values(
                                                                as.numeric(as.factor(label_manual))[order.dendrogram(dend)]
                                                )]
                                
                                # We shall add the trait type to the labels:
                                labels(dend) <- paste(as.character(label_manual)[order.dendrogram(dend)],
                                                      "(",labels(dend),")", 
                                                      sep = "")
                                
                                # We hang the dendrogram a bit:
                                dend <- hang.dendrogram(dend, hang_height = 0.1)
                                # reduce the size of the labels:
                                dend <- set(dend, "labels_cex", c(0.8,1.1))
                                dend <- set(dend, "leaves_pch", 0.4)
                                # And plot:
                                if((prs_matrix == 'PCAresults_p_val_1.rds')|
                                   (prs_matrix == "PCAresults_p_val_1e-05.rds")|
                                   (prs_matrix == "PCAresults_p_val_1e-06.rds")|
                                   (prs_matrix == "PCAresults_p_val_1e-07.rds")|
                                   (prs_matrix == "PCAresults_p_val_5e-08.rds")){
                                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendrograms/Dendrogram_Pheno_",
                                                            scaling_method,'_',gsub("PCAresults_","",gsub(".rds","",prs_matrix)),'_',hc_method,'_',dis_mat_method,'_',agg_method,".jpg"),
                                                     width = width, height = height)
                                                par(mar=c(5,5,10,80),cex=1,font=3)
                                                plot(dend, cex.main=2,
                                                     main = "",
                                                     horiz =  TRUE,
                                                     nodePar = list(cex = .007)
                                                )
                                                mtext(paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                             dim(clusters[[1]])[1],' individuals',"\n",
                                                             'at P-value threshold of ',
                                                             unlist(strsplit(prs_matrix, "[_.]"))[4]), side = 3, line = 1, cex = 2)
                                                #legend("bottomleft", legend = label, fill = rainbow_hcl(5))
                                                dev.off()
                                                
                                } else {
                                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendrograms/Dendrogram_Pheno_",
                                                            scaling_method,'_',gsub("PCAresults_","",gsub(".rds","",prs_matrix)),'_',hc_method,'_',dis_mat_method,'_',agg_method,".jpg"),
                                                     width = width, height = height)
                                                par(mar=c(5,5,10,80),cex=1,font=3)
                                                plot(dend, cex.main=2,
                                                     main = "",
                                                     horiz =  TRUE,
                                                     nodePar = list(cex = .007)
                                                )
                                                mtext(paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                             dim(clusters[[1]])[1],' individuals',"\n",
                                                             'at P-value threshold of ',
                                                             paste0(unlist(strsplit(prs_matrix, "[_.]"))[4], '.',
                                                                    unlist(strsplit(prs_matrix, "[_.]"))[5])), side = 3, line = 1, cex = 2)
                                                #legend("bottomleft", legend = label, fill = rainbow_hcl(5))
                                                dev.off()
                                }
                }
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function generates all of the heatmaps for us:
plot_heat <- function(lmatrix,
                      scaling_method = "Scaled",
                      hc_method = c("Agglomerative","Divisive"),
                      dis_mat_method=c("Manhattan","Euclidean","Correlation"),
                      label_subset = levels(as.factor(meta_pheno$category_manual)),
                      agg_method = "ward.D2",
                      trait_type = levels(as.factor(meta_pheno$trait_type)),
                      width=1800,
                      height=900){ # lmatrix : list of our PRS matrices
                for(prs_matrix in names(lmatrix)){
                                matrix <- lmatrix[[prs_matrix]]$residuals[,-dim(lmatrix[[prs_matrix]]$residuals)[2]]
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
                                # including trait type-----------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
                                #--------------------------------------------------------------------------------------------------------------------------
                                names <- strsplit(colnames(matrix),'_')
                                # label_trait_type <- NULL
                                label_manual <- NULL
                                for(i in names(matrix)){
                                                label_manual <- append(label_manual,meta_pheno$category_manual[which(meta_pheno$new_pheno_annot==i)])
                                }
                                label <- rev(levels(as.factor(label_manual)))

                                if(hc_method=="Agglomerative"){
                                                if(dis_mat_method == "Manhattan"){
                                                                Rowv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(dist(matrix,method="manhattan"),method = agg_method)),k= 4), hang_height = 0.1)
                                                                Colv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(dist(t(matrix),method="manhattan"),method = agg_method)),k= length(label)), hang_height = 0.1)
                                                } else if (dis_mat_method == "Euclidean"){
                                                                Rowv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(dist(matrix,method="euclidean"),method = agg_method)),k= 4), hang_height = 0.1)
                                                                Colv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(dist(t(matrix),method="euclidean"),method = agg_method)),k= length(label)), hang_height = 0.1)
                                                } else if (dis_mat_method == "Correlation"){
                                                                Rowv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(as.dist(1 - cor(t(matrix), use = "pairwise.complete.obs", method = "pearson")),method=agg_method)),k= 4), hang_height = 0.1)
                                                                Colv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(as.dist(1 - cor(matrix, use = "pairwise.complete.obs", method = "pearson")),method=agg_method)),k= length(label)), hang_height = 0.1)
                                                }
                                }
                                
                                if(hc_method=="Divisive"){
                                                if(dis_mat_method == "Manhattan"){
                                                                Rowv  <- hang.dendrogram(color_branches(as.dendrogram(diana(dist(matrix),method="manhattan")),k= 4), hang_height = 0.1)
                                                                Colv  <- hang.dendrogram(color_branches(as.dendrogram(diana(dist(t(matrix)),method="manhattan")),k= length(label)), hang_height = 0.1)
                                                } else if (dis_mat_method == "Euclidean"){
                                                                Rowv  <- hang.dendrogram(color_branches(as.dendrogram(diana(dist(matrix),method="euclidean")),k= 4), hang_height = 0.1)
                                                                Colv  <- hang.dendrogram(color_branches(as.dendrogram(diana(dist(t(matrix)),method="euclidean")),k= length(label)), hang_height = 0.1)
                                                } else if (dis_mat_method == "Correlation"){
                                                                Rowv  <- hang.dendrogram(color_branches(as.dendrogram(diana(as.dist(1 - cor(t(matrix), use = "pairwise.complete.obs", method = "pearson")))),k= 4), hang_height = 0.1)
                                                                Colv  <- hang.dendrogram(color_branches(as.dendrogram(diana(as.dist(1 - cor(matrix, use = "pairwise.complete.obs", method = "pearson")))),k= length(label)), hang_height = 0.1)
                                                }
                                }
                                # reduce the size of the labels:
                                # dend <- set(dend, "labels_cex", c(0.8,1.1))
                                # dend <- set(dend, "leaves_pch", 0.4)
                                # And plot:
                                if((prs_matrix == 'PCAresults_p_val_1.rds')|
                                   (prs_matrix == "PCAresults_p_val_1e-05.rds")|
                                   (prs_matrix == "PCAresults_p_val_1e-06.rds")|
                                   (prs_matrix == "PCAresults_p_val_1e-07.rds")|
                                   (prs_matrix == "PCAresults_p_val_5e-08.rds")){
                                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Heatmaps/Heatmap_Pheno_",
                                                            scaling_method,'_',gsub("PCAresults_","",gsub(".rds","",prs_matrix)),'_',hc_method,'_',dis_mat_method,'_',agg_method,".jpg"),
                                                     width = width, height = height)
                                                par(mar=c(5,5,5,5),cex=1,font=3)
                                                heatmap(as.matrix(matrix), Rowv = Rowv, Colv = Colv,
                                                        scale = "none",
                                                        col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
                                                legend(x="topleft", legend=c("Inevrsely Correlated", "No correlation", "Correled"), 
                                                       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
                                                mtext(paste0('Heatmap of PRSs at phenotypes level',"\n",
                                                             dim(matrix)[1],' individuals and ',dim(matrix)[2],' phenotypes',"\n",
                                                             'at P-value threshold of ',
                                                             unlist(strsplit(prs_matrix, "[_.]"))[4]), side = 3, line = 1, cex = 2)
                                                dev.off()
                                } else {
                                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Heatmaps/Heatmap_Pheno_",
                                                            scaling_method,'_',gsub("PCAresults_","",gsub(".rds","",prs_matrix)),'_',hc_method,'_',dis_mat_method,'_',agg_method,".jpg"),
                                                     width = width, height = height)
                                                par(mar=c(5,5,5,5),cex=1,font=3)
                                                heatmap(as.matrix(matrix), Rowv = Rowv, Colv = Colv,
                                                        scale = "none",
                                                        col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
                                                legend(x="topleft", legend=c("Inevrsely Correlated", "No correlation", "Correled"), 
                                                       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
                                                mtext(paste0('Heatmap of PRSs for ',"\n",
                                                             dim(matrix)[1],' individuals and ',dim(matrix)[2],' phenotypes',"\n",
                                                             'at P-value threshold of ',
                                                             paste0(unlist(strsplit(prs_matrix, "[_.]"))[4], '.',
                                                                    unlist(strsplit(prs_matrix, "[_.]"))[5])), side = 3, line = 1, cex = 2)
                                                dev.off()
                                }
                }
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


                # ------------------------------------- #
                # Agglomerative Hierarchical Clustering #
                # ------------------------------------- #

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function will create the best performing clustering sructure of the matrices:
agglom_coef <- function(label_subset = levels(as.factor(meta_pheno$category_manual)),
                        file_name,
                        trait_type = levels(as.factor(meta_pheno$trait_type))){
                agglom_result_dat <- matrix(ncol = 7)
                for(prs_matrix in names(results.S)){
                                matrix <- results.S[[prs_matrix]]$residuals[,-dim(results.S[[prs_matrix]]$residuals)[2]]
                                                                
                                # including trait type-----------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
                                #--------------------------------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
                                # methods to assess
                                m <- c( "average", "single", "complete", "ward")
                                names(m) <- c( "average", "single", "complete", "ward")
                                for(dis_mat in c('Manhattan', 'Euclidean', 'Correlation')) {
                                                if (dis_mat == "Manhattan") {
                                                                # function to compute coefficient
                                                                ac <- function(x) {agnes(dist(t(matrix),method = "manhattan"),method = x)$ac}
                                                                row <- map_dbl(m, ac)
                                                }
                                                if (dis_mat == "Euclidean") {
                                                                # function to compute coefficient
                                                                ac <- function(x) {agnes(dist(t(matrix),method = "euclidean"),method = x)$ac}
                                                                row <- map_dbl(m, ac)
                                                }
                                                if (dis_mat == "Correlation") {
                                                                # Pairwise correlation between samples (columns)
                                                                cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                                                                # function to compute coefficient
                                                                ac <- function(x) {agnes(as.dist(1-cols.cor),method = x)$ac}
                                                                row <- map_dbl(m, ac)
                                                }
                                                row <- c(prs_matrix,dis_mat,as.vector(row))
                                                agglom_result_dat <- rbind(agglom_result_dat,row)
                                }

                }
                agglom_result_dat1 <- agglom_result_dat[2:dim(agglom_result_dat)[1],]
                agglom_result_dat1 <- as.data.frame(agglom_result_dat1)
                names(agglom_result_dat1)<-c('P-value threshhold','dissimilarity matrix',"average", "single", "complete", "ward")
                write.csv(agglom_result_dat1,paste0('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Manuscript/',file_name,'.csv'),row.names = FALSE)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Now We take a look at the Dendrogram of the phenotypes of the residuals of the PRSs
plot_hc(results.S,"Scaled","Agglomerative","Correlation",agg_method = "ward.D2",width=2500,height=1500)
agglom_coef(file_name = 'Agglomerative_Coef_PRSice_no_MHC')
plot_heat(results.S,"Scaled","Agglomerative","Correlation",agg_method = "ward.D2",width=1500,height=1500)

##
##
##
##
# Let's take a look at Subset 3: only the ICD10

# Choosing P-val < 0.0001 and results.S:

prs_matrix <- "PCAresults_p_val_0_0001.rds"
label_subset=levels(as.factor(meta_pheno[which(meta_pheno$trait_type=='icd10'),]$category_manual))
trait_type <- 'icd10'
matrix <- results.S[[prs_matrix]]$residuals[,-dim(results.S[[prs_matrix]]$residuals)[2]]
matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
# including trait type-----------------------------------------------------------------------------------------------------
matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
#--------------------------------------------------------------------------------------------------------------------------
# label_trait_type <- NULL
cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
clusters <- hclust(as.dist(1-cols.cor),method = "ward.D2")
# Cut tree into 17 groups
k =17
sub_grp <- cutree(clusters, k = k)
# Number of members in each cluster
table(sub_grp)
sub_grp[which((sub_grp == 17)|
                              (sub_grp == 16)|
                              (sub_grp == 15)|
                              (sub_grp == 14)|
                              (sub_grp == 13)|
                              (sub_grp == 12)|
                              (sub_grp == 11)|
                              (sub_grp == 9))]
meta_pheno$category_manual[which(meta_pheno$new_pheno_annot %in% 
                                                 names(sub_grp[which((sub_grp == 17)|
                                                                                     (sub_grp == 16)|
                                                                                     (sub_grp == 15)|
                                                                                     (sub_grp == 14)|
                                                                                     (sub_grp == 13)|
                                                                                     (sub_grp == 12)|
                                                                                     (sub_grp == 11)|
                                                                                     (sub_grp == 9))]))]


# Let's remove some of the categories: lifestyle_environment, environmental_injuries,injuries_poisonings
label_subset=c("abnormal_symptoms","circulatory_system","congenital_anomalies","dermatologic",                 
               "digestive","endocrine_metabolic","genitourinary","immune_blood","infectious_diseases",          
               "mental_disorders_neurological","musculoskeletal","neoplasms","respiratory","sense_organs")
plot_hc(results.S,"Scaled","Agglomerative","Correlation",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500,trait_type = 'icd10')

# Let's take a look at the clusters:
matrix <- results.S[[prs_matrix]]$residuals[,-dim(results.S[[prs_matrix]]$residuals)[2]]
matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
# including trait type-----------------------------------------------------------------------------------------------------
matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
#--------------------------------------------------------------------------------------------------------------------------
# label_trait_type <- NULL
cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
clusters <- hclust(as.dist(1-cols.cor),method = "ward.D2")
# Cut tree into 14 groups
k =14
sub_grp <- cutree(clusters, k = k)
# Number of members in each cluster
table(sub_grp)
meta_pheno$category_manual[which(meta_pheno$new_pheno_annot %in% 
                                                 names(sub_grp[which((sub_grp == 11)|
                                                                                     (sub_grp == 12)|
                                                                                     (sub_grp == 13)|
                                                                                     (sub_grp == 14))]))]
sub_grp[which((sub_grp == 11)|
                              (sub_grp == 12)|
                              (sub_grp == 13)|
                              (sub_grp == 14))]

groups <- as.data.frame(sub_grp)
groups$pheno <- rownames(groups)
groups$pheno_category <- meta_pheno[which(meta_pheno$new_pheno_annot %in% groups$pheno),]$category_manual
table(groups$sub_grp,groups$pheno_category)
write.csv(groups,'/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Manuscript/Hierarchical Clustering/P_val_0_0001_subset.csv')



# Let's remove some of the categories: lifestyle_environment, environmental_injuries,injuries_poisonings
label_subset=c("abnormal_symptoms","circulatory_system","congenital_anomalies","dermatologic",                 
               "digestive","endocrine_metabolic","genitourinary","immune_blood","infectious_diseases",          
               "mental_disorders_neurological","musculoskeletal","neoplasms","respiratory","sense_organs")
plot_hc(results.S,"Scaled","Agglomerative","Correlation",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500,trait_type = 'icd10')

# Let's take a look at the clusters:
matrix <- results.S[[prs_matrix]]$residuals[,-dim(results.S[[prs_matrix]]$residuals)[2]]
matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
# including trait type-----------------------------------------------------------------------------------------------------
matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
#--------------------------------------------------------------------------------------------------------------------------
# label_trait_type <- NULL
cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
clusters <- hclust(as.dist(1-cols.cor),method = "ward.D2")
# Cut tree into 14 groups
k =14
sub_grp <- cutree(clusters, k = k)
# Number of members in each cluster
table(sub_grp)
meta_pheno$category_manual[which(meta_pheno$new_pheno_annot %in% 
                                                 names(sub_grp[which((sub_grp == 11)|
                                                                                     (sub_grp == 12)|
                                                                                     (sub_grp == 13)|
                                                                                     (sub_grp == 14))]))]
sub_grp[which((sub_grp == 11)|
                              (sub_grp == 12)|
                              (sub_grp == 13)|
                              (sub_grp == 14))]

groups <- as.data.frame(sub_grp)
groups$pheno <- rownames(groups)
groups$pheno_category <- meta_pheno[which(meta_pheno$new_pheno_annot %in% groups$pheno),]$category_manual
table(groups$sub_grp,groups$pheno_category)
write.csv(groups,'/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Manuscript/Hierarchical Clustering/P_val_0_0001_subset.csv')

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


                # ------------------------------------- #
                # Divisive Hierarchical Clustering      #
                # ------------------------------------- #

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function will create the best performing clustering sructure of the matrices:
divisive_coef <- function(label_subset = levels(as.factor(meta_pheno$category_manual)),
                          file_name,
                          trait_type = levels(as.factor(meta_pheno$trait_type))){
                divisive_result_dat <- matrix(ncol = 4)
                for(prs_matrix in names(results.S)){
                                matrix <- results.S[[prs_matrix]]$residuals[,-dim(results.S[[prs_matrix]]$residuals)[2]]
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
                                # including trait type-----------------------------------------------------------------------------------------------------
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
                                #--------------------------------------------------------------------------------------------------------------------------
                                for(dis_mat in c('Manhattan', 'Euclidean', 'Correlation')) {
                                                if (dis_mat == "Manhattan") {
                                                                dclust <- diana(dist(t(matrix),method="manhattan"),metric = "manhattan")
                                                                row <- dclust$dc
                                                }
                                                if (dis_mat == "Euclidean") {
                                                                dclust <- diana(dist(t(matrix),method="euclidean"),metric = "euclidean")
                                                                row <- dclust$dc
                                                }
                                                if (dis_mat == "Correlation") {
                                                                cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                                                                dclust <- diana(as.dist(1-cols.cor))
                                                                row <- dclust$dc
                                                }
                                                row <- c(prs_matrix,method,dis_mat,as.vector(row))
                                                divisive_result_dat <- rbind(divisive_result_dat,row)
                                }

                }
                divisive_result_dat1 <- divisive_result_dat[2:dim(divisive_result_dat)[1],]
                divisive_result_dat1 <- as.data.frame(divisive_result_dat1)
                names(divisive_result_dat1)<-c('P-value threshhold','scaling method','dissimilarity matrix',"DC")
                write.csv(divisive_result_dat1,paste0('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Manuscript/Hierarchical Clustering/',file_name,'.csv'),row.names = FALSE)
                
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Now We take a look at the Dendrogram of the phenotypes of the residuals of the PRSs

# file_name = 
plot_hc(results.S,"Scaled","Divisive","Manhattan",agg_method = "",width=2500,height=1500)


# Let's permute the matrix and see the results:

set.seed(1234)

for(prs_matrix in names(results1)){
                matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
                matrix1 <-as.data.frame(lapply(as.data.frame(matrix), sample))
                # Pairwise correlation between samples (columns)
                cols.cor <- cor(matrix1, use = "pairwise.complete.obs", method = "pearson")
                dclust <- diana(as.dist(1-cols.cor)) #   dist(t(matrix1),method="manhattan"),metric = "manhattan"
                print(prs_matrix)
                print(dclust$dc)
}
# meta_pheno$herit_label



