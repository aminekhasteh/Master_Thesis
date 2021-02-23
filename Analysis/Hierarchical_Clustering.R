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

# Reading datasets:
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated_new.csv")

# -----------------------------# Here we can alternate between the PRS with no MHC and with MHC # -----------------------------------#                                                                                       #|
                                                                                                                                     #|                                                                                                     
#prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_Resid_PCA_Clust/                  #|
prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_No_MHC_Resid_PCA_Clust/"           #|
                                                                                                                                     #|
#------------------------------------------------------------------------------------------------------------------------------------#

# read in PRS files
results.S <- list()
results.C <- list()
for(method in c("Scaled","Centered")){
                prs_filenames <- list.files(paste0(prs_path,'/',method),full.names = T)
                if (method == 'Scaled'){
                                for (prsfile in prs_filenames) {
                                                results.C[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
                                }
                }
                if (method == 'Centered'){
                                for (prsfile in prs_filenames) {
                                                results.S[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
                                }
                }
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
# This function generates all of the plots for us:
plot_hc <- function(lmatrix,
                    scaling_method = c("Center","Scaled"),
                    hc_method = c("Agglomerative","Divisive"),
                    dis_mat_method=c("Manhattan","Euclidean","Correlation"),
                    label_subset = levels(as.factor(meta_pheno$category_manual)),
                    agg_method = "ward.D2",
                    width=1800,
                    height=900){ # lmatrix : list of our PRS matrices
                for(prs_matrix in names(lmatrix)){
                                matrix <- lmatrix[[prs_matrix]]$residuals[,-dim(lmatrix[[prs_matrix]]$residuals)[2]]
                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
                                names <- strsplit(colnames(matrix),'_')
                                # label_trait_type <- NULL
                                label_manual <- NULL
                                # label_h <- NULL
                                # for(i in 1:length(names)){
                                #                 label_trait_type <- append(label_trait_type,names[[i]][1])
                                #                 label_h <- append(label_h,as.numeric(strsplit(names[[i]][(length(names[[i]])-1)],'h')[[1]][2]))
                                # }
                                for(i in names(matrix)){
                                                label_manual <- append(label_manual,meta_pheno[which(meta_pheno$new_pheno_annot==i),]$category_manual)
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

                # ------------------------------------- #
                # Agglomerative Hierarchical Clustering #
                # ------------------------------------- #

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function will create the best performing clustering sructure of the matrices:
agglom_coef <- function(label_subset = levels(as.factor(meta_pheno$category_manual)),
                        file_name){
                agglom_result_dat <- matrix(ncol = 7)
                for(prs_matrix in names(results.S)){
                                for(method in c("Scaled","Centered")){
                                                print(paste0(method,'-',prs_matrix))
                                                if(method == "Scaled"){
                                                                matrix <- results.S[[prs_matrix]]$residuals[,-dim(results.S[[prs_matrix]]$residuals)[2]]
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
                                                                                row <- c(prs_matrix,method,dis_mat,as.vector(row))
                                                                                agglom_result_dat <- rbind(agglom_result_dat,row)
                                                                }
                                                                
                                                }
                                                if(method == "Centered"){
                                                                matrix <- results.C[[prs_matrix]]$residuals[,-dim(results.C[[prs_matrix]]$residuals)[2]]
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
                                                                                row <- c(prs_matrix,method,dis_mat,as.vector(row))
                                                                                agglom_result_dat <- rbind(agglom_result_dat,row)
                                                                }
                                                                
                                                }
                                                
                                                
                                }
                }
                agglom_result_dat1 <- agglom_result_dat[2:dim(agglom_result_dat)[1],]
                agglom_result_dat1 <- as.data.frame(agglom_result_dat1)
                names(agglom_result_dat1)<-c('P-value threshhold','scaling method','dissimilarity matrix',"average", "single", "complete", "ward")
                write.csv(agglom_result_dat1,paste0('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Manuscript/Hierarchical Clustering/',file_name,'.csv'),row.names = FALSE)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Now We take a look at the Dendrogram of the phenotypes of the residuals of the PRSs
plot_hc(results.S,"Scaled","Agglomerative","Manhattan",agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",agg_method = "ward.D2",width=2500,height=1500)

plot_hc(results.S,"Scaled","Agglomerative","Manhattan",agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",agg_method = "complete",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",agg_method = "complete",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",agg_method = "complete",width=2500,height=1500)

plot_hc(results.S,"Scaled","Agglomerative","Manhattan",agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",agg_method = "single",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",agg_method = "single",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",agg_method = "single",width=2500,height=1500)
# Now that we have our plots, let's take a look at their coefficients:
agglom_coef(file_name = 'clustering_performance_Agglo_No_MHC')

# Let's exlcude some of the phenotypes based on their categories

                # Subset_1 #

label_subset = c("circulatory_system","neoplasms","digestive","neurological","respiratory","sense_organs",
                 "genitourinary","endocrine_metabolic","mental_disorders","infectious_diseases","musculoskeletal")
plot_hc(results.S,"Scaled","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)

plot_hc(results.S,"Scaled","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)

plot_hc(results.S,"Scaled","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
# Now that we have our plots, let's take a look at their coefficients:
agglom_coef(label_subset=label_subset,file_name = 'clustering_performance_Agglo_No_MHC_subset_1')

                # Subset_2 #

label_subset = c("circulatory_system","neoplasms","digestive","neurological","sense_organs","mental_disorders")
plot_hc(results.S,"Scaled","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",label_subset=label_subset,agg_method = "ward.D2",width=2500,height=1500)

plot_hc(results.S,"Scaled","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",label_subset=label_subset,agg_method = "complete",width=2500,height=1500)

plot_hc(results.S,"Scaled","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Manhattan",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Euclidean",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.S,"Scaled","Agglomerative","Correlation",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
plot_hc(results.C,"Center","Agglomerative","Correlation",label_subset=label_subset,agg_method = "single",width=2500,height=1500)
# Now that we have our plots, let's take a look at their coefficients:
agglom_coef(label_subset=label_subset,file_name = 'clustering_performance_Agglo_No_MHC_subset_2')
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


                # ------------------------------------- #
                # Divisive Hierarchical Clustering      #
                # ------------------------------------- #

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This function will create the best performing clustering sructure of the matrices:
divisive_coef <- function(label_subset = levels(as.factor(meta_pheno$category_manual)),
                          file_name){
                divisive_result_dat <- matrix(ncol = 4)
                for(prs_matrix in names(results.S)){
                                for(method in c("Scaled","Centered")){
                                                print(paste0(method,'-',prs_matrix))
                                                if(method == "Scaled"){
                                                                matrix <- results.S[[prs_matrix]]$residuals[,-dim(results.S[[prs_matrix]]$residuals)[2]]
                                                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
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
                                                if(method == "Centered"){
                                                                matrix <- results.C[[prs_matrix]]$residuals[,-dim(results.C[[prs_matrix]]$residuals)[2]]
                                                                matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
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
plot_hc(results.C,"Center","Divisive","Manhattan",agg_method = "",width=2500,height=1500)
plot_hc(results.S,"Scaled","Divisive","Euclidean",agg_method = "",width=2500,height=1500)
plot_hc(results.C,"Center","Divisive","Euclidean",agg_method = "",width=2500,height=1500)
plot_hc(results.S,"Scaled","Divisive","Correlation",agg_method = "",width=2500,height=1500)
plot_hc(results.C,"Center","Divisive","Correlation",agg_method = "",width=2500,height=1500)
# Now that we have our plots, let's take a look at their coefficients:
divisive_coef(file_name = 'clustering_performance_Divisive_No_MHC')

# Let's exlcude some of the phenotypes based on their categories
                
                # Subset_1 #

label_subset = c("circulatory_system","neoplasms","digestive","neurological","respiratory","sense_organs","genitourinary","endocrine_metabolic","mental_disorders","infectious_diseases","musculoskeletal")
plot_hc(results.S,"Scaled","Divisive","Manhattan",agg_method = "",label_subset=label_subset,width=2500,height=1500)
plot_hc(results.C,"Center","Divisive","Manhattan",agg_method = "",label_subset=label_subset,width=2500,height=1500)
plot_hc(results.S,"Scaled","Divisive","Euclidean",agg_method = "",label_subset=label_subset,width=2500,height=1500)
plot_hc(results.C,"Center","Divisive","Euclidean",agg_method = "",label_subset=label_subset,width=2500,height=1500)
plot_hc(results.S,"Scaled","Divisive","Correlation",agg_method = "",label_subset=label_subset,width=2500,height=1500)
plot_hc(results.C,"Center","Divisive","Correlation",agg_method = "",label_subset=label_subset,width=2500,height=1500)
# Now that we have our plots, let's take a look at their coefficients:
divisive_coef(label_subset=label_subset,file_name = 'clustering_performance_Divisive_No_MHC_subset_1')
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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

##############################################################################################
################################   |                         |   #############################                                    
################################   | Clusters on Individuals |   #############################
################################   |                         |   #############################
##############################################################################################

###We take a look at the Dendogram of the phenotypes of the residuals of PRSs at an individual level.
#for(prs_matrix in names(results1)){
#                 matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
#                 matrix <- t(matrix)
#                 colnames(matrix) <- results1[[prs_matrix]]$residuals$IID
#                 clusters <- hclust(dist(t(matrix)),method = "ward.D")
#                 dend <- as.dendrogram(clusters)
#                 # order it the closest we can to the order of the observations:
#                 dend <- rotate(dend, 1:(length(clusters$height)+1))
#                 # Color the branches based on the clusters:
#                 dend <- color_branches(dend, k = 5) #, groupLabels=iris_species)
#                 # We hang the dendrogram a bit:
#                 dend <- hang.dendrogram(dend, hang_height = 0.1)
#                 # reduce the size of the labels:
#                 dend <- set(dend, "labels_cex", 0.5)
#                 # And plot:
#                 par(mar = c(3, 3, 3, 7))
#                 if((prs_matrix == 'PCAresults_p_val_1.rds')|
#                    (prs_matrix == "PCAresults_p_val_1e-05.rds")|
#                    (prs_matrix == "PCAresults_p_val_1e-06.rds")|
#                    (prs_matrix == "PCAresults_p_val_1e-07.rds")|
#                    (prs_matrix == "PCAresults_p_val_5e-08.rds")){
#                                 jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Ind_",
#                                             gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
#                                      width = 1200, height = 1500)
#                                 plot(dend, cex.main=1.25,
#                                      main =paste0('Dendogram of PRSs at individuals level',"\n",
#                                                   dim(clusters[[1]])[1],' individuals',"\n",
#                                                   'at P-value threshold of ',
#                                                   unlist(strsplit(prs_matrix, "[_.]"))[4]),
#                                      horiz =  TRUE,
#                                      nodePar = list(cex = .007)
#                                 )
#                                 dev.off()
#                                 # circlize_dendrogram(dend)
#                 } else {
#                                 jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Ind_",
#                                             gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
#                                      width = 1200, height = 1500)
#                                 plot(dend, cex.main=1.25,
#                                      main =paste0('Dendogram of PRSs at individuals level',"\n",
#                                                   dim(clusters[[1]])[1],' individuals',"\n",
#                                                   'at P-value threshold of ',
#                                                   paste0(unlist(strsplit(prs_matrix, "[_.]"))[4], '.',
#                                                          unlist(strsplit(prs_matrix, "[_.]"))[5])),
#                                      horiz =  TRUE,
#                                      nodePar = list(cex = .007)
#                                 )
#                                 dev.off()
#                                 # circlize_dendrogram(dend)
#                 }
# }


