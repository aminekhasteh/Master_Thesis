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

geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_Resid_PCA_Clust/",full.names = T)

check_sim_pheno <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_Resid_PCA_Clust/duplicated_phenotypes.rds")

# read in PRS files
results1 <- list()
for (prsfile in prs_filenames[2:11]) {  #[2:11]
                results1[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
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

# We take a look at the Dendogram of the phenotypes of the residuals of the PRSs.

# Agglomerative

for(prs_matrix in names(results1)){
                matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
                # methods to assess
                m <- c( "average", "single", "complete", "ward")
                names(m) <- c( "average", "single", "complete", "ward")
                
                # Pairwise correlation between samples (columns)
                cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                # function to compute coefficient
                ac <- function(x) {
                                agnes(as.dist(1-cols.cor), method = x)$ac  # dist(t(matrix),method="manhattan")
                }
                print(prs_matrix)
                print(map_dbl(m, ac))
}


# The Agglomorative clustering with Ward method performed the best
for(prs_matrix in names(results1)){
                matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
                #matrix <- matrix[,which(!names(matrix) %in% check_sim_pheno[[prs_matrix]])] # removing any duplicate phenotype
                
                names <- strsplit(colnames(matrix),'_')
                label_trait_type <- NULL
                label_manual <- NULL
                label_h <- NULL
                for(i in 1:length(names)){
                                label_trait_type <- append(label_trait_type,names[[i]][1])
                                label_h <- append(label_h,as.numeric(strsplit(names[[i]][(length(names[[i]])-1)],'h')[[1]][2]))
                }
                for(i in names(matrix)){
                                label_manual <- append(label_manual,meta_pheno[which(meta_pheno$new_pheno_annot==i),]$category_manual)
                }
                
                # label <- rev(levels(as.factor(label_trait_type))) # This label didn't work well
                # label <- rev(levels(as.factor(ntile(label_h, 4)))) # This label didn't work as well!
                label <- rev(levels(as.factor(label_manual)))
                #cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                #clusters <- hclust(dist(t(matrix),method="manhattan"),method = "ward.D2") # as.dist(1-cols.cor)
                clusters <- diana(dist(t(matrix),method="euclidean"),metric = 'euclidean') # as.dist(1-cols.cor)
                dend <- as.dendrogram(clusters)
                # order it the closest we can to the order of the observations:
                dend <- rotate(dend, 1:(length(clusters$height)+1)) #length(clusters$labels) 
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
                #dend <- set(dend, "leaves_cex", 1.2)
                # And plot:
                if((prs_matrix == 'PCAresults_p_val_1.rds')|
                   (prs_matrix == "PCAresults_p_val_1e-05.rds")|
                   (prs_matrix == "PCAresults_p_val_1e-06.rds")|
                   (prs_matrix == "PCAresults_p_val_1e-07.rds")|
                   (prs_matrix == "PCAresults_p_val_5e-08.rds")){
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Pheno_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 3500, height = 2500)
                                par(cex=1.6)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  unlist(strsplit(prs_matrix, "[_.]"))[4]),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                legend("bottomleft", legend = label, fill = rainbow_hcl(5))
                                dev.off()
                                # circlize_dendrogram(dend)
                                
                } else {
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Pheno_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 3500, height = 2500)
                                par(cex=1.6)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  paste0(unlist(strsplit(prs_matrix, "[_.]"))[4], '.',
                                                         unlist(strsplit(prs_matrix, "[_.]"))[5])),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                legend("bottomleft", legend = label, fill = rainbow_hcl(5))
                                dev.off()
                                # circlize_dendrogram(dend)
                }
}

# Let's take a look at some of the odd clusters:
  ## p_val_1e-05
  ### k=20
prs_matrix <- "PCAresults_p_val_1e-05.rds"
matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
clusters <- hclust(dist(t(matrix)),method = "ward.D2")
# Cut tree into 11 groups (based on the dendogram, we are selecting 11 groups of phenotypes)
k =20
sub_grp <- cutree(clusters, k = k)
# Number of members in each cluster
table(sub_grp)

sub_grp[which(sub_grp == 20)]
meta_pheno$category_manual[which(names(sub_grp[which(sub_grp == 20)]) %in% meta_pheno$new_pheno_annot)]
sub_grp[which(sub_grp == 19)]
meta_pheno$category_manual[which(names(sub_grp[which(sub_grp == 19)]) %in% meta_pheno$new_pheno_annot)]
sub_grp[which(sub_grp == 12)]
table(meta_pheno$category_manual[which(names(sub_grp[which(sub_grp == 1)]) %in% meta_pheno$new_pheno_annot)])

table(meta_pheno$category_manual[which(names(sub_grp) %in% meta_pheno$new_pheno_annot)])

new_matrix <- matrix
for(i in 1:11){
                print(i)
                # Let's see which method performs better for our new dataset
                m <- c( "average", "single", "complete", "ward")
                names(m) <- c( "average", "single", "complete", "ward")
                # function to compute coefficient
                ac <- function(x) {
                                agnes(t(new_matrix) , method = x)$ac
                }
                print(map_dbl(m, ac))
                
                clusters <- hclust(dist(t(new_matrix)),method = "ward.D2")
                # Cut tree into 11 groups (based on the dendogram, we are selecting 11 groups of phenotypes)
                k = 4
                sub_grp <- cutree(clusters, k = k)
                # Number of members in each cluster
                print(table(sub_grp))
                print("IC_Alzheimer's_disease_h0.14431_n910" %in% names(sub_grp))
                new_matrix <- new_matrix[,names(sub_grp[which(sub_grp == 1)])]
                
}

k = 20
sub_grp <- cutree(clusters, k = k)
# Number of members in each cluster
print(table(sub_grp))



# Divisive

 for(prs_matrix in names(results1)){
                 matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
                 # Pairwise correlation between samples (columns)
                 # cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                 dclust <- diana(dist(t(matrix),method="manhattan")) #   as.dist(1-cols.cor)
                 print(prs_matrix)
                 print(dclust$dc)
 }

for(prs_matrix in names(results1)){
                matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
                #matrix <- matrix[,which(!names(matrix) %in% check_sim_pheno[[prs_matrix]])] # removing any duplicate phenotype
                matrix <-as.data.frame(lapply(as.data.frame(matrix), sample))
                names <- strsplit(colnames(matrix),'_')
                label_trait_type <- NULL
                label_manual <- NULL
                label_h <- NULL
                for(i in 1:length(names)){
                                label_trait_type <- append(label_trait_type,names[[i]][1])
                                label_h <- append(label_h,as.numeric(strsplit(names[[i]][(length(names[[i]])-1)],'h')[[1]][2]))
                }
                for(i in names(matrix)){
                                label_manual <- append(label_manual,meta_pheno[which(meta_pheno$new_pheno_annot==i),]$category_manual)
                }
                
                # label <- rev(levels(as.factor(label_trait_type))) # This label didn't work well
                # label <- rev(levels(as.factor(ntile(label_h, 4)))) # This label didn't work as well!
                label <- rev(levels(as.factor(label_manual)))
                cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                #clusters <- hclust(dist(t(matrix),method="manhattan"),method = "ward.D2") #  dist(t(matrix),method="euclidean")
                clusters <- diana(as.dist(1-cols.cor)) # as.dist(1-cols.cor)  dist(t(matrix),method="manhattan"),metric="euclidean"
                dend <- as.dendrogram(clusters)
                # order it the closest we can to the order of the observations:
                dend <- rotate(dend, 1:(length(clusters$height)+1)) #length(clusters$labels) 
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
                #dend <- set(dend, "leaves_cex", 1.2)
                # And plot:
                if((prs_matrix == 'PCAresults_p_val_1.rds')|
                   (prs_matrix == "PCAresults_p_val_1e-05.rds")|
                   (prs_matrix == "PCAresults_p_val_1e-06.rds")|
                   (prs_matrix == "PCAresults_p_val_1e-07.rds")|
                   (prs_matrix == "PCAresults_p_val_5e-08.rds")){
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Pheno_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 3500, height = 2500)
                                par(cex=1.6)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  unlist(strsplit(prs_matrix, "[_.]"))[4]),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                legend("bottomleft", legend = label, fill = rainbow_hcl(5))
                                dev.off()
                                # circlize_dendrogram(dend)
                                
                } else {
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Pheno_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 3500, height = 2500)
                                par(cex=1.6)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  paste0(unlist(strsplit(prs_matrix, "[_.]"))[4], '.',
                                                         unlist(strsplit(prs_matrix, "[_.]"))[5])),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                legend("bottomleft", legend = label, fill = rainbow_hcl(5))
                                dev.off()
                                # circlize_dendrogram(dend)
                }
}

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


