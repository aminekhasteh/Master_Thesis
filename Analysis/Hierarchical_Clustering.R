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
library(Hmisc)
library(dplyr)

# Reading datasets:
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated.csv")

geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_Resid_PCA_Clust/",full.names = T)

check_sim_pheno <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_Resid_PCA_Clust/duplicated_phenotypes.rds")

# read in PRS files
results1 <- list()
for (prsfile in prs_filenames[2:11]) {
                results1[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
}

##############################################################################
#############################   |            |   #############################                                    
#############################   | Dendograms |   #############################
#############################   |            |   #############################
##############################################################################

# We take a look at the Dendogram of the phenotypes of the residuals of the PRSs.

for(prs_matrix in names(results1)){
                matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
                matrix <- matrix[,which(!names(matrix) %in% check_sim_pheno[[prs_matrix]])] # removing any duplicate phenotype
                
                names <- strsplit(colnames(matrix),'_')
                label_trait_type <- NULL
                label_manual <- NULL
                label_h <- NULL
                for(i in 1:length(names)){
                                label_trait_type <- append(label_trait_type,names[[i]][1])
                                label_h <- append(label_h,as.numeric(strsplit(names[[i]][(length(names[[i]])-1)],'h')[[1]][2]))
                }
                for(i in names(matrix)){
                                label_manual <- append(label_manual,meta_pheno[which(meta_pheno$phenocode_annotate_lst==i),]$category_manual)
                }
                
                # label <- rev(levels(as.factor(label_trait_type))) # This label didn't work well
                # label <- rev(levels(as.factor(ntile(label_h, 4)))) # This label didn't work as well!
                label <- rev(levels(as.factor(label_manual)))
                
                clusters <- hclust(dist(t(matrix)))
                dend <- as.dendrogram(clusters)
                # order it the closest we can to the order of the observations:
                dend <- rotate(dend, 1:length(clusters$labels)) 
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
                dend <- set(dend, "labels_cex", 0.5)
                # And plot:
                if((prs_matrix == 'PCAresults_p_val_1.rds')|
                   (prs_matrix == "PCAresults_p_val_1e-05.rds")|
                   (prs_matrix == "PCAresults_p_val_1e-06.rds")|
                   (prs_matrix == "PCAresults_p_val_1e-07.rds")|
                   (prs_matrix == "PCAresults_p_val_5e-08.rds")){
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Pheno_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 1200, height = 1500)
                                par(cex=1.3)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  unlist(strsplit(prs_matrix, "[_.]"))[4]),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                legend("topleft", legend = label, fill = rainbow_hcl(5))
                                dev.off()
                                # circlize_dendrogram(dend)
                                
                } else {
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Pheno_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 1200, height = 1500)
                                par(cex=1.3, mar=c(5, 8, 4, 1))
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at phenotypes level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  paste0(unlist(strsplit(prs_matrix, "[_.]"))[4], '.',
                                                         unlist(strsplit(prs_matrix, "[_.]"))[5])),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                legend("topleft", legend = label, fill = rainbow_hcl(5))
                                dev.off()
                                # circlize_dendrogram(dend)
                }
}


###We take a look at the Dendogram of the phenotypes of the residuals of PRSs at an individual level.



for(prs_matrix in names(results1)){
                matrix <- results1[[prs_matrix]]$residuals[,-dim(results1[[prs_matrix]]$residuals)[2]]
                matrix <- t(matrix)
                colnames(matrix) <- results1[[prs_matrix]]$residuals$IID
                clusters <- hclust(dist(t(matrix)))
                dend <- as.dendrogram(clusters)
                # order it the closest we can to the order of the observations:
                dend <- rotate(dend, 1:(length(clusters$height)+1))
                # Color the branches based on the clusters:
                dend <- color_branches(dend, k = 5) #, groupLabels=iris_species)
                # We hang the dendrogram a bit:
                dend <- hang.dendrogram(dend, hang_height = 0.1)
                # reduce the size of the labels:
                dend <- set(dend, "labels_cex", 0.5)
                # And plot:
                par(mar = c(3, 3, 3, 7))
                if((prs_matrix == 'PCAresults_p_val_1.rds')|
                   (prs_matrix == "PCAresults_p_val_1e-05.rds")|
                   (prs_matrix == "PCAresults_p_val_1e-06.rds")|
                   (prs_matrix == "PCAresults_p_val_1e-07.rds")|
                   (prs_matrix == "PCAresults_p_val_5e-08.rds")){
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Ind_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 1200, height = 1500)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at individuals level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  unlist(strsplit(prs_matrix, "[_.]"))[4]),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                dev.off()
                                # circlize_dendrogram(dend)
                } else {
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Dendograms/Dendogram_Ind_",
                                            gsub("PCAresults_","",gsub(".rds","",prs_matrix)),".jpg"),
                                     width = 1200, height = 1500)
                                plot(dend, cex.main=1.25,
                                     main =paste0('Dendogram of PRSs at individuals level',"\n",
                                                  dim(clusters[[1]])[1],' individuals',"\n",
                                                  'at P-value threshold of ',
                                                  paste0(unlist(strsplit(prs_matrix, "[_.]"))[4], '.',
                                                         unlist(strsplit(prs_matrix, "[_.]"))[5])),
                                     horiz =  TRUE,
                                     nodePar = list(cex = .007)
                                )
                                dev.off()
                                # circlize_dendrogram(dend)
                }
}


