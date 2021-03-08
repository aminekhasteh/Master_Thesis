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

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#############################################################################################
################################   |                   |   ##################################                                   
################################   | PCA on Phenotypes |   ##################################
################################   |                   |   ##################################
#############################################################################################

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#results.S <- results.S2
# This loop will find the p-values of each PC for each selected phenotype in ROSMAP, and generate plots for them
results.G <- list()
for (pc in names(results.S)){
                res.pcs <- as.data.frame(results.S[[pc]]$pca$x)
                res.pcs$IID <- results.S[[pc]]$residuals$IID
                pheno_pcs <- merge(res.pcs,ROSmaster,by="IID")
                pheno_pcs2 <- merge(pheno_pcs,geno_pcs,by="IID")
                ## identify relationships of PCs with phenotypes of interest
                pheno.list <- c("amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt","neuroticism_12","anxiety_20items","cesdsum_lv","educ","mf3123","it3123","vm3123","pput3123")
                
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
                
                if((pc == 'PCAresults_p_val_1.rds')|
                   (pc == "PCAresults_p_val_1e-05.rds")|
                   (pc == "PCAresults_p_val_1e-06.rds")|
                   (pc == "PCAresults_p_val_1e-07.rds")|
                   (pc == "PCAresults_p_val_5e-08.rds")){
                                g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                                                geom_bar(stat="identity",position="dodge")+
                                                geom_hline(yintercept = -log10(0.05/nrow(assocres)),col="red",lty=2)+ 
                                                geom_text(data=subset(assocres,p<0.05/nrow(assocres)),aes(label=pc))+
                                                theme_minimal()+
                                                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                                                ggtitle(paste0('PRSs at P-value threshhold of',"\n",
                                                               unlist(strsplit(pc, "[_.]"))[4]))
                                
                } else {
                                g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                                                geom_bar(stat="identity",position="dodge")+
                                                geom_hline(yintercept = -log10(0.05/nrow(assocres)),col="red",lty=2)+ 
                                                geom_text(data=subset(assocres,p<0.05/nrow(assocres)),aes(label=pc))+
                                                theme_minimal()+
                                                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                                                ggtitle(paste0('PRSs at P-value threshhold of',"\n",
                                                               paste0(unlist(strsplit(pc, "[_.]"))[4], '.',
                                                                      unlist(strsplit(pc, "[_.]"))[5])))
                }
                results.G[[pc]] <- list(resultsols=assocres,
                                        plot=g)
                
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


################################   |                    |   ##################################                                   
                               #   |     P < 0.0001     |   #
################################   |                    |   ##################################
p_val <- 'PCAresults_p_val_0_0001.rds'
results.G[[p_val]]$plot

# PC8 for cognitive decline and amyloid_sqrt:
pcnum <- 12
### get contributing PRS for given PC
results.S[[p_val]]$res.var$cos2[,pcnum][order(results.S[[p_val]]$res.var$cos2[,pcnum],decreasing = T)][1:20]

fviz_cos2(results.S$`PCAresults_p_val_1e_05.rds`$pca,choice="var",axes=pcnum, top=10)
fviz_contrib(results.S$`PCAresults_p_val_1e_05.rds`$pca,choice="var",axes=pcnum, top=30)


################################   |                    |   ##################################                                   
                               #   |     P < 1e-05      |   #
################################   |                    |   ##################################
p_val <- 'PCAresults_p_val_1e-05.rds'
results.G[[p_val]]$plot

# PC11 for cognitive decline and amyloid_sqrt:
pcnum <- 9
### get contributing PRS for given PC
results.S[[p_val]]$res.var$cos2[,pcnum][order(results.S[[p_val]]$res.var$cos2[,pcnum],decreasing = T)][1:10]

fviz_cos2(results.S$`PCAresults_p_val_1e_05.rds`$pca,choice="var",axes=pcnum, top=10)
fviz_contrib(results.S$`PCAresults_p_val_1e_05.rds`$pca,choice="var",axes=pcnum, top=30)


################################   |                    |   ##################################                                   
                               #   |     P < 0.01      |   #
################################   |                    |   ##################################
p_val <- 'PCAresults_p_val_0_01.rds'
results.G[[p_val]]$plot

# PC11 for cognitive decline and amyloid_sqrt:
pcnum <- 8
### get contributing PRS for given PC
results.S[[p_val]]$res.var$cos2[,pcnum][order(results.S[[p_val]]$res.var$cos2[,pcnum],decreasing = T)][1:10]

fviz_cos2(results.S$`PCAresults_p_val_1e_05.rds`$pca,choice="var",axes=pcnum, top=10)
fviz_contrib(results.S$`PCAresults_p_val_1e_05.rds`$pca,choice="var",axes=pcnum, top=30)


################################   |                    |   ##################################                                   
                               #   |     P < 1e-06      |   #
################################   |                    |   ##################################
p_val <- 'PCAresults_p_val_1e-05.rds'
results.G[[p_val]]$plot

# PC15 for cognitive decline and amyloid_sqrt and tangles:
pcnum <- 6
### get contributing PRS for given PC
results.S[[p_val]]$res.var$cos2[,pcnum][order(results.S[[p_val]]$res.var$cos2[,pcnum],decreasing = T)][1:10]

fviz_cos2(results.S$`PCAresults_p_val_1e-05.rds`$pca,choice="var",axes=pcnum, top=10)
fviz_contrib(results.S$`PCAresults_p_val_1e-05.rds`$pca,choice="var",axes=pcnum, top=10)


################################   |                    |   ##################################                                   
                               #   |     P < 1e-07      |   #
################################   |                    |   ##################################
p_val <- 'PCAresults_p_val_5e-08.rds'
results.G[[p_val]]$plot

# PC11 for cognitive decline and amyloid_sqrt and tangles:
pcnum <- 14
### get contributing PRS for given PC
results.S[[p_val]]$res.var$cos2[,pcnum][order(results.S[[p_val]]$res.var$cos2[,pcnum],decreasing = T)][1:10]

fviz_cos2(results.S$`PCAresults_p_val_1e_05.rds`$pca,choice="var",axes=pcnum, top=10)
fviz_contrib(results.S$`PCAresults_p_val_1e_05.rds`$pca,choice="var",axes=pcnum, top=30)



ggplot(data=pheno_pcs2,aes(x=PC8,y=cogn_global_random_slope))+
                geom_point()+
                #geom_point(data=subset(pheno_pcs2,IID %in% highlightiid),col="red",size=3)+
                geom_smooth(method="lm")+
                theme_minimal()









### get cos2 for top 5 components

corrplot(results.S$PCAresults_p_val_0_0001.rds$res.var$cos2[1:10,1:5],is.corr = F,tl.cex = 0.3)




pheno_pcs3 <- merge(pheno_pcs2,prs_list$p_val_1e_05.txt,by="IID")
summary(lm(data=pheno_pcs3, cogn_global_random_slope ~ PC9))
summary(lm(data=pheno_pcs3, cogn_global_random_slope ~ icd10_C71_C71_Malignant_neoplasm_of_brain))
