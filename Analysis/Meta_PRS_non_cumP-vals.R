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
library(NbClust)
library(ggrepel) # for overlapping labels in ggplot



# Reading datasets:
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated_new.csv")
geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

# -----------------------------# Here we can alternate between the PRS with no MHC and with MHC # ----------------------------------------#                                                                                       #|
#|                                                                                                     
#prs_path <- "/Users/Shadow/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_Resid_PCA_Clust/"        #|
prs_path <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC_Resid_PCA_Clust/"  #|
#|  #|
#|
#-----------------------------------------------------------------------------------------------------------------------------------------#

# read in PRS files
results.S <- list()
prs_filenames <- list.files(prs_path,full.names = T)
for (prsfile in prs_filenames) {
                results.S[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
}

# First set of MEta-PRSs:

# P-val thresh between 0.0001 and 1 #

p_vals <- c("PCAresults_p_val_0_00001.rds","PCAresults_p_val_0_00005.rds",
            "PCAresults_p_val_0_0001.rds","PCAresults_p_val_0_0005.rds",
            "PCAresults_p_val_0_001.rds", "PCAresults_p_val_0_005.rds",
            "PCAresults_p_val_0_01.rds" ,"PCAresults_p_val_0_05.rds")

names <- NULL
for (p in p_vals){
                print(p)
                print(dim(results.S[[p]]$residuals))
                #names <- append(names,colnames(results.S[[p]]$residuals)) # Use for cum P-values
                names <- append(names,unlist(lapply(colnames(results.S[[p]]$residuals), function(x){sub(sprintf("^((.*?_){%d}.*?)_.*",length(str_split(x,"_")[[1]])-3 ), "\\1", x)})))
}

n_occur <- data.frame(table(names))

new_names <- as.character(n_occur[n_occur$Freq ==length(p_vals),]$names)
IID <- new_names[which(new_names=="IID")]
new_names <- new_names[-which(new_names=="IID")]

MetaPRS <- list()
MetaPRS_dat <- matrix(nrow = 2052)
MegaPRS_dat <- matrix(nrow = 2052)
for (name in new_names){
                tmp_dat <- matrix()
                P  <- c("_0_00001","_0_00005","_0_0001","_0_0005",
                        "_0_001","_0_005","_0_01","_0_05")
                i=1
                for (p in p_vals){
                                a <- results.S[[p]]$residuals[paste0(name,P[i])]
                                tmp_dat <- cbind(tmp_dat,a)
                                print(dim(tmp_dat))  
                                i=i+1
                }
                tmp_dat <- tmp_dat[,2:(length(p_vals)+1)]
                pheno_name <- colnames(tmp_dat)[1]
                colnames(tmp_dat) <- c(paste(pheno_name,"_0.00001"),
                                       paste(pheno_name,"_0.00005"),
                                       paste(pheno_name,"_0.0001"),
                                       paste(pheno_name,"_0.0005"),
                                       paste(pheno_name,"_0.001"),
                                       paste(pheno_name,"_0.005"),
                                       paste(pheno_name,"_0.01"),
                                       paste(pheno_name,"_0.05")) #,paste(pheno_name,"_0.05"),paste(pheno_name,"_0.1"),paste(pheno_name,"_1")
                #run PCA
                print(paste("Running PCA"))
                respca <- prcomp(tmp_dat)
                PC1_var_exp <- sum(summary(respca)$importance[2,1]) # variation explained by PC1
                PC2_var_exp <- sum(summary(respca)$importance[2,2]) # variation explained by PC1
                res.var <- get_pca_var(respca)
                
                ### get contributing PRS for PC1
                pc1_top6 <- as.data.frame(cbind(res.var$cos2[,1][order(res.var$cos2[,1],decreasing = T)][1:(length(p_vals))],res.var$contrib[,1][order(res.var$contrib[,1],decreasing = T)][1:(length(p_vals))]))
                pc1_top6$phenotypes <- rownames(pc1_top6)
                names(pc1_top6) <- c("Cos2","Contribution","phenotypes")
                
                # contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
                g1 <- ggplot(pc1_top6, aes(x=reorder(phenotypes, Contribution), y=Contribution)) +
                                geom_bar(stat="identity",position="dodge")+
                                theme_minimal()+
                                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC1'))+
                                coord_flip()
                
                ### get contributing PRS for PC2
                pc2_top6 <- as.data.frame(cbind(res.var$cos2[,2][order(res.var$cos2[,2],decreasing = T)][1:(length(p_vals))],res.var$contrib[,2][order(res.var$contrib[,2],decreasing = T)][1:(length(p_vals))]))
                pc2_top6$phenotypes <- rownames(pc2_top6)
                names(pc2_top6) <- c("Cos2","Contribution","phenotypes")
                
                # contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
                g2 <- ggplot(pc2_top6, aes(x=reorder(phenotypes, Contribution), y=Contribution)) +
                                geom_bar(stat="identity",position="dodge")+
                                theme_minimal()+
                                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC1'))+
                                coord_flip()
                
                MetaPRS[[name]] <- list(pca = respca,
                                        res.var = res.var,
                                        variation_explained_PC1 = PC1_var_exp,
                                        variation_explained_PC2 = PC2_var_exp,
                                        plot_cont_pc1 = g1,
                                        plot_cont_pc2 = g2)
                MetaPRS_dat <- cbind(MetaPRS_dat,respca$x[,1])
                MegaPRS_dat <- cbind(MegaPRS_dat,tmp_dat)
}

# Saving the data files:
# Change for the name:
## Meta-PRS
## Meta-PRS_No_MHC
## Meta-PRS_No_APOE
## Meta-PRS_No_MHC_APOE

# Saving the Meta-PRS as one dataframe (using only PC1 for each phenotype):
MetaPRS_dat <- as.data.frame(MetaPRS_dat[,-1])
colnames(MetaPRS_dat) <- new_names
MetaPRS_dat$IID <- results.S[[p]]$residuals$IID
write.csv(MetaPRS_dat,paste0("/Users/Shadow/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/MetaPRS/Top_6/","Meta-PRS_No_MHC",".csv"))

# Saving each Meta-PRS independently
saveRDS(MetaPRS,file=paste0("/Users/Shadow/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/MetaPRS/Top_6/","Meta-PRS_No_MHC",".rds"))


# Saving the Mega-PRS as one dataframe (using only PC1 for each phenotype):
MegaPRS_dat <- as.data.frame(MegaPRS_dat[,-1])
MegaPRS_dat$IID <- results.S[[p]]$residuals$IID
write.csv(MegaPRS_dat,paste0("/Users/Shadow/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/MetaPRS/Top_6/","Mega-PRS_No_MHC",".csv"))


###
### PC1
###
l=NULL
for(name in names(MetaPRS)){
                l <- append(l,MetaPRS[[name]]$variation_explained_PC1)
}
summary(l)
# Min.    1st Qu. Median  Mean    3rd Qu. Max. 
# 0.3016  0.3250  0.3313  0.3321  0.3380  0.3792 

###
### PC2
###
l1=NULL
for(name in names(MetaPRS)){
                l1 <- append(l1,MetaPRS[[name]]$variation_explained_PC2)
}
summary(l1)
# Min.    1st Qu. Median  Mean    3rd Qu. Max. 
# 0.1631  0.1779  0.1827  0.1825  0.1869  0.2024  

for(name in names(MetaPRS)[50:60]){
                print(MetaPRS$`IC_Alzheimer's_disease_h0.14431_n910`$plot_cont_pc1)
}


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################
################################   |                        |   #############################                                    
################################   | Clusters on Phenotypes |   #############################
################################   |          META          |   #############################
################################   |                        |   #############################
#############################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

matrix <- as.matrix(MetaPRS_dat[,1:(dim(MetaPRS_dat)[2]-1)])

m <- c( "average", "single", "complete", "ward")
ac <- function(x) {agnes(dist(t(matrix),method = "euclidean"),method = x)$ac}
map_dbl(m, ac)
# Euclidean: 0.08731433 0.06191563 0.30533196 0.37906067

# Pairwise correlation between samples (columns)
cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
# function to compute coefficient
ac <- function(x) {agnes(as.dist(1-cols.cor),method = x)$ac}
map_dbl(m, ac)
# Correlation: 0.10738047 0.06258807 0.50514810 0.56570911

#Phenos:
Nclusters1<-NbClust(t(matrix),diss=as.dist(1 - cor(matrix, use = "pairwise.complete.obs", method = "pearson")), distance = NULL, min.nc=2, max.nc=20, 
                    method = "complete", index = "ch")
table(Nclusters1$Best.partition)
# 1  2  3  4  5  6  7  8  9 10 11 
# 33 30 55 41 49 42 32 47 53 34 46
#Individuals:
Nclusters2<-NbClust(matrix,diss=as.dist(1 - cor(t(matrix), use = "pairwise.complete.obs", method = "pearson")), distance = NULL, min.nc=2, max.nc=20, 
                    method = "complete", index = "ch")
table(Nclusters2$Best.partition)
#1    2 
#1069  983 

Colv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(as.dist(1 - cor(matrix, use = "pairwise.complete.obs", method = "pearson")),method = "complete")),k= 11), hang_height = 0.1)
Rowv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(as.dist(1 - cor(t(matrix), use = "pairwise.complete.obs", method = "pearson")),method = "complete")),k= 2, hang_height = 0.1))
jpeg("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Heatmaps/Heatmap_Pheno_Meta-PRS_top8_noncum_complete.jpg",
     width = 1800, height = 1800)
par(mar=c(5,5,5,5),cex=1,font=3)
par(oma=c(15,3,9,3))
heatmap(as.matrix(matrix), Rowv = Rowv, Colv = Colv,
        scale = "none",xlab = NULL, ylab = NULL,
        col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
legend(x="topleft", legend=c("Dissimilar", "No Similarity", "Similar"), 
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
mtext(paste0('Heatmap of PRSs at phenotypes level',"\n",
             dim(matrix)[1],' individuals and ',dim(matrix)[2],' phenotypes',"\n",
             'of Meta-PRSs (top 8 P-values)'), side = 3, line = 1, cex = 2, outer=TRUE)
dev.off()
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################
################################   |                        |   #############################                                    
################################   |          PCA           |   #############################
################################   |                        |   #############################
#############################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

respca <- prcomp(MetaPRS_dat[,1:(dim(MetaPRS_dat)[2]-1)])
res.pcs <- as.data.frame(respca$x)
res.pcs$IID <- MetaPRS_dat$IID
pheno_pcs <- merge(res.pcs,ROSmaster,by="IID")
pheno_pcs2 <- merge(pheno_pcs,geno_pcs,by="IID")
## identify relationships of PCs with phenotypes of interest
pheno.list <- c("amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt",
                "neuroticism_12","anxiety_20items","cesdsum_lv","educ","mf3123","it3123","vm3123","pput3123",
                "stroke_ever","thyroid_ever","headinjrloc_ever","diabetes_sr_rx_ever","cancer_ever","hypertension_ever")

PClist <- paste0("PC",seq(1:50))
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

assocres$colour <- ifelse(assocres$b < 0, "Negative effect","Positive effect")
g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                geom_bar(stat="identity",position="dodge",aes(fill = colour))+
                geom_hline(aes(yintercept = -log10(0.05/nrow(assocres)),color="P-value < 0.0000526"),lty=2)+ 
                geom_hline(aes(yintercept = -log10(0.05),color="P-value < 0.05"),lty=2)+ 
                scale_linetype_manual(name = "limit", values = c(2, 2), 
                                      guide = guide_legend(override.aes = list(color = c("green", "orange"))))+
                geom_text_repel(data=subset(assocres,p<0.05),aes(label=pc))+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                ggtitle(('Meta-PRS'))


res.var <- get_pca_var(respca)


for (pcnum in 1:20){
                var <-sum(summary(respca)$importance[2,pcnum])*100 
                ### get contributing PRS for given PC
                pc_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
                pc_top20$phenotypes <- rownames(pc_top20)
                names(pc_top20) <- c("Cos2","Contribution","phenotypes")
                # represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
                
                # contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
                ggplot(pc_top20, aes(x=reorder(phenotypes, Contribution), y=Contribution)) +
                                geom_bar(stat="identity",position="dodge")+
                                theme_minimal()+
                                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum, 
                                               ' (',var,'% variation explained',')'))+
                                coord_flip()
                ggsave(paste("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/PCA/",
                             "META_top8_noncumPvals_PC",pcnum,"_top20_cont.jpg"), width = 30, height = 27, units = "cm")
}



# # cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
# ggplot(pc_top20, aes(x=reorder(phenotypes, Cos2), y=Cos2)) +
#                 geom_bar(stat="identity",position="dodge")+
#                 theme_minimal()+
#                 theme(axis.text.x=element_text(angle = -90, hjust = 0))+
#                 ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum))+
#                 coord_flip()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################
################################   |                        |   #############################                                    
################################   |       Modelling        |   #############################
################################   |                        |   #############################
#############################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

pca_dat <- as.data.frame(cbind(respca$x))
pca_dat$IID <- MetaPRS_dat$IID
merged_dat <- merge(ROSmaster,pca_dat,by="IID")
merged_dat <- merge(merged_dat,MetaPRS_dat)


############
mod1 <- lm(merged_dat$globcog_random_slope~merged_dat$PC1+merged_dat$PC2+merged_dat$PC3+merged_dat$PC4+merged_dat$PC5+merged_dat$PC6+
                           merged_dat$PC7+merged_dat$PC8+merged_dat$PC9+merged_dat$PC10+merged_dat$PC11+merged_dat$PC12+merged_dat$PC13+
                           merged_dat$PC14+merged_dat$PC15+merged_dat$PC16+merged_dat$PC17+merged_dat$PC18+merged_dat$PC19+merged_dat$PC20)
summary(mod1)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################
################################   |                        |   #############################                                    
################################   | Clusters on Phenotypes |   #############################
################################   |          MEGA          |   #############################
################################   |                        |   #############################
#############################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

matrix <- as.matrix(MegaPRS_dat[,1:(dim(MegaPRS_dat)[2]-1)])

m <- c( "average", "single", "complete", "ward")
ac <- function(x) {agnes(dist(t(matrix),method = "euclidean"),method = x)$ac}
map_dbl(m, ac)
# Euclidean: 0.2535853 0.2243504 0.2953972 0.6859888

# Pairwise correlation between samples (columns)
cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
# function to compute coefficient
ac <- function(x) {agnes(as.dist(1-cols.cor),method = x)$ac}
map_dbl(m, ac)
# Correlation: 0.4371788 0.3923307 0.4984101 0.8115135

#Phenos:
Nclusters1<-NbClust(t(matrix),diss=as.dist(1 - cor(matrix, use = "pairwise.complete.obs", method = "pearson")), distance = NULL, min.nc=2, max.nc=20, 
                    method = "complete", index = "ch")
table(Nclusters1$Best.partition)
#1    2 
#1792 1904
#Individuals:
Nclusters2<-NbClust(matrix,diss=as.dist(1 - cor(t(matrix), use = "pairwise.complete.obs", method = "pearson")), distance = NULL, min.nc=2, max.nc=20, 
                    method = "complete", index = "ch")
table(Nclusters2$Best.partition)
#1    2 
#968 1084

Colv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(as.dist(1 - cor(matrix, use = "pairwise.complete.obs", method = "pearson")),method = "complete")),k= 2), hang_height = 0.1)
Rowv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(as.dist(1 - cor(t(matrix), use = "pairwise.complete.obs", method = "pearson")),method = "complete")),k= 2, hang_height = 0.1))
jpeg("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Heatmaps/Heatmap_Pheno_Mega-PRS_noncumP_top8.jpg",
     width = 1800, height = 1800)
par(mar=c(5,5,5,5),cex=1,font=3)
par(oma=c(15,3,9,3))
heatmap(as.matrix(matrix), Rowv = Rowv, Colv = Colv,
        scale = "none",xlab = NULL, ylab = NULL,
        col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
legend(x="topleft", legend=c("Dissimilar", "No Similarity", "Similar"), 
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
mtext(paste0('Heatmap of PRSs at phenotypes level',"\n",
             dim(matrix)[1],' individuals and ',dim(matrix)[2],' phenotypes',"\n",
             'of Mega-PRSs (top 8 P-values)'), side = 3, line = 1, cex = 2, outer=TRUE)
dev.off()
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################
################################   |                        |   #############################                                    
################################   |          PCA           |   #############################
################################   |                        |   #############################
#############################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

respca <- prcomp(MegaPRS_dat[,1:(dim(MegaPRS_dat)[2]-1)])
res.pcs <- as.data.frame(respca$x)
res.pcs$IID <- MegaPRS_dat$IID
pheno_pcs <- merge(res.pcs,ROSmaster,by="IID")
pheno_pcs2 <- merge(pheno_pcs,geno_pcs,by="IID")
## identify relationships of PCs with phenotypes of interest
pheno.list <- c("amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt",
                "neuroticism_12","anxiety_20items","cesdsum_lv","educ","mf3123","it3123","vm3123","pput3123",
                "stroke_ever","thyroid_ever","headinjrloc_ever","diabetes_sr_rx_ever","cancer_ever","hypertension_ever")

PClist <- paste0("PC",seq(1:50))
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

assocres$colour <- ifelse(assocres$b < 0, "Negative effect","Positive effect")
g <- ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                geom_bar(stat="identity",position="dodge",aes(fill = colour))+
                geom_hline(aes(yintercept = -log10(0.05/nrow(assocres)),color="P-value < 0.0000526"),lty=2)+ 
                geom_hline(aes(yintercept = -log10(0.05),color="P-value < 0.05"),lty=2)+ 
                scale_linetype_manual(name = "limit", values = c(2, 2), 
                                      guide = guide_legend(override.aes = list(color = c("green", "orange"))))+
                geom_text_repel(data=subset(assocres,p<0.05),aes(label=pc))+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                ggtitle(('Meta-PRS'))


res.var <- get_pca_var(respca)

for (pcnum in 1:20){
                var <-sum(summary(respca)$importance[2,pcnum])*100 
                ### get contributing PRS for given PC
                pc_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
                pc_top20$phenotypes <- rownames(pc_top20)
                names(pc_top20) <- c("Cos2","Contribution","phenotypes")
                # represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
                
                # contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
                ggplot(pc_top20, aes(x=reorder(phenotypes, Contribution), y=Contribution)) +
                                geom_bar(stat="identity",position="dodge")+
                                theme_minimal()+
                                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum, 
                                               ' (',var,'% variation explained',')'))+
                                coord_flip()
                ggsave(paste("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/PCA/",
                             "MEgA_top8_nonCumPvals_PC",pcnum,"_top20_cont.jpg"), width = 30, height = 27, units = "cm")
}

# # cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
# ggplot(pc_top20, aes(x=reorder(phenotypes, Cos2), y=Cos2)) +
#                 geom_bar(stat="identity",position="dodge")+
#                 theme_minimal()+
#                 theme(axis.text.x=element_text(angle = -90, hjust = 0))+
#                 ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum))+
#                 coord_flip()

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################
################################   |                        |   #############################                                    
################################   |       Modelling        |   #############################
################################   |                        |   #############################
#############################################################################################
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------

pca_dat <- as.data.frame(cbind(respca$x))
pca_dat$IID <- MegaPRS_dat$IID
merged_dat <- merge(ROSmaster,pca_dat,by="IID")
merged_dat <- merge(merged_dat,MegaPRS_dat)


############