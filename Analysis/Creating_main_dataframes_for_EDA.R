library(data.table)
library(openxlsx)
library(rms)
library(limma)
library(ggplot2)
library(factoextra)
library(diptest)
library(corrplot)
library(fuzzyjoin)
library(harrietr) # dist to long format
library(lsr)
library(dplyr)
library(tibble)

ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated_new.csv")

geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))


# -------------------------------------# Here we can alternate between the PRS with no MHC and with MHC # -------------------------------------------------#                                                                                       #|
                                                                                                                                                           #|                                                                                                     
prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice/",full.names = T) ; no_mhc = FALSE      #|
#prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_No_MHC/",full.names = T) ; no_mhc = TRUE #|
                                                                                                                                                           #|
#----------------------------------------------------------------------------------------------------------------------------------------------------------#


# read in PRS files
prs_list <- list()
for (prsfile in prs_filenames[1:10]) {
                prs_list[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
}

# read the SNP count file
snp_count <- as.data.frame(fread('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice/SNP_count/snp_count_cumulative.txt',header=T))

snp_count$category <- c("p_val_1.txt","p_val_0_1.txt","p_val_0_05.txt","p_val_0_01.txt",
                        "p_val_0_001.txt","p_val_0_0001.txt","p_val_1e-05.txt","p_val_1e-06.txt", 
                        "p_val_1e-07.txt","p_val_5e-08.txt")               
                
## LOOP for PCA analysis
dsets_torun <- names(prs_list)

allres <- list()
for (thresh in dsets_torun) {
                pval <- gsub(".txt","",thresh)
                print(paste("Starting run, p-value threshold=",pval))
                print(paste("number of subjects =",nrow(prs_list[[thresh]])))
                print(paste("number of PRS =",ncol(prs_list[[thresh]])-2))
                
                formod <- merge(prs_list[[thresh]],geno_pcs,by=c("IID","FID"))
                
                # this could be removed from loop, but in case any subjects are removed/added to either geno_pcs or PRS files, the design matrix will need to be recalculated, so left in to prevent error
                designvars <- c(names(geno_pcs)[-c(1:2)])  
                designtext <- paste("model.matrix( ~ ",paste("formod$",designvars,sep="",collapse=" + "),")",sep="")
                design <- eval(parse(text=designtext))
                
                # OLS regression using machinery build for gene expression in limma
                print(paste("Fitting model with genetic PCs"))
                fit <- lmFit(t(formod[,3:ncol(prs_list[[thresh]])]),design)
                
                print(paste("Extracting model residuals"))
                res <- t(residuals.MArrayLM(fit,y = t(formod[,3:ncol(prs_list[[thresh]])])))
                #res <- apply(res,2,scale) 
                
                res <- scale(res,scale=FALSE) # Change this later maybe
                # # SNP count flag:
                # print(paste("Number of phenotypes with more than 100 SNPs used in their PRS is",
                #             (sum(snp_count[which(snp_count$category == thresh),] > 99)-1))) # removing 1 : category column
                # pheno_passed <- names(snp_count[,which(snp_count[which(snp_count$category == thresh),
                #                                                  ] > 99)])[-sum(snp_count[which(snp_count$category == thresh),
                #                                                                           ] > 99)]
                # 
                # res <- as.matrix(as.data.frame(res)[,which(names(as.data.frame(res)) %in% pheno_passed)])
                
                # flag multimodal scores
                print(paste("Flagging scores with multimodal distributions"))
                dts <- apply(res,2,dip.test)
                
                dtsp <- lapply(dts,function(x) { x$p.value })
                multimodal <- which(dtsp < 0.05/ncol(res))
                
                print(paste("Removing",length(multimodal),"multimodal scores"))
                
                if (length(multimodal)>0) {
                                multimodal.names <- colnames(res)[multimodal]
                                res.clean <- res[,-multimodal]
                } else {
                                multimodal.names <- NA
                                res.clean <- res
                }
                
                # Change the phenotype labels
                
                colnames(res.clean) <- meta_pheno[which(meta_pheno$phenocode_annotate_lst %in% colnames(res.clean)),]$new_pheno_annot
                
                # Use Euclidean distance and correlation test to remove duplicated or almost duplicated phenotypes
                
                matrix <- as.data.frame(res.clean[,-dim(res.clean)[2]])
                matrix_dist <- as.matrix(dist(t(matrix)))
                dist_df <- melt_dist(matrix_dist)
                # #dist_df_passed <- dist_df[which(log(dist_df$dist) >= quantile(log(dist_df$dist),0.025)),]
                # check <- dist_df[which(log(dist_df$dist) < quantile(log(dist_df$dist),0.025)),]
                # to_remove = NULL
                # for (i in 1:dim(check)[1]){
                #                 corr <- cor.test(as.numeric(unlist(matrix[check[i,1]])),
                #                                  as.numeric(unlist(matrix[check[i,2]])))
                #                 if((corr$p.value <0.05/dim(matrix[2])) & (abs(corr$estimate)>0.95)){
                #                                 if (grepl( "IC", matrix[check[i,1]], fixed = TRUE)){
                #                                                 to_remove <- append(to_remove,check[i,2])
                #                                 } else if (grepl( "IC", matrix[check[i,2]], fixed = TRUE)){
                #                                                 to_remove <- append(to_remove,check[i,1])
                #                                 } else if(as.numeric(unlist(strsplit(unlist(strsplit(check[i,1],'_'))[length(unlist(strsplit(check[i,1],'_')))],'n'))[2]) > 
                #                                           as.numeric(unlist(strsplit(unlist(strsplit(check[i,2],'_'))[length(unlist(strsplit(check[i,2],'_')))],'n'))[2])){
                #                                                 to_remove <- append(to_remove,check[i,2])
                #                                 } else if (as.numeric(unlist(strsplit(unlist(strsplit(check[i,1],'_'))[length(unlist(strsplit(check[i,1],'_')))],'n'))[2]) < 
                #                                            as.numeric(unlist(strsplit(unlist(strsplit(check[i,2],'_'))[length(unlist(strsplit(check[i,2],'_')))],'n'))[2])){
                #                                                 to_remove <- append(to_remove,check[i,1])
                #                                 } else if (as.numeric(unlist(strsplit(unlist(strsplit(check[i,1],'_'))[(length(unlist(strsplit(check[i,1],'_')))-1)],'h'))[2]) < 
                #                                            as.numeric(unlist(strsplit(unlist(strsplit(check[i,2],'_'))[(length(unlist(strsplit(check[i,2],'_')))-1)],'h'))[2])){
                #                                                 to_remove <- append(to_remove,check[i,1])
                #                                 } else if(as.numeric(unlist(strsplit(unlist(strsplit(check[i,1],'_'))[(length(unlist(strsplit(check[i,1],'_')))-1)],'h'))[2]) > 
                #                                           as.numeric(unlist(strsplit(unlist(strsplit(check[i,2],'_'))[(length(unlist(strsplit(check[i,2],'_')))-1)],'h'))[2])){
                #                                                 to_remove <- append(to_remove,check[i,2])
                #                                 } else {
                #                                                 to_remove <- append(to_remove,sample(c(check[i,1], check[i,2]), 1))
                #                                 }
                #                 }
                # }
                # 
                # print(paste0(length(to_remove)," duplicated phenotypes were found"))
                
                #res.clean <- as.matrix(as.data.frame(res.clean)[,which(!colnames(res.clean)%in%to_remove)])
                
                
                # run PCA
                print(paste("Running PCA"))
                respca <- prcomp(res.clean)
                
                # Calculate subject and variable contributions to each component
                print(paste("Calculating variable and subject contributions to components"))
                res.var <- get_pca_var(respca)
                res.ind <- get_pca_ind(respca)
                
                # Changing res.clean --> dataframe
                res.clean <- as.data.frame(res.clean)
                res.clean$IID <- formod$IID
                
                # Save files in list
                print(paste("Saving data"))
                
                allres <- list(residuals=res.clean,
                               multimodal=multimodal.names,
                               pca=respca,
                               res.var=res.var,
                               res.ind=res.ind)
                if(isTRUE(no_mhc)){
                                saveRDS(allres,file=paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_No_MHC_Resid_PCA_Clust/Scaled/", "PCAresults_",pval,".rds"))
                } else{
                                saveRDS(allres,file=paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice_Resid_PCA_Clust/Scaled/", "PCAresults_",pval,".rds"))         
                }
}

