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

ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated_new.csv")

geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec_no_mhc_new.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

# ------------------------------------------# Here we can alternate between the PRS with no MHC and with MHC # ----------------------------------------------------------- #                                                                                       #|
#|                                                                                                     
#prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice/",full.names = T) ; no_mhc = FALSE        #|
prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_No_MHC/",full.names = T) ; no_mhc = TRUE   #|
#|
  #|
#|
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #


# read in PRS files
prs_list <- list()
for (prsfile in prs_filenames[1:10]) {
                prs_list[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
}

#read the SNP count file
snp_count <- as.data.frame(fread('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_No_MHC/SNP_count/snp_count_cumulative.txt',header=T))

snp_count$category <- c("p_val_1.txt","p_val_0_1.txt","p_val_0_05.txt","p_val_0_01.txt",
                        "p_val_0_001.txt","p_val_0_0001.txt","p_val_1e-05.txt","p_val_1e-06.txt",
                        "p_val_1e-07.txt","p_val_5e-08.txt")
# 
## LOOP for PCA analysis
dsets_torun <- names(prs_list)

allres <- list()

# for (thresh in dsets_torun) {
#                 print(thresh)
#                 print('FID'%in%colnames(prs_list[[thresh]]))
#                 print('IID'%in%colnames(prs_list[[thresh]]))
# }
# 
# # The FID name for "p_val_0_001.txt" (PRSs with MHC with PLINK) is ataFID. Let's change that:
# thresh <- "p_val_0_001.txt"
# colnames(prs_list[[thresh]])[1] <- 'FID'

for (thresh in dsets_torun) {
                pval <- gsub(".txt","",thresh)
                print(paste("Starting run, p-value threshold=",pval))
                print(paste("number of subjects =",nrow(prs_list[[thresh]])))
                print(paste("number of PRS =",ncol(prs_list[[thresh]])-2))
                names(prs_list[[thresh]]) <- gsub(",","",names(prs_list[[thresh]]),fixed=T)
                
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
                
                res <- scale(res,scale=TRUE) # Change this later maybe
                
                # Finding columns with NAs
                # na_cols <- colnames(res)[ apply(res, 2, anyNA) ]
                # print(paste0('number of columns with NAs ',length(na_cols)))
                # res <- res[,which(!colnames(res) %in% na_cols)]
                
                # SNP count flag:
                pheno_failed <- names(snp_count[,which(snp_count[which(snp_count$category == thresh),
                ] < 5)])[-sum(snp_count[which(snp_count$category == thresh),
                ] < 5)]
                print(paste("Number of phenotypes with less than 5 SNPs used in each phenotype is",length(pheno_failed))) # removing 1 : category column
                
                if (length(pheno_failed)>0) {
                                lowsnpcount.names <- pheno_failed
                                res.clean <- res[, !colnames(res) %in% pheno_failed]
                } else {
                                lowsnpcount.names <- NA
                                res.clean <- res
                }
                
                
                # flag multimodal scores - Using dip test
                print(paste("Flagging scores with multimodal distributions"))
                dts <- apply(res.clean,2,dip.test)
                dtsp <- lapply(dts,function(x) { x$p.value })
                multimodal <- which(dtsp < 0.05/ncol(res.clean))
                
                print(paste("Removing",length(multimodal),"multimodal scores"))
                
                if (length(multimodal)>0) {
                                multimodal.names <- colnames(res)[multimodal]
                                res.clean <- res.clean[,-multimodal]
                } else {
                                multimodal.names <- NA
                                res.clean <- res.clean
                }
                
                
                # Change the phenotype labels
                
                #####
                ##
                i1 <- match(colnames(res.clean), meta_pheno$phenocode_annotate_lst)
                i2 <- !is.na(i1) # to take care of non matches which are NA
                colnames(res.clean)[i2] <- meta_pheno$new_pheno_annot[i1[i2]]
                ##
                ####
                
                # Removing phenotypes with acute in them (except for mental health and neurological)
                l <- meta_pheno$new_pheno_annot[which(meta_pheno$category_manual != "mental_disorders_neurological")]
                to_remove_acute <- l[which(grepl( 'Acute', l, fixed = TRUE))]
                acute <- which(colnames(res.clean) %in% to_remove_acute)
                
                if (length(acute)>0) {
                                acute.phenos.names <- to_remove_acute
                                res.clean <- res.clean[,-acute]
                } else {
                                acute.phenos.names <- NA
                                res.clean <- res.clean
                }
                
                print(paste0('removing ',length(to_remove_acute),' acute diseases, except for acute mental or nueroligical phenotypes.'))
                
                # Use correlation test to remove duplicated or almost duplicated phenotypes
                
                matrix <- as.data.frame(res.clean[,-dim(res.clean)[2]])
                
                cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                matrix_dist <- as.matrix(as.dist(cols.cor))
                
                #matrix_dist <- as.matrix(dist(t(matrix)))
                dist_df <- melt_dist(matrix_dist)
                #dist_df_passed <- dist_df[which(log(dist_df$dist) >= quantile(log(dist_df$dist),0.025)),]
                check <- dist_df[which(abs(dist_df$dist) > 0.8 ),]
                to_remove = NULL
                if(dim(check)[1] >0){
                                for (i in 1:dim(check)[1]){
                                                corr <- cor.test(as.numeric(unlist(matrix[check[i,1]])),
                                                                 as.numeric(unlist(matrix[check[i,2]])))
                                                if((corr$p.value <0.05/dim(matrix[2]))){
                                                                if ((grepl( "IC", check[i,1], fixed = TRUE)) & 
                                                                    (grepl( "PH", check[i,2], fixed = TRUE))){
                                                                                to_remove <- append(to_remove,check[i,2])
                                                                } else if ((grepl( "IC", check[i,2], fixed = TRUE)) & 
                                                                           (grepl( "PH", check[i,1], fixed = TRUE))){
                                                                                to_remove <- append(to_remove,check[i,1])
                                                                } else {
                                                                                to_remove <- append(to_remove,sample(c(check[i,1], check[i,2]), 1))
                                                                }
                                                }
                                }
                }
                
                print(paste0(length(to_remove)," duplicated phenotypes were found"))
                
                if (length(to_remove)>0) {
                                res.clean <- res.clean[ , -which(colnames(res.clean) %in% to_remove)]
                                
                } else {
                                res.clean <- res.clean
                }
                
                # run PCA
                # print(paste("Running PCA"))
                # respca <- prcomp(res.clean)
                # 
                # # Calculate subject and variable contributions to each component
                # print(paste("Calculating variable and subject contributions to components"))
                # res.var <- get_pca_var(respca)
                # res.ind <- get_pca_ind(respca)
                
                # Changing res.clean --> dataframe
                res.clean <- as.data.frame(res.clean)
                res.clean$IID <- formod$IID
                
                # Save files in list
                print(paste("Saving data"))
                
                allres <- list(residuals=res.clean,
                               multimodal=multimodal.names,
                               highlycorrelated=check,
                               removed.acute.phenos=acute.phenos.names) 
                # ,
                # pca=respca,
                # res.var=res.var,
                # res.ind=res.ind
                
                if(isTRUE(no_mhc)){
                                
                                saveRDS(allres,file=paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_No_MHC_Resid_PCA_Clust/", "PCAresults_",pval,".rds"))
                } else{
                                saveRDS(allres,file=paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_Resid_PCA_Clust/", "PCAresults_",pval,".rds"))  
                }
}

