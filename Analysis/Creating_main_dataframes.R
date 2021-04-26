library(data.table)
library(rms)
library(limma)
library(ggplot2)
library(diptest)
library(lsr)
library(dplyr)
library(harrietr) # dist to long format

primary_resid_data_GEN <- function(cum=TRUE,
                                   MHC=TRUE,
                                   APOE=TRUE,
                                   use_snp_count=TRUE,
                                   snp_count_n=5,
                                   cor_thresh =0.8) {
                # Reading ROS/MAP phenotype dataset
                ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/ROSMAP_Phenotype/ROSmaster.rds")
                # Reading Filtered PNUKBB manifest
                meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated.csv")
                # Reading phenos we want to remove additionally
                pheno_remove <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/phenos_to_remove.txt")
                names(pheno_remove) <- "pheno"
                # Reading PCA of the ROS/MAP Genotype dataset
                geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
                names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

                # ----------------------------------------------------------- # Reading PRS matrices # ----------------------------------------------------------- #                                                                                       #|
                
                if(isTRUE(cum)){
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC"
                                                snp_categories <- c("p_val_1.txt","p_val_0_1.txt","p_val_0_05.txt","p_val_0_01.txt",
                                                                    "p_val_0_005.txt","p_val_0_001.txt","p_val_0_0005.txt","p_val_0_0001.txt",
                                                                    "p_val_5e-05.txt","p_val_1e-05.txt","p_val_5e-06.txt","p_val_1e-06.txt",
                                                                    "p_val_5e-07.txt","p_val_1e-07.txt","p_val_5e-08.txt","p_val_1e-08.txt",
                                                                    "p_val_5e-09.txt")
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_APOE"
                                                snp_categories <- c("p_val_1.txt","p_val_0_1.txt","p_val_0_05.txt","p_val_0_01.txt",
                                                                    "p_val_0_005.txt","p_val_0_001.txt","p_val_0_0005.txt","p_val_0_0001.txt",
                                                                    "p_val_5e-05.txt","p_val_1e-05.txt","p_val_5e-06.txt","p_val_1e-06.txt",
                                                                    "p_val_5e-07.txt","p_val_1e-07.txt","p_val_5e-08.txt","p_val_1e-08.txt",
                                                                    "p_val_5e-09.txt")
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC_APOE"
                                                snp_categories <- c("p_val_1.txt","p_val_0_1.txt","p_val_0_05.txt","p_val_0_01.txt",
                                                                    "p_val_0_005.txt","p_val_0_001.txt","p_val_0_0005.txt","p_val_0_0001.txt",
                                                                    "p_val_5e-05.txt","p_val_1e-05.txt","p_val_5e-06.txt","p_val_1e-06.txt",
                                                                    "p_val_5e-07.txt","p_val_1e-07.txt","p_val_5e-08.txt","p_val_1e-08.txt",
                                                                    "p_val_5e-09.txt")
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice"
                                                snp_categories <- c("p_val_1.txt","p_val_0_1.txt","p_val_0_05.txt","p_val_0_01.txt",
                                                                    "p_val_0_005.txt","p_val_0_001.txt","p_val_0_0005.txt","p_val_0_0001.txt",
                                                                    "p_val_5e-05.txt","p_val_1e-05.txt","p_val_5e-06.txt","p_val_1e-06.txt",
                                                                    "p_val_5e-07.txt","p_val_1e-07.txt","p_val_5e-08.txt")
                                }
                } else {
                                snp_categories <- c("p_val_0_1.txt","p_val_0_05.txt","p_val_0_01.txt",
                                                    "p_val_0_005.txt","p_val_0_001.txt","p_val_0_0005.txt","p_val_0_0001.txt",
                                                    "p_val_0_00005.txt","p_val_0_00001.txt","p_val_0_000005.txt","p_val_0_000001.txt",
                                                    "p_val_0_0000005.txt","p_val_0_0000001.txt","p_val_0_00000005.txt","p_val_0.txt")
                                if(isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC"
                                }
                                if(!isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_APOE"
                                }
                                if(isTRUE(MHC) & isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC_APOE"
                                }
                                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice"
                                }
                }
                
                prs_filenames <- list.files(path,full.names = T)
                if(isTRUE(use_snp_count)){
                                snp_count <- as.data.frame(fread(paste0(path,'/SNP_count/snp_count.txt'),header=T))
                                snp_count$category <- snp_categories
                }
                #------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
                # read in PRS files
                prs_list <- list()
                for (prsfile in prs_filenames[which(grepl(".txt", prs_filenames))]) {
                                prs_list[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
                }

                ## LOOP for PCA analysis
                dsets_torun <- names(prs_list)
                allres <- list()
                details <- matrix(ncol = 6)
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
                                
                                res <- scale(res,scale=TRUE)

                                #Finding columns with NAs
                                na_cols <- colnames(res)[ apply(res, 2, anyNA) ]
                                print(paste0('number of columns with NAs ',length(na_cols)))
                                res <- res[,which(!colnames(res) %in% na_cols)]
                                
                                #SNP count flag:
                                if (isTRUE(use_snp_count)){
                                                pheno_failed <- names(snp_count[,which(snp_count[which(snp_count$category == thresh),
                                                ] < snp_count_n)])[-sum(snp_count[which(snp_count$category == thresh),
                                                ] < snp_count_n)]
                                                print(paste("Number of phenotypes with less than",snp_count_n,"SNPs used in each phenotype is",length(pheno_failed))) # removing 1 : category column
                                                
                                                if (length(pheno_failed)>0) {
                                                                lowsnpcount.names <- pheno_failed
                                                                res <- res[, !colnames(res) %in% pheno_failed]
                                                } else {
                                                                lowsnpcount.names <- NA
                                                                res <- res
                                                }
                                                
                                }
                                
                                # Change the phenotype labels
                                i1 <- match(colnames(res), meta_pheno$phenocode_annotate_lst)
                                i2 <- !is.na(i1) # to take care of non matches which are NA
                                colnames(res)[i2] <- meta_pheno$new_pheno_annot[i1[i2]]
                                
                                # Removing extra phenotypes we selected manually
                                phenos_to_remove <- which(colnames(res) %in% pheno_remove$pheno)
                                
                                if (length(phenos_to_remove)>0) {
                                                phenos_to_remove.phenos.names <- pheno_remove$pheno[which(pheno_remove$pheno %in% colnames(res))]
                                                res.clean <- res[,-phenos_to_remove]
                                } else {
                                                phenos_to_remove.phenos.names <- NA
                                                res.clean <- res
                                }
                                
                                print(paste0('removing ',length(phenos_to_remove),' acute diseases, except for acute mental or nueroligical phenotypes, and other administrative phenotypes, selected manually.'))


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

                                
                                # Use correlation test to remove duplicated or almost duplicated phenotypes
                                
                                matrix <- as.data.frame(res.clean[,-dim(res.clean)[2]])
                                
                                cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
                                matrix_dist <- as.matrix(as.dist(cols.cor))
                                
                                #matrix_dist <- as.matrix(dist(t(matrix)))
                                dist_df <- melt_dist(matrix_dist)
                                #dist_df_passed <- dist_df[which(log(dist_df$dist) >= quantile(log(dist_df$dist),0.025)),]
                                check <- dist_df[which(abs(dist_df$dist) > cor_thresh ),]
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

                                # Changing res.clean --> dataframe
                                res.clean <- as.data.frame(res.clean)
                                res.clean$IID <- formod$IID
                                
                                # Save files in list
                                print(paste("Saving data"))
                                
                                allres <- list(residuals=res.clean,
                                               multimodal=multimodal.names,
                                               highlycorrelated=check)
                                
                                row <- c(pval,(ncol(prs_list[[thresh]])-2),length(pheno_failed),length(multimodal),length(to_remove),dim(res.clean)[2])
                                details <- rbind(details,row)
                                
                                saveRDS(allres,file=paste0(path,"_Resid/", "RESIDresults_",pval,".rds"))

                }
                details <- details[2:dim(details)[1],]
                details <- as.data.frame(details)
                names(details)<-c('P-value threshhold','Number of PRSs',paste0("Number of phenotypes with less than ", snp_count_n," SNPs"), 
                                  "Number of phenotypes with multimodal distribution", 
                                  paste0("Number of duplicated phenotypes with |corr|<",cor_thresh), "Number of PRSs after QC")
                write.csv(details,paste0(path,"_Resid/",'details.csv'),row.names = FALSE)
}

#############################################################

primary_resid_data_GEN(cum = TRUE, MHC = TRUE, APOE = FALSE, use_snp_count = TRUE, snp_count_n = 5, cor_thresh = 0.8)
primary_resid_data_GEN(cum = TRUE, MHC = TRUE, APOE = TRUE, use_snp_count = TRUE, snp_count_n = 5, cor_thresh = 0.8)
primary_resid_data_GEN(cum = TRUE, MHC = FALSE, APOE = FALSE, use_snp_count = TRUE, snp_count_n = 5, cor_thresh = 0.8)
