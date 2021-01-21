---
title: "Creating the dataframes for EDA"
output: github_document
---
                
## Packages and importing datasets:
                
```{r setup, include=FALSE}
library(data.table)
library(openxlsx)
library(rms)
library(limma)
library(ggplot2)
library(factoextra)
library(diptest)
library(corrplot)

ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/Sum_Stats/Pan_UK_Biobank_phenotype_manifest_h2_more_0.05_both_sex_non_ambigous.csv")

geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PLINK/",full.names = T)

# read in PRS files
prs_list <- list()
for (prsfile in prs_filenames[1:10]) {
                prs_list[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
}
```

## LOOP for PCA analysis

```{r setup, include=FALSE}
dsets_torun <- names(prs_list)

results1 <- list()
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
                res <- apply(res,2,scale)
                
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
                
                # run PCA
                print(paste("Running PCA"))
                respca <- prcomp(res.clean)
                
                # Calculate subject and variable contributions to each component
                print(paste("Calculating variable and subject contributions to components"))
                res.var <- get_pca_var(respca)
                res.ind <- get_pca_ind(respca)
                
                # Save files in list
                print(paste("Saving data"))
                res.clean <- as.data.frame(res.clean)
                res.clean$IID <- formod$IID
                allres <- list(residuals=res.clean,
                               multimodal=multimodal.names,
                               pca=respca,
                               res.var=res.var,
                               res.ind=res.ind)
                
                saveRDS(allres,file=paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PLINK_resid/", "PCAresults_",pval,".rds"))
                results1[[pval]] <- allres
}
```