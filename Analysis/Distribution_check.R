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

geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

# ------------------------------------------# Here we can alternate between the PRS with no MHC and with MHC # ----------------------------------------------------------- #                                                                                       #|
#|                                                                                                     
#prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice/",full.names = T) ; no_mhc = FALSE        #|
prs_filenames <- list.files("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC/",full.names = T) ; no_mhc = TRUE   #|
#|
#|
#|
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# read in PRS files
prs_list <- list()
for (prsfile in prs_filenames[1:18]) {
                prs_list[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
}

##############################
# Plotting histograms of some of the PRSs randomly
dsets_torun <- names(prs_list)

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
                
                res.clean <- scale(res,scale=TRUE) # Change this later maybe
                
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
                set.seed(56372522)
                selected_pheno <- sample(colnames(res.clean), 15)
                
                for(pheno in selected_pheno){
                                jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/PRS_Distribution/",
                                            pval,' _',pheno,".jpg"),
                                     width = 1000, height = 1500)
                                par(mfrow=c(3,1))
                                hist(res.clean[,pheno],main=paste0('scaled residuals PRS','\n',
                                                                   pheno, ' ',pval),
                                     xlab = '',ylab='',col='lightblue')
                                hist(prs_list[[thresh]][,pheno],main=paste0('Original PRS','\n',
                                                                            pheno, ' ',pval),
                                     xlab = '',ylab='',col='brown')
                                hist(scale(prs_list[[thresh]][,pheno]),main=paste0('scaled PRS','\n',
                                                                                   pheno, ' ',pval),
                                     xlab = '',ylab='',col='green')
 
                                dev.off()
                }
                
                
                
}
