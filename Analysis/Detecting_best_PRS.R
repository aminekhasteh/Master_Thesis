######################
prs_path_worked <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PLINK_Worked/"
prs_path_2 <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PLINK/"
prs_path_1 <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PLINK(1.9)/"

# PRS(1.9) and PRS_worked are the same

# Working PRS
prs_list_worked <- list()
prs_filenames <- list.files(prs_path_worked,full.names = T)
for (prsfile in prs_filenames[1:10]) {
                prs_list_worked[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
}

# PLINK2
prs_list_2 <- list()
prs_filenames1 <- list.files(prs_path_2,full.names = T)
for (prsfile in prs_filenames1[1:10]) {
                prs_list_2[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
}

# PLINK1.9
prs_list_1 <- list()
prs_filenames2 <- list.files(prs_path_1,full.names = T)
for (prsfile in prs_filenames2[1:10]) {
                prs_list_1[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- as.data.frame(fread(prsfile,header=T,sep=","))
}


a <- prs_list_worked$p_val_0_0001.txt
b <- prs_list_2$p_val_0_0001.txt
c <- prs_list_1$p_val_0_0001.txt

dim(a)
dim(b)
dim(c)

d <- a$`icd10_G30_G30_Alzheimer's_disease_`
f <- b$`IC_Alzheimer's_disease_h0.14431_n910`
all.equal(a,c)

################################

prs_path_prs_res <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PRSice_Resid_PCA_Clust/"
prs_path_2_res <- "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/CLUMP_500_0.2/PRS_PLINK_Resid_PCA_Clust/"

# Working PRS/PLINK1.9
results.prs <- list()
prs_filenames <- list.files(prs_path_prs_res,full.names = T)
for (prsfile in prs_filenames) {
                results.prs[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
}


# PLINK2
results.S2 <- list()
prs_filenames <- list.files(prs_path_2_res,full.names = T)
for (prsfile in prs_filenames) {
                results.S2[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
}


a <- results.prs$PCAresults_p_val_0_0001.rds
b <- results.S2$PCAresults_p_val_0_0001.rds

f <- a$residuals$`IC_Alzheimer's_disease_h0.14431_n910`
e <- b$residuals$`IC_Alzheimer's_disease_h0.14431_n910`

cor.test(f,e)


############################


