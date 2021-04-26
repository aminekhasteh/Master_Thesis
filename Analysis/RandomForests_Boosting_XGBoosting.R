# Importing libraries we need:
library(data.table)
library(randomForest)
library(randomForestExplainer)
library(metaforest)
library(ggplot2)
library(gbm)


# Reading ROS/MAP phenotype dataset
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/ROSMAP_Phenotype/ROSmaster.rds")
# Reading Filtered PNUKBB manifest
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/code/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated.csv")
# Reading PCA of the ROS/MAP Genotype dataset
geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/PCA_Genotype/geno_qc.eigenvec_new.txt",header=F)
names(geno_pcs) <- c("FID","IID",paste0("genoPC",seq(1:10)))

cum=TRUE
MHC=TRUE
APOE=FALSE

if(isTRUE(cum)){
                if(isTRUE(MHC) & !isTRUE(APOE)){
                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC_Resid/"
                                #path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Cum/No_MHC/"
                }
                if(!isTRUE(MHC) & isTRUE(APOE)){
                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_APOE_Resid/"
                                #path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Cum/No_APOE/"
                }
                if(isTRUE(MHC) & isTRUE(APOE)){
                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_No_MHC_APOE_Resid/"
                                #path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Cum/No_MHC_APOE/"
                }
                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Cum_P-val/PRS_PRSice_Resid/"
                                #path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Cum/With_MHC_APOE/"
                }
} else {
                if(isTRUE(MHC) & !isTRUE(APOE)){
                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC_Resid/"
                                #path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Non_Cum/No_MHC/"
                }
                if(!isTRUE(MHC) & isTRUE(APOE)){
                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_APOE_Resid/"
                                #path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Non_Cum/No_APOE/"
                }
                if(isTRUE(MHC) & isTRUE(APOE)){
                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_No_MHC_APOE_Resid/"
                                #path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Non_Cum/No_MHC_APOE/"
                }
                if(!isTRUE(MHC) & !isTRUE(APOE)){
                                path = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Datasets/CLUMP_500_0.2/Non_Cum_P-val/PRS_PRSice_Resid/"
                                #path_to_save = "/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/PCA/Non_Cum/With_MHC_APOE/"
                }
}

# read in residuals of PRSs files
lmatrix <- list()
prs_filenames <- list.files(path,full.names = T)
for (prsfile in prs_filenames[which(grepl(".rds", prs_filenames))]) {
                lmatrix[[tail(strsplit(prsfile,split="/")[[1]],n=1)]] <- readRDS(prsfile)
}

# Selecting variables we want to use in our model:
rosmap_phenos <- c("FID","IID","amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt",
                   "neuroticism_12","anxiety_20items","cesdsum_lv","educ","mf3123","it3123","vm3123","pput3123",
                   "stroke_ever","thyroid_ever","headinjrloc_ever","diabetes_sr_rx_ever","cancer_ever","hypertension_ever",
                    "age_bl", "msex", "heart_ever", "claudication_ever", "smoking","dcfdx_lv") # "apoe_genotype",
rosmap_response <- c("amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt",
                     "neuroticism_12","anxiety_20items","cesdsum_lv","mf3123","it3123","vm3123","pput3123","dcfdx_lv")
rosmap_predictors <- rosmap_phenos[which(!rosmap_phenos %in% rosmap_response)][3:15]

# Combine all datasets together:


main_data <- matrix(nrow = 2052)
for (prs_matrix in names(lmatrix)){
                tmp_dat <-lmatrix[[prs_matrix]]$residuals
                tmp_dat <- tmp_dat[,which(!names(tmp_dat)%in%c("IID","FID"))]
                if((prs_matrix == 'RESIDresults_p_val_1.rds')|
                   (prs_matrix == "RESIDresults_p_val_1e-05.rds")|
                   (prs_matrix == "RESIDresults_p_val_1e-06.rds")|
                   (prs_matrix == "RESIDresults_p_val_1e-07.rds")|
                   (prs_matrix == "RESIDresults_p_val_5e-08.rds")|
                   (prs_matrix == "RESIDresults_p_val_1e-08.rds")|
                   (prs_matrix == "RESIDresults_p_val_1e-09.rds")|
                   (prs_matrix == "RESIDresults_p_val_5e-05.rds")|
                   (prs_matrix == "RESIDresults_p_val_5e-06.rds")|
                   (prs_matrix == "RESIDresults_p_val_5e-07.rds")|
                   (prs_matrix == "RESIDresults_p_val_5e-09.rds")){
                                p_val <- paste0(unlist(strsplit(prs_matrix, "[_.]"))[4])
                } else {
                                p_val <- paste0(unlist(strsplit(prs_matrix, "[_.]"))[4], '.',
                                                unlist(strsplit(prs_matrix, "[_.]"))[5])
                }
                names(tmp_dat) <- paste0(names(tmp_dat),'_',p_val)
                main_data <- cbind(main_data,tmp_dat)
}
main_data <- as.data.frame(main_data[,-1])
main_data$IID <- lmatrix[[prs_matrix]]$residuals$IID
main_data <- merge(ROSmaster[,rosmap_phenos],main_data,by.x=c("IID"),by.y=c("IID"))

# Changing some of the colnames for modelling:

names(main_data) <- gsub("-", "_", names(main_data))
names(main_data) <- gsub("'", "", names(main_data))
names(main_data) <- gsub(";", "", names(main_data))
names(main_data) <- gsub('[(]', "", names(main_data))
names(main_data) <- gsub(")", "", names(main_data))
names(main_data) <- gsub('[[]', "", names(main_data))
names(main_data) <- gsub("[]]", "", names(main_data))
names(main_data) <- gsub("%", "", names(main_data))
names(main_data) <- gsub("[+]", "", names(main_data))
names(main_data) <- gsub("_", "_", names(main_data))
names(main_data) <- gsub("__", "_", names(main_data))
names(main_data) <- gsub('["]', "", names(main_data))

# dim(main_data): 2046X11921

# First, we need to scale the ROS/MAP predictors we selected:
                ## "educ","stroke_ever", "thyroid_ever", "headinjrloc_ever"
                ## "cancer_ever", "hypertension_ever", "apoe_genotype", "age_bl" "msex"               
                ## "heart_ever", "claudication_ever", "smoking", "diabetes_sr_rx_ever"

# The only continous variables we have for now: educ, age_bl

main_data$age_bl <- scale(main_data$age_bl)
main_data$educ <- scale(main_data$educ)

# Second, we need to remove the response variables we don't need to use:
                ## "amyloid_sqrt","tangles_sqrt","cogn_global_random_slope",
                ## "age_death","parksc_lv_sqrt",
                ## "neuroticism_12","anxiety_20items","cesdsum_lv","mf3123",
                ## "it3123","vm3123","pput3123","dcfdx_lv"
# Sacling response variables:
resp <- "cogn_global_random_slope"
main_data[,resp] <- scale(main_data[,resp])
# removing the rest of the response variables we don't need
main_data <- main_data[,-which(names(main_data) %in% rosmap_response[which(rosmap_response!=resp)])]

# Removing FID and IID:
main_data <- main_data[,-which(names(main_data) %in% c("FID","IID"))]

# removing any NAs:
main_data <- na.omit(main_data)
print(dim(main_data))

# Setting up the training and test data sets:
set.seed(21545476)
train <- sample(1:nrow(main_data), round(0.80*nrow(main_data)))
test <- main_data[-train,]

# Now, we model:
bag.mod1 <- randomForest(main_data[,resp]~.,data=main_data[,-which(names(main_data) %in% resp)],
                         subset=train,mtry=100,ntree=500,importance=TRUE,keep.inbag=TRUE)                
bag.mod1

# Call:
#                 randomForest(formula = main_data[, resp] ~ ., data = main_data[,      -which(names(main_data) %in% resp)], mtry = 100, ntree = 500,      importance = TRUE, subset = train) 
# Type of random forest: regression
# Number of trees: 500
# No. of variables tried at each split: 100
# 
# Mean of squared residuals: 0.9122338
# % Var explained: 0.57

# Call:
#                 randomForest(formula = main_data[, resp] ~ ., data = main_data[,      -which(names(main_data) %in% resp)], mtry = ncol(main_data),      importance = TRUE, subset = train) 
# Type of random forest: regression
# Number of trees: 500
# No. of variables tried at each split: 11943
# 
# Mean of squared residuals: 0.9329572
# % Var explained: -1.69

yhat.mod1 <- predict(bag.mod1,newdata=test)
plot(yhat.mod1,test[,resp])
abline(0,1)
mean((yhat.mod1-test[,resp])^2)
varImpPlot(bag.mod1)

##################
# Boosting
start_time <- Sys.time()
boost.mod1 <-gbm(main_data[,resp]~.,data=main_data[,-which(names(main_data) %in% resp)],
                 distribution = "gaussian",n.trees=50,interaction.depth = 2,cv.folds = 5)
end_time <- Sys.time()
end_time - start_time

summary_boosting_mod1.1 <- summary(boost.mod1)
write.csv(summary_boosting_mod1,"/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/boosting_results_mod1.csv",row.names=FALSE)
saveRDS(boost.mod1,"/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Thesis_Project/Code/Tables/boosting_results_mod1.rds")

###########
# gbm(formula = main_data[, resp] ~ ., distribution = "gaussian", 
#     data = main_data[, -which(names(main_data) %in% resp)], n.trees = 50, 
#     interaction.depth = 2, cv.folds = 5)
# A gradient boosted model with gaussian loss function.
# 50 iterations were performed.
# The best cross-validation iteration was 11.
# There were 11942 predictors of which 18 had non-zero influence.

# Predictions:
yhat.boost.mod1 <- predict(boost.mod1,newdata = test,n.trees=50)
mean((yhat.boost.mod1-test[,resp])^2)
# 0.7743661

plot(yhat.boost.mod1,test[,resp], 
     xlab="predicted response on test set",ylab="Actual response values on test set")
abline(0,1)

test_outliers <- test[which(test$cogn_global_random_slope <=-4),]
test_outliers$IID

hist((yhat.boost.mod1-test[,resp]),xlab="predicted-actual", main='Distribution of residuals on test set')


g <- ggplot(summary_boosting_mod1.1[1:18,], aes(x=reorder(var, rel.inf), y=rel.inf)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle("Variable Importance")+
                coord_flip()
g

plot(boost.mod1,i="PH_Dementias_h0.09265_n2229_5e_08")
plot(boost.mod1,i="CA_nasal_sinus_disorder_h0.12930_n1805_0.005")
plot(boost.mod1,i="cancer_ever")
plot(boost.mod1,i="CA_PRS_alfacalcidol_h0.11230_n443_0.1")
plot(boost.mod1,i="CA_Dilation_of_urethra_NEC_h0.13670_n986_0.0001")
