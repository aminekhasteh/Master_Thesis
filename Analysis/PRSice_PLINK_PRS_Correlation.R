library(corrplot)

plink_prs <- as.data.frame(fread('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PLINK/p_val_0_0001.txt',header=T,sep=","))
prsice_prs <- as.data.frame(fread('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice/p_val_0.0001_PRSice.txt',header=T,sep=","))

sum(names(prsice_prs) %in% names(plink_prs))
which(names(prsice_prs) %in% names(plink_prs))
plink_prs_new <- plink_prs[,which(names(plink_prs)%in%names(prsice_prs))]

cor_lst <- NULL
for(item in names(plink_prs_new)[3:1205]){
                cor <- cor(plink_prs_new[,item],prsice_prs[,item])
                if(cor < 0.2){
                                print(item)
                }
                cor_lst <- append(cor_lst,cor)
}
summary(cor_lst)
plot(cor_lst,xlab = 'phenotype index',ylab='correlation')
sum(cor_lst>0.95)

prsice_prs$phecode_275.53_Disorders_of_phosphorus_metabolism_


snp_count_plink <- as.data.frame(fread('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PLINK/SNP_count/snp_count_p_val_0_0001.txt',header=T,sep=","))
snp_count_prsice <- as.data.frame(fread('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PRS_PRSice/SNP_count/snp_count_cumulative.txt',header=T,sep="\t"))
snp_count_prsice$categorical_6143_Transport_type_for_commuting_to_job_workplace_None_of_the_above
summary(snp_count_plink$categorical_6143_Transport_type_for_commuting_to_job_workplace_None_of_the_)
summary(snp_count_plink$categorical_6143_Transport_type_for_commuting_to_job_workplace_None_of_the_above)
