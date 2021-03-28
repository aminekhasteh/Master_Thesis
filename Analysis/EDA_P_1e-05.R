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
library(RColorBrewer)
library(NbClust)

# Reading datasets:
ROSmaster <- readRDS("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/ROSMAP_Phenotype/ROSmaster.rds")

# Will change this later
meta_pheno <- read.csv("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno_annotated_new.csv")
geno_pcs <- read.table("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Datasets/PCA_Genotype/geno_qc.eigenvec_no_mhc_new.txt",header=F)
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

################################   |                    |   ##################################                                   
                               #   |     P < 1e-05      |   #
################################   |                    |   ##################################

prs_matrix <- 'PCAresults_p_val_1e-05.rds'
matrix <- results.S[[prs_matrix]]$residuals[,-dim(results.S[[prs_matrix]]$residuals)[2]]
Nclusters<-NbClust(t(matrix), distance = "euclidean", min.nc=2, max.nc=10, 
             method = "complete", index = "ch")
table(Nclusters$Best.partition)
#-----------------------------------------------------------------------------------------------------

# subset pheno categories-----------------------------------------------------------------------------------------------------
label_subset = levels(as.factor(meta_pheno$category_manual))  #levels(as.factor(meta_pheno[which(meta_pheno$trait_type=='icd10'),]$category_manual))
matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which(meta_pheno$category_manual %in% label_subset),]$new_pheno_annot)]
# including trait type-----------------------------------------------------------------------------------------------------
trait_type <- levels(as.factor(meta_pheno$trait_type)) #'icd10'
matrix <- matrix[,which(colnames(matrix) %in% meta_pheno[which((meta_pheno$trait_type %in% trait_type)),]$new_pheno_annot)]
#--------------------------------------------------------------------------------------------------------------------------
# label_trait_type <- NULL
cols.cor <- cor(matrix, use = "pairwise.complete.obs", method = "pearson")
clusters <- hclust(as.dist(1-cols.cor),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 968   9 

sub_grp[which((sub_grp == 2))]

# CA_sarcoidosis_h0.22162_n971                                 PH_Chronic_hepatitis_h0.13914_n394 
# IC_Intestinal_malabsorption_h0.15994_n2750                   PR_carbimazole_h0.13762_n1184 
# CA_Medical_information_h0.13412_n509    ------> Doctor diagnosed sarcoidosis                     
# PH_Sicca_syndrome_h0.22589_n786 
# PH_Diffuse_diseases_of_connective_tissue_h0.06227_n1382      PH_Systemic_lupus_erythematosus_h0.12527_n540 
# PH_Lupus_(localized_and_systemic)_h0.23891_n660 


# table(meta_pheno$category_manual[which(meta_pheno$new_pheno_annot %in% 
#                                                        names(sub_grp[which((sub_grp == 2))]))])


###  |
###  | A)
###  |

# Let's take a look at sub group 1:

sub_gp_names_1 <- sub_grp[which((sub_grp == 1))]
matrix_selected_1 <- matrix[,names(sub_gp_names_1)]
cols.cor1 <- cor(matrix_selected_1, use = "pairwise.complete.obs", method = "pearson")
clusters1 <- hclust(as.dist(1-cols.cor1),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters1, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 955  13

sub_grp[which((sub_grp == 2))]

# CA_malignant_melanoma_h0.05001_n3594                                   # CA_Re-excision_of_skin_margins_NEC_h0.08608_n1152 
# CA_Re-excision_of_skin_margins_of_head_or_neck_h0.19527_n446           # IC_Malignant_melanoma_of_skin_h0.06200_n2788 
# IC_Carcinoma_in_situ_of_skin_h0.06237_n993                             # IC_Skin_changes_due_to_chronic_exposure_to_nonionising_radiation_h0.05910_n3472 
# CA_removal_of_rodent_ulcer_basal_cell_carcinoma_(bcc)_h0.17277_n1090   # CA_skin_cancer_h0.07432_n1436 
# IC_Other_malignant_neoplasms_of_skin_h0.08286_n15426                   # IC_Melanoma_in_situ_h0.20869_n1035 
# CA_Excision_of_lesion_of_external_ear_h0.08864_n2223                   # CA_Excision_of_lesion_of_external_nose_h0.06814_n3604 
# CA_basal_cell_carcinoma_h0.11939_n4278 

######  |
######  | B)
######  |

# Let's take a look at sub group 1 again:

sub_gp_names_2 <- sub_grp[which((sub_grp == 1))]
matrix_selected_2 <- matrix_selected_1[,names(sub_gp_names_2)]
cols.cor2 <- cor(matrix_selected_2, use = "pairwise.complete.obs", method = "pearson")
clusters2 <- hclust(as.dist(1-cols.cor2),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters2, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 945  10 

sub_grp[which((sub_grp == 2))]

# CA_PRS_sotalol_h0.24089_n722           # CA_PRS_flecainide_h0.05766_n738 
# PR_flecainide_h0.14622_n982            # PR_amiodarone_h0.07738_n1394 
# PR_anti-arrhythmic_h0.07920_n2231      # CA_Direct_current_cardioversion_h0.16027_n4048 
# CA_atrial_fibrillation_h0.14715_n3526  # PR_digoxin_h0.22051_n1816 
# PR_anti-arrhythmic_non-selective_beta-blocker_beta_blocker_h0.23745_n1084 
# CA_PRS_digoxin_h0.15375_n1192 


#########  |
#########  | C)
#########  |

# Let's take a look at sub group 1 again:

sub_gp_names_3 <- sub_grp[which((sub_grp == 1))]
matrix_selected_3 <- matrix_selected_2[,names(sub_gp_names_3)]
cols.cor3 <- cor(matrix_selected_3, use = "pairwise.complete.obs", method = "pearson")
clusters3 <- hclust(as.dist(1-cols.cor3),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters3, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 937   8 

sub_grp[which((sub_grp == 2))]

# PH_Other_disorders_of_gallbladder_h0.07971_n1937         PH_Calculus_of_bile_duct_h0.06227_n3253 
# CA_Administration_h0.07447_n2223 ------> Hepatobiliary & pancreatic surgery                        
# CA_cholecystitis_h0.23650_n468 
# IC_Other_diseases_of_biliary_tract_h0.06784_n2927        PH_Cholelithiasis_with_acute_cholecystitis_h0.12561_n2145 
# PH_Other_disorders_of_biliary_tract_h0.06595_n2972       IC_Other_diseases_of_gallbladder_h0.06736_n2461 

####**************************************************************####
#### Removing: PH_Other_disorders_of_biliary_tract_h0.06595_n2972 ####
####           CA_cholecystitis_h0.23650_n468                     ####
####**************************************************************####

############  |
############  | D)
############  |

# Let's take a look at sub group 1 again:

sub_gp_names_4 <- sub_grp[which((sub_grp == 1))]
matrix_selected_4 <- matrix_selected_3[,names(sub_gp_names_4)]
cols.cor4 <- cor(matrix_selected_4, use = "pairwise.complete.obs", method = "pearson")
clusters4 <- hclust(as.dist(1-cols.cor4),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters4, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 925  12

sub_grp[which((sub_grp == 2))]

# PH_Other_disorders_of_arteries_and_arterioles_h0.16520_n2182 
# PH_Atherosclerosis_of_the_extremities_h0.22051_n1156 
# PH_Atherosclerosis_h0.17687_n1948 
# IC_Occlusion_and_stenosis_of_precerebral_arteries_not_resulting_in_cerebral_infarction_h0.18076_n1448 
# PH_Arterial_embolism_and_thrombosis_h0.26356_n1271 
# CA_Percutaneous_transluminal_balloon_angioplasty_and_insertion_of_1-2_drug-eluting_stents_into_coronary_artery_h0.09059_n6513 
# IC_Other_disorders_of_arteries_and_arterioles_h0.18741_n2018 
# IC_Subsequent_myocardial_infarction_h0.15252_n859 
# CA_PRS_glyceryl_trinitrate_h0.08286_n1469 
# PH_Occlusion_and_stenosis_of_precerebral_arteries_h0.14965_n1911 
# PH_Congenital_anomalies_of_great_vessels_h0.06979_n2784 
# IC_Arterial_embolism_and_thrombosis_h0.25611_n1259 

####********************************************************************************************####
#### Removing: PH_Arterial_embolism_and_thrombosis_h0.26356_n1271                               ####
####           PH_Other_disorders_of_arteries_and_arterioles_h0.16520_n2182                     ####
####           PH_Occlusion_and_stenosis_of_precerebral_arteries_h0.14965_n1911                 ####
####********************************************************************************************####


###############  |
###############  | E)
###############  |
# Let's take a look at sub group 1 again:

sub_gp_names_5 <- sub_grp[which((sub_grp == 1))]
matrix_selected_5 <- matrix_selected_4[,names(sub_gp_names_5)]
cols.cor5 <- cor(matrix_selected_5, use = "pairwise.complete.obs", method = "pearson")
clusters5 <- hclust(as.dist(1-cols.cor5),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters5, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 920   5 

sub_grp[which((sub_grp == 2))]

# IC_Alzheimer's_disease_h0.14431_n910        PH_Dementias_h0.09265_n2229 
# CA_Administration_h0.07773_n777 ----> Old age psychiatry            
# CA_Family_history_h0.10308_n2171  -----> Alzheimer's disease/dementia                       
# PH_Delirium_due_to_conditions_classified_elsewhere_h0.06327_n1533 


##################  |
##################  | F)
##################  |
# Let's take a look at sub group 1 again:

sub_gp_names_6 <- sub_grp[which((sub_grp == 1))]
matrix_selected_6 <- matrix_selected_5[,names(sub_gp_names_6)]
cols.cor6 <- cor(matrix_selected_6, use = "pairwise.complete.obs", method = "pearson")
clusters6 <- hclust(as.dist(1-cols.cor6),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters6, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 901  19  

sub_grp[which((sub_grp == 2))]

# PH_Myeloproliferative_disease_h0.11202_n1439 
# IC_Other_neoplasms_of_uncertain_or_unknown_behaviour_of_lymphoid_haematopoietic_and_related_tissue_h0.12844_n924 
# PH_Cancer_of_brain_and_nervous_system_h0.10767_n806 
# CA_PRS_warfarin_h0.09064_n4714 
# PR_selective_estrogen_receptor_modulator_SERM_h0.06724_n3239 
# PH_Multiple_sclerosis_h0.14763_n1695 
# PH_Polycythemia_secondary_h0.07781_n439 
# CA_multiple_sclerosis_h0.24039_n1689 
# PH_Manlignant_and_unknown_neoplasms_of_brain_and_nervous_system_h0.10888_n989 
# CA_Administration_h0.05791_n1539 ------> Hepatology
# CA_Venography_h0.13726_n957 
# IC_Malignant_neoplasm_of_brain_h0.14387_n744 
# CA_PRS_migraleve_tablet_h0.11963_n635 
# PR_tamoxifen_h0.09469_n2749 
# PH_Coagulation_defects_h0.25830_n1445 
# CA_hereditary_genetic_haematological_disorder_h0.14141_n412 
# CA_mastectomy_h0.06877_n4841 
# IC_Follow-up_care_involving_plastic_surgery_h0.07468_n3195 
# IC_Other_venous_embolism_and_thrombosis_h0.13398_n578 

####************************************************####
#### Removing: CA_multiple_sclerosis_h0.24039_n1689 ####
####************************************************####


#####################  |
#####################  | G)
#####################  |
# Let's take a look at sub group 1 again:

sub_gp_names_7 <- sub_grp[which((sub_grp == 1))]
matrix_selected_7 <- matrix_selected_6[,names(sub_gp_names_7)]
cols.cor7 <- cor(matrix_selected_7, use = "pairwise.complete.obs", method = "pearson")
clusters7 <- hclust(as.dist(1-cols.cor7),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters7, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 892   9   

sub_grp[which((sub_grp == 2))]

# PH_Cancer_of_bronchus;_lung_h0.14038_n3755                              PR_tiotropium_h0.28428_n3572 
# IC_Malignant_neoplasm_of_bronchus_and_lung_h0.15624_n3423               PH_Aortic_aneurysm_h0.14852_n2228 
# IC_Emphysema_h0.27690_n2651                                             PH_Abdominal_aortic_aneurysm_h0.21482_n1418 
# PH_Other_aneurysm_h0.11643_n2862                                        PH_Cancer_within_the_respiratory_system_h0.14201_n4310 
# IC_Aortic_aneurysm_and_dissection_h0.14630_n2224 

####************************************************####
#### Removing: PH_Aortic_aneurysm_h0.14852_n2228    ####
####************************************************####


########################  |
########################  | H)
########################  |
# Let's take a look at sub group 1 again:

sub_gp_names_8 <- sub_grp[which((sub_grp == 1))]
matrix_selected_8 <- matrix_selected_7[,names(sub_gp_names_8)]
cols.cor8 <- cor(matrix_selected_8, use = "pairwise.complete.obs", method = "pearson")
clusters8 <- hclust(as.dist(1-cols.cor8),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters8, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 879  13   

sub_grp[which((sub_grp == 2))]

# PH_Simple_and_unspecified_goiter_h0.12447_n822 
# PH_Other_disorders_of_tympanic_membrane_h0.13502_n1719 
# IC_Postprocedural_endocrine_and_metabolic_disorders_not_elsewhere_classified_h0.12287_n1975 
# PH_Dermatophytosis__Dermatomycosis_h0.22433_n430 
# CA_Lobectomy_of_thyroid_gland_NEC_h0.12791_n1021 
# CA_thyroid_problem_(not_cancer)_h0.13228_n1268 
# IC_Perforation_of_tympanic_membrane_h0.14721_n1270 
# IC_Dermatophytosis_h0.11636_n356 
# CA_Total_thyroidectomy_h0.13596_n715 
# PH_Secondary_hypothyroidism_h0.13637_n1696 
# PH_Nontoxic_nodular_goiter_h0.14372_n1673 
# PH_Dermatophytosis_h0.12431_n421 
# PH_Perforation_of_tympanic_membrane_h0.14998_n1339 

####**************************************************************####
#### Removing: PH_Perforation_of_tympanic_membrane_h0.14998_n1339 ####
####           PH_Dermatophytosis__Dermatomycosis_h0.22433_n430   ####
####**************************************************************####


###########################  |
###########################  | I)
###########################  |
# Let's take a look at sub group 1 again:

sub_gp_names_9 <- sub_grp[which((sub_grp == 1))]
matrix_selected_9 <- matrix_selected_8[,names(sub_gp_names_9)]
cols.cor9 <- cor(matrix_selected_9, use = "pairwise.complete.obs", method = "pearson")
clusters9 <- hclust(as.dist(1-cols.cor9),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters9, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 870   9  

sub_grp[which((sub_grp == 2))]

# PH_Cancer_of_esophagus_h0.14266_n1068     IC_Disorders_of_autonomic_nervous_system_h0.09325_n363 
# IC_Other_strabismus_h0.13268_n1205        PH_Disorders_of_the_autonomic_nervous_system_h0.05396_n372 
# CA_COVID_19_2_h0.06504_n881               CA_COVID_19_3_h0.05460_n935 
# PH_Strabismus_(not_specified_as_paralytic)_h0.12019_n1210         
# IC_Malignant_neoplasm_of_oesophagus_h0.16215_n1057 
# CA_COVID_19_1_h0.19358_n780 

####**********************************************************************####
#### Removing: PH_Cancer_of_esophagus_h0.14266_n1068                      ####
####           CA_COVID_19_3_h0.05460_n935                                ####
####           CA_COVID_19_1_h0.19358_n780                                ####
####           PH_Disorders_of_the_autonomic_nervous_system_h0.05396_n372 ####
####           PH_Strabismus_(not_specified_as_paralytic)_h0.12019_n1210  ####
####**********************************************************************####


##############################  |
##############################  | J)
##############################  |
# Let's take a look at sub group 1 again:

sub_gp_names_10 <- sub_grp[which((sub_grp == 1))]
matrix_selected_10 <- matrix_selected_9[,names(sub_gp_names_10)]
cols.cor10 <- cor(matrix_selected_10, use = "pairwise.complete.obs", method = "pearson")
clusters10 <- hclust(as.dist(1-cols.cor10),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters10, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 854  16   

sub_grp[which((sub_grp == 2))]

# CA_Polypectomy_of_internal_nose_h0.12774_n2566 
# PR_carbonic_anhydrase_inhibitor_h0.05243_n1268 
# CA_Trabeculectomy_h0.20019_n1088 
# PH_Psoriatic_arthropathy_h0.12362_n999 
# CA_PRS_lumigan_0.3mg_ml_eye_drops_h0.16462_n320 
# PR_calcipotriol_h0.10344_n1615 
# PR_prostaglandin_analog_beta_antagonist_combination_prostaglandin_analog_beta_blocker_h0.19299_n490 
# CA_nasal_polyp_surgery_nasal_polypectomy_h0.23092_n1397 
# CA_PRS_latanoprost_h0.17362_n944 
# IC_Psoriatic_and_enteropathic_arthropathies_h0.11485_n642 
# PR_brimonidine_h0.17904_n550 
# CA_nasal_polyps_h0.13798_n2037 
# CA_PRS_xalatan_0.005%_eye_drops_h0.16593_n1354 
# CA_Employment_history_h0.06346_n6935   ------> Breathing problems responsible for leaving job
# CA_PRS_dovobet_ointment_h0.14678_n549 
# CA_Employment_history_h0.06959_n5629   ------> Breathing problems improved/stopped away from workplace or on holiday

####**********************************************************************####
#### Removing: PH_Psoriatic_arthropathy_h0.12362_n999                     ####
####                                                                      ####
####**********************************************************************####

#################################  |
#################################  | K)
#################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_11 <- sub_grp[which((sub_grp == 1))]
matrix_selected_11 <- matrix_selected_10[,names(sub_gp_names_11)]
cols.cor11 <- cor(matrix_selected_11, use = "pairwise.complete.obs", method = "pearson")
clusters11 <- hclust(as.dist(1-cols.cor11),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters11, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 843  11  

sub_grp[which((sub_grp == 2))]

# IC_Other_aplastic_anaemias_h0.14005_n682 
# CA_Ligation_of_long_saphenous_vein_h0.13943_n5384 
# IC_Congenital_deformities_of_feet_h0.12035_n285 
# PH_Congenital_deformities_of_feet_h0.13789_n293 
# PH_Aplastic_anemia_h0.06452_n696 
# CA_varicose_veins_h0.22910_n1711 
# CA_Combined_operations_on_primary_long_saphenous_vein_h0.09645_n1487 
# PH_Paraproteinemia_h0.13285_n359 
# PH_Disorders_of_plasma_protein_metabolism_h0.22800_n462 
# PH_Disorders_of_protein_plasma_amino-acid_transport_and_metabolism_h0.20387_n526 
# CA_Ligation_of_short_saphenous_vein_h0.13915_n628 


####**********************************************************************####
#### Removing: PH_Congenital_deformities_of_feet_h0.13789_n293            ####
####                                                                      ####
####**********************************************************************####


####################################  |
####################################  | L)
####################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_12 <- sub_grp[which((sub_grp == 1))]
matrix_selected_12 <- matrix_selected_11[,names(sub_gp_names_12)]
cols.cor12 <- cor(matrix_selected_12, use = "pairwise.complete.obs", method = "pearson")
clusters12 <- hclust(as.dist(1-cols.cor12),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters12, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 830  13 

sub_grp[which((sub_grp == 2))]

# PH_Colon_cancer_h0.09863_n4345 
# PH_Torsion_dystonia_h0.09055_n383 
# PH_Malignant_neoplasm_of_rectum_rectosigmoid_junction_and_anus_h0.09327_n3010 
# CA_Endoscopic_snare_resection_of_lesion_of_lower_bowel_using_fibreoptic_sigmoidoscope_h0.10516_n2843 
# IC_Alcoholic_liver_disease_h0.10576_n1211 
# IC_Other_inflammatory_liver_diseases_h0.07584_n1053 
# IC_Malignant_neoplasm_of_colon_h0.11563_n4232 
# PH_Liver_abscess_and_sequelae_of_chronic_liver_disease_h0.06127_n1610 
# IC_Dystonia_h0.09012_n376 
# PH_Portal_hypertension_h0.09393_n845 
# PH_Esophageal_bleeding_(varices_hemorrhage)_h0.13031_n2492 
# CA_colon_cancer_sigmoid_cancer_h0.06907_n1520 
# CA_Peranal_excision_of_lesion_of_rectum_h0.06545_n633 


####*********************************************************####
#### Removing: PH_Torsion_dystonia_h0.09055_n383             ####
####                                                         ####
####*********************************************************####

#######################################  |
#######################################  | M)
#######################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_13 <- sub_grp[which((sub_grp == 1))]
matrix_selected_13 <- matrix_selected_12[,names(sub_gp_names_13)]
cols.cor13 <- cor(matrix_selected_13, use = "pairwise.complete.obs", method = "pearson")
clusters13 <- hclust(as.dist(1-cols.cor13),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters13, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 821   9  

sub_grp[which((sub_grp == 2))]

# PH_Non-Hodgkins_lymphoma_h0.07787_n2552                             IC_Zoster_[herpes_zoster]_h0.12765_n729 
# IC_Diseases_of_salivary_glands_h0.13834_n1330                       PH_Herpes_zoster_h0.10593_n738 
# PH_Lump_or_mass_in_breast_h0.05668_n1917                            IC_Other_and_unspecified_types_of_non-Hodgkin's_lymphoma_h0.10630_n1661 
# PH_Cancer_of_other_lymphoid_histiocytic_tissue_h0.05339_n3267       PH_Diseases_of_the_salivary_glands_h0.12386_n1408 
# PH_Abnormal_findings_on_mammogram_or_breast_exam_h0.06330_n1985 


####***************************************************************####
#### Removing: PH_Herpes_zoster_h0.10593_n738                      ####
####           PH_Diseases_of_the_salivary_glands_h0.12386_n1408   ####
####                                                               ####
####***************************************************************####

#######################################  |
#######################################  | M)
#######################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_14 <- sub_grp[which((sub_grp == 1))]
matrix_selected_14 <- matrix_selected_13[,names(sub_gp_names_14)]
cols.cor14 <- cor(matrix_selected_14, use = "pairwise.complete.obs", method = "pearson")
clusters14 <- hclust(as.dist(1-cols.cor14),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters14, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 807  14  

sub_grp[which((sub_grp == 2))]

# IC_Cardiomyopathy_h0.13293_n1914                    CA_Endoscopic_resection_of_lesion_of_bladder_h0.08989_n2656 
# CA_crohns_disease_h0.23849_n1431                    PH_Subarachnoid_hemorrhage_h0.09891_n1158 
# IC_Subarachnoid_haemorrhage_h0.12467_n1130          PH_Suppurative_and_unspecified_otitis_media_h0.09107_n1073 
# IC_Malignant_neoplasm_of_bladder_h0.06816_n2883     IC_Crohn's_disease_[regional_enteritis]_h0.12674_n2273 
# PH_Cardiomyopathy_h0.11981_n1948                    CA_Endoscopic_cauterisation_of_lesion_of_bladder_h0.15770_n1406 
# CA_bladder_cancer_h0.09682_n1118                    PH_Benign_neoplasm_of_brain_cranial_nerves_meninges_h0.13365_n1118 
# IC_Suppurative_and_unspecified_otitis_media_h0.12996_n1004 
# PH_Benign_neoplasm_of_brain_and_other_parts_of_nervous_system_h0.13048_n1144 

####**********************************************************************####
#### Removing: PH_Cardiomyopathy_h0.11981_n1948                           ####
####           PH_Subarachnoid_hemorrhage_h0.09891_n1158                  ####
####           PH_Suppurative_and_unspecified_otitis_media_h0.09107_n1073 ####
####           CA_crohns_disease_h0.23849_n1431                           ####
####                                                                      ####
####**********************************************************************####

##########################################  |
##########################################  | N)
##########################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_15 <- sub_grp[which((sub_grp == 1))]
matrix_selected_15 <- matrix_selected_14[,names(sub_gp_names_15)]
cols.cor15 <- cor(matrix_selected_15, use = "pairwise.complete.obs", method = "pearson")
clusters15 <- hclust(as.dist(1-cols.cor15),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters15, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 801   6 

sub_grp[which((sub_grp == 2))]

# IC_Peritonitis_h0.14609_n1155                      PH_Unspecified_erythematous_condition_h0.13976_n694 
# IC_Other_erythematous_conditions_h0.14158_n773     PH_Peritonitis_and_retroperitoneal_infections_h0.15826_n1265 
# PH_Mineral_deficiency_NEC_h0.14123_n437            IC_Deficiency_of_other_nutrient_elements_h0.21767_n423 

####**********************************************************************####
#### Removing: PH_Mineral_deficiency_NEC_h0.14123_n437                    ####
####                                                                      ####
####**********************************************************************####

#############################################  |
#############################################  | O)
#############################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_16 <- sub_grp[which((sub_grp == 1))]
matrix_selected_16 <- matrix_selected_15[,names(sub_gp_names_16)]
cols.cor16 <- cor(matrix_selected_16, use = "pairwise.complete.obs", method = "pearson")
clusters16 <- hclust(as.dist(1-cols.cor16),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters16, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 796   5 

sub_grp[which((sub_grp == 2))]

# PH_Suicidal_ideation_or_attempt_h0.10331_n2260                     PH_Postinflammatory_pulmonary_fibrosis_h0.16874_n1438       
# IC_Other_interstitial_pulmonary_diseases_h0.20741_n1711            PH_Seborrheic_keratosis_h0.06244_n4169 
# PH_Suicide_or_self-inflicted_injury_h0.10974_n2149

################################################  |
################################################  | P)
################################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_17 <- sub_grp[which((sub_grp == 1))]
matrix_selected_17 <- matrix_selected_16[,names(sub_gp_names_17)]
cols.cor17 <- cor(matrix_selected_17, use = "pairwise.complete.obs", method = "pearson")
clusters17 <- hclust(as.dist(1-cols.cor17),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters17, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 786  10  

sub_grp[which((sub_grp == 2))]

# PH_Polyneuropathy_in_diabetes_h0.13099_n629                          IC_Other_disorders_of_ear_not_elsewhere_classified_h0.06987_n945 
# PH_Pernicious_anemia_h0.21874_n1045                                  PH_Pneumonitis_due_to_inhalation_of_food_or_vomitus_h0.06481_n1081 
# PH_Diseases_and_other_conditions_of_the_tongue_h0.10725_n1712        IC_Polyneuropathy_in_diseases_classified_elsewhere_h0.13163_n679 
# PH_Tinnitus_h0.06025_n812                                            IC_Vitamin_B12_deficiency_anaemia_h0.22289_n1269 
# IC_Diseases_of_tongue_h0.12034_n1699                                 IC_Pneumonitis_due_to_solids_and_liquids_h0.06263_n1117 

####**********************************************************************************####
#### Removing: PH_Polyneuropathy_in_diabetes_h0.13099_n629                            ####
####           PH_Pneumonitis_due_to_inhalation_of_food_or_vomitus_h0.06481_n1081     ####
####           PH_Diseases_and_other_conditions_of_the_tongue_h0.10725_n1712          ####
####                                                                                  ####
####**********************************************************************************####


###################################################  |
###################################################  | Q)
###################################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_18 <- sub_grp[which((sub_grp == 1))]
matrix_selected_18 <- matrix_selected_17[,names(sub_gp_names_18)]
cols.cor18 <- cor(matrix_selected_18, use = "pairwise.complete.obs", method = "pearson")
clusters18 <- hclust(as.dist(1-cols.cor18),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters18, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 778   8   

sub_grp[which((sub_grp == 2))]

# CA_Digital_fasciectomy_h0.22818_n1334                         PH_Acquired_spondylolisthesis_h0.07726_n2339                  
# IC_Fibroblastic_disorders_h0.19900_n3773                      PH_Fracture_of_neck_of_femur_h0.12761_n2127 
# PH_Other_acquired_musculoskeletal_deformity_h0.05406_n2867    PH_Fracture_of_unspecified_part_of_femur_h0.10289_n3267 
# PH_Chondrocalcinosis_h0.13002_n485                            PH_Crystal_arthropathies_h0.14731_n587 


######################################################  |
######################################################  | R)
######################################################  |
# Let's take a look at sub group 1 again:

sub_gp_names_19 <- sub_grp[which((sub_grp == 1))]
matrix_selected_19 <- matrix_selected_18[,names(sub_gp_names_19)]
cols.cor19 <- cor(matrix_selected_19, use = "pairwise.complete.obs", method = "pearson")
clusters19 <- hclust(as.dist(1-cols.cor19),method = "ward.D2")
# Cut tree into 2 groups
k = 2
sub_grp <- cutree(clusters19, k = k)
# Number of members in each cluster
table(sub_grp)

# sub_grp
# 1     2 
# 415 363    

sub_grp[which((sub_grp == 2))]

cor.test(matrix$CA_Administration_h0.23874_n1506, matrix$CA_type_1_diabetes_h0.14695_n420)

#-----------------------------------------------------------------------------------------------------
to_remove_phenos <- c("PH_Polyneuropathy_in_diabetes_h0.13099_n629",
                      "PH_Pneumonitis_due_to_inhalation_of_food_or_vomitus_h0.06481_n1081",
                      "PH_Diseases_and_other_conditions_of_the_tongue_h0.10725_n1712","PH_Mineral_deficiency_NEC_h0.14123_n437",
                      "PH_Cardiomyopathy_h0.11981_n1948",
                      "PH_Subarachnoid_hemorrhage_h0.09891_n1158",
                      "PH_Suppurative_and_unspecified_otitis_media_h0.09107_n1073",
                      "CA_crohns_disease_h0.23849_n1431",
                      "PH_Herpes_zoster_h0.10593_n738",
                      "PH_Diseases_of_the_salivary_glands_h0.12386_n1408",
                      "PH_Torsion_dystonia_h0.09055_n383",
                      "PH_Congenital_deformities_of_feet_h0.13789_n293",
                      "PH_Psoriatic_arthropathy_h0.12362_n999",
                      "PH_Cancer_of_esophagus_h0.14266_n1068",
                      "CA_COVID_19_3_h0.05460_n935",
                      "CA_COVID_19_1_h0.19358_n780",
                      "PH_Disorders_of_the_autonomic_nervous_system_h0.05396_n372",
                      "PH_Strabismus_(not_specified_as_paralytic)_h0.12019_n1210",
                      "PH_Perforation_of_tympanic_membrane_h0.14998_n1339",
                      "PH_Dermatophytosis__Dermatomycosis_h0.22433_n430",
                      "PH_Aortic_aneurysm_h0.14852_n2228",
                      "CA_multiple_sclerosis_h0.24039_n1689",
                      "PH_Arterial_embolism_and_thrombosis_h0.26356_n1271",
                      "PH_Other_disorders_of_arteries_and_arterioles_h0.16520_n2182",
                      "PH_Occlusion_and_stenosis_of_precerebral_arteries_h0.14965_n1911",
                      "PH_Other_disorders_of_biliary_tract_h0.06595_n2972",
                      "CA_cholecystitis_h0.23650_n468",
                      "CA_malignant_melanoma_h0.05001_n3594")

# Removing the selected phenotypes:
selected_phenos <- which(names(matrix)%in%to_remove_phenos)
matrix_updated <- matrix[,-selected_phenos]

# Let's take a look at the heatmap of the updated PRS matrix
label_manual <- NULL
for(i in names(matrix)){
                label_manual <- append(label_manual,meta_pheno$category_manual[which(meta_pheno$new_pheno_annot==i)])
}
label <- rev(levels(as.factor(label_manual)))

Rowv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(as.dist(1 - cor(t(matrix_updated), use = "pairwise.complete.obs", method = "pearson")),method='ward.D2')),k= 4), hang_height = 0.1)
Colv  <- hang.dendrogram(color_branches(as.dendrogram(hclust(as.dist(1 - cor(matrix_updated, use = "pairwise.complete.obs", method = "pearson")),method='ward.D2')),k= length(label)), hang_height = 0.1)
jpeg(paste0("/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Plots/Heatmaps/Heatmap_Pheno_",
            'Scaled','_',gsub("PCAresults_","",gsub(".rds","",prs_matrix)),'_','Agglomerative','_','Correlation','_','Ward',".jpg"),
     width = 2000, height = 2000)
par(mar=c(5,5,5,5),cex=1,font=3)
heatmap(as.matrix(matrix_updated), Rowv = Rowv, Colv = Colv,
        scale = "none",
        col= colorRampPalette(brewer.pal(8, "Oranges"))(25))
legend(x="topleft", legend=c("Inevrsely Correlated", "No correlation", "Correled"), 
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
mtext(paste0('Heatmap of PRSs at phenotypes level',"\n",
             dim(matrix_updated)[1],' individuals and ',dim(matrix_updated)[2],' phenotypes',"\n",
             'at P-value threshold of ',
             unlist(strsplit(prs_matrix, "[_.]"))[4]), side = 3, line = 1, cex = 2)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------------
#############################################################################################
################################   |                        |   #############################                                    
################################   |      PCA Analysis      |   #############################
################################   |                        |   #############################
#############################################################################################
#-----------------------------------------------------------------------------------------------------------------------------------------

# First, we make the phenotype labels smaller:
names <- names(matrix_updated)
names1 <-gsub("\\_n[0-9]+","",names)
names1 <-gsub("_and_","_",names1)
names1 <-gsub("_or_","_",names1)
names1 <-gsub("_of_","_",names1)
names1 <-gsub("_the_","_",names1)
names1 <-gsub("_to_","_",names1)
names(matrix_updated) <- names1
print(paste("Running PCA"))
respca <- prcomp(matrix_updated)

# Calculate subject and variable contributions to each component
print(paste("Calculating variable and subject contributions to components"))
res.var <- get_pca_var(respca)
res.ind <- get_pca_ind(respca)

res.pcs <- as.data.frame(respca$x)
res.pcs$IID <- results.S[[prs_matrix]]$residuals$IID
pheno_pcs <- merge(res.pcs,ROSmaster,by="IID")
pheno_pcs2 <- merge(pheno_pcs,geno_pcs,by="IID")
## identify relationships of PCs with phenotypes of interest
pheno.list <- c("amyloid_sqrt","tangles_sqrt","cogn_global_random_slope","age_death","parksc_lv_sqrt",
                "neuroticism_12","agreeableness","cesdsum_lv","educ","pput3123","conscientiousness",
                "extraversion","openness","dxpark","alcohol_g_bl","smoking",
                "hypertension_bl","diabetes_sr_rx_bl","heart_bl")

PClist <- paste0("PC",seq(1:10))
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
ggplot(data=assocres, aes(x=pheno,y=-log10(p),group=pc))+
                geom_bar(stat="identity",position="dodge",aes(fill = colour))+
                geom_hline(aes(yintercept = -log10(0.05/nrow(assocres)),color="P-value < 0.00026"),lty=2)+ 
                geom_hline(aes(yintercept = -log10(0.05),color="P-value < 0.05"),lty=2)+ 
                scale_linetype_manual(name = "limit", values = c(2, 2), 
                                      guide = guide_legend(override.aes = list(color = c("green", "orange"))))+
                geom_text_repel(data=subset(assocres,p<0.05),aes(label=pc))+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -45, hjust = 0))+
                ggtitle(paste0('PRSs at P-value threshhold of',' 1e-05'))

# Proportion of variance explained by the first 500 PCs:
sum(summary(respca)$importance[2,1:500]) #0.81475

# First 100 PCs:
fviz_eig(respca,ncp = 100,
         choice = "variance",  # "eigenvalue"
         geom = "line",
         barfill = "lightblue",
         barcolor = "steelblue",
         main='Proportion of variance explained by the first 100 PCs: 26.755%')

###
##### PC1 ---> autoimmune dieases/ Class II HLA #### Check MHC region HLA 27, HLA B27, HLA DQB1/DRB1
###

sum(summary(respca)$importance[2,1]) # 0.00631
pcnum <- 1 
### get contributing PRS for given PC1
pc1_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc1_top20$phenotypes <- rownames(pc1_top20)
names(pc1_top20) <- c("Cos2","Contribution","phenotypes")
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc1_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = 90,hjust=1,vjust=1))+
                ggtitle(paste0('Quality of represenation of top 20 variables','\n', 'to PC',pcnum))+
                coord_flip()

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc1_top20, aes(x=reorder(phenotypes, +Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = 90,hjust=1,vjust=1))+
                ggtitle(paste0('Quality of represenation of top 20 variables','\n', 'to PC',pcnum))+
                coord_flip()


# CA_Medical_information_h0.13412 ----> Doctor diagnosed sarcoidosis


###
##### PC2 ---> arrhythmia PC
###

sum(summary(respca)$importance[2,2]) # 0.0061
pcnum <- 2
### get contributing PRS for given PC2
pc2_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc2_top20$phenotypes <- rownames(pc2_top20)
names(pc2_top20) <- c("Cos2","Contribution","phenotypes")
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc2_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc2_top20, aes(x=reorder(phenotypes, -Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables to PC',pcnum,' 1e-05'))


###
##### PC3 ---> Skin cancer PC
###

sum(summary(respca)$importance[2,3]) # 0.00511
pcnum <- 3
### get contributing PRS for given PC3
pc3_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc3_top20$phenotypes <- rownames(pc3_top20)
names(pc3_top20) <- c("Cos2","Contribution","phenotypes")
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc3_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc3_top20, aes(x=reorder(phenotypes, -Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables to PC',pcnum,' 1e-05'))


###
##### PC4 ---> ??
###

sum(summary(respca)$importance[2,4]) # 0.00391
pcnum <- 4
### get contributing PRS for given PC4
pc4_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc4_top20$phenotypes <- rownames(pc4_top20)
names(pc4_top20) <- c("Cos2","Contribution","phenotypes")
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc4_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc4_top20, aes(x=reorder(phenotypes, -Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables to PC',pcnum,' 1e-05'))


# CA_Administration_h0.07447 ----> Hepatobiliary & pancreatic surgery

###
##### PC5 ---> AD/Dementia PC
###

sum(summary(respca)$importance[2,5]) # 0.00373
pcnum <- 5 
### get contributing PRS for given PC5
pc5_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc5_top20$phenotypes <- rownames(pc5_top20)
names(pc5_top20) <- c("Cos2","Contribution","phenotypes")
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc5_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = 45,hjust=1,vjust=1))+
                ggtitle(paste0('Quality of represenation of top 20 variables','\n', 'to PC',pcnum))+
                coord_flip()

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc5_top20, aes(x=reorder(phenotypes, +Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum))+
                coord_flip()


# CA_Family_history_h0.10308 ----> Alzheimer's disease/dementia
# CA_Administration_h0.07773 ----> Old age psychiatry


###
##### PC6 ---> Lung/respiratory PC
###

sum(summary(respca)$importance[2,6]) # 0.00358
pcnum <- 6
### get contributing PRS for given PC6
pc6_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc6_top20$phenotypes <- rownames(pc6_top20)
names(pc6_top20) <- c("Cos2","Contribution","phenotypes")
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc6_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc6_top20, aes(x=reorder(phenotypes, -Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables to PC',pcnum,' 1e-05'))



# CA_Administration_h0.07447 ----> Hepatobiliary & pancreatic surgery


###
##### PC7 ---> Brain PC
###

sum(summary(respca)$importance[2,7]) # 0.00347
pcnum <- 7 
### get contributing PRS for given PC7
pc7_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc7_top20$phenotypes <- rownames(pc7_top20)
names(pc7_top20) <- c("Cos2","Contribution","phenotypes")#
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc7_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc7_top20, aes(x=reorder(phenotypes, +Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables','\n', 'to PC',pcnum))+
                coord_flip()

# CA_Administration_h0.05791 ----> Hepatology
# CA_Administration_h0.07447 ----> Hepatobiliary & pancreatic surgery


###
##### PC8 
###

sum(summary(respca)$importance[2,8]) # 0.00343
pcnum <- 8
### get contributing PRS for given PC8
pc8_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc8_top20$phenotypes <- rownames(pc8_top20)
names(pc8_top20) <- c("Cos2","Contribution","phenotypes")
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc8_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc8_top20, aes(x=reorder(phenotypes, -Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables to PC',pcnum,' 1e-05'))

###
##### PC9
###

sum(summary(respca)$importance[2,9]) # 0.00326
pcnum <- 9
### get contributing PRS for given PC9
pc9_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc9_top20$phenotypes <- rownames(pc9_top20)
names(pc9_top20) <- c("Cos2","Contribution","phenotypes")
# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc9_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc9_top20, aes(x=reorder(phenotypes, -Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables to PC',pcnum,' 1e-05'))

###
##### PC10 ---> Colon PC
###

sum(summary(respca)$importance[2,10]) # 0.00318
pcnum <- 10
### get contributing PRS for given PC10
pc10_top20 <- as.data.frame(cbind(res.var$cos2[,pcnum][order(res.var$cos2[,pcnum],decreasing = T)][1:20],res.var$contrib[,pcnum][order(res.var$contrib[,pcnum],decreasing = T)][1:20]))
pc10_top20$phenotypes <- rownames(pc10_top20)
names(pc10_top20) <- c("Cos2","Contribution","phenotypes")
# CA_Administration_h0.05791 ----> Hepatology
# represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.

# cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
ggplot(pc10_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))

# contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
ggplot(pc10_top20, aes(x=reorder(phenotypes, -Contribution), y=Contribution)) +
                geom_bar(stat="identity",position="dodge")+
                theme_minimal()+
                theme(axis.text.x=element_text(angle = -90, hjust = 0))+
                ggtitle(paste0('Contributions (%) of top 20 variables to PC',pcnum,' 1e-05'))


# Contribution of PC 1 to 10
fviz_contrib(respca, choice = "var", axes = 1:10,top=20)
fviz_cos2(respca, choice = "var", axes = 1:10,top=20)
# CA_Medical_information_h0.13412_n509 ----> Doctor diagnosed sarcoidosis

# pc1_10_top20 <- as.data.frame(cbind(res.var$cos2[1:20,1:10][order(res.var$cos2[1:20,1:10],decreasing = T)],res.var$contrib[1:20,1:10][order(res.var$contrib[1:20,1:10],decreasing = T)]))
# pc1_10_top20$phenotypes <- rownames(pc1_10_top20)
# names(pc1_10_top20) <- c("Cos2","Contribution","phenotypes")
# # CA_Administration_h0.05791 ----> Hepatology
# # represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
# 
# # cos2 : represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
# ggplot(pc1_10_top20, aes(x=reorder(phenotypes, -Cos2), y=Cos2)) +
#                 geom_bar(stat="identity",position="dodge")+
#                 theme_minimal()+
#                 theme(axis.text.x=element_text(angle = -90, hjust = 0))+
#                 ggtitle(paste0('Quality of represenation of top 20 variables to PC',pcnum,' 1e-05'))
# 
# # contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
# ggplot(pc1_10_top20, aes(x=reorder(phenotypes, -Contribution), y=Contribution)) +
#                 geom_bar(stat="identity",position="dodge")+
#                 theme_minimal()+
#                 theme(axis.text.x=element_text(angle = -90, hjust = 0))+
#                 ggtitle(paste0('Contributions (%) of top 20 variables to PC',pcnum,' 1e-05'))
# 


# META-PRS
pca_dat <- as.data.frame(cbind(respca$x))
pca_dat$IID <- results.S[[prs_matrix]]$residuals$IID
merged_dat <- merge(ROSmaster,pca_dat,by="IID")
merged_dat <- merge(merged_dat,results.S[[prs_matrix]]$residuals)


############
pc5_mod <- lm(merged_dat$cogn_global_random_slope~merged_dat$PC5)
summary(pc5_mod) # R2 = 0.03866, P-val < 2.2e-16
prs_mod <- lm(merged_dat$cogn_global_random_slope~merged_dat$`IC_Alzheimer's_disease_h0.14431_n910`)
summary(prs_mod) # R2 = 0.04586, P-val < 2e-16


pc5_7_mod <- lm(merged_dat$cogn_global_random_slope~merged_dat$PC5+merged_dat$PC7)
summary(pc5_7_mod) # R2 = 0.04234, P-val-PC5 < 2.2e-16, P-val-PC7=0.00671

pc1.10_mod <- lm(merged_dat$cogn_global_random_slope~merged_dat$PC1+merged_dat$PC2+merged_dat$PC3+
                             merged_dat$PC4+merged_dat$PC5+merged_dat$PC6+merged_dat$PC7+
                             merged_dat$PC8+merged_dat$PC9+merged_dat$PC10)
summary(pc1.10_mod) # R2 = 0.05094

#############
pc5_mod <- lm(merged_dat$amyloid_sqrt~merged_dat$PC5)
summary(pc5_mod) # R2 = 0.05296, P-val < 2.2e-16
prs_mod <- lm(merged_dat$amyloid_sqrt~merged_dat$`IC_Alzheimer's_disease_h0.14431_n910`)
summary(prs_mod) # R2 = 0.06859, P-val < 2e-16

pc5_7_mod <- lm(merged_dat$amyloid_sqrt~merged_dat$PC5+merged_dat$PC7)
summary(pc5_7_mod) # R2 = 0.06117, P-val-PC5 < 2.2e-16, P-val-PC7=0.00101

pc1.10_mod <- lm(merged_dat$amyloid_sqrt~merged_dat$PC1+merged_dat$PC2+merged_dat$PC3+
                                 merged_dat$PC4+merged_dat$PC5+merged_dat$PC6+merged_dat$PC7+
                                 merged_dat$PC8+merged_dat$PC9+merged_dat$PC10)
summary(pc1.10_mod) # R2 = 0.07712
#############
pc5_mod <- lm(merged_dat$tangles_sqrt~merged_dat$PC5)
summary(pc5_mod) # R2 = 0.02876, P-val = 1.51e-09
prs_mod <- lm(merged_dat$tangles_sqrt~merged_dat$`IC_Alzheimer's_disease_h0.14431_n910`)
summary(prs_mod) # R2 = 0.04336, P-val < 9.458e-14

pc5_7_mod <- lm(merged_dat$tangles_sqrt~merged_dat$PC5+merged_dat$PC7)
summary(pc5_7_mod) # R2 = 0.03534, P-val-PC5 < 1.3e-09, P-val-PC7=0.00356

pc1.10_mod <- lm(merged_dat$tangles_sqrt~merged_dat$PC1+merged_dat$PC2+merged_dat$PC3+
                                 merged_dat$PC4+merged_dat$PC5+merged_dat$PC6+merged_dat$PC7+
                                 merged_dat$PC8+merged_dat$PC9+merged_dat$PC10)
summary(pc1.10_mod) # R2 = 0.04929
