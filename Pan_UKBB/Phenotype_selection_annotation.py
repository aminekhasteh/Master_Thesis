import pandas as pd
import os
import numpy as np

ukbb_manifest = pd.read_excel('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/Pan_UK_Biobank_phenotype_manifest.xlsx',
                                engine='openpyxl')

##############################################################################
#############################   |            |   #############################                                    
#############################   | SELECTION  |   #############################
#############################   |            |   #############################
##############################################################################

# Selecting only European individuals
ukbb_manifest_EUR = ukbb_manifest[ukbb_manifest['pops']=='EUR']

print(ukbb_manifest_EUR.shape) # (4230, 48)

# Replacing any NA value with an empty str
ukbb_manifest_EUR = ukbb_manifest_EUR.replace(np.nan, '', regex=True)

# Selecting phenotypes with heritability scale of bigger or equal to 5%
ukbb_manifest_EUR_h2_05 = ukbb_manifest_EUR[ukbb_manifest_EUR['saige_heritability_EUR']>=0.05]

print(ukbb_manifest_EUR_h2_05.shape) # (1576, 48)

# We also want to select non-gender specific phenotypes:
ukbb_manifest_EUR_h2_05['pheno_sex'].value_counts()

   # both_sexes    1547 --> we only want these!
   # females         21
   # males            8

ukbb_manifest_EUR_h2_05_both_sex = ukbb_manifest_EUR_h2_05.loc[ukbb_manifest_EUR_h2_05['pheno_sex']=='both_sexes']

# However, there are still phenotypes that are gender specific. We can find them by the number of case/control for each gender.
# If the number full cases for each gender is 0, we remoev them:

ukbb_manifest_EUR_h2_05_both_sex = ukbb_manifest_EUR_h2_05_both_sex.loc[(ukbb_manifest_EUR_h2_05_both_sex['n_cases_full_cohort_females']!=0) & 
                                                                        (ukbb_manifest_EUR_h2_05_both_sex['n_cases_full_cohort_males']!=0)]

print(ukbb_manifest_EUR_h2_05_both_sex.shape) # (1327, 48)

ukbb_manifest_EUR_h2_05_both_sex['n_cases_full_cohort_females'].describe()

   # count      1327.00000
   # mean        874.75584
   # std        6182.01862
   # min          17.00000
   # 25%         217.50000
   # 50%         425.00000
   # 75%         766.50000
   # max      219607.00000

ukbb_manifest_EUR_h2_05_both_sex['n_cases_full_cohort_males'].describe()

   # count      1327.000000
   # mean        791.773173
   # std        5238.677193
   # min           1.000000 --> too low, we need to remove phenotypes with low male/female cases
   # 25%         210.000000
   # 50%         379.000000
   # 75%         709.000000
   # max      186645.000000
   
ukbb_manifest_EUR_h2_05_both_sex = ukbb_manifest_EUR_h2_05_both_sex.loc[(ukbb_manifest_EUR_h2_05_both_sex['n_cases_full_cohort_males'] > 40) & 
                                                                        (ukbb_manifest_EUR_h2_05_both_sex['n_cases_full_cohort_females'] > 40)]

print(ukbb_manifest_EUR_h2_05_both_sex.shape) # (1307, 48)

# Saving this file
ukbb_manifest_EUR_h2_05_both_sex.to_csv('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex.csv')

# Manually going through the descriptions before annotations.

##############################################################################
#############################   |            |   #############################                                    
#############################   | ANNOTATION |   #############################
#############################   |            |   #############################
##############################################################################

# importing the new dataset:

# use once you have th final PRS
#ukbb_manifest_main = pd.read_csv('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/ukbb_manifest_EUR_h2_05_both_sex_selected_pheno.csv')
ukbb_manifest_main = pd.read_csv('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/OLD/Pan_UK_Biobank_phenotype_manifest_h2_more_0.05_both_sex_non_ambigous.csv')

###print(ukbb_manifest_main.shape) # (1272, 49)
print(ukbb_manifest_main.shape) # (1305, 49)
ukbb_manifest_main['trait_type'].value_counts()

phenocode_annotate_lst = []
for rows in range(1305):
   print(rows)
   print(row['trait_type'])
   row = ukbb_manifest_main.iloc[rows,:]
   if row['trait_type'] =='phecode':
      phenocode_annotate = 'PH_'+ row.description.replace(' ','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] =='prescriptions': 
      phenocode_annotate = 'PR_'+ row.phenocode.replace(' ','_').replace('|','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] =='icd10':
      phenocode_annotate = 'IC_'+ row.category.split("|",)[-1][5:].replace(' ','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] =='categorical':
      phenocode_annotate = 'CA_'+ row.category.split(">",)[-1].replace(' ','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   if row['trait_type'] =='continuous':
      phenocode_annotate = 'CO_'+ row.category.split(">",)[-1][5:].replace(' ','_') + '_h' + str("{:.5f}".format(row.saige_heritability_EUR)) + '_n' + row.n_cases_full_cohort_both_sexes.astype(str)
   phenocode_annotate_lst = phenocode_annotate_lst +[phenocode_annotate.replace('__','_')]

ukbb_manifest_main['phenocode_annotate_lst'] = phenocode_annotate_lst

# Saving this file
ukbb_manifest_main.to_csv('/Users/amink/OneDrive/Documents/Current Jobs/Masters Thesis/Code/Master_Thesis/Pan_UKBB/OLD/ukbb_manifest_EUR_h2_05_both_sex_annotated.csv')
