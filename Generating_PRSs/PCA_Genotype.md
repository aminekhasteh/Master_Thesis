Generating PCA from the Genotype Dataset
================

## PCA

Population structure is the principal source of confounding in GWAS and
is usually accounted for by incorporating principal components (PCs) as
covariates.

First we perform prunning:

``` bash
#!/bin/bash --login
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=QC
#SBATCH --output QC_pc_1.out.txt

cd /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/QC_Geno/ROSMAP/ROSMAP_PCA/

module load PLINK2/1.90b3.46

plink \
    --bfile /external/rprshnas01/netdata_kcni/dflab/data/rosmap/genotype/TOPmed_imputed/vcf/merged/merged_overlap_rs \
    --extract /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/QC_Geno/ROSMAP/ROSMAP_QC_rsid/ROSMAP.QC.snplist \
    --keep /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/QC_Geno/ROSMAP/ROSMAP_QC_rsid/ROSMAP.QC.valid.sample \
    --indep-pairwise 200 50 0.25 \
    --out ROSMAP.QC1
```

  - 9329439 variants loaded from .bim file.
  - 7779119 variants and 2052 people pass filters and QC.
  - Pruned 89184 variants from chromosome 22, leaving 8622.
  - Pruning complete. 7113162 of 7779119 variants removed.

Now, we make the bed file to perform PCA:

``` bash
#!/bin/bash --login
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=QC
#SBATCH --output QC_pc_2.out.txt

cd /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/QC_Geno/ROSMAP/ROSMAP_PCA/

module load PLINK2/1.90b3.46

plink \
    --bfile /external/rprshnas01/netdata_kcni/dflab/data/rosmap/genotype/TOPmed_imputed/vcf/merged/merged_overlap_rs \
    --extract ROSMAP.QC1.prune.in \
    --keep /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/QC_Geno/ROSMAP/ROSMAP_QC_rsid/ROSMAP.QC.valid.sample \
    --make-bed \
    --out geno_forPCA
```

Calculating the first 10 PCs:

``` bash
#!/bin/bash --login
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=QC
#SBATCH --output QC_pc_3.out.txt

cd /external/rprshnas01/netdata_kcni/dflab/team/ak/Thesis/QC_Geno/ROSMAP/ROSMAP_PCA/

module load PLINK2/1.90b3.46

plink \
    --bfile geno_forPCA \
    --pca 10\
    --out geno_qc
```

  - 665957 variants loaded from .bim file.
  - 665957 variants and 2052 people pass filters and QC.
