#!/bin/bash
#SBATCH --job-name=plink_parentage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%j.out

set -eou pipefail

#############################################################################################################

echo -e "Initializing file paths & importing relevant modules...\n"
DIR="/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed"
module load PLINK/2.00-20250830-avx2

echo -e "Filtering low-quality SNPs...\n"

# Filtering monomorphic SNPs & SNPs with extreme Mendelian inconsistencies.
        # Doesn't throw away any family trios. Only highly inconsistent SNPs.

plink2 --bfile ${DIR}/ALL \
       --sheep --mac 1 --nonfounders \
       --make-bed \
       --out ${DIR}/ALL_qc1

#############################################################################################################

echo -e "Filtering for ubiquitous SNPs...\n"
plink2 --bfile ${DIR}/ALL_qc1 \
       --sheep --nonfounders \
       --geno 0 \
       --not-chr X \
       --export A \
       --out ${DIR}/genotype_mat

#############################################################################################################

rm ${DIR}/ALL_qc1.{bed,fam,bim} # removing intermediate files
rm ${DIR}/*.log # removing plink log files.

