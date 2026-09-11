#!/bin/bash
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%j.out

set -eou pipefail

#############################################################################################################

echo -e "Initializing file paths & importing relevant modules...\n"

DIR="/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed"
RES="/mnt/gpfs/persist/projects/2026_morrisonz_masters/2_results/aim0/0-PREPROCESSING"
VCF="/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed/beagle_vcf"

module load PLINK/2.00-20250830-avx2

#############################################################################################################

# PLINK2 snapshot is bugged and doesn't permit filtering only on IID. This is a temporary fix for that.
join -1 1 -2 2 -o 2.1,2.2 <(tail -n +2 ${DIR}/keep_ids.txt | sort) <(sort -k2,2 ${DIR}/ALL.fam) > ${DIR}/keep_ids_fidiid.txt

echo -e "Creating pedigree-corrected PLINK files....\n"

plink2 --bfile ${DIR}/ALL \
       --sheep \
       --keep ${DIR}/keep_ids_fidiid.txt \
       --make-bed \
       --out ${DIR}/ALL_qc1

rm ${DIR}/ALL_qc1.{fam,log}
cp ${DIR}/ALL_qc1_temp.fam ${DIR}/ALL_qc1.fam # replaces .fam with resolved pedigree from 03_parentage_preprocessing.

#############################################################################################################
#
#echo -e "Separating BFILE into sample and reference sets...\n"
#
#plink2 --bfile ${DIR}/ALL_qc1 \
#       --sheep \
#       --keep ${DIR}/hd_genotyped_ids.txt \
#       --make-bed \
#       --out ${DIR}/reference
#
##############################################################################################################
#
#echo -e "Calculating genome-wide FST between reference sets...\n"
#
## Specifying factor levels for pairwise FST (FID, IID, LVL)
#awk '{print $1, $2, $1}' ${DIR}/reference.fam > ${DIR}/within_fid.txt
#
#plink2 --bfile ${DIR}/reference \
#       --sheep --nonfounders \
#       --within ${DIR}/within_fid.txt \
#       --fst CATPHENO method=wc \
#       --out ${RES}/reference \
#
#rm ${DIR}/within_fid.txt
#rm ${DIR}/reference.{bed,bim,fam}
#rm ${RES}/reference.{fst.summary,log}
#
#
############################################################################################################

echo -e "partitioning by chr, and converting to VCF format...\n"

for CHR in {1..26}; do

plink2 --bfile ${DIR}/ALL_qc1 \
       --sheep \
       --keep-fam <(echo "low") \
       --chr ${CHR} \
       --export vcf \
       --out ${VCF}/samples_low_chr${CHR}

plink2 --bfile ${DIR}/ALL_qc1 \
       --sheep \
       --keep-fam <(echo "high") \
       --chr ${CHR} \
       --export vcf \
       --out ${VCF}/samples_high_chr${CHR}

plink2 --bfile ${DIR}/ALL_qc1 \
       --sheep \
       --chr ${CHR} \
       --export vcf \
       --out ${VCF}/samples_all_chr${CHR}

done

rm ${VCF}/*.log
