# README!

# This script is for extracting relevant pedigree and genotypic data, and for construction of PLINK .bed, .bim, and .fam file formats.
      # SNP(s) that exhibit high rates of inconsistency between arrays are dropped.
      # SNP(s) that exhibit duplicate physical positions are dropped
      # SNP call(s) with impossible values (X > 2 | X < 0) are dropped.

# Relevant file formats are constructed only from individuals that are pedigreed AND genotyped.
      # If either type of data is missing, these records are dropped.

# Metadata is also extracted in the final step (ASRIDs of high- and low-density samples).


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                             PART 1: PEDIGREE & PHENOTYPE DATA CLEANUP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library(data.table)
library(genio)

cat("\nImporting raw pedigree and phenotype data...\n\n")
sil_ped <- fread("/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/01_raw/pedigree.csv",
                 nrows = Inf, fill = TRUE, strip.white = TRUE, na.strings = "")

sil_phen <- fread("/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/01_raw/phenotypes.csv",
                  nrows = Inf, fill = TRUE, strip.white = TRUE, na.strings = "")

# Check to prevent downstream errors when indexing, filtering, or merging on ASRID.
stopifnot("ASRID(s) in sil_ped have duplicate records!\n" = !anyDuplicated(sil_ped$ASRID))



cat("Extracting flock 3633 (and founders) from raw pedigree...\n")
flk3633_info <- sil_ped[FLOCK == 3633, 
                      .(ASRID, 
                        SireASRID, 
                        DamASRID,
                        SireFLK = sub("^[^.]+\\.([^.]+)\\..*$", "\\1", SireASRID),
                        DamFLK  = sub("^[^.]+\\.([^.]+)\\..*$", "\\1", DamASRID)) ]


founding_sire <- unique(flk3633_info[SireFLK != "3633", SireASRID]) 
founding_dam <- unique(flk3633_info[DamFLK != "3633", DamASRID]) 
flk3633_animals <- unique(flk3633_info[, ASRID])

cat(length(flk3633_animals), "animals in flock 3633, and", length(founding_sire) + length(founding_dam), "founders.\n\n")

keep_ids <- c(founding_dam, founding_sire, flk3633_animals)
sil_ped <- sil_ped[ASRID %in% keep_ids]



cat("Merging line and sex phenotypes...\n")
inconsistent_line <- sil_phen[, .(n = uniqueN(line)), by = ASRID][n > 1]
inconsistent_sex <- sil_phen[, .(n = uniqueN(SIL_sex)), by = ASRID][n > 1]

stopifnot("ASRID(s) in sil_phen have more than one distinct 'line' value - match() will pick the first occurrence!\n" = nrow(inconsistent_line) == 0,
          "ASRID(s) in sil_phen have more than one distinct 'SIL_sex' value - match() will pick the first occurrence!\n" = nrow(inconsistent_sex) == 0 )

sil_ped[, line := sil_phen[match(sil_ped$ASRID, sil_phen$ASRID), line]]
cat(sum(is.na(sil_ped$line)), "animals with unresolved lines remain.\n")

sil_ped[, sex := sil_phen[match(sil_ped$ASRID, sil_phen$ASRID), SIL_sex]]
cat(sum(is.na(sil_ped$sex)), "animals with unresolved sex remain.\n\n")



rm(list = setdiff(ls(), "sil_ped")); invisible(gc())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                    PART 2: GENOTYPE DATA CLEANUP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

sil_geno <- fread("/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/01_raw/genotypes.csv",
                  nrows = Inf, fill = TRUE, strip.white = TRUE, na.strings = "")

sil_snps <- fread("/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/01_raw/HDsnps.csv",
                  nrows = Inf, fill = TRUE, strip.white = TRUE, na.strings = "")

sil_merg <- fread("/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/01_raw/InconsistentSNPs.csv", 
                  nrows = Inf, fill = TRUE, strip.white = TRUE, na.strings = "")

# Checks to prevent downstream errors when indexing, filtering, or merging on SNP and/or animal ID.
stopifnot("ASRID(s) in sil_geno have duplicate records!\n" = !anyDuplicated(sil_geno$ASRID),
          "SNP ID(s) in sil_geno have duplicate records!\n" = !anyDuplicated(names(sil_geno)),
          "SNP ID(s) in sil_snps have duplicate records!\n" = !anyDuplicated(sil_snps$SNP),
          "SNP ID(s) in sil_merg have duplicate records!\n" = !anyDuplicated(sil_merg$SNP))



cat("Dropping SNP(s) with inconsistent calls across arrays...\n")
inconsistent_snps <- sil_merg[inconsistentrate >= 0.01, SNP] 
  
if(length(inconsistent_snps) != 0) {
sil_snps <- sil_snps[SNP %notin% inconsistent_snps]

sil_geno[, (inconsistent_snps) := NULL] 
cat(length(inconsistent_snps), "inconsistent SNP(s) removed.\n")
 }



cat("\nDropping SNP(s) with duplicate positions across arrays...\n")
duplicate_snps <- sil_snps[duplicated(sil_snps, by = c("chrom", "pos")), SNP]

if(length(duplicate_snps) != 0) {
sil_snps <- sil_snps[SNP %notin% duplicate_snps]

sil_geno[, (duplicate_snps) := NULL] 
cat(length(duplicate_snps), "duplicate SNP(s) removed.\n") 
 }



cat("\nDropping SNP calls with impossible values...\n")
sil_geno_snpcols <- setdiff(names(sil_geno), "ASRID")

n_changed <- 0

for (col in sil_geno_snpcols) {
  x <- sil_geno[[col]]
  bad <- x > 2L | x < 0L
  nb <- sum(bad, na.rm = TRUE)
  
  if (nb > 0) {n_changed <- n_changed + nb
  x[bad] <- NA_integer_
  set(sil_geno, j = col, value = x)
  }
  
}

cat(n_changed, "calls removed.\n\n")



rm(list = setdiff(ls(), c("sil_ped", "sil_snps", "sil_geno", "sil_geno_snpcols"))); invisible(gc())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                            PART 3: CREATING PLINK .BIM and .BED FILES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat("Assembling .bim file from HD array positions...\n")
all_bimfile <- sil_snps[, .(chr = chrom, id = SNP, posg = 0L, pos = pos, alt = "B", ref = "A")]

all_bim_snpcols <- all_bimfile$id

intersect_snps <- intersect(sil_geno_snpcols, all_bim_snpcols)
sil_only_snps <- setdiff(sil_geno_snpcols, all_bim_snpcols)
bim_only_snps <- setdiff(all_bim_snpcols, sil_geno_snpcols)

cat(length(intersect_snps), "SNP(s) shared between .bim file and genotype matrix.\n")
all_bimfile <- all_bimfile[id %chin% intersect_snps]

cat(length(sil_only_snps), "SNP(s) present in genotype matrix but not .bim file (dropped from genotype matrix).\n")
cat(length(bim_only_snps), "SNP(s) present in .bim file but not genotype matrix (dropped from .bim file).\n\n")
sil_geno[, (sil_only_snps) := NULL]



cat("Ordering SNP(s) by physical position...\n")
bim_chr_levels <- c(1:26, "X")
all_bimfile[, chr_order := match(chr, bim_chr_levels)]
setorder(all_bimfile, chr_order, pos)
all_bimfile[, chr_order := NULL]



cat("Forcing identical order .bim file and genotype matrix...\n\n")
sil_geno <- sil_geno[, c("ASRID", all_bimfile$id), with = FALSE]



rm(list = setdiff(ls(), c("sil_ped", "sil_snps", "sil_geno", "all_bimfile"))); invisible(gc())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                   PART 4: CREATING PLINK .BED & .FAM FORMATS (FULL SET)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat("Constructing PLINK .bed and .fam files...\n\n")
all_genotyped_animals <- sil_geno$ASRID
all_pedigreed_animals <- sil_ped$ASRID

all_intersect_animals <- intersect(all_genotyped_animals, all_pedigreed_animals)

all_famfile <- sil_ped[ASRID %chin% all_intersect_animals]
all_bedfile <- sil_geno[ASRID %chin% all_intersect_animals]

all_bedfile_ids <- all_bedfile$ASRID



cat("Matching order between .fam and .bed files...\n\n")
all_famfile <- all_famfile[match(all_bedfile_ids, all_famfile$ASRID)]

stopifnot("Misaligned ASRID(s) remain after re-ordering!\n\n" = identical(all_famfile$ASRID, all_bedfile$ASRID))



cat("Checking for missing parents...\n")
unique_sire <- unique(all_famfile$SireASRID)
unique_dam <- unique(all_famfile$DamASRID)

missing_sire <- unique_sire[unique_sire %notin% all_famfile$ASRID]
missing_dam <- unique_dam[unique_dam %notin% all_famfile$ASRID]

cat(length(missing_sire), "sires and", length(missing_dam), "dams missing.\n\n")

all_famfile[SireASRID %chin% missing_sire, SireASRID := "0"]
all_famfile[DamASRID %chin% missing_dam, DamASRID := "0"]

all_famfile <- all_famfile[, .(fam = line, 
                               id = ASRID, 
                               pat = SireASRID,
                               mat = DamASRID, 
                               sex = sex, 
                               pheno = NA)]



cat("Transposing genotype matrix...\n\n")
all_bedfile[, ASRID := NULL]

X <- as.matrix(transpose(all_bedfile)) 

rownames(X) <- all_bimfile$id
colnames(X) <- all_famfile$id



cat("Writing out PLINK2 file formats...\n")
genio::write_plink(file = "/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed/ALL",
                   X = X, bim = all_bimfile, fam = all_famfile)



rm(list = setdiff(ls(), c("X", "all_famfile", "all_bedfile_ids"))); invisible(gc())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                           PART 5: EXTRACTING METADATA FOR REFERENCE SETS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat("\nExtracting IDs of high-density animals to be included in reference sets...\n\n")
sample_density <- colSums(!is.na(X))

hd_genotyped_animals <- all_bedfile_ids[sample_density >= 500000] # > 500,000 SNP threshold to be defined as HD.
ld_genotyped_animals <- all_bedfile_ids[sample_density <  500000] # < 500,000 SNP is LD.

stopifnot(length(hd_genotyped_animals) + length(ld_genotyped_animals) == length(all_bedfile_ids))

fwrite(all_famfile[id %chin% hd_genotyped_animals, .(fam, id)], file = "/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed/hd_genotyped_ids.txt",
       row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

fwrite(all_famfile[id %chin% ld_genotyped_animals, .(fam, id)], file = "/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed/ld_genotyped_ids.txt",
       row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


