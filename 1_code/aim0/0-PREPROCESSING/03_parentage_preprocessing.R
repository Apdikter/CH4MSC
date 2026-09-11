# README!

# This script is for detecting pedigree inconsistencies based on pairwise opposing homozygote (OH) counts.
  # Parents are considered 'inconsistent' if OH counts exceed expectations under a priori error rates, and either 'resolved' or set to missing (0).
    # Candidate parents are those whom satisfy opposing_homozygote_count < error_rate & a year-of-birth lower than that of the individual in question.
      # If >2 candidate parents exist for a single individual, an error is thrown and that animal is not resolved. This is intentional.

# This script is also for detecting line-assignment inconsistencies based on z-scored PC1 distance to to their respective group mean.
  # Individuals with z scores >2 SD are considered 'inconsistent' and are set to missing (NA) for later resolution through the pedigree.

# Crossbred individuals (born of a high- and low-line animal) and all their offspring are also removed as part of this process.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                   PART 1: PEDIGREE QUALITY CONTROL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library("data.table")
library("ggplot2")

cat("Importing genotype and pedigree data...\n\n")

sil_geno <- fread("/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed/genotype_mat.raw", 
                  nrows = Inf, fill = TRUE, strip.white = TRUE, na.strings = c("NA", ""), 
                  drop = c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))

sil_ped <- fread("/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed/genotype_mat.raw", 
                 nrows = Inf, fill = TRUE, strip.white = TRUE, na.strings = c("NA", ""), 
                 select = c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))

stopifnot("Number of records in sil_ped and sil_geno disagree!\n\n" = nrow(sil_geno) == nrow(sil_ped))



cat("Checking integrity of raw genotype matrix...\n\n")
monomorphic_snps <- unname(which(apply(sil_geno, 2, sd, na.rm = TRUE) == 0))
cat(length(monomorphic_snps), "monomorphic SNPs detected.\nRemoving...\n\n")
if (length(monomorphic_snps) > 0) sil_geno[, (monomorphic_snps) := NULL]



cat("Converting to atomic (numeric) matrix...\n\n")
geno_matrix <- as.matrix(sil_geno)
rownames(geno_matrix) <- sil_ped$IID



cat("Generating Opposing Homozygote (OH) matrix...\n\n") # logic from hsphase (Fordosi et al., 2014).
geno_matrix[is.na(geno_matrix)] <- 9 

cm <- geno_matrix - (floor(geno_matrix/9) * 8)
result <- (floor(cm/2) %*% t(round(((cm) - 2)/2))) * (-1)
ohg_matrix <- t(result) + result

stopifnot("Number of records in sil_ped and ohg_mat disagree!\n" = nrow(ohg_matrix) == nrow(sil_ped),
          "Order of records in sil_ped and ohg_mat disagree!\n\n" = all(rownames(ohg_matrix) == sil_ped$IID))



cat("Verifying parentage...\n\n")
err_rate <- ncol(geno_matrix) * (1 - 0.998) # threshold based on a priori knowledge of error rates.

n_sire <- n_sire_confirmed <- 0
n_dam <- n_dam_confirmed <- 0

sire_results <- vector("list", nrow(sil_ped))
dam_results <- vector("list", nrow(sil_ped))

for (i in seq_len(nrow(sil_ped))) {
  
  id <- sil_ped[i, IID]
  sire <- sil_ped[i, PAT]
  dam <- sil_ped[i, MAT]
  
  if (sire != "0") {
    n_sire <- n_sire + 1
    ohg_sire <- ohg_matrix[id, sire]
    sire_ok <- ohg_sire < err_rate
    
      if (sire_ok) {
        n_sire_confirmed <- n_sire_confirmed + 1
      } else sire_results[[i]] <- data.table(IID = id, PAT = sire, OHG_PAT = ohg_sire, STATUS_PAT = "inconsistent")
    } else sire_results[[i]] <- data.table(IID = id, PAT = 0, OHG_PAT = NA, STATUS_PAT = "unassigned")
  
  
  if(dam != "0") {
    n_dam <- n_dam + 1
    ohg_dam <- ohg_matrix[id, dam]
    dam_ok  <- ohg_dam < err_rate

      if(dam_ok) {
        n_dam_confirmed  <- n_dam_confirmed + 1  
      } else dam_results[[i]] <- data.table(IID = id, MAT = dam, OHG_MAT = ohg_dam, STATUS_MAT = "inconsistent")
    } else dam_results[[i]] <- data.table(IID = id, MAT = 0, OHG_MAT = NA, STATUS_MAT = "unassigned") 
}

sire_flags <- rbindlist(sire_results)
dam_flags <- rbindlist(dam_results)

sire_accuracy <- n_sire_confirmed / n_sire
dam_accuracy <- n_dam_confirmed / n_dam

cat("Sire accuracy rate:", sire_accuracy, "\n")
cat("Dam accuracy rate:", dam_accuracy, "\n\n")

cat(sire_flags[STATUS_PAT == "inconsistent", .N] + dam_flags[STATUS_MAT == "inconsistent", .N], "inconsistent parents identified.\n")
cat(sire_flags[STATUS_PAT == "unassigned", .N] + dam_flags[STATUS_MAT == "unassigned", .N], "unassigned parents identified.\n\n")



cat("Checking if any parent-offspring relationships are resolvable...\n\n")
top_n_candidates <- 10  # candidates considered per flagged ID.
max_candidates <- 2   # max 'resolvable' candidates before manual review.

flagged_ids <- union(sire_flags$IID, dam_flags$IID)
flagged_ohg_matrix <- ohg_matrix[flagged_ids, ]

flagged_results <- vector("list", nrow(flagged_ohg_matrix))

for (i in seq_len(nrow(flagged_ohg_matrix))) { # ugly code
  
    id    <- rownames(flagged_ohg_matrix)[i]
    row   <- flagged_ohg_matrix[id, ]
    row   <- row[names(row) != id]
    top_n <- sort(row)[seq_len(top_n_candidates)]
  
  flagged_results[[i]] <- data.table(
    IID           = id,
    YOB_IID       = as.numeric(sub("\\..*", "", id)),
    CANDIDATE     = names(top_n),
    YOB_CANDIDATE = as.numeric(sub("\\..*", "", names(top_n))),
    OHG           = top_n)
}


flagged_results <- rbindlist(flagged_results)

resolvable <- flagged_results[(OHG < err_rate) & (YOB_IID > YOB_CANDIDATE)] # ewes lamb once per year: parent YOB > offspring YOB. 
resolvable[, CANDIDATE_SEX := sil_ped[match(CANDIDATE, IID), SEX]] 



id_counts <- resolvable[, .N, by = IID][N > max_candidates] 

if (nrow(id_counts) > 0) {
  cat("Warning! The following records have >", max_candidates, "potential parents and require manual review:\n", id_counts$IID)
  resolvable <- resolvable[IID %notin% id_counts$IID]
  cat("Removing...\n\n")
}



cat("Resolving relationships...\n\n")
sil_ped[resolvable[CANDIDATE_SEX == 1], PAT := i.CANDIDATE, on = "IID"]
unresolvable_sire_ids <- setdiff(sire_flags$IID, resolvable[CANDIDATE_SEX == 1, IID])
sil_ped[IID %in% unresolvable_sire_ids, PAT := 0]


sil_ped[resolvable[CANDIDATE_SEX == 2], MAT := i.CANDIDATE, on = "IID"]
unresolvable_dam_ids <- setdiff(dam_flags$IID,  resolvable[CANDIDATE_SEX == 2, IID])
sil_ped[IID %in% unresolvable_dam_ids, MAT := 0]



rm(list = setdiff(ls(), c("geno_matrix", "sil_ped"))); invisible(gc())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                       PART 2: OUTLIER DETECTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat("Detecting outliers in each selection line...\n\n")

stopifnot("Order of records in pedigree and genotype matrix disagree.\n" = all(rownames(geno_matrix) == sil_ped$IID))

pca_results <- prcomp(geno_matrix, center = TRUE, scale. = TRUE)

pca_results <- data.table(PC1 = pca_results$x[, 1], 
                          FID = sil_ped$FID,
                          IID = sil_ped$IID,
                          YOB = as.numeric(sub("\\..*", "", sil_ped$IID)))

pca_results[, abs_dist := abs(PC1 - mean(PC1)), by = FID]
pca_results[, z := abs_dist / sd(abs_dist), by = FID]

# If an animal exhibits extreme departure (>2 SD) from its group's mean PC1, FID is set to missing (NA).
  # These animals are reassigned in subsequent steps based on information from the resolved pedigree.

outliers <- pca_results[z >= 2, IID]

sil_ped[IID %chin% outliers, FID := NA]



rm(list = setdiff(ls(), c("geno_matrix", "sil_ped"))); invisible(gc())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                            PART 3: PEDIGREE SORTING & 'IMPUTATION'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat("Topologically sorting pedigree...\n\n") # Logic from Kahn, 1962

done <- rep(FALSE, nrow(sil_ped))
out <- integer(0)

while(length(out) < nrow(sil_ped)) {
  
  sire_ready <- sil_ped$PAT == "0" | sil_ped$PAT %chin% sil_ped$IID[done]
  dam_ready  <- sil_ped$MAT == "0" | sil_ped$MAT %chin% sil_ped$IID[done]
  
  ready <- which(!done & sire_ready & dam_ready)
  
  if(length(ready) == 0) stop("Pedigree contains loops or missing parents\n\n")
  out <- c(out, ready)
  done[ready] <- TRUE
}

sil_ped_sorted <- sil_ped[out]



rm(list = setdiff(ls(), c("geno_matrix", "sil_ped_sorted"))); invisible(gc())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                        PART 4: LINE 'IMPUTATION'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat("Imputation of line assignments...\n\n")

conflicts <- vector("character")
repeat {
  
  missing_before <- sum(is.na(sil_ped_sorted$FID))
  
  for (i in sil_ped_sorted$IID) {
    
    offspring_line <- sil_ped_sorted[PAT == i | MAT == i, FID]
    offspring_line <- unique(na.omit(offspring_line))
    
    sire <- sil_ped_sorted[IID == i, PAT]
    if (sire != "0") {
      sire_line <- sil_ped_sorted[IID == sire, FID]
      sire_offspring_line <- sil_ped_sorted[PAT == sire, FID]
    }   else sire_line <- sire_offspring_line <- NA
     
    dam  <- sil_ped_sorted[IID == i, MAT]
    if (dam  != "0") {
      dam_line <- sil_ped_sorted[IID == dam, FID]
      dam_offspring_line <- sil_ped_sorted[MAT == dam, FID]
    }   else dam_line <- dam_offspring_line <- NA
  
    if (!is.na(sire_line) && !is.na(dam_line) && sire_line != dam_line) {
      conflicts <- c(conflicts, i)
      next  
    }
    
    all_line <- unique(na.omit(c(offspring_line, sire_line, sire_offspring_line, dam_line, dam_offspring_line)))
    if (length(all_line) == 1) sil_ped_sorted[IID == i, FID := ..all_line]
    
  }
  
  missing_after <- sum(is.na(sil_ped_sorted$FID))
  if (missing_after == missing_before) break
  
}

conflicts <- unique(conflicts)
cat(sum(is.na(sil_ped_sorted$FID)), "animals with unresolved lines remain...\n\n")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                         PART 5: CROSSBRED DETECTION & RECURSIVE REMOVAL
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat("Detecting sire-dam-offspring line conflicts...\n\n")

for (i in seq_len(nrow(sil_ped_sorted))) {
  
  id   <- sil_ped_sorted$IID[i]
  fid  <- sil_ped_sorted$FID[i]
  
  offspring_fid <- sil_ped_sorted[PAT == id | MAT == id, FID]
  
  
  sire <- sil_ped_sorted$PAT[i]
    
  if (sire != "0") {
    sire_fid <- sil_ped_sorted[IID == sire, FID]
    sire_offspring_fid <- sil_ped_sorted[PAT == sire | MAT == sire, FID]
  }   else sire_fid <- sire_offspring_fid <- NA
  
  
  dam  <- sil_ped_sorted$MAT[i]
  
  if (dam != "0") {
    dam_fid <- sil_ped_sorted[IID == dam, FID]
    dam_offspring_fid <- sil_ped_sorted[PAT == dam | MAT == dam, FID]
  }   else dam_fid <- dam_offspring_fid <- NA
  
  
  
  all_fid <- na.omit(c(fid, offspring_fid, sire_fid, sire_offspring_fid, dam_fid, dam_offspring_fid))
  
  if (length(unique(all_fid)) > 1) {
    conflicts <- c(conflicts, id)
  }
}

conflicts <- unique(conflicts)
cat(length(conflicts), "conflicts detected.\n")



if (length(conflicts) > 0) { 
  cat("Extracting offspring of crossbred records...\n")
  repeat {
    initial_length <- length(conflicts)
    conflicts <- unique(c(conflicts, sil_ped_sorted[PAT %chin% conflicts, IID], sil_ped_sorted[MAT %chin% conflicts, IID]))
    # Recursively extracts offspring of crossbred individuals and appends to a growing vector.
    # Loop breaks when vector stops growing, e.g., all descendants of 'crossbred' individuals are accounted for.
    
    final_length <- length(conflicts)
    if(initial_length == final_length) break
  }
  
  sil_ped_sorted <- sil_ped_sorted[IID %notin% conflicts]
  cat(length(conflicts), "crossbred records (and descendants) removed.\n\n")
}



rm(list = setdiff(ls(), c("geno_matrix", "sil_ped_sorted"))); invisible(gc())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                       PART 5: WRITING OUT FILES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

fwrite(sil_ped_sorted[!is.na(FID)], file = "/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed/ALL_qc1_temp.fam",
       row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)


fwrite(sil_ped_sorted[!is.na(FID), .(IID)], file = "/mnt/gpfs/persist/projects/2026_morrisonz_masters/0_data/02_processed/keep_ids.txt",
       row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                            PART 6: ILLUSTRATION OF POPULATION STRUCTURE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

library("ggnewscale")

cat("Performing PCA...\n\n")

stopifnot("Pedigree ASRID(s) missing from genotype matrix" = all(sil_ped_sorted$IID %in% rownames(geno_matrix)))
geno_matrix <- geno_matrix[match(sil_ped_sorted$IID, rownames(geno_matrix)), ]

pca <- prcomp(geno_matrix, center = TRUE, scale. = TRUE)

pca_results <- data.table(PC1 = pca$x[, 1], 
                          PC2 = pca$x[, 2], 
                          Line = sil_ped_sorted$FID, ID = sil_ped_sorted$IID, 
                          YOB = as.numeric(sub("\\..*", "", sil_ped_sorted$IID)))


cat("Plotting PC1 VS PC2...\n\n")
ggplot() +
  geom_point(data = pca_results[Line == "high"], aes(x = PC1, y = PC2, color = YOB), size = 2, shape = 1, alpha = 1) +
  scale_color_gradient(low = "#ffcccc", high = "#a80000", name = "YOB (High)") +
  new_scale_color() +
  geom_point(data = pca_results[Line == "low"], aes(x = PC1, y = PC2, color = YOB),size = 2, shape = 1, alpha = 1) +
  scale_color_gradient(low = "#cce5ff", high = "#003d80", name = "YOB (Low)") +
  labs(title = "Principal Component Analysis (PCA)", subtitle = "1702 SNP subset shared by all animals in flock 3633.") +
  theme_minimal()

ggsave("/mnt/gpfs/persist/projects/2026_morrisonz_masters/2_results/aim0/0-PREPROCESSING/PCA.png")






q()

cat("Distribution of OH counts by relationship status...\n\n")
ggplot(data = as.data.frame(parent_ohg), aes(x = parent_ohg)) + 
  geom_density(adjust = 50, fill = alpha("#cce5ff", 0.2), colour = alpha("#13265C", 0.8), linewidth = 0.5) +
  annotate("segment", x = 16, xend = 16, y = 0, yend = Inf, linewidth = 0.5, linetype = 2, colour = "#a80000") + 
  annotate("text", x = 20, y = 0.51, label = "16 SNP Threshold", size = 3.5, colour = "#a80000") +
  labs(title = "Density Distribution of Opposing Homozygote (OH) Counts",
       subtitle = "Parent-Offspring Relationships", 
       x = "Number of OH loci", y = "") +
  theme_minimal() + xlim(0, 100)


