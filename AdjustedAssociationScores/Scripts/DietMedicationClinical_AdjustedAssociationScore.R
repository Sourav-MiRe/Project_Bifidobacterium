# DIET, MEEDICATION AND CLINICAL ADJUSTED ASSOCIATION-SCORE COMPUTATION AND ANALYSIS

library(vegan)
library(ade4)
library(rrcov)
library(ppcor)
library(dplyr)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)
library(tidyr)

load("/data/StrongAssociationScores_new_16s.RData")
load("/data/StrongAssociationScores_new_wgs.RData")

### Functions ###

#to compute PCA and PCoA components 
compute_pca_pcoa <- function(df_temp) {
  set.seed(123) 
  row_ids <- rownames(df_temp)
  
  df_temp <- as.data.frame(apply(df_temp, 2, as.numeric))
  rownames(df_temp) <- row_ids
  df_temp <- na.omit(df_temp)
  
  ## ---- PCoA ----
  dist <- vegdist(df_temp, method = "euclidean")
  pcoa_results <- dudi.pco(dist, scannf = FALSE, nf = 5)
  
  # PCoA eigenvalues & variance explained
  pcoa_eig <- pcoa_results$eig
  pcoa_var_exp <- pcoa_eig / sum(pcoa_eig)
  pcoa_var_first3 <- sum(pcoa_var_exp[1:3])
  
  ## ---- Robust PCA ----
  rob_pca <- PcaHubert(df_temp, k = ncol(df_temp), scale = TRUE, mcd = TRUE)
  
  # PCA eigenvalues & variance explained
  pca_eig <- rob_pca@eigenvalues
  pca_var_exp <- pca_eig / sum(pca_eig)
  pca_var_first3 <- sum(pca_var_exp[1:3])
  
  ## ---- Scores ----
  pca_scores <- as.data.frame(rob_pca@scores[, 1:3])
  colnames(pca_scores) <- paste0("PCA", 1:3)
  
  pcoa_scores <- as.data.frame(pcoa_results$li[, 1:3])
  colnames(pcoa_scores) <- paste0("PCoA", 1:3)
  
  combined <- cbind(pca_scores, pcoa_scores)
  
  ## ---- Return both components + variance explained ----
  attr(combined, "variance_explained") <- list(
    PCA_first3 = pca_var_first3,
    PCoA_first3 = pcoa_var_first3,
    PCA_each = pca_var_exp,
    PCoA_each = pcoa_var_exp
  )
  
  return(combined)
}

extract_variance_summary <- function(obj, study_name) {
  var <- attr(obj, "variance_explained")
  
  data.frame(
    Study = study_name,
    PCA_PC1 = var$PCA_each[1],
    PCA_PC2 = var$PCA_each[2],
    PCA_PC3 = var$PCA_each[3],
    PCA_First3 = var$PCA_first3,
    PCoA_PC1 = var$PCoA_each[1],
    PCoA_PC2 = var$PCoA_each[2],
    PCoA_PC3 = var$PCoA_each[3],
    PCoA_First3 = var$PCoA_first3
  )
}

#compute pcor and make association_df in 16s
compute_partial_corr_associations_16s <- function(spdata,df_association_all_new_16s,bif_name,bif_studies) {
  if (!requireNamespace("ppcor", quietly = TRUE)) {
    stop("Please install the 'ppcor' package: install.packages('ppcor')")
  }
  library(ppcor)
  
  # Filter samples by age category
  studies <- bif_studies #4
  
  non_bif_species <- rownames(df_association_all_new_16s)
  
  corr_mat <- matrix(NA, nrow = length(non_bif_species), ncol = length(studies),
                     dimnames = list(non_bif_species, studies))
  pval_mat <- corr_mat
  dir_mat  <- corr_mat
  
  for (study in studies) {
    study_df <- subset(spdata, study_name == study)
    if (nrow(study_df) < 5) next
    
    for (sp in non_bif_species) {
      if (!all(c(sp, bif_name, "PCA1", "PCA2", "PCA3", "PCoA1", "PCoA2", "PCoA3") %in% colnames(study_df))) next #pca, pcoa 6
      
      x <- study_df[[sp]]
      y <- study_df[[bif_name]]
      covars <- study_df[, c("PCA1", "PCA2", "PCA3", "PCoA1", "PCoA2", "PCoA3")] #pca, pcoa 6
      
      if (sum(complete.cases(x, y, covars)) < 5) next
      
      test_res <- try(ppcor::pcor.test(x, y, covars, method = "spearman"), silent = TRUE)
      if (inherits(test_res, "try-error")) next
      
      corr_mat[sp, study] <- test_res$estimate
      pval_mat[sp, study] <- test_res$p.value
    }
  }
  
  corr_mat[is.na(corr_mat)] <- 0
  pval_mat[is.na(pval_mat)] <- 1
  
  dir_mat[corr_mat > 0 & pval_mat <= 0.05] <- 2
  dir_mat[corr_mat > 0 & pval_mat > 0.05 & pval_mat <= 0.1] <- 1
  dir_mat[corr_mat < 0 & pval_mat <= 0.05] <- -2
  dir_mat[corr_mat < 0 & pval_mat > 0.05 & pval_mat <= 0.1] <- -1
  dir_mat[is.na(dir_mat)] <- 0
  
  positive <- apply(dir_mat, 1, function(x) sum(x == 2, na.rm = TRUE))
  negative <- apply(dir_mat, 1, function(x) sum(x == -2, na.rm = TRUE))
  total <- length(studies)
  
  score <- ((positive - negative) / total) *
    (1 - ((pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
  
  association <- data.frame(
    positive = positive,
    negative = negative,
    total = total,
    score = score,
    row.names = non_bif_species
  )
  
  output_list <- list(
    corr = as.data.frame(corr_mat),
    pval = as.data.frame(pval_mat),
    dir = as.data.frame(dir_mat),
    association = association
  )
  
  return(output_list)
}

#compute pcor and make association_df in wgs
compute_partial_corr_associations_wgs <- function(spdata,df_association_all_new_wgs,bif_name,bif_studies) {
  if (!requireNamespace("ppcor", quietly = TRUE)) {
    stop("Please install the 'ppcor' package: install.packages('ppcor')")
  }
  library(ppcor)
  
  # Filter samples by age category
  studies <- bif_studies #4
  
  non_bif_species <- rownames(df_association_all_new_wgs)
  
  corr_mat <- matrix(NA, nrow = length(non_bif_species), ncol = length(studies),
                     dimnames = list(non_bif_species, studies))
  pval_mat <- corr_mat
  dir_mat  <- corr_mat
  
  for (study in studies) {
    study_df <- subset(spdata, study_name == study)
    if (nrow(study_df) < 5) next
    
    for (sp in non_bif_species) {
      if (!all(c(sp, bif_name, "PCA1", "PCA2", "PCA3", "PCoA1", "PCoA2", "PCoA3") %in% colnames(study_df))) next #pca, pcoa 6
      
      x <- study_df[[sp]]
      y <- study_df[[bif_name]]
      covars <- study_df[, c("PCA1", "PCA2", "PCA3", "PCoA1", "PCoA2", "PCoA3")] #pca, pcoa 6
      
      if (sum(complete.cases(x, y, covars)) < 5) next
      
      test_res <- try(ppcor::pcor.test(x, y, covars, method = "spearman"), silent = TRUE)
      if (inherits(test_res, "try-error")) next
      
      corr_mat[sp, study] <- test_res$estimate
      pval_mat[sp, study] <- test_res$p.value
    }
  }
  
  corr_mat[is.na(corr_mat)] <- 0
  pval_mat[is.na(pval_mat)] <- 1
  
  dir_mat[corr_mat > 0 & pval_mat <= 0.05] <- 2
  dir_mat[corr_mat > 0 & pval_mat > 0.05 & pval_mat <= 0.1] <- 1
  dir_mat[corr_mat < 0 & pval_mat <= 0.05] <- -2
  dir_mat[corr_mat < 0 & pval_mat > 0.05 & pval_mat <= 0.1] <- -1
  dir_mat[is.na(dir_mat)] <- 0
  
  positive <- apply(dir_mat, 1, function(x) sum(x == 2, na.rm = TRUE))
  negative <- apply(dir_mat, 1, function(x) sum(x == -2, na.rm = TRUE))
  total <- length(studies)
  
  score <- ((positive - negative) / total) *
    (1 - ((pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
  
  association <- data.frame(
    positive = positive,
    negative = negative,
    total = total,
    score = score,
    row.names = non_bif_species
  )
  
  output_list <- list(
    corr = as.data.frame(corr_mat),
    pval = as.data.frame(pval_mat),
    dir = as.data.frame(dir_mat),
    association = association
  )
  
  return(output_list)
}


### DIET 16S STUDIES ###

load("/data/CuestaZuluaga_DietData_spProfile.RData")
load("/data/NUAGE_DietData_spProfile.RData")
load("/data/Hernandez_DietData_spProfile.RData")
load("/data/BartonW_DietData_spProfile.RData")

pc_CuestaZuluaga <- compute_pca_pcoa(CuestaZuluaga_metadata_diet)
pc_NUAGE <- compute_pca_pcoa(NUAGE_metadata_diet)
pc_Hernandez <- compute_pca_pcoa(Hernandez_metadata_diet)
pc_BartonW <- compute_pca_pcoa(BartonW_metadata_diet)

summary_table <- rbind(
  extract_variance_summary(pc_CuestaZuluaga, "CuestaZuluaga"),
  extract_variance_summary(pc_NUAGE,         "NUAGE"),
  extract_variance_summary(pc_Hernandez,     "Hernandez"),
  extract_variance_summary(pc_BartonW,       "BartonW")
)

df_diet_16s_variance <- summary_table

df_diet_16s_variance[ , -1] <- round(df_diet_16s_variance[ , -1] * 100, 2)
df_diet_16s_variance

write.xlsx(df_diet_16s_variance, "/results/diet_16s_variance_explained_percentage.xlsx")

#-------------------------------------------------------------------------------------------------
##### Bind PC's + SpProfiles ######

#### 1. CuestaZuluagaJ_2018  #####
identical(rownames(CuestaZuluaga_speciesProfile), rownames(pc_CuestaZuluaga))
CuestaZuluagaJ_combined <- cbind(CuestaZuluaga_speciesProfile, pc_CuestaZuluaga)
#View(CuestaZuluagaJ_combined)

#### 2. NUAGE  #####
identical(rownames(NUAGE_speciesProfile), rownames(pc_NUAGE))
NUAGE_combined <- cbind(NUAGE_speciesProfile, pc_NUAGE)
#View(NUAGE_combined)

#### 3. Hernandez_2023  #####
identical(rownames(hernandez_spProfile), rownames(pc_Hernandez))
Hernandez_combined <- cbind(hernandez_spProfile, pc_Hernandez)
#View(Hernandez_combined)

#### 4. BartonW_2018  #####
identical(rownames(BartonW_speciesProfile), rownames(pc_BartonW))
BartonW_combined <- cbind(BartonW_speciesProfile, pc_BartonW)
#View(BartonW_combined)

##### Overall spdata creation ######
spdata <- bind_rows(CuestaZuluagaJ_combined, NUAGE_combined, Hernandez_combined, BartonW_combined)
dim(spdata)
#View(spdata)

any(is.na(spdata[, !colnames(spdata) %in% "study_name"]))

spdata[, !colnames(spdata) %in% "study_name"] <- 
  replace(spdata[, !colnames(spdata) %in% "study_name"], 
          is.na(spdata[, !colnames(spdata) %in% "study_name"]), 
          0)

sum(rowSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)
sum(colSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)

bif_studies <- unique(spdata$study_name)

#---------------------------------------------------------------------------------------------------------------------------------------
##### compute pcor and make association_df ######

df_all_Bif_detection <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_detection",bif_studies)

df_all_Bif_longum <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_longum",bif_studies)

df_all_Bif_adolescentis <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_adolescentis",bif_studies)

df_all_Bif_pseudocatenulatum <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_pseudocatenulatum",bif_studies)

df_all_Bif_dentium <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_dentium",bif_studies)

df_all_Bif_bifidum <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_bifidum",bif_studies)

df_all_Bif_breve <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_breve",bif_studies)

df_all_Bif_catenulatum <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_catenulatum",bif_studies)

df_all_Bif_animalis <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_animalis",bif_studies)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### extract association df's ######
names(df_all_Bif_detection)
# [1] "corr" "pval" "dir" "association"

df_association_detection_16s_diet       <- df_all_Bif_detection$association
df_association_longum_16s_diet          <- df_all_Bif_longum$association
df_association_adolescentis_16s_diet    <- df_all_Bif_adolescentis$association
df_association_pseudocatenulatum_16s_diet <- df_all_Bif_pseudocatenulatum$association
df_association_dentium_16s_diet         <- df_all_Bif_dentium$association
df_association_bifidum_16s_diet         <- df_all_Bif_bifidum$association
df_association_breve_16s_diet           <- df_all_Bif_breve$association
df_association_catenulatum_16s_diet     <- df_all_Bif_catenulatum$association
df_association_animalis_16s_diet        <- df_all_Bif_animalis$association


####### make one combine df for association score ###########
df_association_16s_diet <- data.frame(
  detection       = df_association_detection_16s_diet$score,
  longum          = df_association_longum_16s_diet$score,
  adolescentis    = df_association_adolescentis_16s_diet$score,
  pseudocatenulatum = df_association_pseudocatenulatum_16s_diet$score,
  dentium         = df_association_dentium_16s_diet$score,
  bifidum         = df_association_bifidum_16s_diet$score,
  breve           = df_association_breve_16s_diet$score,
  catenulatum     = df_association_catenulatum_16s_diet$score,
  animalis        = df_association_animalis_16s_diet$score
)

# Preserve non-bif species names as rownames
rownames(df_association_16s_diet) <- rownames(df_association_detection_16s_diet)

# Quick check
dim(df_association_16s_diet)
#View(df_association_16s_diet)

#-----------------------------------------------------------------------------------------------------------------
###### correlating diet_association with association score ########

### 1. load the new association and make a subset ###
#View(df_association_all_new_16s)

# Subset adult bif columns
df_association_new_adult_16s <- df_association_all_new_16s[, grep("^adult_", colnames(df_association_all_new_16s))]
#View(df_association_new_adult_16s)

# Subset senior bif columns
df_association_new_senior_16s <- df_association_all_new_16s[, grep("^senior_", colnames(df_association_all_new_16s))]
#View(df_association_new_senior_16s)

# Subset adult/senior bif columns
# Make sure both data frames have identical rownames (species)
stopifnot(all(rownames(df_association_new_adult_16s) == rownames(df_association_new_senior_16s)))

# Create mean_association dataframe
mean_association_16s <- data.frame(
  detection = rowMeans(cbind(df_association_new_adult_16s$adult_detection,
                             df_association_new_senior_16s$senior_detection), na.rm = TRUE),
  longum = rowMeans(cbind(df_association_new_adult_16s$adult_longum,
                          df_association_new_senior_16s$senior_longum), na.rm = TRUE),
  adolescentis = rowMeans(cbind(df_association_new_adult_16s$adult_adolescentis,
                                df_association_new_senior_16s$senior_adolescentis), na.rm = TRUE),
  pseudocatenulatum = rowMeans(cbind(df_association_new_adult_16s$adult_pseudocatenulatum,
                                     df_association_new_senior_16s$senior_pseudocatenulatum), na.rm = TRUE),
  dentium = rowMeans(cbind(df_association_new_adult_16s$adult_dentium,
                           df_association_new_senior_16s$senior_dentium), na.rm = TRUE),
  bifidum = rowMeans(cbind(df_association_new_adult_16s$adult_bifidum,
                           df_association_new_senior_16s$senior_bifidum), na.rm = TRUE),
  breve = rowMeans(cbind(df_association_new_adult_16s$adult_breve,
                         df_association_new_senior_16s$senior_breve), na.rm = TRUE),
  catenulatum = rowMeans(cbind(df_association_new_adult_16s$adult_catenulatum,
                               df_association_new_senior_16s$senior_catenulatum), na.rm = TRUE),
  animalis = rowMeans(cbind(df_association_new_adult_16s$adult_animalis,
                            df_association_new_senior_16s$senior_animalis), na.rm = TRUE)
)

rownames(mean_association_16s) <- rownames(df_association_new_adult_16s)

class(mean_association_16s)
#View(mean_association_16s)

#2. Correlation
# Function to format p-values in scientific E notation
fmt_p <- function(x) {
  format(x, scientific = TRUE, digits = 3) |> toupper()
}

# Bifidobacterium species list
bif_species <- c("detection", "longum", "adolescentis", "pseudocatenulatum", "bifidum", "dentium", "breve", "catenulatum", "animalis")

# Create empty list to store results
results <- list()

for (sp in bif_species) {
  
  # === Mean correlation ===
  test_mean <- cor.test(
    df_association_16s_diet[[sp]],
    mean_association_16s[[sp]],
    method = "spearman",
    exact = TRUE
  )
  
  # Store as one row in list
  results[[sp]] <- data.frame(
    Bif_species = sp,
    corr_mean   = round(test_mean$estimate, 3),
    pval_mean   = fmt_p(test_mean$p.value),
    row.names = NULL
  )
}

# Combine into a single data frame
diet_16s_corr_df <- do.call(rbind, results)

# View final table
diet_16s_corr_df

write.csv(diet_16s_corr_df, "/results/diet_16s_corr.csv")

######################################################################################################################################################################
### DIET WGS STUDIES ###

load("/data/GhoshT_DietData_spProfile.RData")
load("/data/KeohaneDM_DietData_spProfile.RData")
load("/data/HeitzBushartA_DietData_spProfile.RData")
load("/data/JefferyIB_DietData_spProfile.RData")

pc_GhoshT_2020 <- compute_pca_pcoa(GhoshT_2020_metadata_diet)
pc_KeohaneDM <- compute_pca_pcoa(KeohaneDM_metadata_diet)
pc_HeitzBushartA_2016 <- compute_pca_pcoa(HeitzBushartA_2016_metadata_diet)
pc_JefferyIB <- compute_pca_pcoa(JefferyIB_metadata_diet)

summary_table <- rbind(
  extract_variance_summary(pc_GhoshT_2020, "GhoshT_2020"),
  extract_variance_summary(pc_KeohaneDM,         "KeohaneDM_2020"),
  extract_variance_summary(pc_HeitzBushartA_2016,     "Heitz-BuschartA_2016"),
  extract_variance_summary(pc_JefferyIB,       "JefferyIB_2020")
)

df_diet_wgs_variance <- summary_table

df_diet_wgs_variance[ , -1] <- round(df_diet_wgs_variance[ , -1] * 100, 2)
#df_diet_wgs_variance

write.xlsx(df_diet_wgs_variance, "/results/diet_wgs_variance_explained_percentage.xlsx")

#-------------------------------------------------------------------------------------------------
##### Bind PC's + SpProfiles ######

#### 1. GhoshT_2020  #####
identical(rownames(GhoshT_2020_speciesProfile), rownames(pc_GhoshT_2020))
GhoshT_2020_combined <- cbind(GhoshT_2020_speciesProfile, pc_GhoshT_2020)
#View(GhoshT_2020_combined)

#### 2. KeohaneDM_2020  #####
identical(rownames(KeohaneDM_speciesProfile), rownames(pc_KeohaneDM))
KeohaneDM_combined <- cbind(KeohaneDM_speciesProfile, pc_KeohaneDM)
#View(KeohaneDM_combined)

#### 3. Heitz-BuschartA_2016  #####
identical(rownames(HeitzBushartA_2016_speciesProfile), rownames(pc_HeitzBushartA_2016))
HeitzBushartA_2016_combined <- cbind(HeitzBushartA_2016_speciesProfile, pc_HeitzBushartA_2016)
#View(HeitzBushartA_2016_combined)

#### 4. JefferyIB_2020  #####
identical(rownames(JefferyIB_speciesProfile), rownames(pc_JefferyIB))
JefferyIB_combined <- cbind(JefferyIB_speciesProfile, pc_JefferyIB)
#View(JefferyIB_combined)

##### Overall spdata creation ######
spdata <- bind_rows(GhoshT_2020_combined, KeohaneDM_combined, HeitzBushartA_2016_combined, JefferyIB_combined)
dim(spdata)
#View(spdata)

any(is.na(spdata[, !colnames(spdata) %in% "study_name"]))
sum(rowSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)
sum(colSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)

bif_studies <- unique(spdata$study_name)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
##### compute pcor and make association_df ######

df_all_Bif_detection <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_detection",bif_studies)

df_all_Bif_longum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_longum",bif_studies)

df_all_Bif_adolescentis <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_adolescentis",bif_studies)

df_all_Bif_pseudocatenulatum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_pseudocatenulatum",bif_studies)

df_all_Bif_dentium <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_dentium",bif_studies)

df_all_Bif_bifidum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_bifidum",bif_studies)

df_all_Bif_breve <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_breve",bif_studies)

df_all_Bif_catenulatum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_catenulatum",bif_studies)

df_all_Bif_animalis <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_animalis",bif_studies)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### extract association df's ######
names(df_all_Bif_detection)
# [1] "corr" "pval" "dir" "association"

df_association_detection_wgs_diet       <- df_all_Bif_detection$association
df_association_longum_wgs_diet          <- df_all_Bif_longum$association
df_association_adolescentis_wgs_diet    <- df_all_Bif_adolescentis$association
df_association_pseudocatenulatum_wgs_diet <- df_all_Bif_pseudocatenulatum$association
df_association_dentium_wgs_diet         <- df_all_Bif_dentium$association
df_association_bifidum_wgs_diet         <- df_all_Bif_bifidum$association
df_association_breve_wgs_diet           <- df_all_Bif_breve$association
df_association_catenulatum_wgs_diet     <- df_all_Bif_catenulatum$association
df_association_animalis_wgs_diet        <- df_all_Bif_animalis$association

####### make one combine df for association score ###########
df_association_wgs_diet <- data.frame(
  detection       = df_association_detection_wgs_diet$score,
  longum          = df_association_longum_wgs_diet$score,
  adolescentis    = df_association_adolescentis_wgs_diet$score,
  pseudocatenulatum = df_association_pseudocatenulatum_wgs_diet$score,
  dentium         = df_association_dentium_wgs_diet$score,
  bifidum         = df_association_bifidum_wgs_diet$score,
  breve           = df_association_breve_wgs_diet$score,
  catenulatum     = df_association_catenulatum_wgs_diet$score,
  animalis        = df_association_animalis_wgs_diet$score
)

# Preserve non-bif species names as rownames
rownames(df_association_wgs_diet) <- rownames(df_association_detection_wgs_diet)

# Quick check
dim(df_association_wgs_diet)
#View(df_association_wgs_diet)

#----------------------------------------------------------------------------------------------------------------------
###### correlating diet_association with association score ########

### 1. load the new association and make a subset ###
#View(df_association_all_new_wgs)

# Subset adult bif columns
df_association_new_adult_wgs <- df_association_all_new_wgs[, grep("^adult_", colnames(df_association_all_new_wgs))]
#View(df_association_new_adult_wgs)

# Subset senior bif columns
df_association_new_senior_wgs <- df_association_all_new_wgs[, grep("^senior_", colnames(df_association_all_new_wgs))]
#View(df_association_new_senior_wgs)

# Subset adult/senior bif columns
# Make sure both data frames have identical rownames (species)
stopifnot(all(rownames(df_association_new_adult_wgs) == rownames(df_association_new_senior_wgs)))

# Create mean_association dataframe
mean_association_wgs <- data.frame(
  detection = rowMeans(cbind(df_association_new_adult_wgs$adult_detection,
                             df_association_new_senior_wgs$senior_detection), na.rm = TRUE),
  longum = rowMeans(cbind(df_association_new_adult_wgs$adult_longum,
                          df_association_new_senior_wgs$senior_longum), na.rm = TRUE),
  adolescentis = rowMeans(cbind(df_association_new_adult_wgs$adult_adolescentis,
                                df_association_new_senior_wgs$senior_adolescentis), na.rm = TRUE),
  pseudocatenulatum = rowMeans(cbind(df_association_new_adult_wgs$adult_pseudocatenulatum,
                                     df_association_new_senior_wgs$senior_pseudocatenulatum), na.rm = TRUE),
  dentium = rowMeans(cbind(df_association_new_adult_wgs$adult_dentium,
                           df_association_new_senior_wgs$senior_dentium), na.rm = TRUE),
  bifidum = rowMeans(cbind(df_association_new_adult_wgs$adult_bifidum,
                           df_association_new_senior_wgs$senior_bifidum), na.rm = TRUE),
  breve = rowMeans(cbind(df_association_new_adult_wgs$adult_breve,
                         df_association_new_senior_wgs$senior_breve), na.rm = TRUE),
  catenulatum = rowMeans(cbind(df_association_new_adult_wgs$adult_catenulatum,
                               df_association_new_senior_wgs$senior_catenulatum), na.rm = TRUE),
  animalis = rowMeans(cbind(df_association_new_adult_wgs$adult_animalis,
                            df_association_new_senior_wgs$senior_animalis), na.rm = TRUE)
)

rownames(mean_association_wgs) <- rownames(df_association_new_adult_wgs)

class(mean_association_wgs)
#View(mean_association_wgs)

#2. Correlation
# Function to format p-values in scientific E notation
fmt_p <- function(x) {
  format(x, scientific = TRUE, digits = 3) |> toupper()
}

# Bifidobacterium species list
bif_species <- c("detection", "longum", "adolescentis", "pseudocatenulatum", "bifidum", "dentium", "breve", "catenulatum", "animalis")

# Create empty list to store results
results <- list()

for (sp in bif_species) {

  # === Mean correlation ===
  test_mean <- cor.test(
    df_association_wgs_diet[[sp]],
    mean_association_wgs[[sp]],
    method = "spearman",
    exact = TRUE
  )
  
  # Store as one row in list
  results[[sp]] <- data.frame(
    Bif_species = sp,
    corr_mean   = round(test_mean$estimate, 3),
    pval_mean   = fmt_p(test_mean$p.value),
    row.names = NULL
  )
}

# Combine into a single data frame
diet_wgs_corr_df <- do.call(rbind, results)

# View final table
diet_wgs_corr_df

write.csv(diet_wgs_corr_df, "/results/diet_wgs_corr.csv")

######################################################################################################################################################################
### MEDICATION 16S STUDIES ###

load("/data/NUAGE_MedData_spProfile.RData")

pc_NUAGE <- compute_pca_pcoa(NUAGE_metadata_med)

summary_table <- extract_variance_summary(pc_NUAGE, "NUAGE")

df_med_16s_variance <- summary_table

df_med_16s_variance[ , -1] <- round(df_med_16s_variance[ , -1] * 100, 2)
df_med_16s_variance

write.xlsx(df_med_16s_variance, "/results/med_16s_variance_explained_percentage.xlsx")

#--------------------------------------------------------------------------------------------------------------------
##### Bind PC's + SpProfiles ######
identical(rownames(NUAGE_speciesProfile), rownames(pc_NUAGE))
NUAGE_combined <- cbind(NUAGE_speciesProfile, pc_NUAGE)
#View(NUAGE_combined)

#-----------------------------------------------------------------
######### pcorr.test function ##############
compute_partial_corr_matrix <- function(df, non_bifs, bifs, covariates) {
  if (!requireNamespace("ppcor", quietly = TRUE)) {
    stop("Please install 'ppcor' package first: install.packages('ppcor')")
  }
  
  # Initialize empty matrices
  corr_mat <- matrix(NA, nrow = length(non_bifs), ncol = length(bifs),
                     dimnames = list(non_bifs, bifs))
  pval_mat <- corr_mat
  
  # Loop over all combinations
  for (x in non_bifs) {
    for (y in bifs) {
      data_sub <- df[, c(x, y, covariates)]
      data_sub <- data_sub[complete.cases(data_sub), ]
      
      if (nrow(data_sub) >= 5) {
        res <- ppcor::pcor.test(data_sub[[x]], data_sub[[y]],
                                data_sub[, covariates, drop = FALSE],
                                method = "spearman")
        
        corr_mat[x, y] <- res$estimate
        pval_mat[x, y] <- res$p.value
      }
    }
  }
  
  # Replace missing values
  corr_mat[is.na(corr_mat)] <- 0
  pval_mat[is.na(pval_mat)] <- 1
  
  # Convert to data frames and rename columns
  corr_df <- as.data.frame(corr_mat)
  names(corr_df) <- paste0(names(corr_df), "_corr")
  
  pval_df <- as.data.frame(pval_mat)
  names(pval_df) <- paste0(names(pval_df), "_pval")
  
  return(list(corr = corr_df, pval = pval_df))
}

#for 16s studies
non_bifs <- rownames(df_association_all_new_16s)

bifs <- c("Bifidobacterium_detection", "Bifidobacterium_longum", "Bifidobacterium_adolescentis", "Bifidobacterium_pseudocatenulatum", "Bifidobacterium_dentium", "Bifidobacterium_bifidum", "Bifidobacterium_breve", "Bifidobacterium_catenulatum", "Bifidobacterium_animalis")

covars <- c("PCA1", "PCA2", "PCA3", "PCoA1", "PCoA2", "PCoA3")

NUAGE_pcor <- compute_partial_corr_matrix(NUAGE_combined, non_bifs, bifs, covars)

#mannual correction
pcor.test(NUAGE_combined$Akkermansia_muciniphila,NUAGE_combined$Bifidobacterium_animalis, NUAGE_combined[, covars],  method = "spearman")

# Extract the correlation dataframe
#NUAGE
NUAGE_corr <- NUAGE_pcor[["corr"]]
dim(NUAGE_corr)
#View(NUAGE_corr)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### calculate the prevalence ######
#View(NUAGE_speciesProfile)

# Remove study_name column
species_only <- NUAGE_speciesProfile[, !colnames(NUAGE_speciesProfile) %in% "study_name"]

# Calculate prevalence (%)
#prevalence_percent <- colSums(species_only > 0, na.rm = TRUE) / nrow(species_only) * 100

# Vector of Bifidobacterium columns
bif_cols <- c(
  "Bifidobacterium_detection",
  "Bifidobacterium_longum",
  "Bifidobacterium_adolescentis",
  "Bifidobacterium_pseudocatenulatum",
  "Bifidobacterium_bifidum",
  "Bifidobacterium_dentium",
  "Bifidobacterium_breve",
  "Bifidobacterium_catenulatum",
  "Bifidobacterium_animalis"
)

# Total samples
total_samples <- nrow(NUAGE_speciesProfile)

# Loop over columns and calculate prevalence
prevalence_df <- data.frame(
  Bif_species = bif_cols,
  prevalence_percent = sapply(bif_cols, function(col) {
    sum(NUAGE_speciesProfile[[col]] > 0, na.rm = TRUE) / total_samples * 100
  })
)

#prevalence_df

#---------------------------------------------------------------------------------------------------------------
####### correlation between medication association scores + association scores ###################
#View(df_association_new_senior_16s)

identical(rownames(NUAGE_corr), rownames(df_association_new_senior_16s))

# Run Spearman correlations and store results
tests <- list(
  detection = cor.test(
    NUAGE_corr$Bifidobacterium_detection_corr,
    df_association_new_senior_16s$senior_detection,
    method = "spearman",
    exact = TRUE
  ),
  longum = cor.test(
    NUAGE_corr$Bifidobacterium_longum_corr,
    df_association_new_senior_16s$senior_longum,
    method = "spearman",
    exact = TRUE
  ),
  adolescentis = cor.test(
    NUAGE_corr$Bifidobacterium_adolescentis_corr,
    df_association_new_senior_16s$senior_adolescentis,
    method = "spearman",
    exact = TRUE
  )
)

med_16s_corr_df <- data.frame(
  Bif_species  = names(tests),
  corr_senior  = sapply(tests, function(x) unname(x$estimate)),
  pval_senior  = sapply(tests, function(x) x$p.value),
  row.names = NULL
)

med_16s_corr_df

write.csv(med_16s_corr_df, "/results/med_16s_corr.csv")

######################################################################################################################################################################
### MEDICATION WGS STUDIES ###

load("/data/GhoshT_MedData_spProfile.RData")
load("/data/MetaCardis_MedData_spProfile.RData")
load("/data/HeitzBushartA_MedData_spProfile.RData")

pc_GhoshT_2020 <- compute_pca_pcoa(GhoshT_2020_metadata_med)
pc_MetaCardis <- compute_pca_pcoa(MetaCardis_metadata_med)

compute_pca_pcoa_Heitz <- function(df_temp) {
  
  ## --- Clean Data ---
  df_temp <- as.data.frame(lapply(df_temp, as.numeric))
  df_temp <- na.omit(df_temp)
  row_ids <- rownames(df_temp)
  
  ## ============================
  ##         SAFE PCA
  ## ============================
  safe_pca <- function(df) {
    out <- tryCatch({
      PcaHubert(df, k = min(3, ncol(df)), scale = TRUE, mcd = TRUE)
    }, error = function(e) NA)
    
    # If fails → classical PCA
    if (is.na(out)) {
      message("HeitzBuschartA: Robust PCA failed — switching to classical PCA.")
      pr <- prcomp(df, center = TRUE, scale. = TRUE)
      eig <- pr$sdev^2
      
      return(list(
        scores = pr$x[, 1:3],
        eigenvalues = eig
      ))
    }
    
    return(list(
      scores = out@scores[, 1:3],
      eigenvalues = out@eigenvalues
    ))
  }
  
  # Run safe PCA
  pca_res <- safe_pca(df_temp)
  pca_eig <- pca_res$eigenvalues
  pca_var_each <- pca_eig / sum(pca_eig)
  pca_first3 <- sum(pca_var_each[1:3])
  
  PCA_df <- as.data.frame(pca_res$scores)
  colnames(PCA_df) <- paste0("PCA", 1:3)
  
  
  ## ============================
  ##            PCoA
  ## ============================
  dist_mat <- dist(df_temp)
  pcoa_res <- cmdscale(dist_mat, eig = TRUE, k = 3)
  
  PCoA_df <- as.data.frame(pcoa_res$points)
  colnames(PCoA_df) <- paste0("PCoA", 1:3)
  
  eig <- pcoa_res$eig[pcoa_res$eig > 0]
  pcoa_var_each <- eig / sum(eig)
  pcoa_first3 <- sum(pcoa_var_each[1:3])
  
  
  ## ============================
  ##         COMBINED SCORES
  ## ============================
  combined <- cbind(PCA_df, PCoA_df)
  
  ## ============================
  ##         STORE ATTRIBUTES
  ## ============================
  attr(combined, "variance_explained") <- list(
    PCA_each = pca_var_each,
    PCA_first3 = pca_first3,
    PCoA_each = pcoa_var_each,
    PCoA_first3 = pcoa_first3
  )
  
  return(combined)
}

pc_HeitzBushartA_2016 <- compute_pca_pcoa_Heitz(HeitzBushartA_2016_metadata_med)

summary_table <- rbind(
  extract_variance_summary(pc_GhoshT_2020, "GhoshT_2020"),
  extract_variance_summary(pc_MetaCardis,         "MetaCardis_2020_a"),
  extract_variance_summary(pc_HeitzBushartA_2016,     "Heitz-BuschartA_2016")
)

df_med_wgs_variance <- summary_table

df_med_wgs_variance[ , -1] <- round(df_med_wgs_variance[ , -1] * 100, 2)
df_med_wgs_variance

write.xlsx(df_med_wgs_variance, "/results/med_wgs_variance_explained_percentage.xlsx")

#---------------------------------------------------------------------------------------------------------
##### Bind PC's + SpProfiles ######

#### 1. GhoshT_2020  #####
identical(rownames(GhoshT_2020_speciesProfile), rownames(pc_GhoshT_2020))
GhoshT_2020_combined <- cbind(GhoshT_2020_speciesProfile, pc_GhoshT_2020)
#View(GhoshT_2020_combined)

#### 2. MetaCardis_2020_a  #####
identical(rownames(MetaCardis_speciesProfile), rownames(pc_MetaCardis))
MetaCardis_combined <- cbind(MetaCardis_speciesProfile, pc_MetaCardis)
#View(MetaCardis_combined)

#### 3. Heitz-BuschartA_2016  #####
identical(rownames(HeitzBushartA_2016_speciesProfile), rownames(pc_HeitzBushartA_2016))
HeitzBushartA_2016_combined <- cbind(HeitzBushartA_2016_speciesProfile, pc_HeitzBushartA_2016)
#View(HeitzBushartA_2016_combined)

##### Overall spdata creation ######
spdata <- bind_rows(GhoshT_2020_combined, MetaCardis_combined, HeitzBushartA_2016_combined)
dim(spdata)
#View(spdata)

any(is.na(spdata[, !colnames(spdata) %in% "study_name"]))
sum(rowSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)
sum(colSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)

bif_studies <- unique(spdata$study_name)

#---------------------------------------------------------------------------------------------------------------------------------
##### compute pcor and make association_df ######

df_all_Bif_detection <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_detection",bif_studies)

df_all_Bif_longum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_longum",bif_studies)

df_all_Bif_adolescentis <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_adolescentis",bif_studies)

df_all_Bif_pseudocatenulatum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_pseudocatenulatum",bif_studies)

df_all_Bif_dentium <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_dentium",bif_studies)

df_all_Bif_bifidum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_bifidum",bif_studies)

df_all_Bif_breve <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_breve",bif_studies)

df_all_Bif_catenulatum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_catenulatum",bif_studies)

df_all_Bif_animalis <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_animalis",bif_studies)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### extract association df's ######
names(df_all_Bif_detection)
# [1] "corr" "pval" "dir" "association"

df_association_detection_wgs_med       <- df_all_Bif_detection$association
df_association_longum_wgs_med          <- df_all_Bif_longum$association
df_association_adolescentis_wgs_med    <- df_all_Bif_adolescentis$association
df_association_pseudocatenulatum_wgs_med <- df_all_Bif_pseudocatenulatum$association
df_association_dentium_wgs_med         <- df_all_Bif_dentium$association
df_association_bifidum_wgs_med         <- df_all_Bif_bifidum$association
df_association_breve_wgs_med           <- df_all_Bif_breve$association
df_association_catenulatum_wgs_med     <- df_all_Bif_catenulatum$association
df_association_animalis_wgs_med        <- df_all_Bif_animalis$association

####### make one combine df for association score ###########
df_association_wgs_med <- data.frame(
  detection       = df_association_detection_wgs_med$score,
  longum          = df_association_longum_wgs_med$score,
  adolescentis    = df_association_adolescentis_wgs_med$score,
  pseudocatenulatum = df_association_pseudocatenulatum_wgs_med$score,
  dentium         = df_association_dentium_wgs_med$score,
  bifidum         = df_association_bifidum_wgs_med$score,
  breve           = df_association_breve_wgs_med$score,
  catenulatum     = df_association_catenulatum_wgs_med$score,
  animalis        = df_association_animalis_wgs_med$score
)

# Preserve non-bif species names as rownames
rownames(df_association_wgs_med) <- rownames(df_association_detection_wgs_med)

# Quick check
dim(df_association_wgs_med)
#View(df_association_wgs_med)

#------------------------------------------------------------------------------------------------------------------------------
###### correlating medication_association with association score ########

### 1. load the new association and make a subset ###
#View(df_association_new_adult_wgs)
#View(df_association_new_senior_wgs)
#View(mean_association_wgs)

#2. Correlation
# Function to format p-values in scientific E notation
fmt_p <- function(x) {
  format(x, scientific = TRUE, digits = 3) |> toupper()
}

# Bifidobacterium species list
bif_species <- c("detection", "longum", "adolescentis", "pseudocatenulatum", "bifidum", "dentium", "breve", "catenulatum", "animalis")

# Create empty list to store results
results <- list()

for (sp in bif_species) {
  
  # === Mean correlation ===
  test_mean <- cor.test(
    df_association_wgs_med[[sp]],
    mean_association_wgs[[sp]],
    method = "spearman",
    exact = TRUE
  )
  
  # Store as one row in list
  results[[sp]] <- data.frame(
    Bif_species = sp,
    corr_mean   = round(test_mean$estimate, 3),
    pval_mean   = fmt_p(test_mean$p.value),
    row.names = NULL
  )
}

# Combine into a single data frame
med_wgs_corr_df <- do.call(rbind, results)

# View final table
#med_wgs_corr_df

write.csv(med_wgs_corr_df, "/results/med_wgs_corr.csv")

######################################################################################################################################################################
### CLINICAL 16S STUDIES ###

load("/data/CuestaZuluaga_BloodData_spProfile.RData")
load("/data/BartonW_BloodData_spProfile.RData")
load("/data/Hernandez_BloodData_spProfile.RData")
load("/data/He_BloodData_spProfile.RData")

pc_CuestaZuluaga <- compute_pca_pcoa(CuestaZuluaga_metadata_blood)
pc_BartonW <- compute_pca_pcoa(BartonW_metadata_blood)
pc_Hernandez <- compute_pca_pcoa(Hernandez_metadata_blood)
pc_He <- compute_pca_pcoa(He_metadata_blood)

summary_table <- rbind(
  extract_variance_summary(pc_CuestaZuluaga, "CuestaZuluagaJ_2018"),
  extract_variance_summary(pc_BartonW,         "BartonW_2018"),
  extract_variance_summary(pc_Hernandez,     "Hernandez_2023"),
  extract_variance_summary(pc_He,       "He")
)

df_clinical_16s_variance <- summary_table

df_clinical_16s_variance[ , -1] <- round(df_clinical_16s_variance[ , -1] * 100, 2)
df_clinical_16s_variance

write.xlsx(df_clinical_16s_variance, "/results/clinical_16s_variance_explained_percentage.xlsx")

#------------------------------------------------------------------------------------------------------------------------------------------
##### Bind PC's + SpProfiles ######

#### 1. CuestaZuluagaJ_2018  #####
identical(rownames(CuestaZuluaga_speciesProfile), rownames(pc_CuestaZuluaga))
CuestaZuluagaJ_combined <- cbind(CuestaZuluaga_speciesProfile, pc_CuestaZuluaga)
#View(CuestaZuluagaJ_combined)

#### 2. BartonW_2018  #####
identical(rownames(BartonW_speciesProfile), rownames(pc_BartonW))
BartonW_combined <- cbind(BartonW_speciesProfile, pc_BartonW)
#View(BartonW_combined)

#### 3. Hernandez_2023  #####
identical(rownames(hernandez_spProfile), rownames(pc_Hernandez))
Hernandez_combined <- cbind(hernandez_spProfile, pc_Hernandez)
#View(Hernandez_combined)

#### 4. He  #####
identical(rownames(He_speciesProfile), rownames(pc_He))
He_combined <- cbind(He_speciesProfile, pc_He)
#View(He_combined)

##### Overall spdata creation ######
spdata <- bind_rows(CuestaZuluagaJ_combined, BartonW_combined, Hernandez_combined, He_combined)
dim(spdata)
#View(spdata)

any(is.na(spdata[, !colnames(spdata) %in% "study_name"]))

spdata[, !colnames(spdata) %in% "study_name"] <- 
  replace(spdata[, !colnames(spdata) %in% "study_name"], 
          is.na(spdata[, !colnames(spdata) %in% "study_name"]), 
          0)

sum(rowSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)
sum(colSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)

bif_studies <- unique(spdata$study_name)

#-------------------------------------------------------------------------------------------------------------------------------------------
##### compute pcor and make association_df ######

df_all_Bif_detection <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_detection",bif_studies)

df_all_Bif_longum <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_longum",bif_studies)

df_all_Bif_adolescentis <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_adolescentis",bif_studies)

df_all_Bif_pseudocatenulatum <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_pseudocatenulatum",bif_studies)

df_all_Bif_dentium <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_dentium",bif_studies)

df_all_Bif_bifidum <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_bifidum",bif_studies)

df_all_Bif_breve <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_breve",bif_studies)

df_all_Bif_catenulatum <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_catenulatum",bif_studies)

df_all_Bif_animalis <- compute_partial_corr_associations_16s(spdata,df_association_all_new_16s,bif_name = "Bifidobacterium_animalis",bif_studies)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### extract association df's ######
names(df_all_Bif_detection)
# [1] "corr" "pval" "dir" "association"

df_association_detection_16s_blood       <- df_all_Bif_detection$association
df_association_longum_16s_blood          <- df_all_Bif_longum$association
df_association_adolescentis_16s_blood    <- df_all_Bif_adolescentis$association
df_association_pseudocatenulatum_16s_blood <- df_all_Bif_pseudocatenulatum$association
df_association_dentium_16s_blood         <- df_all_Bif_dentium$association
df_association_bifidum_16s_blood         <- df_all_Bif_bifidum$association
df_association_breve_16s_blood           <- df_all_Bif_breve$association
df_association_catenulatum_16s_blood     <- df_all_Bif_catenulatum$association
df_association_animalis_16s_blood        <- df_all_Bif_animalis$association

####### make one combine df for association score ###########
df_association_16s_blood <- data.frame(
  detection       = df_association_detection_16s_blood$score,
  longum          = df_association_longum_16s_blood$score,
  adolescentis    = df_association_adolescentis_16s_blood$score,
  pseudocatenulatum = df_association_pseudocatenulatum_16s_blood$score,
  dentium         = df_association_dentium_16s_blood$score,
  bifidum         = df_association_bifidum_16s_blood$score,
  breve           = df_association_breve_16s_blood$score,
  catenulatum     = df_association_catenulatum_16s_blood$score,
  animalis        = df_association_animalis_16s_blood$score
)

# Preserve non-bif species names as rownames
rownames(df_association_16s_blood) <- rownames(df_association_detection_16s_blood)

# Quick check
dim(df_association_16s_blood)
#View(df_association_16s_blood)

#--------------------------------------------------------------------------------------------------------------------------
####### correlation between clinical association scores + association scores ###################
#View(df_association_new_adult_16s)
#View(df_association_new_senior_16s)
#View(mean_association_16s)

#2. Correlation
# Function to format p-values in scientific E notation
fmt_p <- function(x) {
  format(x, scientific = TRUE, digits = 3) |> toupper()
}

# Bifidobacterium species list
bif_species <- c("detection", "longum", "adolescentis", "pseudocatenulatum", "bifidum", "dentium", "breve", "catenulatum", "animalis")

# Create empty list to store results
results <- list()

for (sp in bif_species) {
  
  # === Mean correlation ===
  test_mean <- cor.test(
    df_association_16s_blood[[sp]],
    mean_association_16s[[sp]],
    method = "spearman",
    exact = TRUE
  )
  
  # Store as one row in list
  results[[sp]] <- data.frame(
    Bif_species = sp,
    corr_mean   = round(test_mean$estimate, 3),
    pval_mean   = fmt_p(test_mean$p.value),
    row.names = NULL
  )
}

# Combine into a single data frame
clinical_16s_corr_df <- do.call(rbind, results)

# View final table
clinical_16s_corr_df

write.csv(clinical_16s_corr_df, "/results/clinical_16s_corr.csv")

######################################################################################################################################################################
### CLINICAL WGS STUDIES ###

load("/data/GhoshT_2020_BloodData_spProfile.RData")
load("/data/MetaCardis_BloodData_spProfile.RData")
load("/data/Cronin_BloodData_spProfile.RData")
load("/data/KeohaneDM_BloodData_spProfile.RData")
load("/data/HeitzBushartA_2016_BloodData_spProfile.RData")

pc_GhoshT_2020 <- compute_pca_pcoa(GhoshT_2020_metadata_blood)
pc_MetaCardis <- compute_pca_pcoa(MetaCardis_metadata_blood)
pc_Cronin <- compute_pca_pcoa(Cronin_metadata_blood)
pc_KeohaneDM <- compute_pca_pcoa(KeohaneDM_metadata_blood)
pc_HeitzBushartA_2016 <- compute_pca_pcoa(HeitzBushartA_2016_metadata_blood)

summary_table <- rbind(
  extract_variance_summary(pc_GhoshT_2020, "GhoshT_2020"),
  extract_variance_summary(pc_MetaCardis,         "MetaCardis_2020_a"),
  extract_variance_summary(pc_Cronin,     "Cronin_et_al_2018"),
  extract_variance_summary(pc_KeohaneDM,       "KeohaneDM_2020"),
  extract_variance_summary(pc_HeitzBushartA_2016,       "Heitz-BuschartA_2016")
)

df_clinical_wgs_variance <- summary_table

df_clinical_wgs_variance[ , -1] <- round(df_clinical_wgs_variance[ , -1] * 100, 2)
df_clinical_wgs_variance

write.xlsx(df_clinical_wgs_variance, "/results/clinical_wgs_variance_explained_percentage.xlsx")

#---------------------------------------------------------------------------------------------------------------------------
##### Bind PC's + SpProfiles ######

#### 1. GhoshT_2020  #####
identical(rownames(GhoshT_2020_speciesProfile), rownames(pc_GhoshT_2020))
GhoshT_2020_combined <- cbind(GhoshT_2020_speciesProfile, pc_GhoshT_2020)
#View(GhoshT_2020_combined)

#### 2. MetaCardis_2020_a  #####
identical(rownames(MetaCardis_speciesProfile), rownames(pc_MetaCardis))
MetaCardis_combined <- cbind(MetaCardis_speciesProfile, pc_MetaCardis)
#View(MetaCardis_combined)

#### 3. Cronin_et_al_2018  #####
identical(rownames(Cronin_speciesProfile), rownames(pc_Cronin))
Cronin_combined <- cbind(Cronin_speciesProfile, pc_Cronin)
#View(Cronin_combined)

#### 4. KeohaneDM_2020  #####
identical(rownames(KeohaneDM_speciesProfile), rownames(pc_KeohaneDM))
KeohaneDM_combined <- cbind(KeohaneDM_speciesProfile, pc_KeohaneDM)
#View(KeohaneDM_combined)

#### 5. Heitz-BuschartA_2016  #####
identical(rownames(HeitzBushartA_2016_speciesProfile), rownames(pc_HeitzBushartA_2016))
HeitzBushartA_2016_combined <- cbind(HeitzBushartA_2016_speciesProfile, pc_HeitzBushartA_2016)
#View(HeitzBushartA_2016_combined)

##### Overall spdata creation ######
spdata <- bind_rows(GhoshT_2020_combined, MetaCardis_combined, Cronin_combined, KeohaneDM_combined, HeitzBushartA_2016_combined)
dim(spdata)
#View(spdata)

any(is.na(spdata[, !colnames(spdata) %in% "study_name"]))
sum(rowSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)
sum(colSums(spdata[, !colnames(spdata) %in% "study_name"]) == 0)

bif_studies <- unique(spdata$study_name)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### compute pcor and make association_df ######

df_all_Bif_detection <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_detection",bif_studies)

df_all_Bif_longum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_longum",bif_studies)

df_all_Bif_adolescentis <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_adolescentis",bif_studies)

df_all_Bif_pseudocatenulatum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_pseudocatenulatum",bif_studies)

df_all_Bif_dentium <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_dentium",bif_studies)

df_all_Bif_bifidum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_bifidum",bif_studies)

df_all_Bif_breve <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_breve",bif_studies)

df_all_Bif_catenulatum <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_catenulatum",bif_studies)

df_all_Bif_animalis <- compute_partial_corr_associations_wgs(spdata,df_association_all_new_wgs,bif_name = "Bifidobacterium_animalis",bif_studies)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### extract association df's ######
names(df_all_Bif_detection)
# [1] "corr" "pval" "dir" "association"

df_association_detection_wgs_blood       <- df_all_Bif_detection$association
df_association_longum_wgs_blood          <- df_all_Bif_longum$association
df_association_adolescentis_wgs_blood    <- df_all_Bif_adolescentis$association
df_association_pseudocatenulatum_wgs_blood <- df_all_Bif_pseudocatenulatum$association
df_association_dentium_wgs_blood         <- df_all_Bif_dentium$association
df_association_bifidum_wgs_blood         <- df_all_Bif_bifidum$association
df_association_breve_wgs_blood           <- df_all_Bif_breve$association
df_association_catenulatum_wgs_blood     <- df_all_Bif_catenulatum$association
df_association_animalis_wgs_blood        <- df_all_Bif_animalis$association

####### make one combine df for association score ###########
df_association_wgs_blood <- data.frame(
  detection       = df_association_detection_wgs_blood$score,
  longum          = df_association_longum_wgs_blood$score,
  adolescentis    = df_association_adolescentis_wgs_blood$score,
  pseudocatenulatum = df_association_pseudocatenulatum_wgs_blood$score,
  dentium         = df_association_dentium_wgs_blood$score,
  bifidum         = df_association_bifidum_wgs_blood$score,
  breve           = df_association_breve_wgs_blood$score,
  catenulatum     = df_association_catenulatum_wgs_blood$score,
  animalis        = df_association_animalis_wgs_blood$score
)

# Preserve non-bif species names as rownames
rownames(df_association_wgs_blood) <- rownames(df_association_detection_wgs_blood)

# Quick check
dim(df_association_wgs_blood)
#View(df_association_wgs_blood)

#------------------------------------------------------------------------------------------------------------------------------
###### correlating clinical_association with association score ########

### 1. load the new association and make a subset ###
#View(df_association_new_adult_wgs)
#View(df_association_new_senior_wgs)
#View(mean_association_wgs)

#2. Correlation
# Function to format p-values in scientific E notation
fmt_p <- function(x) {
  format(x, scientific = TRUE, digits = 3) |> toupper()
}

# Bifidobacterium species list
bif_species <- c("detection", "longum", "adolescentis", "pseudocatenulatum", "bifidum", "dentium", "breve", "catenulatum", "animalis")

# Create empty list to store results
results <- list()

for (sp in bif_species) {

  # === Mean correlation ===
  test_mean <- cor.test(
    df_association_wgs_blood[[sp]],
    mean_association_wgs[[sp]],
    method = "spearman",
    exact = TRUE
  )
  
  # Store as one row in list
  results[[sp]] <- data.frame(
    Bif_species = sp,
    corr_mean   = round(test_mean$estimate, 3),
    pval_mean   = fmt_p(test_mean$p.value),
    row.names = NULL
  )
}

# Combine into a single data frame
clinical_wgs_corr_df <- do.call(rbind, results)

# View final table
#clinical_wgs_corr_df

write.csv(clinical_wgs_corr_df, "/results/clinical_wgs_corr.csv")

######################################################################################################################################################################
######################################################################################################################################################################
#diet_clinical barplots

#View(diet_16s_corr_df)
diet_16s_corr_df$pval_mean <- as.numeric(diet_16s_corr_df$pval_mean)
diet_16s_corr_df$significance <- ifelse(diet_16s_corr_df$pval_mean <= 0.05, "*", "")

#View(diet_wgs_corr_df)
diet_wgs_corr_df$pval_mean <- as.numeric(diet_wgs_corr_df$pval_mean)
diet_wgs_corr_df$significance <- ifelse(diet_wgs_corr_df$pval_mean <= 0.05, "*", "")

#View(med_16s_corr_df)
med_16s_corr_df$significance <- ifelse(med_16s_corr_df$pval_senior <= 0.05, "*", "")

#View(med_wgs_corr_df)
med_wgs_corr_df$pval_mean <- as.numeric(med_wgs_corr_df$pval_mean)
med_wgs_corr_df$significance <- ifelse(med_wgs_corr_df$pval_mean <= 0.05, "*", "")

#View(clinical_16s_corr_df)
clinical_16s_corr_df$pval_mean <- as.numeric(clinical_16s_corr_df$pval_mean)
clinical_16s_corr_df$significance <- ifelse(clinical_16s_corr_df$pval_mean <= 0.05, "*", "")

#View(clinical_wgs_corr_df)
clinical_wgs_corr_df$pval_mean <- as.numeric(clinical_wgs_corr_df$pval_mean)
clinical_wgs_corr_df$significance <- ifelse(clinical_wgs_corr_df$pval_mean <= 0.05, "*", "")

#### barplots generation ####

bif_colors <- c(
  detection = "#8DA0CB",
  longum = "#66C2A5",
  adolescentis = "#FC8D62",
  pseudocatenulatum = "#A6D854",
  bifidum = "#FFD92F",
  dentium = "#E78AC3",
  breve = "#A6761D",
  catenulatum = "#999999",
  animalis = "#80B1D3"
)

#outdir <- "/results/correlation_barplots_vertical"
#dir.create(outdir, showWarnings = FALSE)

plot_bif_bar <- function(df, corr_col, pval_col, title, filename) {
  
  df <- df %>%
    mutate(
      Bif_species = factor(Bif_species, levels = Bif_species),
      significance = ifelse(.data[[pval_col]] <= 0.05, "*", "")
    )
  
  p <- ggplot(df, aes(x = Bif_species, y = .data[[corr_col]], fill = Bif_species)) +
    geom_bar(stat = "identity", width = 0.65) +
    geom_text(
      aes(label = significance, y = .data[[corr_col]] + 0.03),
      size = 5
    ) +
    scale_fill_manual(values = bif_colors, name = "Bifidobacterium") +
    #coord_flip() +
    labs(
      title = title,
      y = "Spearman's rho",
      x = ""
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "right",
      axis.line.y = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y = element_line(color = "black", linewidth = 0.4),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 12),
      plot.title = element_text(face = "bold")
    ) +
    ylim(0, max(df[[corr_col]]) + 0.12)
  
  ggsave(
    filename = filename,
    plot = p,
    width = 7,
    height = 4,
    units = "in"
  )
}

## Diet
plot_bif_bar(diet_16s_corr_df, "corr_mean", "pval_mean",
             "Diet (16S)", "/results/diet_16s_correlation.pdf")

plot_bif_bar(diet_wgs_corr_df, "corr_mean", "pval_mean",
             "Diet (WGS)", "/results/diet_wgs_correlation.pdf")

## Medication diet
plot_bif_bar(med_16s_corr_df, "corr_senior", "pval_senior",
             "Medication (16S)", "/results/med_16s_correlation.pdf")

plot_bif_bar(med_wgs_corr_df, "corr_mean", "pval_mean",
             "Medication (WGS)", "/results/med_wgs_correlation.pdf")

## Clinical
plot_bif_bar(clinical_16s_corr_df, "corr_mean", "pval_mean",
             "Clinical (16S)", "/results/clinical_16s_correlation.pdf")

plot_bif_bar(clinical_wgs_corr_df, "corr_mean", "pval_mean",
             "Clinical (WGS)", "/results/clinical_wgs_correlation.pdf")

# to see average correlation values 
all_corr_values <- c(
  diet_16s_corr_df$corr_mean,
  diet_wgs_corr_df$corr_mean,
  med_16s_corr_df$corr_senior,
  med_wgs_corr_df$corr_mean,
  clinical_16s_corr_df$corr_mean,
  clinical_wgs_corr_df$corr_mean
)

avg_corr <- mean(abs(all_corr_values), na.rm = TRUE)
avg_corr

bif_features <- c(
  "longum",
  "adolescentis",
  "breve",
  "dentium",
  "detection"
)

filtered_corr <- c(
  diet_16s_corr_df$corr_mean[diet_16s_corr_df$Bif_species %in% bif_features],
  diet_wgs_corr_df$corr_mean[diet_wgs_corr_df$Bif_species %in% bif_features],
  med_16s_corr_df$corr_senior[med_16s_corr_df$Bif_species %in% bif_features],
  med_wgs_corr_df$corr_mean[med_wgs_corr_df$Bif_species %in% bif_features],
  clinical_16s_corr_df$corr_mean[clinical_16s_corr_df$Bif_species %in% bif_features],
  clinical_wgs_corr_df$corr_mean[clinical_wgs_corr_df$Bif_species %in% bif_features]
)

mean(abs(filtered_corr), na.rm = TRUE)

#------------------------------------------------------------------------------
#variance explained 

variance_dfs <- list(
  df_diet_16s_variance,
  df_diet_wgs_variance,
  df_med_16s_variance,
  df_med_wgs_variance,
  df_clinical_16s_variance,
  df_clinical_wgs_variance
)

variance_dfs <- lapply(variance_dfs, function(x) {
  x$Max_VarianceExplained <- pmax(x$PCA_First3, x$PCoA_First3, na.rm = TRUE)
  x
})


# Reassign back (same order as above)
df_diet_16s_variance     <- variance_dfs[[1]]
df_diet_wgs_variance     <- variance_dfs[[2]]
df_med_16s_variance      <- variance_dfs[[3]]
df_med_wgs_variance      <- variance_dfs[[4]]
df_clinical_16s_variance <- variance_dfs[[5]]
df_clinical_wgs_variance <- variance_dfs[[6]]

#Study name
df_diet_16s_variance$Study <- c("CuestaZuluagaJ_2018","NUAGE","Hernandez_2023","BartonW_2018")
df_diet_wgs_variance$Study
df_med_16s_variance$Study
df_med_wgs_variance$Study
df_clinical_16s_variance$Study
df_clinical_wgs_variance$Study

study_colors <- c(
  CuestaZuluagaJ_2018   = "#8DA0CB",
  NUAGE                = "#66C2A5",
  Hernandez_2023       = "#FC8D62",
  BartonW_2018         = "#A6D854",
  GhoshT_2020          = "#E78AC3",
  KeohaneDM_2020       = "#FFD92F",
  HeitzBuschartA_2016  = "#E9AFAF",
  JefferyIB_2020       = "#999999",
  MetaCardis_2020_a    = "#80B1D3",
  Cronin_et_al_2018    = "#B3B3B3",
  He                   = "#BEBADA"
)

#outdir <- "/results/variance_explained_barplots"


plot_variance_bar <- function(df, title, filename) {
  
  df <- df %>%
    mutate(
      Study = factor(Study, levels = Study)
    )
  
  p <- ggplot(df, aes(x = Study, y = Max_VarianceExplained, fill = Study)) +
    geom_bar(stat = "identity", width = 0.65) +
    scale_fill_manual(values = study_colors, name = "Study") +
    labs(
      title = title,
      y = "Max Variance Explained (%)",
      x = ""
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line.y = element_line(color = "black", linewidth = 0.4),
      axis.ticks.y = element_line(color = "black", linewidth = 0.4),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    ) +
    ylim(0, 110)
  
  ggsave(
    filename = filename,
    plot = p,
    width = 7,
    height = 4,
    units = "in"
  )
}

## Diet
plot_variance_bar(df_diet_16s_variance, "Diet (16S)", "/results/diet_16s_variance.pdf")
plot_variance_bar(df_diet_wgs_variance, "Diet (WGS)", "/results/diet_wgs_variance.pdf")

## Medication
plot_variance_bar(df_med_16s_variance, "Medication (16S)", "/results/med_16s_variance.pdf")
plot_variance_bar(df_med_wgs_variance, "Medication (WGS)", "/results/med_wgs_variance.pdf")

## Clinical
plot_variance_bar(df_clinical_16s_variance, "Clinical (16S)", "/results/clinical_16s_variance.pdf")
plot_variance_bar(df_clinical_wgs_variance, "Clinical (WGS)", "/results/clinical_wgs_variance.pdf")


#
# Combine all Max_VarianceExplained values
all_max_variance <- c(
  df_diet_16s_variance$Max_VarianceExplained,
  df_diet_wgs_variance$Max_VarianceExplained,
  df_med_16s_variance$Max_VarianceExplained,
  df_med_wgs_variance$Max_VarianceExplained,
  df_clinical_16s_variance$Max_VarianceExplained,
  df_clinical_wgs_variance$Max_VarianceExplained
)

# Average
avg_max_variance_explained <- mean(all_max_variance, na.rm = TRUE)

avg_max_variance_explained



######################################################################################################################################################################
######################################################################################################################################################################
### Graded Adjusted Association Scores 

load("/data/GroupOrders.RData")

df_partial_corr_16s_med <- NUAGE_corr
df_association_16s_clinical <- df_association_16s_blood
df_association_wgs_clinical <- df_association_wgs_blood

####  Part 1  ####
#maintain the bif order in all the adjusted assocaition scores df
bif_order <- c("detection","longum","adolescentis","pseudocatenulatum","bifidum","dentium","breve","catenulatum","animalis")

clean_bif_names <- function(x) {
  x <- gsub("^Bifidobacterium_", "", x)
  x <- gsub("_corr$", "", x)
  x
}
colnames(df_partial_corr_16s_med) <- clean_bif_names(colnames(df_partial_corr_16s_med))

reorder_bif_cols <- function(df, order_vec) {
  df[, order_vec, drop = FALSE]
}

df_association_16s_diet       <- reorder_bif_cols(df_association_16s_diet, bif_order)
df_association_16s_clinical  <- reorder_bif_cols(df_association_16s_clinical, bif_order)
df_partial_corr_16s_med      <- reorder_bif_cols(df_partial_corr_16s_med, bif_order)

df_association_wgs_diet      <- reorder_bif_cols(df_association_wgs_diet, bif_order)
df_association_wgs_clinical <- reorder_bif_cols(df_association_wgs_clinical, bif_order)
df_association_wgs_med      <- reorder_bif_cols(df_association_wgs_med, bif_order)

lapply(
  list(
    df_association_16s_diet,
    df_association_16s_clinical,
    df_partial_corr_16s_med,
    df_association_wgs_diet,
    df_association_wgs_clinical,
    df_association_wgs_med
  ),
  colnames)


####  Part 2 a ####
#concatenate all the dataframe into one for wgs
colnames(df_association_wgs_diet) <- paste0("diet_", colnames(df_association_wgs_diet))

colnames(df_association_wgs_med) <- paste0("medication_", colnames(df_association_wgs_med))

colnames(df_association_wgs_clinical) <- paste0("clinical_", colnames(df_association_wgs_clinical))

AdjustedAssociationScores_wgs_combined <- cbind(df_association_wgs_diet,df_association_wgs_med,df_association_wgs_clinical)
identical(rownames(AdjustedAssociationScores_wgs_combined), rownames(df_association_wgs_clinical))

####  Part 2 b ####
#concatenate all the dataframe into one for 16s
colnames(df_association_16s_diet) <- paste0("diet_", colnames(df_association_16s_diet))

colnames(df_partial_corr_16s_med) <- paste0("medication_", colnames(df_partial_corr_16s_med))

colnames(df_association_16s_clinical) <- paste0("clinical_", colnames(df_association_16s_clinical))

AdjustedAssociationScores_16s_combined <- cbind(df_association_16s_diet, df_partial_corr_16s_med, df_association_16s_clinical)
identical(rownames(AdjustedAssociationScores_16s_combined), rownames(df_association_16s_clinical))


####  Part 3 ####
#reorder non-bif rownames based on Group 
AdjustedAssociationScores_wgs_combined <- AdjustedAssociationScores_wgs_combined[rownames(carpet_StrongAssociationScore_Adjusted_wgs), ,drop = FALSE]
identical(rownames(AdjustedAssociationScores_wgs_combined),rownames(carpet_StrongAssociationScore_Adjusted_wgs))

AdjustedAssociationScores_16s_combined <- AdjustedAssociationScores_16s_combined[rownames(carpet_StrongAssociationScore_Adjusted_16s), ,drop = FALSE]
identical(rownames(AdjustedAssociationScores_16s_combined),rownames(carpet_StrongAssociationScore_Adjusted_16s))

save(AdjustedAssociationScores_16s_combined, AdjustedAssociationScores_wgs_combined, file = "/results/AdjustedAssociationScores_combined_DietClinical.RData")


####  Part 4 ####
####### Graded Association Scores #########

### For WGS ###

subject_groups <- c(
  "diet",
  "medication",
  "clinical")

thresholds_for_association_scores <- as.data.frame(
  matrix(NA_real_, nrow = length(subject_groups), ncol = 2,
         dimnames = list(subject_groups, c("Lower", "Higher")))
)


for (grp in subject_groups) {
  
  cols <- grep(paste0("^", grp, "_"), colnames(AdjustedAssociationScores_wgs_combined), value = TRUE)
  
  vals <- unlist(AdjustedAssociationScores_wgs_combined[, cols, drop = FALSE])
  
  qs <- quantile(vals, probs = c(0.1, 0.9), na.rm = TRUE)
  
  thresholds_for_association_scores[grp, "Lower"]  <- qs[1]
  thresholds_for_association_scores[grp, "Higher"] <- qs[2]
}


grade_associations <- function(x, group) {
  ifelse(
    x > thresholds_for_association_scores[group, "Higher"],  2,
    ifelse(
      x < thresholds_for_association_scores[group, "Lower"], -2,
      sign(x)
    )
  )
}

dig_AdjustedAssociationScores_wgs_combined <- do.call(
  cbind,
  lapply(subject_groups, function(grp) {
    
    cols <- grep(paste0("^", grp, "_"), colnames(AdjustedAssociationScores_wgs_combined), value = TRUE)
    
    apply(
      AdjustedAssociationScores_wgs_combined[, cols, drop = FALSE],
      2,
      grade_associations,
      group = grp
    )
  })
)

dig_AdjustedAssociationScores_wgs_combined <- as.data.frame(dig_AdjustedAssociationScores_wgs_combined)

df_overall_pattern_wgs <- data.frame(
  positive = apply(dig_AdjustedAssociationScores_wgs_combined, 1, function(x) sum(x == 2, na.rm = TRUE)),
  negative = apply(dig_AdjustedAssociationScores_wgs_combined, 1, function(x) sum(x == -2, na.rm = TRUE)))




### For 16S ###

subject_groups <- c(
  "diet",
  "medication",
  "clinical")

thresholds_for_association_scores_16s <- as.data.frame(
  matrix(NA_real_, nrow = length(subject_groups), ncol = 2,
         dimnames = list(subject_groups, c("Lower", "Higher")))
)


for (grp in subject_groups) {
  
  cols <- grep(paste0("^", grp, "_"), colnames(AdjustedAssociationScores_16s_combined), value = TRUE)
  
  vals <- unlist(AdjustedAssociationScores_16s_combined[, cols, drop = FALSE])
  
  qs <- quantile(vals, probs = c(0.1, 0.9), na.rm = TRUE)
  
  thresholds_for_association_scores_16s[grp, "Lower"]  <- qs[1]
  thresholds_for_association_scores_16s[grp, "Higher"] <- qs[2]
}


grade_associations <- function(x, group) {
  ifelse(
    x > thresholds_for_association_scores_16s[group, "Higher"],  2,
    ifelse(
      x < thresholds_for_association_scores_16s[group, "Lower"], -2,
      sign(x)
    )
  )
}

dig_AdjustedAssociationScores_16s_combined <- do.call(
  cbind,
  lapply(subject_groups, function(grp) {
    
    cols <- grep(paste0("^", grp, "_"), colnames(AdjustedAssociationScores_16s_combined), value = TRUE)
    
    apply(
      AdjustedAssociationScores_16s_combined[, cols, drop = FALSE],
      2,
      grade_associations,
      group = grp
    )
  })
)

dig_AdjustedAssociationScores_16s_combined <- as.data.frame(dig_AdjustedAssociationScores_16s_combined)

df_overall_pattern_16s <- data.frame(
  positive = apply(dig_AdjustedAssociationScores_16s_combined, 1, function(x) sum(x == 2, na.rm = TRUE)),
  negative = apply(dig_AdjustedAssociationScores_16s_combined, 1, function(x) sum(x == -2, na.rm = TRUE)))


############ Re-ordering Columns of WGS #########

species_order <- c(
  "detection",
  "longum",
  "adolescentis",
  "pseudocatenulatum",
  "bifidum",
  "dentium",
  "breve",
  "catenulatum",
  "animalis")

group_order <- c(
  "diet",
  "medication",
  "clinical")

cn <- colnames(dig_AdjustedAssociationScores_wgs_combined)

## split column names into group + species
parts <- do.call(rbind, strsplit(cn, "_", fixed = TRUE))
col_df <- data.frame(
  group   = parts[, 1],
  species = parts[, 2],
  colname = cn,
  stringsAsFactors = FALSE
)

## convert to ordered factors
col_df$species <- factor(col_df$species, levels = species_order)
col_df$group   <- factor(col_df$group, levels = group_order)

## order: species first, then group
col_df <- col_df[order(col_df$species, col_df$group), ]

## reorder dataframe
dig_AdjustedAssociationScores_wgs_combined <-
  dig_AdjustedAssociationScores_wgs_combined[, col_df$colname]


############ Re-ordering Columns of 16S #########

species_order <- c(
  "detection",
  "longum",
  "adolescentis",
  "pseudocatenulatum",
  "bifidum",
  "dentium",
  "breve",
  "catenulatum",
  "animalis")

group_order <- c(
  "diet",
  "medication",
  "clinical")

cn <- colnames(dig_AdjustedAssociationScores_16s_combined)

## split column names into group + species
parts <- do.call(rbind, strsplit(cn, "_", fixed = TRUE))
col_df <- data.frame(
  group   = parts[, 1],
  species = parts[, 2],
  colname = cn,
  stringsAsFactors = FALSE
)

## convert to ordered factors
col_df$species <- factor(col_df$species, levels = species_order)
col_df$group   <- factor(col_df$group, levels = group_order)

## order: species first, then group
col_df <- col_df[order(col_df$species, col_df$group), ]

## reorder dataframe
dig_AdjustedAssociationScores_16s_combined <-
  dig_AdjustedAssociationScores_16s_combined[, col_df$colname]

########## Heat map of the Graded Associations in WGS ######

mat <- as.matrix(dig_AdjustedAssociationScores_wgs_combined)

col_fun <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("deeppink4","pink","gray95","skyblue","royalblue4")
)


ht <- Heatmap(
  mat,
  name = "Association Scores",
  col = col_fun,
  cluster_rows = FALSE,                   
  cluster_columns = FALSE,               
  show_row_dend = FALSE,                 
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7.7),
  column_names_gp = gpar(fontsize = 7.7),
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title = "Association Scores"
  ),
  width = unit(ncol(mat) * 3, "mm"),
  height = unit(nrow(mat) * 3, "mm"),
  
  # Use cell_fun instead of layer_fun
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width = width, height = height, 
              gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
  }
)

pdf("/results/dig_AdjustedAssociationScores_wgs_combined_new.pdf", 
    height = 20, width = 16)
draw(ht)
dev.off()

##### Carpet #####

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_dig_AdjustedAssociationScores_wgs_combined <- as.data.frame(carpet_df)

#write.xlsx(carpet_dig_AdjustedAssociationScores_wgs_combined, "carpet_dig_AdjustedAssociationScores_wgs_combined_NEW.xlsx", rowNames = TRUE)

#####  Line Plot #####

association_long <- df_overall_pattern_wgs %>%
  mutate(Species = rownames(.)) %>%
  pivot_longer(cols = c("positive", "negative"),
               names_to = "Association",
               values_to = "Count")

# Retain rowname order
association_long$Species <- factor(
  association_long$Species,
  levels = rownames(df_overall_pattern_wgs)
)

pdf("/results/lineplot_AdjustedAssociationScores_wgs_NEW.pdf",
    height = 1, width = 10)

ggplot(association_long, aes(x = Species, y = Count, group = Association, color = Association)) +
  geom_line(linewidth = 0.3) +
  geom_point(shape = 16, size = 0.7) +
  scale_color_manual(values = c("positive" = "blue", "negative" = "red")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "Species", y = "Count", color = "Association Type")

dev.off()


########## Heat map of the Graded Associations in 16S ######


mat <- as.matrix(dig_AdjustedAssociationScores_16s_combined)

col_fun <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("deeppink4","pink","gray95","skyblue","royalblue4")
)


ht <- Heatmap(
  mat,
  name = "Association Scores",
  col = col_fun,
  cluster_rows = FALSE,                   
  cluster_columns = FALSE,               
  show_row_dend = FALSE,                 
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7.7),
  column_names_gp = gpar(fontsize = 7.7),
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title = "Association Scores"
  ),
  width = unit(ncol(mat) * 3, "mm"),
  height = unit(nrow(mat) * 3, "mm"),
  
  # Use cell_fun instead of layer_fun
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width = width, height = height, 
              gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
  }
)

pdf("/results/dig_AdjustedAssociationScores_16s_combined_new.pdf", 
    height = 25, width = 16)
draw(ht)
dev.off()

##### Carpet #####

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_dig_AdjustedAssociationScores_16s_combined <- as.data.frame(carpet_df)

#write.xlsx(carpet_dig_AdjustedAssociationScores_16s_combined, "carpet_dig_AdjustedAssociationScores_16s_combined_NEW.xlsx", rowNames = TRUE)

#####  Line Plot #####

association_long <- df_overall_pattern_16s %>%
  mutate(Species = rownames(.)) %>%
  pivot_longer(cols = c("positive", "negative"),
               names_to = "Association",
               values_to = "Count")

# Retain rowname order
association_long$Species <- factor(
  association_long$Species,
  levels = rownames(df_overall_pattern_16s)
)

pdf("/results/lineplot_AdjustedAssociationScores_16s_new.pdf",
    height = 1, width = 10)

ggplot(association_long, aes(x = Species, y = Count, group = Association, color = Association)) +
  geom_line(linewidth = 0.3) +
  geom_point(shape = 16, size = 0.7) +
  scale_color_manual(values = c("positive" = "blue", "negative" = "red")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "Species", y = "Count", color = "Association Type")

dev.off()

######################################################################################################################################################################
######################################################################################################################################################################



























































































