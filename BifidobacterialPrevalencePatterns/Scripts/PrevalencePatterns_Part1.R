
load("/data/bifido_45809_new_2014_not_normalized_Harnandez_added.RData")

library(dplyr)
library(pcaPP)
library(psych)
library(vegan)

#shannon45809 <- diversity(SpeciesProfile_All, index = "shannon")

# calculate_pielou <- function(abundance) {
#   H <- diversity(abundance, index = "shannon")
#   S <- sum(abundance > 0)
#   if (S > 1) {
#     E <- H / log(S)
#   } else {
#     E <- NA
#   }
#   return(E)
# }

SpeciesProfile_All[is.na(SpeciesProfile_All)] <- 0

names(SpeciesProfile_All)[1265] <- "Intestinibacter_bartlettii"
names(SpeciesProfile_All)[1362] <- "Erysipelatoclostridium_ramosum"

spdata_45809 <- SpeciesProfile_All
metadata_45809 <- metadata_updated

AllBifidoSpecies <- grep("Bifidobacterium",colnames(spdata_45809),value=TRUE)

spdata_Bifs_rel_45809 <- apply(spdata_45809[,AllBifidoSpecies],2,function(x)(ifelse(is.na(x),0,x)))

spdata_Bifs_detect_45809 <- apply(spdata_Bifs_rel_45809,2,function(x)(ifelse(x>0.0001,1,0)))

studies_45809 <- unique(metadata_45809$study_name)

spdata_45809 <- spdata_45809[intersect(rownames(spdata_45809),rownames(metadata_45809)),]

metadata_45809 <- metadata_45809[intersect(rownames(spdata_45809),rownames(metadata_45809)),]

infant <- rownames(metadata_45809[metadata_45809$age_category == "infant",])

adult <- rownames(metadata_45809[metadata_45809$age_category == "adult",])

senior <- rownames(metadata_45809[metadata_45809$age_category == "senior",])

AllSpecies <- colnames(spdata_45809)

spdata_45809$study_name <- metadata_45809$study_name

infant_studies <- unique(metadata_45809[infant,"study_name"])

adult_studies <- unique(metadata_45809[adult,"study_name"])

senior_studies <- unique(metadata_45809[senior,"study_name"])



compute_detection <- function(data,var1_list,grouping_variable,grouping_list)
{
  detection_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
  rownames(detection_matrix) <- var1_list
  colnames(detection_matrix) <- grouping_list
  for(i in 1:length(var1_list))
  {
    var1 <- var1_list[i]
    for(j in 1:length(grouping_list))
    {
      group <- grouping_list[j]
      detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]>0))/length(data[data[,grouping_variable]==group,var1])
    }
  }
  return(detection_matrix)
}

### To calculate the prevalence please run the following codes ###

# print("Infants Overall Species Detection")
# infant_species_detection <- compute_detection(spdata_45809[infant,],AllSpecies,"study_name",infant_studies)
# 
# print("Adult Overall Species Detection")
# adult_species_detection <- compute_detection(spdata_45809[adult,],AllSpecies,"study_name",adult_studies)
# 
# print("Senior Overall Species Detection")
# senior_species_detection <- compute_detection(spdata_45809[senior,],AllSpecies,"study_name",senior_studies)

load("G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\MANUSCRIPT_ITEMS\\REVISED_FINAL\\CodeOcean_GitHub\\BifidobacterialPrevalencePatterns\\SpeciesDetection_AgeCategories.RData")

infant_select_species <- names(which(apply(infant_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(infant_species_detection)>=0.33))
adult_select_species <- names(which(apply(adult_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(adult_species_detection)>=0.33))
senior_select_species <- names(which(apply(senior_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(senior_species_detection)>=0.33))


infant_select_non_bif_species <- setdiff(infant_select_species,AllBifidoSpecies)
adult_select_non_bif_species <- setdiff(adult_select_species,AllBifidoSpecies)
senior_select_non_bif_species <- setdiff(senior_select_species,AllBifidoSpecies)

infant_studies_sample_numbers <- table(metadata_45809[infant,"study_name"])
infant_select_studies <- names(which(infant_studies_sample_numbers>=20))

adult_studies_sample_numbers <- table(metadata_45809[adult,"study_name"])
adult_select_studies <- names(which(adult_studies_sample_numbers>=20))

senior_studies_sample_numbers <- table(metadata_45809[senior,"study_name"])
senior_select_studies <- names(which(senior_studies_sample_numbers>=20))

MajorBifidoSpecies <- c("Bifidobacterium_adolescentis","Bifidobacterium_animalis","Bifidobacterium_bifidum","Bifidobacterium_breve","Bifidobacterium_catenulatum","Bifidobacterium_dentium","Bifidobacterium_longum","Bifidobacterium_pseudocatenulatum")



###################################################################
                      ### Infant   ####

df_infant_bifido_new <- as.data.frame(matrix(NA, length(infant_select_studies), 30))
rownames(df_infant_bifido_new) <- infant_select_studies
colnames(df_infant_bifido_new) <- c(
  MajorBifidoSpecies,
  "BifPrevalence",
  "BrayAssocAnimalis",
  "KendallAssocAnimalis",
  "BrayAssocCatenulatum",
  "KendallAssocCatenulatum",
  "BrayAssocBreve",
  "KendallAssocBreve",
  "BrayAssocDentium",
  "KendallAssocDentium",
  "BrayAssocBifidum",
  "KendallAssocBifidum",
  "BrayAssocPseudocatenulatum",
  "KendallAssocPseudocatenulatum",
  "BrayAssocAdolescentis",
  "KendallAssocAdolescentis",
  "BrayAssocLongum",
  "KendallAssocLongum",
  "BrayAssocBifDetection",
  "KendallAssocBifDetection",
  "SequencingType",
  "CohortType",
  "AgeCategory")

for (i in seq_along(infant_select_studies)) {
  study_name <- infant_select_studies[i]
  print(study_name)
  
  ## ---- Get samples for this study ----
  study_rows <- intersect(
    infant,
    rownames(metadata_45809)[metadata_45809$study_name == study_name]
  )
  
  if (length(study_rows) == 0) next
  
  ## ---- Prevalence of major Bifidobacterium species ----
  df_infant_bifido_new[study_name, 1:8] <- apply(
    spdata_45809[study_rows, MajorBifidoSpecies, drop = FALSE],
    2,
    function(x) length(x[(x > 0.0001) & !is.na(x)])
  ) / length(study_rows)
  
  ## ---- Overall Bifidobacterium detection prevalence ----
  df_infant_bifido_new[study_name, "BifPrevalence"] <- length(which(
    apply(
      spdata_45809[study_rows, AllBifidoSpecies, drop = FALSE],
      1,
      function(x) length(x[(x > 0.0001) & !is.na(x)])
    ) > 0
  )) / length(study_rows)
  
  ## ---- Subset to non-bifidobacterium species ----
  temp_spdata <- spdata_45809[study_rows, infant_select_non_bif_species, drop = FALSE]
  temp_spdata <- apply(temp_spdata, 2, function(x) ifelse(is.na(x), 0, x))
  temp_spdata <- temp_spdata[rowSums(temp_spdata) > 0, , drop = FALSE]
  temp_spdata <- temp_spdata[, colSums(temp_spdata) > 0, drop = FALSE]
  
  if (nrow(temp_spdata) < 2 || ncol(temp_spdata) < 2) {
    warning(paste("Not enough data for PERMANOVA in study:", study_name))
    next
  }
  
  temp_spdata <- temp_spdata / rowSums(temp_spdata)
  
  ## ---- Bifidobacterium detection ----
  BifidoDetection <- apply(
    spdata_45809[rownames(temp_spdata), AllBifidoSpecies, drop = FALSE],
    1,
    function(x) length(x[x >= 0.0001])
  )
  
  ## ---- Metadata for adonis2: ONLY BifidoDetection ----
  temp_meta <- data.frame(
    BifidoDetection = BifidoDetection
  )
  
  ## ---- PERMANOVA: no covariates, only BifidoDetection ----
  print("Bray influence")
  tempAdonisBray <- adonis2(
    vegdist(temp_spdata, method = "bray") ~ BifidoDetection,
    data = temp_meta
  )
  df_infant_bifido_new[study_name, "BrayAssocBifDetection"] <- tempAdonisBray$Pr[1]
  
  print("Kendall influence")
  tempAdonisKendall <- adonis2(
    as.dist(1 - cor.fk(t(temp_spdata)) / 2) ~ BifidoDetection,
    data = temp_meta
  )
  df_infant_bifido_new[study_name, "KendallAssocBifDetection"] <- tempAdonisKendall$Pr[1]
  
  ## ---- Fill SequencingType, CohortType, AgeCategory from metadata_45809 ----
  # Sequencing type (seq_type)
  seq_vals <- unique(metadata_45809[study_rows, "seq_type"])
  seq_vals <- seq_vals[!is.na(seq_vals)]
  if (length(seq_vals) > 1) {
    warning(paste("Multiple seq_type values for study:", study_name, "– using first."))
  }
  if (length(seq_vals) >= 1) {
    df_infant_bifido_new[study_name, "SequencingType"] <- seq_vals[1]
  }
  
  # Cohort type (cohort)
  cohort_vals <- unique(metadata_45809[study_rows, "cohort"])
  cohort_vals <- cohort_vals[!is.na(cohort_vals)]
  if (length(cohort_vals) > 1) {
    warning(paste("Multiple cohort values for study:", study_name, "– using first."))
  }
  if (length(cohort_vals) >= 1) {
    df_infant_bifido_new[study_name, "CohortType"] <- cohort_vals[1]
  }
  
  # Age category (age_category)
  age_vals <- unique(metadata_45809[study_rows, "age_category"])
  age_vals <- age_vals[!is.na(age_vals)]
  if (length(age_vals) > 1) {
    warning(paste("Multiple age_category values for study:", study_name, "– using first."))
  }
  if (length(age_vals) >= 1) {
    df_infant_bifido_new[study_name, "AgeCategory"] <- age_vals[1]
  }
}

####### For Eight different Bifidobacterium – example for B. pseudocatenulatum #######

# for (i in seq_along(infant_select_studies)) {
#   study_name <- infant_select_studies[i]
#   print(study_name)
#   
#   # ---- Get samples for this study ----
#   study_rows <- intersect(
#     infant,
#     rownames(metadata_45809)[metadata_45809$study_name == study_name]
#   )
#   
#   if (length(study_rows) == 0) {
#     warning(study_name, ": no samples found - skipping")
#     next
#   }
#   
#   # ---- Subset to non-bifidobacterium species ----
#   temp_spdata <- spdata_45809[study_rows, infant_select_non_bif_species, drop = FALSE]
#   temp_spdata <- apply(temp_spdata, 2, function(x) ifelse(is.na(x), 0, x))
#   temp_spdata <- temp_spdata[rowSums(temp_spdata) > 0, , drop = FALSE]
#   temp_spdata <- temp_spdata[, colSums(temp_spdata) > 0, drop = FALSE]
#   
#   # skip if too few samples or species after filtering
#   if (nrow(temp_spdata) < 3 || ncol(temp_spdata) < 2) {
#     warning(study_name, ": insufficient data after filtering - skipping")
#     next
#   }
#   
#   temp_spdata <- temp_spdata / rowSums(temp_spdata)
#   
#   # ---- Bifidobacterium pseudocatenulatum abundance ----
#   temp_pseudocatenulatum <- spdata_45809[rownames(temp_spdata), "Bifidobacterium_pseudocatenulatum"]
#   
#   # ---- Metadata for adonis2: ONLY temp_longum ----
#   temp_meta <- data.frame(
#     temp_pseudocatenulatum = temp_pseudocatenulatum
#   )
#   
#   # ---- Bray-Curtis PERMANOVA: no covariates ----
#   print("Bray influence (no covariates)")
#   tempAdonisBray <- adonis2(
#     vegdist(temp_spdata, method = "bray") ~ temp_pseudocatenulatum,
#     data = temp_meta
#   )
#   df_infant_bifido_new[study_name, "BrayAssocPseudocatenulatum"] <- tempAdonisBray$Pr[1]
#   
#   # ---- Kendall PERMANOVA: no covariates ----
#   print("Kendall influence (no covariates)")
#   tempAdonisKendall <- adonis2(
#     as.dist(1 - cor.fk(t(temp_spdata)) / 2) ~ temp_pseudocatenulatum,
#     data = temp_meta
#   )
#   df_infant_bifido_new[study_name, "KendallAssocPseudocatenulatum"] <- tempAdonisKendall$Pr[1]
# }



####################################################################
#####     For Overall Detection in ADULT Study-Cohorts   #######

# df_adult_bifido_new <- as.data.frame(matrix(NA, length(adult_select_studies), 30))
# rownames(df_adult_bifido_new) <- adult_select_studies
# colnames(df_adult_bifido_new) <- c(
#   MajorBifidoSpecies,
#   "BifPrevalence",
#   "BrayAssocAnimalis",
#   "KendallAssocAnimalis",
#   "BrayAssocCatenulatum",
#   "KendallAssocCatenulatum",
#   "BrayAssocBreve",
#   "KendallAssocBreve",
#   "BrayAssocDentium",
#   "KendallAssocDentium",
#   "BrayAssocBifidum",
#   "KendallAssocBifidum",
#   "BrayAssocPseudocatenulatum",
#   "KendallAssocPseudocatenulatum",
#   "BrayAssocAdolescentis",
#   "KendallAssocAdolescentis",
#   "BrayAssocLongum",
#   "KendallAssocLongum",
#   "BrayAssocBifDetection",
#   "KendallAssocBifDetection",
#   "SequencingType",
#   "CohortType",
#   "AgeCategory"
# )
# 
# for (i in seq_along(adult_select_studies)) {
#   study_name <- adult_select_studies[i]
#   print(study_name)
#   
#   ## ---- Get samples for this study ----
#   study_rows <- intersect(
#     adult,
#     rownames(metadata_45809)[metadata_45809$study_name == study_name]
#   )
#   
#   if (length(study_rows) == 0) {
#     warning(study_name, ": no samples found - skipping")
#     next
#   }
#   
#   ## ---- Prevalence of major Bifidobacterium species ----
#   df_adult_bifido_new[study_name, 1:8] <- apply(
#     spdata_45809[study_rows, MajorBifidoSpecies, drop = FALSE],
#     2,
#     function(x) length(x[(x > 0.0001) & !is.na(x)])
#   ) / length(study_rows)
#   
#   ## ---- Overall Bifidobacterium detection prevalence ----
#   df_adult_bifido_new[study_name, "BifPrevalence"] <- length(which(
#     apply(
#       spdata_45809[study_rows, AllBifidoSpecies, drop = FALSE],
#       1,
#       function(x) length(x[(x > 0.0001) & !is.na(x)])
#     ) > 0
#   )) / length(study_rows)
#   
#   ## ---- Subset to non-bifidobacterium species ----
#   temp_spdata <- spdata_45809[study_rows, adult_select_non_bif_species, drop = FALSE]
#   temp_spdata <- apply(temp_spdata, 2, function(x) ifelse(is.na(x), 0, x))
#   temp_spdata <- temp_spdata[rowSums(temp_spdata) > 0, , drop = FALSE]
#   temp_spdata <- temp_spdata[, colSums(temp_spdata) > 0, drop = FALSE]
#   
#   # skip if too few samples or species after filtering
#   if (nrow(temp_spdata) < 2 || ncol(temp_spdata) < 2) {
#     warning(study_name, ": insufficient data after filtering - skipping")
#     next
#   }
#   
#   temp_spdata <- temp_spdata / rowSums(temp_spdata)
#   
#   ## ---- Bifidobacterium detection ----
#   BifidoDetection <- apply(
#     spdata_45809[rownames(temp_spdata), AllBifidoSpecies, drop = FALSE],
#     1,
#     function(x) length(x[x >= 0.0001])
#   )
#   
#   ## ---- Metadata for adonis2: ONLY BifidoDetection ----
#   temp_meta <- data.frame(
#     BifidoDetection = BifidoDetection
#   )
#   
#   ## ---- PERMANOVA: no covariates, only BifidoDetection ----
#   print("Bray influence (no covariates)")
#   tempAdonisBray <- adonis2(
#     vegdist(temp_spdata, method = "bray") ~ BifidoDetection,
#     data = temp_meta
#   )
#   df_adult_bifido_new[study_name, "BrayAssocBifDetection"] <- tempAdonisBray$Pr[1]
#   
#   print("Kendall influence (no covariates)")
#   tempAdonisKendall <- adonis2(
#     as.dist(1 - cor.fk(t(temp_spdata)) / 2) ~ BifidoDetection,
#     data = temp_meta
#   )
#   df_adult_bifido_new[study_name, "KendallAssocBifDetection"] <- tempAdonisKendall$Pr[1]
#   
#   ## ---- Fill SequencingType, CohortType, AgeCategory from metadata_45809 ----
#   # Sequencing type (seq_type)
#   seq_vals <- unique(metadata_45809[study_rows, "seq_type"])
#   seq_vals <- seq_vals[!is.na(seq_vals)]
#   if (length(seq_vals) > 1) {
#     warning(paste("Multiple seq_type values for study:", study_name, "- using first."))
#   }
#   if (length(seq_vals) >= 1) {
#     df_adult_bifido_new[study_name, "SequencingType"] <- seq_vals[1]
#   }
#   
#   # Cohort type (cohort)
#   cohort_vals <- unique(metadata_45809[study_rows, "cohort"])
#   cohort_vals <- cohort_vals[!is.na(cohort_vals)]
#   if (length(cohort_vals) > 1) {
#     warning(paste("Multiple cohort values for study:", study_name, "- using first."))
#   }
#   if (length(cohort_vals) >= 1) {
#     df_adult_bifido_new[study_name, "CohortType"] <- cohort_vals[1]
#   }
#   
#   # Age category (age_category)
#   age_vals <- unique(metadata_45809[study_rows, "age_category"])
#   age_vals <- age_vals[!is.na(age_vals)]
#   if (length(age_vals) > 1) {
#     warning(paste("Multiple age_category values for study:", study_name, "- using first."))
#   }
#   if (length(age_vals) >= 1) {
#     df_adult_bifido_new[study_name, "AgeCategory"] <- age_vals[1]
#   }
# }

####### For Eight different Bifs in Adult  ######

# for (i in seq_along(adult_select_studies)) {
#   study_name <- adult_select_studies[i]
#   print(study_name)
#   
#   # ---- Get samples for this study ----
#   study_rows <- intersect(
#     adult,
#     rownames(metadata_45809)[metadata_45809$study_name == study_name]
#   )
#   
#   if (length(study_rows) == 0) {
#     warning(study_name, ": no samples found - skipping")
#     next
#   }
#   
#   # ---- Subset to non-bifidobacterium species ----
#   temp_spdata <- spdata_45809[study_rows, adult_select_non_bif_species, drop = FALSE]
#   temp_spdata <- apply(temp_spdata, 2, function(x) ifelse(is.na(x), 0, x))
#   temp_spdata <- temp_spdata[rowSums(temp_spdata) > 0, , drop = FALSE]
#   temp_spdata <- temp_spdata[, colSums(temp_spdata) > 0, drop = FALSE]
#   
#   # skip if too few samples or species after filtering
#   if (nrow(temp_spdata) < 3 || ncol(temp_spdata) < 2) {
#     warning(study_name, ": insufficient data after filtering - skipping")
#     next
#   }
#   
#   temp_spdata <- temp_spdata / rowSums(temp_spdata)
#   
#   # ---- B. longum abundance ----
#   temp_animalis <- spdata_45809[rownames(temp_spdata), "Bifidobacterium_animalis"]
#   
#   # ---- Metadata for adonis2: ONLY temp_longum ----
#   temp_meta <- data.frame(
#     temp_animalis = temp_animalis
#   )
#   
#   # ---- Bray-Curtis PERMANOVA (no covariates) ----
#   print("Bray influence (no covariates)")
#   tempAdonisBray <- adonis2(
#     vegdist(temp_spdata, method = "bray") ~ temp_animalis,
#     data = temp_meta
#   )
#   df_adult_bifido_new[study_name, "BrayAssocAnimalis"] <- tempAdonisBray$Pr[1]
#   
#   # ---- Kendall PERMANOVA (no covariates) ----
#   print("Kendall influence (no covariates)")
#   tempAdonisKendall <- adonis2(
#     as.dist(1 - cor.fk(t(temp_spdata)) / 2) ~ temp_animalis,
#     data = temp_meta
#   )
#   df_adult_bifido_new[study_name, "KendallAssocAnimalis"] <- tempAdonisKendall$Pr[1]
# }
# 
# save(df_adult_bifido_new, file = "/data/df_adult_bifido_new.RData")




#######################################################################
################### SENIOR ####################

# df_senior_bifido_new <- as.data.frame(matrix(NA, length(senior_select_studies), 30))
# rownames(df_senior_bifido_new) <- senior_select_studies
# colnames(df_senior_bifido_new) <- c(
#   MajorBifidoSpecies,
#   "BifPrevalence",
#   "BrayAssocAnimalis",
#   "KendallAssocAnimalis",
#   "BrayAssocCatenulatum",
#   "KendallAssocCatenulatum",
#   "BrayAssocBreve",
#   "KendallAssocBreve",
#   "BrayAssocDentium",
#   "KendallAssocDentium",
#   "BrayAssocBifidum",
#   "KendallAssocBifidum",
#   "BrayAssocPseudocatenulatum",
#   "KendallAssocPseudocatenulatum",
#   "BrayAssocAdolescentis",
#   "KendallAssocAdolescentis",
#   "BrayAssocLongum",
#   "KendallAssocLongum",
#   "BrayAssocBifDetection",
#   "KendallAssocBifDetection",
#   "SequencingType",
#   "CohortType",
#   "AgeCategory"
# )
# 
# for (i in seq_along(senior_select_studies)) {
#   study_name <- senior_select_studies[i]
#   print(study_name)
#   
#   ## ---- Get samples for this study ----
#   study_rows <- intersect(
#     senior,
#     rownames(metadata_45809)[metadata_45809$study_name == study_name]
#   )
#   
#   if (length(study_rows) == 0) {
#     warning(study_name, ": no samples found - skipping")
#     next
#   }
#   
#   ## ---- Prevalence of major Bifidobacterium species ----
#   df_senior_bifido_new[study_name, 1:8] <- apply(
#     spdata_45809[study_rows, MajorBifidoSpecies, drop = FALSE],
#     2,
#     function(x) length(x[(x > 0.0001) & !is.na(x)])
#   ) / length(study_rows)
#   
#   ## ---- Overall Bifidobacterium detection prevalence ----
#   df_senior_bifido_new[study_name, "BifPrevalence"] <- length(which(
#     apply(
#       spdata_45809[study_rows, AllBifidoSpecies, drop = FALSE],
#       1,
#       function(x) length(x[(x > 0.0001) & !is.na(x)])
#     ) > 0
#   )) / length(study_rows)
#   
#   ## ---- Subset to non-bifidobacterium species ----
#   temp_spdata <- spdata_45809[study_rows, senior_select_non_bif_species, drop = FALSE]
#   temp_spdata <- apply(temp_spdata, 2, function(x) ifelse(is.na(x), 0, x))
#   temp_spdata <- temp_spdata[rowSums(temp_spdata) > 0, , drop = FALSE]
#   temp_spdata <- temp_spdata[, colSums(temp_spdata) > 0, drop = FALSE]
#   
#   # skip if too few samples or species after filtering
#   if (nrow(temp_spdata) < 2 || ncol(temp_spdata) < 2) {
#     warning(study_name, ": insufficient data after filtering - skipping")
#     next
#   }
#   
#   temp_spdata <- temp_spdata / rowSums(temp_spdata)
#   
#   ## ---- Bifidobacterium detection ----
#   BifidoDetection <- apply(
#     spdata_45809[rownames(temp_spdata), AllBifidoSpecies, drop = FALSE],
#     1,
#     function(x) length(x[x >= 0.0001])
#   )
#   
#   ## ---- Metadata for adonis2: ONLY BifidoDetection ----
#   temp_meta <- data.frame(
#     BifidoDetection = BifidoDetection
#   )
#   
#   ## ---- PERMANOVA: no covariates, only BifidoDetection ----
#   print("Bray influence (no covariates)")
#   tempAdonisBray <- adonis2(
#     vegdist(temp_spdata, method = "bray") ~ BifidoDetection,
#     data = temp_meta
#   )
#   df_senior_bifido_new[study_name, "BrayAssocBifDetection"] <- tempAdonisBray$Pr[1]
#   
#   print("Kendall influence (no covariates)")
#   tempAdonisKendall <- adonis2(
#     as.dist(1 - cor.fk(t(temp_spdata)) / 2) ~ BifidoDetection,
#     data = temp_meta
#   )
#   df_senior_bifido_new[study_name, "KendallAssocBifDetection"] <- tempAdonisKendall$Pr[1]
#   
#   ## ---- Fill SequencingType, CohortType, AgeCategory from metadata_45809 ----
#   # Sequencing type (seq_type)
#   seq_vals <- unique(metadata_45809[study_rows, "seq_type"])
#   seq_vals <- seq_vals[!is.na(seq_vals)]
#   if (length(seq_vals) > 1) {
#     warning(paste("Multiple seq_type values for study:", study_name, "- using first."))
#   }
#   if (length(seq_vals) >= 1) {
#     df_senior_bifido_new[study_name, "SequencingType"] <- seq_vals[1]
#   }
#   
#   # Cohort type (cohort)
#   cohort_vals <- unique(metadata_45809[study_rows, "cohort"])
#   cohort_vals <- cohort_vals[!is.na(cohort_vals)]
#   if (length(cohort_vals) > 1) {
#     warning(paste("Multiple cohort values for study:", study_name, "- using first."))
#   }
#   if (length(cohort_vals) >= 1) {
#     df_senior_bifido_new[study_name, "CohortType"] <- cohort_vals[1]
#   }
#   
#   # Age category (age_category)
#   age_vals <- unique(metadata_45809[study_rows, "age_category"])
#   age_vals <- age_vals[!is.na(age_vals)]
#   if (length(age_vals) > 1) {
#     warning(paste("Multiple age_category values for study:", study_name, "- using first."))
#   }
#   if (length(age_vals) >= 1) {
#     df_senior_bifido_new[study_name, "AgeCategory"] <- age_vals[1]
#   }
# }
# 
# 
# 
# ####### Eight different Bifs ######
# 
# for (i in seq_along(senior_select_studies)) {
#   study_name <- senior_select_studies[i]
#   print(study_name)
#   
#   # ---- Get samples for this study ----
#   study_rows <- intersect(
#     senior,
#     rownames(metadata_45809)[metadata_45809$study_name == study_name]
#   )
#   
#   if (length(study_rows) == 0) {
#     warning(study_name, ": no samples found - skipping")
#     next
#   }
#   
#   # ---- Subset to non-bifidobacterium species ----
#   temp_spdata <- spdata_45809[study_rows, senior_select_non_bif_species, drop = FALSE]
#   temp_spdata <- apply(temp_spdata, 2, function(x) ifelse(is.na(x), 0, x))
#   temp_spdata <- temp_spdata[rowSums(temp_spdata) > 0, , drop = FALSE]
#   temp_spdata <- temp_spdata[, colSums(temp_spdata) > 0, drop = FALSE]
#   
#   # skip if too few samples or species after filtering
#   if (nrow(temp_spdata) < 3 || ncol(temp_spdata) < 2) {
#     warning(study_name, ": insufficient data after filtering - skipping")
#     next
#   }
#   
#   temp_spdata <- temp_spdata / rowSums(temp_spdata)
#   
#   # ---- B. longum abundance ----
#   temp_longum <- spdata_45809[rownames(temp_spdata), "Bifidobacterium_longum"]
#   
#   # ---- Metadata for adonis2: ONLY temp_longum ----
#   temp_meta <- data.frame(
#     temp_longum = temp_longum
#   )
#   
#   # ---- Bray-Curtis PERMANOVA (no covariates) ----
#   print("Bray influence (no covariates)")
#   tempAdonisBray <- adonis2(
#     vegdist(temp_spdata, method = "bray") ~ temp_longum,
#     data = temp_meta
#   )
#   df_senior_bifido_new[study_name, "BrayAssocLongum"] <- tempAdonisBray$Pr[1]
#   
#   # ---- Kendall PERMANOVA (no covariates) ----
#   print("Kendall influence (no covariates)")
#   tempAdonisKendall <- adonis2(
#     as.dist(1 - cor.fk(t(temp_spdata)) / 2) ~ temp_longum,
#     data = temp_meta
#   )
#   df_senior_bifido_new[study_name, "KendallAssocLongum"] <- tempAdonisKendall$Pr[1]
# }


#################################################################
#################################################################

load("/data/df_infant_bifido_new.RData")
load("/data/df_adult_bifido_new.RData")
load("/data/df_senior_bifido_new.RData")


df_infant_bifido_new[,"AgeCategory"] <- "Infant"
df_adult_bifido_new[,"AgeCategory"] <- "Adult"
df_senior_bifido_new[,"AgeCategory"] <- "Senior"

rownames(df_infant_bifido_new) <- paste0("infant:",rownames(df_infant_bifido_new))
rownames(df_adult_bifido_new) <- paste0("adult:",rownames(df_adult_bifido_new))
rownames(df_senior_bifido_new) <- paste0("senior:",rownames(df_senior_bifido_new))

df_combined_bifido <- as.data.frame(rbind(df_infant_bifido_new,df_adult_bifido_new,df_senior_bifido_new))

#save(df_combined_bifido, file="/data/df_combined_AgeCategories_AllAdjusted_new.RData")














