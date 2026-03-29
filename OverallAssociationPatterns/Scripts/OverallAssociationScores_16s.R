load("G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\QMP_analysis\\QMP_workspace.RData")
load("/data/bifido_45809_new_2014_not_normalized_Harnandez_added.RData")
load("/data/full_AssociationScoresAgeCategory_16s_original.RData")
load("/data/full_AssociationScoresCohortType_16s_original.RData")

library(dplyr)
library(pcaPP)
library(psych)
library(RColorBrewer)
library(metap)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ppcor)
library(vegan)

grep("_16s$",ls(),value = T)

SpeciesProfile_All[is.na(SpeciesProfile_All)] <- 0

spdata_16s <- SpeciesProfile_All[rownames(metadata_updated[metadata_updated$seq_type == '16s',]),]
metadata_16s <- metadata_updated[rownames(spdata_16s),]
metadata_16s$cohort <- ifelse(metadata_16s$cohort == "IndustrializedUrban","IndustrializedUrban","Others")

spdata_16s$Bifidobacterium_detection <- apply(spdata_16s[, AllBifidoSpecies], 1, function(x) length(x[x >= 0.0001]))
all(rownames(spdata_16s) == rownames(metadata_16s))
spdata_16s$age_category <- metadata_16s$age_category
spdata_16s$cohort <- metadata_16s$cohort
spdata_16s$study_name <- metadata_16s$study_name

names(spdata_16s)[1265] <- "Intestinibacter_bartlettii"
names(spdata_16s)[1362] <- "Erysipelatoclostridium_ramosum"


##################### Study Filtration ################

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

refined_data_infant <- spdata_16s[spdata_16s$age_category == 'infant',]
refined_data_adult <- spdata_16s[spdata_16s$age_category == 'adult',]
refined_data_senior <- spdata_16s[spdata_16s$age_category == 'senior',]
bifido_species <- c("Bifidobacterium_detection","Bifidobacterium_longum","Bifidobacterium_adolescentis","Bifidobacterium_pseudocatenulatum","Bifidobacterium_dentium","Bifidobacterium_bifidum","Bifidobacterium_breve","Bifidobacterium_catenulatum","Bifidobacterium_animalis")

MajorBifsDetectionInfant <- compute_detection(refined_data_infant,bifido_species,"study_name",unique(refined_data_infant$study_name))
Bifidobacterium_catenulatum_infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_catenulatum",]>=0.10)]
Bifidobacterium_breve_infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_breve",]>=0.10)]
Bifidobacterium_dentium_infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_dentium",]>=0.10)]
Bifidobacterium_bifidum_infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_bifidum",]>=0.10)]
Bifidobacterium_pseudocatenulatum_infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_pseudocatenulatum",]>=0.10)]
Bifidobacterium_adolescentis_infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_adolescentis",]>=0.10)]
Bifidobacterium_longum_infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_longum",]>=0.10)]
Bifidobacterium_detection_infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_detection",]>=0.10)]

MajorBifsDetectionAdult <- compute_detection(refined_data_adult,bifido_species,"study_name",unique(refined_data_adult$study_name))
Bifidobacterium_catenulatum_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_catenulatum",]>=0.10)]
Bifidobacterium_breve_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_breve",]>=0.10)]
Bifidobacterium_dentium_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_dentium",]>=0.10)]
Bifidobacterium_bifidum_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_bifidum",]>=0.10)]
Bifidobacterium_pseudocatenulatum_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_pseudocatenulatum",]>=0.10)]
Bifidobacterium_adolescentis_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_adolescentis",]>=0.10)]
Bifidobacterium_longum_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_longum",]>=0.10)]
Bifidobacterium_animalis_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_animalis",]>=0.10)]
Bifidobacterium_detection_adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_detection",]>=0.10)]

MajorBifsDetectionSenior <- compute_detection(refined_data_senior,bifido_species,"study_name",unique(refined_data_senior$study_name))
Bifidobacterium_catenulatum_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_catenulatum",]>=0.10)]
Bifidobacterium_breve_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_breve",]>=0.10)]
Bifidobacterium_dentium_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_dentium",]>=0.10)]
Bifidobacterium_bifidum_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_bifidum",]>=0.10)]
Bifidobacterium_pseudocatenulatum_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_pseudocatenulatum",]>=0.10)]
Bifidobacterium_adolescentis_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_adolescentis",]>=0.10)]
Bifidobacterium_longum_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_longum",]>=0.10)]
Bifidobacterium_animalis_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_animalis",]>=0.10)]
Bifidobacterium_detection_senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_detection",]>=0.10)]



######################################################################
###########    Association Calculation    ###############

############## AGE CATEGORY Wise ############

compute_corr_associations_by_AgeCategory <- function(spdata, df_association, age_category, bif_name, bif_studies) {
  
  spdata_age <- spdata[spdata$age_category == age_category & spdata$study_name %in% bif_studies, ]
  
  studies <- intersect(unique(spdata_age$study_name), bif_studies)
  
  non_bif_species <- rownames(df_association)
  
  corr_mat <- matrix(NA, nrow = length(non_bif_species), ncol = length(studies),
                     dimnames = list(non_bif_species, studies))
  pval_mat <- corr_mat
  dir_mat  <- corr_mat
  
  for (study in studies) {
    study_df <- subset(spdata_age, study_name == study)
    if (nrow(study_df) < 5) next
    
    for (sp in non_bif_species) {
      if (!all(c(sp, bif_name) %in% colnames(study_df))) next
      
      x <- study_df[[sp]]
      y <- study_df[[bif_name]]
      
      valid_idx <- complete.cases(x, y)
      if (sum(valid_idx) < 5) next
      
      test_res <- try(cor.test(x[valid_idx], y[valid_idx], method = "spearman", exact = TRUE), silent = TRUE)
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
  
  positive <- apply(dir_mat, 1, function(x) sum(x %in% c(1, 2)))
  negative <- apply(dir_mat, 1, function(x) sum(x %in% c(-1, -2)))
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


df_association_infant_detection_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_infant_detection_16s,"infant","Bifidobacterium_detection",Bifidobacterium_detection_infant_studies)
df_association_infant_longum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_infant_longum_16s,"infant","Bifidobacterium_longum",Bifidobacterium_longum_infant_studies)
df_association_infant_adolescentis_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_infant_adolescentis_16s,"infant","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_infant_studies)
df_association_infant_pseudocatenulatum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_infant_pseudocatenulatum_16s,"infant","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_infant_studies)
df_association_infant_dentium_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_infant_dentium_16s,"infant","Bifidobacterium_dentium",Bifidobacterium_dentium_infant_studies)
df_association_infant_bifidum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_infant_bifidum_16s,"infant","Bifidobacterium_bifidum",Bifidobacterium_bifidum_infant_studies)
df_association_infant_breve_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_infant_breve_16s,"infant","Bifidobacterium_breve",Bifidobacterium_breve_infant_studies)
df_association_infant_catenulatum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_infant_catenulatum_16s,"infant","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_infant_studies)


# df_association_adult_detection_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_detection_16s,"adult","Bifidobacterium_detection",Bifidobacterium_detection_adult_studies)
# df_association_adult_longum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_longum_16s,"adult","Bifidobacterium_longum",Bifidobacterium_longum_adult_studies)
# df_association_adult_adolescentis_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_adolescentis_16s,"adult","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_adult_studies)
# df_association_adult_pseudocatenulatum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_pseudocatenulatum_16s,"adult","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_adult_studies)
# df_association_adult_dentium_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_dentium_16s,"adult","Bifidobacterium_dentium",Bifidobacterium_dentium_adult_studies)
# df_association_adult_bifidum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_bifidum_16s,"adult","Bifidobacterium_bifidum",Bifidobacterium_bifidum_adult_studies)
# df_association_adult_breve_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_breve_16s,"adult","Bifidobacterium_breve",Bifidobacterium_breve_adult_studies)
# df_association_adult_catenulatum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_catenulatum_16s,"adult","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_adult_studies)
# df_association_adult_animalis_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_adult_animalis_16s,"adult","Bifidobacterium_animalis",Bifidobacterium_animalis_adult_studies)
# 
# 
# df_association_senior_detection_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_detection_16s,"senior","Bifidobacterium_detection",Bifidobacterium_detection_senior_studies)
# df_association_senior_longum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_longum_16s,"senior","Bifidobacterium_longum",Bifidobacterium_longum_senior_studies)
# df_association_senior_adolescentis_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_adolescentis_16s,"senior","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_senior_studies)
# df_association_senior_pseudocatenulatum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_pseudocatenulatum_16s,"senior","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_senior_studies)
# df_association_senior_dentium_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_dentium_16s,"senior","Bifidobacterium_dentium",Bifidobacterium_dentium_senior_studies)
# df_association_senior_bifidum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_bifidum_16s,"senior","Bifidobacterium_bifidum",Bifidobacterium_bifidum_senior_studies)
# df_association_senior_breve_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_breve_16s,"senior","Bifidobacterium_breve",Bifidobacterium_breve_senior_studies)
# df_association_senior_catenulatum_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_catenulatum_16s,"senior","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_senior_studies)
# df_association_senior_animalis_16s_new <- compute_corr_associations_by_AgeCategory(spdata_16s, df_association_senior_animalis_16s,"senior","Bifidobacterium_animalis",Bifidobacterium_animalis_senior_studies)
# 
# 
# ###########  Cohort Lifestyle Wise #########################
# 
# refined_data_IndustrializedUrban <- spdata_16s[spdata_16s$age_category == 'adult' & spdata_16s$cohort == 'IndustrializedUrban',]
# refined_data_Others <- spdata_16s[spdata_16s$age_category == 'adult' & spdata_16s$cohort == 'Others',]
# 
# MajorBifsDetectionIndustrializedUrban <- compute_detection(refined_data_IndustrializedUrban,bifido_species,"study_name",unique(refined_data_IndustrializedUrban$study_name))
# Bifidobacterium_catenulatum_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_catenulatum",]>=0.10)]
# Bifidobacterium_breve_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_breve",]>=0.10)]
# Bifidobacterium_dentium_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_dentium",]>=0.10)]
# Bifidobacterium_bifidum_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_bifidum",]>=0.10)]
# Bifidobacterium_pseudocatenulatum_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_pseudocatenulatum",]>=0.10)]
# Bifidobacterium_adolescentis_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_adolescentis",]>=0.10)]
# Bifidobacterium_longum_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_longum",]>=0.10)]
# Bifidobacterium_detection_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_detection",]>=0.10)]
# Bifidobacterium_animalis_IndustrializedUrban_studies <- colnames(MajorBifsDetectionIndustrializedUrban)[which(MajorBifsDetectionIndustrializedUrban["Bifidobacterium_animalis",]>=0.10)]
# 
# MajorBifsDetectionOthers <- compute_detection(refined_data_Others,bifido_species,"study_name",unique(refined_data_Others$study_name))
# Bifidobacterium_catenulatum_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_catenulatum",]>=0.10)]
# Bifidobacterium_breve_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_breve",]>=0.10)]
# Bifidobacterium_dentium_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_dentium",]>=0.10)]
# Bifidobacterium_bifidum_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_bifidum",]>=0.10)]
# Bifidobacterium_pseudocatenulatum_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_pseudocatenulatum",]>=0.10)]
# Bifidobacterium_adolescentis_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_adolescentis",]>=0.10)]
# Bifidobacterium_longum_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_longum",]>=0.10)]
# Bifidobacterium_animalis_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_animalis",]>=0.10)]
# Bifidobacterium_detection_Others_studies <- colnames(MajorBifsDetectionOthers)[which(MajorBifsDetectionOthers["Bifidobacterium_detection",]>=0.10)]
# 
# 
# ###################################################################################
# 
# compute_corr_associations_by_cohort <- function(spdata, df_association, cohort_name, bif_name, bif_studies) {
#   
#   spdata_cohort <- spdata[spdata$age_category == "adult" & spdata$cohort %in% cohort_name, ]
#   
#   if (nrow(spdata_cohort) == 0) {
#     stop(paste("No samples found for cohort:", cohort_name))
#   }
#   
#   studies <- bif_studies
#   
#   non_bif_species <- rownames(df_association)
#   
#   corr_mat <- matrix(NA, nrow = length(non_bif_species), ncol = length(studies),
#                      dimnames = list(non_bif_species, studies))
#   pval_mat <- corr_mat
#   dir_mat  <- corr_mat
#   
#   for (study in studies) {
#     study_df <- subset(spdata_cohort, study_name == study)
#     if (nrow(study_df) < 5) next
#     
#     for (sp in non_bif_species) {
#       if (!all(c(sp, bif_name) %in% colnames(study_df))) next
#       
#       x <- study_df[[sp]]
#       y <- study_df[[bif_name]]
#       
#       if (sum(complete.cases(x, y)) < 5) next
#       
#       test_res <- try(cor.test(x, y, method = "spearman", exact = TRUE), silent = TRUE)
#       if (inherits(test_res, "try-error")) next
#       
#       corr_mat[sp, study] <- test_res$estimate
#       pval_mat[sp, study] <- test_res$p.value
#     }
#   }
#   
#   corr_mat[is.na(corr_mat)] <- 0
#   pval_mat[is.na(pval_mat)] <- 1
#   
#   dir_mat[corr_mat > 0 & pval_mat <= 0.05] <- 2
#   dir_mat[corr_mat > 0 & pval_mat > 0.05 & pval_mat <= 0.1] <- 1
#   dir_mat[corr_mat < 0 & pval_mat <= 0.05] <- -2
#   dir_mat[corr_mat < 0 & pval_mat > 0.05 & pval_mat <= 0.1] <- -1
#   dir_mat[is.na(dir_mat)] <- 0
#   
#   positive <- apply(dir_mat, 1, function(x) sum(x %in% c(1, 2)))
#   negative <- apply(dir_mat, 1, function(x) sum(x %in% c(-1, -2)))
#   total <- length(studies)
#   
#   score <- ((positive - negative) / total) *
#     (1 - ((pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
#   
#   association <- data.frame(
#     positive = positive,
#     negative = negative,
#     total = total,
#     score = score,
#     row.names = non_bif_species
#   )
#   
#   output_list <- list(
#     corr = as.data.frame(corr_mat),
#     pval = as.data.frame(pval_mat),
#     dir = as.data.frame(dir_mat),
#     association = association
#   )
#   
#   return(output_list)
# }
# 
# 
# df_association_IndustrializedUrban_detection_16s_new <- compute_corr_associations_by_cohort(spdata_16s,df_association_IndustrializedUrban_detection_16s,"IndustrializedUrban","Bifidobacterium_detection",Bifidobacterium_detection_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_longum_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_IndustrializedUrban_longum_16s,"IndustrializedUrban","Bifidobacterium_longum",Bifidobacterium_longum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_adolescentis_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_IndustrializedUrban_adolescentis_16s,"IndustrializedUrban","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_pseudocatenulatum_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_IndustrializedUrban_pseudocatenulatum_16s,"IndustrializedUrban","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_dentium_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_IndustrializedUrban_dentium_16s,"IndustrializedUrban","Bifidobacterium_dentium",Bifidobacterium_dentium_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_bifidum_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_IndustrializedUrban_bifidum_16s,"IndustrializedUrban","Bifidobacterium_bifidum",Bifidobacterium_bifidum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_breve_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_IndustrializedUrban_breve_16s,"IndustrializedUrban","Bifidobacterium_breve",Bifidobacterium_breve_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_catenulatum_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_IndustrializedUrban_catenulatum_16s,"IndustrializedUrban","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_animalis_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_IndustrializedUrban_animalis_16s,"IndustrializedUrban","Bifidobacterium_animalis",Bifidobacterium_animalis_IndustrializedUrban_studies)
# 
# df_association_Others_detection_16s_new <- compute_corr_associations_by_cohort(spdata_16s,df_association_Others_detection_16s,"Others","Bifidobacterium_detection",Bifidobacterium_detection_Others_studies)
# df_association_Others_longum_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_Others_longum_16s,"Others","Bifidobacterium_longum",Bifidobacterium_longum_Others_studies)
# df_association_Others_adolescentis_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_Others_adolescentis_16s,"Others","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_Others_studies)
# df_association_Others_pseudocatenulatum_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_Others_pseudocatenulatum_16s,"Others","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_Others_studies)
# df_association_Others_dentium_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_Others_dentium_16s,"Others","Bifidobacterium_dentium",Bifidobacterium_dentium_Others_studies)
# df_association_Others_bifidum_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_Others_bifidum_16s,"Others","Bifidobacterium_bifidum",Bifidobacterium_bifidum_Others_studies)
# df_association_Others_breve_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_Others_breve_16s,"Others","Bifidobacterium_breve",Bifidobacterium_breve_Others_studies)
# df_association_Others_catenulatum_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_Others_catenulatum_16s,"Others","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_Others_studies)
# df_association_Others_animalis_16s_new <- compute_corr_associations_by_cohort(spdata_16s, df_association_Others_animalis_16s,"Others","Bifidobacterium_animalis",Bifidobacterium_animalis_Others_studies)
# 
# 
# ####################################################################
# 
# adjusted_lists <- ls(pattern = "_new$")
# save(list = adjusted_lists, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\AssociationScores_16s\\Overall_AssociationScores_16s_new.RData")
# 
# 
# list_names <- ls()[sapply(ls(), function(x) is.list(get(x)))]
# 
# for (lname in list_names) {
#   lst <- get(lname)
#   
#   if (!is.null(lst$association)) { 
#     assign(lname, lst$association)  
#   } else {
#     warning(paste("List", lname, "does not contain 'association'"))
#   }
# }
# 
# save(list = list_names, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\AssociationScores_16s\\df_AssociationScores_new_16s.RData")

##########################################################################
# 
# load("G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\AssociationScores_16s\\df_AssociationScores_original_16s.RData")
# load("G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\AssociationScores_16s\\df_AssociationScores_new_16s.RData")
# 
# corr_result <- data.frame(
#   pair = character(),
#   corr_orig_new = numeric(),
#   pval_orig_new = numeric(),
#   stringsAsFactors = FALSE)
# 
# 
# original_dfs <- ls(pattern = "_original$")
# 
# for (orig_name in original_dfs) {
#   new_name <- sub("_original$", "_new", orig_name)
#   
#   # Check if the corresponding _new version exists
#   if (!exists(new_name)) next
#   
#   df_orig <- get(orig_name)
#   df_new  <- get(new_name)
#   
#   # Ensure "score" column exists in both
#   if (!("score" %in% colnames(df_orig)) || !("score" %in% colnames(df_new))) {
#     warning(paste("Missing 'score' column in", orig_name, "or", new_name))
#     next
#   }
#   
#   # Safe correlation computation
#   test_res <- try(cor.test(df_orig$score, df_new$score, method = "spearman", exact = TRUE), silent = TRUE)
#   
#   if (inherits(test_res, "try-error")) next
#   
#   pair_name <- gsub("^df_association_|_wgs.*$", "", orig_name)  # extract readable pair name
#   
#   corr_result <- rbind(corr_result, data.frame(
#     pair = pair_name,
#     corr_orig_new = test_res$estimate,
#     pval_orig_new = test_res$p.value,
#     stringsAsFactors = FALSE
#   ))
# }
# 
# rownames(corr_result) <- corr_result$pair
# corr_result$pair <- NULL
# 
# df <- corr_result
# df <- as.data.frame(df)
# rownames(df) <- gsub("_16s_original","",rownames(df))
# 
# # Define desired order
# categories <- c("infant", "adult", "senior", 
#                 "IndustrializedUrban", "Others")
# bifnames <- c("detection", "longum", "adolescentis", "pseudocatenulatum",
#               "dentium", "bifidum", "breve", "catenulatum", "animalis")
# 
# # Generate all desired combinations
# desired_order <- unlist(lapply(categories, function(cat) {
#   paste0(cat, "_", bifnames)
# }))
# 
# # Remove infant_animalis (since it doesn’t exist)
# desired_order <- desired_order[!desired_order %in% "infant_animalis"]
# 
# # Now reorder based on this sequence (keeping only those that exist in your data)
# available_order <- desired_order[desired_order %in% rownames(df)]
# 
# # Reorder
# df_reordered <- df[available_order, ]
# 
# 
# library(xlsx)
# write.xlsx(df_reordered,"G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\AssociationScores_16s\\Original_vs_New_Corr_16s.xlsx")
# 

####################################################################
### Strong Association Score Calculation for 16S ####

####### Considering Strong Associations ######

load("G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\MANUSCRIPT_ITEMS\\REVISED_FINAL\\CodeOcean_GitHub\\OverallAssociationPatterns\\DATA\\Overall_AssociationScores_16s_new.RData")

obj_names <- ls(pattern = "^df_association_")

for (obj in obj_names) {
  
  # Extract the list
  lst <- get(obj)
  
  # Extract the dir matrix
  dir_mat <- as.matrix(lst$dir)
  
  # Count strong positives (only 2)
  positive <- apply(dir_mat, 1, function(x) sum(x == 2, na.rm = TRUE))
  
  # Count strong negatives (only -2)
  negative <- apply(dir_mat, 1, function(x) sum(x == -2, na.rm = TRUE))
  
  # Total number of studies (columns)
  total <- ncol(dir_mat)
  
  # Association score
  score <- ((positive - negative) / total) *
    (1 - ((pmin(positive, negative) + 0.00001) /
            (pmax(positive, negative) + 0.00001)))
  
  # Create strong_association dataframe
  strong_association <- data.frame(
    positive = positive,
    negative = negative,
    total = total,
    score = score,
    row.names = rownames(dir_mat)
  )
  
  # Add to list
  lst$strong_association <- strong_association
  
  # Save back to environment
  assign(obj, lst)
}

adjusted_lists <- ls(pattern = "_new$")
save(list = adjusted_lists, file = "/results/Overall_Strong_AssociationScores_16s_new.RData")

#####################################################################

load("/data/Overall_Strong_AssociationScores_16s_new.RData")

list_names <- ls()[sapply(ls(), function(x) is.list(get(x)))]

for (lname in list_names) {
  lst <- get(lname)
  
  if (!is.null(lst$strong_association)) { 
    assign(lname, lst$strong_association)  
  } else {
    warning(paste("List", lname, "does not contain 'association'"))
  }
}

save(list = list_names, file = "/data/df_Strong_AssociationScores_new_16s.RData")

##################################################################
### Combined Association Score df for 16S ###

df_strong_association_all_16s <- read.xlsx("df_strong_association_all_16s_filtered.xlsx", sheet = 'Sheet1')
rownames(df_strong_association_all_16s) <- df_strong_association_all_16s$c..Acetanaerobacterium_elongatum....Acidaminococcus_fermentans...
df_strong_association_all_16s$c..Acetanaerobacterium_elongatum....Acidaminococcus_fermentans... <- NULL

col_fun <- colorRamp2(
  c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
  c("#8B0053", "#CC3377", "#F6B6D2", "#FFFFFF", "#CDE9B6", "#77C679", "#004D00")
)

# Convert to numeric matrix
mat <- as.matrix(df_strong_association_all_16s)

ht <- Heatmap(
  mat,
  name = "Association Scores",
  col = col_fun,
  cluster_rows = TRUE,                   
  cluster_columns = FALSE,               
  show_row_dend = FALSE,                 
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7.7),
  column_names_gp = gpar(fontsize = 7.7),
  heatmap_legend_param = list(
    at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
    labels = c("-1", "-0.75", "-0.5", "-0.25", "0", 
               "0.25", "0.5", "0.75", "1"),
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

pdf("/results/StrongAssociationScores_16s.pdf", 
    height = 25, width = 16)
draw(ht)
dev.off()

###### CARPET ######

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_AssociationScores_Heatmap <- as.data.frame(carpet_df)

#write.xlsx(carpet_AssociationScores_Heatmap, "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\AssociationScores_16s\\carpet_StrongAssociationScores_16s.xlsx")

save(carpet_AssociationScores_Heatmap,file="/results/df_strong_association_all_new_16s.RData")



###################################################################
######## Shannon Association Scores across Age groups ###########

Bif_list <- rev(MajorBifsOrdered)
SpeciesProfile_All[is.na(SpeciesProfile_All)] <- 0
spdata_16s <- SpeciesProfile_All[rownames(SpeciesProfile_All) %in% rownames(metadata_updated[metadata_updated$seq_type == '16s',]),]

Bif_SpeciesProf <- spdata_16s[,Bif_list]
Bif_SpeciesProf$Bifidobacterium_detection <- apply(SpeciesProfile_All[rownames(Bif_SpeciesProf), AllBifidoSpecies], 1, function(x) length(x[x >= 0.0001]))

load("/data/df_strong_association_all_new_16s.RData")

metadata_16s <- metadata_updated[rownames(spdata_16s),]
metadata_16s$cohort <- ifelse(metadata_16s$cohort == "IndustrializedUrban","IndustrializedUrban","Others")

#spdata_16s$Bifidobacterium_detection <- apply(spdata_16s[, AllBifidoSpecies], 1, function(x) length(x[x >= 0.0001]))
all(rownames(spdata_16s) == rownames(metadata_16s))
#spdata_16s$age_category <- metadata_16s$age_category
#spdata_16s$cohort <- metadata_16s$cohort
#spdata_16s$study_name <- metadata_16s$study_name

names(spdata_16s)[1265] <- "Intestinibacter_bartlettii"
names(spdata_16s)[1362] <- "Erysipelatoclostridium_ramosum"

nonbif_SpeciesProfile_16s <- spdata_16s[,names(spdata_16s) %in% rownames(carpet_AssociationScores_Heatmap)]

Bif_SpeciesProf$shannon_index <- diversity(nonbif_SpeciesProfile_16s, index = 'shannon')

calculate_pielou <- function(abundance) {
  H <- diversity(abundance, index = "shannon") 
  S <- sum(abundance > 0)                     
  if (S > 1) {
    E <- H / log(S)                           
  } else {
    E <- NA                                 
  }
  return(E)
}


Bif_SpeciesProf$pielou <- calculate_pielou(nonbif_SpeciesProfile_16s)

all(rownames(Bif_SpeciesProf) == rownames(metadata_16s))
Bif_SpeciesProf$age_category <- metadata_16s$age_category
Bif_SpeciesProf$study_name <- metadata_16s$study_name
Bif_SpeciesProf$cohort <- metadata_16s$cohort

shannon_pielou_Bif_df_16s <- as.data.frame(Bif_SpeciesProf)

# Create list to store results for each age category
age_category_list_shannon <- list()

# Loop through each unique age category
for (age_cat in unique(shannon_pielou_Bif_df_16s$age_category)) {
  
  temp_df <- subset(shannon_pielou_Bif_df_16s, age_category == age_cat)
  bifido_species <- colnames(temp_df)[1:9]
  studies <- unique(temp_df$study_name)
  
  corr_df <- matrix(NA, nrow = length(studies), ncol = length(bifido_species),
                    dimnames = list(studies, bifido_species))
  pval_df <- corr_df
  
  # Loop through each study
  for (study in studies) {
    study_df <- subset(temp_df, study_name == study)
    
    for (sp in bifido_species) {
      if (length(unique(study_df[[sp]])) > 1 && length(unique(study_df$shannon_index)) > 1) {
        res <- suppressWarnings(cor.test(study_df[[sp]], study_df$shannon_index,
                                         method = "spearman", exact = TRUE))
        corr_df[study, sp] <- res$estimate
        pval_df[study, sp] <- res$p.value
      } else {
        corr_df[study, sp] <- 0
        pval_df[study, sp] <- 1
      }
    }
  }
  
  # Replace NAs with defaults (corr=0, pval=1)
  corr_df[is.na(corr_df)] <- 0
  pval_df[is.na(pval_df)] <- 1
  
  # Compute direction matrix
  dir <- matrix(0, nrow = nrow(corr_df), ncol = ncol(corr_df),
                dimnames = dimnames(corr_df))
  
  dir[corr_df > 0 & pval_df <= 0.05] <- 2
  dir[corr_df > 0 & pval_df > 0.05 & pval_df <= 0.1] <- 1
  dir[corr_df < 0 & pval_df <= 0.05] <- -2
  dir[corr_df < 0 & pval_df > 0.05 & pval_df <= 0.1] <- -1
  
  # Build association summary
  association_df <- data.frame(
    Positive = apply(dir, 2, function(x) sum(x == 2, na.rm = TRUE)),
    Negative = apply(dir, 2, function(x) sum(x == -2, na.rm = TRUE)),
    Total = nrow(dir)
  )
  
  # Add Score column
  association_df$Score <- apply(association_df, 1, function(x) {
    pos <- x["Positive"]; neg <- x["Negative"]; tot <- x["Total"]
    ((pos - neg) / tot) * (1 - ((min(pos, neg) + 0.00001) / (max(pos, neg) + 0.00001)))
  })
  
  # Store results in list
  age_category_list_shannon[[age_cat]] <- list(
    corr_df = corr_df,
    pval_df = pval_df,
    dir = dir,
    association_df = association_df
  )
}

#View(age_category_list_shannon[["infant"]]$association_df)


## Across Lifestyle Association Scores between Load and different Bifs #####

# Create list to store results for each cohort
cohort_list_shannon <- list()

# Loop through each unique cohort
for (cohort_name in unique(shannon_pielou_Bif_df_16s$cohort)) {
  
  temp_df <- subset(shannon_pielou_Bif_df_16s, cohort == cohort_name)
  bifido_species <- colnames(temp_df)[1:9]
  studies <- unique(temp_df$study_name)
  
  corr_df <- matrix(NA, nrow = length(studies), ncol = length(bifido_species),
                    dimnames = list(studies, bifido_species))
  pval_df <- corr_df
  
  # Loop through each study
  for (study in studies) {
    study_df <- subset(temp_df, study_name == study)
    
    for (sp in bifido_species) {
      if (length(unique(study_df[[sp]])) > 1 && length(unique(study_df$shannon_index)) > 1) {
        res <- suppressWarnings(cor.test(study_df[[sp]], study_df$shannon_index,
                                         method = "spearman", exact = TRUE))
        corr_df[study, sp] <- res$estimate
        pval_df[study, sp] <- res$p.value
      } else {
        corr_df[study, sp] <- 0
        pval_df[study, sp] <- 1
      }
    }
  }
  
  # Replace NAs with defaults (corr=0, pval=1)
  corr_df[is.na(corr_df)] <- 0
  pval_df[is.na(pval_df)] <- 1
  
  # Compute direction matrix
  dir <- matrix(0, nrow = nrow(corr_df), ncol = ncol(corr_df),
                dimnames = dimnames(corr_df))
  
  dir[corr_df > 0 & pval_df <= 0.05] <- 2
  dir[corr_df > 0 & pval_df > 0.05 & pval_df <= 0.1] <- 1
  dir[corr_df < 0 & pval_df <= 0.05] <- -2
  dir[corr_df < 0 & pval_df > 0.05 & pval_df <= 0.1] <- -1
  
  # Build association summary
  association_df <- data.frame(
    Positive = apply(dir, 2, function(x) sum(x == 2, na.rm = TRUE)),
    Negative = apply(dir, 2, function(x) sum(x == -2, na.rm = TRUE)),
    Total = nrow(dir)
  )
  
  # Add Score column
  association_df$Score <- apply(association_df, 1, function(x) {
    pos <- x["Positive"]; neg <- x["Negative"]; tot <- x["Total"]
    ((pos - neg) / tot) * (1 - ((min(pos, neg) + 0.00001) / (max(pos, neg) + 0.00001)))
  })
  
  # Store results in list
  cohort_list_shannon[[cohort_name]] <- list(
    corr_df = corr_df,
    pval_df = pval_df,
    dir = dir,
    association_df = association_df
  )
}

#View(cohort_list_shannon[["IndustrializedUrban"]]$association_df)

############ Combined Association Score df #######

combined_AssociationScores_shannon_df <- data.frame(infant = age_category_list_shannon[["infant"]][["association_df"]]$Score,
                                                    adult = age_category_list_shannon[["adult"]][["association_df"]]$Score,
                                                    senior = age_category_list_shannon[["senior"]][["association_df"]]$Score,
                                                    IndustrializedUrban = cohort_list_shannon[["IndustrializedUrban"]][["association_df"]]$Score,
                                                    Others = cohort_list_shannon[["Others"]][["association_df"]]$Score,
                                                    row.names = rownames(cohort_list_shannon[["Others"]][["association_df"]]))

combined_AssociationScores_shannon_df <- combined_AssociationScores_shannon_df[c(9,1:8),]

# Define desired order
categories <- c("infant", "adult", "senior", 
                "IndustrializedUrban", "Others")
bifnames <- c("detection", "longum", "adolescentis", "pseudocatenulatum",
              "dentium", "bifidum", "breve", "catenulatum", "animalis")

desired_order <- unlist(lapply(categories, function(cat) {
  paste0(cat, "_", bifnames)
}))

# Remove infant_animalis (since it doesn’t exist)
desired_order <- desired_order[!desired_order %in% "infant_animalis"]


df <- combined_AssociationScores_shannon_df
df_long <- as.data.frame(as.table(as.matrix(df)))
df_long$Var1 <- sub("^Bifidobacterium_", "", df_long$Var1)
df_long <- df_long[!(df_long$Var1 == "animalis" & df_long$Var2 == "infant"), ]
df_long$combined <- paste(df_long$Var2, df_long$Var1, sep = "_")
df_association_shannon <- as.data.frame(t(df_long$Freq))
colnames(df_association_shannon) <- df_long$combined
rownames(df_association_shannon) <- "shannon"
df_association_shannon <- df_association_shannon[, intersect(desired_order, colnames(df_association_shannon))]

##################################################################
######## Pielou Association Scores across Age groups ###########

# Create list to store results for each age category
age_category_list_pielou <- list()

# Loop through each unique age category
for (age_cat in unique(shannon_pielou_Bif_df_16s$age_category)) {
  
  temp_df <- subset(shannon_pielou_Bif_df_16s, age_category == age_cat)
  bifido_species <- colnames(temp_df)[1:9]
  studies <- unique(temp_df$study_name)
  
  corr_df <- matrix(NA, nrow = length(studies), ncol = length(bifido_species),
                    dimnames = list(studies, bifido_species))
  pval_df <- corr_df
  
  # Loop through each study
  for (study in studies) {
    study_df <- subset(temp_df, study_name == study)
    
    for (sp in bifido_species) {
      if (length(unique(study_df[[sp]])) > 1 && length(unique(study_df$pielou)) > 1) {
        res <- suppressWarnings(cor.test(study_df[[sp]], study_df$pielou,
                                         method = "spearman", exact = TRUE))
        corr_df[study, sp] <- res$estimate
        pval_df[study, sp] <- res$p.value
      } else {
        corr_df[study, sp] <- 0
        pval_df[study, sp] <- 1
      }
    }
  }
  
  # Replace NAs with defaults (corr=0, pval=1)
  corr_df[is.na(corr_df)] <- 0
  pval_df[is.na(pval_df)] <- 1
  
  # Compute direction matrix
  dir <- matrix(0, nrow = nrow(corr_df), ncol = ncol(corr_df),
                dimnames = dimnames(corr_df))
  
  dir[corr_df > 0 & pval_df <= 0.05] <- 2
  dir[corr_df > 0 & pval_df > 0.05 & pval_df <= 0.1] <- 1
  dir[corr_df < 0 & pval_df <= 0.05] <- -2
  dir[corr_df < 0 & pval_df > 0.05 & pval_df <= 0.1] <- -1
  
  # Build association summary
  association_df <- data.frame(
    Positive = apply(dir, 2, function(x) sum(x == 2, na.rm = TRUE)),
    Negative = apply(dir, 2, function(x) sum(x == -2, na.rm = TRUE)),
    Total = nrow(dir)
  )
  
  # Add Score column
  association_df$Score <- apply(association_df, 1, function(x) {
    pos <- x["Positive"]; neg <- x["Negative"]; tot <- x["Total"]
    ((pos - neg) / tot) * (1 - ((min(pos, neg) + 0.00001) / (max(pos, neg) + 0.00001)))
  })
  
  # Store results in list
  age_category_list_pielou[[age_cat]] <- list(
    corr_df = corr_df,
    pval_df = pval_df,
    dir = dir,
    association_df = association_df
  )
}

#View(age_category_list_pielou[["infant"]]$association_df)


## Across Lifestyle Association Scores between Load and different Bifs #####

# Create list to store results for each cohort
cohort_list_pielou <- list()

# Loop through each unique cohort
for (cohort_name in unique(shannon_pielou_Bif_df_16s$cohort)) {
  
  temp_df <- subset(shannon_pielou_Bif_df_16s, cohort == cohort_name)
  bifido_species <- colnames(temp_df)[1:9]
  studies <- unique(temp_df$study_name)
  
  corr_df <- matrix(NA, nrow = length(studies), ncol = length(bifido_species),
                    dimnames = list(studies, bifido_species))
  pval_df <- corr_df
  
  # Loop through each study
  for (study in studies) {
    study_df <- subset(temp_df, study_name == study)
    
    for (sp in bifido_species) {
      if (length(unique(study_df[[sp]])) > 1 && length(unique(study_df$pielou)) > 1) {
        res <- suppressWarnings(cor.test(study_df[[sp]], study_df$pielou,
                                         method = "spearman", exact = TRUE))
        corr_df[study, sp] <- res$estimate
        pval_df[study, sp] <- res$p.value
      } else {
        corr_df[study, sp] <- 0
        pval_df[study, sp] <- 1
      }
    }
  }
  
  # Replace NAs with defaults (corr=0, pval=1)
  corr_df[is.na(corr_df)] <- 0
  pval_df[is.na(pval_df)] <- 1
  
  # Compute direction matrix
  dir <- matrix(0, nrow = nrow(corr_df), ncol = ncol(corr_df),
                dimnames = dimnames(corr_df))
  
  dir[corr_df > 0 & pval_df <= 0.05] <- 2
  dir[corr_df > 0 & pval_df > 0.05 & pval_df <= 0.1] <- 1
  dir[corr_df < 0 & pval_df <= 0.05] <- -2
  dir[corr_df < 0 & pval_df > 0.05 & pval_df <= 0.1] <- -1
  
  # Build association summary
  association_df <- data.frame(
    Positive = apply(dir, 2, function(x) sum(x == 2, na.rm = TRUE)),
    Negative = apply(dir, 2, function(x) sum(x == -2, na.rm = TRUE)),
    Total = nrow(dir)
  )
  
  # Add Score column
  association_df$Score <- apply(association_df, 1, function(x) {
    pos <- x["Positive"]; neg <- x["Negative"]; tot <- x["Total"]
    ((pos - neg) / tot) * (1 - ((min(pos, neg) + 0.00001) / (max(pos, neg) + 0.00001)))
  })
  
  # Store results in list
  cohort_list_pielou[[cohort_name]] <- list(
    corr_df = corr_df,
    pval_df = pval_df,
    dir = dir,
    association_df = association_df
  )
}

#View(cohort_list_pielou[["IndustrializedUrban"]]$association_df)

############ Combined Association Score df #######

combined_AssociationScores_pielou_df <- data.frame(infant = age_category_list_pielou[["infant"]][["association_df"]]$Score,
                                                   adult = age_category_list_pielou[["adult"]][["association_df"]]$Score,
                                                   senior = age_category_list_pielou[["senior"]][["association_df"]]$Score,
                                                   IndustrializedUrban = cohort_list_pielou[["IndustrializedUrban"]][["association_df"]]$Score,
                                                   Others = cohort_list_pielou[["Others"]][["association_df"]]$Score,
                                                   row.names = rownames(cohort_list_pielou[["Others"]][["association_df"]]))

combined_AssociationScores_pielou_df <- combined_AssociationScores_pielou_df[c(9,1:8),]

df <- combined_AssociationScores_pielou_df
df_long <- as.data.frame(as.table(as.matrix(df)))
df_long$Var1 <- sub("^Bifidobacterium_", "", df_long$Var1)
df_long <- df_long[!(df_long$Var1 == "animalis" & df_long$Var2 == "infant"), ]
df_long$combined <- paste(df_long$Var2, df_long$Var1, sep = "_")
df_association_pielou <- as.data.frame(t(df_long$Freq))
colnames(df_association_pielou) <- df_long$combined
rownames(df_association_pielou) <- "pielou"
df_association_pielou <- df_association_pielou[, intersect(desired_order, colnames(df_association_pielou))]

df_strong_association_all_16s_expanded <- as.data.frame(rbind(carpet_AssociationScores_Heatmap,df_association_shannon,df_association_pielou))

save(df_strong_association_all_16s_expanded,file="/results/StrongAssociationScores_exapanded_16s.RData")

##################################################################









