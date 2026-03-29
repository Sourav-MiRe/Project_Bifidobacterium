  
load("/data/QMP_workspace.RData")
load("/data/bifido_45809_new_2014_not_normalized_Harnandez_added.RData")
load("/data/Shannon_Pielou_Load_df.RData")
load("/data/full_AssociationScoresAgeCategory_wgs_original.RData")
load("/data/full_AssociationScoresCohortType_wgs_original.RData")

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
library(xlsx)


############# Association Scores in WGS #############


grep("_wgs$",ls(),value = T)

SpeciesProfile_All[is.na(SpeciesProfile_All)] <- 0

spdata_45809_wgs <- SpeciesProfile_All[rownames(wgs_SpeciesProfile_filtered_norm),]
shannon_pielou_load_wgs <- shannon_pielou_load_df[rownames(wgs_SpeciesProfile_filtered_norm),]
spdata_45809_wgs$Bifidobacterium_detection <- apply(spdata_45809_wgs[, AllBifidoSpecies], 1, function(x) length(x[x >= 0.0001]))
spdata_45809_wgs$age_category <- shannon_pielou_load_wgs$age_category
spdata_45809_wgs$cohort <- shannon_pielou_load_wgs$cohort
spdata_45809_wgs$study_name <- shannon_pielou_load_wgs$study_name

names(spdata_45809_wgs)[1265] <- "Intestinibacter_bartlettii"
names(spdata_45809_wgs)[1362] <- "Erysipelatoclostridium_ramosum"


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

refined_data_infant <- spdata_45809_wgs[spdata_45809_wgs$age_category == 'infant',]
refined_data_adult <- spdata_45809_wgs[spdata_45809_wgs$age_category == 'adult',]
refined_data_senior <- spdata_45809_wgs[spdata_45809_wgs$age_category == 'senior',]

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


df_association_infant_detection_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_detection_wgs,"infant","Bifidobacterium_detection",Bifidobacterium_detection_infant_studies)
df_association_infant_longum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_longum_wgs,"infant","Bifidobacterium_longum",Bifidobacterium_longum_infant_studies)
df_association_infant_adolescentis_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_adolescentis_wgs,"infant","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_infant_studies)
df_association_infant_pseudocatenulatum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_pseudocatenulatum_wgs,"infant","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_infant_studies)
df_association_infant_dentium_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_dentium_wgs,"infant","Bifidobacterium_dentium",Bifidobacterium_dentium_infant_studies)
df_association_infant_bifidum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_bifidum_wgs,"infant","Bifidobacterium_bifidum",Bifidobacterium_bifidum_infant_studies)
df_association_infant_breve_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_breve_wgs,"infant","Bifidobacterium_breve",Bifidobacterium_breve_infant_studies)
df_association_infant_catenulatum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_catenulatum_wgs,"infant","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_infant_studies)

#### Similarly for Adult and Senior, run the following code

# df_association_adult_detection_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_detection_wgs,"adult","Bifidobacterium_detection",Bifidobacterium_detection_adult_studies)
# df_association_adult_longum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_longum_wgs,"adult","Bifidobacterium_longum",Bifidobacterium_longum_adult_studies)
# df_association_adult_adolescentis_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_adolescentis_wgs,"adult","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_adult_studies)
# df_association_adult_pseudocatenulatum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_pseudocatenulatum_wgs,"adult","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_adult_studies)
# df_association_adult_dentium_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_dentium_wgs,"adult","Bifidobacterium_dentium",Bifidobacterium_dentium_adult_studies)
# df_association_adult_bifidum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_bifidum_wgs,"adult","Bifidobacterium_bifidum",Bifidobacterium_bifidum_adult_studies)
# df_association_adult_breve_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_breve_wgs,"adult","Bifidobacterium_breve",Bifidobacterium_breve_adult_studies)
# df_association_adult_catenulatum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_catenulatum_wgs,"adult","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_adult_studies)
# df_association_adult_animalis_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_animalis_wgs,"adult","Bifidobacterium_animalis",Bifidobacterium_animalis_adult_studies)
# 
# 
# df_association_senior_detection_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_detection_wgs,"senior","Bifidobacterium_detection",Bifidobacterium_detection_senior_studies)
# df_association_senior_longum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_longum_wgs,"senior","Bifidobacterium_longum",Bifidobacterium_longum_senior_studies)
# df_association_senior_adolescentis_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_adolescentis_wgs,"senior","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_senior_studies)
# df_association_senior_pseudocatenulatum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_pseudocatenulatum_wgs,"senior","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_senior_studies)
# df_association_senior_dentium_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_dentium_wgs,"senior","Bifidobacterium_dentium",Bifidobacterium_dentium_senior_studies)
# df_association_senior_bifidum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_bifidum_wgs,"senior","Bifidobacterium_bifidum",Bifidobacterium_bifidum_senior_studies)
# df_association_senior_breve_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_breve_wgs,"senior","Bifidobacterium_breve",Bifidobacterium_breve_senior_studies)
# df_association_senior_catenulatum_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_catenulatum_wgs,"senior","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_senior_studies)
# df_association_senior_animalis_wgs_new <- compute_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_animalis_wgs,"senior","Bifidobacterium_animalis",Bifidobacterium_animalis_senior_studies)
# 
# 
# ###########  Cohort Lifestyle Wise #########################
# Association-Scores stratified by Industrialization status  

# refined_data_IndustrializedUrban <- spdata_45809_wgs[spdata_45809_wgs$age_category == 'adult' & spdata_45809_wgs$cohort == 'IndustrializedUrban',]
# refined_data_UrbanRuralMixed <- spdata_45809_wgs[spdata_45809_wgs$age_category == 'adult' & spdata_45809_wgs$cohort == 'UrbanRuralMixed',]
# refined_data_RuralTribal <- spdata_45809_wgs[spdata_45809_wgs$age_category == 'adult' & spdata_45809_wgs$cohort == 'RuralTribal',]
# 
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
# MajorBifsDetectionUrbanRuralMixed <- compute_detection(refined_data_UrbanRuralMixed,bifido_species,"study_name",unique(refined_data_UrbanRuralMixed$study_name))
# Bifidobacterium_catenulatum_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_catenulatum",]>=0.10)]
# Bifidobacterium_breve_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_breve",]>=0.10)]
# Bifidobacterium_dentium_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_dentium",]>=0.10)]
# Bifidobacterium_bifidum_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_bifidum",]>=0.10)]
# Bifidobacterium_pseudocatenulatum_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_pseudocatenulatum",]>=0.10)]
# Bifidobacterium_adolescentis_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_adolescentis",]>=0.10)]
# Bifidobacterium_longum_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_longum",]>=0.10)]
# Bifidobacterium_animalis_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_animalis",]>=0.10)]
# Bifidobacterium_detection_UrbanRuralMixed_studies <- colnames(MajorBifsDetectionUrbanRuralMixed)[which(MajorBifsDetectionUrbanRuralMixed["Bifidobacterium_detection",]>=0.10)]
# 
# MajorBifsDetectionRuralTribal <- compute_detection(refined_data_RuralTribal,bifido_species,"study_name",unique(refined_data_RuralTribal$study_name))
# Bifidobacterium_catenulatum_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_catenulatum",]>=0.10)]
# Bifidobacterium_breve_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_breve",]>=0.10)]
# Bifidobacterium_dentium_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_dentium",]>=0.10)]
# Bifidobacterium_bifidum_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_bifidum",]>=0.10)]
# Bifidobacterium_pseudocatenulatum_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_pseudocatenulatum",]>=0.10)]
# Bifidobacterium_adolescentis_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_adolescentis",]>=0.10)]
# Bifidobacterium_longum_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_longum",]>=0.10)]
# Bifidobacterium_animalis_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_animalis",]>=0.10)]
# Bifidobacterium_detection_RuralTribal_studies <- colnames(MajorBifsDetectionRuralTribal)[which(MajorBifsDetectionRuralTribal["Bifidobacterium_detection",]>=0.10)]
# 

###################################################################################

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
# df_association_IndustrializedUrban_detection_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs,df_association_IndustrializedUrban_detection_wgs,"IndustrializedUrban","Bifidobacterium_detection",Bifidobacterium_detection_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_longum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_longum_wgs,"IndustrializedUrban","Bifidobacterium_longum",Bifidobacterium_longum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_adolescentis_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_adolescentis_wgs,"IndustrializedUrban","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_pseudocatenulatum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_pseudocatenulatum_wgs,"IndustrializedUrban","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_dentium_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_dentium_wgs,"IndustrializedUrban","Bifidobacterium_dentium",Bifidobacterium_dentium_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_bifidum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_bifidum_wgs,"IndustrializedUrban","Bifidobacterium_bifidum",Bifidobacterium_bifidum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_breve_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_breve_wgs,"IndustrializedUrban","Bifidobacterium_breve",Bifidobacterium_breve_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_catenulatum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_catenulatum_wgs,"IndustrializedUrban","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_animalis_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_animalis_wgs,"IndustrializedUrban","Bifidobacterium_animalis",Bifidobacterium_animalis_IndustrializedUrban_studies)
# 
# df_association_UrbanRuralMixed_detection_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs,df_association_UrbanRuralMixed_detection_wgs,"UrbanRuralMixed","Bifidobacterium_detection",Bifidobacterium_detection_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_longum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_longum_wgs,"UrbanRuralMixed","Bifidobacterium_longum",Bifidobacterium_longum_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_adolescentis_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_adolescentis_wgs,"UrbanRuralMixed","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_pseudocatenulatum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_pseudocatenulatum_wgs,"UrbanRuralMixed","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_dentium_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_dentium_wgs,"UrbanRuralMixed","Bifidobacterium_dentium",Bifidobacterium_dentium_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_bifidum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_bifidum_wgs,"UrbanRuralMixed","Bifidobacterium_bifidum",Bifidobacterium_bifidum_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_breve_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_breve_wgs,"UrbanRuralMixed","Bifidobacterium_breve",Bifidobacterium_breve_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_catenulatum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_catenulatum_wgs,"UrbanRuralMixed","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_animalis_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_animalis_wgs,"UrbanRuralMixed","Bifidobacterium_animalis",Bifidobacterium_animalis_UrbanRuralMixed_studies)
# 
# df_association_RuralTribal_detection_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs,df_association_RuralTribal_detection_wgs,"RuralTribal","Bifidobacterium_detection",Bifidobacterium_detection_RuralTribal_studies)
# df_association_RuralTribal_longum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_longum_wgs,"RuralTribal","Bifidobacterium_longum",Bifidobacterium_longum_RuralTribal_studies)
# df_association_RuralTribal_adolescentis_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_adolescentis_wgs,"RuralTribal","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_RuralTribal_studies)
# df_association_RuralTribal_pseudocatenulatum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_pseudocatenulatum_wgs,"RuralTribal","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_RuralTribal_studies)
# df_association_RuralTribal_dentium_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_dentium_wgs,"RuralTribal","Bifidobacterium_dentium",Bifidobacterium_dentium_RuralTribal_studies)
# df_association_RuralTribal_bifidum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_bifidum_wgs,"RuralTribal","Bifidobacterium_bifidum",Bifidobacterium_bifidum_RuralTribal_studies)
# df_association_RuralTribal_breve_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_breve_wgs,"RuralTribal","Bifidobacterium_breve",Bifidobacterium_breve_RuralTribal_studies)
# df_association_RuralTribal_catenulatum_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_catenulatum_wgs,"RuralTribal","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_RuralTribal_studies)
# df_association_RuralTribal_animalis_wgs_new <- compute_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_animalis_wgs,"RuralTribal","Bifidobacterium_animalis",Bifidobacterium_animalis_RuralTribal_studies)

####################################################################

# overall_lists <- ls(pattern = "_new$")
# save(list = overall_lists, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\Overall_AssociationScores_new.RData")
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
# save(list = list_names, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\df_AssociationScores_new_wgs.RData")

#######################################################################
# ####### Considering Strong Associations (p <= 0.05) ######
# 
# load("G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\Overall_AssociationScores_new.RData")
# 
# obj_names <- ls(pattern = "^df_association_")
# 
# for (obj in obj_names) {
#   
#   # Extract the list
#   lst <- get(obj)
#   
#   # Extract the dir matrix
#   dir_mat <- as.matrix(lst$dir)
#   
#   # Count strong positives (only 2)
#   positive <- apply(dir_mat, 1, function(x) sum(x == 2, na.rm = TRUE))
#   
#   # Count strong negatives (only -2)
#   negative <- apply(dir_mat, 1, function(x) sum(x == -2, na.rm = TRUE))
#   
#   # Total number of studies (columns)
#   total <- ncol(dir_mat)
#   
#   # Association score
#   score <- ((positive - negative) / total) *
#     (1 - ((pmin(positive, negative) + 0.00001) /
#             (pmax(positive, negative) + 0.00001)))
#   
#   # Create strong_association dataframe
#   strong_association <- data.frame(
#     positive = positive,
#     negative = negative,
#     total = total,
#     score = score,
#     row.names = rownames(dir_mat)
#   )
#   
#   # Add to list
#   lst$strong_association <- strong_association
#   
#   # Save back to environment
#   assign(obj, lst)
# }
# 
# adjusted_lists <- ls(pattern = "_new$")
# save(list = adjusted_lists, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\Overall_Strong_AssociationScores_new.RData")

#####################################################################

load("/data/Overall_Strong_AssociationScores_new.RData")

list_names <- ls()[sapply(ls(), function(x) is.list(get(x)))]

for (lname in list_names) {
  lst <- get(lname)
  
  if (!is.null(lst$strong_association)) { 
    assign(lname, lst$strong_association)  
  } else {
    warning(paste("List", lname, "does not contain 'association'"))
  }
}

save(list = list_names, file = "/results/df_Strong_AssociationScores_new_wgs.RData")

####################################################################
### Combined Association Scores df ######

load("/data/df_Strong_AssociationScores_new_wgs.RData")
df_names <- ls(pattern = "^df_association_.*_wgs_new$")
df_list <- lapply(df_names, function(obj) {
  df <- get(obj)                      # get the data frame
  
  # short name for the column: strip prefix and suffix
  short_name <- sub("^df_association_(.*)_wgs_new$", "\\1", obj)
  
  # make a 2-column data frame: species + score
  out <- data.frame(
    species = rownames(df),
    score   = df[,"score"],
    stringsAsFactors = FALSE
  )
  
  # rename 'score' column to the short name
  colnames(out)[2] <- short_name
  
  out
})
combined <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE),
                   df_list)
rownames(combined) <- combined$species
combined$species <- NULL

df_strong_association_all_new_wgs <- combined
df_strong_association_all_new_wgs[is.na(df_strong_association_all_new_wgs)] <- 0

categories <- c("infant", "adult", "senior", 
                "IndustrializedUrban", "UrbanRuralMixed", "RuralTribal")
bifnames <- c("detection", "longum", "adolescentis", "pseudocatenulatum",
              "dentium", "bifidum", "breve", "catenulatum", "animalis")
desired_order <- unlist(lapply(categories, function(cat) {
  paste0(cat, "_", bifnames)
}))
desired_order <- desired_order[desired_order != "infant_animalis"]

df_strong_association_all_new_wgs <- df_strong_association_all_new_wgs[,desired_order]

# Score >= 0.25 in at least 5 columns
keep_rows <- rowSums(abs(df_strong_association_all_new_wgs) >= 0.25, na.rm = TRUE) >= 5

######################################################################
# For each category, get indices of columns that belong to it
# cat_cols <- lapply(categories, function(cat) {
#   grep(paste0("^", cat), colnames(df_strong_association_all_new_wgs))
# })
# 
# names(cat_cols) <- categories
# 
# # Logical matrix: rows = species, columns = categories
# category_hit <- sapply(cat_cols, function(idx) {
#   
#   if (length(idx) == 0) {
#     return(rep(FALSE, nrow(df_strong_association_all_new_wgs)))
#   }
#   
#   # TRUE if any column in that category >= 0.25
#   apply(abs(df_strong_association_all_new_wgs[, idx, drop = FALSE]) >= 0.25,
#         1, any, na.rm = TRUE)
#   
# })
# 
# keep <- rowSums(category_hit) >= 4

###########################################################################

df_strong_association_all_new_wgs_filtered <- df_strong_association_all_new_wgs[keep_rows, ]


#####################################################################
#####################################################################
#####################################################################
# Across Age-Category Association Scores between Load and different Bifs #####

Bif_list <- rev(MajorBifsOrdered)
Bif_SpeciesProf <- wgs_SpeciesProfile_filtered_norm[,Bif_list]
Bif_SpeciesProf$Bifidobacterium_detection <- apply(wgs_SpeciesProfile_filtered_norm[,names(wgs_SpeciesProfile_filtered_norm) %in% AllBifidoSpecies], 1, function(x) length(x[x >= 0.0001]))
wgs_shannon_pielou_load_df <- shannon_pielou_load_df[rownames(Bif_SpeciesProf),]

wgs_shannon_pielou_load_Bif_df <- as.data.frame(cbind(Bif_SpeciesProf, wgs_shannon_pielou_load_df))


# Create list to store results for each age category
age_category_list <- list()

# Loop through each unique age category
for (age_cat in unique(wgs_shannon_pielou_load_Bif_df$age_category)) {
  
  temp_df <- subset(wgs_shannon_pielou_load_Bif_df, age_category == age_cat)
  bifido_species <- colnames(temp_df)[1:9]
  studies <- unique(temp_df$study_name)
  
  corr_df <- matrix(NA, nrow = length(studies), ncol = length(bifido_species),
                    dimnames = list(studies, bifido_species))
  pval_df <- corr_df
  
  # Loop through each study
  for (study in studies) {
    study_df <- subset(temp_df, study_name == study)
    
    for (sp in bifido_species) {
      if (length(unique(study_df[[sp]])) > 1 && length(unique(study_df$load)) > 1) {
        res <- suppressWarnings(cor.test(study_df[[sp]], study_df$load,
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
  age_category_list[[age_cat]] <- list(
    corr_df = corr_df,
    pval_df = pval_df,
    dir = dir,
    association_df = association_df
  )
}


# View(age_category_list[["adult"]]$corr_df)
# View(age_category_list[["adult"]]$pval_df)
# View(age_category_list[["adult"]]$dir)
# View(age_category_list[["adult"]]$association_df)


## Across Lifestyle Association Scores between Load and different Bifs #####

# Create list to store results for each cohort
cohort_list <- list()

# Loop through each unique cohort
for (cohort_name in unique(wgs_shannon_pielou_load_Bif_df$cohort)) {
  
  temp_df <- subset(wgs_shannon_pielou_load_Bif_df, cohort == cohort_name)
  bifido_species <- colnames(temp_df)[1:9]
  studies <- unique(temp_df$study_name)
  
  corr_df <- matrix(NA, nrow = length(studies), ncol = length(bifido_species),
                    dimnames = list(studies, bifido_species))
  pval_df <- corr_df
  
  # Loop through each study
  for (study in studies) {
    study_df <- subset(temp_df, study_name == study)
    
    for (sp in bifido_species) {
      if (length(unique(study_df[[sp]])) > 1 && length(unique(study_df$load)) > 1) {
        res <- suppressWarnings(cor.test(study_df[[sp]], study_df$load,
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
  cohort_list[[cohort_name]] <- list(
    corr_df = corr_df,
    pval_df = pval_df,
    dir = dir,
    association_df = association_df
  )
}


View(cohort_list[["UrbanRuralMixed"]]$corr_df)
View(cohort_list[["UrbanRuralMixed"]]$association_df)

############ Combined Association Score df #######

combined_AssociationScores_df <- data.frame(infant = age_category_list[["infant"]][["association_df"]]$Score,
                                            adult = age_category_list[["adult"]][["association_df"]]$Score,
                                            senior = age_category_list[["senior"]][["association_df"]]$Score,
                                            IndustrializedUrban = cohort_list[["IndustrializedUrban"]][["association_df"]]$Score,
                                            UrbanRuralMixed = cohort_list[["UrbanRuralMixed"]][["association_df"]]$Score,
                                            RuralTribal = cohort_list[["RuralTribal"]][["association_df"]]$Score,
                                            row.names = rownames(cohort_list[["RuralTribal"]][["association_df"]]))

combined_AssociationScores_df <- combined_AssociationScores_df[c(9,1:8),]


df <- combined_AssociationScores_df
df_long <- as.data.frame(as.table(as.matrix(df)))
df_long$Var1 <- sub("^Bifidobacterium_", "", df_long$Var1)
df_long <- df_long[!(df_long$Var1 == "animalis" & df_long$Var2 == "infant"), ]
df_long$combined <- paste(df_long$Var2, df_long$Var1, sep = "_")
df_association_load <- as.data.frame(t(df_long$Freq))
colnames(df_association_load) <- df_long$combined
rownames(df_association_load) <- "load"
desired_order <- c(
  "infant_detection","infant_longum","infant_adolescentis","infant_pseudocatenulatum",
  "infant_dentium","infant_bifidum","infant_breve","infant_catenulatum",
  "adult_detection","adult_longum","adult_adolescentis","adult_pseudocatenulatum",
  "adult_dentium","adult_bifidum","adult_breve","adult_catenulatum","adult_animalis",
  "senior_detection","senior_longum","senior_adolescentis","senior_pseudocatenulatum",
  "senior_dentium","senior_bifidum","senior_breve","senior_catenulatum","senior_animalis",
  "IndustrializedUrban_detection","IndustrializedUrban_longum","IndustrializedUrban_adolescentis",
  "IndustrializedUrban_pseudocatenulatum","IndustrializedUrban_dentium","IndustrializedUrban_bifidum",
  "IndustrializedUrban_breve","IndustrializedUrban_catenulatum","IndustrializedUrban_animalis",
  "UrbanRuralMixed_detection","UrbanRuralMixed_longum","UrbanRuralMixed_adolescentis",
  "UrbanRuralMixed_pseudocatenulatum","UrbanRuralMixed_dentium","UrbanRuralMixed_bifidum",
  "UrbanRuralMixed_breve","UrbanRuralMixed_catenulatum","UrbanRuralMixed_animalis",
  "RuralTribal_detection","RuralTribal_longum","RuralTribal_adolescentis",
  "RuralTribal_pseudocatenulatum","RuralTribal_dentium","RuralTribal_bifidum",
  "RuralTribal_breve","RuralTribal_catenulatum","RuralTribal_animalis"
)

df_association_load <- df_association_load[, intersect(desired_order, colnames(df_association_load))]

##################################################################
######## Shannon Association Scores across Age groups ###########

# Create list to store results for each age category
age_category_list_shannon <- list()

# Loop through each unique age category
for (age_cat in unique(wgs_shannon_pielou_load_Bif_df$age_category)) {
  
  temp_df <- subset(wgs_shannon_pielou_load_Bif_df, age_category == age_cat)
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

View(age_category_list_shannon[["infant"]]$association_df)


## Across Lifestyle Association Scores between Load and different Bifs #####

# Create list to store results for each cohort
cohort_list_shannon <- list()

# Loop through each unique cohort
for (cohort_name in unique(wgs_shannon_pielou_load_Bif_df$cohort)) {
  
  temp_df <- subset(wgs_shannon_pielou_load_Bif_df, cohort == cohort_name)
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

View(cohort_list_shannon[["IndustrializedUrban"]]$association_df)

############ Combined Association Score df #######

combined_AssociationScores_shannon_df <- data.frame(infant = age_category_list_shannon[["infant"]][["association_df"]]$Score,
                                                    adult = age_category_list_shannon[["adult"]][["association_df"]]$Score,
                                                    senior = age_category_list_shannon[["senior"]][["association_df"]]$Score,
                                                    IndustrializedUrban = cohort_list_shannon[["IndustrializedUrban"]][["association_df"]]$Score,
                                                    UrbanRuralMixed = cohort_list_shannon[["UrbanRuralMixed"]][["association_df"]]$Score,
                                                    RuralTribal = cohort_list_shannon[["RuralTribal"]][["association_df"]]$Score,
                                                    row.names = rownames(cohort_list_shannon[["RuralTribal"]][["association_df"]]))

combined_AssociationScores_shannon_df <- combined_AssociationScores_shannon_df[c(9,1:8),]


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
for (age_cat in unique(wgs_shannon_pielou_load_Bif_df$age_category)) {
  
  temp_df <- subset(wgs_shannon_pielou_load_Bif_df, age_category == age_cat)
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

View(age_category_list_pielou[["infant"]]$association_df)


## Across Lifestyle Association Scores between Load and different Bifs #####

# Create list to store results for each cohort
cohort_list_pielou <- list()

# Loop through each unique cohort
for (cohort_name in unique(wgs_shannon_pielou_load_Bif_df$cohort)) {
  
  temp_df <- subset(wgs_shannon_pielou_load_Bif_df, cohort == cohort_name)
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

View(cohort_list_pielou[["IndustrializedUrban"]]$association_df)

############ Combined Association Score df #######

combined_AssociationScores_pielou_df <- data.frame(infant = age_category_list_pielou[["infant"]][["association_df"]]$Score,
                                                   adult = age_category_list_pielou[["adult"]][["association_df"]]$Score,
                                                   senior = age_category_list_pielou[["senior"]][["association_df"]]$Score,
                                                   IndustrializedUrban = cohort_list_pielou[["IndustrializedUrban"]][["association_df"]]$Score,
                                                   UrbanRuralMixed = cohort_list_pielou[["UrbanRuralMixed"]][["association_df"]]$Score,
                                                   RuralTribal = cohort_list_pielou[["RuralTribal"]][["association_df"]]$Score,
                                                   row.names = rownames(cohort_list_pielou[["RuralTribal"]][["association_df"]]))

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

df_strong_association_all_wgs_expanded <- as.data.frame(rbind(df_strong_association_all_new_wgs_filtered,df_association_shannon,df_association_pielou,df_association_load))

#save(df_strong_association_all_wgs_expanded,file="StrongAssociationScores_exapanded_Fig3_Raw.RData")


#####################################################################
#####################################################################
### Strong Association Scores Expanded ######

load("/data/StrongAssociationScores_exapanded_Fig3_Raw.RData")

categories <- c("infant", "adult", "senior", 
                "IndustrializedUrban", "UrbanRuralMixed", "RuralTribal")
bifnames <- c("detection", "longum", "adolescentis", "pseudocatenulatum",
              "bifidum", "dentium", "breve", "catenulatum", "animalis")
desired_order <- unlist(lapply(categories, function(cat) {
  paste0(cat, "_", bifnames)
}))
desired_order <- desired_order[desired_order != "infant_animalis"]

df_strong_association_all_wgs_expanded <- df_strong_association_all_wgs_expanded[,desired_order]


### Heatmap ####

col_fun <- colorRamp2(
  c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
  c("#8B0053", "#CC3377", "#F6B6D2", "#FFFFFF", "#CDE9B6", "#77C679", "#004D00")
)

# Convert to numeric matrix
mat <- as.matrix(df_strong_association_all_wgs_expanded)

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

pdf("/results/StrongAssociationScores_wgs_expanded.pdf", 
    height = 20, width = 16)
draw(ht)
dev.off()

###### CARPET ######

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_AssociationScores_Heatmap <- as.data.frame(carpet_df)

write.xlsx(carpet_AssociationScores_Heatmap, "/results/carpet_StrongAssociationScores_Fig3_expanded.xlsx")

save(carpet_AssociationScores_Heatmap,file="/results/df_strong_association_all_new_wgs_expanded.RData")


#####################################################################
######################################################################

load("/data/dig_carpet_overall.RData")

mat <- as.matrix(dig_carpet_overall)

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

pdf("/results/StrongAssociationScores_wgs_expanded_latest.pdf", 
    height = 20, width = 16)
draw(ht)
dev.off()

##### Carpet #####

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_DiseaseAssociation_Heatmap <- as.data.frame(carpet_df)

write.xlsx(carpet_DiseaseAssociation_Heatmap, "/results/carpet_StrongAssociationScores_Fig3_expanded_latest.xlsx")



