load("/data/Shannon_Pielou_Load_df.RData")
load("/data/all_16s_WGS_combined.RData")
load("/data/bifido_45809_new_2014_not_normalized_Harnandez_added.RData")
load("/data/Mismatched_names_NCBI_taxonomy.RData")
load("/data/SpeciesDetection_AgeCategories.RData")

library(dplyr)
library(randomForest)
library(pcaPP)
library(psych)
library(RColorBrewer)
library(metap)
library(ggplot2)
library(vegan)

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

oob_validation <- function(data,metadata,features_to_check,grouping_column,metadata_column,groups_list)
{
  df_Perform <- as.data.frame(matrix(NA,length(groups_list),2))
  rownames(df_Perform) <- groups_list
  colnames(df_Perform) <- c("oob_corr","oob_p")
  
  pairwise_matrix_r <- as.data.frame(matrix(NA,length(groups_list),length(groups_list)))
  rownames(pairwise_matrix_r) <- groups_list
  colnames(pairwise_matrix_r) <- groups_list
  
  pairwise_matrix_p <- as.data.frame(matrix(NA,length(groups_list),length(groups_list)))
  rownames(pairwise_matrix_p) <- groups_list
  colnames(pairwise_matrix_p) <- groups_list
  
  mat_importance <- as.data.frame(matrix(0,length(groups_list),length(features_to_check)))
  rownames(mat_importance) <- groups_list
  colnames(mat_importance) <- features_to_check
  
  group_list_rows <- rownames(data[data[,grouping_column] %in% groups_list,])
  rows_with_metadata <- rownames(metadata[!is.na(metadata[,metadata_column]),])
  df_pred <- data[intersect(rows_with_metadata,group_list_rows),features_to_check]
  
  for(i in 1:length(groups_list))
  {
    print(i)
    group_name <- groups_list[i]
    train_rows <- intersect(rownames(df_pred),rownames(metadata)[metadata[,grouping_column] == group_name])
    print(train_rows)
    df_train_pred <- df_pred[train_rows,]
    df_train_pred <- apply(df_train_pred,2,function(x)(ifelse(is.nan(x),0,x)))
    df_train_pred <- apply(df_train_pred,2,function(x)(ifelse(is.na(x),0,x)))
    df_train_pred <- df_train_pred[rowSums(df_train_pred)>0,colSums(df_train_pred)>0]
    train_rows <- rownames(df_train_pred)
    
    df_train_response <- metadata[train_rows,metadata_column]
    
    rf_oob_temp <- randomForest(df_train_response~.,df_train_pred)
    
    df_oob_corr <- cor.test(rf_oob_temp$predicted,rf_oob_temp$y,method="spearman",use="pairwise.complete")
    
    print(df_oob_corr)
    
    mat_importance[group_name,rownames(rf_oob_temp$importance)] <- rf_oob_temp$importance
    
    for(j in 1:length(groups_list))
    {
      test_study <- groups_list[j]
      if(i == j)
      {
        pairwise_matrix_r[i,j] <- df_oob_corr$estimate
        pairwise_matrix_p[i,j] <- df_oob_corr$p.value
      }
      else
      {
        test_study_rows <- intersect(rownames(df_pred),rownames(metadata)[metadata[,grouping_column] == test_study])
        test_study_rows <- intersect(rows_with_metadata,test_study_rows)
        df_pred_test_study <- data[test_study_rows,colnames(df_train_pred)]
        response_test_study <- metadata[test_study_rows,metadata_column]
        df_test_study_corr <- cor.test(predict(rf_oob_temp,df_pred_test_study),response_test_study,method="spearman",use="pairwise.complete")
        pairwise_matrix_r[i,j] <- df_test_study_corr$estimate
        pairwise_matrix_p[i,j] <- df_test_study_corr$p.value
      }
    }
    
    df_Perform[i,1] <- df_oob_corr$estimate
    df_Perform[i,2] <- df_oob_corr$p.value
    
    i <- i + 1
  }
  
  df_Correlation <- list("r"=pairwise_matrix_r,"p"=pairwise_matrix_p)
  
  returnlist <- list("Perform"=df_Perform,"Correlation"=df_Correlation,"FeatureImportance"=mat_importance)
  
  return(returnlist)
}

rank_scale=function(x)
{
  x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

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

update_species_names_vec <- function(species_names, df_name) {
  
  # Function to check if a name contains numbers
  contains_numbers <- function(name) grepl("[0-9]", name)
  
  original_names <- species_names  # Keep copy for logging
  
  # Step 1: Apply mismatched_14_24
  name_mapping_14_24 <- mismatched_14_24 %>%
    mutate(names_24 = gsub(" ", "_", names_24)) %>%
    filter(!contains_numbers(names_14) & !contains_numbers(names_24))
  
  name_map_14_24 <- setNames(name_mapping_14_24$names_14, name_mapping_14_24$names_24)
  
  species_names <- ifelse(species_names %in% names(name_map_14_24) & !contains_numbers(species_names), 
                          name_map_14_24[species_names], species_names)
  
  # Step 2: Apply mismatched_18_24
  name_mapping_18_24 <- mismatched_18_24 %>%
    mutate(names_24 = gsub(" ", "_", names_24)) %>%
    filter(!contains_numbers(names_18) & !contains_numbers(names_24))
  
  name_map_18_24 <- setNames(name_mapping_18_24$names_18, name_mapping_18_24$names_24)
  
  changed_indices <- species_names %in% names(name_map_18_24) & !contains_numbers(species_names)
  
  species_names[changed_indices] <- name_map_18_24[species_names[changed_indices]]
  species_names[changed_indices] <- paste0(species_names[changed_indices], "%%%")
  
  # Step 3: Apply mismatched_14_18
  name_mapping_14_18 <- mismatched_14_18 %>%
    mutate(names_18 = gsub(" ", "_", names_18)) %>%
    filter(!contains_numbers(names_14) & !contains_numbers(names_18))
  
  name_map_14_18 <- setNames(name_mapping_14_18$names_14, name_mapping_14_18$names_18)
  
  species_names <- ifelse(species_names %in% names(name_map_14_18) & !contains_numbers(species_names), 
                          name_map_14_18[species_names], species_names)
  
  # Final cleanup: replace spaces, remove brackets
  species_names <- gsub(" ", "_", species_names)
  species_names <- gsub("\\[|\\]", "", species_names)
  
  # Log changes
  changes <- data.frame(
    Original = original_names,
    Updated = species_names
  ) %>%
    filter(Original != Updated)
  
  log_file <- paste0(df_name, "_name_changes_log.txt")
  #write.table(changes, file = log_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\n***********************\n")
  cat("Changes for", df_name, ":\n")
  print(changes)
  cat("\n***********************\n")
  
  return(species_names)
}

SpeciesProfile_All[is.na(SpeciesProfile_All)] <- 0

spdata_46218 <- SpeciesProfile_All
metadata_46218 <- metadata_updated

AllBifidoSpecies <- grep("Bifidobacterium",colnames(spdata_46218),value=TRUE)

spdata_Bifs_rel_46218 <- apply(spdata_46218[,AllBifidoSpecies],2,function(x)(ifelse(is.na(x),0,x)))

spdata_Bifs_detect_46218 <- apply(spdata_Bifs_rel_46218,2,function(x)(ifelse(x>0.0001,1,0)))

studies_46218 <- unique(metadata_46218$study_name)

spdata_46218 <- spdata_46218[intersect(rownames(spdata_46218),rownames(metadata_46218)),]

metadata_46218 <- metadata_46218[intersect(rownames(spdata_46218),rownames(metadata_46218)),]

infant <- rownames(metadata_46218[metadata_46218$age_category == "infant",])

adult <- rownames(metadata_46218[metadata_46218$age_category == "adult",])

senior <- rownames(metadata_46218[metadata_46218$age_category == "senior",])

AllSpecies <- colnames(spdata_46218)

spdata_46218$study_name <- metadata_46218$study_name

infant_studies <- unique(metadata_46218[infant,"study_name"])

adult_studies <- unique(metadata_46218[adult,"study_name"])

senior_studies <- unique(metadata_46218[senior,"study_name"])


spdata_Bifs_detect_46218 <- ifelse(spdata_Bifs_detect_46218 > 0,1,0)
spdata_Bifs_detect_46218 <- as.data.frame(spdata_Bifs_detect_46218)

Bifidobacterium_detection <- rowSums(spdata_Bifs_detect_46218)

SpeciesProfile_All$Bifidobacterium_detection <- Bifidobacterium_detection[rownames(spdata_Bifs_detect_46218)]
SpeciesProfile_All$study_name <- metadata_updated$study_name

names(SpeciesProfile_All)[1265] <- "Intestinibacter_bartlettii"
names(SpeciesProfile_All)[1362] <- "Erysipelatoclostridium_ramosum"

infant_select_species <- names(which(apply(infant_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(infant_species_detection)>=0.33))
adult_select_species <- names(which(apply(adult_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(adult_species_detection)>=0.33))
senior_select_species <- names(which(apply(senior_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(senior_species_detection)>=0.33))


infant_select_non_bif_species <- setdiff(infant_select_species,AllBifidoSpecies)
adult_select_non_bif_species <- setdiff(adult_select_species,AllBifidoSpecies)
senior_select_non_bif_species <- setdiff(senior_select_species,AllBifidoSpecies)

infant_studies_sample_numbers <- table(metadata_46218[infant,"study_name"])
infant_select_studies <- names(which(infant_studies_sample_numbers>=20))

adult_studies_sample_numbers <- table(metadata_46218[adult,"study_name"])
adult_select_studies <- names(which(adult_studies_sample_numbers>=20))

senior_studies_sample_numbers <- table(metadata_46218[senior,"study_name"])
senior_select_studies <- names(which(senior_studies_sample_numbers>=20))

MajorBifidoSpecies <- c("Bifidobacterium_adolescentis","Bifidobacterium_animalis","Bifidobacterium_bifidum","Bifidobacterium_breve","Bifidobacterium_catenulatum","Bifidobacterium_dentium","Bifidobacterium_longum","Bifidobacterium_pseudocatenulatum")

#################################################################
#### Function to run the entire pipeline ####

rf_correlation_pipeline <- function(refined_data,refined_species,species_name,SpeciesProfile_All,shannon_pielou_load_df,adolescentis_infant_combined,oob_validation,rank_scale,groups_list,metadata_updated) {
  message(">>> Normalizing data...")
  refined_data <- as.data.frame(t(apply(refined_data, 1, function(x) x / sum(x))))
  refined_data[is.na(refined_data)] <- 0
  
  message(">>> Adding metadata and feature columns...")
  shannon <- diversity(refined_data, index = "shannon")
  pielou  <- calculate_pielou(refined_data)
  refined_data$shannon <- shannon
  refined_data$pielou <- pielou
  refined_data$load <- shannon_pielou_load_df[rownames(refined_data), "load"]
  refined_data$load[is.na(refined_data$load)] <- 0
  refined_data$study_name <- SpeciesProfile_All[rownames(refined_data), "study_name"]
  refined_data[[species_name]] <- SpeciesProfile_All[rownames(refined_data), species_name]
  
  
  # Add the 3 new features to refined species list
  refined_species <- c(refined_species, "shannon", "pielou", "load")
  
  message(">>> Running Random Forest OOB validation...")
  RF_result <- oob_validation(
    refined_data,
    refined_data,
    refined_species,
    "study_name",
    species_name,
    groups_list
  )
  
  FI_scaled <- as.data.frame(t(apply(RF_result$FeatureImportance, 1, rank_scale)))
  
  message(">>> Running Spearman correlations for ALL refined features...")
  
  study_order <- rownames(RF_result[["Perform"]])
  target_cols <- refined_species
  
  estimate_df <- matrix(NA, nrow = length(study_order), ncol = length(target_cols),
                        dimnames = list(study_order, target_cols))
  pval_df <- estimate_df
  dir_df <- estimate_df
  
  for (study in study_order) {
    study_data <- refined_data[refined_data$study_name == study, ]
    
    if (nrow(study_data) >= 3) {
      for (col in target_cols) {
        
        x <- study_data[[species_name]]
        y <- study_data[[col]]
        
        valid_idx <- complete.cases(x, y)
        
        if (sum(valid_idx) >= 3) {
          test <- suppressWarnings(
            cor.test(x[valid_idx], y[valid_idx], method = "spearman", exact = TRUE)
          )
          
          corr <- test$estimate
          pval <- test$p.value
          
          estimate_df[study, col] <- corr
          pval_df[study, col] <- pval
          
          # Direction coding
          dir_val <- 0
          if (!is.na(corr) && !is.na(pval)) {
            if (corr > 0 && pval <= 0.05) dir_val <- 2
            else if (corr > 0 && pval > 0.05 && pval <= 0.1) dir_val <- 1
            else if (corr < 0 && pval <= 0.05) dir_val <- -2
            else if (corr < 0 && pval > 0.05 && pval <= 0.1) dir_val <- -1
          }
          
          dir_df[study, col] <- dir_val
        }
      }
    }
  }
  
  corr_list <- list(
    estimate = as.data.frame(estimate_df),
    pval = as.data.frame(pval_df),
    dir = as.data.frame(dir_df)
  )
  
  message(">>> Building association summary for ALL refined features...")
  
  association_subset <- data.frame(
    positive = rep(0, length(target_cols)),
    negative = rep(0, length(target_cols)),
    total = rep(0, length(target_cols)),
    row.names = target_cols,
    stringsAsFactors = FALSE
  )
  
  for (col in target_cols) {
    vals <- corr_list$dir[[col]]
    association_subset[col, "positive"] <- sum(vals %in% c(1, 2), na.rm = TRUE)
    association_subset[col, "negative"] <- sum(vals %in% c(-1, -2), na.rm = TRUE)
    association_subset[col, "total"]    <- length(vals)
  }
  
  # Replace NA for load if any
  corr_list$estimate$load[is.na(corr_list$estimate$load)] <- 0
  corr_list$pval$load[is.na(corr_list$pval$load)] <- 1
  
  message(">>> Final assembly of combined output list...")
  
  study_names <- rownames(corr_list$dir)
  
  # Extract corresponding metadata from metadata_updated
  meta_temp <- metadata_updated[match(study_names, metadata_updated$study_name), ]
  
  # Create final metadata with correct rownames and selected columns
  metadata_final <- data.frame(
    seq_type  = meta_temp$seq_type,
    lifestyle = meta_temp$cohort,
    row.names = study_names,
    stringsAsFactors = FALSE)
  
  combined_new <- list(
    metadata = metadata_final,
    dir = corr_list$dir,
    oob = RF_result$Perform,
    pval = corr_list$pval,
    FI = RF_result$FeatureImportance,
    FI_scaled = FI_scaled,
    estimate = corr_list$estimate,
    association = association_subset
  )
  
  message("Pipeline completed for: ", species_name)
  return(combined_new)
}



###################    INFANT    #####################

refined_data <- SpeciesProfile_All[infant,]

MajorBifs <- c("Bifidobacterium_animalis","Bifidobacterium_catenulatum","Bifidobacterium_breve","Bifidobacterium_dentium","Bifidobacterium_bifidum","Bifidobacterium_pseudocatenulatum","Bifidobacterium_adolescentis","Bifidobacterium_longum","Bifidobacterium_detection")

MajorBifsDetectionInfant <- compute_detection(refined_data,MajorBifs,"study_name",infant_select_studies)
Bifidobacterium_animalis_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_animalis",]>=0.10)]
Bifidobacterium_catenulatum_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_catenulatum",]>=0.10)]
Bifidobacterium_breve_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_breve",]>=0.10)]
Bifidobacterium_dentium_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_dentium",]>=0.10)]
Bifidobacterium_bifidum_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_bifidum",]>=0.10)]
Bifidobacterium_pseudocatenulatum_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_pseudocatenulatum",]>=0.10)]
Bifidobacterium_adolescentis_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_adolescentis",]>=0.10)]
Bifidobacterium_longum_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_longum",]>=0.10)]
Bifidobacterium_detection_Infant_studies <- colnames(MajorBifsDetectionInfant)[which(MajorBifsDetectionInfant["Bifidobacterium_detection",]>=0.10)]


refined_species_infant_detection <- names(detection_infant_combined[["FI_scaled"]])
refined_data_infant_detection <- refined_data[infant,refined_species_infant_detection]
detection_infant_combined_new <- rf_correlation_pipeline(refined_data_infant_detection,refined_species_infant_detection,"Bifidobacterium_detection",SpeciesProfile_All,shannon_pielou_load_df,detection_infant_combined,oob_validation,rank_scale,Bifidobacterium_detection_Infant_studies,metadata_updated)

##### To generate for other Bifidobacterium in different age categories you can run the following codes

# refined_species_infant_longum <- names(longum_infant_combined[["FI_scaled"]])
# refined_data_infant_longum <- SpeciesProfile_All[infant,refined_species_infant_longum]
# longum_infant_combined_new <- rf_correlation_pipeline(refined_data_infant_longum,refined_species_infant_longum,"Bifidobacterium_longum",SpeciesProfile_All,shannon_pielou_load_df,longum_infant_combined,oob_validation,rank_scale,Bifidobacterium_longum_Infant_studies,metadata_updated)
# 
# refined_species_infant_adolescentis <- names(adolescentis_infant_combined[["FI_scaled"]])
# refined_data_infant_adolescentis <- SpeciesProfile_All[infant,refined_species_infant_adolescentis]
# adolescentis_infant_combined_new <- rf_correlation_pipeline(refined_data_infant_adolescentis,refined_species_infant_adolescentis,"Bifidobacterium_adolescentis",SpeciesProfile_All,shannon_pielou_load_df,adolescentis_infant_combined,oob_validation,rank_scale,Bifidobacterium_adolescentis_Infant_studies,metadata_updated)
# 
# refined_species_infant_pseudocatenulatum <- names(pseudocatenulatum_infant_combined[["FI_scaled"]])
# refined_data_infant_pseudocatenulatum <- SpeciesProfile_All[infant,refined_species_infant_pseudocatenulatum]
# pseudocatenulatum_infant_combined_new <- rf_correlation_pipeline(refined_data_infant_pseudocatenulatum,refined_species_infant_pseudocatenulatum,"Bifidobacterium_pseudocatenulatum",SpeciesProfile_All,shannon_pielou_load_df,pseudocatenulatum_infant_combined,oob_validation,rank_scale,Bifidobacterium_pseudocatenulatum_Infant_studies,metadata_updated)
# 
# refined_species_infant_bifidum <- names(bifidum_infant_combined[["FI_scaled"]])
# refined_data_infant_bifidum <- SpeciesProfile_All[infant,refined_species_infant_bifidum]
# bifidum_infant_combined_new <- rf_correlation_pipeline(refined_data_infant_bifidum,refined_species_infant_bifidum,"Bifidobacterium_bifidum",SpeciesProfile_All,shannon_pielou_load_df,bifidum_infant_combined,oob_validation,rank_scale,Bifidobacterium_bifidum_Infant_studies,metadata_updated)
# 
# refined_species_infant_dentium <- names(dentium_infant_combined[["FI_scaled"]])
# refined_data_infant_dentium <- SpeciesProfile_All[infant,refined_species_infant_dentium]
# dentium_infant_combined_new <- rf_correlation_pipeline(refined_data_infant_dentium,refined_species_infant_dentium,"Bifidobacterium_dentium",SpeciesProfile_All,shannon_pielou_load_df,dentium_infant_combined,oob_validation,rank_scale,Bifidobacterium_dentium_Infant_studies,metadata_updated)
# 
# refined_species_infant_breve <- names(breve_infant_combined[["FI_scaled"]])
# refined_data_infant_breve <- SpeciesProfile_All[infant,refined_species_infant_breve]
# breve_infant_combined_new <- rf_correlation_pipeline(refined_data_infant_breve,refined_species_infant_breve,"Bifidobacterium_breve",SpeciesProfile_All,shannon_pielou_load_df,breve_infant_combined,oob_validation,rank_scale,Bifidobacterium_breve_Infant_studies,metadata_updated)
# 
# refined_species_infant_catenulatum <- names(catenulatum_infant_combined[["FI_scaled"]])
# refined_data_infant_catenulatum <- SpeciesProfile_All[infant,refined_species_infant_catenulatum]
# catenulatum_infant_combined_new <- rf_correlation_pipeline(refined_data_infant_catenulatum,refined_species_infant_catenulatum,"Bifidobacterium_catenulatum",SpeciesProfile_All,shannon_pielou_load_df,catenulatum_infant_combined,oob_validation,rank_scale,Bifidobacterium_catenulatum_Infant_studies,metadata_updated)
# 
# 
# #################################################################
# ###########          ADULT          ################
# 
# refined_data <- SpeciesProfile_All[adult,]
# 
# MajorBifs <- c("Bifidobacterium_animalis","Bifidobacterium_catenulatum","Bifidobacterium_breve","Bifidobacterium_dentium","Bifidobacterium_bifidum","Bifidobacterium_pseudocatenulatum","Bifidobacterium_adolescentis","Bifidobacterium_longum","Bifidobacterium_detection")
# 
# MajorBifsDetectionAdult <- compute_detection(refined_data,MajorBifs,"study_name",c(adult_select_studies,"Hernandez_2023"))
# Bifidobacterium_animalis_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_animalis",]>=0.10)]
# Bifidobacterium_catenulatum_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_catenulatum",]>=0.10)]
# Bifidobacterium_breve_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_breve",]>=0.10)]
# Bifidobacterium_dentium_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_dentium",]>=0.10)]
# Bifidobacterium_bifidum_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_bifidum",]>=0.10)]
# Bifidobacterium_pseudocatenulatum_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_pseudocatenulatum",]>=0.10)]
# Bifidobacterium_adolescentis_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_adolescentis",]>=0.10)]
# Bifidobacterium_longum_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_longum",]>=0.10)]
# Bifidobacterium_detection_Adult_studies <- colnames(MajorBifsDetectionAdult)[which(MajorBifsDetectionAdult["Bifidobacterium_detection",]>=0.10)]
# 
# 
# rf_correlation_pipeline <- function(refined_data,refined_species,species_name,SpeciesProfile_All,shannon_pielou_load_df,adolescentis_infant_combined,oob_validation,rank_scale,groups_list,metadata_updated) {
#   message(">>> Normalizing data...")
#   refined_data <- as.data.frame(t(apply(refined_data, 1, function(x) x / sum(x))))
#   refined_data[is.na(refined_data)] <- 0
#   
#   message(">>> Adding metadata and feature columns...")
#   shannon <- diversity(refined_data, index = "shannon")
#   pielou  <- calculate_pielou(refined_data)
#   refined_data$shannon <- shannon
#   refined_data$pielou <- pielou
#   refined_data$load <- shannon_pielou_load_df[rownames(refined_data), "load"]
#   refined_data$load[is.na(refined_data$load)] <- 0
#   refined_data$study_name <- SpeciesProfile_All[rownames(refined_data), "study_name"]
#   refined_data[[species_name]] <- SpeciesProfile_All[rownames(refined_data), species_name]
#   
#   
#   # Add the 3 new features to refined species list
#   refined_species <- c(refined_species, "shannon", "pielou", "load")
#   
#   message(">>> Running Random Forest OOB validation...")
#   RF_result <- oob_validation(
#     refined_data,
#     refined_data,
#     refined_species,
#     "study_name",
#     species_name,
#     groups_list
#   )
#   
#   FI_scaled <- as.data.frame(t(apply(RF_result$FeatureImportance, 1, rank_scale)))
#   
#   message(">>> Running Spearman correlations for ALL refined features...")
#   
#   study_order <- rownames(RF_result[["Perform"]])
#   target_cols <- refined_species
#   
#   estimate_df <- matrix(NA, nrow = length(study_order), ncol = length(target_cols),
#                         dimnames = list(study_order, target_cols))
#   pval_df <- estimate_df
#   dir_df <- estimate_df
#   
#   for (study in study_order) {
#     study_data <- refined_data[refined_data$study_name == study, ]
#     
#     if (nrow(study_data) >= 3) {
#       for (col in target_cols) {
#         
#         x <- study_data[[species_name]]
#         y <- study_data[[col]]
#         
#         valid_idx <- complete.cases(x, y)
#         
#         if (sum(valid_idx) >= 3) {
#           test <- suppressWarnings(
#             cor.test(x[valid_idx], y[valid_idx], method = "spearman", exact = TRUE)
#           )
#           
#           corr <- test$estimate
#           pval <- test$p.value
#           
#           estimate_df[study, col] <- corr
#           pval_df[study, col] <- pval
#           
#           # Direction coding
#           dir_val <- 0
#           if (!is.na(corr) && !is.na(pval)) {
#             if (corr > 0 && pval <= 0.05) dir_val <- 2
#             else if (corr > 0 && pval > 0.05 && pval <= 0.1) dir_val <- 1
#             else if (corr < 0 && pval <= 0.05) dir_val <- -2
#             else if (corr < 0 && pval > 0.05 && pval <= 0.1) dir_val <- -1
#           }
#           
#           dir_df[study, col] <- dir_val
#         }
#       }
#     }
#   }
#   
#   corr_list <- list(
#     estimate = as.data.frame(estimate_df),
#     pval = as.data.frame(pval_df),
#     dir = as.data.frame(dir_df)
#   )
#   
#   message(">>> Building association summary for ALL refined features...")
#   
#   association_subset <- data.frame(
#     positive = rep(0, length(target_cols)),
#     negative = rep(0, length(target_cols)),
#     total = rep(0, length(target_cols)),
#     row.names = target_cols,
#     stringsAsFactors = FALSE
#   )
#   
#   for (col in target_cols) {
#     vals <- corr_list$dir[[col]]
#     association_subset[col, "positive"] <- sum(vals %in% c(1, 2), na.rm = TRUE)
#     association_subset[col, "negative"] <- sum(vals %in% c(-1, -2), na.rm = TRUE)
#     association_subset[col, "total"]    <- length(vals)
#   }
#   
#   # Replace NA for load if any
#   corr_list$estimate$load[is.na(corr_list$estimate$load)] <- 0
#   corr_list$pval$load[is.na(corr_list$pval$load)] <- 1
#   
#   message(">>> Final assembly of combined output list...")
#   
#   study_names <- rownames(corr_list$dir)
#   
#   # Extract corresponding metadata from metadata_updated
#   meta_temp <- metadata_updated[match(study_names, metadata_updated$study_name), ]
#   
#   # Create final metadata with correct rownames and selected columns
#   metadata_final <- data.frame(
#     seq_type  = meta_temp$seq_type,
#     lifestyle = meta_temp$cohort,
#     row.names = study_names,
#     stringsAsFactors = FALSE)
#   
#   combined_new <- list(
#     metadata = metadata_final,
#     dir = corr_list$dir,
#     oob = RF_result$Perform,
#     pval = corr_list$pval,
#     FI = RF_result$FeatureImportance,
#     FI_scaled = FI_scaled,
#     estimate = corr_list$estimate,
#     association = association_subset
#   )
#   
#   message("Pipeline completed for: ", species_name)
#   return(combined_new)
# }
# 
# 
# refined_species_adult_detection <- names(detection_adult_combined[["FI_scaled"]])
# refined_data_adult_detection <- refined_data[adult,refined_species_adult_detection]
# detection_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_detection,refined_species_adult_detection,"Bifidobacterium_detection",SpeciesProfile_All,shannon_pielou_load_df,detection_adult_combined,oob_validation,rank_scale,Bifidobacterium_detection_Adult_studies,metadata_updated)
# 
# refined_species_adult_longum <- names(longum_adult_combined[["FI_scaled"]])
# refined_data_adult_longum <- SpeciesProfile_All[adult,refined_species_adult_longum]
# longum_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_longum,refined_species_adult_longum,"Bifidobacterium_longum",SpeciesProfile_All,shannon_pielou_load_df,longum_adult_combined,oob_validation,rank_scale,Bifidobacterium_longum_Adult_studies,metadata_updated)
# 
# refined_species_adult_adolescentis <- names(adolescentis_adult_combined[["FI_scaled"]])
# refined_data_adult_adolescentis <- SpeciesProfile_All[adult,refined_species_adult_adolescentis]
# adolescentis_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_adolescentis,refined_species_adult_adolescentis,"Bifidobacterium_adolescentis",SpeciesProfile_All,shannon_pielou_load_df,adolescentis_adult_combined,oob_validation,rank_scale,Bifidobacterium_adolescentis_Adult_studies,metadata_updated)
# 
# refined_species_adult_pseudocatenulatum <- names(pseudocatenulatum_adult_combined[["FI_scaled"]])
# refined_data_adult_pseudocatenulatum <- SpeciesProfile_All[adult,refined_species_adult_pseudocatenulatum]
# pseudocatenulatum_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_pseudocatenulatum,refined_species_adult_pseudocatenulatum,"Bifidobacterium_pseudocatenulatum",SpeciesProfile_All,shannon_pielou_load_df,pseudocatenulatum_adult_combined,oob_validation,rank_scale,Bifidobacterium_pseudocatenulatum_Adult_studies,metadata_updated)
# 
# refined_species_adult_bifidum <- names(bifidum_adult_combined[["FI_scaled"]])
# refined_data_adult_bifidum <- SpeciesProfile_All[adult,refined_species_adult_bifidum]
# bifidum_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_bifidum,refined_species_adult_bifidum,"Bifidobacterium_bifidum",SpeciesProfile_All,shannon_pielou_load_df,bifidum_adult_combined,oob_validation,rank_scale,Bifidobacterium_bifidum_Adult_studies,metadata_updated)
# 
# refined_species_adult_dentium <- names(dentium_adult_combined[["FI_scaled"]])
# refined_data_adult_dentium <- SpeciesProfile_All[adult,refined_species_adult_dentium]
# dentium_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_dentium,refined_species_adult_dentium,"Bifidobacterium_dentium",SpeciesProfile_All,shannon_pielou_load_df,dentium_adult_combined,oob_validation,rank_scale,Bifidobacterium_dentium_Adult_studies,metadata_updated)
# 
# refined_species_adult_breve <- names(breve_adult_combined[["FI_scaled"]])
# refined_data_adult_breve <- SpeciesProfile_All[adult,refined_species_adult_breve]
# breve_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_breve,refined_species_adult_breve,"Bifidobacterium_breve",SpeciesProfile_All,shannon_pielou_load_df,breve_adult_combined,oob_validation,rank_scale,Bifidobacterium_breve_Adult_studies,metadata_updated)
# 
# refined_species_adult_catenulatum <- names(catenulatum_adult_combined[["FI_scaled"]])
# refined_data_adult_catenulatum <- SpeciesProfile_All[adult,refined_species_adult_catenulatum]
# catenulatum_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_catenulatum,refined_species_adult_catenulatum,"Bifidobacterium_catenulatum",SpeciesProfile_All,shannon_pielou_load_df,catenulatum_adult_combined,oob_validation,rank_scale,Bifidobacterium_catenulatum_Adult_studies,metadata_updated)
# 
# refined_species_adult_animalis <- names(animalis_adult_combined[["FI_scaled"]])
# refined_data_adult_animalis <- SpeciesProfile_All[adult,refined_species_adult_animalis]
# animalis_adult_combined_new <- rf_correlation_pipeline(refined_data_adult_animalis,refined_species_adult_animalis,"Bifidobacterium_animalis",SpeciesProfile_All,shannon_pielou_load_df,animalis_adult_combined,oob_validation,rank_scale,Bifidobacterium_animalis_Adult_studies,metadata_updated)
# 
# 
# 
# 
# 
# #################################################################
# ###########          SENIOR          ################
# 
# 
# refined_data <- SpeciesProfile_All[senior,]
# 
# MajorBifs <- c("Bifidobacterium_animalis","Bifidobacterium_catenulatum","Bifidobacterium_breve","Bifidobacterium_dentium","Bifidobacterium_bifidum","Bifidobacterium_pseudocatenulatum","Bifidobacterium_adolescentis","Bifidobacterium_longum","Bifidobacterium_detection")
# 
# MajorBifsDetectionSenior <- compute_detection(refined_data,MajorBifs,"study_name",c(senior_select_studies,"Hernandez_2023"))
# Bifidobacterium_animalis_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_animalis",]>=0.10)]
# Bifidobacterium_catenulatum_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_catenulatum",]>=0.10)]
# Bifidobacterium_breve_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_breve",]>=0.10)]
# Bifidobacterium_dentium_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_dentium",]>=0.10)]
# Bifidobacterium_bifidum_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_bifidum",]>=0.10)]
# Bifidobacterium_pseudocatenulatum_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_pseudocatenulatum",]>=0.10)]
# Bifidobacterium_adolescentis_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_adolescentis",]>=0.10)]
# Bifidobacterium_longum_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_longum",]>=0.10)]
# Bifidobacterium_detection_Senior_studies <- colnames(MajorBifsDetectionSenior)[which(MajorBifsDetectionSenior["Bifidobacterium_detection",]>=0.10)]
# 
# 
# rf_correlation_pipeline <- function(refined_data,refined_species,species_name,SpeciesProfile_All,shannon_pielou_load_df,adolescentis_infant_combined,oob_validation,rank_scale,groups_list,metadata_updated) {
#   message(">>> Normalizing data...")
#   refined_data <- as.data.frame(t(apply(refined_data, 1, function(x) x / sum(x))))
#   refined_data[is.na(refined_data)] <- 0
#   
#   message(">>> Adding metadata and feature columns...")
#   shannon <- diversity(refined_data, index = "shannon")
#   pielou  <- calculate_pielou(refined_data)
#   refined_data$shannon <- shannon
#   refined_data$pielou <- pielou
#   refined_data$load <- shannon_pielou_load_df[rownames(refined_data), "load"]
#   refined_data$load[is.na(refined_data$load)] <- 0
#   refined_data$study_name <- SpeciesProfile_All[rownames(refined_data), "study_name"]
#   refined_data[[species_name]] <- SpeciesProfile_All[rownames(refined_data), species_name]
#   
#   
#   # Add the 3 new features to refined species list
#   refined_species <- c(refined_species, "shannon", "pielou", "load")
# 
#   message(">>> Running Random Forest OOB validation...")
#   RF_result <- oob_validation(
#     refined_data,
#     refined_data,
#     refined_species,
#     "study_name",
#     species_name,
#     groups_list
#   )
#   
#   FI_scaled <- as.data.frame(t(apply(RF_result$FeatureImportance, 1, rank_scale)))
#   
#   message(">>> Running Spearman correlations for ALL refined features...")
#   
#   study_order <- rownames(RF_result[["Perform"]])
#   target_cols <- refined_species
#   
#   estimate_df <- matrix(NA, nrow = length(study_order), ncol = length(target_cols),
#                         dimnames = list(study_order, target_cols))
#   pval_df <- estimate_df
#   dir_df <- estimate_df
#   
#   for (study in study_order) {
#     study_data <- refined_data[refined_data$study_name == study, ]
#     
#     if (nrow(study_data) >= 3) {
#       for (col in target_cols) {
#         
#         x <- study_data[[species_name]]
#         y <- study_data[[col]]
#         
#         valid_idx <- complete.cases(x, y)
#         
#         if (sum(valid_idx) >= 3) {
#           test <- suppressWarnings(
#             cor.test(x[valid_idx], y[valid_idx], method = "spearman", exact = TRUE)
#           )
#           
#           corr <- test$estimate
#           pval <- test$p.value
#           
#           estimate_df[study, col] <- corr
#           pval_df[study, col] <- pval
#           
#           # Direction coding
#           dir_val <- 0
#           if (!is.na(corr) && !is.na(pval)) {
#             if (corr > 0 && pval <= 0.05) dir_val <- 2
#             else if (corr > 0 && pval > 0.05 && pval <= 0.1) dir_val <- 1
#             else if (corr < 0 && pval <= 0.05) dir_val <- -2
#             else if (corr < 0 && pval > 0.05 && pval <= 0.1) dir_val <- -1
#           }
#           
#           dir_df[study, col] <- dir_val
#         }
#       }
#     }
#   }
#   
#   corr_list <- list(
#     estimate = as.data.frame(estimate_df),
#     pval = as.data.frame(pval_df),
#     dir = as.data.frame(dir_df)
#   )
#   
#   message(">>> Building association summary for ALL refined features...")
#   
#   association_subset <- data.frame(
#     positive = rep(0, length(target_cols)),
#     negative = rep(0, length(target_cols)),
#     total = rep(0, length(target_cols)),
#     row.names = target_cols,
#     stringsAsFactors = FALSE
#   )
#   
#   for (col in target_cols) {
#     vals <- corr_list$dir[[col]]
#     association_subset[col, "positive"] <- sum(vals %in% c(1, 2), na.rm = TRUE)
#     association_subset[col, "negative"] <- sum(vals %in% c(-1, -2), na.rm = TRUE)
#     association_subset[col, "total"]    <- length(vals)
#   }
#   
#   # Replace NA for load if any
#   corr_list$estimate$load[is.na(corr_list$estimate$load)] <- 0
#   corr_list$pval$load[is.na(corr_list$pval$load)] <- 1
#   
#   message(">>> Final assembly of combined output list...")
#   
#   study_names <- rownames(corr_list$dir)
#   
#   # Extract corresponding metadata from metadata_updated
#   meta_temp <- metadata_updated[match(study_names, metadata_updated$study_name), ]
#   
#   # Create final metadata with correct rownames and selected columns
#   metadata_final <- data.frame(
#     seq_type  = meta_temp$seq_type,
#     lifestyle = meta_temp$cohort,
#     row.names = study_names,
#     stringsAsFactors = FALSE)
#   
#   combined_new <- list(
#     metadata = metadata_final,
#     dir = corr_list$dir,
#     oob = RF_result$Perform,
#     pval = corr_list$pval,
#     FI = RF_result$FeatureImportance,
#     FI_scaled = FI_scaled,
#     estimate = corr_list$estimate,
#     association = association_subset
#   )
#   
#   message("Pipeline completed for: ", species_name)
#   return(combined_new)
# }
# 
# 
# refined_species_senior_detection <- names(detection_senior_combined[["FI_scaled"]])
# refined_data_senior_detection <- refined_data[,refined_species_senior_detection]
# detection_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_detection,refined_species_senior_detection,"Bifidobacterium_detection",SpeciesProfile_All,shannon_pielou_load_df,detection_senior_combined,oob_validation,rank_scale,Bifidobacterium_detection_Senior_studies,metadata_updated)
# 
# refined_species_senior_longum <- names(longum_senior_combined[["FI_scaled"]])
# refined_data_senior_longum <- SpeciesProfile_All[senior,refined_species_senior_longum]
# longum_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_longum,refined_species_senior_longum,"Bifidobacterium_longum",SpeciesProfile_All,shannon_pielou_load_df,longum_senior_combined,oob_validation,rank_scale,Bifidobacterium_longum_Senior_studies,metadata_updated)
# 
# refined_species_senior_adolescentis <- names(adolescentis_senior_combined[["FI_scaled"]])
# refined_data_senior_adolescentis <- SpeciesProfile_All[senior,refined_species_senior_adolescentis]
# adolescentis_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_adolescentis,refined_species_senior_adolescentis,"Bifidobacterium_adolescentis",SpeciesProfile_All,shannon_pielou_load_df,adolescentis_senior_combined,oob_validation,rank_scale,Bifidobacterium_adolescentis_Senior_studies,metadata_updated)
# 
# refined_species_senior_pseudocatenulatum <- names(pseudocatenulatum_senior_combined[["FI_scaled"]])
# refined_data_senior_pseudocatenulatum <- SpeciesProfile_All[senior,refined_species_senior_pseudocatenulatum]
# pseudocatenulatum_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_pseudocatenulatum,refined_species_senior_pseudocatenulatum,"Bifidobacterium_pseudocatenulatum",SpeciesProfile_All,shannon_pielou_load_df,pseudocatenulatum_senior_combined,oob_validation,rank_scale,Bifidobacterium_pseudocatenulatum_Senior_studies,metadata_updated)
# 
# refined_species_senior_bifidum <- names(bifidum_senior_combined[["FI_scaled"]])
# refined_data_senior_bifidum <- SpeciesProfile_All[senior,refined_species_senior_bifidum]
# bifidum_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_bifidum,refined_species_senior_bifidum,"Bifidobacterium_bifidum",SpeciesProfile_All,shannon_pielou_load_df,bifidum_senior_combined, oob_validation,rank_scale,Bifidobacterium_bifidum_Senior_studies,metadata_updated)
# 
# refined_species_senior_dentium <- names(dentium_senior_combined[["FI_scaled"]])
# refined_data_senior_dentium <- SpeciesProfile_All[senior,refined_species_senior_dentium]
# dentium_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_dentium,refined_species_senior_dentium,"Bifidobacterium_dentium",SpeciesProfile_All,shannon_pielou_load_df,dentium_senior_combined,oob_validation,rank_scale,Bifidobacterium_dentium_Senior_studies,metadata_updated)
# 
# refined_species_senior_breve <- names(breve_senior_combined[["FI_scaled"]])
# refined_data_senior_breve <- SpeciesProfile_All[senior,refined_species_senior_breve]
# breve_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_breve,refined_species_senior_breve,"Bifidobacterium_breve",SpeciesProfile_All,shannon_pielou_load_df,breve_senior_combined, oob_validation,rank_scale,Bifidobacterium_breve_Senior_studies,metadata_updated)
# 
# refined_species_senior_catenulatum <- names(catenulatum_senior_combined[["FI_scaled"]])
# refined_data_senior_catenulatum <- SpeciesProfile_All[senior,refined_species_senior_catenulatum]
# catenulatum_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_catenulatum,refined_species_senior_catenulatum,"Bifidobacterium_catenulatum",SpeciesProfile_All,shannon_pielou_load_df,catenulatum_senior_combined,oob_validation,rank_scale,Bifidobacterium_catenulatum_Senior_studies,metadata_updated)
# 
# refined_species_senior_animalis <- names(animalis_senior_combined[["FI_scaled"]])
# refined_data_senior_animalis <- SpeciesProfile_All[senior,refined_species_senior_animalis]
# animalis_senior_combined_new <- rf_correlation_pipeline(refined_data_senior_animalis,refined_species_senior_animalis,"Bifidobacterium_animalis",SpeciesProfile_All,shannon_pielou_load_df,animalis_senior_combined,oob_validation,rank_scale,Bifidobacterium_animalis_Senior_studies,metadata_updated)
# 
# 
# save(detection_infant_combined_new,detection_adult_combined_new,detection_senior_combined_new,adolescentis_infant_combined_new, adolescentis_adult_combined_new,adolescentis_senior_combined_new,animalis_adult_combined_new,animalis_senior_combined_new,bifidum_infant_combined_new,bifidum_adult_combined_new,bifidum_senior_combined_new,catenulatum_infant_combined_new,catenulatum_adult_combined_new,catenulatum_senior_combined_new,dentium_infant_combined_new,dentium_adult_combined_new,dentium_senior_combined_new,breve_infant_combined_new,breve_adult_combined_new,breve_senior_combined_new, longum_infant_combined_new,longum_adult_combined_new,longum_senior_combined_new,pseudocatenulatum_infant_combined_new,pseudocatenulatum_adult_combined_new,pseudocatenulatum_senior_combined_new, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\Fig2_ReComputation\\all_16s_WGS_combined_NEW.RData")
# 
# 







