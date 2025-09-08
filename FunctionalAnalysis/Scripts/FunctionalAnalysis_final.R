library(dplyr)
library(ggplot2)
library(randomForest)
library(Metrics)
library(xlsx)
library(pheatmap)
library(cluster)      
library(factoextra)   
library(ggplot2)
library(ggrepel)
library(ade4)
library(vegan)
library(stats)
library(pcaPP)
library(ComplexHeatmap)
library(circlize)

load("G:\\My Drive\\Bifido_Project\\df_association_CommonSpecies_SeniorAdultMean.RData")
load("G:\\My Drive\\Bifido_Project\\mean_FunctionalProfileAll_groups.RData")

####################################################################
# Random Forest on Mean Functional Profile
#####################################################################

rf_models <- list()

correlation_results <- data.frame(
  Variable = character(),
  Spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)


for (var_name in colnames(df_association_CommonSpecies_modified)) {
  
  response <- df_association_CommonSpecies_modified[[var_name]]
  
  rf_model <- randomForest(x = mean_FunctionaProfile_groups, y = response)
  
  rf_models[[var_name]] <- rf_model
  
  predicted <- rf_model$predicted
  
  cor_test <- cor.test(predicted, response, method = "spearman")
  
  correlation_results <- rbind(correlation_results, data.frame(
    Variable = var_name,
    Spearman_rho = cor_test$estimate,
    p_value = cor_test$p.value,
    stringsAsFactors = FALSE
  ))
}

rm(rf_model)

rownames(correlation_results) <- correlation_results$Variable
correlation_results$Variable <- NULL

rownames(correlation_results) <- gsub("Bifidobacterium_","",rownames(correlation_results))
rownames(correlation_results) <- paste0(rownames(correlation_results),"_score")

write.xlsx(correlation_results,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/correlation_results_MeanFunctionalProfile.xlsx")


##### Extracting Important features (90% cumulative FI) ######

top_features_rf_MeanFunctionalProfile <- list()

for (var_name in names(rf_models)) {
  
  importance_vals <- importance(rf_models[[var_name]])
  
  importance_df <- data.frame(
    Feature = rownames(importance_vals),
    Importance = importance_vals[, 1],
    stringsAsFactors = FALSE
  )
  
  importance_df <- importance_df[order(-importance_df$Importance), ]
  
  importance_df$Cumulative <- cumsum(importance_df$Importance) / sum(importance_df$Importance)
  
  top_features <- importance_df[importance_df$Cumulative <= 0.9, c("Feature", "Importance")]
  
  top_features_rf_MeanFunctionalProfile[[var_name]] <- top_features
}


#all(diff(top_features_rf_MeanFunctionalProfile[["Bifidobacterium_longum"]]$Importance) <= 0)


for (var_name in names(top_features_rf_MeanFunctionalProfile)) {
  
  top_features <- top_features_rf_MeanFunctionalProfile[[var_name]]
  
  # Merge with annotation metadata
  annotated <- merge(
    top_features,
    annotation_df_mean_FunctionalProfile,
    by.x = "Feature",
    by.y = "Groups",
    all.x = TRUE
  )
  
  annotated <- annotated[, c("Feature", "Importance", "Total_features", "Feature_ids", "Annotations")]
  
  annotated <- annotated[order(-annotated$Importance), ]
  
  top_features_rf_MeanFunctionalProfile[[var_name]] <- annotated
}

##########################################################################
### Filtering features based on certain criteria ###
##########################################################################

top_features_rf_MeanFunctionalProfile$detection$bif_score <- "detection_score"
top_features_rf_MeanFunctionalProfile$detection$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$detection)

top_features_rf_MeanFunctionalProfile$Bifidobacterium_longum$bif_score <- "longum_score"
top_features_rf_MeanFunctionalProfile$Bifidobacterium_longum$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$Bifidobacterium_longum)

top_features_rf_MeanFunctionalProfile$Bifidobacterium_adolescentis$bif_score <- "adolescentis_score"
top_features_rf_MeanFunctionalProfile$Bifidobacterium_adolescentis$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$Bifidobacterium_adolescentis)

top_features_rf_MeanFunctionalProfile$Bifidobacterium_pseudocatenulatum$bif_score <- "pseudocatenulatum_score"
top_features_rf_MeanFunctionalProfile$Bifidobacterium_pseudocatenulatum$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$Bifidobacterium_pseudocatenulatum)

top_features_rf_MeanFunctionalProfile$Bifidobacterium_dentium$bif_score <- "dentium_score"
top_features_rf_MeanFunctionalProfile$Bifidobacterium_dentium$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$Bifidobacterium_dentium)

top_features_rf_MeanFunctionalProfile$Bifidobacterium_bifidum$bif_score <- "bifidum_score"
top_features_rf_MeanFunctionalProfile$Bifidobacterium_bifidum$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$Bifidobacterium_bifidum)

top_features_rf_MeanFunctionalProfile$Bifidobacterium_breve$bif_score <- "breve_score"
top_features_rf_MeanFunctionalProfile$Bifidobacterium_breve$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$Bifidobacterium_breve)

top_features_rf_MeanFunctionalProfile$Bifidobacterium_catenulatum$bif_score <- "catenulatum_score"
top_features_rf_MeanFunctionalProfile$Bifidobacterium_catenulatum$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$Bifidobacterium_catenulatum)

top_features_rf_MeanFunctionalProfile$Bifidobacterium_animalis$bif_score <- "animalis_score"
top_features_rf_MeanFunctionalProfile$Bifidobacterium_animalis$rank <- 1:nrow(top_features_rf_MeanFunctionalProfile$Bifidobacterium_animalis)


combined_df_FI_all <- as.data.frame(rbind(top_features_rf_MeanFunctionalProfile$detection,
                                          top_features_rf_MeanFunctionalProfile$Bifidobacterium_longum,
                                          top_features_rf_MeanFunctionalProfile$Bifidobacterium_adolescentis,
                                          top_features_rf_MeanFunctionalProfile$Bifidobacterium_pseudocatenulatum,
                                          top_features_rf_MeanFunctionalProfile$Bifidobacterium_dentium,
                                          top_features_rf_MeanFunctionalProfile$Bifidobacterium_bifidum,
                                          top_features_rf_MeanFunctionalProfile$Bifidobacterium_breve,
                                          top_features_rf_MeanFunctionalProfile$Bifidobacterium_catenulatum,
                                          top_features_rf_MeanFunctionalProfile$Bifidobacterium_animalis))
all_feaures_detection <- data.frame(
  Group = colnames(mean_FunctionaProfile_groups),
  Detection = colSums(mean_FunctionaProfile_groups > 0))

combined_df_FI_all$Detection <- all_feaures_detection$Detection[match(combined_df_FI_all$Feature, all_feaures_detection$Group)]

####### Applying filtration criteria #####

library(dplyr)

feature_category_counts <- combined_df_FI_all %>%
  group_by(Feature) %>%
  summarise(
    n_categories = n_distinct(bif_score),
    has_rank_1_to_50 = any(rank >= 1 & rank <= 50),
    has_detection_2_or_more = any(Detection >= 2),
    .groups = "drop"
  ) %>%
  filter(n_categories > 1, has_rank_1_to_50, has_detection_2_or_more)

combined_df_FI_filtered <- combined_df_FI_all %>%
  filter(Feature %in% feature_category_counts$Feature)

names(combined_df_FI_filtered)[8] <- "Detection_NonBif"

sum(table(combined_df_FI_filtered$Feature) > 3)
names(table(combined_df_FI_filtered$Feature) > 3)[table(combined_df_FI_filtered$Feature) > 3]
combined_df_FI_filtered[combined_df_FI_filtered$Feature == 'Group10788',]

combined_df_FI_filtered <- combined_df_FI_filtered[,c(1,2,6,7,8,3,4,5)]

model_wise_important_features <- data.frame(model = names(table(combined_df_FI_filtered$bif_score)), no_of_features = as.numeric(table(combined_df_FI_filtered$bif_score)))

write.xlsx(model_wise_important_features,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/model_wise_important_features.xlsx")
write.xlsx(combined_df_FI_filtered,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/combined_df_FI_filtered.xlsx")

######################################################################
## Correlation with Important features 
#######################################################################

# Load required packages
library(dplyr)
library(readxl)
library(openxlsx)

#Data loaded 
load("mean_FunctionaProfile_groups.RData")
#load("df_association_CommonSpecies_modified.RData")
annotation <- read_excel("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/combined_df_FI_filtered.xlsx")

# Step 1: Collapse duplicated Feature (group) entries
annotation_unique <- annotation %>%
  group_by(Feature) %>%
  summarise(
    Importance = paste(Importance, collapse = ", "),
    Bif_score = paste(Bif_score, collapse = ", "),
    Rank = paste(Rank, collapse = ", "),
    Detection_NonBif = first(Detection_NonBif),  # All values are the same per group
    Total_features = first(Total_features),
    Feature_ids = first(Feature_ids),
    Annotations = first(Annotations)
  ) %>%
  ungroup()

# Step 2: Create Detection_NonBif_list using mean_FunctionaProfile_groups
# mean_FunctionaProfile_groups has species as rownames and group names as columns

# Transpose for easier access (species in columns)
t_profile <- as.data.frame(mean_FunctionaProfile_groups)
t_profile$Species <- rownames(mean_FunctionaProfile_groups)

# Function to get species names with >0 value for a group
get_nonbif_species <- function(group_name) {
  if (!group_name %in% colnames(mean_FunctionaProfile_groups)) {
    return(NA)
  }
  species_values <- mean_FunctionaProfile_groups[, group_name]
  species_names <- rownames(mean_FunctionaProfile_groups)[species_values > 0]
  return(paste(species_names, collapse = ", "))
}

# Apply to each group in annotation_unique
annotation_unique$Detection_NonBif_list <- sapply(annotation_unique$Feature, get_nonbif_species)

write.xlsx(annotation_unique, file = "G:/My Drive/Bifido_Project/New_FunctionalAnalysis/combined_df_FI_filtered_unique.xlsx", rowNames = FALSE)

#_____________________________________________________________________________________________________________________________

#cor.test(mean_FunctionaProfile_groups$Group1, df_association_CommonSpecies_modified$Bifidobacterium_pseudocatenulatum, method = 'spearman', exact = T)

bif_list <- c("detection","longum", "adolescentis","pseudocatenulatum", "bifidum", "breve", "catenulatum", "dentium", "animalis")
names(df_association_CommonSpecies_modified)[9] <- "Bifidobacterium_detection"

unique_groups_bif_score_corr <- data.frame(matrix(NA, nrow = nrow(annotation_unique), ncol = length(bif_list) * 2))
rownames(unique_groups_bif_score_corr) <- annotation_unique$Feature
colnames(unique_groups_bif_score_corr) <- as.vector(sapply(bif_list, function(b) c(paste0(b, "_score_corr"), paste0(b, "_score_pval"))))

for (i in seq_len(nrow(annotation_unique))) {
  
  group <- annotation_unique$Feature[i]
  bif_scores_str <- annotation_unique$Bif_score[i]
  
  if (!group %in% colnames(mean_FunctionaProfile_groups)) next  # Skip if group column not present
  
  functional_values <- mean_FunctionaProfile_groups[, group]
  
  # Extract bif species list (remove _score, split by comma, trim whitespace)
  bif_species <- trimws(unlist(strsplit(gsub("_score", "", bif_scores_str), ",")))
  
  for (bif in bif_species) {
    bif_colname <- paste0("Bifidobacterium_", bif)
    if (!bif_colname %in% colnames(df_association_CommonSpecies_modified)) next 
    bif_scores <- df_association_CommonSpecies_modified[[bif_colname]]
    
    test_result <- suppressWarnings(cor.test(functional_values, bif_scores, method = "spearman", exact = TRUE))
    unique_groups_bif_score_corr[group, paste0(bif, "_score_corr")] <- test_result$estimate
    unique_groups_bif_score_corr[group, paste0(bif, "_score_pval")] <- test_result$p.value
  }
}

#check for sanity 
cor.test(mean_FunctionaProfile_groups$Group10405, df_association_CommonSpecies_modified$Bifidobacterium_longum, method = 'spearman', exact = T)

#write.xlsx(unique_groups_bif_score_corr, file = "unique_groups_bif_score_corr.xlsx", rowNames = TRUE)

#directinality
uniGroup_BifScore_corr_direct <- data.frame(matrix(0, nrow = nrow(unique_groups_bif_score_corr), ncol = length(bif_list)))
rownames(uniGroup_BifScore_corr_direct) <- rownames(unique_groups_bif_score_corr)
colnames(uniGroup_BifScore_corr_direct) <- paste0(bif_list, "_corr")

for (bif in bif_list) {
  corr_col <- paste0(bif, "_score_corr")
  pval_col <- paste0(bif, "_score_pval")
  
  corr <- unique_groups_bif_score_corr[[corr_col]]
  pval <- unique_groups_bif_score_corr[[pval_col]]
  
  score <- ifelse(is.na(corr) | is.na(pval), 0,
                  ifelse(corr < 0 & pval < 0.05, -1,
                         ifelse(corr > 0 & pval < 0.05, 1, 0)))
  
  uniGroup_BifScore_corr_direct[[paste0(bif, "_corr")]] <- score
}


#for (bif in c("longum", "adolescentis", "animals", "breve", "bifidum", "catenulatum", "pseudocatenulatum", "dentium")) print(unique_groups_bif_score_corr[unique_groups_bif_score_corr[[paste0(bif, "_score_pval")]] > 0.05, c(paste0(bif, "_score_corr"), paste0(bif, "_score_pval"))])

#Group10788 has longum non-significant

save(unique_groups_bif_score_corr, uniGroup_BifScore_corr_direct, file = "G:/My Drive/Bifido_Project/New_FunctionalAnalysis/features_bifs_correlation.RData")

#########################################################################
#### Correlation heatmap 
#########################################################################

corr_df <- unique_groups_bif_score_corr[, grep("_score_corr$", colnames(unique_groups_bif_score_corr))]
pval_df <- unique_groups_bif_score_corr[, grep("_score_pval$", colnames(unique_groups_bif_score_corr))]
species <- unique(gsub("_score_.*", "", colnames(unique_groups_bif_score_corr)))

direction_matrix <- matrix(0,nrow = nrow(unique_groups_bif_score_corr), ncol = 9)
rownames(direction_matrix) <- rownames(unique_groups_bif_score_corr)
colnames(direction_matrix) <- species

for (i in 1:nrow(direction_matrix)) {
  for (j in 1:ncol(direction_matrix)) {
    direction_matrix[i,j] <- ifelse(corr_df[i,j] > 0 & pval_df[i,j]<=0.05,+2,ifelse(corr_df[i,j] > 0 & 0.05<pval_df[i,j] & pval_df[i,j]<=0.01,1,ifelse(corr_df[i,j] < 0 & pval_df[i,j]<=0.05,-2, ifelse(corr_df[i,j] < 0 & 0.05<pval_df[i,j] & pval_df[i,j]<=0.01,-1,0))))
    direction_matrix[is.na(direction_matrix)] <- 0
  }
}
direction_matrix <- as.data.frame(direction_matrix)

write.xlsx(direction_matrix,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/corr_direction_matrix.xlsx", rowNames = T)

#################################################################
## Plotting ##

load("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/features_bifs_correlation.RData")

imp_features_corr_df <- uniGroup_BifScore_corr_direct[rowSums(abs(uniGroup_BifScore_corr_direct)) >= 2,]

imp_features_corr_df_t <- t(imp_features_corr_df)
rownames(imp_features_corr_df_t) <- gsub("corr","score", rownames(imp_features_corr_df_t))

imp_features_corr_df_annotation <- read.xlsx("G:\\My Drive\\Bifido_Project\\New_FunctionalAnalysis\\final_feature_bif_association_NEW_modified.xlsx", sheetName = 'Sheet1', check.names =F)
imp_features_corr_df_annotation <- imp_features_corr_df_annotation[,-4]
rownames(imp_features_corr_df_annotation) <- imp_features_corr_df_annotation$Groups
imp_features_corr_df_annotation <- imp_features_corr_df_annotation[colnames(imp_features_corr_df_t),]
imp_features_corr_df_annotation$positive_association <- rowSums(imp_features_corr_df_annotation[, 5:13] == 1, na.rm = TRUE)
imp_features_corr_df_annotation$negative_association <- rowSums(imp_features_corr_df_annotation[, 5:13] == -1, na.rm = TRUE)
imp_features_corr_df_annotation$total_associations <- rowSums(imp_features_corr_df_annotation[,c(14,15)])


imp_features_corr_df_annotation_unique <- imp_features_corr_df_annotation %>%
  group_by(Additional_Annotation) %>%
  slice_max(total_associations, with_ties = FALSE) %>%
  ungroup()

write.xlsx(imp_features_corr_df_annotation_unique,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/imp_features_corr_df_annotation_unique.xlsx")

imp_features_corr_df_t <- as.data.frame(imp_features_corr_df_t)
imp_features_corr_df_t_filtered <- imp_features_corr_df_t[,names(imp_features_corr_df_t) %in% imp_features_corr_df_annotation_unique$Groups]
imp_features_corr_df_t_filtered <- imp_features_corr_df_t_filtered[,imp_features_corr_df_annotation_unique$Groups]
colnames(imp_features_corr_df_t_filtered) <- imp_features_corr_df_annotation_unique$Additional_Annotation

### Now 'Response to toxic substance' and 'Lactate dehydrogenase' have two types of naming so keeping one based on the associations (positive associations)
## So we need to remove column101 (Response to toxic substance) and column71 (Lactate dehydrogenase)

imp_features_corr_df_t_filtered <- imp_features_corr_df_t_filtered[,-c(71,101)]

custom_colors <- c("#e31a1c", "white","#377eb8")

breaks <- seq(-1.5, 1.5, length.out = 4)

pdf("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/features_corr_heatmap.pdf", height = 8, width = 20)
p <- pheatmap::pheatmap(
  imp_features_corr_df_t_filtered,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,       
  treeheight_col = 0,       
  color = custom_colors,
  breaks = breaks,
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("-1", "0", "1"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cellheight = 9,
  cellwidth = 7,
  fontsize_row = 8,
  fontsize_col = 7,
  main = "Associating feaures with Bif scores"
)
dev.off()

ordered_rows <- rownames(imp_features_corr_df_t_filtered)[p$tree_row$order]
ordered_cols <- colnames(imp_features_corr_df_t_filtered)[p$tree_col$order]

carpet_df <- imp_features_corr_df_t_filtered[ordered_rows, ordered_cols]
write.xlsx(carpet_df,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/annotation_heatmap_carpet.xlsx")

imp_features_corr_df_t_groups <- as.data.frame(imp_features_corr_df_t)
imp_features_corr_df_t_groups <- imp_features_corr_df_t_groups[,names(imp_features_corr_df_t_groups) %in% imp_features_corr_df_annotation_unique$Groups]
imp_features_corr_df_t_groups <- imp_features_corr_df_t_groups[,imp_features_corr_df_annotation_unique$Groups]
imp_features_corr_df_t_groups <- imp_features_corr_df_t_groups[,-c(71,101)]

p <- pheatmap::pheatmap(
  imp_features_corr_df_t_groups,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,       
  treeheight_col = 0,       
  color = custom_colors,
  breaks = breaks,
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("-1", "0", "1"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cellheight = 9,
  cellwidth = 7,
  fontsize_row = 8,
  fontsize_col = 7,
  main = "Associating feaures with Bif scores"
)
dev.off()

ordered_rows <- rownames(imp_features_corr_df_t_groups)[p$tree_row$order]
ordered_cols <- colnames(imp_features_corr_df_t_groups)[p$tree_col$order]

carpet_df_groups <- imp_features_corr_df_t_groups[ordered_rows, ordered_cols]
write.xlsx(carpet_df_groups,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/annotation_heatmap_carpet_groups.xlsx")

imp_features_corr_df_annotation_unique <- imp_features_corr_df_annotation_unique[-c(71,101),]

annotation_unique_features142 <- annotation_unique[annotation_unique$Feature %in% imp_features_corr_df_annotation_unique$Groups,]
annotation_unique_features142 <- as.data.frame(annotation_unique_features142)
rownames(annotation_unique_features142) <- annotation_unique_features142$Feature

annotation_unique_features142 <- annotation_unique_features142[imp_features_corr_df_annotation_unique$Groups,]


#################################################################
## Annotation homogenization of important features ##

load("G:/My Drive/Bifido_Project/mean_FunctionalProfileAll_groups.RData")

annotation_imp_features_corr_df <- annotation_df_mean_FunctionalProfile[annotation_df_mean_FunctionalProfile$Groups %in% rownames(imp_features_corr_df),]
rownames(annotation_imp_features_corr_df) <- annotation_imp_features_corr_df$Groups
annotation_imp_features_corr_df <- annotation_imp_features_corr_df[rownames(imp_features_corr_df),]

#################################################################
## Revised Annotations 
#################################################################

load("subset66_127_annotation_imp_features_corr_df.RData")
names(subset_annotation_imp_features_corr_df)[5] <- "Additional_Annotation"

library(xlsx)
subset1_65 <- read.xlsx("FeatureAnnotation65.xlsx")
rownames(subset1_65) <- subset1_65$Groups

overall_unique_feature_annotation <- as.data.frame(rbind(subset1_65,subset_annotation_imp_features_corr_df))
overall_unique_feature_annotation_modified <- overall_unique_feature_annotation[,-4]

feature_bif_associations <- imp_features_corr_df
#feature_bif_associations <- feature_bif_associations[,c(1,2,7,5,4,8,6,3)]
names(feature_bif_associations) <- gsub("_corr","",names(feature_bif_associations))
all(rownames(feature_bif_associations) == rownames(overall_unique_feature_annotation_modified))

final_feature_bif_association <- as.data.frame(cbind(overall_unique_feature_annotation_modified,feature_bif_associations))

write.xlsx(final_feature_bif_association,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/final_feature_bif_association.xlsx")

final_feature_bif_association_old <- read.xlsx("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/final_feature_bif_association_old.xlsx", sheet = 'Sheet1', check.names = F)

rownames(final_feature_bif_association_old) <- final_feature_bif_association_old$Groups
length(intersect(rownames(final_feature_bif_association_old),rownames(final_feature_bif_association_NEW)))

final_feature_bif_association_NEW <- as.data.frame(cbind(annotation_imp_features_corr_df,imp_features_corr_df))

final_feature_bif_association_NEW$Additional_Annotation <- final_feature_bif_association_old[rownames(final_feature_bif_association_NEW), "Additional.Annotation"]
final_feature_bif_association_NEW <- final_feature_bif_association_NEW[,c(1,2,3,4,14,5,6,7,8,9,10,11,12,13)]

write.xlsx(imp_features_corr_df_annotation_unique,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/final_feature_bif_association_NEW.xlsx")

######### Adding the Detection and re-run ###############

load("df_association_all.RData")

detection_association_score <- df_association_all_wgs[,c('adult_detection','senior_detection')]
detection_association_score_selected <- detection_association_score[rownames(df_association_CommonSpecies_modified),]
all(rownames(df_association_CommonSpecies_modified) == rownames(detection_association_score_selected))
detection_association_score_selected$detection <- rowMeans(detection_association_score_selected)

df_association_CommonSpecies_modified$detection <- detection_association_score_selected$detection

######################################################################
## Adding additional columns

final_feature_bif_association_NEW <- read.xlsx("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/final_feature_bif_association_NEW.xlsx", sheetName = 'Sheet1')
all(annotation_unique_features142$Feature == final_feature_bif_association_NEW$Groups)
final_feature_bif_association_NEW$Bif_Association <- annotation_unique_features142$Bif_score
final_feature_bif_association_NEW$Detection_NonBif <- annotation_unique_features142$Detection_NonBif
final_feature_bif_association_NEW$Detection_NonBif_list <- annotation_unique_features142$Detection_NonBif_list

write.xlsx(final_feature_bif_association_NEW, "G:/My Drive/Bifido_Project/New_FunctionalAnalysis/final_feature_bif_association_NEW_selected.xlsx")

###################################################################
## Running RF with features covering 90% Cumulative FI 

rf_models_new <- list()

correlation_results_top <- data.frame(
  Variable = character(),
  Spearman_rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE)

for (var_name in colnames(df_association_CommonSpecies_modified)) {
  
  species_name <- gsub("Bifidobacterium_", "", var_name)
  species_key <- paste0("Bifidobacterium_", species_name)

  top_features <- top_features_rf_MeanFunctionalProfile[[species_key]]$Feature
  
  selected_features_df <- mean_FunctionaProfile_groups[, top_features, drop = FALSE]
  
  response <- df_association_CommonSpecies_modified[[var_name]]
  
  rf_model <- randomForest(x = selected_features_df, y = response)
  
  rf_models_new[[var_name]] <- rf_model
  
  predicted <- rf_model$predicted
  
  cor_test <- cor.test(predicted, response, method = "spearman", exact = T)
  
  correlation_results_top <- rbind(correlation_results_top, data.frame(
    Variable = var_name,
    Spearman_rho = cor_test$estimate,
    p_value = cor_test$p.value,
    stringsAsFactors = FALSE
  ))
}

rownames(correlation_results_top) <- gsub("Bifidobacterium_", "", correlation_results_top$Variable)
rownames(correlation_results_top) <- paste0(rownames(correlation_results_top), "_score")
correlation_results_top$Variable <- NULL

write.xlsx(correlation_results_top,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/correlation_results_MeanFunctionalProfile_topFeatures.xlsx")

#### RF corr scatter plots ###

library(ggplot2)

output_dir <- "G:/My Drive/Bifido_Project/New_FunctionalAnalysis/"

for (var_name in names(rf_models_new)) {
  
  species_label <- gsub("Bifidobacterium_", "", var_name)
  plot_filename <- paste0(output_dir, species_label, "_score_corrplot.pdf")
  
  rf_model <- rf_models_new[[var_name]]
  predicted <- rf_model$predicted
  actual <- df_association_CommonSpecies_modified[[var_name]]
  
  plot_df <- data.frame(
    Predicted = predicted,
    Actual = actual
  )
  
  p <- ggplot(plot_df, aes(x = Actual, y = Predicted)) +
    geom_point(size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "blue3", linewidth = 1.1) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 16),   
      axis.text.x = element_text(size = 16),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.05, "inches")
    )
  
  ggsave(plot_filename, plot = p, height = 4, width = 4)
}


####################### New Association Heatmap #######
input_data <- read.xlsx("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/temp_input.xlsx", sheetName = 'Sheet1')
input_df <- input_data[,6:14]
rownames(input_df) <- input_data$Annotation

input_df_t <- as.data.frame(t(input_df))



custom_colors <- c("#e31a1c", "white","#377eb8")

breaks <- seq(-1.5, 1.5, length.out = 4)

pdf("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/features_corr_heatmap_NEW.pdf", height = 20, width = 8)
grid::grid.newpage()
p <- pheatmap::pheatmap(
  input_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,       
  treeheight_col = 0,       
  color = custom_colors,
  breaks = breaks,
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("-1", "0", "1"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cellheight = 9,
  cellwidth = 7,
  fontsize_row = 8,
  fontsize_col = 7,
  angle_col = 90,
  main = "Associating feaures with Bif scores"
)
dev.off()




custom_colors <- c("#e31a1c", "white","#377eb8")
breaks <- seq(-1, 1, length.out = 3)
col_fun <- colorRamp2(breaks, custom_colors)

ht <- Heatmap(
  as.matrix(input_df),
  name = "score",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 90,
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),  
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 0.5), 
  heatmap_legend_param = list(
    at = c(-1, 0, 1),
    labels = c("-1", "0", "1")), 
  column_title = "Associating features with Bif scores",
  width = unit(ncol(as.matrix(input_df)) * 3, "mm"),
  height = unit(nrow(as.matrix(input_df)) * 3, "mm"))

pdf("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/features_corr_heatmap_NEW.pdf",
    height = 20, width = 8)
draw(ht)
dev.off()


#ordered_rows <- rownames(input_df_t)[p$tree_row$order]
#ordered_cols <- colnames(input_df_t)[p$tree_col$order]

clustered_functions <- rownames(as.matrix(input_df))[row_order(draw(ht))]
clustered_Bifs <- colnames(as.matrix(input_df))[column_order(draw(ht))]

carpet_df <- input_df[clustered_functions,clustered_Bifs]

carpet_df <- input_df_t[ordered_rows, ordered_cols]
write.xlsx(carpet_df,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/annotation_heatmap_carpet_NEW.xlsx")

rownames(input_data) <- input_data$Annotation
input_data_ordered <- input_data[rownames(carpet_df),]
input_data_ordered[,6:14] <- carpet_df[,1:9]
names(input_data_ordered)[6:14] <- names(carpet_df)
write.xlsx(input_data_ordered,"G:/My Drive/Bifido_Project/New_FunctionalAnalysis/input_data_ordered.xlsx")

###############################################################

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridisLite)
library(colorspace)
library(scales)
library(purrr)

broad_category <- read.xlsx("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/input_data_ordered.xlsx", sheetName = 'Sheet1')
broad_category_positive <- broad_category[broad_category$Groups %in% rownames(consistent_positive),]
broad_category_negative <- broad_category[broad_category$Groups %in% rownames(consistent_negative),]

############ Prepare Positive Data ############
df <- as.data.frame(table(broad_category_positive$Broad_process_category))
colnames(df) <- c("Category", "Count")

df <- df %>%
  mutate(Percentage = Count / sum(Count) * 100,
         Label = paste0(round(Percentage, 1), "%"))

############ Prepare Negative Data ############
df_neg <- as.data.frame(table(broad_category_negative$Broad_process_category))
colnames(df_neg) <- c("Category", "Count")

df_neg <- df_neg %>%
  mutate(Percentage = Count / sum(Count) * 100,
         Label = paste0(round(Percentage, 1), "%"))

############ Color Assignment ############
# Base colors for positive categories

base_pal_pos <- brewer.pal(8, "Pastel2")
cats_pos <- unique(df$Category)
colors_pos <- setNames(colorRampPalette(base_pal_pos)(length(cats_pos)), cats_pos)

# Find new categories in negative that are not in positive
cats_neg <- unique(df_neg$Category)
new_cats <- setdiff(cats_neg, cats_pos)

# Assign extra colors for new categories
base_pal_extra <- c(brewer.pal(9, "Pastel1"), brewer.pal(12, "Set3"))
extra_colors <- setNames(colorRampPalette(base_pal_extra)(length(new_cats)), new_cats)

# Final color mapping for negative categories
colors_neg <- c(colors_pos, extra_colors)

# Attach color mapping to data frames
df$Color <- colors_pos[df$Category]
df_neg$Color <- colors_neg[df_neg$Category]

############ Plot Positive Pie Chart ############
pdf("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/Positive_functions.pdf", height = 3.5, width = 7)
ggplot(df, aes(x = "", y = Count, fill = Category)) +
  geom_col(color = "gray40", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            color = "black", size = 3.5) +
  scale_fill_manual(values = colors_pos) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 11)
  )
dev.off()

############ Plot Negative Pie Chart ############
pdf("G:/My Drive/Bifido_Project/New_FunctionalAnalysis/Negative_functions.pdf", height = 3.5, width = 7)
ggplot(df_neg, aes(x = "", y = Count, fill = Category)) +
  geom_col(color = "gray40", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            color = "black", size = 3.5) +
  scale_fill_manual(values = colors_neg) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 11)
  )
dev.off()



library(tidyr)
input_data_ordered_expanded <- input_data_ordered_NEW %>%
  separate_rows(Feature_ids, sep = ";")





