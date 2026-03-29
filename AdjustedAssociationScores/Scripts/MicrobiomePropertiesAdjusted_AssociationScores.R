############ Analysis of Microbial Load with Metacardis Data #########

library(openxlsx)
library(psych)
library(dplyr)


load("/data/bifido_45809_new_2014_not_normalized_Harnandez_added.RData")
load("/data/spdata_Bifs_detect_46218.RData")
load("/data/df_Strong_AssociationScores_new_wgs.RData")

SpeciesProfile_All[is.na(SpeciesProfile_All)] <- 0
SpeciesProfile_All <- SpeciesProfile_All/rowSums(SpeciesProfile_All)

names(SpeciesProfile_All)[1265] <- "Intestinibacter_bartlettii"
names(SpeciesProfile_All)[1362] <- "Erysipelatoclostridium_ramosum"


metacardis_load <- read.xlsx("/data/MetaCardis_Experimental_Load.xlsx")

metacardis_SpeciesProfile <- SpeciesProfile_All[rownames(metadata_updated[metadata_updated$study_name == "MetaCardis_2020_a",]),]
metacardis_metadata <- metadata_updated[rownames(metacardis_SpeciesProfile),]

metacardis_SpeciesProfile_common <- metacardis_SpeciesProfile[intersect(metacardis_load$ID,rownames(metacardis_SpeciesProfile)),]
metacardis_metadata_common <- metacardis_metadata[rownames(metacardis_SpeciesProfile_common),]
metacardis_load_common <- metacardis_load[metacardis_load$ID %in% rownames(metacardis_SpeciesProfile_common),]

all(metacardis_load_common$ID == rownames(metacardis_SpeciesProfile_common))

Bifidobacterium_detection <- rowSums(spdata_Bifs_detect_46218)
metacardis_SpeciesProfile_common$Bifidobacterium_detection <- Bifidobacterium_detection[rownames(metacardis_SpeciesProfile_common)]

MajorBifs <- c("Bifidobacterium_animalis","Bifidobacterium_catenulatum","Bifidobacterium_breve","Bifidobacterium_dentium","Bifidobacterium_bifidum","Bifidobacterium_pseudocatenulatum","Bifidobacterium_adolescentis","Bifidobacterium_longum","Bifidobacterium_detection")

metacardis_BifsProfile <- metacardis_SpeciesProfile_common[,names(metacardis_SpeciesProfile_common) %in% MajorBifs]

metacardis_SpeciesProfile_common_adult <- metacardis_SpeciesProfile_common[rownames(metacardis_SpeciesProfile_common) %in% rownames(metacardis_metadata_common[metacardis_metadata_common$age_category == 'adult',]),]
metacardis_SpeciesProfile_common_senior <- metacardis_SpeciesProfile_common[rownames(metacardis_SpeciesProfile_common) %in% rownames(metacardis_metadata_common[metacardis_metadata_common$age_category == 'senior',]),]

metacardis_BifsProfile_adult <- metacardis_BifsProfile[rownames(metacardis_SpeciesProfile_common_adult),]
metacardis_BifsProfile_senior <- metacardis_BifsProfile[rownames(metacardis_SpeciesProfile_common_senior),]

metacardis_SpeciesProfile_common_adult_selected <- metacardis_SpeciesProfile_common_adult[,names(metacardis_SpeciesProfile_common_adult) %in% rownames(df_association_adult_detection_wgs_new)]
metacardis_SpeciesProfile_common_senior_selected <- metacardis_SpeciesProfile_common_senior[,names(metacardis_SpeciesProfile_common_senior) %in% rownames(df_association_senior_detection_wgs_new)]

metacardis_load_common_adult <- metacardis_load_common[metacardis_load_common$ID %in% rownames(metacardis_SpeciesProfile_common_adult_selected),]
metacardis_load_common_senior <- metacardis_load_common[metacardis_load_common$ID %in% rownames(metacardis_SpeciesProfile_common_senior_selected),]

qmp_metacardis_common_adult <- metacardis_SpeciesProfile_common_adult_selected * metacardis_load_common_adult$Load
qmp_metacardis_common_senior <- metacardis_SpeciesProfile_common_senior_selected * metacardis_load_common_senior$Load

qmp_Bifs_adult <- metacardis_BifsProfile_adult[,1:8] * metacardis_load_common_adult$Load
qmp_Bifs_adult$Bifidobacterium_detection <- metacardis_BifsProfile_adult$Bifidobacterium_detection

qmp_Bifs_senior <- metacardis_BifsProfile_senior[,1:8] * metacardis_load_common_senior$Load
qmp_Bifs_senior$Bifidobacterium_detection <- metacardis_BifsProfile_senior$Bifidobacterium_detection

AssociationScores_adult_wgs <- data.frame(Bifidobacterium_detection = df_association_adult_detection_wgs_new$score,
                                          Bifidobacterium_longum = df_association_adult_longum_wgs_new$score,
                                          Bifidobacterium_adolescentis = df_association_adult_adolescentis_wgs_new$score,
                                          Bifidobacterium_pseudocatenulatum = df_association_adult_pseudocatenulatum_wgs_new$score,
                                          Bifidobacterium_dentium = df_association_adult_dentium_wgs_new$score,
                                          Bifidobacterium_bifidum = df_association_adult_bifidum_wgs_new$score,
                                          Bifidobacterium_breve = df_association_adult_breve_wgs_new$score,
                                          Bifidobacterium_catenulatum = df_association_adult_catenulatum_wgs_new$score,
                                          Bifidobacterium_animalis = df_association_adult_animalis_wgs_new$score,
                                          row.names = rownames(df_association_adult_detection_wgs_new))

AssociationScores_senior_wgs <- data.frame(Bifidobacterium_detection = df_association_senior_detection_wgs_new$score,
                                           Bifidobacterium_longum = df_association_senior_longum_wgs_new$score,
                                           Bifidobacterium_adolescentis = df_association_senior_adolescentis_wgs_new$score,
                                           Bifidobacterium_pseudocatenulatum = df_association_senior_pseudocatenulatum_wgs_new$score,
                                           Bifidobacterium_dentium = df_association_senior_dentium_wgs_new$score,
                                           Bifidobacterium_bifidum = df_association_senior_bifidum_wgs_new$score,
                                           Bifidobacterium_breve = df_association_senior_breve_wgs_new$score,
                                           Bifidobacterium_catenulatum = df_association_senior_catenulatum_wgs_new$score,
                                           Bifidobacterium_animalis = df_association_senior_animalis_wgs_new$score,
                                           row.names = rownames(df_association_senior_detection_wgs_new))

qmp_metacardis_common_adult <- qmp_metacardis_common_adult[,rownames(AssociationScores_adult_wgs)]
qmp_metacardis_common_senior <- qmp_metacardis_common_senior[,rownames(AssociationScores_senior_wgs)]

###################################################################
#### Adult Microbiomes #####

qmp_association_df_adult <- matrix(
  NA,
  nrow = ncol(qmp_metacardis_common_adult),
  ncol = ncol(qmp_Bifs_adult),
  dimnames = list(
    colnames(qmp_metacardis_common_adult),
    colnames(qmp_Bifs_adult)
  )
)


for (i in seq_along(qmp_metacardis_common_adult)) {
  for (j in seq_along(qmp_Bifs_adult)) {
    
    test <- cor.test(qmp_metacardis_common_adult[, i],
                     qmp_Bifs_adult[, j],
                     method = "spearman",
                     exact = TRUE)
    
    qmp_association_df_adult[i, j] <- unname(test$estimate)
  }
}

qmp_association_df_adult <- as.data.frame(qmp_association_df_adult)

qmp_association_df_adult <- qmp_association_df_adult[,names(AssociationScores_adult_wgs)]


##### Correlation between QMP Association and Association Scores

qmp_corr_result_df_adult <- data.frame(
  corr = rep(NA, ncol(qmp_association_df_adult)),
  pvalue = rep(NA, ncol(qmp_association_df_adult)),
  row.names = colnames(qmp_association_df_adult)
)


for (i in seq_along(qmp_association_df_adult)) {
  
  test <- cor.test(
    AssociationScores_adult_wgs[, i],
    qmp_association_df_adult[, i],
    method = "spearman",
    exact = TRUE)
  
  qmp_corr_result_df_adult[i, "corr"]   <- unname(test$estimate)
  qmp_corr_result_df_adult[i, "pvalue"] <- test$p.value
}

###################################################################
#### Senior Microbiomes #####

qmp_association_df_senior <- matrix(
  NA,
  nrow = ncol(qmp_metacardis_common_senior),
  ncol = ncol(qmp_Bifs_senior),
  dimnames = list(
    colnames(qmp_metacardis_common_senior),
    colnames(qmp_Bifs_senior)
  )
)


for (i in seq_along(qmp_metacardis_common_senior)) {
  for (j in seq_along(qmp_Bifs_senior)) {
    
    test <- cor.test(qmp_metacardis_common_senior[, i],
                     qmp_Bifs_senior[, j],
                     method = "spearman",
                     exact = TRUE)
    
    qmp_association_df_senior[i, j] <- unname(test$estimate)
  }
}

qmp_association_df_senior <- as.data.frame(qmp_association_df_senior)

qmp_association_df_senior <- qmp_association_df_senior[,names(AssociationScores_senior_wgs)]


##### Correlation between QMP Association and Association Scores

qmp_corr_result_df_senior <- data.frame(
  corr = rep(NA, ncol(qmp_association_df_senior)),
  pvalue = rep(NA, ncol(qmp_association_df_senior)),
  row.names = colnames(qmp_association_df_senior)
)


for (i in seq_along(qmp_association_df_senior)) {
  
  test <- cor.test(
    AssociationScores_senior_wgs[, i],
    qmp_association_df_senior[, i],
    method = "spearman",
    exact = TRUE)
  
  qmp_corr_result_df_senior[i, "corr"]   <- unname(test$estimate)
  qmp_corr_result_df_senior[i, "pvalue"] <- test$p.value
}

##################################################################
################# Using Relative Abundance ##########
###################################################################
#### Adult Microbiomes #####

metacardis_rel_abundance_common_adult <- metacardis_SpeciesProfile_common_adult_selected[,rownames(AssociationScores_adult_wgs)]
metacardis_rel_abundance_common_senior <- metacardis_SpeciesProfile_common_senior_selected[,rownames(AssociationScores_adult_wgs)]

rel_ab_association_df_adult <- matrix(
  NA,
  nrow = ncol(metacardis_rel_abundance_common_adult),
  ncol = ncol(metacardis_BifsProfile_adult),
  dimnames = list(
    colnames(metacardis_rel_abundance_common_adult),
    colnames(metacardis_BifsProfile_adult)
  )
)


for (i in seq_along(metacardis_rel_abundance_common_adult)) {
  for (j in seq_along(metacardis_BifsProfile_adult)) {
    
    test <- cor.test(metacardis_rel_abundance_common_adult[, i],
                     metacardis_BifsProfile_adult[, j],
                     method = "spearman",
                     exact = TRUE)
    
    rel_ab_association_df_adult[i, j] <- unname(test$estimate)
  }
}

rel_ab_association_df_adult <- as.data.frame(rel_ab_association_df_adult)

rel_ab_association_df_adult <- rel_ab_association_df_adult[,names(AssociationScores_adult_wgs)]


##### Correlation between QMP Association and Association Scores

rel_ab_corr_result_df_adult <- data.frame(
  corr = rep(NA, ncol(rel_ab_association_df_adult)),
  pvalue = rep(NA, ncol(rel_ab_association_df_adult)),
  row.names = colnames(rel_ab_association_df_adult)
)


for (i in seq_along(rel_ab_association_df_adult)) {
  
  test <- cor.test(
    AssociationScores_adult_wgs[, i],
    rel_ab_association_df_adult[, i],
    method = "spearman",
    exact = TRUE)
  
  rel_ab_corr_result_df_adult[i, "corr"]   <- unname(test$estimate)
  rel_ab_corr_result_df_adult[i, "pvalue"] <- test$p.value
}

###################################################################
#### Senior Microbiomes #####

metacardis_rel_abundance_common_senior <- metacardis_SpeciesProfile_common_senior_selected[,rownames(AssociationScores_senior_wgs)]
metacardis_rel_abundance_common_senior <- metacardis_SpeciesProfile_common_senior_selected[,rownames(AssociationScores_senior_wgs)]

rel_ab_association_df_senior <- matrix(
  NA,
  nrow = ncol(metacardis_rel_abundance_common_senior),
  ncol = ncol(metacardis_BifsProfile_senior),
  dimnames = list(
    colnames(metacardis_rel_abundance_common_senior),
    colnames(metacardis_BifsProfile_senior)
  )
)


for (i in seq_along(metacardis_rel_abundance_common_senior)) {
  for (j in seq_along(metacardis_BifsProfile_senior)) {
    
    test <- cor.test(metacardis_rel_abundance_common_senior[, i],
                     metacardis_BifsProfile_senior[, j],
                     method = "spearman",
                     exact = TRUE)
    
    rel_ab_association_df_senior[i, j] <- unname(test$estimate)
  }
}

rel_ab_association_df_senior <- as.data.frame(rel_ab_association_df_senior)

rel_ab_association_df_senior <- rel_ab_association_df_senior[,names(AssociationScores_senior_wgs)]


##### Correlation between QMP Association and Association Scores

rel_ab_corr_result_df_senior <- data.frame(
  corr = rep(NA, ncol(rel_ab_association_df_senior)),
  pvalue = rep(NA, ncol(rel_ab_association_df_senior)),
  row.names = colnames(rel_ab_association_df_senior)
)


for (i in seq_along(rel_ab_association_df_senior)) {
  
  test <- cor.test(
    AssociationScores_senior_wgs[, i],
    rel_ab_association_df_senior[, i],
    method = "spearman",
    exact = TRUE)
  
  rel_ab_corr_result_df_senior[i, "corr"]   <- unname(test$estimate)
  rel_ab_corr_result_df_senior[i, "pvalue"] <- test$p.value
}


#####################################################################

# library(xlsx)
# write.xlsx(qmp_corr_result_df_adult,"G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\QMP_analysis\\qmp_corr_result_df_adult.xlsx")
# write.xlsx(qmp_corr_result_df_senior,"G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\QMP_analysis\\qmp_corr_result_df_senior.xlsx")

rownames(qmp_corr_result_df_adult) <- gsub("Bifidobacterium","",rownames(qmp_corr_result_df_adult))
rownames(qmp_corr_result_df_senior) <- gsub("Bifidobacterium","",rownames(qmp_corr_result_df_senior))

rownames(qmp_corr_result_df_adult) <- gsub("_","",rownames(qmp_corr_result_df_adult))
rownames(qmp_corr_result_df_senior) <- gsub("_","",rownames(qmp_corr_result_df_senior))

qmp_corr_result_df_adult <- qmp_corr_result_df_adult[c(1,2,3,4,6,5,7,8,9),]
qmp_corr_result_df_senior <- qmp_corr_result_df_senior[c(1,2,3,4,6,5,7,8,9),]


library(dplyr)
library(ggplot2)
library(tibble)

## =========================
## Color scheme (same as before)
## =========================
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

## =========================
## Prepare dataframe
## =========================
plot_df <- qmp_corr_result_df_adult %>%
  rownames_to_column("bif") %>%
  mutate(
    bif = factor(bif, levels = bif),   # preserve order left → right
    star_label = ifelse(pvalue <= 0.05, "*", ""),
    star_y = corr + 0.02 * diff(range(corr, na.rm = TRUE))
  )

## =========================
## Vertical bar plot
## =========================
ggplot(plot_df, aes(x = bif, y = corr, fill = bif)) +
  geom_col(width = 0.7) +
  
  geom_text(
    aes(y = star_y, label = star_label),
    color = "black",
    size = 5
  ) +
  
  scale_fill_manual(values = bif_colors, drop = TRUE) +
  
  labs(
    x = NULL,
    y = "Correlation",
    title = "Adult QMP Correlations"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

## =========================
## Prepare dataframe
## =========================
plot_df <- qmp_corr_result_df_senior %>%
  rownames_to_column("bif") %>%
  mutate(
    bif = factor(bif, levels = bif),   # preserve order left → right
    star_label = ifelse(pvalue <= 0.05, "*", ""),
    star_y = corr + 0.02 * diff(range(corr, na.rm = TRUE))
  )

## =========================
## Vertical bar plot
## =========================
ggplot(plot_df, aes(x = bif, y = corr, fill = bif)) +
  geom_col(width = 0.7) +
  
  geom_text(
    aes(y = star_y, label = star_label),
    color = "black",
    size = 5
  ) +
  
  scale_fill_manual(values = bif_colors, drop = TRUE) +
  
  labs(
    x = NULL,
    y = "Correlation",
    title = "Senior QMP Correlations"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


####################################################################
################  Adjusted Association Scores ####################

load("/data/QMP_workspace.RData")
load("/data/Shannon_Pielou_Load_df.RData")
load("/data/full_AssociationScoresAgeCategory_wgs_original.RData")
load("/data/full_AssociationScoresCohortType_wgs_original.RData")
load("/data/bifido_45809_new_2014_not_normalized_Harnandez_added.RData")

library(dplyr)
library(pcaPP)
library(psych)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ppcor)

grep("_wgs$",ls(),value = T)

SpeciesProfile_All[is.na(SpeciesProfile_All)] <- 0

spdata_45809_wgs <- SpeciesProfile_All[rownames(wgs_SpeciesProfile_filtered_norm),]
shannon_pielou_load_wgs <- shannon_pielou_load_df[rownames(wgs_SpeciesProfile_filtered_norm),]
spdata_45809_wgs$Bifidobacterium_detection <- apply(spdata_45809_wgs[, AllBifidoSpecies], 1, function(x) length(x[x >= 0.0001]))
spdata_45809_wgs$load <- shannon_pielou_load_wgs$load
spdata_45809_wgs$shannon <- shannon_pielou_load_wgs$shannon_index
spdata_45809_wgs$pielou <- shannon_pielou_load_wgs$pielou
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

bifido_species <- c("Bifidobacterium_detection","Bifidobacterium_longum","Bifidobacterium_adolescentis","Bifidobacterium_pseudocatenulatum","Bifidobacterium_bifidum","Bifidobacterium_dentium","Bifidobacterium_breve","Bifidobacterium_catenulatum","Bifidobacterium_animalis")

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
###########    PARTIAL CORRELATION    ###############

############## AGE CATEGORY Wise ############

compute_partial_corr_associations_by_AgeCategory <- function(spdata,df_association,age_category,bif_name,bif_studies) {
  if (!requireNamespace("ppcor", quietly = TRUE)) {
    stop("Please install the 'ppcor' package: install.packages('ppcor')")
  }
  library(ppcor)
  
  # Filter samples by age category
  spdata_age <- spdata[spdata$age_category == age_category & spdata$study_name %in% bif_studies, ]
  
  
  studies <- bif_studies
  
  non_bif_species <- rownames(df_association)
  
  corr_mat <- matrix(NA, nrow = length(non_bif_species), ncol = length(studies),
                     dimnames = list(non_bif_species, studies))
  pval_mat <- corr_mat
  dir_mat  <- corr_mat
  
  for (study in studies) {
    study_df <- subset(spdata_age, study_name == study)
    if (nrow(study_df) < 5) next
    
    for (sp in non_bif_species) {
      if (!all(c(sp, bif_name, "load", "shannon", "pielou") %in% colnames(study_df))) next
      
      x <- study_df[[sp]]
      y <- study_df[[bif_name]]
      covars <- study_df[, c("load", "shannon", "pielou")]
      
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

df_association_infant_detection_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_detection_wgs,"infant","Bifidobacterium_detection",Bifidobacterium_detection_infant_studies)
df_association_infant_longum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_longum_wgs,"infant","Bifidobacterium_longum",Bifidobacterium_longum_infant_studies)
df_association_infant_adolescentis_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_adolescentis_wgs,"infant","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_infant_studies)
df_association_infant_pseudocatenulatum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_pseudocatenulatum_wgs,"infant","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_infant_studies)
df_association_infant_dentium_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_dentium_wgs,"infant","Bifidobacterium_dentium",Bifidobacterium_dentium_infant_studies)
df_association_infant_bifidum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_bifidum_wgs,"infant","Bifidobacterium_bifidum",Bifidobacterium_bifidum_infant_studies)
df_association_infant_breve_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_breve_wgs,"infant","Bifidobacterium_breve",Bifidobacterium_breve_infant_studies)
df_association_infant_catenulatum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_infant_catenulatum_wgs,"infant","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_infant_studies)


# df_association_adult_detection_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_detection_wgs,"adult","Bifidobacterium_detection",Bifidobacterium_detection_adult_studies)
# df_association_adult_longum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_longum_wgs,"adult","Bifidobacterium_longum",Bifidobacterium_longum_adult_studies)
# df_association_adult_adolescentis_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_adolescentis_wgs,"adult","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_adult_studies)
# df_association_adult_pseudocatenulatum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_pseudocatenulatum_wgs,"adult","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_adult_studies)
# df_association_adult_dentium_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_dentium_wgs,"adult","Bifidobacterium_dentium",Bifidobacterium_dentium_adult_studies)
# df_association_adult_bifidum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_bifidum_wgs,"adult","Bifidobacterium_bifidum",Bifidobacterium_bifidum_adult_studies)
# df_association_adult_breve_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_breve_wgs,"adult","Bifidobacterium_breve",Bifidobacterium_breve_adult_studies)
# df_association_adult_catenulatum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_catenulatum_wgs,"adult","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_adult_studies)
# df_association_adult_animalis_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_adult_animalis_wgs,"adult","Bifidobacterium_animalis",Bifidobacterium_animalis_adult_studies)
# 
# 
# df_association_senior_detection_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_detection_wgs,"senior","Bifidobacterium_detection",Bifidobacterium_detection_senior_studies)
# df_association_senior_longum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_longum_wgs,"senior","Bifidobacterium_longum",Bifidobacterium_longum_senior_studies)
# df_association_senior_adolescentis_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_adolescentis_wgs,"senior","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_senior_studies)
# df_association_senior_pseudocatenulatum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_pseudocatenulatum_wgs,"senior","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_senior_studies)
# df_association_senior_dentium_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_dentium_wgs,"senior","Bifidobacterium_dentium",Bifidobacterium_dentium_senior_studies)
# df_association_senior_bifidum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_bifidum_wgs,"senior","Bifidobacterium_bifidum",Bifidobacterium_bifidum_senior_studies)
# df_association_senior_breve_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_breve_wgs,"senior","Bifidobacterium_breve",Bifidobacterium_breve_senior_studies)
# df_association_senior_catenulatum_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_catenulatum_wgs,"senior","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_senior_studies)
# df_association_senior_animalis_wgs_new <- compute_partial_corr_associations_by_AgeCategory(spdata_45809_wgs, df_association_senior_animalis_wgs,"senior","Bifidobacterium_animalis",Bifidobacterium_animalis_senior_studies)
# 
# 
# ############### Formatting Name of objects ###################
# 
# df_association_infant_detection_wgs_adjusted <- df_association_infant_detection_wgs_new
# df_association_infant_longum_wgs_adjusted <- df_association_infant_longum_wgs_new
# df_association_infant_adolescentis_wgs_adjusted <- df_association_infant_adolescentis_wgs_new
# df_association_infant_pseudocatenulatum_wgs_adjusted <- df_association_infant_pseudocatenulatum_wgs_new
# df_association_infant_dentium_wgs_adjusted <- df_association_infant_dentium_wgs_new
# df_association_infant_bifidum_wgs_adjusted <- df_association_infant_bifidum_wgs_new
# df_association_infant_breve_wgs_adjusted <- df_association_infant_breve_wgs_new
# df_association_infant_catenulatum_wgs_adjusted <- df_association_infant_catenulatum_wgs_new
# 
# df_association_adult_detection_wgs_adjusted <- df_association_adult_detection_wgs_new
# df_association_adult_longum_wgs_adjusted <- df_association_adult_longum_wgs_new
# df_association_adult_adolescentis_wgs_adjusted <- df_association_adult_adolescentis_wgs_new
# df_association_adult_pseudocatenulatum_wgs_adjusted <- df_association_adult_pseudocatenulatum_wgs_new
# df_association_adult_dentium_wgs_adjusted <- df_association_adult_dentium_wgs_new
# df_association_adult_bifidum_wgs_adjusted <- df_association_adult_bifidum_wgs_new
# df_association_adult_breve_wgs_adjusted <- df_association_adult_breve_wgs_new
# df_association_adult_catenulatum_wgs_adjusted <- df_association_adult_catenulatum_wgs_new
# df_association_adult_animalis_wgs_adjusted <- df_association_adult_animalis_wgs_new
# 
# df_association_senior_detection_wgs_adjusted <- df_association_senior_detection_wgs_new
# df_association_senior_longum_wgs_adjusted <- df_association_senior_longum_wgs_new
# df_association_senior_adolescentis_wgs_adjusted <- df_association_senior_adolescentis_wgs_new
# df_association_senior_pseudocatenulatum_wgs_adjusted <- df_association_senior_pseudocatenulatum_wgs_new
# df_association_senior_dentium_wgs_adjusted <- df_association_senior_dentium_wgs_new
# df_association_senior_bifidum_wgs_adjusted <- df_association_senior_bifidum_wgs_new
# df_association_senior_breve_wgs_adjusted <- df_association_senior_breve_wgs_new
# df_association_senior_catenulatum_wgs_adjusted <- df_association_senior_catenulatum_wgs_new
# df_association_senior_animalis_wgs_adjusted <- df_association_senior_animalis_wgs_new
# 
# 
# ###########  Cohort Lifestyle Wise #########################
# 
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
# 
# 
# compute_partial_corr_associations_by_cohort <- function(spdata, df_association, cohort_name, bif_name, bif_studies) {
#   if (!requireNamespace("ppcor", quietly = TRUE)) {
#     stop("Please install the 'ppcor' package: install.packages('ppcor')")
#   }
#   library(ppcor)
#   
#   # Filter samples: only adults and matching cohort
#   spdata_cohort <- spdata[spdata$age_category == "adult" & spdata$cohort %in% cohort_name, ]
#   
#   # Stop if no matching samples
#   if (nrow(spdata_cohort) == 0) {
#     stop(paste("No samples found for cohort:", cohort_name))
#   }
#   
#   # Get list of unique studies
#   studies <- bif_studies
#   
#   # Get all non-bif species
#   non_bif_species <- rownames(df_association)
#   
#   # Initialize matrices
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
#       if (!all(c(sp, bif_name, "load", "shannon", "pielou") %in% colnames(study_df))) next
#       
#       x <- study_df[[sp]]
#       y <- study_df[[bif_name]]
#       covars <- study_df[, c("load", "shannon", "pielou")]
#       
#       if (sum(complete.cases(x, y, covars)) < 5) next
#       
#       test_res <- try(ppcor::pcor.test(x, y, covars, method = "spearman"), silent = TRUE)
#       if (inherits(test_res, "try-error")) next
#       
#       corr_mat[sp, study] <- test_res$estimate
#       pval_mat[sp, study] <- test_res$p.value
#     }
#   }
#   
#   # Replace NA with defaults
#   corr_mat[is.na(corr_mat)] <- 0
#   pval_mat[is.na(pval_mat)] <- 1
#   
#   # Direction matrix
#   dir_mat[corr_mat > 0 & pval_mat <= 0.05] <- 2
#   dir_mat[corr_mat > 0 & pval_mat > 0.05 & pval_mat <= 0.1] <- 1
#   dir_mat[corr_mat < 0 & pval_mat <= 0.05] <- -2
#   dir_mat[corr_mat < 0 & pval_mat > 0.05 & pval_mat <= 0.1] <- -1
#   dir_mat[is.na(dir_mat)] <- 0
#   
#   # Count positive and negative directions
#   positive <- apply(dir_mat, 1, function(x) sum(x %in% c(1, 2)))
#   negative <- apply(dir_mat, 1, function(x) sum(x %in% c(-1, -2)))
#   total <- length(studies)
#   
#   # Score formula
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
#   # Output list
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
# df_association_IndustrializedUrban_detection_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs,df_association_IndustrializedUrban_detection_wgs,"IndustrializedUrban","Bifidobacterium_detection",Bifidobacterium_detection_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_longum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_longum_wgs,"IndustrializedUrban","Bifidobacterium_longum",Bifidobacterium_longum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_adolescentis_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_adolescentis_wgs,"IndustrializedUrban","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_pseudocatenulatum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_pseudocatenulatum_wgs,"IndustrializedUrban","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_dentium_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_dentium_wgs,"IndustrializedUrban","Bifidobacterium_dentium",Bifidobacterium_dentium_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_bifidum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_bifidum_wgs,"IndustrializedUrban","Bifidobacterium_bifidum",Bifidobacterium_bifidum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_breve_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_breve_wgs,"IndustrializedUrban","Bifidobacterium_breve",Bifidobacterium_breve_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_catenulatum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_catenulatum_wgs,"IndustrializedUrban","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_IndustrializedUrban_studies)
# df_association_IndustrializedUrban_animalis_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_IndustrializedUrban_animalis_wgs,"IndustrializedUrban","Bifidobacterium_animalis",Bifidobacterium_animalis_IndustrializedUrban_studies)
# 
# df_association_UrbanRuralMixed_detection_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs,df_association_UrbanRuralMixed_detection_wgs,"UrbanRuralMixed","Bifidobacterium_detection",Bifidobacterium_detection_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_longum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_longum_wgs,"UrbanRuralMixed","Bifidobacterium_longum",Bifidobacterium_longum_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_adolescentis_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_adolescentis_wgs,"UrbanRuralMixed","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_pseudocatenulatum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_pseudocatenulatum_wgs,"UrbanRuralMixed","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_dentium_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_dentium_wgs,"UrbanRuralMixed","Bifidobacterium_dentium",Bifidobacterium_dentium_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_bifidum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_bifidum_wgs,"UrbanRuralMixed","Bifidobacterium_bifidum",Bifidobacterium_bifidum_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_breve_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_breve_wgs,"UrbanRuralMixed","Bifidobacterium_breve",Bifidobacterium_breve_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_catenulatum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_catenulatum_wgs,"UrbanRuralMixed","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_UrbanRuralMixed_studies)
# df_association_UrbanRuralMixed_animalis_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_UrbanRuralMixed_animalis_wgs,"UrbanRuralMixed","Bifidobacterium_animalis",Bifidobacterium_animalis_UrbanRuralMixed_studies)
# 
# df_association_RuralTribal_detection_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs,df_association_RuralTribal_detection_wgs,"RuralTribal","Bifidobacterium_detection",Bifidobacterium_detection_RuralTribal_studies)
# df_association_RuralTribal_longum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_longum_wgs,"RuralTribal","Bifidobacterium_longum",Bifidobacterium_longum_RuralTribal_studies)
# df_association_RuralTribal_adolescentis_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_adolescentis_wgs,"RuralTribal","Bifidobacterium_adolescentis",Bifidobacterium_adolescentis_RuralTribal_studies)
# df_association_RuralTribal_pseudocatenulatum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_pseudocatenulatum_wgs,"RuralTribal","Bifidobacterium_pseudocatenulatum",Bifidobacterium_pseudocatenulatum_RuralTribal_studies)
# df_association_RuralTribal_dentium_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_dentium_wgs,"RuralTribal","Bifidobacterium_dentium",Bifidobacterium_dentium_RuralTribal_studies)
# df_association_RuralTribal_bifidum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_bifidum_wgs,"RuralTribal","Bifidobacterium_bifidum",Bifidobacterium_bifidum_RuralTribal_studies)
# df_association_RuralTribal_breve_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_breve_wgs,"RuralTribal","Bifidobacterium_breve",Bifidobacterium_breve_RuralTribal_studies)
# df_association_RuralTribal_catenulatum_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_catenulatum_wgs,"RuralTribal","Bifidobacterium_catenulatum",Bifidobacterium_catenulatum_RuralTribal_studies)
# df_association_RuralTribal_animalis_wgs_new <- compute_partial_corr_associations_by_cohort(spdata_45809_wgs, df_association_RuralTribal_animalis_wgs,"RuralTribal","Bifidobacterium_animalis",Bifidobacterium_animalis_RuralTribal_studies)
# 
# ############### Formatting Name of objects ###################
# 
# df_association_IndustrializedUrban_detection_wgs_adjusted <- df_association_IndustrializedUrban_detection_wgs_new
# df_association_IndustrializedUrban_longum_wgs_adjusted <- df_association_IndustrializedUrban_longum_wgs_new
# df_association_IndustrializedUrban_adolescentis_wgs_adjusted <- df_association_IndustrializedUrban_adolescentis_wgs_new
# df_association_IndustrializedUrban_pseudocatenulatum_wgs_adjusted <- df_association_IndustrializedUrban_pseudocatenulatum_wgs_new
# df_association_IndustrializedUrban_dentium_wgs_adjusted <- df_association_IndustrializedUrban_dentium_wgs_new
# df_association_IndustrializedUrban_bifidum_wgs_adjusted <- df_association_IndustrializedUrban_bifidum_wgs_new
# df_association_IndustrializedUrban_breve_wgs_adjusted <- df_association_IndustrializedUrban_breve_wgs_new
# df_association_IndustrializedUrban_catenulatum_wgs_adjusted <- df_association_IndustrializedUrban_catenulatum_wgs_new
# df_association_IndustrializedUrban_animalis_wgs_adjusted <- df_association_IndustrializedUrban_animalis_wgs_new
# 
# df_association_UrbanRuralMixed_detection_wgs_adjusted <- df_association_UrbanRuralMixed_detection_wgs_new
# df_association_UrbanRuralMixed_longum_wgs_adjusted <- df_association_UrbanRuralMixed_longum_wgs_new
# df_association_UrbanRuralMixed_adolescentis_wgs_adjusted <- df_association_UrbanRuralMixed_adolescentis_wgs_new
# df_association_UrbanRuralMixed_pseudocatenulatum_wgs_adjusted <- df_association_UrbanRuralMixed_pseudocatenulatum_wgs_new
# df_association_UrbanRuralMixed_dentium_wgs_adjusted <- df_association_UrbanRuralMixed_dentium_wgs_new
# df_association_UrbanRuralMixed_bifidum_wgs_adjusted <- df_association_UrbanRuralMixed_bifidum_wgs_new
# df_association_UrbanRuralMixed_breve_wgs_adjusted <- df_association_UrbanRuralMixed_breve_wgs_new
# df_association_UrbanRuralMixed_catenulatum_wgs_adjusted <- df_association_UrbanRuralMixed_catenulatum_wgs_new
# df_association_UrbanRuralMixed_animalis_wgs_adjusted <- df_association_UrbanRuralMixed_animalis_wgs_new
# 
# df_association_RuralTribal_detection_wgs_adjusted <- df_association_RuralTribal_detection_wgs_new
# df_association_RuralTribal_longum_wgs_adjusted <- df_association_RuralTribal_longum_wgs_new
# df_association_RuralTribal_adolescentis_wgs_adjusted <- df_association_RuralTribal_adolescentis_wgs_new
# df_association_RuralTribal_pseudocatenulatum_wgs_adjusted <- df_association_RuralTribal_pseudocatenulatum_wgs_new
# df_association_RuralTribal_dentium_wgs_adjusted <- df_association_RuralTribal_dentium_wgs_new
# df_association_RuralTribal_bifidum_wgs_adjusted <- df_association_RuralTribal_bifidum_wgs_new
# df_association_RuralTribal_breve_wgs_adjusted <- df_association_RuralTribal_breve_wgs_new
# df_association_RuralTribal_catenulatum_wgs_adjusted <- df_association_RuralTribal_catenulatum_wgs_new
# df_association_RuralTribal_animalis_wgs_adjusted <- df_association_RuralTribal_animalis_wgs_new
# 
# 
# #################################################################
# 
# adjusted_lists <- ls(pattern = "_adjusted$")
# save(list = adjusted_lists, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\Overall_AssociationScores_adjusted_NEW.RData")
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
# save(list = list_names, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\df_AssociationScores_adjusted_wgs.RData")
# 

###################################################################
###################################################################

# ## Strong Association Scores after adjustment
# 
# load("/data/Overall_AssociationScores_adjusted_NEW.RData")
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
# adjusted_lists <- ls(pattern = "_adjusted$")
# save(list = adjusted_lists, file = "/results/Overall_Strong_AssociationScores_adjusted_wgs.RData")

#####################################################################

load("/data/Overall_Strong_AssociationScores_adjusted_wgs.RData")

list_names <- ls()[sapply(ls(), function(x) is.list(get(x)))]

for (lname in list_names) {
  lst <- get(lname)
  
  if (!is.null(lst$strong_association)) { 
    assign(lname, lst$strong_association)  
  } else {
    warning(paste("List", lname, "does not contain 'association'"))
  }
}

save(list = list_names, file = "/results/df_Strong_AssociationScores_adjusted_wgs.RData")

####################################################################
### Combined Association Scores df ######

load("/data/df_Strong_AssociationScores_adjusted_wgs.RData")
df_names <- ls(pattern = "^df_association_.*_wgs_adjusted$")
df_list <- lapply(df_names, function(obj) {
  df <- get(obj)                      # get the data frame
  
  # short name for the column: strip prefix and suffix
  short_name <- sub("^df_association_(.*)_wgs_adjusted$", "\\1", obj)
  
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

df_strong_association_all_adjusted_wgs <- combined
df_strong_association_all_adjusted_wgs[is.na(df_strong_association_all_adjusted_wgs)] <- 0

categories <- c("infant", "adult", "senior", 
                "IndustrializedUrban", "UrbanRuralMixed", "RuralTribal")
bifnames <- c("detection", "longum", "adolescentis", "pseudocatenulatum",
              "dentium", "bifidum", "breve", "catenulatum", "animalis")
desired_order <- unlist(lapply(categories, function(cat) {
  paste0(cat, "_", bifnames)
}))
desired_order <- desired_order[desired_order != "infant_animalis"]

df_strong_association_all_adjusted_wgs <- df_strong_association_all_adjusted_wgs[,desired_order]


load("/data/df_strong_association_all_new_wgs.RData")

df_strong_association_all_adjusted_wgs_filtered <- df_strong_association_all_adjusted_wgs[rownames(df_strong_association_all_new_wgs_filtered), names(df_strong_association_all_new_wgs_filtered)]



### Heatmap ####

col_fun <- colorRamp2(
  c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
  c("#8B0053", "#CC3377", "#F6B6D2", "#FFFFFF", "#CDE9B6", "#77C679", "#004D00")
)

# Convert to numeric matrix
mat <- as.matrix(df_strong_association_all_adjusted_wgs_filtered)

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

pdf("/results/StrongAssociationScores_adjusted_wgs.pdf", 
    height = 20, width = 16)
draw(ht)
dev.off()

###### CARPET ######

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_AssociationScores_Heatmap <- as.data.frame(carpet_df)

#write.xlsx(carpet_AssociationScores_Heatmap, "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\carpet_StrongAssociationScores_adjusted_wgs.xlsx")

save(df_strong_association_all_adjusted_wgs_filtered,file="/results/df_strong_association_all_adjusted_wgs.RData")

#####################################################################
####################################################################
## Original vs Adjusted Correlation (Strong Associations)

load("/data/df_Strong_AssociationScores_new_wgs.RData")
load("/data/df_Strong_AssociationScores_adjusted_wgs.RData")

corr_result <- data.frame(
  pair = character(),
  corr_new_adj = numeric(),
  pval_new_adj = numeric(),
  stringsAsFactors = FALSE
)

# Get all *_new data frames
new_dfs <- ls(pattern = "_new$")

# Helper function for safe correlation test
safe_corr <- function(x, y) {
  res <- try(cor.test(x, y, method = "spearman", exact = TRUE), silent = TRUE)
  if (inherits(res, "try-error")) return(c(NA, NA))
  return(c(res$estimate, res$p.value))
}

for (new_name in new_dfs) {
  adj_name <- sub("_new$", "_adjusted", new_name)
  
  # Check if corresponding adjusted version exists
  if (!exists(adj_name)) next
  
  df_new <- get(new_name)
  df_adj <- get(adj_name)
  
  # Ensure "score" column exists in both
  if (!("score" %in% colnames(df_new)) || !("score" %in% colnames(df_adj))) {
    warning(paste("Missing 'score' column in one of:", new_name, adj_name))
    next
  }
  
  # Compute correlation: new vs adjusted
  c2 <- safe_corr(df_new$score, df_adj$score)
  
  pair_name <- gsub("^df_association_|_wgs.*$", "", new_name)
  
  corr_result <- rbind(corr_result, data.frame(
    pair = pair_name,
    corr_new_adj = c2[1],
    pval_new_adj = c2[2],
    stringsAsFactors = FALSE
  ))
}

rownames(corr_result) <- corr_result$pair
corr_result$pair <- NULL


# library(xlsx)
# write.xlsx(corr_result,"G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\original_vs_adjusted_AssociationScores\\AdjustedVsNew_StrongAssociation_wgs_Corr.xlsx")





















