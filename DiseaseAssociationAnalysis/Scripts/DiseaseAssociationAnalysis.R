library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

load("/data/DiseaseAssociation_InputData.RData")

Metadata_disease_15437_filtered <- Metadata_disease_15437[!Metadata_disease_15437$study_condition %in% c('HIV','Liver_Disease','Prediabetes','Polyps'),]
Nonbif_SpeciesProfile_disease_15437_filtered <- Nonbif_SpeciesProfile_disease_15437[rownames(Nonbif_SpeciesProfile_disease_15437) %in% rownames(Metadata_disease_15437_filtered),]
Bif_SpeciesProfile_disease_15437_filtered <- Bif_SpeciesProfile_disease_15437[rownames(Bif_SpeciesProfile_disease_15437) %in% rownames(Metadata_disease_15437_filtered),]

disease_seqType_summary_df <- Metadata_disease_15437_filtered %>%
  group_by(study_condition, seq_type, study_name) %>%
  summarise(samples = n(), .groups = "drop")
disease_seqType_summary_df <- as.data.frame(disease_seqType_summary_df)


disease_seqType_summary_df_2 <- Metadata_disease_15437_filtered %>%
  group_by(study_condition, seq_type) %>%
  summarise(studies = n_distinct(study_name), .groups = "drop") %>%
  as.data.frame()

disease_seqType_summary_df_3 <- Metadata_disease_15437_filtered %>%
  group_by(study_condition) %>%
  summarise(studies = n_distinct(study_name), .groups = "drop") %>%
  as.data.frame()


get_corr_lists_final <- function(data, bif_SpeciesProfile_filtered) {
  
  # Extract metadata
  study_names <- data$study_name
  study_conditions <- data$study_condition
  
  # Remove metadata from `data` (keep only non-bifs)
  non_bifs <- data[, !(colnames(data) %in% c("study_name", "study_condition")), drop = FALSE]
  
  # Identify bif taxa
  bif_taxa <- colnames(bif_SpeciesProfile_filtered)
  
  # Unique diseases
  diseases <- unique(study_conditions)
  
  # Initialize final output
  bif_results <- vector("list", length(bif_taxa))
  names(bif_results) <- bif_taxa
  
  # Loop over bif taxa
  for (bif in bif_taxa) {
    
    disease_list <- list()
    
    # Loop over diseases
    for (disease in diseases) {
      
      # Subset samples for this disease
      idx_disease <- which(study_conditions == disease)
      non_bifs_disease <- non_bifs[idx_disease, , drop = FALSE]
      bif_disease <- bif_SpeciesProfile_filtered[idx_disease, bif, drop = FALSE]
      study_names_disease <- study_names[idx_disease]
      
      # Get unique studies
      unique_studies <- unique(study_names_disease)
      
      # Create matrices (always 2D)
      corr_mat <- matrix(NA, nrow = length(unique_studies), ncol = ncol(non_bifs_disease),
                         dimnames = list(unique_studies, colnames(non_bifs_disease)))
      pval_mat <- corr_mat
      dir_mat  <- corr_mat
      
      # Studywise correlation
      for (study in unique_studies) {
        idx_study <- which(study_names_disease == study)
        non_bifs_study <- non_bifs_disease[idx_study, , drop = FALSE]
        bif_study <- bif_disease[idx_study, , drop = FALSE]
        
        for (nb in colnames(non_bifs_study)) {
          if (length(unique(non_bifs_study[[nb]])) > 1 && length(unique(bif_study[[1]])) > 1) {
            ct <- suppressWarnings(cor.test(non_bifs_study[[nb]], bif_study[[1]], method = "spearman"))
            corr_mat[study, nb] <- ct$estimate
            pval_mat[study, nb]  <- ct$p.value
            
            # Direction encoding
            if (ct$estimate > 0 && ct$p.value <= 0.05) {
              dir_mat[study, nb] <- 2
            } else if (ct$estimate > 0 && ct$p.value <= 0.1) {
              dir_mat[study, nb] <- 1
            } else if (ct$estimate < 0 && ct$p.value <= 0.05) {
              dir_mat[study, nb] <- -2
            } else if (ct$estimate < 0 && ct$p.value <= 0.1) {
              dir_mat[study, nb] <- -1
            } else {
              dir_mat[study, nb] <- 0
            }
          } else {
            corr_mat[study, nb] <- NA
            pval_mat[study, nb]  <- NA
            dir_mat[study, nb]   <- 0
          }
        }
      }
      
      # Convert to data.frames
      corr_df <- as.data.frame(corr_mat, stringsAsFactors = FALSE)
      pval_df <- as.data.frame(pval_mat, stringsAsFactors = FALSE)
      dir_df  <- as.data.frame(dir_mat, stringsAsFactors = FALSE)
      
      # Association summary (use apply to handle 1-col edge case)
      association_df <- data.frame(
        Positive = apply(as.matrix(dir_df), 2, function(x) sum(x == 2, na.rm = TRUE)),
        Negative = apply(as.matrix(dir_df), 2, function(x) sum(x == -2, na.rm = TRUE)),
        Total    = nrow(dir_df)
      )
      rownames(association_df) <- colnames(dir_df)
      
      # Add scores column
      association_df$scores <- apply(association_df, 1, function(x) {
        score1 <- (x[1] - x[2]) / x[3]
        score2 <- 1 - ((min(x[1:2]) + 0.0001) / (max(x[1:2]) + 0.0001))
        score1 * score2
      })
      
      # Save four dfs for this disease
      disease_list[[disease]] <- list(
        Corr = corr_df,
        Pval = pval_df,
        Dir = dir_df,
        Association = association_df
      )
    }
    
    # Save disease list under bif
    bif_results[[bif]] <- disease_list
  }
  
  return(bif_results)
}

ControlDiseaseAssociation_results_NEW <- get_corr_lists_final(Nonbif_SpeciesProfile_disease_15437_filtered, Bif_SpeciesProfile_disease_15437_filtered)

#save(ControlDiseaseAssociation_results_NEW, file = "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\Fig4_ReComputation\\ControlDiseaseAssociation_results_Latest.RData")

all_assoc_scores <- list()

diseases <- unique(unlist(lapply(ControlDiseaseAssociation_results_NEW, names)))
bifs <- names(ControlDiseaseAssociation_results_NEW)

for (disease in diseases) {
  for (bif in bifs) {
    assoc_df <- ControlDiseaseAssociation_results_NEW[[bif]][[disease]]$Association
    
    scores <- assoc_df$scores
    colname <- paste0(bif, "__", disease)
    
    all_assoc_scores[[colname]] <- scores
  }
}

nonbif_bif_disease_control_association_df_NEW <- do.call(cbind, all_assoc_scores)
rownames(nonbif_bif_disease_control_association_df_NEW) <- rownames(ControlDiseaseAssociation_results_NEW[[1]][[1]]$Association)
nonbif_bif_disease_control_association_df_NEW <- as.data.frame(nonbif_bif_disease_control_association_df_NEW)

nonbif_bif_disease_control_association_df_NEW <- nonbif_bif_disease_control_association_df_NEW[!rowSums(nonbif_bif_disease_control_association_df_NEW) == 0,]

#write.xlsx(nonbif_bif_disease_control_association_df_NEW,"G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\Fig4_ReComputation\\DiseaseAssociationScores_AllSpecies.xlsx")

###### Filtration ######## (For more than one disease a non-bif shows abs(score)>=0.25)

disease_names <- sapply(strsplit(colnames(nonbif_bif_disease_control_association_df_NEW), "__"), `[`, 2)

keep_rows <- apply(nonbif_bif_disease_control_association_df_NEW, 1, function(x) {
  
  disease_hits <- tapply(x, disease_names, function(vals) any(abs(vals) >= 0.25, na.rm = TRUE))
  
  sum(disease_hits, na.rm = TRUE) > 2
})

DiseaseControl_Association_df_filtered_NEW <- nonbif_bif_disease_control_association_df_NEW[keep_rows, ]

names(DiseaseControl_Association_df_filtered_NEW) <- sub("IBD_GutInflammation","IBD", names(DiseaseControl_Association_df_filtered_NEW))
names(DiseaseControl_Association_df_filtered_NEW) <- sub("Bifidobacterium_","", names(DiseaseControl_Association_df_filtered_NEW))

species_order <- c(
  "detection", "longum", "adolescentis", "pseudocatenulatum",
  "dentium", "bifidum", "breve", "catenulatum", "animalis"
)

disease_order <- c("Control", "CRC", "Polyps", "Prediabetes", "T2D", 
                   "CVD", "IBD", "Covid", "Parkinsons", "Liver_Disease")

current_cols <- colnames(DiseaseControl_Association_df_filtered_NEW)

split_names <- strsplit(current_cols, "__")
species_vec  <- sapply(split_names, `[`, 1)
disease_vec  <- sapply(split_names, `[`, 2)

disease_factor <- factor(disease_vec, levels = disease_order)
species_factor <- factor(species_vec, levels = species_order)

new_order <- order(disease_factor, species_factor, na.last = TRUE)

DiseaseControl_Association_df_filtered_NEW <- DiseaseControl_Association_df_filtered_NEW[, new_order]

### Heatmap ####

col_fun <- colorRamp2(
  c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
  c("#8B0053", "#CC3377", "#F6B6D2", "#FFFFFF", "#CDE9B6", "#77C679", "#004D00")
)

# Convert to numeric matrix
mat <- as.matrix(DiseaseControl_Association_df_filtered_NEW)

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

pdf("/results/DiseaseControl_Association_heatmap.pdf", 
    height = 20, width = 16)
draw(ht)
dev.off()

###### CARPET ####

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_DiseaseAssociation_Heatmap <- as.data.frame(carpet_df)

#write.xlsx(carpet_DiseaseAssociation_Heatmap, "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\Fig4_ReComputation\\carpet_DiseaseAssociation_Heatmap.xlsx")


###################################################################
###################################################################
######## Shannon Association Scores across Diseases ###########

Bif_SpeciesProfile_disease_15437$Bifidobacterium_detection <- apply(Bif_SpeciesProfile_disease_15437, 1, function(x) length(x[x >= 0.0001]))
disease_shannon <- diversity(Nonbif_SpeciesProfile_disease_15437[,1:355], index = "shannon")

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
disease_pielou <- calculate_pielou(Nonbif_SpeciesProfile_disease_15437[,1:355])

disease_shannon_pielou_df <- data.frame(shannon = disease_shannon, pielou = disease_pielou)

combined_df <- as.data.frame(cbind(Nonbif_SpeciesProfile_disease_15437[,1:355],Bif_SpeciesProfile_disease_15437,disease_shannon_pielou_df,Nonbif_SpeciesProfile_disease_15437[,356:357]))

combined_df_filtered <- combined_df[!combined_df$study_condition == 'HIV',]

# Create list to store results for each age category
disease_list_shannon <- list()

# Loop through each unique age category
for (disease in unique(combined_df_filtered$study_condition)) {
  
  temp_df <- subset(combined_df_filtered, study_condition %in% disease)
  bifido_species <- colnames(temp_df)[356:364]
  studies <- unique(temp_df$study_name)
  
  corr_df <- matrix(NA, nrow = length(studies), ncol = length(bifido_species),
                    dimnames = list(studies, bifido_species))
  pval_df <- corr_df
  
  # Loop through each study
  for (study in studies) {
    study_df <- subset(temp_df, study_name == study)
    
    for (sp in bifido_species) {
      if (length(unique(study_df[[sp]])) > 1 && length(unique(study_df$shannon)) > 1) {
        res <- suppressWarnings(cor.test(study_df[[sp]], study_df$shannon,
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
  disease_list_shannon[[disease]] <- list(
    corr_df = corr_df,
    pval_df = pval_df,
    dir = dir,
    association_df = association_df
  )
}

################ making the Shannon df ###########

disease_order <- c("Control", "CRC", "Polyps", "Prediabetes", "T2D",
                   "CVD", "IBD_GutInflammation", "Covid", "Parkinsons", "Liver_Disease")

bif_order <- c("detection", "longum", "adolescentis", "pseudocatenulatum",
               "dentium", "bifidum", "breve", "catenulatum", "animalis")

score_vector <- c()

for (disease in disease_order) {
  # Extract the association_df for the disease
  assoc_df <- disease_list_shannon[[disease]]$association_df
  
  # Match the bifidobacteria order (row names assumed to be like "Bifidobacterium_xxx")
  bif_scores <- sapply(bif_order, function(bif) {
    row_name <- paste0("Bifidobacterium_", bif)
    if (row_name %in% rownames(assoc_df)) {
      return(assoc_df[row_name, "Score"])
    } else {
      return(NA)  # if missing
    }
  })
  
  # Name the vector elements as e.g., Control_detection, Control_longum, etc.
  names(bif_scores) <- paste0(disease, "_", bif_order)
  
  # Append to main vector
  score_vector <- c(score_vector, bif_scores)
}

shannon_df <- as.data.frame(t(score_vector))
rownames(shannon_df) <- "shannon"

###################################################################
########## Calculating Scores for Pielou ##########

# Create list to store results for each age category
disease_list_pielou <- list()

# Loop through each unique age category
for (disease in unique(combined_df_filtered$study_condition)) {
  
  temp_df <- subset(combined_df_filtered, study_condition %in% disease)
  bifido_species <- colnames(temp_df)[356:364]
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
  disease_list_pielou[[disease]] <- list(
    corr_df = corr_df,
    pval_df = pval_df,
    dir = dir,
    association_df = association_df
  )
}

################ making the Shannon df ###########

disease_order <- c("Control", "CRC", "Polyps", "Prediabetes", "T2D",
                   "CVD", "IBD_GutInflammation", "Covid", "Parkinsons", "Liver_Disease")

bif_order <- c("detection", "longum", "adolescentis", "pseudocatenulatum",
               "dentium", "bifidum", "breve", "catenulatum", "animalis")

score_vector <- c()

for (disease in disease_order) {
  # Extract the association_df for the disease
  assoc_df <- disease_list_pielou[[disease]]$association_df
  
  # Match the bifidobacteria order (row names assumed to be like "Bifidobacterium_xxx")
  bif_scores <- sapply(bif_order, function(bif) {
    row_name <- paste0("Bifidobacterium_", bif)
    if (row_name %in% rownames(assoc_df)) {
      return(assoc_df[row_name, "Score"])
    } else {
      return(NA)  # if missing
    }
  })
  
  # Name the vector elements as e.g., Control_detection, Control_longum, etc.
  names(bif_scores) <- paste0(disease, "_", bif_order)
  
  # Append to main vector
  score_vector <- c(score_vector, bif_scores)
}

pielou_df <- as.data.frame(t(score_vector))
rownames(pielou_df) <- "pielou"

load("/data/carpet_DiseaseAssociation_Heatmap_reordered.RData")

names(pielou_df) = names(shannon_df) = names(carpet_DiseaseAssociation_Heatmap_reordered)

combined_shannon_pielou_disease_df <- as.data.frame(rbind(carpet_DiseaseAssociation_Heatmap_reordered, shannon_df, pielou_df))



load("/data/carpet_DiseaseAssociation_Heatmap_expanded.RData")

species_order <- c("detection", "longum", "adolescentis", "pseudocatenulatum", "bifidum", "dentium", "breve", "catenulatum", "animalis")

disease_order <- c("Control", "CRC", "Polyps", "Prediabetes", "T2D", 
                   "CVD", "IBD", "Covid", "Parkinsons", "Liver_Disease")

current_cols <- colnames(carpet_DiseaseAssociation_Heatmap_expanded)

split_names <- strsplit(current_cols, "__")
species_vec  <- sapply(split_names, `[`, 1)
disease_vec  <- sapply(split_names, `[`, 2)

disease_factor <- factor(disease_vec, levels = disease_order)
species_factor <- factor(species_vec, levels = species_order)

new_order <- order(disease_factor, species_factor, na.last = TRUE)

combined_shannon_pielou_disease_df <- carpet_DiseaseAssociation_Heatmap_expanded[, new_order]

#####################################################################

load("/data/carpet_DiseaseAssociation_Heatmap_expanded.RData")

carpet_DiseaseAssociation_Heatmap_expanded <- carpet_DiseaseAssociation_Heatmap_expanded[,grep("Liver_Disease|Polyps|Prediabetes",colnames(carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE,invert=TRUE)]

df_association_patterns_diseases <- data.frame("Positive"=apply(carpet_DiseaseAssociation_Heatmap_expanded,1,function(x)(length(x[x>=0.20]))),"Negative"=apply(carpet_DiseaseAssociation_Heatmap_expanded,1,function(x)(length(x[x<=(-0.20)]))))

ordered_species <- rownames(df_association_patterns_diseases)[rev(order(df_association_patterns_diseases[,1]-df_association_patterns_diseases[,2]))]

ordered_carpet_DiseaseAssociation_Heatmap_expanded <- carpet_DiseaseAssociation_Heatmap_expanded[ordered_species,]


subject_groups <- c("Control","CRC","T2D","CVD","IBD","Covid","Parkinsons")

thresholds_for_association_scores <- as.data.frame(matrix(NA,7,2))
rownames(thresholds_for_association_scores) <- c("Control","CRC","T2D","CVD","IBD","Covid","Parkinsons")
colnames(thresholds_for_association_scores) <- c("Lower","Higher")

for(i in 1:7)
{
  group <- subject_groups[i]
  thresholds_for_association_scores[i,1] <- as.numeric(quantile(unlist(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(group,colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)]),prob=c(0,0.1,0.90,1))[2])
  thresholds_for_association_scores[i,2] <- as.numeric(quantile(unlist(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(group,colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)]),prob=c(0,0.1,0.90,1))[3])
}

grade_associations <- function(x,group)
{
  y <- ifelse((x>=thresholds_for_association_scores[group,2]),2,ifelse((x<=thresholds_for_association_scores[group,1]),-2,sign(x)))
  return(y)
}

dig_ordered_carpet_groups <- as.data.frame(cbind(apply(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(subject_groups[1],colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)],2,function(x)(grade_associations(x,subject_groups[1]))),apply(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(subject_groups[2],colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)],2,function(x)(grade_associations(x,subject_groups[2]))),apply(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(subject_groups[3],colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)],2,function(x)(grade_associations(x,subject_groups[3]))),apply(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(subject_groups[4],colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)],2,function(x)(grade_associations(x,subject_groups[4]))),apply(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(subject_groups[5],colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)],2,function(x)(grade_associations(x,subject_groups[5]))),apply(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(subject_groups[6],colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)],2,function(x)(grade_associations(x,subject_groups[6]))),apply(ordered_carpet_DiseaseAssociation_Heatmap_expanded[,grep(subject_groups[7],colnames(ordered_carpet_DiseaseAssociation_Heatmap_expanded),value=TRUE)],2,function(x)(grade_associations(x,subject_groups[7])))))

df_overall_pattern <- data.frame("positive"=apply(dig_ordered_carpet_groups,1,function(x)(length(x[x==2]))),"negative"=apply(dig_ordered_carpet_groups,1,function(x)(length(x[x==(-2)]))))

####################################################################
####################################################################

load("/data/251227_Disease_Association_Analysis.RData")

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

###### Heat map #####

mat <- as.matrix(dig_ordered_carpet_groups[,1:63])

col_fun <- colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("royalblue4", "skyblue", "gray95","pink","deeppink4")
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

pdf("/results/Dig_DiseaseControl_StrongAssociation_ordered.pdf", 
    height = 20, width = 16)
draw(ht)
dev.off()

##### Carpet #####

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_DiseaseAssociation_Heatmap <- as.data.frame(carpet_df)

#write.xlsx(carpet_DiseaseAssociation_Heatmap, "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\Fig4_ReComputation\\Dig_carpet_DiseaseAssociation_Strong_ordered.xlsx")

#####  Line Plot #####

library(tidyr)
library(dplyr)
library(ggplot2)

association_long <- df_overall_pattern %>%
  mutate(Species = rownames(.)) %>%
  pivot_longer(cols = c("positive", "negative"),
               names_to = "Association",
               values_to = "Count")

# Retain rowname order
association_long$Species <- factor(
  association_long$Species,
  levels = rownames(df_overall_pattern)
)

pdf("/results/dig_disease_association_lineplot.pdf",
    height = 1, width = 10)

ggplot(association_long, aes(x = Species, y = Count, group = Association, color = Association)) +
  geom_line(linewidth = 0.3) +
  geom_point(shape = 16, size = 0.7) +
  scale_color_manual(values = c("positive" = "red", "negative" = "blue")) +
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

#####################################################################
########### Re-ordering ############

cn <- colnames(carpet_DiseaseAssociation_Heatmap)

# Split into feature and disease
tmp <- do.call(rbind, strsplit(cn, "__"))
df <- data.frame(
  colname = cn,
  feature = tmp[, 1],
  disease = tmp[, 2],
  stringsAsFactors = FALSE
)

# Desired orders
disease_order <- c(
  "Control", "CRC", "T2D", "CVD", "IBD", "Covid", "Parkinsons"
)

feature_order <- c(
  "detection",
  "longum",
  "adolescentis",
  "pseudocatenulatum",
  "bifidum",
  "dentium",
  "breve",
  "catenulatum",
  "animalis"
)

df$feature <- factor(df$feature, levels = feature_order)
df$disease <- factor(df$disease, levels = disease_order)

# Order first by feature, then by disease
df <- df[order(df$feature, df$disease), ]

# Reorder the matrix/dataframe
carpet_DiseaseAssociation_Heatmap_new <-
  carpet_DiseaseAssociation_Heatmap[, df$colname]

mat <- as.matrix(carpet_DiseaseAssociation_Heatmap_new)

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

pdf("/results/Dig_DiseaseControl_StrongAssociation_ordered_NEW.pdf", 
    height = 20, width = 16)
draw(ht)
dev.off()

##### Carpet #####

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_DiseaseAssociation_Heatmap <- as.data.frame(carpet_df)

#write.xlsx(carpet_DiseaseAssociation_Heatmap, "G:\\My Drive\\Bifido_Project\\NatComm_SubmissionV1\\MajorRevision1\\Fig4_ReComputation\\Dig_carpet_DiseaseAssociation_Strong_ordered_NEW.xlsx")

#####  Line Plot #####

library(tidyr)
library(dplyr)
library(ggplot2)

association_long <- df_overall_pattern %>%
  mutate(Species = rownames(.)) %>%
  pivot_longer(cols = c("positive", "negative"),
               names_to = "Association",
               values_to = "Count")

# Retain rowname order
association_long$Species <- factor(
  association_long$Species,
  levels = rownames(df_overall_pattern)
)

pdf("/results/dig_disease_association_lineplot.pdf",
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





