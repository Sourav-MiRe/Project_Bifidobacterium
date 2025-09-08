library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

load("G:\\My Drive\\Bifido_Project\\BioRxiv_submission\\GitHub_folders\\DiseaseAssociationPatterns\\DATA\\DiseaseAssociation_InputData.RData")

Metadata_disease_15437_filtered <- Metadata_disease_15437[!Metadata_disease_15437$study_condition == 'HIV',]
Nonbif_SpeciesProfile_disease_15437_filtered <- Nonbif_SpeciesProfile_disease_15437[rownames(Nonbif_SpeciesProfile_disease_15437) %in% rownames(Metadata_disease_15437_filtered),]
Bif_SpeciesProfile_disease_15437_filtered <- Bif_SpeciesProfile_disease_15437[rownames(Bif_SpeciesProfile_disease_15437) %in% rownames(Metadata_disease_15437_filtered),]


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
        Positive = apply(as.matrix(dir_df), 2, function(x) sum(x %in% c(1, 2), na.rm = TRUE)),
        Negative = apply(as.matrix(dir_df), 2, function(x) sum(x %in% c(-1, -2), na.rm = TRUE)),
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

diseases <- unique(sub(".*__", "", colnames(DiseaseControl_Association_df_filtered_NEW)))

new_order <- unlist(lapply(diseases, function(d) {
  paste0(species_order, "__", d)
}))

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

pdf("G:/My Drive/Bifido_Project/DiseaseWise_AssociationScores/DiseaseControl_Association_heatmap_NEW.pdf", 
    height = 20, width = 16)
draw(ht)
dev.off()

###### CARPET ####

ht_drawn <- draw(ht)

row_order <- row_order(ht_drawn)
col_order <- column_order(ht_drawn)

carpet_df <- mat[row_order, col_order]

carpet_DiseaseAssociation_Heatmap <- as.data.frame(carpet_df)

save(carpet_DiseaseAssociation_Heatmap, file = "G:/My Drive/Bifido_Project/DiseaseWise_AssociationScores/carpet_DiseaseAssociation_Heatmap.RData")


#### Re-order the columns acoording to diseases #####

disease_order <- c("Control", "CRC", "Polyps", "Prediabetes", "T2D", 
                   "CVD", "IBD", "Covid", "Parkinsons", "Liver_Disease")

disease <- sub(".*__", "", names(carpet_DiseaseAssociation_Heatmap))

new_col_order <- order(match(disease, disease_order))

carpet_DiseaseAssociation_Heatmap_reordered <- carpet_DiseaseAssociation_Heatmap[, new_col_order]

save(carpet_DiseaseAssociation_Heatmap_reordered, file = "G:/My Drive/Bifido_Project/DiseaseWise_AssociationScores/carpet_DiseaseAssociation_Heatmap_reordered.RData")



#################################################################
### Sequence-Type wise Analysis ###

disease_seqType_summary_df <- Metadata_disease_15437 %>%
  group_by(study_condition, seq_type, study_name) %>%
  summarise(samples = n(), .groups = "drop")
disease_seqType_summary_df <- as.data.frame(disease_seqType_summary_df)


disease_seqType_summary_df_2 <- Metadata_disease_15437 %>%
  group_by(study_condition, seq_type) %>%
  summarise(studies = n_distinct(study_name), .groups = "drop") %>%
  as.data.frame()

############## WGS Analysis ###############

metadata_diseases3_wgs <- Metadata_disease_15437[Metadata_disease_15437$study_condition %in% c("Control","CVD","CRC","IBD_GutInflammation") & Metadata_disease_15437$seq_type == 'WGS',]
Nonbif_SpeciesProfile_disease3_wgs <- Nonbif_SpeciesProfile_disease_15437[rownames(Nonbif_SpeciesProfile_disease_15437) %in% rownames(metadata_diseases3_wgs),]
Bif_SpeciesProfile_disease3_wgs <- Bif_SpeciesProfile_disease_15437[rownames(Bif_SpeciesProfile_disease_15437) %in% rownames(metadata_diseases3_wgs),]


wgs_ControlDisease3_AssociationScores <- get_corr_lists_final(Nonbif_SpeciesProfile_disease3_wgs, Bif_SpeciesProfile_disease3_wgs)

all_assoc_scores_wgs <- list()

diseases_wgs <- unique(unlist(lapply(wgs_ControlDisease3_AssociationScores, names)))
bifs_wgs <- names(wgs_ControlDisease3_AssociationScores)

for (diseases in diseases_wgs) {
  for (bifs in bifs_wgs) {
    assoc_df <- wgs_ControlDisease3_AssociationScores[[bifs]][[diseases]]$Association
    
    scores <- assoc_df$scores
    colname <- paste0(bifs, "__", diseases)
    
    all_assoc_scores_wgs[[colname]] <- scores
  }
}

nonbif_bif_3diseases_control_association_wgs <- do.call(cbind, all_assoc_scores_wgs)
rownames(nonbif_bif_3diseases_control_association_wgs) <- rownames(wgs_ControlDisease3_AssociationScores[[1]][[1]]$Association)
nonbif_bif_3diseases_control_association_wgs <- as.data.frame(nonbif_bif_3diseases_control_association_wgs)

##################### 16S Analysis #####################

metadata_diseases3_16s <- Metadata_disease_15437[Metadata_disease_15437$study_condition %in% c("Control","CVD","CRC","IBD_GutInflammation") & Metadata_disease_15437$seq_type == '16s',]
Nonbif_SpeciesProfile_disease3_16s <- Nonbif_SpeciesProfile_disease_15437[rownames(Nonbif_SpeciesProfile_disease_15437) %in% rownames(metadata_diseases3_16s),]
Bif_SpeciesProfile_disease3_16s <- Bif_SpeciesProfile_disease_15437[rownames(Bif_SpeciesProfile_disease_15437) %in% rownames(metadata_diseases3_16s),]


ControlDisease3_AssociationScores_16s <- get_corr_lists_final(Nonbif_SpeciesProfile_disease3_16s, Bif_SpeciesProfile_disease3_16s)

all_assoc_scores_16s <- list()

diseases_16s <- unique(unlist(lapply(ControlDisease3_AssociationScores_16s, names)))
bifs_16s <- names(ControlDisease3_AssociationScores_16s)

for (diseases in diseases_16s) {
  for (bifs in bifs_16s) {
    assoc_df <- ControlDisease3_AssociationScores_16s[[bifs]][[diseases]]$Association
    
    scores <- assoc_df$scores
    colname <- paste0(bifs, "__", diseases)
    
    all_assoc_scores_16s[[colname]] <- scores
  }
}

nonbif_bif_3diseases_control_association_16s <- do.call(cbind, all_assoc_scores_16s)
rownames(nonbif_bif_3diseases_control_association_16s) <- rownames(ControlDisease3_AssociationScores_16s[[1]][[1]]$Association)
nonbif_bif_3diseases_control_association_16s <- as.data.frame(nonbif_bif_3diseases_control_association_16s)

#####################################################################
### 16S WGS Combined Heatmap ###

nonbif_bif_3diseases_control_association_wgs <- nonbif_bif_3diseases_control_association_wgs[,names(nonbif_bif_3diseases_control_association_16s)]
all(names(nonbif_bif_3diseases_control_association_wgs) == names(nonbif_bif_3diseases_control_association_16s))
all(rownames(nonbif_bif_3diseases_control_association_wgs) == rownames(nonbif_bif_3diseases_control_association_16s))

names(nonbif_bif_3diseases_control_association_wgs) <- paste0(names(nonbif_bif_3diseases_control_association_wgs),"__wgs")
names(nonbif_bif_3diseases_control_association_16s) <- paste0(names(nonbif_bif_3diseases_control_association_16s),"__16s")

names(nonbif_bif_3diseases_control_association_wgs) <- gsub("IBD_GutInflammation","IBD",names(nonbif_bif_3diseases_control_association_wgs))
names(nonbif_bif_3diseases_control_association_16s) <- gsub("IBD_GutInflammation","IBD",names(nonbif_bif_3diseases_control_association_16s))

names(nonbif_bif_3diseases_control_association_wgs) <- gsub("Bifidobacterium_","",names(nonbif_bif_3diseases_control_association_wgs))
names(nonbif_bif_3diseases_control_association_16s) <- gsub("Bifidobacterium_","",names(nonbif_bif_3diseases_control_association_16s))


species_order <- c("detection", "longum", "adolescentis", 
                   "pseudocatenulatum", "dentium", 
                   "bifidum", "breve", "catenulatum", "animalis")

diseases <- c("Control", "CRC", "CVD", "IBD")

ordered_colnames <- function(suffix) {
  unlist(lapply(diseases, function(d) {
    paste0(species_order, "__", d, "__", suffix)
  }))
}

nonbif_bif_3diseases_control_association_16s <- 
  nonbif_bif_3diseases_control_association_16s[, 
                                               ordered_colnames("16s")]

nonbif_bif_3diseases_control_association_wgs <- 
  nonbif_bif_3diseases_control_association_wgs[, 
                                               ordered_colnames("wgs")]


groups <- c("Control", "CRC", "CVD", "IBD")


combined_list <- list()


for (g in groups) {
  cols_16s <- grep(g, colnames(nonbif_bif_3diseases_control_association_16s), value = TRUE)
  cols_wgs <- grep(g, colnames(nonbif_bif_3diseases_control_association_wgs), value = TRUE)
  
  combined_list[[g]] <- cbind(
    nonbif_bif_3diseases_control_association_16s[, cols_16s, drop = FALSE],
    nonbif_bif_3diseases_control_association_wgs[, cols_wgs, drop = FALSE]
  )
}

DiseaseAssociation_16s_wgs_df <- do.call(cbind, combined_list)

names(DiseaseAssociation_16s_wgs_df) <- gsub("^(Control|CRC|CVD|IBD)\\.", "", names(DiseaseAssociation_16s_wgs_df))

DiseaseAssociation_16s_wgs_df_filtered <- DiseaseAssociation_16s_wgs_df[!rowSums(DiseaseAssociation_16s_wgs_df) == 0,]

###### Filtration ######## (For more than one disease group a non-bif shows abs(score)>=0.25)

library(stringr)

df <- DiseaseAssociation_16s_wgs_df_filtered  

col_info <- do.call(rbind, strsplit(colnames(df), "__"))
colnames(col_info) <- c("Bif", "Disease", "SeqType")

diseases <- unique(col_info[, "Disease"])


check_row <- function(row_vals) {
  
  disease_hits <- sapply(diseases, function(d) {
    
    vals_16s <- row_vals[col_info[, "Disease"] == d & col_info[, "SeqType"] == "16s"]
    vals_wgs <- row_vals[col_info[, "Disease"] == d & col_info[, "SeqType"] == "wgs"]
    
    
    any(abs(vals_16s) >= 0.25) && any(abs(vals_wgs) >= 0.25)
  })
  
  sum(disease_hits) >= 1
}

keep_rows <- apply(df, 1, check_row)

DiseaseAssociation_filtered <- df[keep_rows, ]

DiseaseAssociation_16s_wgs_df_final <- DiseaseAssociation_filtered


### Heatmap ###

library(ComplexHeatmap)
library(circlize)
library(grid)

col_fun <- colorRamp2(
  c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1),
  c("#8B0053", "#CC3377", "#F6B6D2", "#FFFFFF", "#CDE9B6", "#77C679","#004D00")
)

mat <- as.matrix(DiseaseAssociation_16s_wgs_df_final)

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
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
    labels = c("-1", "-0.75", "-0.5", "-0.25", "0", 
               "0.25", "0.5", "0.75", "1"),
    title = "Association Scores"
  ),
  width = unit(ncol(mat) * 3.3, "mm"),
  height = unit(nrow(mat) * 3.3, "mm"),
  
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width = width, height = height, 
              gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
    
  }
)

pdf("G:/My Drive/Bifido_Project/DiseaseWise_AssociationScores/DiseaseControl_Association_16s_wgs_heatmap.pdf", 
    height = 14, width = 13)
draw(ht)
dev.off()









