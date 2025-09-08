#Figure S13
#PCoA for seq type

# Load required libraries
library(tidyverse)
library(vegan)
library(ggplot2)

# Load the RData file
load("OriginalAssociationScores_Modified.RData")

make_seqtype_pcoa_plot <- function(species, df_list) {
  # Extract adult/senior featureImportance for both 16S and WGS
  adult_16s <- df_list[[paste0("df_association_adult_", species, "_16s")]][["featureImportance"]]
  senior_16s <- df_list[[paste0("df_association_senior_", species, "_16s")]][["featureImportance"]]
  adult_wgs <- df_list[[paste0("df_association_adult_", species, "_wgs")]][["featureImportance"]]
  senior_wgs <- df_list[[paste0("df_association_senior_", species, "_wgs")]][["featureImportance"]]
  
  # Check column consistency
  if (!identical(colnames(adult_16s), colnames(senior_16s)) || !identical(colnames(adult_wgs), colnames(senior_wgs))) {
    stop("Column names don't match between adult and senior matrices for either 16S or WGS")
  }
  
  # Define a helper function to merge adult/senior using row-means logic
  merge_studies <- function(adult_df, senior_df) {
    species_cols <- colnames(adult_df)
    all_studies <- union(rownames(adult_df), rownames(senior_df))
    merged_df <- matrix(NA, nrow = length(all_studies), ncol = length(species_cols),
                        dimnames = list(all_studies, species_cols))
    for (study in all_studies) {
      in_adult <- study %in% rownames(adult_df)
      in_senior <- study %in% rownames(senior_df)
      if (in_adult && !in_senior) {
        merged_df[study, ] <- as.numeric(adult_df[study, ])
      } else if (!in_adult && in_senior) {
        merged_df[study, ] <- as.numeric(senior_df[study, ])
      } else if (in_adult && in_senior) {
        merged_df[study, ] <- as.numeric((adult_df[study, ] + senior_df[study, ]) / 2)
      }
    }
    return(as.data.frame(merged_df))
  }
  
  # Merge 16S and WGS separately
  merged_16s <- merge_studies(adult_16s, senior_16s)
  merged_wgs <- merge_studies(adult_wgs, senior_wgs)
  
  # Add seq_type labels to rownames
  rownames(merged_16s) <- paste0(rownames(merged_16s), "_16S")
  rownames(merged_wgs) <- paste0(rownames(merged_wgs), "_WGS")
  
  # Combine and clean
  combined_df <- bind_rows(merged_16s, merged_wgs)
  combined_df[is.na(combined_df)] <- 0
  
  # Rank scaling per row
  scaled_df <- t(apply(combined_df, 1, function(x) {
    if (max(x) != min(x)) {
      (x - min(x)) / (max(x) - min(x))
    } else {
      rep(0, length(x))
    }
  }))
  scaled_df <- as.data.frame(scaled_df)
  
  # PCoA
  dist_mat <- dist(scaled_df)
  pcoa <- cmdscale(dist_mat, k = 2)
  pcoa_df <- as.data.frame(pcoa)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  pcoa_df$SeqType <- ifelse(grepl("_16S$", rownames(pcoa_df)), "16S", "WGS")
  pcoa_df$SeqType <- factor(pcoa_df$SeqType, levels = c("16S", "WGS"))
  
  # PERMANOVA
  permanova <- adonis2(dist_mat ~ SeqType, data = pcoa_df, permutations = 999, method = "euclidean")
  r2 <- formatC(permanova$R2[1], format = "f", digits = 3)
  pval <- formatC(permanova$`Pr(>F)`[1], format = "e", digits = 2)
  permanova_text <- paste0("PERMANOVA\nR² = ", r2, "\np = ", pval)
  
  # Centroids
  centroids <- pcoa_df %>%
    group_by(SeqType) %>%
    summarise(PCoA1_centroid = mean(PCoA1), PCoA2_centroid = mean(PCoA2), .groups = "drop")
  pcoa_df <- left_join(pcoa_df, centroids, by = "SeqType")
  
  # Plot
  plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = SeqType)) +
    geom_point(size = 3) +
    geom_segment(aes(xend = PCoA1_centroid, yend = PCoA2_centroid), linetype = "solid") +
    scale_color_manual(values = c("16S" = "darkorange", "WGS" = "steelblue")) +
    theme_minimal() +
    labs(title = paste("PCoA: 16S vs WGS for", species),
         x = "Principal Coordinate 1",
         y = "Principal Coordinate 2") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 12)
          ,legend.position = "none"
    )
  
  # Save plot
  pdf_name <- paste0("PCoA_", species, "_seqtype.pdf")
  ggsave(pdf_name, plot = plot, width = 8, height = 5)
  cat("Saved plot:", pdf_name, "\n")
  cat("PERMANOVA result for", species, ":\n")
  print(permanova_text)
  
  return(scaled_df)
}

# List of your 8 Bifidobacterium species (use exact column names as in your featureImportance matrices)
bif_species_list <- c("adolescentis", "animalis", "bifidum", "breve", "catenulatum", "dentium", "longum", "pseudocatenulatum", "detection")

scaled_dfs <- list()
for (sp in bif_species_list) {
  cat("Running PCoA for:", sp, "\n")
  scaled_df <- make_seqtype_pcoa_plot(sp, df_list = mget(ls(pattern = "df_association_")))
  scaled_dfs[[sp]] <- scaled_df
}

#____________________________________________________________________________________________________
#____________________________________________________________________________________________________
#Figure S13
#PCoA for age-category

# Load required libraries
library(tidyverse)
library(vegan)
library(ggplot2)

# Load the RData file
load("OriginalAssociationScores_Modified.RData")

# Get all objects from the loaded file
all_dfs <- ls()

# List of species to loop over
species_list <- c("adolescentis", "bifidum", "breve", "catenulatum", "dentium", "longum", "pseudocatenulatum", "detection")

# Define prefixes and group labels
age_groups <- c("infant", "adult", "senior")
prefixes <- c("infant" = "I__", "adult" = "A__", "senior" = "S__")
group_labels <- c("I" = "Infant", "A" = "Adult", "S" = "Senior")

# Function to generate scaled matrix + PCoA matrix and save plot
make_pcoa_plot <- function(species, method, input_dfs) {
  feature_matrices <- list()
  
  for (age in age_groups) {
    df_name <- paste0("df_association_", age, "_", species, "_", method)
    if (df_name %in% input_dfs) {
      df <- get(df_name)
      feature_data <- df[["featureImportance"]]
      rownames(feature_data) <- paste0(prefixes[age], rownames(feature_data))
      feature_matrices[[age]] <- feature_data
    }
  }
  
  combined <- bind_rows(feature_matrices)
  combined[is.na(combined)] <- 0
  
  # Row-wise min-max scaling
  scaled <- t(apply(combined, 1, function(x) {
    if (max(x) != min(x)) {
      (x - min(x)) / (max(x) - min(x))
    } else {
      rep(0, length(x))
    }
  }))
  rownames(scaled) <- rownames(combined)
  scaled_df <- as.data.frame(scaled)
  
  # PCoA using Euclidean distance
  dist_mat <- dist(scaled_df)
  pcoa <- cmdscale(dist_mat, k = 2)
  pcoa_df <- as.data.frame(pcoa)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  pcoa_df$Group <- substr(rownames(pcoa_df), 1, 1)
  pcoa_df$Group <- factor(group_labels[pcoa_df$Group], levels = c("Infant", "Adult", "Senior"))
  
  # PERMANOVA
  permanova <- adonis2(dist_mat ~ Group, data = pcoa_df, permutations = 999, method = 'euclidean')
  r2 <- formatC(permanova$R2[1], format = "f", digits = 3)
  pval <- formatC(permanova$`Pr(>F)`[1], format = "e", digits = 2)
  permanova_text <- paste0("PERMANOVA\nR² = ", r2, "\np = ", pval)
  
  # Centroids for arrow visualization
  centroids <- pcoa_df %>%
    group_by(Group) %>%
    summarise(PCoA1_centroid = mean(PCoA1), PCoA2_centroid = mean(PCoA2), .groups = "drop")
  
  pcoa_df <- left_join(pcoa_df, centroids, by = "Group")
  
  # Plot
  plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 4) +
    geom_segment(aes(xend = PCoA1_centroid, yend = PCoA2_centroid), linetype = "solid") +
    scale_color_manual(values = c("Infant" ="forestgreen", "Adult" = "orangered2", "Senior" = "royalblue3")) +
    #annotate("text", x = text_x, y = text_y, label = permanova_text,
    #hjust = 1, vjust = 1, size = 4.5, fontface = "italic") +
    theme_minimal() +
    labs(
      title = paste("PCoA Plot for", species, method),
      x = "Principal Coordinate 1",
      y = "Principal Coordinate 2",
    ) +
    theme(
      panel.grid = element_blank(),         
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
  
  # Save plot
  pdf_name <- paste0("PCoA_", species, "_", method, ".pdf")
  ggsave(pdf_name, plot = plot, width = 5, height = 5)
  cat("Saved plot:", pdf_name, "\n")
  
  # Print PERMANOVA stats to console
  cat("PERMANOVA result for", species, method, ":\n")
  print(permanova_text)
  
  # Save dataframes
  assign(paste0("df_association_", species, "_", method), scaled_df, envir = .GlobalEnv)
  assign(paste0("pcoa_", species, "_", method), pcoa_df, envir = .GlobalEnv)
}

# Loop through species and method
for (sp in species_list) {
  for (method in c("16s", "wgs")) {
    relevant_dfs <- all_dfs[grepl(paste0("df_association_.*_", sp, "_", method), all_dfs)]
    if (length(relevant_dfs) > 0) {
      make_pcoa_plot(sp, method, relevant_dfs)
    }
  }
}

#____________________________________________________________________________________________________
#____________________________________________________________________________________________________
#____________________________________________________________________________________________________
#Figure S14
#PCoA for cohort lifestyle

# Load required libraries
library(tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)

# Load the RData file
load("OriginalAssociationScores_Modified.RData")
load("cohort_metadata.RData")

# Filter all objects in the environment
df_list <- mget(ls(pattern = "^df_association_"), envir = .GlobalEnv)


make_pcoa_plot <- function(species, method = "wgs", df_list = mget(ls(pattern = "df_association_")), cohort_metadata) {
  # Extract relevant adult and senior dataframes
  adult_df <- df_list[[paste0("df_association_adult_", species, "_", method)]]$featureImportance
  senior_df <- df_list[[paste0("df_association_senior_", species, "_", method)]]$featureImportance
  
  # Check column name consistency
  if (!identical(colnames(adult_df), colnames(senior_df))) {
    stop("Adult and Senior featureImportance column names do not match.")
  }
  
  species_cols <- colnames(adult_df)
  
  # Merge studies: take unique union, then intersect with cohort metadata
  all_studies <- union(rownames(adult_df), rownames(senior_df))
  valid_studies <- intersect(all_studies, rownames(cohort_metadata))
  
  # Initialize empty matrix with valid study names
  merged_df <- matrix(NA, nrow = length(valid_studies), ncol = length(species_cols))
  rownames(merged_df) <- valid_studies
  colnames(merged_df) <- species_cols
  
  # Fill merged_df with values from adult/senior or row mean of both
  for (study in valid_studies) {
    in_adult <- study %in% rownames(adult_df)
    in_senior <- study %in% rownames(senior_df)
    
    if (in_adult && !in_senior) {
      merged_df[study, ] <- as.numeric(adult_df[study, ])
    } else if (!in_adult && in_senior) {
      merged_df[study, ] <- as.numeric(senior_df[study, ])
    } else if (in_adult && in_senior) {
      merged_df[study, ] <- as.numeric((adult_df[study, ] + senior_df[study, ]) / 2)
    }
  }
  
  # Add lifestyle prefix from cohort metadata
  lifestyle_prefix <- c(
    "IndustrializedUrban" = "IU__",
    "RuralTribal" = "RT__",
    "UrbanRuralMixed" = "URM__"
  )
  study_types <- cohort_metadata[rownames(merged_df), , drop = FALSE]
  prefixed_rownames <- paste0(lifestyle_prefix[study_types$`Cohort-Type`], rownames(merged_df))
  rownames(merged_df) <- prefixed_rownames
  
  # Rank scaling across rows
  scaled <- t(apply(merged_df, 1, function(x) {
    if (max(x) != min(x)) {
      (x - min(x)) / (max(x) - min(x))
    } else {
      rep(0, length(x))
    }
  }))
  rownames(scaled) <- rownames(merged_df)
  scaled_df <- as.data.frame(scaled)
  scaled_df[is.na(scaled_df)] <- 0
  
  # Compute Euclidean distance and PCoA
  dist_mat <- dist(scaled_df)
  pcoa <- cmdscale(dist_mat, k = 2)
  pcoa_df <- as.data.frame(pcoa)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  
  # Extract group from prefixed rowname
  pcoa_df$Group <- substr(rownames(pcoa_df), 1, regexpr("__", rownames(pcoa_df)) - 1)
  group_labels <- c("IU" = "IndustrializedUrban", "RT" = "RuralTribal", "URM" = "UrbanRuralMixed")
  pcoa_df$Group <- factor(group_labels[pcoa_df$Group], levels = c("IndustrializedUrban", "RuralTribal", "UrbanRuralMixed"))
  
  # PERMANOVA
  permanova <- adonis2(dist_mat ~ Group, data = pcoa_df, permutations = 999, method = "euclidean")
  r2 <- formatC(permanova$R2[1], format = "f", digits = 3)
  pval <- formatC(permanova$`Pr(>F)`[1], format = "e", digits = 2)
  permanova_text <- paste0("PERMANOVA\nR² = ", r2, "\np = ", pval)
  
  # Compute centroids
  centroids <- pcoa_df %>%
    group_by(Group) %>%
    summarise(PCoA1_centroid = mean(PCoA1), PCoA2_centroid = mean(PCoA2), .groups = "drop")
  pcoa_df <- left_join(pcoa_df, centroids, by = "Group")
  
  # Plot
  plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 4) +
    geom_segment(aes(xend = PCoA1_centroid, yend = PCoA2_centroid), linetype = "solid") +
    scale_color_manual(values = c("IndustrializedUrban" = "mediumpurple", "RuralTribal" = "seagreen", "UrbanRuralMixed" = "darkorange")) +
    theme_minimal() +
    labs(title = paste("PCoA Plot for", species, method),
         x = "Principal Coordinate 1",
         y = "Principal Coordinate 2") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 12),
          legend.position = "none")
  
  pdf_name <- paste0("PCoA_", species, "_lifestyle_", method, ".pdf")
  ggsave(pdf_name, plot = plot, width = 8, height = 5)
  cat("Saved plot:", pdf_name, "\n")
  cat("PERMANOVA result for", species, method, ":\n")
  print(permanova_text)
  
  return(scaled_df)
}

# List of your 8 Bifidobacterium species (use exact column names as in your featureImportance matrices)
bif_species_list <- c("detection", "longum", "adolescentis", "pseudocatenulatum","bifidum", "breve", "dentium", "catenulatum")

# Run the function for each species
for (species in bif_species_list) {
  cat("Running PCoA for:", species, "\n")
  make_pcoa_plot(species, "wgs", df_list, cohort_metadata)
}

scaled_dfs <- list()
for (species in bif_species_list) {
  cat("Running PCoA for:", species, "\n")
  scaled_dfs[[species]] <- make_pcoa_plot(species, "wgs", df_list, cohort_metadata)
}

#_________________________________________________________________________________________________________________
#Figure S14
#PCoA for IndustrializedUrban vs Others

# Load required libraries
library(tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)

# Load the RData file
load("OriginalAssociationScores_Modified.RData")
load("cohort_metadata.RData")

# Filter all objects in the environment
df_list <- mget(ls(pattern = "^df_association_"), envir = .GlobalEnv)


make_pcoa_plot <- function(species, method = "wgs", df_list = mget(ls(pattern = "df_association_")), cohort_metadata, output_dir = "PCoA_lifestyle_2") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  # Extract relevant adult and senior dataframes
  adult_df <- df_list[[paste0("df_association_adult_", species, "_", method)]]$featureImportance
  senior_df <- df_list[[paste0("df_association_senior_", species, "_", method)]]$featureImportance
  
  # Check column name consistency
  if (!identical(colnames(adult_df), colnames(senior_df))) {
    stop("Adult and Senior featureImportance column names do not match.")
  }
  
  species_cols <- colnames(adult_df)
  
  # Merge studies: take unique union, then intersect with cohort metadata
  all_studies <- union(rownames(adult_df), rownames(senior_df))
  valid_studies <- intersect(all_studies, rownames(cohort_metadata))
  
  # Initialize empty matrix with valid study names
  merged_df <- matrix(NA, nrow = length(valid_studies), ncol = length(species_cols))
  rownames(merged_df) <- valid_studies
  colnames(merged_df) <- species_cols
  
  # Fill merged_df with values from adult/senior or row mean of both
  for (study in valid_studies) {
    in_adult <- study %in% rownames(adult_df)
    in_senior <- study %in% rownames(senior_df)
    
    if (in_adult && !in_senior) {
      merged_df[study, ] <- as.numeric(adult_df[study, ])
    } else if (!in_adult && in_senior) {
      merged_df[study, ] <- as.numeric(senior_df[study, ])
    } else if (in_adult && in_senior) {
      merged_df[study, ] <- as.numeric((adult_df[study, ] + senior_df[study, ]) / 2)
    }
  }
  
  # Add lifestyle prefix from cohort metadata
  lifestyle_prefix <- c(
    "IndustrializedUrban" = "IU__",
    "RuralTribal" = "RT__",
    "UrbanRuralMixed" = "URM__"
  )
  study_types <- cohort_metadata[rownames(merged_df), , drop = FALSE]
  prefixed_rownames <- paste0(lifestyle_prefix[study_types$`Cohort-Type`], rownames(merged_df))
  rownames(merged_df) <- prefixed_rownames
  
  # Rank scaling across rows
  scaled <- t(apply(merged_df, 1, function(x) {
    if (max(x) != min(x)) {
      (x - min(x)) / (max(x) - min(x))
    } else {
      rep(0, length(x))
    }
  }))
  rownames(scaled) <- rownames(merged_df)
  scaled_df <- as.data.frame(scaled)
  scaled_df[is.na(scaled_df)] <- 0
  
  # Compute Euclidean distance and PCoA
  dist_mat <- dist(scaled_df)
  pcoa <- cmdscale(dist_mat, k = 2)
  pcoa_df <- as.data.frame(pcoa)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  
  # Extract group from prefixed rowname
  #pcoa_df$Group <- substr(rownames(pcoa_df), 1, regexpr("__", rownames(pcoa_df)) - 1)
  #group_labels <- c("IU" = "IndustrializedUrban", "RT" = "RuralTribal", "URM" = "UrbanRuralMixed")
  #pcoa_df$Group <- factor(group_labels[pcoa_df$Group], levels = c("IndustrializedUrban", "RuralTribal", "UrbanRuralMixed"))
  group_raw <- substr(rownames(pcoa_df), 1, regexpr("__", rownames(pcoa_df)) - 1)
  pcoa_df$Group <- ifelse(group_raw == "IU", "IndustrializedUrban", "Others")
  pcoa_df$Group <- factor(pcoa_df$Group, levels = c("IndustrializedUrban", "Others"))
  
  # PERMANOVA
  permanova <- adonis2(dist_mat ~ Group, data = pcoa_df, permutations = 999, method = "euclidean")
  r2 <- formatC(permanova$R2[1], format = "f", digits = 3)
  pval <- formatC(permanova$`Pr(>F)`[1], format = "e", digits = 2)
  permanova_text <- paste0("PERMANOVA\nR² = ", r2, "\np = ", pval)
  
  # Compute centroids
  centroids <- pcoa_df %>%
    group_by(Group) %>%
    summarise(PCoA1_centroid = mean(PCoA1), PCoA2_centroid = mean(PCoA2), .groups = "drop")
  pcoa_df <- left_join(pcoa_df, centroids, by = "Group")
  
  # Plot
  pdf(file.path(output_dir, paste0("PCoA_", species, "_IU_vs_Others_", method, ".pdf")),
      width = 8, height = 5)
  
  plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 4) +
    geom_segment(aes(xend = PCoA1_centroid, yend = PCoA2_centroid), linetype = "solid") +
    scale_color_manual(values = c("IndustrializedUrban" = "mediumpurple", "Others" = "hotpink")) +
    theme_minimal() +
    labs(title = paste("PCoA:IU vs Others-", species, method),
         x = "Principal Coordinate 1",
         y = "Principal Coordinate 2") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 12),
          legend.position = "none")
  print(plot)
  dev.off()
  
  #pdf_name <- paste0("PCoA_", species, "_IU_vs_Others_", method, ".pdf")
  #ggsave(pdf_name, plot = plot, width = 8, height = 5)
  #cat("Saved plot:", pdf_name, "\n")
  #print(plot)
  cat("PERMANOVA result for", species, method, ":\n")
  print(permanova_text)
  
  return(scaled_df)
}

#make_pcoa_plot("longum", "wgs", df_list, cohort_metadata)

# List of your 8 Bifidobacterium species (use exact column names as in your featureImportance matrices)
bif_species_list <- c("detection", "longum", "adolescentis", "pseudocatenulatum","bifidum", "breve", "dentium", "catenulatum")

# Run the function for each species
scaled_dfs <- list()
for (species in bif_species_list) {
  cat("Running IU vs Others PCoA for:", species, "\n")
  scaled_dfs[[species]] <- make_pcoa_plot(species, "wgs", df_list, cohort_metadata)
}

#____________________________________________________________________________________________________
#____________________________________________________________________________________________________


