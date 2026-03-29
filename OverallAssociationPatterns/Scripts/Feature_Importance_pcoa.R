#Figure S12
#PCoA for seq type

# Load required libraries
library(tidyverse)
library(vegan)
library(ggplot2)

# Load the file
load("/data/all_16s_WGS_combined_NEW.RData")

association_dfs <- list()

obj_names <- ls()

for (obj in obj_names) {
  
  x <- get(obj)
  
  if (!is.list(x)) next
  if (!all(c("FI_scaled", "metadata") %in% names(x))) next
  
  parts <- strsplit(obj, "_")[[1]]
  bifname <- parts[1]
  age     <- parts[2]
  
  # initialize levels if needed
  if (!age %in% names(association_dfs)) {
    association_dfs[[age]] <- list()
  }
  if (!bifname %in% names(association_dfs[[age]])) {
    association_dfs[[age]][[bifname]] <- list()
  }
  
  FI_scaled <- x$FI_scaled
  metadata  <- x$metadata
  
  common_studies <- intersect(rownames(FI_scaled), rownames(metadata))
  FI_scaled <- FI_scaled[common_studies, , drop = FALSE]
  metadata  <- metadata[common_studies, , drop = FALSE]
  
  # split by sequencing type
  studies_16s <- rownames(metadata)[metadata$seq_type == "16s"]
  studies_wgs <- rownames(metadata)[metadata$seq_type == "WGS"]
  
  association_dfs[[age]][[bifname]][["16S"]] <-
    FI_scaled[studies_16s, , drop = FALSE]
  
  association_dfs[[age]][[bifname]][["WGS"]] <-
    FI_scaled[studies_wgs, , drop = FALSE]
}


names(association_dfs)
names(association_dfs$adult)
names(association_dfs$senior)
names(association_dfs$infant)

outdir <- "/results/"

#pcoa plots
make_seqtype_pcoa_plot <- function(species, association_dfs) {
  
  # -----------------------------
  # Extract adult/senior FI_scaled
  # -----------------------------
  adult_16s  <- association_dfs$adult[[species]][["16S"]]
  senior_16s <- association_dfs$senior[[species]][["16S"]]
  
  adult_wgs  <- association_dfs$adult[[species]][["WGS"]]
  senior_wgs <- association_dfs$senior[[species]][["WGS"]]
  
  # -----------------------------
  # Sanity checks
  # -----------------------------
  # -----------------------------
  # Harmonise columns (adult vs senior)
  # -----------------------------
  common_cols_16s <- intersect(colnames(adult_16s), colnames(senior_16s))
  common_cols_wgs <- intersect(colnames(adult_wgs), colnames(senior_wgs))
  
  adult_16s  <- adult_16s[, common_cols_16s, drop = FALSE]
  senior_16s <- senior_16s[, common_cols_16s, drop = FALSE]
  
  adult_wgs  <- adult_wgs[, common_cols_wgs, drop = FALSE]
  senior_wgs <- senior_wgs[, common_cols_wgs, drop = FALSE]
  
  if (length(common_cols_16s) == 0 || length(common_cols_wgs) == 0) {
    stop(paste("No common features for species:", species))
  }

  # -----------------------------
  # Helper: merge adult & senior
  # -----------------------------
  merge_studies <- function(adult_df, senior_df) {
    
    species_cols <- colnames(adult_df)
    all_studies <- union(rownames(adult_df), rownames(senior_df))
    
    merged_df <- matrix(
      NA,
      nrow = length(all_studies),
      ncol = length(species_cols),
      dimnames = list(all_studies, species_cols)
    )
    
    for (study in all_studies) {
      in_adult  <- study %in% rownames(adult_df)
      in_senior <- study %in% rownames(senior_df)
      
      if (in_adult && !in_senior) {
        merged_df[study, ] <- as.numeric(adult_df[study, ])
      } else if (!in_adult && in_senior) {
        merged_df[study, ] <- as.numeric(senior_df[study, ])
      } else if (in_adult && in_senior) {
        merged_df[study, ] <- as.numeric((adult_df[study, ] + senior_df[study, ]) / 2)
      }
    }
    
    as.data.frame(merged_df)
  }
  
  # -----------------------------
  # Merge 16S and WGS separately
  # -----------------------------
  merged_16s <- merge_studies(adult_16s, senior_16s)
  merged_wgs <- merge_studies(adult_wgs, senior_wgs)
  
  rownames(merged_16s) <- paste0(rownames(merged_16s), "_16S")
  rownames(merged_wgs) <- paste0(rownames(merged_wgs), "_WGS")
  
  combined_df <- dplyr::bind_rows(merged_16s, merged_wgs)
  combined_df[is.na(combined_df)] <- 0
  
  # -----------------------------
  # Rank scaling per row
  # -----------------------------
  scaled_df <- t(apply(combined_df, 1, function(x) {
    if (max(x) != min(x)) {
      (x - min(x)) / (max(x) - min(x))
    } else {
      rep(0, length(x))
    }
  }))
  scaled_df <- as.data.frame(scaled_df)
  
  # -----------------------------
  # PCoA
  # -----------------------------
  dist_mat <- dist(scaled_df)
  pcoa <- cmdscale(dist_mat, k = 2)
  
  pcoa_df <- as.data.frame(pcoa)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  
  pcoa_df$SeqType <- ifelse(grepl("_16S$", rownames(pcoa_df)), "16S", "WGS")
  pcoa_df$SeqType <- factor(pcoa_df$SeqType, levels = c("16S", "WGS"))
  
  # -----------------------------
  # PERMANOVA
  # -----------------------------
  permanova <- vegan::adonis2(
    dist_mat ~ SeqType,
    data = pcoa_df,
    permutations = 999,
    method = "euclidean"
  )
  
  r2   <- formatC(permanova$R2[1], format = "f", digits = 3)
  pval <- formatC(permanova$`Pr(>F)`[1], format = "e", digits = 2)
  permanova_text <- paste0("PERMANOVA\nR² = ", r2, "\np = ", pval)
  
  # -----------------------------
  # Centroids
  # -----------------------------
  centroids <- pcoa_df %>%
    dplyr::group_by(SeqType) %>%
    dplyr::summarise(
      PCoA1_centroid = mean(PCoA1),
      PCoA2_centroid = mean(PCoA2),
      .groups = "drop"
    )
  
  pcoa_df <- dplyr::left_join(pcoa_df, centroids, by = "SeqType")
  
  # -----------------------------
  # Plot
  # -----------------------------
  plot <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color = SeqType)) +
    geom_point(size = 3) +
    geom_segment(
      aes(xend = PCoA1_centroid, yend = PCoA2_centroid),
      linetype = "solid"
    ) +
    scale_color_manual(values = c("16S" = "darkorange", "WGS" = "steelblue")) +
    theme_minimal() +
    labs(
      title = paste("PCoA: 16S vs WGS for", species),
      x = "Principal Coordinate 1",
      y = "Principal Coordinate 2"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
  
  pdf_name <- file.path(outdir, paste0("PCoA_", species, "_seqtype.pdf"))
  ggsave(pdf_name, plot = plot, width = 8, height = 5)
  
  cat("Saved plot:", pdf_name, "\n")
  cat("PERMANOVA result for", species, ":\n")
  cat(permanova_text, "\n")
  
  return(scaled_df)
}

bif_species_list <- c("detection", "longum", "adolescentis", "pseudocatenulatum", "bifidum", "dentium", "breve", "catenulatum", "animalis")

scaled_dfs <- list()

for (sp in bif_species_list) {
  cat("Running PCoA for:", sp, "\n")
  scaled_dfs[[sp]] <- make_seqtype_pcoa_plot(
    species = sp,
    association_dfs = association_dfs
  )
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
load("/data/all_16s_WGS_combined_NEW.RData")

species_list <- c(
  "detection", "longum", "adolescentis", "pseudocatenulatum",
  "bifidum", "dentium", "breve", "catenulatum", "animalis"
)

age_groups <- c("infant", "adult", "senior")
prefixes <- c("infant" = "I__", "adult" = "A__", "senior" = "S__")
group_labels <- c("I" = "Infant", "A" = "Adult", "S" = "Senior")

methods <- c("16s", "wgs")

dir.create("PCoA_AgeWise", showWarnings = FALSE)

# --------------------------------------------------
# Function: Age-wise PCoA (per sequencing type)
# --------------------------------------------------
make_agewise_pcoa_plot <- function(species, method) {
  
  feature_matrices <- list()
  
  for (age in age_groups) {
    
    obj_name <- paste0(species, "_", age, "_combined_new")
    if (!exists(obj_name)) next
    
    obj <- get(obj_name)
    if (!all(c("FI_scaled", "metadata") %in% names(obj))) next
    
    fi <- obj$FI_scaled
    meta <- obj$metadata
    
    # ensure matching order
    common <- intersect(rownames(fi), rownames(meta))
    fi   <- fi[common, , drop = FALSE]
    meta <- meta[common, , drop = FALSE]
    
    # --------------------------------------------------
    # Subset by sequencing type
    # --------------------------------------------------
    if (method == "16s") {
      keep <- meta$seq_type == "16s"
    } else {
      keep <- meta$seq_type == "WGS"
    }
    
    fi <- fi[keep, , drop = FALSE]
    if (nrow(fi) == 0) next
    
    rownames(fi) <- paste0(prefixes[age], rownames(fi))
    feature_matrices[[age]] <- fi
  }
  
  if (length(feature_matrices) < 2) {
    message("Skipping ", species, " (", method, "): insufficient age groups")
    return(NULL)
  }
  
  # --------------------------------------------------
  # Harmonize features
  # --------------------------------------------------
  common_cols <- Reduce(intersect, lapply(feature_matrices, colnames))
  feature_matrices <- lapply(feature_matrices, function(x) x[, common_cols, drop = FALSE])
  
  combined <- bind_rows(feature_matrices)
  combined[is.na(combined)] <- 0
  
  # --------------------------------------------------
  # Row-wise min–max scaling
  # --------------------------------------------------
  scaled <- t(apply(combined, 1, function(x) {
    if (max(x) != min(x)) {
      (x - min(x)) / (max(x) - min(x))
    } else {
      rep(0, length(x))
    }
  }))
  
  rownames(scaled) <- rownames(combined)
  scaled_df <- as.data.frame(scaled)
  
  # --------------------------------------------------
  # PCoA
  # --------------------------------------------------
  dist_mat <- dist(scaled_df)
  pcoa <- cmdscale(dist_mat, k = 2)
  
  pcoa_df <- as.data.frame(pcoa)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  pcoa_df$Group <- substr(rownames(pcoa_df), 1, 1)
  pcoa_df$Group <- factor(group_labels[pcoa_df$Group],
                          levels = c("Infant", "Adult", "Senior"))
  
  # --------------------------------------------------
  # PERMANOVA
  # --------------------------------------------------
  permanova <- adonis2(dist_mat ~ Group, data = pcoa_df, permutations = 999)
  r2 <- formatC(permanova$R2[1], digits = 3, format = "f")
  pval <- formatC(permanova$`Pr(>F)`[1], digits = 2, format = "e")
  
  # --------------------------------------------------
  # Centroids
  # --------------------------------------------------
  centroids <- pcoa_df %>%
    group_by(Group) %>%
    summarise(
      PCoA1_centroid = mean(PCoA1),
      PCoA2_centroid = mean(PCoA2),
      .groups = "drop"
    )
  
  pcoa_df <- left_join(pcoa_df, centroids, by = "Group")
  
  # --------------------------------------------------
  # Plot
  # --------------------------------------------------
  plot <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color = Group)) +
    geom_point(size = 4) +
    geom_segment(
      aes(xend = PCoA1_centroid, yend = PCoA2_centroid),
      linewidth = 0.6
    ) +
    scale_color_manual(values = c(
      "Infant" = "forestgreen",
      "Adult" = "orangered2",
      "Senior" = "royalblue3"
    )) +
    theme_minimal() +
    labs(
      title = paste("Age-wise PCoA –", species, "(", toupper(method), ")"),
      x = "PCoA1",
      y = "PCoA2"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
  
  # --------------------------------------------------
  # Save
  # --------------------------------------------------
  pdf_name <- file.path(
    "/results/",
    paste0("PCoA_AgeWise_", species, "_", method, ".pdf")
  )
  
  ggsave(pdf_name, plot, width = 5, height = 5)
  
  cat("Saved:", pdf_name,
      "| R² =", r2,
      "| p =", pval, "\n")
}

# --------------------------------------------------
# Run for all species × seq-type
# --------------------------------------------------
for (sp in species_list) {
  if (sp == "animalis") next  # skip animalis as infant group is not present
  for (m in methods) {
    make_agewise_pcoa_plot(sp, m)
  }
}

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
load("/data/all_16s_WGS_combined_NEW.RData")

# List of species for PCoA
bif_species_list <- c(
  "detection", "longum", "adolescentis", "pseudocatenulatum",
  "bifidum", "dentium", "breve", "catenulatum"
)

# Prefixes for lifestyle groups
lifestyle_prefix <- c(
  "IndustrializedUrban" = "IU__",
  "RuralTribal" = "RT__",
  "UrbanRuralMixed" = "URM__"
)

# Age groups to use
age_groups <- c("adult", "senior")

# Create output directory
dir.create("PCoA_Lifestyle", showWarnings = FALSE)

# --------------------------------------------------
# Function: PCoA by lifestyle (adult + senior)
# --------------------------------------------------
make_lifestyle_pcoa <- function(species) {
  
  feature_matrices <- list()
  metadata_list <- list()
  
  # Loop over adult + senior
  for (age in age_groups) {
    obj_name <- paste0(species, "_", age, "_combined_new")
    if (!exists(obj_name)) next
    
    obj <- get(obj_name)
    if (!all(c("FI_scaled", "metadata") %in% names(obj))) next
    
    fi <- obj$FI_scaled
    meta <- obj$metadata
    
    # Ensure matching rownames
    common_rows <- intersect(rownames(fi), rownames(meta))
    fi <- fi[common_rows, , drop = FALSE]
    meta <- meta[common_rows, , drop = FALSE]
    
    feature_matrices[[age]] <- fi
    metadata_list[[age]] <- meta
  }
  
  if (length(feature_matrices) < 1) {
    message("Skipping ", species, ": no adult/senior data")
    return(NULL)
  }
  
  # Combine matrices
  combined <- bind_rows(feature_matrices)
  combined[is.na(combined)] <- 0
  
  # Combine metadata
  combined_meta <- bind_rows(metadata_list)
  
  # Add lifestyle prefix to rownames
  prefixed_rownames <- paste0(lifestyle_prefix[combined_meta$lifestyle], rownames(combined))
  rownames(combined) <- prefixed_rownames
  
  # Row-wise min-max scaling
  scaled <- t(apply(combined, 1, function(x) {
    if (max(x) != min(x)) (x - min(x)) / (max(x) - min(x)) else rep(0, length(x))
  }))
  rownames(scaled) <- rownames(combined)
  scaled_df <- as.data.frame(scaled)
  scaled_df[is.na(scaled_df)] <- 0
  
  # Compute distance & PCoA
  dist_mat <- dist(scaled_df)
  pcoa <- cmdscale(dist_mat, k = 2)
  pcoa_df <- as.data.frame(pcoa)
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  
  # Extract lifestyle group
  pcoa_df$Group <- substr(rownames(pcoa_df), 1, regexpr("__", rownames(pcoa_df)) - 1)
  group_labels <- c("IU" = "IndustrializedUrban",
                    "RT" = "RuralTribal",
                    "URM" = "UrbanRuralMixed")
  pcoa_df$Group <- factor(group_labels[pcoa_df$Group],
                          levels = c("IndustrializedUrban", "RuralTribal", "UrbanRuralMixed"))
  
  # PERMANOVA
  permanova <- adonis2(dist_mat ~ Group, data = pcoa_df, permutations = 999)
  r2 <- formatC(permanova$R2[1], digits = 3, format = "f")
  pval <- formatC(permanova$`Pr(>F)`[1], digits = 2, format = "e")
  cat("Saved:", species, "| R² =", r2, "| p =", pval, "\n")
  
  # Centroids
  centroids <- pcoa_df %>%
    group_by(Group) %>%
    summarise(PCoA1_centroid = mean(PCoA1),
              PCoA2_centroid = mean(PCoA2),
              .groups = "drop")
  pcoa_df <- left_join(pcoa_df, centroids, by = "Group")
  
  # Plot
  plot <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color = Group)) +
    geom_point(size = 4) +
    geom_segment(aes(xend = PCoA1_centroid, yend = PCoA2_centroid), linetype = "solid") +
    scale_color_manual(values = c(
      "IndustrializedUrban" = "mediumpurple",
      "RuralTribal" = "seagreen",
      "UrbanRuralMixed" = "darkorange"
    )) +
    theme_minimal() +
    labs(title = paste("PCoA Plot –", species),
         x = "PCoA1", y = "PCoA2") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 12),
          legend.position = "none")
  
  # Save PDF
  pdf_name <- file.path("/results/", paste0("PCoA_", species, "_lifestyle.pdf"))
  ggsave(pdf_name, plot = plot, width = 8, height = 5)
  
  return(scaled_df)
}

# --------------------------------------------------
# Run for all 9 species
# --------------------------------------------------
scaled_dfs <- list()
for (species in bif_species_list) {
  cat("Running PCoA for:", species, "\n")
  scaled_dfs[[species]] <- make_lifestyle_pcoa(species)
}

#____________________________________________________________________________________________________
