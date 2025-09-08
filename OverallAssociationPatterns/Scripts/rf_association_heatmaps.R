
library(pheatmap)
library(dplyr)
library(xlsx)
library(ggplot2)
library(ggrepel)

load("G:/My Drive/Bifido_Project/AssociationHeatmaps/OriginalAssociationScores_Modified.RData")
load("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\cohort_metadata.RData")

adolescentis_infant_16s_FI <- df_association_infant_adolescentis_16s$featureImportance[,df_association_infant_adolescentis_16s$topFeatures]
adolescentis_infant_16s_pval <- df_association_infant_adolescentis_16s$pvalue[,df_association_infant_adolescentis_16s$topFeatures]
adolescentis_infant_16s_est <- df_association_infant_adolescentis_16s$est[,df_association_infant_adolescentis_16s$topFeatures]
adolescentis_infant_16s_association <- df_association_infant_adolescentis_16s$association[df_association_infant_adolescentis_16s$topFeatures,]
all(names(adolescentis_infant_16s_est) == rownames(adolescentis_infant_16s_association ))


adolescentis_infant_wgs_FI <- df_association_infant_adolescentis_wgs$featureImportance[,df_association_infant_adolescentis_wgs$topFeatures]
adolescentis_infant_wgs_pval <- df_association_infant_adolescentis_wgs$pvalue[,df_association_infant_adolescentis_wgs$topFeatures]
adolescentis_infant_wgs_est <- df_association_infant_adolescentis_wgs$est[,df_association_infant_adolescentis_wgs$topFeatures]
adolescentis_infant_wgs_association <- df_association_infant_adolescentis_wgs$association[df_association_infant_adolescentis_wgs$topFeatures,]
all(names(adolescentis_infant_wgs_est) == rownames(adolescentis_infant_wgs_association ))

df1 <- adolescentis_infant_16s_association[, 1:3]
df2 <- adolescentis_infant_wgs_association[, 1:3]

df1$Species <- rownames(df1)
df2$Species <- rownames(df2)
merged_df <- merge(df1, df2, by = "Species", all = TRUE)

merged_df[is.na(merged_df)] <- 0

merged_df$positive <- merged_df$positive.x + merged_df$positive.y
merged_df$negative <- merged_df$negative.x + merged_df$negative.y
merged_df$total    <- merged_df$total.x + merged_df$total.y

adolescentis_infant_association_combined <- merged_df[, c("Species", "positive", "negative", "total")]
rownames(adolescentis_infant_association_combined) <- adolescentis_infant_association_combined$Species
adolescentis_infant_association_combined$Species <- NULL
rm(merged_df)

adolescentis_infant_est_combined <- bind_rows(df_association_infant_adolescentis_16s$est,df_association_infant_adolescentis_wgs$est)
adolescentis_infant_est_combined <- adolescentis_infant_est_combined[,rownames(adolescentis_infant_association_combined)]
adolescentis_infant_est_combined[is.na(adolescentis_infant_est_combined)] <- 0

adolescentis_infant_FI_combined <- bind_rows(df_association_infant_adolescentis_16s$featureImportance,df_association_infant_adolescentis_wgs$featureImportance)
adolescentis_infant_FI_combined <- adolescentis_infant_FI_combined[,rownames(adolescentis_infant_association_combined)]
adolescentis_infant_FI_combined[is.na(adolescentis_infant_FI_combined)] <- 0

adolescentis_infant_pval_combined <- bind_rows(df_association_infant_adolescentis_16s$pvalue,df_association_infant_adolescentis_wgs$pvalue)
adolescentis_infant_pval_combined <- adolescentis_infant_pval_combined[,rownames(adolescentis_infant_association_combined)]
adolescentis_infant_pval_combined[is.na(adolescentis_infant_pval_combined)] <- 1

adolescentis_infant_oob_combined <- as.data.frame(rbind(df_association_infant_adolescentis_16s$oob,df_association_infant_adolescentis_wgs$oob))

est_mat <- as.matrix(adolescentis_infant_est_combined)
pval_mat <- as.matrix(adolescentis_infant_pval_combined)
dir_mat <- matrix(0, nrow = nrow(est_mat), ncol = ncol(est_mat))
rownames(dir_mat) <- rownames(est_mat)
colnames(dir_mat) <- colnames(est_mat)
dir_mat[pval_mat <= 0.05] <- 2 * sign(est_mat[pval_mat <= 0.05])
dir_mat[pval_mat > 0.05 & pval_mat <= 0.1] <- 1 * sign(est_mat[pval_mat > 0.05 & pval_mat <= 0.1])
adolescentis_infant_dir_combined <- as.data.frame(dir_mat)

adolescentis_infant_metadata_combined <- as.data.frame(matrix(NA,nrow = nrow(adolescentis_infant_oob_combined),ncol = 2))
rownames(adolescentis_infant_metadata_combined) <- rownames(adolescentis_infant_oob_combined)
colnames(adolescentis_infant_metadata_combined) <- c("seq_type","lifestyle")

adolescentis_infant_metadata_combined$seq_type <- ifelse(1:nrow(adolescentis_infant_metadata_combined) <= nrow(df_association_infant_adolescentis_16s$oob),"16s","wgs")

adolescentis_infant_metadata_combined$lifestyle <- cohort_metadata$`Cohort-Type`[match(rownames(adolescentis_infant_metadata_combined), rownames(cohort_metadata))]

adolescentis_infant_combined <- list(
  metadata = adolescentis_infant_metadata_combined,
  dir = adolescentis_infant_dir_combined,
  oob = adolescentis_infant_oob_combined,
  pval = adolescentis_infant_pval_combined,
  FI = adolescentis_infant_FI_combined,
  estimate = adolescentis_infant_est_combined,
  association = adolescentis_infant_association_combined)

###################################################################
## Actual Data generation function for heat maps
##################################################################

generate_combined_rf_association_list <- function(species_name, group_name, cohort_metadata) {
  # Construct object names dynamically
  obj_prefix_16s <- paste0("df_association_", group_name, "_", species_name, "_16s")
  obj_prefix_wgs <- paste0("df_association_", group_name, "_", species_name, "_wgs")
  
  # Get the actual data objects
  df_16s <- get(obj_prefix_16s)
  df_wgs <- get(obj_prefix_wgs)
  
  # Extract components from 16S
  est_16s <- df_16s$est[, df_16s$topFeatures]
  pval_16s <- df_16s$pvalue[, df_16s$topFeatures]
  FI_16s <- df_16s$featureImportance[, df_16s$topFeatures]
  assoc_16s <- df_16s$association[df_16s$topFeatures, ]
  
  # Extract components from WGS
  est_wgs <- df_wgs$est[, df_wgs$topFeatures]
  pval_wgs <- df_wgs$pvalue[, df_wgs$topFeatures]
  FI_wgs <- df_wgs$featureImportance[, df_wgs$topFeatures]
  assoc_wgs <- df_wgs$association[df_wgs$topFeatures, ]
  
  # Combine association data
  df1 <- assoc_16s[, 1:3]; df1$Species <- rownames(df1)
  df2 <- assoc_wgs[, 1:3]; df2$Species <- rownames(df2)
  
  merged_df <- merge(df1, df2, by = "Species", all = TRUE)
  merged_df[is.na(merged_df)] <- 0
  merged_df$positive <- merged_df$positive.x + merged_df$positive.y
  merged_df$negative <- merged_df$negative.x + merged_df$negative.y
  merged_df$total    <- merged_df$total.x + merged_df$total.y
  
  association_combined <- merged_df[, c("Species", "positive", "negative", "total")]
  rownames(association_combined) <- association_combined$Species
  association_combined$Species <- NULL
  
  # Combine est, FI, pval
  est_combined <- bind_rows(df_16s$est, df_wgs$est)
  est_combined <- est_combined[, rownames(association_combined)]
  est_combined[is.na(est_combined)] <- 0
  
  FI_combined <- bind_rows(df_16s$featureImportance, df_wgs$featureImportance)
  FI_combined <- FI_combined[, rownames(association_combined)]
  FI_combined[is.na(FI_combined)] <- 0
  
  rank_scale=function(x)
  {
    x <- rank(x);
    y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
    y <- ifelse(is.nan(y),0,y)
    return(y);
  }
  
  FI_combined_scaled <- as.data.frame(t(apply(FI_combined,1,rank_scale)))
  
  pval_combined <- bind_rows(df_16s$pvalue, df_wgs$pvalue)
  pval_combined <- pval_combined[, rownames(association_combined)]
  pval_combined[is.na(pval_combined)] <- 1
  
  oob_combined <- as.data.frame(rbind(df_16s$oob, df_wgs$oob))
  
  # Direction matrix
  est_mat <- as.matrix(est_combined)
  pval_mat <- as.matrix(pval_combined)
  dir_mat <- matrix(0, nrow = nrow(est_mat), ncol = ncol(est_mat),
                    dimnames = list(rownames(est_mat), colnames(est_mat)))
  dir_mat[pval_mat <= 0.05] <- 2 * sign(est_mat[pval_mat <= 0.05])
  dir_mat[pval_mat > 0.05 & pval_mat <= 0.1] <- 1 * sign(est_mat[pval_mat > 0.05 & pval_mat <= 0.1])
  dir_combined <- as.data.frame(dir_mat)
  
  # Metadata
  metadata_combined <- as.data.frame(matrix(NA, nrow = nrow(oob_combined), ncol = 2))
  rownames(metadata_combined) <- rownames(oob_combined)
  colnames(metadata_combined) <- c("seq_type", "lifestyle")
  
  metadata_combined$seq_type <- ifelse(seq_len(nrow(metadata_combined)) <= nrow(df_16s$oob), "16s", "wgs")
  metadata_combined$lifestyle <- cohort_metadata$`Cohort-Type`[match(rownames(metadata_combined), rownames(cohort_metadata))]
  
  # Return as a named list
  combined_list <- list(
    metadata = metadata_combined,
    dir = dir_combined,
    oob = oob_combined,
    pval = pval_combined,
    FI = FI_combined,
    FI_scaled = FI_combined_scaled,
    estimate = est_combined,
    association = association_combined
  )
  
  return(combined_list)
}

detection_infant_combined <- generate_combined_rf_association_list("detection", "infant", cohort_metadata)
detection_adult_combined <- generate_combined_rf_association_list("detection", "adult", cohort_metadata)
detection_senior_combined <- generate_combined_rf_association_list("detection", "senior", cohort_metadata)

adolescentis_infant_combined <- generate_combined_rf_association_list("adolescentis", "infant", cohort_metadata)
adolescentis_adult_combined <- generate_combined_rf_association_list("adolescentis", "adult", cohort_metadata)
adolescentis_senior_combined <- generate_combined_rf_association_list("adolescentis", "senior", cohort_metadata)

animalis_adult_combined <- generate_combined_rf_association_list("animalis", "adult", cohort_metadata)
animalis_senior_combined <- generate_combined_rf_association_list("animalis", "senior", cohort_metadata)

bifidum_infant_combined <- generate_combined_rf_association_list("bifidum", "infant", cohort_metadata)
bifidum_adult_combined <- generate_combined_rf_association_list("bifidum", "adult", cohort_metadata)
bifidum_senior_combined <- generate_combined_rf_association_list("bifidum", "senior", cohort_metadata)

catenulatum_infant_combined <- generate_combined_rf_association_list("catenulatum", "infant", cohort_metadata)
catenulatum_adult_combined <- generate_combined_rf_association_list("catenulatum", "adult", cohort_metadata)
catenulatum_senior_combined <- generate_combined_rf_association_list("catenulatum", "senior", cohort_metadata)

dentium_infant_combined <- generate_combined_rf_association_list("dentium", "infant", cohort_metadata)
dentium_adult_combined <- generate_combined_rf_association_list("dentium", "adult", cohort_metadata)
dentium_senior_combined <- generate_combined_rf_association_list("dentium", "senior", cohort_metadata)

breve_infant_combined <- generate_combined_rf_association_list("breve", "infant", cohort_metadata)
breve_adult_combined <- generate_combined_rf_association_list("breve", "adult", cohort_metadata)
breve_senior_combined <- generate_combined_rf_association_list("breve", "senior", cohort_metadata)

longum_infant_combined <- generate_combined_rf_association_list("longum", "infant", cohort_metadata)
longum_adult_combined <- generate_combined_rf_association_list("longum", "adult", cohort_metadata)
longum_senior_combined <- generate_combined_rf_association_list("longum", "senior", cohort_metadata)

pseudocatenulatum_infant_combined <- generate_combined_rf_association_list("pseudocatenulatum", "infant", cohort_metadata)
pseudocatenulatum_adult_combined <- generate_combined_rf_association_list("pseudocatenulatum", "adult", cohort_metadata)
pseudocatenulatum_senior_combined <- generate_combined_rf_association_list("pseudocatenulatum", "senior", cohort_metadata)

save(detection_infant_combined,detection_adult_combined,detection_senior_combined,adolescentis_infant_combined, adolescentis_adult_combined,adolescentis_senior_combined,animalis_adult_combined,animalis_senior_combined,bifidum_infant_combined,bifidum_adult_combined,bifidum_senior_combined,catenulatum_infant_combined,catenulatum_adult_combined,catenulatum_senior_combined,dentium_infant_combined,dentium_adult_combined,dentium_senior_combined,breve_infant_combined,breve_adult_combined,breve_senior_combined, longum_infant_combined,longum_adult_combined,longum_senior_combined,pseudocatenulatum_infant_combined,pseudocatenulatum_adult_combined,pseudocatenulatum_senior_combined, file = "all_16s_WGS_combined.RData")


############## Heat map function ################

library(ComplexHeatmap)
library(circlize)
library(grid)
library(circlize)

FI_scaled <- as.matrix(detection_infant_combined[["FI_scaled"]])
dir_matrix <- as.matrix(detection_infant_combined[["dir"]])

col_fun <- colorRamp2(
  c(min(FI_scaled, na.rm = TRUE),
    quantile(FI_scaled, 0.25, na.rm = TRUE),
    0,
    quantile(FI_scaled, 0.75, na.rm = TRUE),
    max(FI_scaled, na.rm = TRUE)),
  
  c("#ffffcc",  # soft yellow
    "#ffe0b2",  # light peach
    "white",    # center
    "#f4a582",  # light coral
    "#d6604d"   # coral
  )
)

ht <- Heatmap(
  FI_scaled,
  name = "FI_scaled",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",  
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(title = "Feature Importance"),
  width = unit(ncol(FI_scaled) * 1.7, "mm"),
  height = unit(nrow(FI_scaled) * 1.7, "mm"),
  border = TRUE,
  layer_fun = function(j_index, i_index, x, y, width, height, fill) {
    grid.rect(x, y, width = width, height = height,
              gp = gpar(col = "grey60", fill = NA, lwd = 0.4))
    
    val <- dir_matrix[cbind(i_index, j_index)]
    shape_size <- unit(0.4, "mm")
    
    blue_circles <- which(val == 2)
    red_circles <- which(val == -2)
    blue_triangles <- which(val == 1)
    red_triangles <- which(val == -1)
    
    if (length(blue_circles) > 0)
      grid.circle(x[blue_circles], y[blue_circles], r = shape_size, gp = gpar(fill = "blue", col = NA))
    if (length(red_circles) > 0)
      grid.circle(x[red_circles], y[red_circles], r = shape_size, gp = gpar(fill = "darkred", col = NA))
    
    draw_triangles <- function(indices, color) {
      for (idx in indices) {
        xi <- x[idx]
        yi <- y[idx]
        triangle_x <- unit.c(xi, xi - shape_size, xi + shape_size)
        triangle_y <- unit.c(yi + shape_size, yi - shape_size, yi - shape_size)
        grid.polygon(triangle_x, triangle_y, gp = gpar(fill = color, col = NA))
      }
    }
    
    draw_triangles(blue_triangles, "blue")
    draw_triangles(red_triangles, "darkred")
  }
)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\detection_infant_FI_heatmap.pdf", height = 4, width = 13)
draw(ht)
dev.off()

### association line plots ###

clustered_studies <- rownames(FI_scaled)[row_order(draw(ht))]
clustered_species <- colnames(FI_scaled)[column_order(draw(ht))]

association_ordered <- detection_infant_combined$association[clustered_species,]

association_long <- association_ordered %>%
  tibble::rownames_to_column("Species") %>%
  select(Species, positive, negative) %>%
  pivot_longer(cols = c("positive", "negative"), names_to = "Association", values_to = "Count")

association_long$Species <- factor(association_long$Species, levels = rownames(association_ordered))

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\detection_infant_lineplot.pdf", height = 1, width = 10)
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

metadata_ordered <- detection_infant_combined$metadata[clustered_studies,]

### Seq-type heatmap strip ### 

seq_type_vec <- metadata_ordered$seq_type

seq_type_colors <- structure(
  c("#B2DF8A", "#CAB2D6"),
  names = unique(seq_type_vec)
)

seq_type_mat <- matrix(seq_type_vec, ncol = 1)
rownames(seq_type_mat) <- rownames(metadata_ordered)

ht_seq_type <- Heatmap(
  seq_type_mat,
  name = "Seq Type",
  col = seq_type_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(2, "mm"),
  height = unit(nrow(seq_type_mat) * 1.7, "mm"),
  border = FALSE,  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = "Seq Type"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(
      x = x, y = y, width = width, height = height,
      gp = gpar(col = "grey30", fill = NA, lwd = 0.5)
    )
  }
)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\detection_infant_SeqType_strip.pdf", height = 4, width = 3)
draw(ht_seq_type)
dev.off()

### lifestyle heat map strip ###

lifestyle_vec <- metadata_ordered$lifestyle

lifestyle_color_map <- c(
  "IndustrializedUrban" = "hotpink1",
  "UrbanRuralMixed" = "lightpink",
  "RuralTribal" = "violetred4"
)

lifestyle_colors <- lifestyle_color_map[unique(lifestyle_vec)]

lifestyle_mat <- matrix(lifestyle_vec, ncol = 1)
rownames(lifestyle_mat) <- rownames(metadata_ordered)

ht_lifestyle <- Heatmap(
  lifestyle_mat,
  name = "Lifestyle",
  col = lifestyle_colors,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(2, "mm"),
  height = unit(nrow(lifestyle_mat) * 1.7, "mm"),
  border = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = "Lifestyle"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(
      x = x, y = y, width = width, height = height,
      gp = gpar(col = "grey30", fill = NA, lwd = 0.5)
    )
  }
)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\detection_infant_Lifestyle_strip.pdf", height = 4, width = 3)
draw(ht_lifestyle)
dev.off()

######## OOB Correlation heat map strip #####

oob_ordered <- detection_infant_combined$oob[clustered_studies,]

oob_corr_vec <- oob_ordered$oob_corr
names(oob_corr_vec) <- rownames(oob_ordered)

oob_corr_mat <- matrix(oob_corr_vec, ncol = 1)
rownames(oob_corr_mat) <- rownames(oob_ordered)

col_fun_corr <- colorRamp2(c(-1, 0, 1), c("mediumpurple4", "white", "orange"))

oob_p_vec <- oob_ordered$oob_p

ht_oob_corr <- Heatmap(
  oob_corr_mat,
  name = "OOB Corr",
  col = col_fun_corr,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(3.5, "mm"),
  height = unit(nrow(oob_corr_mat) * 1.7, "mm"),
  border = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = "OOB Corr"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    grid.rect(
      x = x, y = y, width = width, height = height,
      gp = gpar(col = "grey30", fill = NA, lwd = 0.5)
    )
    
    
    if (!is.na(oob_p_vec[i]) && oob_p_vec[i] <= 0.05) {
      grid.text("*", x = x, y = y, gp = gpar(col = "black", fontsize = 8, fontface = "bold"))
    }
  }
)

pdf("G:/My Drive/Bifido_Project/AssociationHeatmaps/detection_infant_oob_strip.pdf", height = 4, width = 3)
draw(ht_oob_corr)
dev.off()

##################################################################
### Function to generate all plots together
##################################################################

library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)
library(ggrepel)
library(tidyverse)

load("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\all_16s_WGS_combined.RData")

generate_association_heatmaps <- function(data_list, prefix, age_group, output_dir = "G:/My Drive/Bifido_Project/AssociationHeatmaps") {
  
  ## === 1. Feature Importance Heatmap ===
  FI_scaled <- as.matrix(data_list[["FI_scaled"]])
  dir_matrix <- as.matrix(data_list[["dir"]])
  
  col_fun <- colorRamp2(
    c(min(FI_scaled, na.rm = TRUE),
      quantile(FI_scaled, 0.25, na.rm = TRUE),
      0,
      quantile(FI_scaled, 0.75, na.rm = TRUE),
      max(FI_scaled, na.rm = TRUE)),
    c("#ffffcc", "#ffe0b2", "white", "#f4a582", "#d6604d")
  )
  
  ht <- Heatmap(
    FI_scaled,
    name = "FI_scaled",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 4.6),
    column_names_gp = gpar(fontsize = 4.6),
    heatmap_legend_param = list(title = "Feature Importance"),
    width = unit(ncol(FI_scaled) * 1.5, "mm"),
    height = unit(nrow(FI_scaled) * 1.5, "mm"),
    border = TRUE,
    layer_fun = function(j_index, i_index, x, y, width, height, fill) {
      grid.rect(x, y, width = width, height = height,
                gp = gpar(col = "grey60", fill = NA, lwd = 0.4))
      val <- dir_matrix[cbind(i_index, j_index)]
      shape_size <- unit(0.4, "mm")
      
      draw_triangles <- function(indices, color) {
        for (idx in indices) {
          triangle_x <- unit.c(x[idx], x[idx] - shape_size, x[idx] + shape_size)
          triangle_y <- unit.c(y[idx] + shape_size, y[idx] - shape_size, y[idx] - shape_size)
          grid.polygon(triangle_x, triangle_y, gp = gpar(fill = color, col = NA))
        }
      }
      
      grid.circle(x[val == 2], y[val == 2], r = shape_size, gp = gpar(fill = "blue", col = NA))
      grid.circle(x[val == -2], y[val == -2], r = shape_size, gp = gpar(fill = "darkred", col = NA))
      draw_triangles(which(val == 1), "blue")
      draw_triangles(which(val == -1), "darkred")
    }
  )
  
  pdf(file.path(output_dir, paste0(prefix, "_", age_group, "_FI_heatmap.pdf")), height = 8, width = 16)
  draw(ht)
  dev.off()
  
  clustered_studies <- rownames(FI_scaled)[row_order(draw(ht))]
  clustered_species <- colnames(FI_scaled)[column_order(draw(ht))]
  
  ## === 2. Association Line Plot ===
  association_ordered <- data_list$association[clustered_species, ]
  association_long <- association_ordered %>%
    rownames_to_column("Species") %>%
    select(Species, positive, negative) %>%
    pivot_longer(cols = c("positive", "negative"), names_to = "Association", values_to = "Count")
  
  association_long$Species <- factor(association_long$Species, levels = rownames(association_ordered))
  
  p <- ggplot(association_long, aes(x = Species, y = Count, group = Association, color = Association)) +
    geom_line(linewidth = 0.6) +
    geom_point(shape = 16, size = 1.8) +
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
  ggsave(filename = file.path(output_dir, paste0(prefix, "_", age_group, "_lineplot.pdf")),plot = p, width = 16, height = 1.5, units = "in")
  
  metadata_ordered <- data_list$metadata[clustered_studies, ]
  
  ## === 3. Seq Type Strip ===
  seq_type_vec <- metadata_ordered$seq_type
  seq_type_colors <- structure(
    c("#B2DF8A", "#CAB2D6"),
    names = unique(seq_type_vec)
  )
  seq_type_mat <- matrix(seq_type_vec, ncol = 1)
  rownames(seq_type_mat) <- rownames(metadata_ordered)
  
  ht_seq_type <- Heatmap(
    seq_type_mat,
    name = "Seq Type",
    col = seq_type_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    width = unit(2, "mm"),
    height = unit(nrow(seq_type_mat) * 1.7, "mm"),
    border = FALSE,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = "Seq Type"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width = width, height = height, gp = gpar(col = "grey30", fill = NA, lwd = 0.5))
    }
  )
  
  pdf(file.path(output_dir, paste0(prefix, "_", age_group, "_SeqType_strip.pdf")), height = 10, width = 3)
  draw(ht_seq_type)
  dev.off()
  
  ## === 4. Lifestyle Strip ===
  lifestyle_vec <- metadata_ordered$lifestyle
  lifestyle_color_map <- c(
    "IndustrializedUrban" = "hotpink1",
    "UrbanRuralMixed" = "lightpink",
    "RuralTribal" = "violetred4"
  )
  lifestyle_colors <- lifestyle_color_map[names(lifestyle_color_map) %in% unique(lifestyle_vec)]
  
  lifestyle_mat <- matrix(lifestyle_vec, ncol = 1)
  rownames(lifestyle_mat) <- rownames(metadata_ordered)
  
  ht_lifestyle <- Heatmap(
    lifestyle_mat,
    name = "Lifestyle",
    col = lifestyle_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    width = unit(2, "mm"),
    height = unit(nrow(lifestyle_mat) * 1.7, "mm"),
    border = FALSE,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = "Lifestyle"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width = width, height = height, gp = gpar(col = "grey30", fill = NA, lwd = 0.5))
    }
  )
  
  pdf(file.path(output_dir, paste0(prefix, "_", age_group, "_Lifestyle_strip.pdf")), height = 10, width = 3)
  draw(ht_lifestyle)
  dev.off()
  
  ## === 5. OOB Corr Strip ===
  oob_ordered <- data_list$oob[clustered_studies, ]
  
  oob_corr_vec <- oob_ordered$oob_corr
  names(oob_corr_vec) <- rownames(oob_ordered)
  
  oob_corr_mat <- matrix(oob_corr_vec, ncol = 1)
  rownames(oob_corr_mat) <- rownames(oob_ordered)
  
  oob_p_vec <- oob_ordered$oob_p
  
  col_fun_corr <- colorRamp2(c(-1, 0, 1), c("mediumpurple4", "white", "orange"))
  
  ht_oob_corr <- Heatmap(
    oob_corr_mat,
    name = "OOB Corr",
    col = col_fun_corr,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    width = unit(3.5, "mm"),
    height = unit(nrow(oob_corr_mat) * 1.7, "mm"),
    border = FALSE,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = "OOB Corr"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width = width, height = height, gp = gpar(col = "grey30", fill = NA, lwd = 0.5))
      if (!is.na(oob_p_vec[i]) && oob_p_vec[i] <= 0.05) {
        grid.text("*", x = x, y = y, gp = gpar(col = "black", fontsize = 8, fontface = "bold"))
      }
    }
  )
  
  pdf(file.path(output_dir, paste0(prefix, "_", age_group, "_oob_strip.pdf")), height = 10, width = 3)
  draw(ht_oob_corr)
  dev.off()
}

generate_association_heatmaps(detection_infant_combined, "detection", "infant")
generate_association_heatmaps(detection_adult_combined, "detection", "adult")
generate_association_heatmaps(detection_senior_combined, "detection", "senior")

generate_association_heatmaps(longum_infant_combined, "longum", "infant")
generate_association_heatmaps(longum_adult_combined, "longum", "adult")
generate_association_heatmaps(longum_senior_combined, "longum", "senior")

generate_association_heatmaps(adolescentis_infant_combined, "adolescentis", "infant")
generate_association_heatmaps(adolescentis_adult_combined, "adolescentis", "adult")
generate_association_heatmaps(adolescentis_senior_combined, "adolescentis", "senior")

generate_association_heatmaps(pseudocatenulatum_infant_combined, "pseudocatenulatum", "infant")
generate_association_heatmaps(pseudocatenulatum_adult_combined, "pseudocatenulatum", "adult")
generate_association_heatmaps(pseudocatenulatum_senior_combined, "pseudocatenulatum", "senior")

generate_association_heatmaps(bifidum_infant_combined, "bifidum", "infant")
generate_association_heatmaps(bifidum_adult_combined, "bifidum", "adult")
generate_association_heatmaps(bifidum_senior_combined, "bifidum", "senior")

generate_association_heatmaps(breve_infant_combined, "breve", "infant")
generate_association_heatmaps(breve_adult_combined, "breve", "adult")
generate_association_heatmaps(breve_senior_combined, "breve", "senior")

generate_association_heatmaps(dentium_infant_combined, "dentium", "infant")
generate_association_heatmaps(dentium_adult_combined, "dentium", "adult")
generate_association_heatmaps(dentium_senior_combined, "dentium", "senior")

generate_association_heatmaps(catenulatum_infant_combined, "catenulatum", "infant")
generate_association_heatmaps(catenulatum_adult_combined, "catenulatum", "adult")
generate_association_heatmaps(catenulatum_senior_combined, "catenulatum", "senior")

generate_association_heatmaps(animalis_adult_combined, "animalis", "adult")
generate_association_heatmaps(animalis_senior_combined, "animalis", "senior")


############ Supplementary Table Feature Importance Senior (TableS4)  #######

senior_detection <- as.data.frame(cbind(detection_senior_combined$metadata,detection_senior_combined$FI_scaled,detection_senior_combined$oob))
write.xlsx(senior_detection,"G:\\My Drive\\Bifido_Project\\senior_detection.xlsx")

senior_longum <- as.data.frame(cbind(longum_senior_combined$metadata,longum_senior_combined$FI_scaled,longum_senior_combined$oob))
write.xlsx(senior_longum,"G:\\My Drive\\Bifido_Project\\senior_longum.xlsx")

senior_adolescentis <- as.data.frame(cbind(adolescentis_senior_combined$metadata,adolescentis_senior_combined$FI_scaled,adolescentis_senior_combined$oob))
write.xlsx(senior_adolescentis,"G:\\My Drive\\Bifido_Project\\senior_adolescentis.xlsx")

senior_pseudocatenulatum <- as.data.frame(cbind(pseudocatenulatum_senior_combined$metadata,pseudocatenulatum_senior_combined$FI_scaled,pseudocatenulatum_senior_combined$oob))
write.xlsx(senior_pseudocatenulatum,"G:\\My Drive\\Bifido_Project\\senior_pseudocatenulatum.xlsx")

senior_bifidum <- as.data.frame(cbind(bifidum_senior_combined$metadata,bifidum_senior_combined$FI_scaled,bifidum_senior_combined$oob))
write.xlsx(senior_bifidum,"G:\\My Drive\\Bifido_Project\\senior_bifidum.xlsx")

senior_breve <- as.data.frame(cbind(breve_senior_combined$metadata,breve_senior_combined$FI_scaled,breve_senior_combined$oob))
write.xlsx(senior_breve,"G:\\My Drive\\Bifido_Project\\senior_breve.xlsx")

senior_dentium <- as.data.frame(cbind(dentium_senior_combined$metadata,dentium_senior_combined$FI_scaled,dentium_senior_combined$oob))
write.xlsx(senior_dentium,"G:\\My Drive\\Bifido_Project\\senior_dentium.xlsx")

senior_catenulatum <- as.data.frame(cbind(catenulatum_senior_combined$metadata,catenulatum_senior_combined$FI_scaled,catenulatum_senior_combined$oob))
write.xlsx(senior_catenulatum,"G:\\My Drive\\Bifido_Project\\senior_catenulatum.xlsx")

senior_animalis <- as.data.frame(cbind(animalis_senior_combined$metadata,animalis_senior_combined$FI_scaled,animalis_senior_combined$oob))
write.xlsx(senior_animalis,"G:\\My Drive\\Bifido_Project\\senior_animalis.xlsx")

####################################################################
## ==== Supplementary tables Dir TableS5 to TableS10 =====

adolescentis_infant_16s_association <- df_association_infant_adolescentis_16s$association[df_association_infant_adolescentis_16s$topFeatures,]
est_mat <- as.matrix(df_association_infant_adolescentis_16s$est[,df_association_infant_adolescentis_16s$topFeatures])
pval_mat <- as.matrix(df_association_infant_adolescentis_16s$pvalue[,df_association_infant_adolescentis_16s$topFeatures])
dir_mat <- matrix(0, nrow = nrow(est_mat), ncol = ncol(est_mat))
rownames(dir_mat) <- rownames(est_mat)
colnames(dir_mat) <- colnames(est_mat)
dir_mat[pval_mat <= 0.05] <- 2 * sign(est_mat[pval_mat <= 0.05])
dir_mat[pval_mat > 0.05 & pval_mat <= 0.1] <- 1 * sign(est_mat[pval_mat > 0.05 & pval_mat <= 0.1])
adolescentis_infant_16s_dir <- as.data.frame(dir_mat)
adolescentis_infant_16s_metadata <- cohort_metadata[rownames(adolescentis_infant_16s_dir),,drop = F]
adolescentis_infant_16s_suppl <- as.data.frame(rbind(t(adolescentis_infant_16s_metadata),t(adolescentis_infant_16s_dir)))
adolescentis_infant_16s_suppl$total_pos[2:nrow(adolescentis_infant_16s_suppl)] <- adolescentis_infant_16s_association$positive
adolescentis_infant_16s_suppl$total_neg[2:nrow(adolescentis_infant_16s_suppl)] <- adolescentis_infant_16s_association$negative

#### Function to run ###

generate_suppl_df <- function(association_list, cohort_metadata) {
  
  # Extract association, estimate, p-value, and top features
  association <- association_list$association[association_list$topFeatures, ]
  est_mat <- as.matrix(association_list$est[, association_list$topFeatures])
  pval_mat <- as.matrix(association_list$pvalue[, association_list$topFeatures])
  
  # Initialize direction matrix
  dir_mat <- matrix(0, nrow = nrow(est_mat), ncol = ncol(est_mat))
  rownames(dir_mat) <- rownames(est_mat)
  colnames(dir_mat) <- colnames(est_mat)
  
  # Fill direction matrix
  dir_mat[pval_mat <= 0.05] <- 2 * sign(est_mat[pval_mat <= 0.05])
  dir_mat[pval_mat > 0.05 & pval_mat <= 0.1] <- 1 * sign(est_mat[pval_mat > 0.05 & pval_mat <= 0.1])
  
  # Convert to data frame
  dir_df <- as.data.frame(dir_mat)
  
  # Match metadata for the rows
  metadata_matched <- cohort_metadata[rownames(dir_df), , drop = FALSE]
  
  # Combine metadata and directional matrix
  combined_suppl <- as.data.frame(rbind(t(metadata_matched), t(dir_df)))
  
  # Add total positive and negative associations
  combined_suppl$total_pos[2:nrow(combined_suppl)] <- association$positive
  combined_suppl$total_neg[2:nrow(combined_suppl)] <- association$negative
  
  return(combined_suppl)
}

##=====  INFANT 16S ======##

detection_infant_16s_suppl <- generate_suppl_df(df_association_infant_detection_16s, cohort_metadata)
write.xlsx(detection_infant_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\detection_infant_16s_suppl.xlsx")

longum_infant_16s_suppl <- generate_suppl_df(df_association_infant_longum_16s, cohort_metadata)
write.xlsx(longum_infant_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\longum_infant_16s_suppl.xlsx")

adolescentis_infant_16s_suppl <- generate_suppl_df(df_association_infant_adolescentis_16s, cohort_metadata)
write.xlsx(adolescentis_infant_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\adolescentis_infant_16s_suppl.xlsx")

pseudocatenulatum_infant_16s_suppl <- generate_suppl_df(df_association_infant_pseudocatenulatum_16s, cohort_metadata)
write.xlsx(pseudocatenulatum_infant_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\pseudocatenulatum_infant_16s_suppl.xlsx")

bifidum_infant_16s_suppl <- generate_suppl_df(df_association_infant_bifidum_16s, cohort_metadata)
write.xlsx(bifidum_infant_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\bifidum_infant_16s_suppl.xlsx")

breve_infant_16s_suppl <- generate_suppl_df(df_association_infant_breve_16s, cohort_metadata)
write.xlsx(breve_infant_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\breve_infant_16s_suppl.xlsx")

dentium_infant_16s_suppl <- generate_suppl_df(df_association_infant_dentium_16s, cohort_metadata)
write.xlsx(dentium_infant_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\dentium_infant_16s_suppl.xlsx")

catenulatum_infant_16s_suppl <- generate_suppl_df(df_association_infant_catenulatum_16s, cohort_metadata)
write.xlsx(catenulatum_infant_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\catenulatum_infant_16s_suppl.xlsx")



##=====  Adult 16S ======##

detection_adult_16s_suppl <- generate_suppl_df(df_association_adult_detection_16s, cohort_metadata)
write.xlsx(detection_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\detection_adult_16s_suppl.xlsx")

longum_adult_16s_suppl <- generate_suppl_df(df_association_adult_longum_16s, cohort_metadata)
write.xlsx(longum_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\longum_adult_16s_suppl.xlsx")

adolescentis_adult_16s_suppl <- generate_suppl_df(df_association_adult_adolescentis_16s, cohort_metadata)
write.xlsx(adolescentis_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\adolescentis_adult_16s_suppl.xlsx")

pseudocatenulatum_adult_16s_suppl <- generate_suppl_df(df_association_adult_pseudocatenulatum_16s, cohort_metadata)
write.xlsx(pseudocatenulatum_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\pseudocatenulatum_adult_16s_suppl.xlsx")

bifidum_adult_16s_suppl <- generate_suppl_df(df_association_adult_bifidum_16s, cohort_metadata)
write.xlsx(bifidum_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\bifidum_adult_16s_suppl.xlsx")

breve_adult_16s_suppl <- generate_suppl_df(df_association_adult_breve_16s, cohort_metadata)
write.xlsx(breve_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\breve_adult_16s_suppl.xlsx")

dentium_adult_16s_suppl <- generate_suppl_df(df_association_adult_dentium_16s, cohort_metadata)
write.xlsx(dentium_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\dentium_adult_16s_suppl.xlsx")

catenulatum_adult_16s_suppl <- generate_suppl_df(df_association_adult_catenulatum_16s, cohort_metadata)
write.xlsx(catenulatum_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\catenulatum_adult_16s_suppl.xlsx")

animalis_adult_16s_suppl <- generate_suppl_df(df_association_adult_animalis_16s, cohort_metadata)
write.xlsx(animalis_adult_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\animalis_adult_16s_suppl.xlsx")


##=====  Senior 16S ======##

detection_senior_16s_suppl <- generate_suppl_df(df_association_senior_detection_16s, cohort_metadata)
write.xlsx(detection_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\detection_senior_16s_suppl.xlsx")

longum_senior_16s_suppl <- generate_suppl_df(df_association_senior_longum_16s, cohort_metadata)
write.xlsx(longum_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\longum_senior_16s_suppl.xlsx")

adolescentis_senior_16s_suppl <- generate_suppl_df(df_association_senior_adolescentis_16s, cohort_metadata)
write.xlsx(adolescentis_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\adolescentis_senior_16s_suppl.xlsx")

pseudocatenulatum_senior_16s_suppl <- generate_suppl_df(df_association_senior_pseudocatenulatum_16s, cohort_metadata)
write.xlsx(pseudocatenulatum_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\pseudocatenulatum_senior_16s_suppl.xlsx")

bifidum_senior_16s_suppl <- generate_suppl_df(df_association_senior_bifidum_16s, cohort_metadata)
write.xlsx(bifidum_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\bifidum_senior_16s_suppl.xlsx")

breve_senior_16s_suppl <- generate_suppl_df(df_association_senior_breve_16s, cohort_metadata)
write.xlsx(breve_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\breve_senior_16s_suppl.xlsx")

dentium_senior_16s_suppl <- generate_suppl_df(df_association_senior_dentium_16s, cohort_metadata)
write.xlsx(dentium_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\dentium_senior_16s_suppl.xlsx")

catenulatum_senior_16s_suppl <- generate_suppl_df(df_association_senior_catenulatum_16s, cohort_metadata)
write.xlsx(catenulatum_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\catenulatum_senior_16s_suppl.xlsx")

animalis_senior_16s_suppl <- generate_suppl_df(df_association_senior_animalis_16s, cohort_metadata)
write.xlsx(animalis_senior_16s_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\animalis_senior_16s_suppl.xlsx")


##=====  INFANT WGS ======##

detection_infant_wgs_suppl <- generate_suppl_df(df_association_infant_detection_wgs, cohort_metadata)
write.xlsx(detection_infant_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\detection_infant_wgs_suppl.xlsx")

longum_infant_wgs_suppl <- generate_suppl_df(df_association_infant_longum_wgs, cohort_metadata)
write.xlsx(longum_infant_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\longum_infant_wgs_suppl.xlsx")

adolescentis_infant_wgs_suppl <- generate_suppl_df(df_association_infant_adolescentis_wgs, cohort_metadata)
write.xlsx(adolescentis_infant_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\adolescentis_infant_wgs_suppl.xlsx")

pseudocatenulatum_infant_wgs_suppl <- generate_suppl_df(df_association_infant_pseudocatenulatum_wgs, cohort_metadata)
write.xlsx(pseudocatenulatum_infant_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\pseudocatenulatum_infant_wgs_suppl.xlsx")

bifidum_infant_wgs_suppl <- generate_suppl_df(df_association_infant_bifidum_wgs, cohort_metadata)
write.xlsx(bifidum_infant_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\bifidum_infant_wgs_suppl.xlsx")

breve_infant_wgs_suppl <- generate_suppl_df(df_association_infant_breve_wgs, cohort_metadata)
write.xlsx(breve_infant_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\breve_infant_wgs_suppl.xlsx")

dentium_infant_wgs_suppl <- generate_suppl_df(df_association_infant_dentium_wgs, cohort_metadata)
write.xlsx(dentium_infant_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\dentium_infant_wgs_suppl.xlsx")

catenulatum_infant_wgs_suppl <- generate_suppl_df(df_association_infant_catenulatum_wgs, cohort_metadata)
write.xlsx(catenulatum_infant_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\catenulatum_infant_wgs_suppl.xlsx")


##=====  Adult WGS ======##

detection_adult_wgs_suppl <- generate_suppl_df(df_association_adult_detection_wgs, cohort_metadata)
write.xlsx(detection_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\detection_adult_wgs_suppl.xlsx")

longum_adult_wgs_suppl <- generate_suppl_df(df_association_adult_longum_wgs, cohort_metadata)
write.xlsx(longum_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\longum_adult_wgs_suppl.xlsx")

adolescentis_adult_wgs_suppl <- generate_suppl_df(df_association_adult_adolescentis_wgs, cohort_metadata)
write.xlsx(adolescentis_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\adolescentis_adult_wgs_suppl.xlsx")

pseudocatenulatum_adult_wgs_suppl <- generate_suppl_df(df_association_adult_pseudocatenulatum_wgs, cohort_metadata)
write.xlsx(pseudocatenulatum_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\pseudocatenulatum_adult_wgs_suppl.xlsx")

bifidum_adult_wgs_suppl <- generate_suppl_df(df_association_adult_bifidum_wgs, cohort_metadata)
write.xlsx(bifidum_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\bifidum_adult_wgs_suppl.xlsx")

breve_adult_wgs_suppl <- generate_suppl_df(df_association_adult_breve_wgs, cohort_metadata)
write.xlsx(breve_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\breve_adult_wgs_suppl.xlsx")

dentium_adult_wgs_suppl <- generate_suppl_df(df_association_adult_dentium_wgs, cohort_metadata)
write.xlsx(dentium_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\dentium_adult_wgs_suppl.xlsx")

catenulatum_adult_wgs_suppl <- generate_suppl_df(df_association_adult_catenulatum_wgs, cohort_metadata)
write.xlsx(catenulatum_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\catenulatum_adult_wgs_suppl.xlsx")

animalis_adult_wgs_suppl <- generate_suppl_df(df_association_adult_animalis_wgs, cohort_metadata)
write.xlsx(animalis_adult_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\animalis_adult_wgs_suppl.xlsx")


##=====  Senior WGS ======##

detection_senior_wgs_suppl <- generate_suppl_df(df_association_senior_detection_wgs, cohort_metadata)
write.xlsx(detection_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\detection_senior_wgs_suppl.xlsx")

longum_senior_wgs_suppl <- generate_suppl_df(df_association_senior_longum_wgs, cohort_metadata)
write.xlsx(longum_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\longum_senior_wgs_suppl.xlsx")

adolescentis_senior_wgs_suppl <- generate_suppl_df(df_association_senior_adolescentis_wgs, cohort_metadata)
write.xlsx(adolescentis_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\adolescentis_senior_wgs_suppl.xlsx")

pseudocatenulatum_senior_wgs_suppl <- generate_suppl_df(df_association_senior_pseudocatenulatum_wgs, cohort_metadata)
write.xlsx(pseudocatenulatum_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\pseudocatenulatum_senior_wgs_suppl.xlsx")

bifidum_senior_wgs_suppl <- generate_suppl_df(df_association_senior_bifidum_wgs, cohort_metadata)
write.xlsx(bifidum_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\bifidum_senior_wgs_suppl.xlsx")

breve_senior_wgs_suppl <- generate_suppl_df(df_association_senior_breve_wgs, cohort_metadata)
write.xlsx(breve_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\breve_senior_wgs_suppl.xlsx")

dentium_senior_wgs_suppl <- generate_suppl_df(df_association_senior_dentium_wgs, cohort_metadata)
write.xlsx(dentium_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\dentium_senior_wgs_suppl.xlsx")

catenulatum_senior_wgs_suppl <- generate_suppl_df(df_association_senior_catenulatum_wgs, cohort_metadata)
write.xlsx(catenulatum_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\catenulatum_senior_wgs_suppl.xlsx")

animalis_senior_wgs_suppl <- generate_suppl_df(df_association_senior_animalis_wgs, cohort_metadata)
write.xlsx(animalis_senior_wgs_suppl,"G:\\My Drive\\Bifido_Project\\Supplementary_tables_Dir\\animalis_senior_wgs_suppl.xlsx")

####################################################################
## WGS vs 16S association scores heatmap
####################################################################

load("G:\\.shortcut-targets-by-id\\1mJ9DmgZE4NfNNbH72zZNyUWw4tVWKsRf\\Bif_Manuscript\\df_association_all.RData")

common_species_16s_wgs <- intersect(rownames(df_association_all_16s),rownames(df_association_all_wgs))


####################### INFANT ###################
library(psych)

infant_16s_commonSP <- df_association_all_16s[common_species_16s_wgs,grep("infant",names(df_association_all_16s))]
infant_wgs_commonSP <- df_association_all_wgs[common_species_16s_wgs,grep("infant",names(df_association_all_wgs))]

names(infant_16s_commonSP) <- paste0(names(infant_16s_commonSP),"_16s")
names(infant_wgs_commonSP) <- paste0(names(infant_wgs_commonSP),"_wgs")

n <- ncol(infant_wgs_commonSP)
corr_matrix <- matrix(NA, nrow = n, ncol = n)
pval_matrix <- matrix(NA, nrow = n, ncol = n)

colnames(corr_matrix) <- colnames(infant_16s_commonSP)
rownames(corr_matrix) <- colnames(infant_wgs_commonSP)

colnames(pval_matrix) <- colnames(infant_16s_commonSP)
rownames(pval_matrix) <- colnames(infant_wgs_commonSP)

result <- corr.test(infant_16s_commonSP,infant_wgs_commonSP, method = 'spearman')

corr_df_infant <- as.data.frame(result$r)
pval_df_infant <- as.data.frame(result$p)

dir_infant <- matrix(0, nrow = nrow(corr_df_infant), ncol = ncol(corr_df_infant))
rownames(dir_infant) <- rownames(corr_df_infant)
colnames(dir_infant) <- colnames(corr_df_infant)

for (i in 1:nrow(corr_df_infant)) {
  for (j in 1:ncol(corr_df_infant)) {
    corr_val <- corr_df_infant[i, j]
    p_val <- pval_df_infant[i, j]
    
    if (!is.na(corr_val) && !is.na(p_val)) {
      if (corr_val > 0 && p_val <= 0.05) {
        dir_infant[i, j] <- 2
      } else if (corr_val > 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_infant[i, j] <- 1
      } else if (corr_val < 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_infant[i, j] <- -1
      } else if (corr_val < 0 && p_val <= 0.05) {
        dir_infant[i, j] <- -2
      }
    }
  }
}

dir_infant <- as.data.frame(dir_infant)

corr_diag <- diag(as.matrix(corr_df_infant))
pval_diag <- diag(as.matrix(pval_df_infant))

infant_corr_16s_wgs <- data.frame(
  corr = corr_diag,
  pval = pval_diag,
  row.names = rownames(corr_df_infant)
)


library(ggplot2)

infant_corr_16s_wgs$signif <- ifelse(infant_corr_16s_wgs$pval <= 0.05, "*", "")
infant_corr_16s_wgs <- infant_corr_16s_wgs[-8,]
infant_corr_16s_wgs$species <- gsub("infant_","",rownames(infant_corr_16s_wgs))
infant_corr_16s_wgs$species <- gsub("_16s","",infant_corr_16s_wgs$species)

infant_corr_16s_wgs$species <- factor(infant_corr_16s_wgs$species, levels = infant_corr_16s_wgs$species)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\infant_corr_16s_wgs.pdf", height = 2.5, width = 6)
ggplot(infant_corr_16s_wgs, aes(x = species, y = corr)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(aes(label = signif, y = corr + 0.02), size = 6) +
  coord_flip() +
  theme_minimal() +
  labs(x = "", y = "Correlation", title = "16S vs WGS Correlation in Infant") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),     
        panel.grid = element_blank(),              
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

##### Making the heatmap ####

colnames(dir_infant) <- gsub("infant_", "", colnames(dir_infant))
rownames(dir_infant) <- gsub("infant_", "", rownames(dir_infant))
colnames(dir_infant) <- gsub("_wgs","", names(dir_infant))
rownames(dir_infant) <- gsub("_16s","", rownames(dir_infant))

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  breaks = c(-2, -1, 0, 1, 2),
  colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B") 
)


ht <- Heatmap(
  as.matrix(dir_infant),
  name = "Association",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 90,
  column_names_side = "bottom",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "darkgray", fill = NA, lwd = 0.5))
    
    r_val <- round(corr_df_infant[i, j], 2)
    if (!is.na(r_val)) {
      grid.text(label = r_val, x = x, y = y, gp = gpar(fontsize = 7))
    }
  },
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title = "Association"
  ),
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  column_title = "Association scores between 16S and WGS in infant"
)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\corr_16S_WGS_infant.pdf", height = 6, width = 6)
draw(ht)
dev.off()



##################### ADULT #############################


adult_16s_commonSP <- df_association_all_16s[common_species_16s_wgs,grep("adult",names(df_association_all_16s))]
adult_wgs_commonSP <- df_association_all_wgs[common_species_16s_wgs,grep("adult",names(df_association_all_wgs))]

names(adult_16s_commonSP) <- paste0(names(adult_16s_commonSP),"_16s")
names(adult_wgs_commonSP) <- paste0(names(adult_wgs_commonSP),"_wgs")

n <- ncol(adult_wgs_commonSP)
corr_matrix <- matrix(NA, nrow = n, ncol = n)
pval_matrix <- matrix(NA, nrow = n, ncol = n)

colnames(corr_matrix) <- colnames(adult_16s_commonSP)
rownames(corr_matrix) <- colnames(adult_wgs_commonSP)

colnames(pval_matrix) <- colnames(adult_16s_commonSP)
rownames(pval_matrix) <- colnames(adult_wgs_commonSP)

result <- corr.test(adult_16s_commonSP,adult_wgs_commonSP, method = 'spearman')

corr_df_adult <- as.data.frame(result$r)
pval_df_adult <- as.data.frame(result$p)

dir_adult <- matrix(0, nrow = nrow(corr_df_adult), ncol = ncol(corr_df_adult))
rownames(dir_adult) <- rownames(corr_df_adult)
colnames(dir_adult) <- colnames(corr_df_adult)

for (i in 1:nrow(corr_df_adult)) {
  for (j in 1:ncol(corr_df_adult)) {
    corr_val <- corr_df_adult[i, j]
    p_val <- pval_df_adult[i, j]
    
    if (!is.na(corr_val) && !is.na(p_val)) {
      if (corr_val > 0 && p_val <= 0.05) {
        dir_adult[i, j] <- 2
      } else if (corr_val > 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_adult[i, j] <- 1
      } else if (corr_val < 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_adult[i, j] <- -1
      } else if (corr_val < 0 && p_val <= 0.05) {
        dir_adult[i, j] <- -2
      }
    }
  }
}

dir_adult <- as.data.frame(dir_adult)

##### Making the heatmap ####

colnames(dir_adult) <- gsub("adult_", "", colnames(dir_adult))
rownames(dir_adult) <- gsub("adult_", "", rownames(dir_adult))
colnames(dir_adult) <- gsub("_wgs","", names(dir_adult))
rownames(dir_adult) <- gsub("_16s","", rownames(dir_adult))

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  breaks = c(-2, -1, 0, 1, 2),
  colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B") 
)


ht <- Heatmap(
  as.matrix(dir_adult),
  name = "Association",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 90,
  column_names_side = "bottom",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "darkgray", fill = NA, lwd = 0.5))
    
    r_val <- round(corr_df_adult[i, j], 2)
    if (!is.na(r_val)) {
      grid.text(label = r_val, x = x, y = y, gp = gpar(fontsize = 7))
    }
  },
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title = "Association"
  ),
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  column_title = "Association scores between 16S and WGS in adult"
)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\corr_16S_WGS_adult.pdf", height = 6, width = 6)
draw(ht)
dev.off()


##################### SENIOR #############################


senior_16s_commonSP <- df_association_all_16s[common_species_16s_wgs,grep("senior",names(df_association_all_16s))]
senior_wgs_commonSP <- df_association_all_wgs[common_species_16s_wgs,grep("senior",names(df_association_all_wgs))]

names(senior_16s_commonSP) <- paste0(names(senior_16s_commonSP),"_16s")
names(senior_wgs_commonSP) <- paste0(names(senior_wgs_commonSP),"_wgs")

n <- ncol(senior_wgs_commonSP)
corr_matrix <- matrix(NA, nrow = n, ncol = n)
pval_matrix <- matrix(NA, nrow = n, ncol = n)

colnames(corr_matrix) <- colnames(senior_16s_commonSP)
rownames(corr_matrix) <- colnames(senior_wgs_commonSP)

colnames(pval_matrix) <- colnames(senior_16s_commonSP)
rownames(pval_matrix) <- colnames(senior_wgs_commonSP)

result <- corr.test(senior_16s_commonSP,senior_wgs_commonSP, method = 'spearman')

corr_df_senior <- as.data.frame(result$r)
pval_df_senior <- as.data.frame(result$p)

dir_senior <- matrix(0, nrow = nrow(corr_df_senior), ncol = ncol(corr_df_senior))
rownames(dir_senior) <- rownames(corr_df_senior)
colnames(dir_senior) <- colnames(corr_df_senior)

for (i in 1:nrow(corr_df_senior)) {
  for (j in 1:ncol(corr_df_senior)) {
    corr_val <- corr_df_senior[i, j]
    p_val <- pval_df_senior[i, j]
    
    if (!is.na(corr_val) && !is.na(p_val)) {
      if (corr_val > 0 && p_val <= 0.05) {
        dir_senior[i, j] <- 2
      } else if (corr_val > 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_senior[i, j] <- 1
      } else if (corr_val < 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_senior[i, j] <- -1
      } else if (corr_val < 0 && p_val <= 0.05) {
        dir_senior[i, j] <- -2
      }
    }
  }
}

dir_senior <- as.data.frame(dir_senior)

##### Making the heatmap ####

colnames(dir_senior) <- gsub("senior_", "", colnames(dir_senior))
rownames(dir_senior) <- gsub("senior_", "", rownames(dir_senior))
colnames(dir_senior) <- gsub("_wgs","", names(dir_senior))
rownames(dir_senior) <- gsub("_16s","", rownames(dir_senior))

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  breaks = c(-2, -1, 0, 1, 2),
  colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B") 
)


ht <- Heatmap(
  as.matrix(dir_senior),
  name = "Association",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 90,
  column_names_side = "bottom",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "darkgray", fill = NA, lwd = 0.5))
    
    r_val <- round(corr_df_senior[i, j], 2)
    if (!is.na(r_val)) {
      grid.text(label = r_val, x = x, y = y, gp = gpar(fontsize = 7))
    }
  },
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title = "Association"
  ),
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  column_title = "Association scores between 16S and WGS in senior"
)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\corr_16S_WGS_senior.pdf", height = 6, width = 6)
draw(ht)
dev.off()


##################### ADULT + SENIOR Merged  #############################

mean_adult_senior_16s <- (adult_16s_commonSP + senior_16s_commonSP) / 2
mean_adult_senior_wgs <- (adult_wgs_commonSP + senior_wgs_commonSP) / 2


names(mean_adult_senior_16s) <- gsub("adult","adult_senior",names(mean_adult_senior_16s))
names(mean_adult_senior_wgs) <- gsub("adult","adult_senior",names(mean_adult_senior_wgs))

n <- ncol(mean_adult_senior_16s)
corr_matrix <- matrix(NA, nrow = n, ncol = n)
pval_matrix <- matrix(NA, nrow = n, ncol = n)

colnames(corr_matrix) <- colnames(mean_adult_senior_16s)
rownames(corr_matrix) <- colnames(mean_adult_senior_16s)

colnames(pval_matrix) <- colnames(mean_adult_senior_16s)
rownames(pval_matrix) <- colnames(mean_adult_senior_16s)

result <- corr.test(mean_adult_senior_16s,mean_adult_senior_wgs, method = 'spearman')

corr_df_adult_senior <- as.data.frame(result$r)
pval_df_adult_senior <- as.data.frame(result$p)

dir_adult_senior <- matrix(0, nrow = nrow(corr_df_adult_senior), ncol = ncol(corr_df_adult_senior))
rownames(dir_adult_senior) <- rownames(corr_df_adult_senior)
colnames(dir_adult_senior) <- colnames(corr_df_adult_senior)

for (i in 1:nrow(corr_df_adult_senior)) {
  for (j in 1:ncol(corr_df_adult_senior)) {
    corr_val <- corr_df_adult_senior[i, j]
    p_val <- pval_df_adult_senior[i, j]
    
    if (!is.na(corr_val) && !is.na(p_val)) {
      if (corr_val > 0 && p_val <= 0.05) {
        dir_adult_senior[i, j] <- 2
      } else if (corr_val > 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_adult_senior[i, j] <- 1
      } else if (corr_val < 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_adult_senior[i, j] <- -1
      } else if (corr_val < 0 && p_val <= 0.05) {
        dir_adult_senior[i, j] <- -2
      }
    }
  }
}

dir_adult_senior <- as.data.frame(dir_adult_senior)

##### Making the heatmap ####

colnames(dir_adult_senior) <- gsub("adult_senior_", "", colnames(dir_adult_senior))
rownames(dir_adult_senior) <- gsub("adult_senior_", "", rownames(dir_adult_senior))
colnames(dir_adult_senior) <- gsub("_wgs","", names(dir_adult_senior))
rownames(dir_adult_senior) <- gsub("_16s","", rownames(dir_adult_senior))

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  breaks = c(-2, -1, 0, 1, 2),
  colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B") 
)


ht <- Heatmap(
  as.matrix(dir_adult_senior),
  name = "Association",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 90,
  column_names_side = "bottom",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height,
              gp = gpar(col = "darkgray", fill = NA, lwd = 0.5))
    
    r_val <- round(corr_df_adult_senior[i, j], 2)
    if (!is.na(r_val)) {
      grid.text(label = r_val, x = x, y = y, gp = gpar(fontsize = 7))
    }
  },
  heatmap_legend_param = list(
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2"),
    title = "Association"
  ),
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  column_title = "Association scores between 16S and WGS in adult+senior"
)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\corr_16S_WGS_adult_senior_merged.pdf", height = 6, width = 6)
draw(ht)
dev.off()

corr_diag <- diag(as.matrix(corr_df_adult_senior))
pval_diag <- diag(as.matrix(pval_df_adult_senior))

adult_senior_corr_16s_wgs <- data.frame(
  corr = corr_diag,
  pval = pval_diag,
  row.names = rownames(corr_df_adult_senior)
)


adult_senior_corr_16s_wgs$signif <- ifelse(adult_senior_corr_16s_wgs$pval <= 0.05, "*", "")
adult_senior_corr_16s_wgs$species <- gsub("adult_senior_","",rownames(adult_senior_corr_16s_wgs))
adult_senior_corr_16s_wgs$species <- gsub("_16s","",adult_senior_corr_16s_wgs$species)

adult_senior_corr_16s_wgs$species <- factor(adult_senior_corr_16s_wgs$species, levels = adult_senior_corr_16s_wgs$species)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\adult_senior_corr_16s_wgs.pdf", height = 2.5, width = 6)
ggplot(adult_senior_corr_16s_wgs, aes(x = species, y = corr)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(aes(label = signif, y = corr + 0.02), size = 6) +
  coord_flip() +
  theme_minimal() +
  labs(x = "", y = "Correlation", title = "16S vs WGS Correlation in adult+senior") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),     
        panel.grid = element_blank(),              
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

################### Industrialized Urban ##############

IndustrializedUrban_16s_commonSP <- df_association_all_16s[common_species_16s_wgs,grep("IndustrializedUrban",names(df_association_all_16s))]
IndustrializedUrban_wgs_commonSP <- df_association_all_wgs[common_species_16s_wgs,grep("IndustrializedUrban",names(df_association_all_wgs))]

names(IndustrializedUrban_16s_commonSP) <- paste0(names(IndustrializedUrban_16s_commonSP),"_16s")
names(IndustrializedUrban_wgs_commonSP) <- paste0(names(IndustrializedUrban_wgs_commonSP),"_wgs")

n <- ncol(IndustrializedUrban_wgs_commonSP)
corr_matrix <- matrix(NA, nrow = n, ncol = n)
pval_matrix <- matrix(NA, nrow = n, ncol = n)

colnames(corr_matrix) <- colnames(IndustrializedUrban_16s_commonSP)
rownames(corr_matrix) <- colnames(IndustrializedUrban_wgs_commonSP)

colnames(pval_matrix) <- colnames(IndustrializedUrban_16s_commonSP)
rownames(pval_matrix) <- colnames(IndustrializedUrban_wgs_commonSP)

result <- corr.test(IndustrializedUrban_16s_commonSP,IndustrializedUrban_wgs_commonSP, method = 'spearman')

corr_df_IndustrializedUrban <- as.data.frame(result$r)
pval_df_IndustrializedUrban <- as.data.frame(result$p)

dir_IndustrializedUrban <- matrix(0, nrow = nrow(corr_df_IndustrializedUrban), ncol = ncol(corr_df_IndustrializedUrban))
rownames(dir_IndustrializedUrban) <- rownames(corr_df_IndustrializedUrban)
colnames(dir_IndustrializedUrban) <- colnames(corr_df_IndustrializedUrban)

for (i in 1:nrow(corr_df_IndustrializedUrban)) {
  for (j in 1:ncol(corr_df_IndustrializedUrban)) {
    corr_val <- corr_df_IndustrializedUrban[i, j]
    p_val <- pval_df_IndustrializedUrban[i, j]
    
    if (!is.na(corr_val) && !is.na(p_val)) {
      if (corr_val > 0 && p_val <= 0.05) {
        dir_IndustrializedUrban[i, j] <- 2
      } else if (corr_val > 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_IndustrializedUrban[i, j] <- 1
      } else if (corr_val < 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_IndustrializedUrban[i, j] <- -1
      } else if (corr_val < 0 && p_val <= 0.05) {
        dir_IndustrializedUrban[i, j] <- -2
      }
    }
  }
}

dir_IndustrializedUrban <- as.data.frame(dir_IndustrializedUrban)

corr_diag <- diag(as.matrix(corr_df_IndustrializedUrban))
pval_diag <- diag(as.matrix(pval_df_IndustrializedUrban))

IndustrializedUrban_corr_16s_wgs <- data.frame(
  corr = corr_diag,
  pval = pval_diag,
  row.names = rownames(corr_df_IndustrializedUrban)
)

library(ggplot2)

IndustrializedUrban_corr_16s_wgs$signif <- ifelse(IndustrializedUrban_corr_16s_wgs$pval <= 0.05, "*", "")
IndustrializedUrban_corr_16s_wgs <- IndustrializedUrban_corr_16s_wgs[-8,]
IndustrializedUrban_corr_16s_wgs$species <- gsub("IndustrializedUrban_","",rownames(IndustrializedUrban_corr_16s_wgs))
IndustrializedUrban_corr_16s_wgs$species <- gsub("_16s","",IndustrializedUrban_corr_16s_wgs$species)

IndustrializedUrban_corr_16s_wgs$species <- factor(IndustrializedUrban_corr_16s_wgs$species, levels = IndustrializedUrban_corr_16s_wgs$species)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\IndustrializedUrban_corr_16s_wgs.pdf", height = 2.5, width = 6)
ggplot(IndustrializedUrban_corr_16s_wgs, aes(x = species, y = corr)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = signif, y = corr + 0.02), size = 6) +
  coord_flip() +
  theme_minimal() +
  labs(x = "", y = "Correlation", title = "16S vs WGS Correlation in IndustrializedUrban") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),     
        panel.grid = element_blank(),              
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()


################### Others ##############

Others_16s_commonSP <- df_association_all_16s[common_species_16s_wgs,grep("Others",names(df_association_all_16s))]

UrbanRuralMixed_wgs_commonSP <- df_association_all_wgs[common_species_16s_wgs,grep("UrbanRuralMixed",names(df_association_all_wgs))]
RuralTribal_wgs_commonSP <- df_association_all_wgs[common_species_16s_wgs,grep("RuralTribal",names(df_association_all_wgs))]

Others_wgs_commonSP <- (UrbanRuralMixed_wgs_commonSP + RuralTribal_wgs_commonSP)/2

names(Others_wgs_commonSP) <- gsub("UrbanRuralMixed","Others",names(Others_wgs_commonSP))
Others_wgs_commonSP <- Others_wgs_commonSP[,1:5]

names(Others_16s_commonSP) <- paste0(names(Others_16s_commonSP),"_16s")
names(Others_wgs_commonSP) <- paste0(names(Others_wgs_commonSP),"_wgs")

n <- ncol(Others_wgs_commonSP)
corr_matrix <- matrix(NA, nrow = n, ncol = n)
pval_matrix <- matrix(NA, nrow = n, ncol = n)

colnames(corr_matrix) <- colnames(Others_16s_commonSP)
rownames(corr_matrix) <- colnames(Others_wgs_commonSP)

colnames(pval_matrix) <- colnames(Others_16s_commonSP)
rownames(pval_matrix) <- colnames(Others_wgs_commonSP)

result <- corr.test(Others_16s_commonSP,Others_wgs_commonSP, method = 'spearman')

corr_df_Others <- as.data.frame(result$r)
pval_df_Others <- as.data.frame(result$p)

dir_Others <- matrix(0, nrow = nrow(corr_df_Others), ncol = ncol(corr_df_Others))
rownames(dir_Others) <- rownames(corr_df_Others)
colnames(dir_Others) <- colnames(corr_df_Others)

for (i in 1:nrow(corr_df_Others)) {
  for (j in 1:ncol(corr_df_Others)) {
    corr_val <- corr_df_Others[i, j]
    p_val <- pval_df_Others[i, j]
    
    if (!is.na(corr_val) && !is.na(p_val)) {
      if (corr_val > 0 && p_val <= 0.05) {
        dir_Others[i, j] <- 2
      } else if (corr_val > 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_Others[i, j] <- 1
      } else if (corr_val < 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_Others[i, j] <- -1
      } else if (corr_val < 0 && p_val <= 0.05) {
        dir_Others[i, j] <- -2
      }
    }
  }
}

dir_Others <- as.data.frame(dir_Others)

corr_diag <- diag(as.matrix(corr_df_Others))
pval_diag <- diag(as.matrix(pval_df_Others))

Others_corr_16s_wgs <- data.frame(
  corr = corr_diag,
  pval = pval_diag,
  row.names = rownames(corr_df_Others)
)

Others_corr_16s_wgs$signif <- ifelse(Others_corr_16s_wgs$pval <= 0.05, "*", "")
Others_corr_16s_wgs$species <- gsub("Others_","",rownames(Others_corr_16s_wgs))
Others_corr_16s_wgs$species <- gsub("_16s","",Others_corr_16s_wgs$species)

Others_corr_16s_wgs$species <- factor(Others_corr_16s_wgs$species, levels = Others_corr_16s_wgs$species)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\Others_corr_16s_wgs.pdf", height = 2.5, width = 6)
ggplot(Others_corr_16s_wgs, aes(x = species, y = corr)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(aes(label = signif, y = corr + 0.02), size = 6) +
  coord_flip() +
  theme_minimal() +
  labs(x = "", y = "Correlation", title = "16S vs WGS Correlation in Others") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),     
        panel.grid = element_blank(),              
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

####################################################################

load("G:\\.shortcut-targets-by-id\\1mJ9DmgZE4NfNNbH72zZNyUWw4tVWKsRf\\Bif_Manuscript\\AssociationScores_16s.RData")
load("G:\\.shortcut-targets-by-id\\1mJ9DmgZE4NfNNbH72zZNyUWw4tVWKsRf\\Bif_Manuscript\\AssociationScores.RData")

common_species_IU_16s_wgs <- intersect(rownames(df_association_adult_16s_IndustrializedUrban),rownames(df_association_adult_wgs_IndustrializedUrban))

IndustrializedUrban_16s_commonSP <- df_association_adult_16s_IndustrializedUrban[common_species_IU_16s_wgs,]
IndustrializedUrban_wgs_commonSP <- df_association_adult_wgs_IndustrializedUrban[common_species_IU_16s_wgs,]

names(IndustrializedUrban_16s_commonSP) <- paste0(names(IndustrializedUrban_16s_commonSP),"_16s")
names(IndustrializedUrban_wgs_commonSP) <- paste0(names(IndustrializedUrban_wgs_commonSP),"_wgs")

n <- ncol(IndustrializedUrban_wgs_commonSP)
corr_matrix <- matrix(NA, nrow = n, ncol = n)
pval_matrix <- matrix(NA, nrow = n, ncol = n)

colnames(corr_matrix) <- colnames(IndustrializedUrban_16s_commonSP)
rownames(corr_matrix) <- colnames(IndustrializedUrban_wgs_commonSP)

colnames(pval_matrix) <- colnames(IndustrializedUrban_16s_commonSP)
rownames(pval_matrix) <- colnames(IndustrializedUrban_wgs_commonSP)

result <- corr.test(IndustrializedUrban_16s_commonSP,IndustrializedUrban_wgs_commonSP, method = 'spearman')

corr_df_IndustrializedUrban <- as.data.frame(result$r)
pval_df_IndustrializedUrban <- as.data.frame(result$p)

dir_IndustrializedUrban <- matrix(0, nrow = nrow(corr_df_IndustrializedUrban), ncol = ncol(corr_df_IndustrializedUrban))
rownames(dir_IndustrializedUrban) <- rownames(corr_df_IndustrializedUrban)
colnames(dir_IndustrializedUrban) <- colnames(corr_df_IndustrializedUrban)

for (i in 1:nrow(corr_df_IndustrializedUrban)) {
  for (j in 1:ncol(corr_df_IndustrializedUrban)) {
    corr_val <- corr_df_IndustrializedUrban[i, j]
    p_val <- pval_df_IndustrializedUrban[i, j]
    
    if (!is.na(corr_val) && !is.na(p_val)) {
      if (corr_val > 0 && p_val <= 0.05) {
        dir_IndustrializedUrban[i, j] <- 2
      } else if (corr_val > 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_IndustrializedUrban[i, j] <- 1
      } else if (corr_val < 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_IndustrializedUrban[i, j] <- -1
      } else if (corr_val < 0 && p_val <= 0.05) {
        dir_IndustrializedUrban[i, j] <- -2
      }
    }
  }
}

dir_IndustrializedUrban <- as.data.frame(dir_IndustrializedUrban)

corr_diag <- diag(as.matrix(corr_df_IndustrializedUrban))
pval_diag <- diag(as.matrix(pval_df_IndustrializedUrban))

IndustrializedUrban_corr_16s_wgs <- data.frame(
  corr = corr_diag,
  pval = pval_diag,
  row.names = rownames(corr_df_IndustrializedUrban)
)

library(ggplot2)

IndustrializedUrban_corr_16s_wgs$signif <- ifelse(IndustrializedUrban_corr_16s_wgs$pval <= 0.05, "*", "")
IndustrializedUrban_corr_16s_wgs <- IndustrializedUrban_corr_16s_wgs[-8,]
IndustrializedUrban_corr_16s_wgs$species <- gsub("IndustrializedUrban_","",rownames(IndustrializedUrban_corr_16s_wgs))
IndustrializedUrban_corr_16s_wgs$species <- gsub("_16s","",IndustrializedUrban_corr_16s_wgs$species)

IndustrializedUrban_corr_16s_wgs$species <- factor(IndustrializedUrban_corr_16s_wgs$species, levels = IndustrializedUrban_corr_16s_wgs$species)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\IndustrializedUrban_corr_16s_wgs.pdf", height = 2.5, width = 6)
ggplot(IndustrializedUrban_corr_16s_wgs, aes(x = species, y = corr)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(aes(label = signif, y = corr + 0.02), size = 6) +
  coord_flip() +
  theme_minimal() +
  labs(x = "", y = "Correlation", title = "16S vs WGS Correlation in IndustrializedUrban") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),     
        panel.grid = element_blank(),              
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

################### Others ##############

others_wgs_commonSP_URM_RT <- intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),rownames(df_association_adult_wgs_RuralTribal))
others_wgs_URM_commonSP <- df_association_adult_wgs_UrbanRuralMixed[others_wgs_commonSP_URM_RT,]
others_wgs_RT_commonSP <- df_association_adult_wgs_RuralTribal[others_wgs_commonSP_URM_RT,]
df_association_adult_wgs_Others <- (others_wgs_URM_commonSP + others_wgs_RT_commonSP)/2

others_16s_wgs_commonSP <- intersect(rownames(df_association_adult_wgs_Others),rownames(df_association_adult_16s_Others))

Others_16s_commonSP <- df_association_adult_16s_Others[others_16s_wgs_commonSP,]
Others_wgs_commonSP <- df_association_adult_wgs_Others[others_16s_wgs_commonSP,]

Others_wgs_commonSP <- Others_wgs_commonSP[,1:5]

names(Others_16s_commonSP) <- paste0(names(Others_16s_commonSP),"_16s")
names(Others_wgs_commonSP) <- paste0(names(Others_wgs_commonSP),"_wgs")

n <- ncol(Others_wgs_commonSP)
corr_matrix <- matrix(NA, nrow = n, ncol = n)
pval_matrix <- matrix(NA, nrow = n, ncol = n)

colnames(corr_matrix) <- colnames(Others_16s_commonSP)
rownames(corr_matrix) <- colnames(Others_wgs_commonSP)

colnames(pval_matrix) <- colnames(Others_16s_commonSP)
rownames(pval_matrix) <- colnames(Others_wgs_commonSP)

result <- corr.test(Others_16s_commonSP,Others_wgs_commonSP, method = 'spearman')

corr_df_Others <- as.data.frame(result$r)
pval_df_Others <- as.data.frame(result$p)

dir_Others <- matrix(0, nrow = nrow(corr_df_Others), ncol = ncol(corr_df_Others))
rownames(dir_Others) <- rownames(corr_df_Others)
colnames(dir_Others) <- colnames(corr_df_Others)

for (i in 1:nrow(corr_df_Others)) {
  for (j in 1:ncol(corr_df_Others)) {
    corr_val <- corr_df_Others[i, j]
    p_val <- pval_df_Others[i, j]
    
    if (!is.na(corr_val) && !is.na(p_val)) {
      if (corr_val > 0 && p_val <= 0.05) {
        dir_Others[i, j] <- 2
      } else if (corr_val > 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_Others[i, j] <- 1
      } else if (corr_val < 0 && p_val > 0.05 && p_val <= 0.1) {
        dir_Others[i, j] <- -1
      } else if (corr_val < 0 && p_val <= 0.05) {
        dir_Others[i, j] <- -2
      }
    }
  }
}

dir_Others <- as.data.frame(dir_Others)

corr_diag <- diag(as.matrix(corr_df_Others))
pval_diag <- diag(as.matrix(pval_df_Others))

Others_corr_16s_wgs <- data.frame(
  corr = corr_diag,
  pval = pval_diag,
  row.names = rownames(corr_df_Others)
)

Others_corr_16s_wgs$signif <- ifelse(Others_corr_16s_wgs$pval <= 0.05, "*", "")
Others_corr_16s_wgs$species <- gsub("Others_","",rownames(Others_corr_16s_wgs))
Others_corr_16s_wgs$species <- gsub("_16s","",Others_corr_16s_wgs$species)

Others_corr_16s_wgs$species <- factor(Others_corr_16s_wgs$species, levels = Others_corr_16s_wgs$species)

pdf("G:\\My Drive\\Bifido_Project\\AssociationHeatmaps\\Others_corr_16s_wgs.pdf", height = 2.5, width = 6)
ggplot(Others_corr_16s_wgs, aes(x = species, y = corr)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(aes(label = signif, y = corr + 0.02), size = 6) +
  coord_flip() +
  theme_minimal() +
  labs(x = "", y = "Correlation", title = "16S vs WGS Correlation in Others") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),     
        panel.grid = element_blank(),              
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

#######################################################################
### Supplementary Table 11 and 12 ####

library(openxlsx)
wb <- createWorkbook()

# Get all list names that start with "df_association_infant"
infant_list_names <- ls(pattern = "^df_association_infant")

# Loop through each list name
for (name in infant_list_names) {
  # Get the list object
  list_obj <- get(name)
  
  # Extract the 'association' dataframe
  if ("association" %in% names(list_obj)) {
    assoc_df <- list_obj$association
    
    # Add a worksheet with the same name as the list
    addWorksheet(wb, sheetName = gsub("df_association_infant_","",name))
    
    # Write the dataframe to the sheet
    writeData(wb, sheet = gsub("df_association_infant_","",name), x = assoc_df, rowNames = TRUE)
  }
}

saveWorkbook(wb, file = "G:/My Drive/Bifido_Project/AssociationHeatmaps/infant_association_scores.xlsx", overwrite = TRUE)


wb <- createWorkbook()
target_names <- ls(pattern = "adult_wgs_IndustrializedUrban")
target_names <- target_names[-c(3,11,12)]

for (name in target_names) {
  df <- get(name)
  
  if (is.data.frame(df)) {
    
    addWorksheet(wb, sheetName = gsub("_adult_wgs_IndustrializedUrban","",name))
    writeData(wb, sheet = gsub("_adult_wgs_IndustrializedUrban","",name), x = df, rowNames = TRUE)
  }
}

saveWorkbook(wb, file = "adult_wgs_IndustrializedUrban_AssociationScores.xlsx", overwrite = TRUE)

#### Urban Rural Mixed in WGS ###

wb <- createWorkbook()
target_names <- ls(pattern = "adult_wgs_UrbanRuralMixed")
target_names <- target_names[-c(3,11,12)]

for (name in target_names) {
  df <- get(name)
  
  if (is.data.frame(df)) {
    
    addWorksheet(wb, sheetName = gsub("_adult_wgs_UrbanRuralMixed","",name))
    writeData(wb, sheet = gsub("_adult_wgs_UrbanRuralMixed","",name), x = df, rowNames = TRUE)
  }
}

saveWorkbook(wb, file = "adult_wgs_UrbanRuralMixed_AssociationScores.xlsx", overwrite = TRUE)

wb <- createWorkbook()
target_names <- ls(pattern = "adult_16s_Others")
target_names <- target_names[-c(3,11,12)]

for (name in target_names) {
  df <- get(name)
  
  if (is.data.frame(df)) {
    
    addWorksheet(wb, sheetName = gsub("_adult_16s_Others","",name))
    writeData(wb, sheet = gsub("_adult_16s_Others","",name), x = df, rowNames = TRUE)
  }
}

saveWorkbook(wb, file = "adult_16s_Others_AssociationScores.xlsx", overwrite = TRUE)



