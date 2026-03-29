
library(pheatmap)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

load("/data/all_16s_WGS_combined_NEW.RData")
load("/data/cohort_metadata.RData")

##################################################################
### Function to generate all plots together
##################################################################


generate_association_heatmaps <- function(data_list, prefix, age_group, output_dir = "/results/") {
  
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
    row_names_gp = gpar(fontsize = 4.4),
    column_names_gp = gpar(fontsize = 4.4),
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

generate_association_heatmaps(detection_infant_combined_new, "detection", "infant")
generate_association_heatmaps(detection_adult_combined_new, "detection", "adult")
generate_association_heatmaps(detection_senior_combined_new, "detection", "senior")

# generate_association_heatmaps(longum_infant_combined_new, "longum", "infant")
# generate_association_heatmaps(longum_adult_combined_new, "longum", "adult")
# generate_association_heatmaps(longum_senior_combined_new, "longum", "senior")
# 
# generate_association_heatmaps(adolescentis_infant_combined_new, "adolescentis", "infant")
# generate_association_heatmaps(adolescentis_adult_combined_new, "adolescentis", "adult")
# generate_association_heatmaps(adolescentis_senior_combined_new, "adolescentis", "senior")
# 
# generate_association_heatmaps(pseudocatenulatum_infant_combined_new, "pseudocatenulatum", "infant")
# generate_association_heatmaps(pseudocatenulatum_adult_combined_new, "pseudocatenulatum", "adult")
# generate_association_heatmaps(pseudocatenulatum_senior_combined_new, "pseudocatenulatum", "senior")
# 
# generate_association_heatmaps(bifidum_infant_combined_new, "bifidum", "infant")
# generate_association_heatmaps(bifidum_adult_combined_new, "bifidum", "adult")
# generate_association_heatmaps(bifidum_senior_combined_new, "bifidum", "senior")
# 
# generate_association_heatmaps(breve_infant_combined_new, "breve", "infant")
# generate_association_heatmaps(breve_adult_combined_new, "breve", "adult")
# generate_association_heatmaps(breve_senior_combined_new, "breve", "senior")
# 
# generate_association_heatmaps(dentium_infant_combined_new, "dentium", "infant")
# generate_association_heatmaps(dentium_adult_combined_new, "dentium", "adult")
# generate_association_heatmaps(dentium_senior_combined_new, "dentium", "senior")
# 
# generate_association_heatmaps(catenulatum_infant_combined_new, "catenulatum", "infant")
# generate_association_heatmaps(catenulatum_adult_combined_new, "catenulatum", "adult")
# generate_association_heatmaps(catenulatum_senior_combined_new, "catenulatum", "senior")
# 
# generate_association_heatmaps(animalis_adult_combined_new, "animalis", "adult")
# generate_association_heatmaps(animalis_senior_combined_new, "animalis", "senior")
# 
# 
