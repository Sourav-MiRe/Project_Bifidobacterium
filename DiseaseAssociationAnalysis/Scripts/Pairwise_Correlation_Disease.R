#Correlations between the Association-Scores obtained for the different non-Bifidobacterial taxa in the study-sub-cohorts for the six diseases and their controls with respect to the nine Bifidobacterium features
#Figure S16 heatmaps 

library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)

df <- read.xlsx("/data/ordered_carpet_DiseaseAssociation_Heatmap_expanded.xlsx", rowNames = TRUE)
df <- df[!rownames(df) %in% c("pielou", "shannon"), ]

# --- Parse colnames ---
parts <- do.call(rbind, strsplit(colnames(df), "__"))
bif_vec <- parts[,1]
disease_vec <- parts[,2]
col_info <- data.frame(col = colnames(df), bif = bif_vec, disease = disease_vec)

bif_list <- unique(bif_vec)
disease_list <- unique(disease_vec)   # includes Control

# --- Function to compute full disease–disease correlation matrix for one bif ---
make_corr_matrix <- function(bif_name) {
  cols_bif <- col_info$col[col_info$bif == bif_name]
  sub_df <- df[, cols_bif, drop = FALSE]
  
  n <- ncol(sub_df)
  corr_mat <- matrix(NA, n, n, dimnames = list(colnames(sub_df), colnames(sub_df)))
  pval_mat <- matrix(NA, n, n, dimnames = list(colnames(sub_df), colnames(sub_df)))
  
  for(i in 1:n) {
    for(j in 1:n) {
      x <- sub_df[, i]
      y <- sub_df[, j]
      ok <- complete.cases(x, y)
      if(sum(ok) >= 3) {
        t <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
        corr_mat[i, j] <- as.numeric(t$estimate)
        pval_mat[i, j] <- t$p.value
      } else {
        corr_mat[i, j] <- 0     # neutral correlation instead of NA
        pval_mat[i, j] <- 1     # not significant
      }
    }
  }
  
  # Ensure diagonals are always 1 (self-correlation)
  diag(corr_mat) <- 1
  diag(pval_mat) <- 0
  
  # Final cleanup: replace any remaining NA
  corr_mat[is.na(corr_mat)] <- 0
  pval_mat[is.na(pval_mat)] <- 1
  
  list(corr = corr_mat, pval = pval_mat)
}

# --- Color scale ---
col_fun <- colorRamp2(c(-1, 0, 1), c("#D73027", "#FFFFFF", "#1A9850"))

# --- Loop over bifs and make heatmaps ---
for(bif in bif_list) {
  res <- make_corr_matrix(bif)
  corr_mat <- res$corr
  pval_mat <- res$pval
  
  # prepare labels (2 decimals, * if significant)
  labels_mat <- matrix(
    ifelse(!is.na(corr_mat),
           ifelse(pval_mat <= 0.05,
                  sprintf("%.2f*", corr_mat),
                  sprintf("%.2f", corr_mat)),
           ""),
    nrow = nrow(corr_mat),
    dimnames = dimnames(corr_mat)
  )
  
  # --- SAVE OUTPUT TABLES FOR THIS BIF AS XLSX ---
  out_df <- as.data.frame(corr_mat)
  out_df_with_labels <- as.data.frame(labels_mat)
  
  write.xlsx(
    list(
      Correlation_Matrix = out_df,
      Labeled_Matrix = out_df_with_labels,
      Pvalue_Matrix = as.data.frame(pval_mat)
    ),
    paste0("/results/heatmap_table_", bif, ".xlsx"),
    rowNames = TRUE
  )
  
  ht <- Heatmap(
    corr_mat,
    name = "Spearman rho",
    col = col_fun,
    cluster_rows = TRUE,      
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    row_names_side = "left", 
    column_names_side = "top",
    width = unit(ncol(corr_mat) * 13, "mm"),   
    height = unit(nrow(corr_mat) * 13, "mm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width, height,
                gp = gpar(fill = fill, col = "black", lwd = 0.3))
      lab <- labels_mat[i, j]
      if(lab != "") {
        grid.text(lab, x, y, gp = gpar(fontsize = 14, col = "black"))
      }
    }
  )
  
  pdf(paste0("/results/heatmap_", bif, "_disease_corr.pdf"), width = 9, height = 8)
  draw(ht)
  dev.off()
}

