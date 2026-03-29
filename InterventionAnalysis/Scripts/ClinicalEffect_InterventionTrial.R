# Clinical Effect on Receptive-Score in Intervention Trial

library(openxlsx)
library(ggplot2)

print("Trial: RamakrishnanM_2025 (adolescentis)") 
load("/data/RamakrishnanM_adolescentis_df.RData") #df_adolescentis

clinical_before <- read.xlsx("/data/ramakrishnan_clinicaldata.xlsx", rowNames = T, sheet = 1)
clinical_during <- read.xlsx("/data/ramakrishnan_clinicaldata.xlsx", rowNames = T, sheet = 2)
clinical_after <- read.xlsx("/data/ramakrishnan_clinicaldata.xlsx", rowNames = T, sheet = 3)

df_clinical_3_1 <- clinical_after - clinical_before

df_clinical_3_1$Receptive_Score <- NULL
df_clinical_3_1$Receptive_Score <- clinical_before$Receptive_Score

# Step 1: sort based on Receptive_Score
df_sorted <- df_clinical_3_1[order(df_clinical_3_1$Receptive_Score), ]
df_filtered <- df_sorted

# Step 2: Identify clinical columns (1 to 6), excluding Receptive_Score
clinical_cols <- c("AbdominalPain", "Bloating", "Diarrhea", "FecalUrgency")

# Step 3: Compute Reduction count per subject
df_filtered$Reduction <- rowSums(df_filtered[, clinical_cols] < 0)

# remove first and last column based on receptive score
cor.test(df_filtered$Reduction, df_filtered$Receptive_Score)

cor.test(df_filtered[2:9, "Reduction"], df_filtered[2:9, "Receptive_Score"])

cor.test(df_filtered[2:9, "AbdominalPain"], df_filtered[2:9, "Receptive_Score"])

cor.test(df_filtered[2:9, "Diarrhea"], df_filtered[2:9, "Receptive_Score"])

cor.test(df_filtered[2:9, "FecalUrgency"], df_filtered[2:9, "Receptive_Score"])

# Create output directory if it doesn't exist
#out_dir <- "RamakrishnanM_ReceptiveScore_ScatterPlots"
#if (!dir.exists(out_dir)) dir.create(out_dir)

# Variables to plot
vars <- c("Reduction", "AbdominalPain", "Bloating", "Diarrhea", "FecalUrgency")

# Plotting function
plot_scatter <- function(var) {
  p <- ggplot(df_filtered[2:9, ], aes_string(x = var, y = "Receptive_Score")) +
    geom_point(
      size = 3.2,
      alpha = 0.8
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 1
    ) +
    labs(
      x = var,
      y = "Receptive Score"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
  
  ggsave(
    filename = paste0("/results/", var, "_vs_ReceptiveScore.pdf"),
    plot = p,
    width = 5,
    height = 4
  )
  
  return(p)
}

# Generate and save plots
plots <- lapply(vars, plot_scatter)

write.xlsx(df_sorted, "/results/df_ramakrishnan_sorted_new.xlsx", rowNames = TRUE)

############################################################################

print("Trial: BaZ_2021 (animalis)") 

supp_data_filt <- read.xlsx("/data/BaZ_supp_data_filt.xlsx", rowNames = T, sheet = 1)

# Extract subject IDs (remove the suffix after "_")
supp_data_filt$Subject_ID <- sub("_(Baseline|PRE|POST|Capsule)$", "", rownames(supp_data_filt))

# Count how many times each subject appears
subject_counts <- table(supp_data_filt$Subject_ID)

# Keep only subjects having all 4 treatments
valid_subjects <- names(subject_counts[subject_counts == 4])

# Filter data for only complete subjects
df_complete <- supp_data_filt[supp_data_filt$Subject_ID %in% valid_subjects, ]

# Split into four separate dfs
df_baseline <- df_complete[grep("_Baseline$", rownames(df_complete)), ]
df_pre      <- df_complete[grep("_PRE$", rownames(df_complete)), ]
df_post     <- df_complete[grep("_POST$", rownames(df_complete)), ]
df_capsule  <- df_complete[grep("_Capsule$", rownames(df_complete)), ]

# Ensure same subject order across all dfs
ordering <- sub("_(Baseline|PRE|POST|Capsule)$", "", rownames(df_baseline))

df_pre     <- df_pre[match(ordering, sub("_(Baseline|PRE|POST|Capsule)$", "", rownames(df_pre))), ]
df_post    <- df_post[match(ordering, sub("_(Baseline|PRE|POST|Capsule)$", "", rownames(df_post))), ]
df_capsule <- df_capsule[match(ordering, sub("_(Baseline|PRE|POST|Capsule)$", "", rownames(df_capsule))), ]

#Convert H to 1 and L to 0 in all 4 dataframes
convert_HL <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.character(x) | is.factor(x)) {
      x <- as.character(x)
      x[x == "H"] <- 1
      x[x == "L"] <- 0
      return(as.numeric(x))
    } else {
      return(x)
    }
  })
  return(df)
}

df_baseline <- convert_HL(df_baseline)
df_pre      <- convert_HL(df_pre)
df_post     <- convert_HL(df_post)
df_capsule  <- convert_HL(df_capsule)
#if you get warning of NAs, relax just Subject_ID has become NA; we are going to remove it

df_baseline$Subject_ID <- NULL
df_pre$Subject_ID <- NULL
df_post$Subject_ID <- NULL
df_capsule$Subject_ID <- NULL

#Clean Rownames Prefix for All Dataframes
clean_rownames <- function(df) {
  rownames(df) <- sub("_(Baseline|PRE|POST|Capsule)$", "", rownames(df))
  return(df)
}

df_baseline <- clean_rownames(df_baseline)
df_pre      <- clean_rownames(df_pre)
df_post     <- clean_rownames(df_post)
df_capsule  <- clean_rownames(df_capsule)


load("/data/BaZ_animalis_df.RData")
#View(df_animalis_3_bppc)
rownames(df_animalis_3_bppc) <- sub("_b$", "", rownames(df_animalis_3_bppc))

intersect(rownames(df_animalis_3_bppc), rownames(df_pre)) #20

common_subjects <- intersect(rownames(df_pre), rownames(df_animalis_3_bppc))

df_animalis_3_bppc_sub <- df_animalis_3_bppc[common_subjects, ]


# do correlation between base and pre based on 20 subjects
#View(df_animalis_3_bppc_sub)
#View(df_baseline)
#View(df_pre)

identical(rownames(df_animalis_3_bppc_sub), rownames(df_baseline))
identical(rownames(df_baseline), rownames(df_pre))

pre_baseline_change <- df_pre - df_baseline
colnames(pre_baseline_change) <- paste0(colnames(pre_baseline_change), "_change")

identical(rownames(df_animalis_3_bppc_sub), rownames(pre_baseline_change))

# Clinical change data
clinical_df <- pre_baseline_change

# Animalis metrics
animalis_df <- df_animalis_3_bppc_sub[, c("animalis_receptive_score", 
                                          "animalis_increase", 
                                          "animalis_base", 
                                          "animalis_ppc")]

# Create empty results dataframe
results <- data.frame(Clinical_Variable = character(),
                      Animalis_Variable = character(),
                      Correlation = numeric(),
                      P_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each clinical column and each animalis metric
for (clinical_col in colnames(clinical_df)) {
  for (animalis_col in colnames(animalis_df)) {
    
    test <- cor.test(clinical_df[[clinical_col]],
                     animalis_df[[animalis_col]], 
                     method = "spearman",
                     exact = TRUE)
    
    results <- rbind(results, data.frame(
      Clinical_Variable = clinical_col,
      Animalis_Variable = animalis_col,
      Correlation = round(test$estimate, 4),
      P_value = round(test$p.value, 4)
    ))
  }
}

# View results
#results

significant_results <- results[complete.cases(results$Correlation, results$P_value) &results$P_value <= 0.05, ]
significant_results

write.xlsx(results, file = "cor_baseline_pre.xlsx", rowNames = TRUE)
write.xlsx(significant_results, file = "sig_cor_baseline_pre.xlsx", rowNames = TRUE)

#scatterplots
# Create output directory
#out_dir <- "BaZ_Animalis_ReceptiveScore_ScatterPlots"
#if (!dir.exists(out_dir)) dir.create(out_dir)

# Clinical variables to plot
clinical_vars <- c(
  "CD3CD8T_change",
  "IL2_change",
  "BB12TNF_change",
  "BB12IL6_change"
)

# Combine into one plotting dataframe (row-aligned)
plot_df <- data.frame(
  animalis_receptive_score = animalis_df$animalis_receptive_score,
  clinical_df[, clinical_vars, drop = FALSE]
)

# Plotting function
plot_scatter <- function(var) {
  p <- ggplot(plot_df, aes_string(
    x = var,
    y = "animalis_receptive_score"
  )) +
    geom_point(
      size = 3.2,
      alpha = 0.8
    ) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 1
    ) +
    labs(
      x = gsub("_change", " change", var),
      y = "*B. animalis* Receptive Score"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
  
  ggsave(
    filename = paste0("/results/", var, "_vs_AnimalisReceptiveScore.pdf"),
    plot = p,
    width = 5,
    height = 4
  )
  
  return(p)
}

# Generate and save all plots
plots <- lapply(clinical_vars, plot_scatter)



















