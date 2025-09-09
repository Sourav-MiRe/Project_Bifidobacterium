#Intervention_analysis_Part2
#Association Correlations between Discovery and Validation Cohorts

##==== wgs adult/senior correlation between Association-Scores of discovery cohort and Association-Scores of validation cohort ====## 
load("all_species_bifs_associations.RData")
load("bif_ml_analysis_scores.RData")

# Function to create combined correlation data frame for any Bifidobacterium species
create_bif_corr_df <- function(df_list, bif_name) {
  corr_col <- paste0(bif_name, "_corr")
  all_species <- unique(unlist(lapply(df_list, rownames)))
  result_list <- list()
  
  for (study in names(df_list)) {
    df <- df_list[[study]]
    if (corr_col %in% colnames(df)) {
      vec <- df[[corr_col]]
      names(vec) <- rownames(df)
    } else {
      vec <- setNames(rep(0, length(rownames(df))), rownames(df))
    }
    full_vec <- setNames(rep(0, length(all_species)), all_species)
    full_vec[names(vec)] <- vec
    result_list[[paste0(bif_name, "_association_", study)]] <- full_vec
  }
  result_df <- as.data.frame(result_list)
  rownames(result_df) <- all_species
  
  return(result_df)
}

df_list <- list(
  GomezM = df_species_association_longum_1,
  ZhangQ = df_species_association_animalis_1,
  LooijesteijnE = df_species_association_prebio,
  SunB = df_species_association_animalis_2,
  GronbaekI = df_species_association_breve)

df_longum_association_intervention <- create_bif_corr_df(df_list, "longum")
df_breve_association_intervention <- create_bif_corr_df(df_list, "breve")
df_bifidum_association_intervention <- create_bif_corr_df(df_list, "bifidum")
df_adolescentis_association_intervention <- create_bif_corr_df(df_list, "adolescentis")
df_animalis_association_intervention <- create_bif_corr_df(df_list, "animalis")
df_dentium_association_intervention <- create_bif_corr_df(df_list, "dentium")
df_catenulatum_association_intervention <- create_bif_corr_df(df_list, "catenulatum")
df_pseudocatenulatum_association_intervention <- create_bif_corr_df(df_list, "pseudocatenulatum")

#Compute Association-Scores for validation cohort 
# Function to compute positive, negative, and total
summarize_bif_association <- function(intervention_df) {
  positive <- rowSums(intervention_df == 1)
  negative <- rowSums(intervention_df == -1)
  total <- ncol(intervention_df)
  
  summary_df <- data.frame(positive = positive, negative = negative,total = total)
  
  return(summary_df)
}

df_longum_association <- summarize_bif_association(df_longum_association_intervention)
df_breve_association <- summarize_bif_association(df_breve_association_intervention)
df_bifidum_association <- summarize_bif_association(df_bifidum_association_intervention)
df_adolescentis_association <- summarize_bif_association(df_adolescentis_association_intervention)
df_animalis_association <- summarize_bif_association(df_animalis_association_intervention)
df_dentium_association <- summarize_bif_association(df_dentium_association_intervention)
df_catenulatum_association <- summarize_bif_association(df_catenulatum_association_intervention)
df_pseudocatenulatum_association <- summarize_bif_association(df_pseudocatenulatum_association_intervention)

# add scores column
#((Positive-Negative)/(Total))*(1-(Min(Positive,Negative)+0.00001)/(Max(Positive,Negative)+0.00001))
df_longum_association$scores <- with(df_longum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_breve_association$scores <- with(df_breve_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_bifidum_association$scores <- with(df_bifidum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_adolescentis_association$scores <- with(df_adolescentis_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_animalis_association$scores <- with(df_animalis_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_dentium_association$scores <- with(df_dentium_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_catenulatum_association$scores <- with(df_catenulatum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_pseudocatenulatum_association$scores <- with(df_pseudocatenulatum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))

#Function to combine adult/senior Association-Scores of Discovery Cohort
combine_mean_scores <- function(adult_df, senior_df) {
  scores_adult <- adult_df[["association"]][, "scores", drop = FALSE]
  scores_senior <- senior_df[["association"]][, "scores", drop = FALSE]
  common_species <- intersect(rownames(scores_adult), rownames(scores_senior))
  scores_adult_common <- scores_adult[common_species, , drop = FALSE]
  scores_senior_common <- scores_senior[common_species, , drop = FALSE]
  merged <- merge(scores_adult_common, scores_senior_common,
                  by = "row.names", suffixes = c("_adult", "_senior"))
  colnames(merged) <- c("species", "adult_score", "senior_score")
  merged$mean_score <- rowMeans(merged[, c("adult_score", "senior_score")], na.rm = TRUE)
  rownames(merged) <- merged$species
  merged <- merged[, setdiff(colnames(merged), "species")]
  
  return(merged)
}

df_association_longum_wgs <- combine_mean_scores(df_association_adult_longum_wgs, df_association_senior_longum_wgs)
df_association_adolescentis_wgs <- combine_mean_scores(df_association_adult_adolescentis_wgs, df_association_senior_adolescentis_wgs)
df_association_animalis_wgs <- combine_mean_scores(df_association_adult_animalis_wgs, df_association_senior_animalis_wgs)
df_association_breve_wgs <- combine_mean_scores(df_association_adult_breve_wgs, df_association_senior_breve_wgs)
df_association_bifidum_wgs <- combine_mean_scores(df_association_adult_bifidum_wgs, df_association_senior_bifidum_wgs)
df_association_catenulatum_wgs <- combine_mean_scores(df_association_adult_catenulatum_wgs, df_association_senior_catenulatum_wgs)
df_association_pseudocatenulatum_wgs <- combine_mean_scores(df_association_adult_pseudocatenulatum_wgs, df_association_senior_pseudocatenulatum_wgs)
df_association_dentium_wgs <- combine_mean_scores(df_association_adult_dentium_wgs, df_association_senior_dentium_wgs)

library(psych)

#Correlation between discovery cohort and validation cohort
compare_discovery_validation <- function(discovery_df, validation_df) {
  common_species <- intersect(rownames(discovery_df), rownames(validation_df))
  discovery_scores <- discovery_df[common_species, "mean_score"]
  validation_scores <- validation_df[common_species, "scores"]
  corr_result <- corr.test(discovery_scores, validation_scores, method = "spearman")
  result_df <- data.frame(discovery_scores = discovery_scores, validation_scores = validation_scores, coefficient = corr_result$r, p_value = corr_result$p)
  rownames(result_df) <- common_species
  
  return(result_df)
}

df_association_longum_cor <- compare_discovery_validation(df_association_longum_wgs, df_longum_association)
df_association_adolescentis_cor <- compare_discovery_validation(df_association_adolescentis_wgs, df_adolescentis_association)
df_association_animalis_cor <- compare_discovery_validation(df_association_animalis_wgs, df_animalis_association)
df_association_breve_cor <- compare_discovery_validation(df_association_breve_wgs, df_breve_association)
df_association_bifidum_cor <- compare_discovery_validation(df_association_bifidum_wgs, df_bifidum_association)
df_association_catenulatum_cor <- compare_discovery_validation(df_association_catenulatum_wgs, df_catenulatum_association)
df_association_pseudocatenulatum_cor <- compare_discovery_validation(df_association_pseudocatenulatum_wgs, df_pseudocatenulatum_association)
df_association_dentium_cor <- compare_discovery_validation(df_association_dentium_wgs, df_dentium_association)

#Bar plot overall coefficient scores
coef_df <- data.frame(
  Species = c("B. longum", "B. adolescentis", "B. pseudocatenulatum", "B. bifidum", "B. breve", "B. dentium", "B. catenulatum", "B. animalis"),
  Coefficient = c(0.4374053, 0.6484365, 0.4395163, 0.3960654, 0.2914569, 0.6695716, 0.4805518, 0.2215258),
  pval = c(4.44E-05, 5.95E-11, 4.04E-05, 0.000252153, 0.008292533, 8.24E-12, 5.61E-06, 0.04686574))

coef_df$signif = ifelse(coef_df$pval < 0.05, "*", "")
coef_df$Species <- factor(coef_df$Species, levels = c("B. longum", "B. adolescentis", "B. pseudocatenulatum", "B. bifidum", "B. breve", "B. dentium", "B. catenulatum", "B. animalis"))

library(ggplot2)
p <-  ggplot(coef_df, aes(x = Species, y = Coefficient)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(aes(label = signif, y = Coefficient + 0.02), size = 6) +
  labs(
    title = "Coefficient Scores per Bifidobacterium Species",
    y = "Spearman's rho", x = ""
  ) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))
ggsave(filename = "overall_correlation.pdf", plot = p, width = 6, height = 2.5, units = "in")

#_________________________________________________________________________________________________________________________
##==== infant monica correlation between association score and bif correlation ====## 
load("MonicaB_2017_SpeciesProfile.RData")
load("bif_ml_analysis_scores.RData")

library(psych)

common_species <- intersect(names(which(apply(MonicaB_2017_SpeciesProfile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_16s))

df_species_association_monica <- data.frame("longum_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$r,
                                            "adolescentis_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r,
                                            "pseudocatenulatum_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r,
                                            "bifidum_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$r,
                                            "breve_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$r,
                                            "catenulatum_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r,
                                            "dentium_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$r)

common_species <- intersect(rownames(df_species_association_monica), rownames(df_association_infant_adolescentis_16s[["association"]]))

x <- df_species_association_monica[common_species, ]

y <- df_association_infant_longum_16s[["association"]][common_species, ]
cor.test(x$longum_corr, y$scores, method = "spearman", exact = T)

y <- df_association_infant_adolescentis_16s[["association"]][common_species, ]
cor.test(x$adolescentis_corr, y$scores, method = "spearman", exact = T)

y <- df_association_infant_pseudocatenulatum_16s[["association"]][common_species, ]
cor.test(x$pseudocatenulatum_corr, y$scores, method = "spearman", exact = T)

y <- df_association_infant_bifidum_16s[["association"]][common_species, ]
cor.test(x$bifidum_corr, y$scores, method = "spearman", exact = T)

y <- df_association_infant_breve_16s[["association"]][common_species, ]
cor.test(x$breve_corr, y$scores, method = "spearman", exact = T)

y <- df_association_infant_dentium_16s[["association"]][common_species, ]
cor.test(x$dentium_corr, y$scores, method = "spearman", exact = T)

y <- df_association_infant_catenulatum_16s[["association"]][common_species, ]
cor.test(x$catenulatum_corr, y$scores, method = "spearman", exact = T)

#infant bar plot overall coefficient scores
coef_df <- data.frame(
  Species = c("B. longum", "B. adolescentis", "B. pseudocatenulatum", "B. bifidum", "B. breve", "B. dentium", "B. catenulatum"),
  Coefficient = c(0.4738098, 0.2752522, 0.2677653, 0.2728982, 0.5162325, 0.5546554, 0.4811543),
  pval = c(1.61E-06, 0.007578, 0.009462, 0.008131, 1.18E-07, 7.99E-09, 1.05E-06))

coef_df$signif = ifelse(coef_df$pval < 0.05, "*", "")
coef_df$Species <- factor(coef_df$Species, levels = c("B. longum", "B. adolescentis", "B. pseudocatenulatum", "B. bifidum", "B. breve", "B. dentium", "B. catenulatum"))

library(ggplot2)
p <- ggplot(coef_df, aes(x = Species, y = Coefficient)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(aes(label = signif, y = Coefficient + 0.02), size = 6) +
  labs(
    title = "Coefficient Scores per Bifidobacterium Species",
    y = "Spearman's rho", x = ""
  ) +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14))
ggsave(filename = "infant_overall_correlation.pdf", plot = p, width = 7, height = 2.5, units = "in")









