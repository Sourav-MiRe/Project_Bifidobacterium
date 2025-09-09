library(psych)
library(dunn.test)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(MASS)
library(sfsmisc)

load("Intervention_Bifido_trials_data.RData")
load("bif_ml_analysis_scores.RData")

rank_scale=function(x)
{
  x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

train_test_lm <- function(train_df, test_df, target_var, feature1, feature2) {
  formula_str <- as.formula(paste0(target_var, " ~ ", feature1, " + ", feature2))
  lm_model <- lm(formula_str, data = train_df)
  predicted_values <- predict(lm_model, newdata = test_df)
  result_df <- data.frame(Predicted = predicted_values, Actual = test_df[[target_var]], row.names = rownames(test_df))

  corr_result <- corr.test(result_df$Predicted, result_df$Actual, method = "spearman")
  spearman_r <- corr_result$r
  spearman_p <- corr_result$p

  return(list(loocv = result_df, model = lm_model, spearman_r = spearman_r, spearman_p = spearman_p))
}

print("LOOCV Function Read")

#------------------------------------------------------------------------------------
print("Trial 1: Longum Intervention Trial 1")

longum_t0_new <- GomezM_2016_SpeciesProfile[grepl("_t0$", rownames(GomezM_2016_SpeciesProfile)), ]
common_species <- intersect(names(which(apply(longum_t0_new,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_wgs))
df_longum_1 <- data.frame("longum_score"=as.matrix(longum_t0_new[,common_species])%*% ifelse(abs(df_association_all_wgs[common_species,"adult_longum"])>=0,sign(df_association_all_wgs[common_species,"adult_longum"]),0)/sum(ifelse(abs(df_association_all_wgs[common_species,"adult_longum"])>0,sign(df_association_all_wgs[common_species,"adult_longum"]),0))
                          ,"t0"=longum_t0_new[,"Bifidobacterium_longum"])

rownames(df_longum_1) <- gsub("_t0", "", rownames(df_longum_1))

df_longum_1$receptive_score_longum <- rank_scale(df_longum_1[,"longum_score"]) - rank_scale(df_longum_1[,"t0"])

Persisters <- c("A","B","G","H","I","J")
NonPersisters <- setdiff(rownames(df_longum_1),c("A","B","G","H","I","J"))

df_species_association_longum_1 <- data.frame("longum_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_longum"],method="spearman")$r),0)
                                              ,"adolescentis_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),"pseudocatenulatum_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),"bifidum_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_bifidum"],method="spearman")$r),0),"breve_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_breve"],method="spearman")$r),0),"catenulatum_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),"animalis_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_animalis"],method="spearman")$r),0),"dentium_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_dentium"],method="spearman")$r),0),"detection_corr"=ifelse(corr.test(longum_t0_new[,common_species],apply(longum_t0_new[,grep("Bifidobacterium",colnames(longum_t0_new),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],apply(longum_t0_new[,grep("Bifidobacterium",colnames(longum_t0_new),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0))

pdf(file = "Gomez_Boxplot.pdf", width = 5, height = 6)
boxplot(df_longum_1[Persisters,"receptive_score_longum"],df_longum_1[NonPersisters,"receptive_score_longum"], names = c("Persisters", "NonPersisters"), outline=FALSE, col=c("skyblue","bisque"), cex.axis = 1.5)              
dev.off()

wilcox.test(df_longum_1[Persisters,"receptive_score_longum"],df_longum_1[NonPersisters,"receptive_score_longum"])

print("Longum Trial 1 ends")

#------------------------------------------------------------------------------------
print("Trial 2: Animalis Intervention Trial 1")

animalis_1_combined_species_profile <- ZhangQ_2024_SpeciesProfile

common_species <- intersect(names(which(apply(animalis_1_combined_species_profile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_wgs)[df_association_all_wgs[,"adult_animalis"]!=0])

df_species_association_animalis_1  <- data.frame("longum_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_longum"],method="spearman")$r),0),"adolescentis_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),"pseudocatenulatum_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),"bifidum_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),"breve_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_breve"],method="spearman")$r),0),"catenulatum_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),"animalis_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_animalis"],method="spearman")$r),0),"dentium_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_dentium"],method="spearman")$r),0),"detection_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],apply(animalis_1_combined_species_profile[,grep("Bifidobacterium",colnames(animalis_1_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],apply(animalis_1_combined_species_profile[,grep("Bifidobacterium",colnames(animalis_1_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0))

animalis_t0 <- grep("t0",rownames(animalis_1_combined_species_profile),value=TRUE)
animalis_t1 <- grep("t1",rownames(animalis_1_combined_species_profile),value=TRUE)

df_animalis_1 <- data.frame("animalis_score"=as.matrix(animalis_1_combined_species_profile[animalis_t0,common_species])%*% ifelse(abs(df_association_all_wgs[common_species,"adult_animalis"])>0,sign(df_association_all_wgs[common_species,"adult_animalis"]),0),"detection_score"=as.matrix(animalis_1_combined_species_profile[animalis_t0,common_species])%*% ifelse(abs(df_association_all_wgs[common_species,"adult_detection"])>0,sign(df_association_all_wgs[common_species,"adult_detection"]),0),"animalis_t0"=animalis_1_combined_species_profile[animalis_t0,"Bifidobacterium_animalis"],"animalis_t1"=animalis_1_combined_species_profile[animalis_t1,"Bifidobacterium_animalis"],"detection_t0"=apply(animalis_1_combined_species_profile[animalis_t0,grep("Bifidobacterium",colnames(animalis_1_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0.0001]))),"detection_t1"=apply(animalis_1_combined_species_profile[animalis_t1,grep("Bifidobacterium",colnames(animalis_1_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0.0001]))))

df_animalis_1$animalis_receptive_score <- rank_scale(df_animalis_1$animalis_score) - rank_scale(df_animalis_1$animalis_t0)

df_animalis_1$animalis_increase <- rank_scale(df_animalis_1$animalis_t1)-rank_scale(df_animalis_1$animalis_t0)

df_animalis_1_filt <- df_animalis_1[(df_animalis_1$animalis_t0>0)|(df_animalis_1$animalis_t1>0),]

animalis_1_rlm <- rlm(animalis_increase~animalis_receptive_score+animalis_t0,df_animalis_1_filt)
f.robftest(animalis_1_rlm,var="animalis_receptive_score")
animalis_1_rlm[["coefficients"]][["animalis_receptive_score"]]

print("Animalis Trial 1 ends")

#------------------------------------------------------------------------------------
print("Trial 3. Prebiotic Trial \n")

prebio_combined_species_profile <- LooijesteijnE_2024_SpeciesProfile

common_species <- intersect(names(which(apply(prebio_combined_species_profile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_wgs))

prebio_t1 <- grep("_t1",rownames(prebio_combined_species_profile),value=TRUE)
prebio_t2 <- grep("_t2",rownames(prebio_combined_species_profile),value=TRUE)
prebio_t3 <- grep("_t3",rownames(prebio_combined_species_profile),value=TRUE)

df_species_association_prebio <- data.frame("longum_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_longum"],method="spearman")$r),0),"adolescentis_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),"pseudocatenulatum_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),"bifidum_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),"breve_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_breve"],method="spearman")$r),0),"catenulatum_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),"animalis_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_animalis"],method="spearman")$r),0),"dentium_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],prebio_combined_species_profile[,"Bifidobacterium_dentium"],method="spearman")$r),0),"detection_corr"=ifelse(corr.test(prebio_combined_species_profile[,common_species],apply(prebio_combined_species_profile[,grep("Bifidobacterium",colnames(prebio_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(prebio_combined_species_profile[,common_species],apply(prebio_combined_species_profile[,grep("Bifidobacterium",colnames(prebio_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0))

common_species_detection <- intersect(names(which(apply(prebio_combined_species_profile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_wgs)[df_association_all_wgs[,"adult_detection"]>0])

df_prebio_detection <- data.frame("detection_score"=as.matrix(prebio_combined_species_profile[prebio_t2,common_species_detection])%*% ifelse(abs(rowMeans(df_association_all_wgs[common_species_detection,c("adult_detection","senior_detection")]))>0,sign(rowMeans(df_association_all_wgs[common_species_detection,c("adult_detection","senior_detection")])),0),"detection_t0"=apply(prebio_combined_species_profile[prebio_t2,grep("Bifidobacterium",colnames(prebio_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0]))),"detection_t1"=apply(prebio_combined_species_profile[prebio_t3,grep("Bifidobacterium",colnames(prebio_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0]))))

df_prebio_detection$receptive_score <- rank_scale(df_prebio_detection$detection_score) - rank_scale(df_prebio_detection$detection_t0)

df_prebio_detection$increase <- rank_scale(df_prebio_detection$detection_t1) - rank_scale(df_prebio_detection$detection_t0)

prebiotic_rlm <- rlm(increase~receptive_score+detection_t0,df_prebio_detection)
f.robftest(prebiotic_rlm,var="receptive_score")
prebiotic_rlm[["coefficients"]][["receptive_score"]]


print("Prebiotic Trial ends")

#------------------------------------------------------------------------------------
print("Trial 4. Animalis Trial 2\n");

t1 <- grep("Sample_A",grep("_1",rownames(SunB_2022_SpeciesProfile),value=TRUE),value=TRUE)
t2 <- grep("Sample_A",grep("_2",rownames(SunB_2022_SpeciesProfile),value=TRUE),value=TRUE)
t3 <- grep("Sample_A",grep("_3",rownames(SunB_2022_SpeciesProfile),value=TRUE),value=TRUE)

common_species <- intersect(names(which(apply(SunB_2022_SpeciesProfile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_wgs))

SunB_2022_SpeciesProfile[,"Bifidobacterium_catenulatum"] <- 0

df_species_association_animalis_2 <- data.frame("longum_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$r),0),"adolescentis_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),"pseudocatenulatum_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),"bifidum_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),"breve_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$r),0),"catenulatum_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),"animalis_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$r),0),"dentium_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],SunB_2022_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$r),0),"detection_corr"=ifelse(corr.test(SunB_2022_SpeciesProfile[,common_species],apply(SunB_2022_SpeciesProfile[,grep("Bifidobacterium",colnames(SunB_2022_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(SunB_2022_SpeciesProfile[,common_species],apply(SunB_2022_SpeciesProfile[,grep("Bifidobacterium",colnames(SunB_2022_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0))

df_animalis_2 <- data.frame("animalis_score_1"=as.matrix(SunB_2022_SpeciesProfile[t1,common_species])%*% ifelse(abs(df_association_all_wgs[common_species,"senior_animalis"])>0,sign(df_association_all_wgs[common_species,"senior_animalis"]),0),"animalis_score_2"=as.matrix(SunB_2022_SpeciesProfile[t2,common_species])%*% ifelse(abs(df_association_all_wgs[common_species,"senior_animalis"])>0,sign(df_association_all_wgs[common_species,"senior_animalis"]),0),"detection_score_1"=as.matrix(SunB_2022_SpeciesProfile[t1,common_species])%*% df_association_all_wgs[common_species,"senior_detection"],"detection_score_2"=as.matrix(SunB_2022_SpeciesProfile[t2,common_species])%*% df_association_all_wgs[common_species,"senior_detection"],"animalis_t1"=SunB_2022_SpeciesProfile[t1,"Bifidobacterium_animalis"],"animalis_t2"=SunB_2022_SpeciesProfile[t2,"Bifidobacterium_animalis"],"animalis_t3"=SunB_2022_SpeciesProfile[t3,"Bifidobacterium_animalis"],"detection_t1"=apply(SunB_2022_SpeciesProfile[t1,grep("Bifidobacterium",colnames(SunB_2022_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))),"detection_t2"=apply(SunB_2022_SpeciesProfile[t2,grep("Bifidobacterium",colnames(SunB_2022_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))),"detection_t3"=apply(SunB_2022_SpeciesProfile[t3,grep("Bifidobacterium",colnames(SunB_2022_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))))

df_animalis_2$animalis_receptive_score_t1 <- rank_scale(df_animalis_2$animalis_score_1)-rank_scale(df_animalis_2$animalis_t1)

df_animalis_2$animalis_receptive_score_t2 <- rank_scale(df_animalis_2$animalis_score_2)-rank_scale(df_animalis_2$animalis_t2)

df_animalis_2$animalis_increase_t2_t1 <- rank_scale(df_animalis_2$animalis_t2)-rank_scale(df_animalis_2$animalis_t1)

df_animalis_2$animalis_increase_t3_t1 <- rank_scale(df_animalis_2$animalis_t3)-rank_scale(df_animalis_2$animalis_t1)

df_animalis_2$animalis_increase_t3_t2 <- df_animalis_2$animalis_t3-df_animalis_2$animalis_t2

df_animalis_2$detection_receptive_score_t1 <- rank_scale(df_animalis_2$detection_score_1)-rank_scale(df_animalis_2$detection_t1)

df_animalis_2$detection_receptive_score_t2 <- rank_scale(df_animalis_2$detection_score_2)-rank_scale(df_animalis_2$detection_t2)

df_animalis_2$detection_increase_t2_t1 <- df_animalis_2$detection_t2-df_animalis_2$detection_t1

df_animalis_2$detection_increase_t3_t1 <- df_animalis_2$detection_t3-df_animalis_2$detection_t1

df_animalis_2$detection_increase_t3_t2 <- df_animalis_2$detection_t3-df_animalis_2$detection_t2

df_animalis_2_filt <- df_animalis_2[(df_animalis_2$animalis_t1>0)|(df_animalis_2$animalis_t2>0)|(df_animalis_2$animalis_t3>0),]

animalis_2_rlm <- rlm(animalis_increase_t2_t1~animalis_receptive_score_t1+animalis_t1,df_animalis_2_filt)
f.robftest(animalis_2_rlm,var="animalis_receptive_score_t1")
animalis_2_rlm[["coefficients"]][["animalis_receptive_score_t1"]]

print("Trial 4: Animalis Trial 2 ends")

#------------------------------------------------------------------------------------
print("Trial 5: Bifidobacterium Combination MonicaB 2017")

f_fed_m1 <- grep("M1$",grep("\\.ot",grep("\\.1",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
f_fed_m7 <- grep("M7$",grep("\\.ot",grep("\\.1",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
f_fed_m12 <- grep("M12$",grep("\\.ot",grep("\\.1",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
f_fed_m24 <- grep("M24$",grep("\\.ot",grep("\\.1",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)

nf_fed_m1 <- grep("M1$",grep("\\.ot",grep("\\.0",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
nf_fed_m7 <- grep("M7$",grep("\\.ot",grep("\\.0",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
nf_fed_m12 <- grep("M12$",grep("\\.ot",grep("\\.0",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
nf_fed_m24 <- grep("M24$",grep("\\.ot",grep("\\.0",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)

common_species <- intersect(names(which(apply(MonicaB_2017_SpeciesProfile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_16s))

df_species_association_mix_1 <- data.frame("longum_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$r),0),"adolescentis_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),"pseudocatenulatum_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),"bifidum_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),"breve_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$r),0),"catenulatum_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),"animalis_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$r),0),"dentium_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$r),0),"detection_corr"=ifelse(corr.test(MonicaB_2017_SpeciesProfile[,common_species],apply(MonicaB_2017_SpeciesProfile[,grep("Bifidobacterium",colnames(MonicaB_2017_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(MonicaB_2017_SpeciesProfile[,common_species],apply(MonicaB_2017_SpeciesProfile[,grep("Bifidobacterium",colnames(MonicaB_2017_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0))

print("Comparing M1 with M7")
# Taking common subjects across both time-points
f_fed_m1_comp_m1_m7 <- c("100.1.C.B.M1","101.1.N.B.M1","102.1.C.F.M1","103.1.C.BF.M1","105.1.C.F.M1","11.1.C.B.M1","113.1.C.B.M1","114.1.N.F.M1","115.1.N.B.M1","120.1.C.BF.M1","121.1.C.BF.M1","124.1.C.B.M1","125.1.C.B.M1","13.1.N.B.M1","14.1.N.BF.M1","15.1.N.B.M1","19.1.N.B.M1","2.1.C.B.M1","20.1.N.B.M1","25.1.N.F.M1","3.1.N.B.M1","31.1.C.B.M1","36.1.C.F.M1","39.1.C.F.M1","40.1.C.B.M1","41.1.N.B.M1","43.1.C.F.M1","44.1.N.F.M1","46.1.N.F.M1","53.1.N.BF.M1","54.1.N.B.M1","56.1.N.B.M1","58.1.C.B.M1","59.1.N.B.M1","62.1.C.F.M1","65.1.N.B.M1","68.1.N.B.M1","70.1.N.B.M1","73.1.C.B.M1","80.1.C.B.M1","81.1.C.BF.M1","84.1.N.B.M1","89.1.C.B.M1","90.1.N.B.M1")

f_fed_m7_comp_m1_m7 <- c("100.1.C.B.M7","101.1.N.F.M7","102.1.C.F.M7","103.1.C.F.M7","105.1.C.F.M7","11.1.C.B.M7","113.1.C.BF.M7","114.1.N.F.M7","115.1.N.B.M7","120.1.C.BF.M7","121.1.C.F.M7","124.1.C.B.M7","125.1.C.B.M7","13.1.N.F.M7","14.1.N.F.M7","15.1.N.BF.M7","19.1.N.BF.M7","2.1.C.BF.M7","20.1.N.F.M7","25.1.N.F.M7","3.1.N.F.M7","31.1.C.B.M7","36.1.C.F.M7","39.1.C.F.M7","40.1.C.F.M7","41.1.N.BF.M7","43.1.C.F.M7","44.1.N.F.M7","46.1.N.F.M7","53.1.N.B.M7","54.1.N.F.M7","56.1.N.F.M7","58.1.C.F.M7","59.1.N.B.M7","62.1.C.F.M7","65.1.N.F.M7","68.1.N.F.M7","70.1.N.BF.M7","73.1.C.F.M7","80.1.C.B.M7","81.1.C.F.M7","84.1.N.BF.M7","89.1.C.F.M7","90.1.N.B.M7")

df_change_m1_m7 <- data.frame("score_longum_1"=as.matrix(MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,common_species])%*% ifelse(abs(df_association_all_16s[common_species,"infant_longum"])>0,sign(df_association_all_16s[common_species,"infant_longum"]),0),"longum_1"=MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,"Bifidobacterium_longum"],"longum_7"=MonicaB_2017_SpeciesProfile[f_fed_m7_comp_m1_m7,"Bifidobacterium_longum"],"score_breve_1"=as.matrix(MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,common_species])%*% ifelse(abs(df_association_all_16s[common_species,"infant_breve"])>0,sign(df_association_all_16s[common_species,"infant_breve"]),0),"breve_1"=MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,"Bifidobacterium_breve"],"breve_7"=MonicaB_2017_SpeciesProfile[f_fed_m7_comp_m1_m7,"Bifidobacterium_breve"],"score_bifidum_1"=as.matrix(MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,common_species])%*% ifelse(abs(df_association_all_16s[common_species,"infant_bifidum"])>0,sign(df_association_all_16s[common_species,"infant_bifidum"]),0),"bifidum_1"=MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,"Bifidobacterium_bifidum"],"bifidum_7"=MonicaB_2017_SpeciesProfile[f_fed_m7_comp_m1_m7,"Bifidobacterium_bifidum"],"score_detection_1"=as.matrix(MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,common_species])%*% ifelse(abs(df_association_all_16s[common_species,"infant_detection"])>0,sign(df_association_all_16s[common_species,"infant_bifidum"]),0),"combined_1"=apply(MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,c("Bifidobacterium_longum","Bifidobacterium_breve","Bifidobacterium_bifidum")],1,sum),"combined_7"=apply(MonicaB_2017_SpeciesProfile[f_fed_m7_comp_m1_m7,c("Bifidobacterium_longum","Bifidobacterium_breve","Bifidobacterium_bifidum")],1,sum))

df_change_m1_m7$longum_receptive_score <- rank_scale(df_change_m1_m7$score_longum_1) - rank_scale(df_change_m1_m7$longum_1)

df_change_m1_m7$longum_increase <- rank_scale(df_change_m1_m7$longum_7) - rank_scale(df_change_m1_m7$longum_1)

df_change_m1_m7$bifidum_receptive_score <- rank_scale(df_change_m1_m7$score_bifidum_1) - rank_scale(df_change_m1_m7$bifidum_1)

df_change_m1_m7$bifidum_increase <- rank_scale(df_change_m1_m7$bifidum_7) - rank_scale(df_change_m1_m7$bifidum_1)

df_change_m1_m7$breve_receptive_score <- rank_scale(df_change_m1_m7$score_breve_1) - rank_scale(df_change_m1_m7$breve_1)

df_change_m1_m7$breve_increase <- rank_scale(df_change_m1_m7$breve_7) - rank_scale(df_change_m1_m7$breve_1)

df_change_m1_m7$combined_receptive_score <- (df_change_m1_m7$longum_receptive_score + df_change_m1_m7$bifidum_receptive_score + df_change_m1_m7$breve_receptive_score)/3

df_change_m1_m7$combined_increase <- rank_scale(df_change_m1_m7$combined_7) - rank_scale(df_change_m1_m7$combined_1)

mix1_m1_m7_rlm <- rlm(longum_increase~longum_receptive_score+longum_1,df_change_m1_m7)
f.robftest(mix1_m1_m7_rlm,var="longum_receptive_score")
mix1_m1_m7_rlm[["coefficients"]][["longum_receptive_score"]]


mix1_m1_m7_rlm <- rlm(bifidum_increase~bifidum_receptive_score+bifidum_1,df_change_m1_m7)
f.robftest(mix1_m1_m7_rlm,var="bifidum_receptive_score")
mix1_m1_m7_rlm[["coefficients"]][["bifidum_receptive_score"]]


mix1_m1_m7_rlm <- rlm(breve_increase~breve_receptive_score+breve_1,df_change_m1_m7)
f.robftest(mix1_m1_m7_rlm,var="breve_receptive_score")
mix1_m1_m7_rlm[["coefficients"]][["breve_receptive_score"]]

print("Trial 5: Bifidobacterium Combination trial ends")

#------------------------------------------------------------------------------------
print("Trial 6: Animalis Trial 3")

common_species <- intersect(names(which(apply(Baz_2021_SpeciesProfile, 2, function(x) length(x[x > 0])) > 2)),rownames(df_association_all_16s))

df_species_association_animalis_3  <- data.frame(
  "longum_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$r),0),
  "adolescentis_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
  "pseudocatenulatum_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
  "bifidum_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
  "breve_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$r),0),
  "catenulatum_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
  "animalis_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$r),0),
  "dentium_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],Baz_2021_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$r),0),
  "detection_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile[,common_species],apply(Baz_2021_SpeciesProfile[,grep("Bifidobacterium",colnames(Baz_2021_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile[,common_species],apply(Baz_2021_SpeciesProfile[,grep("Bifidobacterium",colnames(Baz_2021_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0))

all_samples <- rownames(Baz_2021_SpeciesProfile)

baseline_samples <- grep("Baseline", all_samples, value = TRUE, ignore.case = TRUE)
pre_samples      <- grep("PRE",      all_samples, value = TRUE, ignore.case = TRUE)
post_samples     <- grep("POST",     all_samples, value = TRUE, ignore.case = TRUE)
capsule_samples  <- grep("Capsule",  all_samples, value = TRUE, ignore.case = TRUE)

baseline_subjects <- gsub("_Baseline", "", baseline_samples, ignore.case = TRUE)
pre_subjects      <- gsub("_PRE",      "", pre_samples,      ignore.case = TRUE)
post_subjects     <- gsub("_POST",     "", post_samples,     ignore.case = TRUE)
capsule_subjects  <- gsub("_Capsule",  "", capsule_samples,  ignore.case = TRUE)

# === 1) Baseline and PRE ===
common_pre <- intersect(baseline_subjects, pre_subjects)
baseline_pre <- Baz_2021_SpeciesProfile[baseline_samples[baseline_subjects %in% common_pre], ]
pre_mat      <- Baz_2021_SpeciesProfile[pre_samples[pre_subjects %in% common_pre], ]
rownames(baseline_pre) <- common_pre
rownames(pre_mat)      <- common_pre

# Baseline and PRE
df_animalis_3_bpre <- data.frame(
  "score_animalis_1" = as.matrix(baseline_pre[, common_species]) %*% ifelse(abs(df_association_all_16s[common_species, "adult_animalis"]) > 0,sign(df_association_all_16s[common_species, "adult_animalis"]), 0),
  "score_detection_1" = as.matrix(baseline_pre[, common_species]) %*% ifelse(abs(df_association_all_16s[common_species, "adult_detection"]) > 0,sign(df_association_all_16s[common_species, "adult_detection"]), 0),
  "animalis_base" = baseline_pre[, "Bifidobacterium_animalis"],
  "animalis_pre"  = pre_mat[, "Bifidobacterium_animalis"],
  "detection_base" = apply(baseline_pre[, grep("Bifidobacterium", colnames(Baz_2021_SpeciesProfile), value = TRUE)],1,function(x) length(x[x > 0.0001])),
  "detection_pre" = apply(pre_mat[, grep("Bifidobacterium", colnames(Baz_2021_SpeciesProfile), value = TRUE)],1,function(x) length(x[x > 0.0001])))

df_animalis_3_bpre$animalis_receptive_score <- rank_scale(df_animalis_3_bpre$score_animalis_1) - rank_scale(df_animalis_3_bpre$animalis_base)

df_animalis_3_bpre$animalis_increase <- rank_scale(df_animalis_3_bpre$animalis_pre)-rank_scale(df_animalis_3_bpre$animalis_base)

df_animalis_3_bpre_filt <- df_animalis_3_bpre[(df_animalis_3_bpre$animalis_base>0)|(df_animalis_3_bpre$animalis_pre>0),]

animalis_3_rlm <- rlm(animalis_increase~animalis_receptive_score+animalis_base,df_animalis_3_bpre_filt)
f.robftest(animalis_3_rlm,var="animalis_receptive_score")
animalis_3_rlm[["coefficients"]][["animalis_receptive_score"]]

# === 2) Baseline and POST ===
common_post <- intersect(baseline_subjects, post_subjects)
baseline_post <- Baz_2021_SpeciesProfile[baseline_samples[baseline_subjects %in% common_post], ]
post_mat      <- Baz_2021_SpeciesProfile[post_samples[post_subjects %in% common_post], ]
rownames(baseline_post) <- common_post
rownames(post_mat)      <- common_post

# Baseline and POST
df_animalis_3_bpost <- data.frame(
  "score_animalis_2" = as.matrix(baseline_post[, common_species]) %*% ifelse(abs(df_association_all_16s[common_species, "adult_animalis"]) > 0,sign(df_association_all_16s[common_species, "adult_animalis"]), 0),
  "score_detection_2" = as.matrix(baseline_post[, common_species]) %*% ifelse(abs(df_association_all_16s[common_species, "adult_detection"]) > 0,sign(df_association_all_16s[common_species, "adult_detection"]), 0),
  "animalis_base" = baseline_post[, "Bifidobacterium_animalis"],
  "animalis_post"  = post_mat[, "Bifidobacterium_animalis"],
  "detection_base" = apply(baseline_post[, grep("Bifidobacterium", colnames(Baz_2021_SpeciesProfile), value = TRUE)],1,function(x) length(x[x > 0.0001])),
  "detection_post" = apply(post_mat[, grep("Bifidobacterium", colnames(Baz_2021_SpeciesProfile), value = TRUE)],1,function(x) length(x[x > 0.0001])))

df_animalis_3_bpost$animalis_receptive_score <- rank_scale(df_animalis_3_bpost$score_animalis_2) - rank_scale(df_animalis_3_bpost$animalis_base)

df_animalis_3_bpost$animalis_increase <- rank_scale(df_animalis_3_bpost$animalis_post)-rank_scale(df_animalis_3_bpost$animalis_base)

df_animalis_3_bpost_filt <- df_animalis_3_bpost[(df_animalis_3_bpost$animalis_base>0)|(df_animalis_3_bpost$animalis_post>0),]

animalis_3_rlm <- rlm(animalis_increase~animalis_receptive_score+animalis_base,df_animalis_3_bpost_filt)
f.robftest(animalis_3_rlm,var="animalis_receptive_score")
animalis_3_rlm[["coefficients"]][["animalis_receptive_score"]]

# === 3) Baseline and Capsule ===
common_cap <- intersect(baseline_subjects, capsule_subjects)
baseline_cap <- Baz_2021_SpeciesProfile[baseline_samples[baseline_subjects %in% common_cap], ]
capsule_mat  <- Baz_2021_SpeciesProfile[capsule_samples[capsule_subjects %in% common_cap], ]
rownames(baseline_cap) <- common_cap
rownames(capsule_mat)  <- common_cap

# Baseline and Capsule
df_animalis_3_bcap <- data.frame(
  "score_animalis_3" = as.matrix(baseline_cap[, common_species]) %*% ifelse(abs(df_association_all_16s[common_species, "adult_animalis"]) > 0,sign(df_association_all_16s[common_species, "adult_animalis"]), 0),
  "score_detection_3" = as.matrix(baseline_cap[, common_species]) %*% ifelse(abs(df_association_all_16s[common_species, "adult_detection"]) > 0,sign(df_association_all_16s[common_species, "adult_detection"]), 0),
  "animalis_base" = baseline_cap[, "Bifidobacterium_animalis"],
  "animalis_cap"  = capsule_mat[, "Bifidobacterium_animalis"],
  "detection_base" = apply(baseline_cap[, grep("Bifidobacterium", colnames(Baz_2021_SpeciesProfile), value = TRUE)],1,function(x) length(x[x > 0.0001])),
  "detection_cap" = apply(capsule_mat[, grep("Bifidobacterium", colnames(Baz_2021_SpeciesProfile), value = TRUE)],1,function(x) length(x[x > 0.0001])))

df_animalis_3_bcap$animalis_receptive_score <- rank_scale(df_animalis_3_bcap$score_animalis_3) - rank_scale(df_animalis_3_bcap$animalis_base)

df_animalis_3_bcap$animalis_increase <- rank_scale(df_animalis_3_bcap$animalis_cap)-rank_scale(df_animalis_3_bcap$animalis_base)

df_animalis_3_bcap_filt <- df_animalis_3_bcap[(df_animalis_3_bcap$animalis_base>0)|(df_animalis_3_bcap$animalis_cap>0),]

animalis_3_rlm <- rlm(animalis_increase~animalis_receptive_score+animalis_base,df_animalis_3_bcap_filt)
f.robftest(animalis_3_rlm,var="animalis_receptive_score")
animalis_3_rlm[["coefficients"]][["animalis_receptive_score"]]

# === 4) Baseline and Mean(PRE, POST, Capsule) ===
common_all <- Reduce(intersect, list(baseline_subjects, pre_subjects, post_subjects, capsule_subjects))

baseline_keep <- baseline_samples[baseline_subjects %in% common_all]
pre_keep      <- pre_samples[pre_subjects %in% common_all]
post_keep     <- post_samples[post_subjects %in% common_all]
capsule_keep  <- capsule_samples[capsule_subjects %in% common_all]

treatment_matrix <- Baz_2021_SpeciesProfile[c(pre_keep, post_keep, capsule_keep), ]

subject_names <- gsub("_(PRE|POST|Capsule)", "", rownames(treatment_matrix), ignore.case = TRUE)

mean_PrePostCap <- aggregate(treatment_matrix, by = list(Subject = subject_names), FUN = median)
anyNA(mean_PrePostCap)

rownames(mean_PrePostCap) <- mean_PrePostCap$Subject
mean_PrePostCap <- mean_PrePostCap[, -1]

baseline_all <- Baz_2021_SpeciesProfile[baseline_keep, ]

baseline_subject_ids <- gsub("_Baseline", "", rownames(baseline_all), ignore.case = TRUE)
rownames(baseline_all) <- baseline_subject_ids

baseline_all <- baseline_all[rownames(mean_PrePostCap), ]

rownames(baseline_all) <- paste0(rownames(baseline_all), "_b")
rownames(mean_PrePostCap) <- paste0(rownames(mean_PrePostCap), "_ppc")

# Baseline and Mean(PRE, POST, Capsule)
df_animalis_3_bppc <- data.frame(
  "score_animalis_4" = as.matrix(baseline_all[, common_species]) %*% ifelse(abs(df_association_all_16s[common_species, "adult_animalis"]) > 0,sign(df_association_all_16s[common_species, "adult_animalis"]), 0),
  "score_detection_4" = as.matrix(baseline_all[, common_species]) %*% ifelse(abs(df_association_all_16s[common_species, "adult_detection"]) > 0,sign(df_association_all_16s[common_species, "adult_detection"]), 0),
  "animalis_base" = baseline_all[, "Bifidobacterium_animalis"],
  "animalis_ppc"  = mean_PrePostCap[, "Bifidobacterium_animalis"],
  "detection_base" = apply(baseline_all[, grep("Bifidobacterium", colnames(Baz_2021_SpeciesProfile), value = TRUE)],1,function(x) length(x[x > 0.0001])),
  "detection_ppc" = apply(mean_PrePostCap[, grep("Bifidobacterium", colnames(Baz_2021_SpeciesProfile), value = TRUE)],1,function(x) length(x[x > 0.0001])))

df_animalis_3_bppc$animalis_receptive_score <- rank_scale(df_animalis_3_bppc$score_animalis_4) - rank_scale(df_animalis_3_bppc$animalis_base)

df_animalis_3_bppc$animalis_increase <- rank_scale(df_animalis_3_bppc$animalis_ppc)-rank_scale(df_animalis_3_bppc$animalis_base)

df_animalis_3_bppc_filt <- df_animalis_3_bppc[(df_animalis_3_bppc$animalis_base>0)|(df_animalis_3_bppc$animalis_ppc>0),]

print("Trial 6: Animalis Trial 3 ends")

#------------------------------------------------------------------------------------
print("Trial 7: Bifidobacterium Breve Trial")

common_species <- intersect(names(which(apply(GronbaekI_2025_SpeciesProfile, 2, function(x) length(x[x > 0])) > 2)),rownames(df_association_all_wgs))

df_species_association_breve  <- data.frame(
  "longum_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$r),0),
  "adolescentis_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
  "pseudocatenulatum_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
  "bifidum_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
  "breve_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$r),0),
  "catenulatum_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
  "animalis_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$r),0),
  "dentium_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],GronbaekI_2025_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$r),0),
  "detection_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],apply(GronbaekI_2025_SpeciesProfile[,grep("Bifidobacterium",colnames(GronbaekI_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile[,common_species],apply(GronbaekI_2025_SpeciesProfile[,grep("Bifidobacterium",colnames(GronbaekI_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0))

#remove Patient38 and Patient61 data
GronbaekI_2025_SpeciesProfile <- GronbaekI_2025_SpeciesProfile[!grepl("Patient38|Patient61", rownames(GronbaekI_2025_SpeciesProfile)), ]

# Define time point sample sets only for Bif195 group
t_week0 <- grep("Bif195__Week0",rownames(GronbaekI_2025_SpeciesProfile),value=TRUE)
t_week4 <- grep("Bif195__Week4",rownames(GronbaekI_2025_SpeciesProfile),value=TRUE)
t_week8 <- grep("Bif195__Week8",rownames(GronbaekI_2025_SpeciesProfile),value=TRUE)
t_week16 <- grep("Bif195__Week16",rownames(GronbaekI_2025_SpeciesProfile),value=TRUE)

df_breve <- data.frame(
  "breve_score_0"=as.matrix(GronbaekI_2025_SpeciesProfile[t_week0,common_species])%*% ifelse(abs(df_association_all_wgs[common_species,"adult_breve"])>0,sign(df_association_all_wgs[common_species,"adult_breve"]),0),
  "breve_score_8"=as.matrix(GronbaekI_2025_SpeciesProfile[t_week8,common_species])%*% ifelse(abs(df_association_all_wgs[common_species,"adult_breve"])>0,sign(df_association_all_wgs[common_species,"adult_breve"]),0),
  "breve_score_16"=as.matrix(GronbaekI_2025_SpeciesProfile[t_week16,common_species])%*% ifelse(abs(df_association_all_wgs[common_species,"adult_breve"])>0,sign(df_association_all_wgs[common_species,"adult_breve"]),0),
  "detection_score_0"=as.matrix(GronbaekI_2025_SpeciesProfile[t_week0,common_species])%*% df_association_all_wgs[common_species,"adult_detection"],
  "detection_score_8"=as.matrix(GronbaekI_2025_SpeciesProfile[t_week8,common_species])%*% df_association_all_wgs[common_species,"adult_detection"],
  "detection_score_16"=as.matrix(GronbaekI_2025_SpeciesProfile[t_week16,common_species])%*% df_association_all_wgs[common_species,"adult_detection"],
  "breve_week0"=GronbaekI_2025_SpeciesProfile[t_week0,"Bifidobacterium_breve"],
  "breve_week4"=GronbaekI_2025_SpeciesProfile[t_week4,"Bifidobacterium_breve"],
  "breve_week8"=GronbaekI_2025_SpeciesProfile[t_week8,"Bifidobacterium_breve"],
  "breve_week16"=GronbaekI_2025_SpeciesProfile[t_week16,"Bifidobacterium_breve"],
  "detection_week0"=apply(GronbaekI_2025_SpeciesProfile[t_week0,grep("Bifidobacterium",colnames(GronbaekI_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))),
  "detection_week4"=apply(GronbaekI_2025_SpeciesProfile[t_week4,grep("Bifidobacterium",colnames(GronbaekI_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))),
  "detection_week8"=apply(GronbaekI_2025_SpeciesProfile[t_week8,grep("Bifidobacterium",colnames(GronbaekI_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))),
  "detection_week16"=apply(GronbaekI_2025_SpeciesProfile[t_week16,grep("Bifidobacterium",colnames(GronbaekI_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))))

df_breve$breve_receptive_score_0 <- rank_scale(df_breve$breve_score_0)-rank_scale(df_breve$breve_week0)

df_breve$breve_receptive_score_8 <- rank_scale(df_breve$breve_score_8)-rank_scale(df_breve$breve_week8)

df_breve$breve_increase_4_0 <- rank_scale(df_breve$breve_week4)-rank_scale(df_breve$breve_week0)

df_breve$breve_increase_8_0 <- rank_scale(df_breve$breve_week8)-rank_scale(df_breve$breve_week0)

df_breve$breve_increase_16_0 <- rank_scale(df_breve$breve_week16)-rank_scale(df_breve$breve_week0)

df_breve$breve_increase_8_4 <- rank_scale(df_breve$breve_week8)-rank_scale(df_breve$breve_week4)

df_breve$breve_increase_16_4 <- rank_scale(df_breve$breve_week16)-rank_scale(df_breve$breve_week4)

df_breve$breve_increase_16_8 <- rank_scale(df_breve$breve_week16)-rank_scale(df_breve$breve_week8)

df_breve_filt <- df_breve[(df_breve$breve_week0>0)|(df_breve$breve_week4>0)|(df_breve$breve_week8>0)|(df_breve$breve_week16>0),]


breve_rlm <- rlm(breve_increase_16_8~breve_receptive_score_8+breve_week8, df_breve_filt)
f.robftest(breve_rlm,var="breve_receptive_score_8")
breve_rlm[["coefficients"]][["breve_receptive_score_8"]]

print("Trial 7: Breve Trial ends")

#------------------------------------------------------------------------------------
print("Trial 8. Adolescentis Trial\n")

probiotic_ids <- rownames(RamakrishnanM_2025_Metadata[RamakrishnanM_2025_Metadata$Treatment == "Probiotic", ])

probiotic_samples <- rownames(RamakrishnanM_2025_SpeciesProfile)[rownames(RamakrishnanM_2025_SpeciesProfile) %in% probiotic_ids]

t1 <- grep("_1$", probiotic_samples, value = TRUE)
t2 <- grep("_2$", probiotic_samples, value = TRUE)
t3 <- grep("_3$", probiotic_samples, value = TRUE)

common_species <- intersect(names(which(apply(RamakrishnanM_2025_SpeciesProfile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_16s))

df_species_association_adolescentis <- data.frame("longum_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$r),0),
                                                  "adolescentis_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
                                                  "pseudocatenulatum_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
                                                  "bifidum_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
                                                  "breve_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$r),0),
                                                  "catenulatum_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
                                                  "animalis_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$r),0),
                                                  "dentium_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],RamakrishnanM_2025_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$r),0),
                                                  "detection_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],apply(RamakrishnanM_2025_SpeciesProfile[,grep("Bifidobacterium",colnames(RamakrishnanM_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile[,common_species],apply(RamakrishnanM_2025_SpeciesProfile[,grep("Bifidobacterium",colnames(RamakrishnanM_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0))

df_adolescentis <- data.frame("adolescentis_score_1"=as.matrix(RamakrishnanM_2025_SpeciesProfile[t1,common_species])%*% ifelse(abs(df_association_all_16s[common_species,"adult_adolescentis"])>0,sign(df_association_all_16s[common_species,"adult_adolescentis"]),0),
                                      "adolescentis_score_2"=as.matrix(RamakrishnanM_2025_SpeciesProfile[t2,common_species])%*% ifelse(abs(df_association_all_16s[common_species,"adult_adolescentis"])>0,sign(df_association_all_16s[common_species,"adult_adolescentis"]),0),
                                      "detection_score_1"=as.matrix(RamakrishnanM_2025_SpeciesProfile[t1,common_species])%*% df_association_all_16s[common_species,"adult_detection"],
                                      "detection_score_2"=as.matrix(RamakrishnanM_2025_SpeciesProfile[t2,common_species])%*% df_association_all_16s[common_species,"adult_detection"],
                                      "adolescentis_t1"=RamakrishnanM_2025_SpeciesProfile[t1,"Bifidobacterium_adolescentis"],
                                      "adolescentis_t2"=RamakrishnanM_2025_SpeciesProfile[t2,"Bifidobacterium_adolescentis"],
                                      "adolescentis_t3"=RamakrishnanM_2025_SpeciesProfile[t3,"Bifidobacterium_adolescentis"],
                                      "detection_t1"=apply(RamakrishnanM_2025_SpeciesProfile[t1,grep("Bifidobacterium",colnames(RamakrishnanM_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))),
                                      "detection_t2"=apply(RamakrishnanM_2025_SpeciesProfile[t2,grep("Bifidobacterium",colnames(RamakrishnanM_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))),
                                      "detection_t3"=apply(RamakrishnanM_2025_SpeciesProfile[t3,grep("Bifidobacterium",colnames(RamakrishnanM_2025_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>=0.0001]))))


df_adolescentis$adolescentis_receptive_score_t1 <- rank_scale(df_adolescentis$adolescentis_score_1)-rank_scale(df_adolescentis$adolescentis_t1)

df_adolescentis$adolescentis_receptive_score_t2 <- rank_scale(df_adolescentis$adolescentis_score_2)-rank_scale(df_adolescentis$adolescentis_t2)

df_adolescentis$adolescentis_increase_t2_t1 <- rank_scale(df_adolescentis$adolescentis_t2)-rank_scale(df_adolescentis$adolescentis_t1)

df_adolescentis$adolescentis_increase_t3_t1 <- rank_scale(df_adolescentis$adolescentis_t3)-rank_scale(df_adolescentis$adolescentis_t1)

df_adolescentis$adolescentis_increase_t3_t2 <- rank_scale(df_adolescentis$adolescentis_t3)-rank_scale(df_adolescentis$adolescentis_t2)

df_adolescentis$adolescentis_increase_t2_t3 <- rank_scale(df_adolescentis$adolescentis_t2)-rank_scale(df_adolescentis$adolescentis_t3)

df_adolescentis$detection_receptive_score_t1 <- rank_scale(df_adolescentis$detection_score_1)-rank_scale(df_adolescentis$detection_t1)

df_adolescentis$detection_receptive_score_t2 <- rank_scale(df_adolescentis$detection_score_2)-rank_scale(df_adolescentis$detection_t2)

df_adolescentis$detection_increase_t2_t1 <- rank_scale(df_adolescentis$detection_t2)-rank_scale(df_adolescentis$detection_t1)

df_adolescentis$detection_increase_t3_t1 <- rank_scale(df_adolescentis$detection_t3)-rank_scale(df_adolescentis$detection_t1)

df_adolescentis$detection_increase_t3_t2 <- rank_scale(df_adolescentis$detection_t3)-rank_scale(df_adolescentis$detection_t2)

df_adolescentis_filt <- df_adolescentis[(df_adolescentis$adolescentis_t1>0)|(df_adolescentis$adolescentis_t2>0)|(df_adolescentis$adolescentis_t3>0),]

adolescentis_rlm <- rlm(adolescentis_increase_t3_t2~adolescentis_receptive_score_t2+adolescentis_t2,df_adolescentis_filt)
f.robftest(adolescentis_rlm,var="adolescentis_receptive_score_t2")
adolescentis_rlm[["coefficients"]][["adolescentis_receptive_score_t2"]]

print("Trial 8: Adolescentis Trial ends")

#------------------------------------------------------------------------------------

#Cross-cohort Validation for B. animalis

subset_df_animalis_1 <- df_animalis_1_filt[, c("animalis_t0", "animalis_receptive_score", "animalis_increase")]
colnames(subset_df_animalis_1) <- c("animalis_baseline", "animalis_receptive_score", "animalis_increase")

subset_df_animalis_2 <- df_animalis_2_filt[, c("animalis_t1", "animalis_receptive_score_t1", "animalis_increase_t2_t1")]
colnames(subset_df_animalis_2) <- c("animalis_baseline", "animalis_receptive_score", "animalis_increase")

subset_df_animalis_3_mean <- df_animalis_3_bppc_filt[, c("animalis_base", "animalis_receptive_score", "animalis_increase")]
colnames(subset_df_animalis_3_mean) <- c("animalis_baseline", "animalis_receptive_score", "animalis_increase")


animalis_ZhangQ_SunB <- train_test_lm(
  subset_df_animalis_1, subset_df_animalis_2,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_SunB_ZhangQ <- train_test_lm(
  subset_df_animalis_2, subset_df_animalis_1,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_ZhangQ_BaZ <- train_test_lm(
  subset_df_animalis_1, subset_df_animalis_3_mean,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_BaZ_ZhangQ <- train_test_lm(
  subset_df_animalis_3_mean, subset_df_animalis_1,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_BaZ_SunB <- train_test_lm(
  subset_df_animalis_3_mean, subset_df_animalis_2,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_SunB_BaZ <- train_test_lm(
  subset_df_animalis_2, subset_df_animalis_3_mean,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

# Scatter plot: actual vs predicted (similarly make all the six-combinations)
ggplot(animalis_SunB_BaZ[["loocv"]], aes(x = Actual, y = Predicted)) +
  geom_point(color = "cornflowerblue", size = 2.5) +
  geom_smooth(method = "lm", color = "black") + 
  labs(
    title = "Actual vs Predicted Increase (SunB -> BaZ)",
    x = "Actual Increase",
    y = "Predicted Increase"
  ) +
  theme_minimal(base_size = 12)









