#Intervention Analysis

library(psych)
library(dunn.test)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(MASS)
library(sfsmisc)
library(ggplot2)
library(psych)
library(dplyr)
library(purrr)
library(openxlsx)

load("/data/Intervention_Bifido_trials_data.RData")
load("/data/StrongAssociationScores_new_16s.RData") #df_association_all_new_16s
load("/data/StrongAssociationScores_new_wgs.RData") #df_association_all_new_wgs

#rank scale function
rank_scale=function(x)
{
  x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  y <- ifelse(is.nan(y),0,y)
  return(y);
}

#loocv function
loocv_lm <- function(data,var1,var2,var3)
{
  fitted_data <- as.data.frame(matrix(NA,nrow(data),2))
  rownames(fitted_data) <- rownames(data)
  colnames(fitted_data) <- c("Predicted","Actual")
  for(i in 1:nrow(data))
  {
    test_sample <- rownames(data)[i]
    train_data <- data[setdiff(rownames(data),test_sample),]
    temp_lm <- lm(as.formula(paste0(var1,"~",var2,"+",var3)),data=train_data)
    fitted_data[test_sample,"Predicted"] <- as.numeric(predict(temp_lm,data[test_sample,]))
    fitted_data[test_sample,"Actual"] <- data[test_sample,var1]
  }
  lm_model <- lm(as.formula(paste0(var1,"~",var2,"+",var3)),data=data)
  return_list <- list("loocv"=fitted_data,"model"=lm_model)
}


#cross-validation function (training with one dataset)
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


#------------------------------------------------------------------------------------
print("Trial: GomezM_2016")

longum_t0_new <- GomezM_2016_SpeciesProfile[grepl("_t0$", rownames(GomezM_2016_SpeciesProfile)), ]
common_species <- intersect(names(which(apply(longum_t0_new,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_new_wgs))

df_species_association_longum_1 <- data.frame("detection_corr"=ifelse(corr.test(longum_t0_new[,common_species],apply(longum_t0_new[,grep("Bifidobacterium",colnames(longum_t0_new),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],apply(longum_t0_new[,grep("Bifidobacterium",colnames(longum_t0_new),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0),
                                              "longum_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_longum"],method="spearman")$r),0),
                                              "adolescentis_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
                                              "pseudocatenulatum_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
                                              "bifidum_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
                                              "dentium_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_dentium"],method="spearman")$r),0),
                                              "breve_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_breve"],method="spearman")$r),0),
                                              "catenulatum_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
                                              "animalis_corr"=ifelse(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(longum_t0_new[,common_species],longum_t0_new[,"Bifidobacterium_animalis"],method="spearman")$r),0))

df_longum_1 <- data.frame("longum_score"=as.matrix(longum_t0_new[,common_species])%*% ifelse(abs(df_association_all_new_wgs[common_species,"adult_longum"])>=0,sign(df_association_all_new_wgs[common_species,"adult_longum"]),0)/sum(ifelse(abs(df_association_all_new_wgs[common_species,"adult_longum"])>0,sign(df_association_all_new_wgs[common_species,"adult_longum"]),0))
                          ,"t0"=longum_t0_new[,"Bifidobacterium_longum"])

rownames(df_longum_1) <- gsub("_t0", "", rownames(df_longum_1))

df_longum_1$receptive_score_longum <- rank_scale(df_longum_1[,"longum_score"]) - rank_scale(df_longum_1[,"t0"])

Persisters <- c("A","B","G","H","I","J")
NonPersisters <- setdiff(rownames(df_longum_1),c("A","B","G","H","I","J"))

pdf(file = "Gomez_Boxplot.pdf", width = 5, height = 6)
boxplot(df_longum_1[Persisters,"receptive_score_longum"],df_longum_1[NonPersisters,"receptive_score_longum"], names = c("Persisters", "NonPersisters"), outline=FALSE, col=c("skyblue","bisque"), cex.axis = 1.5)              
dev.off()

wilcox.test(df_longum_1[Persisters,"receptive_score_longum"],df_longum_1[NonPersisters,"receptive_score_longum"])

df_longum_1$Persisters <- ifelse(rownames(df_longum_1) %in% Persisters,  1,  0)

#rlm (receptive score)
longum_1_rlm <- rlm(Persisters~receptive_score_longum+t0,df_longum_1)
f.robftest(longum_1_rlm,var="receptive_score_longum")
longum_1_rlm[["coefficients"]][["receptive_score_longum"]]

#loocv (receptive score)
longum_1_post_pred_longum <- loocv_lm(df_longum_1,"Persisters","receptive_score_longum","t0")
corr.test(longum_1_post_pred_longum[["loocv"]][["Predicted"]],longum_1_post_pred_longum[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(longum_1_post_pred_longum[["loocv"]][["Predicted"]],longum_1_post_pred_longum[["loocv"]][["Actual"]],method = "spearman")$p

#------------------------------------------------------------------------------------
print("Trial: ZhangQ_2024 (animalis)")

animalis_1_combined_species_profile <- ZhangQ_2024_SpeciesProfile

common_species <- intersect(names(which(apply(animalis_1_combined_species_profile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_new_wgs)[df_association_all_new_wgs[,"adult_animalis"]!=0])

df_species_association_animalis_1  <- data.frame("detection_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],apply(animalis_1_combined_species_profile[,grep("Bifidobacterium",colnames(animalis_1_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],apply(animalis_1_combined_species_profile[,grep("Bifidobacterium",colnames(animalis_1_combined_species_profile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0),
                                                 "longum_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_longum"],method="spearman")$r),0),
                                                 "adolescentis_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
                                                 "pseudocatenulatum_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
                                                 "bifidum_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
                                                 "dentium_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_dentium"],method="spearman")$r),0),
                                                 "breve_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_breve"],method="spearman")$r),0),
                                                 "catenulatum_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
                                                 "animalis_corr"=ifelse(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(animalis_1_combined_species_profile[,common_species],animalis_1_combined_species_profile[,"Bifidobacterium_animalis"],method="spearman")$r),0))

animalis_t0 <- grep("t0",rownames(animalis_1_combined_species_profile),value=TRUE)
animalis_t1 <- grep("t1",rownames(animalis_1_combined_species_profile),value=TRUE)

df_animalis_1 <- data.frame("animalis_score"=as.matrix(animalis_1_combined_species_profile[animalis_t0,common_species])%*% ifelse(abs(df_association_all_new_wgs[common_species,"adult_animalis"])>0,sign(df_association_all_new_wgs[common_species,"adult_animalis"]),0),
                            "animalis_t0"=animalis_1_combined_species_profile[animalis_t0,"Bifidobacterium_animalis"],
                            "animalis_t1"=animalis_1_combined_species_profile[animalis_t1,"Bifidobacterium_animalis"])

df_animalis_1$animalis_receptive_score <- rank_scale(df_animalis_1$animalis_score) - rank_scale(df_animalis_1$animalis_t0)

df_animalis_1$animalis_increase <- rank_scale(df_animalis_1$animalis_t1)-rank_scale(df_animalis_1$animalis_t0)

#removing if animalis_t0 or animalis_t1 is 0
df_animalis_1_filt <- df_animalis_1[(df_animalis_1$animalis_t0>0)|(df_animalis_1$animalis_t1>0),]

#rlm (receptive score)
animalis_1_rlm <- rlm(animalis_increase~animalis_receptive_score+animalis_t0,df_animalis_1_filt)
f.robftest(animalis_1_rlm,var="animalis_receptive_score")
animalis_1_rlm[["coefficients"]][["animalis_receptive_score"]]

#loocv (receptive score)
animalis_1_post_pred_animalis <- loocv_lm(df_animalis_1_filt,"animalis_increase","animalis_receptive_score","animalis_t0")
corr.test(animalis_1_post_pred_animalis[["loocv"]][["Predicted"]],animalis_1_post_pred_animalis[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(animalis_1_post_pred_animalis[["loocv"]][["Predicted"]],animalis_1_post_pred_animalis[["loocv"]][["Actual"]],method = "spearman")$p

#------------------------------------------------------------------------------------
print("Trial: SunB_2022 (animalis)")

SunB_2022_SpeciesProfile[,"Bifidobacterium_catenulatum"] <- 0

probiotic_SunB_2022_SpeciesProfile <- SunB_2022_SpeciesProfile[grepl("Sample_A", rownames(SunB_2022_SpeciesProfile)), ]

common_species <- intersect(names(which(apply(probiotic_SunB_2022_SpeciesProfile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_new_wgs))

df_species_association_animalis_2_new <- data.frame("detection_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],apply(probiotic_SunB_2022_SpeciesProfile[,grep("Bifidobacterium",colnames(probiotic_SunB_2022_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],apply(probiotic_SunB_2022_SpeciesProfile[,grep("Bifidobacterium",colnames(probiotic_SunB_2022_SpeciesProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0),
                                                    "longum_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$r),0),
                                                    "adolescentis_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
                                                    "pseudocatenulatum_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
                                                    "bifidum_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
                                                    "dentium_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$r),0),
                                                    "breve_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$r),0),
                                                    "catenulatum_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
                                                    "animalis_corr"=ifelse(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(probiotic_SunB_2022_SpeciesProfile[,common_species],probiotic_SunB_2022_SpeciesProfile[,"Bifidobacterium_animalis"],method="spearman")$r),0))

t1 <- grep("Sample_A",grep("_1",rownames(SunB_2022_SpeciesProfile),value=TRUE),value=TRUE)
t2 <- grep("Sample_A",grep("_2",rownames(SunB_2022_SpeciesProfile),value=TRUE),value=TRUE)
t3 <- grep("Sample_A",grep("_3",rownames(SunB_2022_SpeciesProfile),value=TRUE),value=TRUE)

df_animalis_2 <- data.frame("animalis_score_1"=as.matrix(SunB_2022_SpeciesProfile[t1,common_species])%*% ifelse(abs(rowMeans(df_association_all_new_wgs[common_species, c("adult_animalis","senior_animalis")]))>0,sign(rowMeans(df_association_all_new_wgs[common_species,c("adult_animalis","senior_animalis")])),0),
                            "animalis_t1"=SunB_2022_SpeciesProfile[t1,"Bifidobacterium_animalis"],
                            "animalis_t2"=SunB_2022_SpeciesProfile[t2,"Bifidobacterium_animalis"],
                            "animalis_t3"=SunB_2022_SpeciesProfile[t3,"Bifidobacterium_animalis"])

df_animalis_2$animalis_receptive_score_t1 <- rank_scale(df_animalis_2$animalis_score_1)-rank_scale(df_animalis_2$animalis_t1)

df_animalis_2$animalis_increase_t2_t1 <- rank_scale(df_animalis_2$animalis_t2)-rank_scale(df_animalis_2$animalis_t1)

df_animalis_2_filt <- df_animalis_2[(df_animalis_2$animalis_t1>0)|(df_animalis_2$animalis_t2>0)|(df_animalis_2$animalis_t3>0),]

#rlm (receptive score)
animalis_2_rlm <- rlm(animalis_increase_t2_t1~animalis_receptive_score_t1+animalis_t1,df_animalis_2)
f.robftest(animalis_2_rlm,var="animalis_receptive_score_t1")
animalis_2_rlm[["coefficients"]][["animalis_receptive_score_t1"]]

#loocv (receptive score)
animalis_2_post_pred_animalis_t2_t1 <- loocv_lm(df_animalis_2_filt,"animalis_increase_t2_t1","animalis_receptive_score_t1","animalis_t1")
corr.test(animalis_2_post_pred_animalis_t2_t1[["loocv"]][["Predicted"]],animalis_2_post_pred_animalis_t2_t1[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(animalis_2_post_pred_animalis_t2_t1[["loocv"]][["Predicted"]],animalis_2_post_pred_animalis_t2_t1[["loocv"]][["Actual"]],method = "spearman")$p

#------------------------------------------------------------------------------------
print("Trial: MonicaB_2017 (Longum, Bifidum, Breve)")

f_fed_m1 <- grep("M1$",grep("\\.ot",grep("\\.1",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
f_fed_m7 <- grep("M7$",grep("\\.ot",grep("\\.1",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
f_fed_m12 <- grep("M12$",grep("\\.ot",grep("\\.1",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
f_fed_m24 <- grep("M24$",grep("\\.ot",grep("\\.1",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)

nf_fed_m1 <- grep("M1$",grep("\\.ot",grep("\\.0",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
nf_fed_m7 <- grep("M7$",grep("\\.ot",grep("\\.0",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
nf_fed_m12 <- grep("M12$",grep("\\.ot",grep("\\.0",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)
nf_fed_m24 <- grep("M24$",grep("\\.ot",grep("\\.0",rownames(MonicaB_2017_Metadata),value=TRUE),value=TRUE,invert=TRUE),value=TRUE)

length(unique(sub("\\..*$", "", rownames(MonicaB_2017_SpeciesProfile)))) #106

common_species <- intersect(names(which(apply(MonicaB_2017_SpeciesProfile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_new_16s))

df_species_association_mix_1 <- data.frame("longum_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_longum"],method="spearman")$r,
                                           "adolescentis_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r,
                                           "pseudocatenulatum_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r,
                                           "bifidum_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_bifidum"],method="spearman")$r,
                                           "dentium_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_dentium"],method="spearman")$r,
                                           "breve_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_breve"],method="spearman")$r,
                                           "catenulatum_corr"=corr.test(MonicaB_2017_SpeciesProfile[,common_species],MonicaB_2017_SpeciesProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r)

print("Comparing M1 with M7")
# Taking common subjects across both time-points
f_fed_m1_comp_m1_m7 <- c("100.1.C.B.M1","101.1.N.B.M1","102.1.C.F.M1","103.1.C.BF.M1","105.1.C.F.M1","11.1.C.B.M1","113.1.C.B.M1","114.1.N.F.M1","115.1.N.B.M1","120.1.C.BF.M1","121.1.C.BF.M1","124.1.C.B.M1","125.1.C.B.M1","13.1.N.B.M1","14.1.N.BF.M1","15.1.N.B.M1","19.1.N.B.M1","2.1.C.B.M1","20.1.N.B.M1","25.1.N.F.M1","3.1.N.B.M1","31.1.C.B.M1","36.1.C.F.M1","39.1.C.F.M1","40.1.C.B.M1","41.1.N.B.M1","43.1.C.F.M1","44.1.N.F.M1","46.1.N.F.M1","53.1.N.BF.M1","54.1.N.B.M1","56.1.N.B.M1","58.1.C.B.M1","59.1.N.B.M1","62.1.C.F.M1","65.1.N.B.M1","68.1.N.B.M1","70.1.N.B.M1","73.1.C.B.M1","80.1.C.B.M1","81.1.C.BF.M1","84.1.N.B.M1","89.1.C.B.M1","90.1.N.B.M1")

f_fed_m7_comp_m1_m7 <- c("100.1.C.B.M7","101.1.N.F.M7","102.1.C.F.M7","103.1.C.F.M7","105.1.C.F.M7","11.1.C.B.M7","113.1.C.BF.M7","114.1.N.F.M7","115.1.N.B.M7","120.1.C.BF.M7","121.1.C.F.M7","124.1.C.B.M7","125.1.C.B.M7","13.1.N.F.M7","14.1.N.F.M7","15.1.N.BF.M7","19.1.N.BF.M7","2.1.C.BF.M7","20.1.N.F.M7","25.1.N.F.M7","3.1.N.F.M7","31.1.C.B.M7","36.1.C.F.M7","39.1.C.F.M7","40.1.C.F.M7","41.1.N.BF.M7","43.1.C.F.M7","44.1.N.F.M7","46.1.N.F.M7","53.1.N.B.M7","54.1.N.F.M7","56.1.N.F.M7","58.1.C.F.M7","59.1.N.B.M7","62.1.C.F.M7","65.1.N.F.M7","68.1.N.F.M7","70.1.N.BF.M7","73.1.C.F.M7","80.1.C.B.M7","81.1.C.F.M7","84.1.N.BF.M7","89.1.C.F.M7","90.1.N.B.M7")

df_change_m1_m7 <- data.frame("score_longum_1"=as.matrix(MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,common_species])%*% ifelse(abs(df_association_all_new_16s[common_species,"infant_longum"])>0,sign(df_association_all_new_16s[common_species,"infant_longum"]),0),
                              "longum_1"=MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,"Bifidobacterium_longum"],
                              "longum_7"=MonicaB_2017_SpeciesProfile[f_fed_m7_comp_m1_m7,"Bifidobacterium_longum"],
                              "score_breve_1"=as.matrix(MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,common_species])%*% ifelse(abs(df_association_all_new_16s[common_species,"infant_breve"])>0,sign(df_association_all_new_16s[common_species,"infant_breve"]),0),
                              "breve_1"=MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,"Bifidobacterium_breve"],
                              "breve_7"=MonicaB_2017_SpeciesProfile[f_fed_m7_comp_m1_m7,"Bifidobacterium_breve"],
                              "score_bifidum_1"=as.matrix(MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,common_species])%*% ifelse(abs(df_association_all_new_16s[common_species,"infant_bifidum"])>0,sign(df_association_all_new_16s[common_species,"infant_bifidum"]),0),
                              "bifidum_1"=MonicaB_2017_SpeciesProfile[f_fed_m1_comp_m1_m7,"Bifidobacterium_bifidum"],
                              "bifidum_7"=MonicaB_2017_SpeciesProfile[f_fed_m7_comp_m1_m7,"Bifidobacterium_bifidum"])

df_change_m1_m7$longum_receptive_score <- rank_scale(df_change_m1_m7$score_longum_1) - rank_scale(df_change_m1_m7$longum_1)

df_change_m1_m7$longum_increase <- rank_scale(df_change_m1_m7$longum_7) - rank_scale(df_change_m1_m7$longum_1)

df_change_m1_m7$bifidum_receptive_score <- rank_scale(df_change_m1_m7$score_bifidum_1) - rank_scale(df_change_m1_m7$bifidum_1)

df_change_m1_m7$bifidum_increase <- rank_scale(df_change_m1_m7$bifidum_7) - rank_scale(df_change_m1_m7$bifidum_1)

df_change_m1_m7$breve_receptive_score <- rank_scale(df_change_m1_m7$score_breve_1) - rank_scale(df_change_m1_m7$breve_1)

df_change_m1_m7$breve_increase <- rank_scale(df_change_m1_m7$breve_7) - rank_scale(df_change_m1_m7$breve_1)

## longum
#rlm (recptive score)
mix1_m1_m7_rlm <- rlm(longum_increase~longum_receptive_score+longum_1,df_change_m1_m7)
f.robftest(mix1_m1_m7_rlm,var="longum_receptive_score")
mix1_m1_m7_rlm[["coefficients"]][["longum_receptive_score"]]

#loocv (receptive score)
mix_1_post_pred_longum_m1_m7 <- loocv_lm(df_change_m1_m7,"longum_increase","longum_receptive_score","longum_1")
corr.test(mix_1_post_pred_longum_m1_m7[["loocv"]][["Predicted"]],mix_1_post_pred_longum_m1_m7[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(mix_1_post_pred_longum_m1_m7[["loocv"]][["Predicted"]],mix_1_post_pred_longum_m1_m7[["loocv"]][["Actual"]],method = "spearman")$p

## bifidum
#rlm (receptive score)
mix1_m1_m7_rlm <- rlm(bifidum_increase~bifidum_receptive_score+bifidum_1,df_change_m1_m7)
f.robftest(mix1_m1_m7_rlm,var="bifidum_receptive_score")
mix1_m1_m7_rlm[["coefficients"]][["bifidum_receptive_score"]]

#loocv (receptive score)
mix_1_post_pred_bifidum_m1_m7 <- loocv_lm(df_change_m1_m7,"bifidum_increase","bifidum_receptive_score","bifidum_1")
corr.test(mix_1_post_pred_bifidum_m1_m7[["loocv"]][["Predicted"]],mix_1_post_pred_bifidum_m1_m7[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(mix_1_post_pred_bifidum_m1_m7[["loocv"]][["Predicted"]],mix_1_post_pred_bifidum_m1_m7[["loocv"]][["Actual"]],method = "spearman")$p

## breve
#rlm (receptice score)
mix1_m1_m7_rlm <- rlm(breve_increase~breve_receptive_score+breve_1,df_change_m1_m7)
f.robftest(mix1_m1_m7_rlm,var="breve_receptive_score")
mix1_m1_m7_rlm[["coefficients"]][["breve_receptive_score"]]

#loocv (receptive score)
mix_1_post_pred_breve_m1_m7 <- loocv_lm(df_change_m1_m7,"breve_increase","breve_receptive_score","breve_1")
corr.test(mix_1_post_pred_breve_m1_m7[["loocv"]][["Predicted"]],mix_1_post_pred_breve_m1_m7[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(mix_1_post_pred_breve_m1_m7[["loocv"]][["Predicted"]],mix_1_post_pred_breve_m1_m7[["loocv"]][["Actual"]],method = "spearman")$p

#------------------------------------------------------------------------------------
print("Trial: BaZ_2021 (animalis)")
table(gsub(".*(Baseline|PRE|POST|Capsule).*", "\\1", rownames(Baz_2021_SpeciesProfile)))
Baz_2021_SpeciesProfile_sub <- Baz_2021_SpeciesProfile[grepl("Baseline|PRE|POST|Capsule", rownames(Baz_2021_SpeciesProfile)),]
table(gsub(".*(Baseline|PRE|POST|Capsule).*", "\\1", rownames(Baz_2021_SpeciesProfile_sub)))

common_species <- intersect(names(which(apply(Baz_2021_SpeciesProfile_sub, 2, function(x) length(x[x > 0])) > 2)),rownames(df_association_all_new_16s))

df_species_association_animalis_3_new  <- data.frame(
  "detection_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],apply(Baz_2021_SpeciesProfile_sub[,grep("Bifidobacterium",colnames(Baz_2021_SpeciesProfile_sub),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],apply(Baz_2021_SpeciesProfile_sub[,grep("Bifidobacterium",colnames(Baz_2021_SpeciesProfile_sub),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0),
  "longum_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_longum"],method="spearman")$r),0),
  "adolescentis_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
  "pseudocatenulatum_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
  "bifidum_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
  "dentium_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_dentium"],method="spearman")$r),0),
  "breve_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_breve"],method="spearman")$r),0),
  "catenulatum_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
  "animalis_corr"=ifelse(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(Baz_2021_SpeciesProfile_sub[,common_species],Baz_2021_SpeciesProfile_sub[,"Bifidobacterium_animalis"],method="spearman")$r),0))

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
  "score_animalis_1" = as.matrix(baseline_pre[, common_species]) %*% ifelse(abs(df_association_all_new_16s[common_species, "adult_animalis"]) > 0,sign(df_association_all_new_16s[common_species, "adult_animalis"]), 0),
  "animalis_base" = baseline_pre[, "Bifidobacterium_animalis"],
  "animalis_pre"  = pre_mat[, "Bifidobacterium_animalis"])

df_animalis_3_bpre$animalis_receptive_score <- rank_scale(df_animalis_3_bpre$score_animalis_1) - rank_scale(df_animalis_3_bpre$animalis_base)

df_animalis_3_bpre$animalis_increase <- rank_scale(df_animalis_3_bpre$animalis_pre)-rank_scale(df_animalis_3_bpre$animalis_base)

df_animalis_3_bpre_filt <- df_animalis_3_bpre[(df_animalis_3_bpre$animalis_base>0)|(df_animalis_3_bpre$animalis_pre>0),]

#rlm (receptive score)
animalis_3_rlm <- rlm(animalis_increase~animalis_receptive_score+animalis_base,df_animalis_3_bpre_filt)
f.robftest(animalis_3_rlm,var="animalis_receptive_score")
animalis_3_rlm[["coefficients"]][["animalis_receptive_score"]]

#loocv (receptive score)
animalis_3_post_pred_animalis_bpre <- loocv_lm(df_animalis_3_bpre_filt,"animalis_increase","animalis_receptive_score","animalis_base")
corr.test(animalis_3_post_pred_animalis_bpre[["loocv"]][["Predicted"]],animalis_3_post_pred_animalis_bpre[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(animalis_3_post_pred_animalis_bpre[["loocv"]][["Predicted"]],animalis_3_post_pred_animalis_bpre[["loocv"]][["Actual"]],method = "spearman")$p

# === 2) Baseline and POST ===
common_post <- intersect(baseline_subjects, post_subjects)
baseline_post <- Baz_2021_SpeciesProfile[baseline_samples[baseline_subjects %in% common_post], ]
post_mat      <- Baz_2021_SpeciesProfile[post_samples[post_subjects %in% common_post], ]
rownames(baseline_post) <- common_post
rownames(post_mat)      <- common_post

# Baseline and POST
df_animalis_3_bpost <- data.frame(
  "score_animalis_2" = as.matrix(baseline_post[, common_species]) %*% ifelse(abs(df_association_all_new_16s[common_species, "adult_animalis"]) > 0,sign(df_association_all_new_16s[common_species, "adult_animalis"]), 0),
  "animalis_base" = baseline_post[, "Bifidobacterium_animalis"],
  "animalis_post"  = post_mat[, "Bifidobacterium_animalis"])

df_animalis_3_bpost$animalis_receptive_score <- rank_scale(df_animalis_3_bpost$score_animalis_2) - rank_scale(df_animalis_3_bpost$animalis_base)

df_animalis_3_bpost$animalis_increase <- rank_scale(df_animalis_3_bpost$animalis_post)-rank_scale(df_animalis_3_bpost$animalis_base)

df_animalis_3_bpost_filt <- df_animalis_3_bpost[(df_animalis_3_bpost$animalis_base>0)|(df_animalis_3_bpost$animalis_post>0),]

#rlm (receptive score)
animalis_3_rlm <- rlm(animalis_increase~animalis_receptive_score+animalis_base,df_animalis_3_bpost_filt)
f.robftest(animalis_3_rlm,var="animalis_receptive_score")
animalis_3_rlm[["coefficients"]][["animalis_receptive_score"]]

#loocv (receptive score)
animalis_3_post_pred_animalis_bpost <- loocv_lm(df_animalis_3_bpost_filt,"animalis_increase","animalis_receptive_score","animalis_base")
corr.test(animalis_3_post_pred_animalis_bpost[["loocv"]][["Predicted"]],animalis_3_post_pred_animalis_bpost[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(animalis_3_post_pred_animalis_bpost[["loocv"]][["Predicted"]],animalis_3_post_pred_animalis_bpost[["loocv"]][["Actual"]],method = "spearman")$p

# === 3) Baseline and Capsule ===
common_cap <- intersect(baseline_subjects, capsule_subjects)
baseline_cap <- Baz_2021_SpeciesProfile[baseline_samples[baseline_subjects %in% common_cap], ]
capsule_mat  <- Baz_2021_SpeciesProfile[capsule_samples[capsule_subjects %in% common_cap], ]
rownames(baseline_cap) <- common_cap
rownames(capsule_mat)  <- common_cap

# Baseline and Capsule
df_animalis_3_bcap <- data.frame(
  "score_animalis_3" = as.matrix(baseline_cap[, common_species]) %*% ifelse(abs(df_association_all_new_16s[common_species, "adult_animalis"]) > 0,sign(df_association_all_new_16s[common_species, "adult_animalis"]), 0),
  "animalis_base" = baseline_cap[, "Bifidobacterium_animalis"],
  "animalis_cap"  = capsule_mat[, "Bifidobacterium_animalis"])

df_animalis_3_bcap$animalis_receptive_score <- rank_scale(df_animalis_3_bcap$score_animalis_3) - rank_scale(df_animalis_3_bcap$animalis_base)

df_animalis_3_bcap$animalis_increase <- rank_scale(df_animalis_3_bcap$animalis_cap)-rank_scale(df_animalis_3_bcap$animalis_base)

df_animalis_3_bcap_filt <- df_animalis_3_bcap[(df_animalis_3_bcap$animalis_base>0)|(df_animalis_3_bcap$animalis_cap>0),]

#rlm (receptive score)
animalis_3_rlm <- rlm(animalis_increase~animalis_receptive_score+animalis_base,df_animalis_3_bcap_filt)
f.robftest(animalis_3_rlm,var="animalis_receptive_score")
animalis_3_rlm[["coefficients"]][["animalis_receptive_score"]]

#loocv (receptive score)
animalis_3_post_pred_animalis_bcap <- loocv_lm(df_animalis_3_bcap_filt,"animalis_increase","animalis_receptive_score","animalis_base")
corr.test(animalis_3_post_pred_animalis_bcap[["loocv"]][["Predicted"]],animalis_3_post_pred_animalis_bcap[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(animalis_3_post_pred_animalis_bcap[["loocv"]][["Predicted"]],animalis_3_post_pred_animalis_bcap[["loocv"]][["Actual"]],method = "spearman")$p

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
  "score_animalis_4" = as.matrix(baseline_all[, common_species]) %*% ifelse(abs(df_association_all_new_16s[common_species, "adult_animalis"]) > 0,sign(df_association_all_new_16s[common_species, "adult_animalis"]), 0),
  "animalis_base" = baseline_all[, "Bifidobacterium_animalis"],
  "animalis_ppc"  = mean_PrePostCap[, "Bifidobacterium_animalis"])

df_animalis_3_bppc$animalis_receptive_score <- rank_scale(df_animalis_3_bppc$score_animalis_4) - rank_scale(df_animalis_3_bppc$animalis_base)

df_animalis_3_bppc$animalis_increase <- rank_scale(df_animalis_3_bppc$animalis_ppc)-rank_scale(df_animalis_3_bppc$animalis_base)

df_animalis_3_bppc_filt <- df_animalis_3_bppc[(df_animalis_3_bppc$animalis_base>0)|(df_animalis_3_bppc$animalis_ppc>0),]

#------------------------------------------------------------------------------------
print("Trial: GronbaekI_2025 (breve)")

GronbaekI_2025_SpeciesProfile <- GronbaekI_2025_SpeciesProfile[!grepl("Patient38|Patient61", rownames(GronbaekI_2025_SpeciesProfile)), ]

table(gsub(".*(Week0|Week4|Week8|Week16).*", "\\1", rownames(GronbaekI_2025_SpeciesProfile)))
GronbaekI_2025_SpeciesProfile_sub <- GronbaekI_2025_SpeciesProfile[grepl("Week0|Week4|Week8|Week16", rownames(GronbaekI_2025_SpeciesProfile)),]
table(gsub(".*(Week0|Week4|Week8|Week16).*", "\\1", rownames(GronbaekI_2025_SpeciesProfile_sub)))

common_species <- intersect(names(which(apply(GronbaekI_2025_SpeciesProfile_sub, 2, function(x) length(x[x > 0])) > 2)),rownames(df_association_all_new_wgs))

df_species_association_breve_1_new  <- data.frame(
  "detection_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],apply(GronbaekI_2025_SpeciesProfile_sub[,grep("Bifidobacterium",colnames(GronbaekI_2025_SpeciesProfile_sub),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],apply(GronbaekI_2025_SpeciesProfile_sub[,grep("Bifidobacterium",colnames(GronbaekI_2025_SpeciesProfile_sub),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0),
  "longum_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_longum"],method="spearman")$r),0),
  "adolescentis_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
  "pseudocatenulatum_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
  "bifidum_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
  "dentium_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_dentium"],method="spearman")$r),0),
  "breve_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_breve"],method="spearman")$r),0),
  "catenulatum_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
  "animalis_corr"=ifelse(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(GronbaekI_2025_SpeciesProfile_sub[,common_species],GronbaekI_2025_SpeciesProfile_sub[,"Bifidobacterium_animalis"],method="spearman")$r),0))

# Define time point sample sets only for Bif195 group
t_week0 <- grep("Bif195__Week0",rownames(GronbaekI_2025_SpeciesProfile),value=TRUE)
t_week4 <- grep("Bif195__Week4",rownames(GronbaekI_2025_SpeciesProfile),value=TRUE)
t_week8 <- grep("Bif195__Week8",rownames(GronbaekI_2025_SpeciesProfile),value=TRUE)
t_week16 <- grep("Bif195__Week16",rownames(GronbaekI_2025_SpeciesProfile),value=TRUE)

df_breve <- data.frame(
  "breve_score_8"=as.matrix(GronbaekI_2025_SpeciesProfile[t_week8,common_species])%*% ifelse(abs(rowMeans(df_association_all_new_wgs[common_species, c("adult_breve", "senior_breve")]))>0,sign(rowMeans(df_association_all_new_wgs[common_species,c("adult_breve", "senior_breve")])),0),
  "breve_week0"=GronbaekI_2025_SpeciesProfile[t_week0,"Bifidobacterium_breve"],
  "breve_week4"=GronbaekI_2025_SpeciesProfile[t_week4,"Bifidobacterium_breve"],
  "breve_week8"=GronbaekI_2025_SpeciesProfile[t_week8,"Bifidobacterium_breve"],
  "breve_week16"=GronbaekI_2025_SpeciesProfile[t_week16,"Bifidobacterium_breve"])

df_breve$breve_receptive_score_8 <- rank_scale(df_breve$breve_score_8)-rank_scale(df_breve$breve_week8)

df_breve$breve_increase_16_8 <- rank_scale(df_breve$breve_week16)-rank_scale(df_breve$breve_week8)

df_breve_filt <- df_breve[(df_breve$breve_week0>0)|(df_breve$breve_week4>0)|(df_breve$breve_week8>0)|(df_breve$breve_week16>0),]

#rlm (receptive score)
breve_rlm <- rlm(breve_increase_16_8~breve_receptive_score_8+breve_week8, df_breve_filt)
f.robftest(breve_rlm,var="breve_receptive_score_8")
breve_rlm[["coefficients"]][["breve_receptive_score_8"]]

#loocv (receptive score)
breve_post_pred_breve_16_8 <- loocv_lm(df_breve,"breve_increase_16_8","breve_receptive_score_8","breve_week8")
corr.test(breve_post_pred_breve_16_8[["loocv"]][["Predicted"]],breve_post_pred_breve_16_8[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(breve_post_pred_breve_16_8[["loocv"]][["Predicted"]],breve_post_pred_breve_16_8[["loocv"]][["Actual"]],method = "spearman")$p

#------------------------------------------------------------------------------------
print("Trial: RamakrishnanM_2025 (adolescentis)")

table(RamakrishnanM_2025_Metadata$Treatment)
RamakrishnanM_2025_SpeciesProfile_probiotic <- RamakrishnanM_2025_SpeciesProfile[rownames(RamakrishnanM_2025_SpeciesProfile) %in% rownames(RamakrishnanM_2025_Metadata[RamakrishnanM_2025_Metadata$Treatment == "Probiotic", ]),]

common_species <- intersect(names(which(apply(RamakrishnanM_2025_SpeciesProfile_probiotic,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_new_16s))

df_species_association_adolescentis_1_new <- data.frame("detection_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],apply(RamakrishnanM_2025_SpeciesProfile_probiotic[,grep("Bifidobacterium",colnames(RamakrishnanM_2025_SpeciesProfile_probiotic),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],apply(RamakrishnanM_2025_SpeciesProfile_probiotic[,grep("Bifidobacterium",colnames(RamakrishnanM_2025_SpeciesProfile_probiotic),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0),
                                                        "longum_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_longum"],method="spearman")$r),0),
                                                        "adolescentis_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
                                                        "pseudocatenulatum_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
                                                        "bifidum_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
                                                        "dentium_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_dentium"],method="spearman")$r),0),
                                                        "breve_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_breve"],method="spearman")$r),0),
                                                        "catenulatum_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
                                                        "animalis_corr"=ifelse(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(RamakrishnanM_2025_SpeciesProfile_probiotic[,common_species],RamakrishnanM_2025_SpeciesProfile_probiotic[,"Bifidobacterium_animalis"],method="spearman")$r),0))

probiotic_ids <- rownames(RamakrishnanM_2025_Metadata[RamakrishnanM_2025_Metadata$Treatment == "Probiotic", ])

probiotic_samples <- rownames(RamakrishnanM_2025_SpeciesProfile)[rownames(RamakrishnanM_2025_SpeciesProfile) %in% probiotic_ids]

t1 <- grep("_1$", probiotic_samples, value = TRUE)
t2 <- grep("_2$", probiotic_samples, value = TRUE)
t3 <- grep("_3$", probiotic_samples, value = TRUE)

df_adolescentis <- data.frame("adolescentis_score_2"=as.matrix(RamakrishnanM_2025_SpeciesProfile[t2,common_species])%*% ifelse(abs(df_association_all_new_16s[common_species,"adult_adolescentis"])>0,sign(df_association_all_new_16s[common_species,"adult_adolescentis"]),0),
                              "adolescentis_t1"=RamakrishnanM_2025_SpeciesProfile[t1,"Bifidobacterium_adolescentis"],
                              "adolescentis_t2"=RamakrishnanM_2025_SpeciesProfile[t2,"Bifidobacterium_adolescentis"],
                              "adolescentis_t3"=RamakrishnanM_2025_SpeciesProfile[t3,"Bifidobacterium_adolescentis"])

df_adolescentis$adolescentis_receptive_score_t2 <- rank_scale(df_adolescentis$adolescentis_score_2)-rank_scale(df_adolescentis$adolescentis_t2)

df_adolescentis$adolescentis_increase_t3_t2 <- rank_scale(df_adolescentis$adolescentis_t3)-rank_scale(df_adolescentis$adolescentis_t2)

df_adolescentis_filt <- df_adolescentis[(df_adolescentis$adolescentis_t1>0)|(df_adolescentis$adolescentis_t2>0)|(df_adolescentis$adolescentis_t3>0),]

#rlm (receptive score)
adolescentis_rlm <- rlm(adolescentis_increase_t3_t2~adolescentis_receptive_score_t2+adolescentis_t2,df_adolescentis)
f.robftest(adolescentis_rlm,var="adolescentis_receptive_score_t2")
adolescentis_rlm[["coefficients"]][["adolescentis_receptive_score_t2"]]

#loocv (receptive score)
adolescentis_post_pred_breve_16_8 <- loocv_lm(df_adolescentis_filt,"adolescentis_increase_t3_t2","adolescentis_receptive_score_t2","adolescentis_t2")
corr.test(adolescentis_post_pred_breve_16_8[["loocv"]][["Predicted"]],adolescentis_post_pred_breve_16_8[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(adolescentis_post_pred_breve_16_8[["loocv"]][["Predicted"]],adolescentis_post_pred_breve_16_8[["loocv"]][["Actual"]],method = "spearman")$p

#------------------------------------------------------------------------------------
print("Trial: KrumbeckJ_2018 (adolescentis and detection)")

#Bifidobacterium_detection

t0 <- grep("_baseline$", rownames(KrumbeckJ_2018_SpProfile), value = TRUE)
t1 <- grep("_treatment$", rownames(KrumbeckJ_2018_SpProfile), value = TRUE)

common_species <- intersect(names(which(apply(KrumbeckJ_2018_SpProfile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_new_16s))

df_species_association_adolescentis_2 <- data.frame("detection_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],apply(KrumbeckJ_2018_SpProfile[,grep("Bifidobacterium",colnames(KrumbeckJ_2018_SpProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],apply(KrumbeckJ_2018_SpProfile[,grep("Bifidobacterium",colnames(KrumbeckJ_2018_SpProfile),value=TRUE)],1,function(x)(length(x[x>0]))),method="spearman")$r),0),
                                                    "longum_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_longum"],method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_longum"],method="spearman")$r),0),
                                                    "adolescentis_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_adolescentis"],method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_adolescentis"],method="spearman")$r),0),
                                                    "pseudocatenulatum_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_pseudocatenulatum"],method="spearman")$r),0),
                                                    "bifidum_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_bifidum"],method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_bifidum"],method="spearman")$r),0),
                                                    "dentium_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_dentium"],method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_dentium"],method="spearman")$r),0),
                                                    "breve_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_breve"],method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_breve"],method="spearman")$r),0),
                                                    "catenulatum_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_catenulatum"],method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_catenulatum"],method="spearman")$r),0),
                                                    "animalis_corr"=ifelse(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_animalis"],method="spearman")$p<=0.1,sign(corr.test(KrumbeckJ_2018_SpProfile[,common_species],KrumbeckJ_2018_SpProfile[,"Bifidobacterium_animalis"],method="spearman")$r),0))

df_mix_detection <- data.frame("detection_score"=as.matrix(KrumbeckJ_2018_SpProfile[t0,common_species])%*% ifelse(abs(rowMeans(df_association_all_new_16s[common_species,c("adult_detection","senior_detection")]))>0,sign(rowMeans(df_association_all_new_16s[common_species,c("adult_detection","senior_detection")])),0),
                               "detection_t0"=apply(KrumbeckJ_2018_SpProfile[t0,grep("Bifidobacterium",colnames(KrumbeckJ_2018_SpProfile),value=TRUE)],1,function(x)(length(x[x>0]))),
                               "detection_t1"=apply(KrumbeckJ_2018_SpProfile[t1,grep("Bifidobacterium",colnames(KrumbeckJ_2018_SpProfile),value=TRUE)],1,function(x)(length(x[x>0]))))

df_mix_detection$receptive_score <- rank_scale(df_mix_detection$detection_score) - rank_scale(df_mix_detection$detection_t0)

df_mix_detection$increase <- rank_scale(df_mix_detection$detection_t1) - rank_scale(df_mix_detection$detection_t0)

#rlm (recpetive score)
mix_rlm <- rlm(increase~receptive_score+detection_t0,df_mix_detection)
f.robftest(mix_rlm,var="receptive_score")
mix_rlm[["coefficients"]][["receptive_score"]]

#loocv (receptive score)
mix_1_post_pred_mix <- loocv_lm(df_mix_detection,"increase","receptive_score","detection_t0")
corr.test(mix_1_post_pred_mix[["loocv"]][["Predicted"]],mix_1_post_pred_mix[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(mix_1_post_pred_mix[["loocv"]][["Predicted"]],mix_1_post_pred_mix[["loocv"]][["Actual"]],method = "spearman")$p

#Bifidobacterium_adolescentis
metadata_adolescentis <- KrumbeckJ_2018_Metadata %>% filter(Treatment_Group %in% c("IVS-1+GOS", "IVS-1"))

t0 <- grep("_baseline$", rownames(metadata_adolescentis), value = TRUE)
t1 <- grep("_treatment$", rownames(metadata_adolescentis), value = TRUE)

common_species <- intersect(names(which(apply(KrumbeckJ_2018_SpProfile,2,function(x)(length(x[x>0])))>2)),rownames(df_association_all_new_16s))

df_adolescentis_2 <- data.frame("adolescentis_score"=as.matrix(KrumbeckJ_2018_SpProfile[t0,common_species])%*% ifelse(abs(rowMeans(df_association_all_new_16s[common_species,c("adult_adolescentis","senior_adolescentis")]))>0,sign(rowMeans(df_association_all_new_16s[common_species,c("adult_adolescentis","senior_adolescentis")])),0),
                                "adolescentis_t0"=KrumbeckJ_2018_SpProfile[t0, "Bifidobacterium_adolescentis"],
                                "adolescentis_t1"=KrumbeckJ_2018_SpProfile[t1, "Bifidobacterium_adolescentis"])

df_adolescentis_2$receptive_score <- rank_scale(df_adolescentis_2$adolescentis_score) - rank_scale(df_adolescentis_2$adolescentis_t0)

df_adolescentis_2$increase <- rank_scale(df_adolescentis_2$adolescentis_t1) - rank_scale(df_adolescentis_2$adolescentis_t0)

df_adolescentis_2_filt <- df_adolescentis_2[(df_adolescentis_2$adolescentis_t0>0)|(df_adolescentis_2$adolescentis_t1>0),]

#rlm (recpetive score)
adolescentis_2_rlm <- rlm(increase~receptive_score+adolescentis_t0,df_adolescentis_2)
f.robftest(adolescentis_2_rlm,var="receptive_score")
adolescentis_2_rlm[["coefficients"]][["receptive_score"]]

#loocv (receptive score)
adolescentis_2_post_pred <- loocv_lm(df_adolescentis_2_filt,"increase","receptive_score","adolescentis_t0")
corr.test(adolescentis_2_post_pred[["loocv"]][["Predicted"]],adolescentis_2_post_pred[["loocv"]][["Actual"]],method = "spearman")$r
corr.test(adolescentis_2_post_pred[["loocv"]][["Predicted"]],adolescentis_2_post_pred[["loocv"]][["Actual"]],method = "spearman")$p

#------------------------------------------------------------------------------------
##########################################################################################

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

animalis_ZhangQ_SunB[["spearman_r"]]
animalis_ZhangQ_SunB[["spearman_p"]]

animalis_SunB_ZhangQ <- train_test_lm(
  subset_df_animalis_2, subset_df_animalis_1,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_SunB_ZhangQ[["spearman_r"]]
animalis_SunB_ZhangQ[["spearman_p"]]

animalis_ZhangQ_BaZ <- train_test_lm(
  subset_df_animalis_1, subset_df_animalis_3_mean,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_ZhangQ_BaZ[["spearman_r"]]
animalis_ZhangQ_BaZ[["spearman_p"]]

animalis_BaZ_ZhangQ <- train_test_lm(
  subset_df_animalis_3_mean, subset_df_animalis_1,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_BaZ_ZhangQ[["spearman_r"]]
animalis_BaZ_ZhangQ[["spearman_p"]]

animalis_BaZ_SunB <- train_test_lm(
  subset_df_animalis_3_mean, subset_df_animalis_2,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_BaZ_SunB[["spearman_r"]]
animalis_BaZ_SunB[["spearman_p"]]

animalis_SunB_BaZ <- train_test_lm(
  subset_df_animalis_2, subset_df_animalis_3_mean,
  "animalis_increase", "animalis_receptive_score", "animalis_baseline")

animalis_SunB_BaZ[["spearman_r"]]
animalis_SunB_BaZ[["spearman_p"]]

# Scatter plot: actual vs predicted (similarly make all the six-combinations)

ggplot(animalis_ZhangQ_SunB[["loocv"]], aes(x = Actual, y = Predicted)) +
  geom_point(color = "cornflowerblue", size = 5) +
  geom_smooth(method = "lm", color = "black") + 
  labs(
    title = "Actual vs Predicted Increase (ZhangQ -> SunB)",
    x = "Actual Increase",
    y = "Predicted Increase"
  ) +
  theme_minimal(base_size = 12)
#------------------------------------------------------------------------------------

#Cross-cohort Validation for B. adolescentis

subset_df_adolescentis <- df_adolescentis_filt[, c("adolescentis_t2", "adolescentis_receptive_score_t2", "adolescentis_increase_t3_t2")]
colnames(subset_df_adolescentis) <- c("adolescentis_baseline", "adolescentis_receptive_score", "adolescentis_increase")

subset_df_adolescentis_2 <- df_adolescentis_2_filt[, c("adolescentis_t0", "receptive_score", "increase")]
colnames(subset_df_adolescentis_2) <- c("adolescentis_baseline", "adolescentis_receptive_score", "adolescentis_increase")

adolescentis_RamakrishnanM_KrumbeckJ <- train_test_lm(
  subset_df_adolescentis, subset_df_adolescentis_2,
  "adolescentis_increase", "adolescentis_receptive_score", "adolescentis_baseline")

adolescentis_RamakrishnanM_KrumbeckJ[["spearman_r"]]
adolescentis_RamakrishnanM_KrumbeckJ[["spearman_p"]]

adolescentis_KrumbeckJ_RamakrishnanM <- train_test_lm(
  subset_df_adolescentis_2, subset_df_adolescentis,
  "adolescentis_increase", "adolescentis_receptive_score", "adolescentis_baseline")

adolescentis_KrumbeckJ_RamakrishnanM[["spearman_r"]]
adolescentis_KrumbeckJ_RamakrishnanM[["spearman_p"]]
#------------------------------------------------------------------------------------

#making boxplot for the rest of the scenarios 
training_testing_df <- data.frame(
  training_cohort = c("SunB_2022","BaZ_2021","RamakrishnanM_2025","KrumbeckJ_2018"),
  testing_cohort = c("ZhangQ_2024","ZhangQ_2024","KrumbeckJ_2018","RamakrishnanM_2025"),
  coefficient = c(0.35,0.35,0.76,0.56),
  p_value = c(0.001,0.001,1.17e-07,0.0897))

training_testing_df

training_testing_df$significance <- ifelse(training_testing_df$p_value <= 0.05, "*", "")

## Preserve row order
training_testing_df$pair <- factor(
  paste(training_testing_df$training_cohort, "→", training_testing_df$testing_cohort),
  levels = paste(training_testing_df$training_cohort, "→", training_testing_df$testing_cohort)
)

p <- ggplot(training_testing_df, aes(x = pair, y = coefficient)) +
  geom_bar(
    stat = "identity",
    fill = "#80E5FF",
    width = 0.6   # <-- controls spacing between bars
  ) +
  geom_text(
    aes(label = ifelse(p_value <= 0.05, "*", "")),
    hjust = -0.2,
    size = 5
  ) +
  labs(
    x = "",
    y = "Spearman's rho",
    title = "Correlation between Training and Testing Cohorts"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  ) +
  ylim(0, max(training_testing_df$coefficient) + 0.1)

p

ggsave(filename = "/results/training_testing_correlation_barplot.pdf",plot = p, width = 7, height = 3, units = "in")
#------------------------------------------------------------------------------
##########################################################################################
#Association Correlations between Discovery and Validation Cohorts

##==== wgs adult/senior correlation between Association-Scores of discovery cohort and Association-Scores of validation cohort ====## 
load("/data/df_Strong_AssociationScores_new_wgs.RData")

############### Functions #################
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

# Function to compute positive, negative, and total
summarize_bif_association <- function(intervention_df) {
  positive <- rowSums(intervention_df == 1)
  negative <- rowSums(intervention_df == -1)
  total <- ncol(intervention_df)
  
  summary_df <- data.frame(positive = positive, negative = negative,total = total)
  
  return(summary_df)
}

#Function to combine adult/senior Association-Scores of Discovery Cohort
combine_mean_scores <- function(adult_df, senior_df) {
  scores_adult <- adult_df[, "score", drop = FALSE]
  scores_senior <- senior_df[, "score", drop = FALSE]
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

#Function for correlation between discovery and validation scores
compare_discovery_validation <- function(discovery_df, validation_df,
                                         discovery_col = "mean_score",
                                         validation_col = "scores",
                                         method = "spearman") {
  
  common_species <- intersect(rownames(discovery_df), rownames(validation_df))
  
  discovery_scores <- discovery_df[common_species, discovery_col]
  validation_scores <- validation_df[common_species, validation_col]
  
  corr_res <- corr.test(discovery_scores, validation_scores, method = method)
  
  result_df <- data.frame(
    correlation = corr_res$r,
    p_value = corr_res$p,
    significance = ifelse(corr_res$p <= 0.05, "*", "")
  )
}

df_list <- list(
  GomezM = df_species_association_longum_1,
  ZhangQ = df_species_association_animalis_1,
  SunB = df_species_association_animalis_2_new,
  GronbaekI = df_species_association_breve_1_new)

#create combined correlation data frame for any Bifidobacterium species
df_longum_association_intervention <- create_bif_corr_df(df_list, "longum")
df_adolescentis_association_intervention <- create_bif_corr_df(df_list, "adolescentis")
df_pseudocatenulatum_association_intervention <- create_bif_corr_df(df_list, "pseudocatenulatum")
df_bifidum_association_intervention <- create_bif_corr_df(df_list, "bifidum")
df_dentium_association_intervention <- create_bif_corr_df(df_list, "dentium")
df_breve_association_intervention <- create_bif_corr_df(df_list, "breve")
df_catenulatum_association_intervention <- create_bif_corr_df(df_list, "catenulatum")
df_animalis_association_intervention <- create_bif_corr_df(df_list, "animalis")

#Compute Association-Scores for validation cohort 
df_longum_association <- summarize_bif_association(df_longum_association_intervention)
df_adolescentis_association <- summarize_bif_association(df_adolescentis_association_intervention)
df_pseudocatenulatum_association <- summarize_bif_association(df_pseudocatenulatum_association_intervention)
df_bifidum_association <- summarize_bif_association(df_bifidum_association_intervention)
df_dentium_association <- summarize_bif_association(df_dentium_association_intervention)
df_breve_association <- summarize_bif_association(df_breve_association_intervention)
df_catenulatum_association <- summarize_bif_association(df_catenulatum_association_intervention)
df_animalis_association <- summarize_bif_association(df_animalis_association_intervention)

# add scores column
#((Positive-Negative)/(Total))*(1-(Min(Positive,Negative)+0.00001)/(Max(Positive,Negative)+0.00001))
df_longum_association$scores <- with(df_longum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_adolescentis_association$scores <- with(df_adolescentis_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_pseudocatenulatum_association$scores <- with(df_pseudocatenulatum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_bifidum_association$scores <- with(df_bifidum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_dentium_association$scores <- with(df_dentium_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_breve_association$scores <- with(df_breve_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_catenulatum_association$scores <- with(df_catenulatum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_animalis_association$scores <- with(df_animalis_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))

#combine adult/senior Association-Scores of Discovery Cohort
df_association_longum_wgs <- combine_mean_scores(df_association_adult_longum_wgs_new, df_association_senior_longum_wgs_new)
df_association_adolescentis_wgs <- combine_mean_scores(df_association_adult_adolescentis_wgs_new, df_association_senior_adolescentis_wgs_new)
df_association_pseudocatenulatum_wgs <- combine_mean_scores(df_association_adult_pseudocatenulatum_wgs_new, df_association_senior_pseudocatenulatum_wgs_new)
df_association_bifidum_wgs <- combine_mean_scores(df_association_adult_bifidum_wgs_new, df_association_senior_bifidum_wgs_new)
df_association_dentium_wgs <- combine_mean_scores(df_association_adult_dentium_wgs_new, df_association_senior_dentium_wgs_new)
df_association_breve_wgs <- combine_mean_scores(df_association_adult_breve_wgs_new, df_association_senior_breve_wgs_new)
df_association_catenulatum_wgs <- combine_mean_scores(df_association_adult_catenulatum_wgs_new, df_association_senior_catenulatum_wgs_new)
df_association_animalis_wgs <- combine_mean_scores(df_association_adult_animalis_wgs_new, df_association_senior_animalis_wgs_new)

#correlation 
species_pairs <- list(
  longum = list(df_association_longum_wgs, df_longum_association),
  adolescentis = list(df_association_adolescentis_wgs, df_adolescentis_association),
  pseudocatenulatum = list(df_association_pseudocatenulatum_wgs, df_pseudocatenulatum_association),
  bifidum = list(df_association_bifidum_wgs, df_bifidum_association),
  dentium = list(df_association_dentium_wgs, df_dentium_association),
  breve = list(df_association_breve_wgs, df_breve_association),
  catenulatum = list(df_association_catenulatum_wgs, df_catenulatum_association),
  animalis = list(df_association_animalis_wgs, df_animalis_association))

correlation_adult_senior_wgs <- imap_dfr(
  species_pairs,
  ~ compare_discovery_validation(.x[[1]], .x[[2]]) %>%
    mutate(species = .y)
) %>%
  dplyr::select(species, everything())

#Bar plot overall coefficient scores
correlation_adult_senior_wgs$species <- factor(correlation_adult_senior_wgs$species, levels = correlation_adult_senior_wgs$species)

p <- ggplot(correlation_adult_senior_wgs, aes(x = species, y = correlation)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(
    aes(label = significance, y = correlation + 0.02), size = 6) +
  labs(
    title = "Correlation between Discovery and Validation Cohorts",
    y = "Spearman's rho",
    x = "") +
  coord_flip() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12))

## Save plot
ggsave(filename = "/results/correlation_adult_senior_wgs.pdf", plot = p, width = 6, height = 2.5, units = "in")

#_________________________________________________________________________________________________________________________
#Table S18A
#merge the validation score into one df
bif_list <- list(
  longum       = df_longum_association,
  adolescentis = df_adolescentis_association,
  pseudocatenulatum = df_pseudocatenulatum_association,
  bifidum      = df_bifidum_association,
  dentium      = df_dentium_association,
  breve        = df_breve_association,
  catenulatum  = df_catenulatum_association,
  animalis     = df_animalis_association
)
all_species <- lapply(bif_list, rownames)
Reduce(intersect, all_species) -> common_species

# Check consistency
sapply(all_species, function(x) all(common_species %in% x))

score_list <- lapply(names(bif_list), function(name) {
  df <- bif_list[[name]]
  
  out <- df[common_species, "scores", drop = FALSE]
  colnames(out) <- paste0(name, "_scores")
  return(out)
})

names(score_list) <- names(bif_list)

bif_scores <- Reduce(function(x, y) cbind(x, y), score_list)

write.xlsx(bif_scores, file = "/results/Table_AssociationScore_Intervention_AdultSenior_wgs.xlsx", rowNames = TRUE)

#_________________________________________________________________________________________________________________________

##==== 16s adult/senior correlation between Association-Scores of discovery cohort and Association-Scores of validation cohort ====## 
load("/data/df_Strong_AssociationScores_new_16s.RData")

df_list <- list(
  BaZ = df_species_association_animalis_3_new,
  RamakrishnanM = df_species_association_adolescentis_1_new,
  KrumbeckJ = df_species_association_adolescentis_2)

#create combined correlation data frame for any Bifidobacterium species
df_longum_association_intervention <- create_bif_corr_df(df_list, "longum")
df_adolescentis_association_intervention <- create_bif_corr_df(df_list, "adolescentis")
df_pseudocatenulatum_association_intervention <- create_bif_corr_df(df_list, "pseudocatenulatum")
df_bifidum_association_intervention <- create_bif_corr_df(df_list, "bifidum")
df_dentium_association_intervention <- create_bif_corr_df(df_list, "dentium")
df_breve_association_intervention <- create_bif_corr_df(df_list, "breve")
df_catenulatum_association_intervention <- create_bif_corr_df(df_list, "catenulatum")
df_animalis_association_intervention <- create_bif_corr_df(df_list, "animalis")


#Compute Association-Scores for validation cohort 
df_longum_association <- summarize_bif_association(df_longum_association_intervention)
df_adolescentis_association <- summarize_bif_association(df_adolescentis_association_intervention)
df_pseudocatenulatum_association <- summarize_bif_association(df_pseudocatenulatum_association_intervention)
df_bifidum_association <- summarize_bif_association(df_bifidum_association_intervention)
df_dentium_association <- summarize_bif_association(df_dentium_association_intervention)
df_breve_association <- summarize_bif_association(df_breve_association_intervention)
df_catenulatum_association <- summarize_bif_association(df_catenulatum_association_intervention)
df_animalis_association <- summarize_bif_association(df_animalis_association_intervention)

# add scores column
#((Positive-Negative)/(Total))*(1-(Min(Positive,Negative)+0.00001)/(Max(Positive,Negative)+0.00001))
df_longum_association$scores <- with(df_longum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_adolescentis_association$scores <- with(df_adolescentis_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_pseudocatenulatum_association$scores <- with(df_pseudocatenulatum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_bifidum_association$scores <- with(df_bifidum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_dentium_association$scores <- with(df_dentium_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_breve_association$scores <- with(df_breve_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_catenulatum_association$scores <- with(df_catenulatum_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))
df_animalis_association$scores <- with(df_animalis_association, ((positive - negative) / total) * (1 - (pmin(positive, negative) + 0.00001) / (pmax(positive, negative) + 0.00001)))

#combine adult/senior Association-Scores of Discovery Cohort
df_association_longum_16s <- combine_mean_scores(df_association_adult_longum_16s_new, df_association_senior_longum_16s_new)
df_association_adolescentis_16s <- combine_mean_scores(df_association_adult_adolescentis_16s_new, df_association_senior_adolescentis_16s_new)
df_association_pseudocatenulatum_16s <- combine_mean_scores(df_association_adult_pseudocatenulatum_16s_new, df_association_senior_pseudocatenulatum_16s_new)
df_association_bifidum_16s <- combine_mean_scores(df_association_adult_bifidum_16s_new, df_association_senior_bifidum_16s_new)
df_association_dentium_16s <- combine_mean_scores(df_association_adult_dentium_16s_new, df_association_senior_dentium_16s_new)
df_association_breve_16s <- combine_mean_scores(df_association_adult_breve_16s_new, df_association_senior_breve_16s_new)
df_association_catenulatum_16s <- combine_mean_scores(df_association_adult_catenulatum_16s_new, df_association_senior_catenulatum_16s_new)
df_association_animalis_16s <- combine_mean_scores(df_association_adult_animalis_16s_new, df_association_senior_animalis_16s_new)

#correlation 
species_pairs <- list(
  longum = list(df_association_longum_16s, df_longum_association),
  adolescentis = list(df_association_adolescentis_16s, df_adolescentis_association),
  pseudocatenulatum = list(df_association_pseudocatenulatum_16s, df_pseudocatenulatum_association),
  bifidum = list(df_association_bifidum_16s, df_bifidum_association),
  dentium = list(df_association_dentium_16s, df_dentium_association),
  breve = list(df_association_breve_16s, df_breve_association),
  catenulatum = list(df_association_catenulatum_16s, df_catenulatum_association),
  animalis = list(df_association_animalis_16s, df_animalis_association))

correlation_adult_senior_16s <- imap_dfr(
  species_pairs,
  ~ compare_discovery_validation(.x[[1]], .x[[2]]) %>%
    mutate(species = .y)
) %>%
  dplyr::select(species, everything())

#Bar plot overall coefficient scores
correlation_adult_senior_16s$species <- factor(correlation_adult_senior_16s$species, levels = correlation_adult_senior_16s$species)

p <- ggplot(correlation_adult_senior_16s, aes(x = species, y = correlation)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(
    aes(label = significance, y = correlation + 0.02), size = 6) +
  labs(
    title = "Correlation between Discovery and Validation Cohorts",
    y = "Spearman's rho",
    x = "") +
  #coord_flip() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12))

## Save plot
ggsave(filename = "/results/correlation_adult_senior_16s_horizontal.pdf", plot = p, width = 6, height = 2.5, units = "in")

#_________________________________________________________________________________________________________________________
#Table S18B
#merge the validation score into one df
bif_list <- list(
  longum       = df_longum_association,
  adolescentis = df_adolescentis_association,
  pseudocatenulatum = df_pseudocatenulatum_association,
  bifidum      = df_bifidum_association,
  dentium      = df_dentium_association,
  breve        = df_breve_association,
  catenulatum  = df_catenulatum_association,
  animalis     = df_animalis_association
)
all_species <- lapply(bif_list, rownames)
Reduce(intersect, all_species) -> common_species

# Check consistency
sapply(all_species, function(x) all(common_species %in% x))

score_list <- lapply(names(bif_list), function(name) {
  df <- bif_list[[name]]
  
  out <- df[common_species, "scores", drop = FALSE]
  colnames(out) <- paste0(name, "_scores")
  return(out)
})

names(score_list) <- names(bif_list)

bif_scores <- Reduce(function(x, y) cbind(x, y), score_list)

write.xlsx(bif_scores, file = "/results/Table_AssociationScore_Intervention_AdultSenior_16s.xlsx", rowNames = TRUE)

#_________________________________________________________________________________________________________________________
##==== infant monica correlation between association score and bif correlation ====## 

load("/data/df_Strong_AssociationScores_new_16s.RData")

###### Correlation Function ######
#Function to run correlation 
run_species_correlation <- function(x_df, y_df, x_col, y_col = "score",
                                    method = "spearman", exact = TRUE) {
  
  common_species <- intersect(rownames(x_df), rownames(y_df))
  
  x_vals <- x_df[common_species, x_col]
  y_vals <- y_df[common_species, y_col]
  
  cor_res <- cor.test(x_vals, y_vals, method = method, exact = exact)
  
  data.frame(
    correlation = unname(cor_res$estimate),
    p_value = cor_res$p.value,
    significance = ifelse(cor_res$p.value <= 0.05, "*", "")
  )
}

species_map <- list(
  longum = list("longum_corr", df_association_infant_longum_16s_new),
  adolescentis = list("adolescentis_corr", df_association_infant_adolescentis_16s_new),
  pseudocatenulatum = list("pseudocatenulatum_corr", df_association_infant_pseudocatenulatum_16s_new),
  bifidum = list("bifidum_corr", df_association_infant_bifidum_16s_new),
  dentium = list("dentium_corr", df_association_infant_dentium_16s_new),
  breve = list("breve_corr", df_association_infant_breve_16s_new),
  catenulatum = list("catenulatum_corr", df_association_infant_catenulatum_16s_new)
)


#Table S19
#View(df_species_association_mix_1)
write.xlsx(df_species_association_mix_1, file = "/results/Table_Correlation_Intervention_Infant_16s.xlsx", rowNames = TRUE)

correlation_BazanellaM_infant_16s <- imap_dfr(
  species_map,
  ~ run_species_correlation(
    x_df = df_species_association_mix_1,
    y_df = .x[[2]],
    x_col = .x[[1]]
  ) %>%
    mutate(species = .y)
) %>%
  dplyr::select(species, everything())

#Bar plot overall coefficient scores
correlation_BazanellaM_infant_16s$species <- factor(correlation_BazanellaM_infant_16s$species, levels = correlation_BazanellaM_infant_16s$species)

p <- ggplot(correlation_BazanellaM_infant_16s, aes(x = species, y = correlation)) +
  geom_bar(stat = "identity", fill = "burlywood") +
  geom_text(
    aes(label = significance, y = correlation + 0.02), size = 6) +
  labs(
    title = "Coefficient Scores per Bifidobacterium Species",
    y = "Spearman's rho",
    x = "") +
  coord_flip() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12))

## Save plot
ggsave(filename = "/results/correlation_BazanellaM_infant_16s.pdf", plot = p, width = 6, height = 2.5, units = "in")

#------------------------------------------------------------------------------------


