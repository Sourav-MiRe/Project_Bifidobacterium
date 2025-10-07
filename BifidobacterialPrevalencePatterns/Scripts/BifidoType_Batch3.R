
load("/data/Bifidotypes_Batch2.RData")
library(dplyr)
library(pcaPP)
library(vegan)
library(gplots)
library(RColorBrewer)

compute_detection <- function(data,var1_list,grouping_variable,grouping_list)
{
  detection_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
  rownames(detection_matrix) <- var1_list
  colnames(detection_matrix) <- grouping_list
  for(i in 1:length(var1_list))
  {
    var1 <- var1_list[i]
    for(j in 1:length(grouping_list))
    {
      group <- grouping_list[j]
      detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]>0))/length(data[data[,grouping_variable]==group,var1])
    }
  }
  return(detection_matrix)
}

batch_corr_grouped <- function(data,metadata,variable_group,metadata_feature,grouping_feature,grouping_list)
{
  df_est <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
  rownames(df_est) <- grouping_list
  colnames(df_est) <- variable_group
  df_p_val <- as.data.frame(matrix(1,length(grouping_list),length(variable_group)))
  rownames(df_p_val) <- grouping_list
  colnames(df_p_val) <- variable_group
  print(grouping_list)
  for(i in 1:length(grouping_list))
  {
    study_name <- grouping_list[i]
    print(study_name)
    study_samples <- rownames(metadata[metadata$study_name == study_name,])
    for(j in 1:length(variable_group))
    {
      species_name <- variable_group[j]
      species_val <- data[study_samples,species_name]
      metadata_val <- metadata[study_samples,metadata_feature]
      species_detect <- length(species_val[species_val>0])/length(species_val)
      metadata_detect <- length(metadata_val[metadata_val>0])/length(metadata_val)
      print(length(species_val[species_val>0]))
      tryCatch(               
        expr = {                     
          if((species_detect > 0.20)&(metadata_detect > 0.20))
          {
            print(paste0(study_name,":",species_name))
            temp_corr <- cor.test(species_val,metadata_val,method="kendall") 
            df_est[i,j] <- as.numeric(temp_corr$estimate)
            df_p_val[i,j] <- as.numeric(temp_corr$p.value)
          }
        },
        error = function(e){ 
          print(e)
          print("Error observed. Moving to next")
        },
        finally = {            
          print("finally Executed")
        }
      )
    }
    
  }
  df_q_val <- apply(df_p_val,2,p.adjust)
  l_fisher <- p.adjust(apply(df_q_val,2,function(x)(sumlog(x)$p)),method="fdr")
  print("l_fisher generated")
  df_dir <- as.data.frame(matrix(0,length(grouping_list),length(variable_group)))
  rownames(df_dir) <- grouping_list
  colnames(df_dir) <- variable_group
  for(i in 1:length(grouping_list))
  {
    for(j in 1:length(variable_group))
    {
      df_dir[i,j] <- ifelse(df_q_val[i,j]<=0.10,3*sign(df_est[i,j]),ifelse(df_p_val[i,j]<=0.05,2*sign(df_est[i,j]),1*sign(df_est[i,j])))
    }
  }
  return_list <- list("est"=df_est,"p.value"=df_p_val,"q.value"=df_q_val,"fisher"=l_fisher,"dir"=df_dir)
  return(return_list)
}

oob_validation <- function(data,metadata,features_to_check,grouping_column,metadata_column,groups_list)
{
  df_Perform <- as.data.frame(matrix(NA,length(groups_list),2))
  rownames(df_Perform) <- groups_list
  colnames(df_Perform) <- c("oob_corr","oob_p")
  
  pairwise_matrix_r <- as.data.frame(matrix(NA,length(groups_list),length(groups_list)))
  rownames(pairwise_matrix_r) <- groups_list
  colnames(pairwise_matrix_r) <- groups_list
  
  pairwise_matrix_p <- as.data.frame(matrix(NA,length(groups_list),length(groups_list)))
  rownames(pairwise_matrix_p) <- groups_list
  colnames(pairwise_matrix_p) <- groups_list
  
  mat_importance <- as.data.frame(matrix(0,length(groups_list),length(features_to_check)))
  rownames(mat_importance) <- groups_list
  colnames(mat_importance) <- features_to_check
  
  group_list_rows <- rownames(data[data[,grouping_column] %in% groups_list,])
  rows_with_metadata <- rownames(metadata[!is.na(metadata[,metadata_column]),])
  df_pred <- data[intersect(rows_with_metadata,group_list_rows),features_to_check]
  
  for(i in 1:length(groups_list))
  {
    print(i)
    group_name <- groups_list[i]
    train_rows <- intersect(rownames(df_pred),rownames(metadata)[metadata[,grouping_column] == group_name])
    print(train_rows)
    df_train_pred <- df_pred[train_rows,]
    df_train_pred <- apply(df_train_pred,2,function(x)(ifelse(is.nan(x),0,x)))
    df_train_pred <- apply(df_train_pred,2,function(x)(ifelse(is.na(x),0,x)))
    df_train_pred <- df_train_pred[rowSums(df_train_pred)>0,colSums(df_train_pred)>0]
    train_rows <- rownames(df_train_pred)
    
    df_train_response <- metadata[train_rows,metadata_column]
    
    rf_oob_temp <- randomForest(df_train_response~.,df_train_pred)
    
    df_oob_corr <- cor.test(rf_oob_temp$predicted,rf_oob_temp$y,method="spearman",use="pairwise.complete")
    
    print(df_oob_corr)
    
    mat_importance[group_name,rownames(rf_oob_temp$importance)] <- rf_oob_temp$importance
    
    for(j in 1:length(groups_list))
    {
      test_study <- groups_list[j]
      if(i == j)
      {
        pairwise_matrix_r[i,j] <- df_oob_corr$estimate
        pairwise_matrix_p[i,j] <- df_oob_corr$p.value
      }
      else
      {
        test_study_rows <- intersect(rownames(df_pred),rownames(metadata)[metadata[,grouping_column] == test_study])
        test_study_rows <- intersect(rows_with_metadata,test_study_rows)
        df_pred_test_study <- data[test_study_rows,colnames(df_train_pred)]
        response_test_study <- metadata[test_study_rows,metadata_column]
        df_test_study_corr <- cor.test(predict(rf_oob_temp,df_pred_test_study),response_test_study,method="spearman",use="pairwise.complete")
        pairwise_matrix_r[i,j] <- df_test_study_corr$estimate
        pairwise_matrix_p[i,j] <- df_test_study_corr$p.value
      }
    }
    
    df_Perform[i,1] <- df_oob_corr$estimate
    df_Perform[i,2] <- df_oob_corr$p.value
    
    i <- i + 1
  }
  
  df_Correlation <- list("r"=pairwise_matrix_r,"p"=pairwise_matrix_p)
  
  returnlist <- list("Perform"=df_Perform,"Correlation"=df_Correlation,"FeatureImportance"=mat_importance)
  
  return(returnlist)
}

TotalBifsDetection <- rowSums(spdata_Bifs_detect_45809)

spdata_45809[is.na(spdata_45809)] <- 0

# Finding the different bifidobacterium detected, considering all the bifidobacterium
spdata_45809$TotalBifsDetection <- TotalBifsDetection[rownames(spdata_Bifs_detect_45809)]
# Taking union of non-bif selected species in adult and in senior
refined_species_list <- union(adult_select_non_bif_species,senior_select_non_bif_species)

refined_species_list <- grep("CAG",refined_species_list,value=TRUE,invert=TRUE)
# Making another dataframe where the rows are of both adult and senior samples and columns are of selected species from both age cateogory
refined_data <- spdata_45809[c(adult,senior),refined_species_list]
# Renormalize the data
refined_data <- as.data.frame(t(apply(refined_data,1,function(x)(x/sum(x)))))

refined_data$study_name <- spdata_45809[rownames(refined_data),"study_name"]
refined_data$TotalBifsDetection <- spdata_45809[rownames(refined_data),"TotalBifsDetection"]
# Adding data to refined data
refined_data$Bifidobacterium_adolescentis <- spdata_Bifs_rel_45809[rownames(refined_data),"Bifidobacterium_adolescentis"]
refined_data$Bifidobacterium_longum <- spdata_Bifs_rel_45809[rownames(refined_data),"Bifidobacterium_longum"]
refined_data$Bifidobacterium_pseudocatenulatum <- spdata_Bifs_rel_45809[rownames(refined_data),"Bifidobacterium_pseudocatenulatum"]
refined_data$Bifidobacterium_bifidum <- spdata_Bifs_rel_45809[rownames(refined_data),"Bifidobacterium_bifidum"]
refined_data$Bifidobacterium_dentium <- spdata_Bifs_rel_45809[rownames(refined_data),"Bifidobacterium_dentium"]
refined_data$Bifidobacterium_breve <- spdata_Bifs_rel_45809[rownames(refined_data),"Bifidobacterium_breve"]
refined_data$Bifidobacterium_catenulatum <- spdata_Bifs_rel_45809[rownames(refined_data),"Bifidobacterium_catenulatum"]
refined_data$Bifidobacterium_animalis <- spdata_Bifs_rel_45809[rownames(refined_data),"Bifidobacterium_animalis"]

MajorBifsDetectionAdultSenior <- compute_detection(refined_data,MajorBifsOrdered,"study_name",unique(adult_select_studies,senior_select_studies))
Bifidobacterium_animalis_AS_studies <- colnames(MajorBifsDetectionAdultSenior)[which(MajorBifsDetectionAdultSenior["Bifidobacterium_animalis",]>=0.10)]
Bifidobacterium_catenulatum_AS_studies <- colnames(MajorBifsDetectionAdultSenior)[which(MajorBifsDetectionAdultSenior["Bifidobacterium_catenulatum",]>=0.10)]
Bifidobacterium_breve_AS_studies <- colnames(MajorBifsDetectionAdultSenior)[which(MajorBifsDetectionAdultSenior["Bifidobacterium_breve",]>=0.10)]
Bifidobacterium_dentium_AS_studies <- colnames(MajorBifsDetectionAdultSenior)[which(MajorBifsDetectionAdultSenior["Bifidobacterium_dentium",]>=0.10)]
Bifidobacterium_bifidum_AS_studies <- colnames(MajorBifsDetectionAdultSenior)[which(MajorBifsDetectionAdultSenior["Bifidobacterium_bifidum",]>=0.10)]
Bifidobacterium_pseudocatenulatum_AS_studies <- colnames(MajorBifsDetectionAdultSenior)[which(MajorBifsDetectionAdultSenior["Bifidobacterium_pseudocatenulatum",]>=0.10)]
Bifidobacterium_adolescentis_AS_studies <- colnames(MajorBifsDetectionAdultSenior)[which(MajorBifsDetectionAdultSenior["Bifidobacterium_adolescentis",]>=0.10)]
Bifidobacterium_longum_AS_studies <- colnames(MajorBifsDetectionAdultSenior)[which(MajorBifsDetectionAdultSenior["Bifidobacterium_longum",]>=0.10)]

RLM_AdultSenior_Bifidobacterium_longum <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"Bifidobacterium_longum","study_name",Bifidobacterium_longum_AS_studies)
RLM_AdultSenior_Bifidobacterium_adolescentis <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"Bifidobacterium_adolescentis","study_name",Bifidobacterium_adolescentis_AS_studies)
RLM_AdultSenior_Bifidobacterium_pseudocatenulatum <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"Bifidobacterium_pseudocatenulatum","study_name",Bifidobacterium_pseudocatenulatum_AS_studies)
RLM_AdultSenior_Bifidobacterium_bifidum <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"Bifidobacterium_bifidum","study_name",Bifidobacterium_bifidum_AS_studies)
RLM_AdultSenior_Bifidobacterium_dentium <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"Bifidobacterium_dentium","study_name",Bifidobacterium_dentium_AS_studies)
RLM_AdultSenior_Bifidobacterium_breve <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"Bifidobacterium_breve","study_name",Bifidobacterium_breve_AS_studies)
RLM_AdultSenior_Bifidobacterium_catenulatum <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"Bifidobacterium_catenulatum","study_name",Bifidobacterium_catenulatum_AS_studies)
RLM_AdultSenior_Bifidobacterium_animalis <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"Bifidobacterium_animalis","study_name",Bifidobacterium_animalis_AS_studies)
RLM_AdultSenior_TotalBifsDetection <- batch_corr_grouped(refined_data,refined_data,refined_species_list,"TotalBifsDetection","study_name",unique(adult_select_studies,senior_select_studies))


NonBifsAssociationAdultSenior <- 
data.frame("Bifidobacterium_longum"=(apply(RLM_AdultSenior_Bifidobacterium_longum$dir,2,function(x)(length(x[!is.na(x)&(x>=2)]))) - apply(RLM_AdultSenior_Bifidobacterium_longum$dir,2,function(x)(length(x[!is.na(x)&(x<=-2)]))))/length(Bifidobacterium_longum_AS_studies),"Bifidobacterium_adolescentis"=(apply(RLM_AdultSenior_Bifidobacterium_adolescentis$dir,2,function(x)(length(x[!is.na(x)&(x>=2)]))) - apply(RLM_AdultSenior_Bifidobacterium_adolescentis$dir,2,function(x)(length(x[!is.na(x)&(x<=-2)]))))/length(Bifidobacterium_adolescentis_AS_studies),"Bifidobacterium_pseudocatenulatum"=(apply(RLM_AdultSenior_Bifidobacterium_pseudocatenulatum$dir,2,function(x)(length(x[!is.na(x)&(x>=2)]))) - apply(RLM_AdultSenior_Bifidobacterium_pseudocatenulatum$dir,2,function(x)(length(x[!is.na(x)&(x<=-2)]))))/length(Bifidobacterium_pseudocatenulatum_AS_studies),"Bifidobacterium_bifidum"=(apply(RLM_AdultSenior_Bifidobacterium_bifidum$dir,2,function(x)(length(x[!is.na(x)&(x>=2)]))) - apply(RLM_AdultSenior_Bifidobacterium_bifidum$dir,2,function(x)(length(x[!is.na(x)&(x<=-2)]))))/length(Bifidobacterium_bifidum_AS_studies),"Bifidobacterium_breve"=(apply(RLM_AdultSenior_Bifidobacterium_breve$dir,2,function(x)(length(x[!is.na(x)&(x>=2)]))) - apply(RLM_AdultSenior_Bifidobacterium_breve$dir,2,function(x)(length(x[!is.na(x)&(x<=-2)]))))/length(Bifidobacterium_breve_AS_studies),"Bifidobacterium_dentium"=(apply(RLM_AdultSenior_Bifidobacterium_dentium$dir,2,function(x)(length(x[!is.na(x)&(x>=2)]))) - apply(RLM_AdultSenior_Bifidobacterium_dentium$dir,2,function(x)(length(x[!is.na(x)&(x<=-2)]))))/length(Bifidobacterium_dentium_AS_studies),"Bifidobacterium_catenulatum"=(apply(RLM_AdultSenior_Bifidobacterium_catenulatum$dir,2,function(x)(length(x[!is.na(x)&(x>=2)]))) - apply(RLM_AdultSenior_Bifidobacterium_catenulatum$dir,2,function(x)(length(x[!is.na(x)&(x<=-2)]))))/length(Bifidobacterium_catenulatum_AS_studies),"Bifidobacterium_animalis"=(apply(RLM_AdultSenior_Bifidobacterium_animalis$dir,2,function(x)(length(x[!is.na(x)&(x>=2)]))) - apply(RLM_AdultSenior_Bifidobacterium_animalis$dir,2,function(x)(length(x[!is.na(x)&(x<=-2)]))))/length(Bifidobacterium_animalis_AS_studies))


Direction_TotalBifsDetection <- data.frame("Positive"=apply(RLM_AdultSenior_TotalBifsDetection$dir,2,function(x)(length(x[x>= 2]))),"Negative"=apply(RLM_AdultSenior_TotalBifsDetection$dir,2,function(x)(length(x[x<= -2]))))

Direction_TotalBifsDetection_Filtered <- Direction_TotalBifsDetection[((Direction_TotalBifsDetection[,1]>=10)&(Direction_TotalBifsDetection[,2]<=5))|((Direction_TotalBifsDetection[,2]>=10)&(Direction_TotalBifsDetection[,1]<=5)),]

select_species <- intersect(names(RLM_AdultSenior_TotalBifsDetection$dir),rownames(Direction_TotalBifsDetection_Filtered))

df_combined_bifido$CohortType <- factor(df_combined_bifido$CohortType,levels=c("IndustrializedUrban","UrbanRuralMixed","RuralTribal"))

#df_combined_bifido$AgeCategory <- factor(df_combined_bifido$CohortType,levels=c("Infant","Adult","S"))

LM_MajorBifsDetection_AgeCategory <- list("r"=matrix(NA,9,7),"p"=matrix(NA,9,7))
rownames(LM_MajorBifsDetection_AgeCategory$r) <- colnames(df_combined_bifido)[1:9]
colnames(LM_MajorBifsDetection_AgeCategory$r) <- c("Infant","Senior","Senior/I","WGS","UrbanRuralMixed","RuralTribal","RuralTribal/M")

rownames(LM_MajorBifsDetection_AgeCategory$p) <- colnames(df_combined_bifido)[1:9]
colnames(LM_MajorBifsDetection_AgeCategory$p) <- c("Infant","Senior","Senior/I","WGS","UrbanRuralMixed","RuralTribal","RuralTribal/M")

for(i in 1:9)
{
	species <- colnames(df_combined_bifido)[i]
	summary_lm <- summary(lm(as.formula(paste0(species,"~","AgeCategory+SequencingType+CohortType")),data=df_combined_bifido))
	summary_lm_1 <- summary(lm(as.formula(paste0(species,"~","AgeCategory+SequencingType")),data=df_combined_bifido[df_combined_bifido$AgeCategory!="Adult",]))
	summary_lm_2 <- summary(lm(as.formula(paste0(species,"~","SequencingType+CohortType")),data=df_combined_bifido[(df_combined_bifido$AgeCategory=="Adult")&(df_combined_bifido$CohortType!="IndustrializedUrban"),]))
	
	LM_MajorBifsDetection_AgeCategory$r[species,1] <- summary_lm$coefficients[2,1]
	LM_MajorBifsDetection_AgeCategory$r[species,2] <- summary_lm$coefficients[3,1]
	LM_MajorBifsDetection_AgeCategory$r[species,3] <- summary_lm_1$coefficients[2,1]
	LM_MajorBifsDetection_AgeCategory$r[species,4] <- summary_lm$coefficients[4,1]
	LM_MajorBifsDetection_AgeCategory$r[species,5] <- summary_lm$coefficients[5,1]
	LM_MajorBifsDetection_AgeCategory$r[species,6] <- summary_lm$coefficients[6,1]
	LM_MajorBifsDetection_AgeCategory$r[species,7] <- summary_lm_2$coefficients[3,1]
	
	LM_MajorBifsDetection_AgeCategory$p[species,1] <- summary_lm$coefficients[2,4]
	LM_MajorBifsDetection_AgeCategory$p[species,2] <- summary_lm$coefficients[3,4]
	LM_MajorBifsDetection_AgeCategory$p[species,3] <- summary_lm_1$coefficients[2,4]
	LM_MajorBifsDetection_AgeCategory$p[species,4] <- summary_lm$coefficients[4,4]
	LM_MajorBifsDetection_AgeCategory$p[species,5] <- summary_lm$coefficients[5,4]
	LM_MajorBifsDetection_AgeCategory$p[species,6] <- summary_lm$coefficients[6,4]
	LM_MajorBifsDetection_AgeCategory$p[species,7] <- summary_lm_2$coefficients[3,4]
}
	
temp_mat <- t(ifelse(LM_MajorBifsDetection_AgeCategory$p<=0.10,sign(LM_MajorBifsDetection_AgeCategory$r),0))

temp_p <- apply(LM_MajorBifsDetection_AgeCategory$p,2,function(x)(ifelse(x<0.001,"***",ifelse(x<0.01,"**",ifelse(x<0.05,"*",ifelse(x<0.1,"@",""))))))

heatmap.2(temp_mat,density="none",trace="none",Rowv=FALSE,Colv=FALSE,col=c("hotpink","white","cornflowerblue"),margins=c(15,10),sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)))

colorInfant <- rgb(169,209,142,max=255)
colorAdult <- rgb(157,195,230,max=255)
colorSenior <- rgb(255,217,102,max=255)

s.class(pcoStudyPrevEuclidean$li,as.factor(df_combined_bifido$AgeCategory),col=c(colorInfant,colorAdult,colorSenior),plabels.col="black",plabels.cex=1.5)

colorWGS <- rgb(100,149,237,max=255)
color16S <- rgb(184,134,11,max=255)

s.class(pcoStudyPrevEuclidean$li,factor(df_combined_bifido$SequencingType,levels=c("wgs","16s")),col=c(colorWGS,color16S),plabels.col="black",plabels.cex=1.5)

save.image("../results/Bifidotypes_Batch3.RData")

RLM_AdultSenior_Bifidobacterium_longum$dir_1 <- sign(RLM_AdultSenior_Bifidobacterium_longum$est) * apply(RLM_AdultSenior_Bifidobacterium_longum$p.value,2,function(x)(ifelse(x<=0.1,1,0)))

RLM_AdultSenior_Bifidobacterium_longum$p.notation <- apply(RLM_AdultSenior_Bifidobacterium_longum$p.value,2,function(x)(ifelse(x<=0.05,"*",ifelse(x<=0.1,"^",""))))

RF_Bifidobacterium_longum <- oob_validation(refined_data,refined_data,select_species,"study_name","Bifidobacterium_longum",Bifidobacterium_longum_AS_studies)

hmp_RF_Bifidobacterium_longum <- heatmap.2(t(apply(RF_Bifidobacterium_longum$FeatureImportance,1,rank_scale)),density="none",trace="none")

temp_mat <- t(apply(RF_Bifidobacterium_longum$FeatureImportance,1,rank_scale))

temp_noteCol <- apply(RLM_AdultSenior_Bifidobacterium_longum$dir_1[rownames(temp_mat),colnames(temp_mat)],2,function(x)(ifelse(x== -1,"red4",ifelse(x==1,"blue4","white"))))

temp_cellNote <- RLM_AdultSenior_Bifidobacterium_longum$p.notation[rownames(temp_mat),colnames(temp_mat)]

temp_mat <- temp_mat[colnames(hmp_RF_Bifidobacterium_longum$carpet),rownames(hmp_RF_Bifidobacterium_longum$carpet)]

temp_noteCol <- temp_noteCol[colnames(hmp_RF_Bifidobacterium_longum$carpet),rownames(hmp_RF_Bifidobacterium_longum$carpet)]

temp_cellNote <- temp_cellNote[colnames(hmp_RF_Bifidobacterium_longum$carpet),rownames(hmp_RF_Bifidobacterium_longum$carpet)]

heatmap.2(temp_mat,density="none",trace="none",Rowv=FALSE,Colv=FALSE,lhei=c(0.5,5),lwid=c(0.5,5),srtCol=90,srtRow=0,cexRow=0.5,cexCol=0.68,margins=c(10,8),cellnote=temp_cellNote,notecol=temp_noteCol,notecex=0.5)