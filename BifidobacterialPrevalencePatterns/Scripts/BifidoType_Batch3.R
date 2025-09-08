
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
# Finding studies for each species where the detection rate in the study is greater than or equal to 10 percent
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


REM_AS_Bifidobacterium_longum <- compute_meta_corr_group(refined_data,refined_species_list,"Bifidobacterium_longum","study_name",Bifidobacterium_longum_AS_studies)
REM_AS_Bifidobacterium_adolescentis <- compute_meta_corr_group(refined_data,refined_species_list,"Bifidobacterium_adolescentis","study_name",Bifidobacterium_adolescentis_AS_studies)
REM_AS_Bifidobacterium_pseudocatenulatum <- compute_meta_corr_group(refined_data,refined_species_list,"Bifidobacterium_pseudocatenulatum","study_name",Bifidobacterium_longum_AS_studies)
REM_AS_Bifidobacterium_breve <- compute_meta_corr_group(refined_data,refined_species_list,"Bifidobacterium_breve","study_name",Bifidobacterium_breve_AS_studies)
REM_AS_Bifidobacterium_dentium <- compute_meta_corr_group(refined_data,refined_species_list,"Bifidobacterium_dentium","study_name",Bifidobacterium_longum_AS_studies)
REM_AS_Bifidobacterium_bifidum <- compute_meta_corr_group(refined_data,refined_species_list,"Bifidobacterium_bifidum","study_name",Bifidobacterium_longum_AS_studies)
REM_AS_Bifidobacterium_catenulatum <- compute_meta_corr_group(refined_data,refined_species_list,"Bifidobacterium_catenulatum","study_name",Bifidobacterium_longum_AS_studies)
REM_AS_Bifidobacterium_animalis <- compute_meta_corr_group(refined_data,refined_species_list,"Bifidobacterium_animalis","study_name",Bifidobacterium_longum_AS_studies)
REM_AS_TotalBifsDetection <- compute_meta_corr_group(refined_data,refined_species_list,"TotalBifsDetection","study_name",unique(adult_select_studies,senior_select_studies))

Direction_AS_TotalBifsDetection <- data.frame("Positive"=apply(RLM_AdultSenior_TotalBifsDetection$dir,2,function(x)(length(x[x>= 2]))),"Negative"=apply(RLM_AdultSenior_TotalBifsDetection$dir,2,function(x)(length(x[x<= -2]))))

Direction_TotalBifsDetection_Filtered <- Direction_TotalBifsDetection[((Direction_TotalBifsDetection[,1]>=10)&(Direction_TotalBifsDetection[,2]<=5))|((Direction_TotalBifsDetection[,2]>=10)&(Direction_TotalBifsDetection[,1]<=5)),]

REM_Filtered_AS_TotalBifsDetection <- REM_AS_TotalBifsDetection[(abs(REM_AS_TotalBifsDetection$dir)>=3)&(REM_AS_TotalBifsDetection$consistency>0.5),]

select_species <- intersect(rownames(REM_Filtered_AS_TotalBifsDetection),rownames(Direction_TotalBifsDetection_Filtered))

REM_Filtered_AS_TotalBifsDetection <- REM_Filtered_AS_TotalBifsDetection[select_species,]

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

save.image("G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\Bifidotypes_Batch3.RData")

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