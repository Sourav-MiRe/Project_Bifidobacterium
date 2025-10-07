
load("../Data/Bifidotypes_Batch1.RData")
library(dplyr)
library(ade4)
library(vegan)
library(adegraphics)
library(gplots)
library(RColorBrewer)

MajorBifsOrdered <- colnames(df_infant_bifido)[order(apply(rbind(df_infant_bifido[,1:8],df_adult_bifido[,1:8],df_senior_bifido[,1:8]),2,function(x)(length(x[x>=0.3]))))]

df_infant_bifido$AgeCategory <- NA
df_infant_bifido[,"AgeCategory"] <- "Infant"

df_adult_bifido$AgeCategory <- NA
df_adult_bifido[,"AgeCategory"] <- "Adult"

df_senior_bifido$AgeCategory <- NA
df_senior_bifido[,"AgeCategory"] <- "Senior"

summary_df_studies <- metadata_45809 %>%
  group_by(study_name, age_category,seq_type, cohort) %>%
  summarise(samples = n(), .groups = "drop")
summary_df_studies <- as.data.frame(summary_df_studies)

summary_df_studies[101,"cohort"] <- "UrbanRuralMixed"
summary_df_studies[114,"cohort"] <- "RuralTribal"

infant_study_details <- summary_df_studies[summary_df_studies$age_category == 'infant',]
df_infant_bifido[,"SequencingType"] <- infant_study_details$seq_type
df_infant_bifido[,"CohortType"] <- infant_study_details$cohort

adult_study_details <- summary_df_studies[summary_df_studies$age_category == 'adult',]
adult_study_details <- adult_study_details[-c(62,70),]
df_adult_bifido <- df_adult_bifido[adult_study_details$study_name,]
df_adult_bifido[,"SequencingType"] <- adult_study_details$seq_type
df_adult_bifido[,"CohortType"] <- adult_study_details$cohort

senior_study_details <- summary_df_studies[summary_df_studies$age_category == 'senior',]
df_senior_bifido <- df_senior_bifido[senior_study_details$study_name,]
df_senior_bifido[,"SequencingType"] <- senior_study_details$seq_type
df_senior_bifido[,"CohortType"] <- senior_study_details$cohort

rownames(df_infant_bifido) <- paste0("infant:",rownames(df_infant_bifido))
rownames(df_adult_bifido) <- paste0("adult:",rownames(df_adult_bifido))
rownames(df_senior_bifido) <- paste0("senior:",rownames(df_senior_bifido))

df_combined_bifido <- as.data.frame(rbind(df_infant_bifido[,c(MajorBifsOrdered,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType","AgeCategory")],df_adult_bifido[,c(MajorBifsOrdered,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType","AgeCategory")],df_senior_bifido[,c(MajorBifsOrdered,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType","AgeCategory")]))


pcoStudyPrevEuclidean <- dudi.pco(vegdist(df_combined_bifido[,MajorBifsOrdered],method="euclidean"),scannf=FALSE,nf=3)

colorInfant <- rgb(169,209,142,max=255)
colorAdult <- rgb(157,195,230,max=255)
colorSenior <- rgb(255,217,102,max=255)

s.class(pcoStudyPrevEuclidean$li,as.factor(df_combined_bifido$AgeCategory),col=c(colorInfant,colorAdult,colorSenior),plabels.col="black",plabels.cex=1.5)

s.class(pcoStudyPrevEuclidean$li,as.factor(df_combined_bifido$SequencingType),col=c("Pink","Orange"),plabels.col="black",plabels.cex=1.5)

s.class(pcoStudyPrevEuclidean$li[paste0("adult:",adult_study_details$study_name),],as.factor(df_combined_bifido[paste0("adult:",adult_study_details$study_name),"CohortType"]),col=c("darkolivegreen1","lightslateblue","hotpink1"),plabels.col="black",plabels.cex=1.2)

#Infant Bifido Plot
temp_mat <- as.matrix(df_infant_bifido[,MajorBifsOrdered])
hmpInfantStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row",col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantImagePrevalence.pdf",height=24,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/InfantImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

#Adult Bifido Plot
temp_mat <- as.matrix(df_adult_bifido[,MajorBifsOrdered])
hmpAdultStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row",col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultImagePrevalence.pdf",height=155.2,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/AdultImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

#Senior Bifido Plot
temp_mat <- as.matrix(df_senior_bifido[,MajorBifsOrdered])
hmpSeniorStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row", col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorImagePrevalence.pdf",height=59.2,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/SeniorImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

# Infant Bifido Prevalence
temp_mat <- as.matrix(df_infant_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/InfantOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/InfantSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantCohortType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/InfantCohortType.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/InfantNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido Prevalence
temp_mat <- as.matrix(df_adult_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/AdultOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/AdultSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultCohortType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/AdultCohortType.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/AdultNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido Prevalence
temp_mat <- as.matrix(df_senior_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/SeniorOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/SeniorSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorCohortType.pdf")
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/SeniorCohortType.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="../Data/SeniorNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

infant_study_vector <- rev(colnames(hmpInfantStudyBifido$carpet))
adult_study_vector <- rev(colnames(hmpAdultStudyBifido$carpet))
senior_study_vector <- rev(colnames(hmpSeniorStudyBifido$carpet))

save(infant_study_vector,adult_study_vector,senior_study_vector,file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\StudyVector.RData")

save.image("../Data/Bifidotypes_Batch2.RData")
