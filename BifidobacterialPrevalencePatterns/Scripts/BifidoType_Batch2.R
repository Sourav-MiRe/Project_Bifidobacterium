

MajorBifsOrdered <- colnames(df_infant_bifido)[order(apply(rbind(df_infant_bifido[,1:8],df_adult_bifido[,1:8],df_senior_bifido[,1:8]),2,function(x)(length(x[x>=0.3]))))]

rownames(df_infant_bifido) <- paste0("infant:",rownames(df_infant_bifido))
rownames(df_adult_bifido) <- paste0("adult:",rownames(df_adult_bifido))
rownames(df_senior_bifido) <- paste0("senior:",rownames(df_senior_bifido))

df_combined_bifido <- as.data.frame(rbind(df_infant_bifido[,c(MajorBifsOrdered,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType")],df_adult_bifido[,c(MajorBifsOrdered,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType")],df_senior_bifido[,c(MajorBifsOrdered,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType")]))

df_combined_bifido$AgeCategory <- NA
df_combined_bifido[rownames(df_infant_bifido),"AgeCategory"] <- "Infant"
df_combined_bifido[rownames(df_adult_bifido),"AgeCategory"] <- "Adult"
df_combined_bifido[rownames(df_senior_bifido),"AgeCategory"] <- "Senior"

#adult_study_details <- read.table("G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\study_list_adult.tsv",sep="\t",row.names=1,header=TRUE)
adult_study_details <- read.csv("notebooks\\figure1\\adult_study_details.csv", row.names = "study_name")
rownames(adult_study_details) <- paste0("adult:",rownames(adult_study_details))
 
#senior_study_details <- read.table("G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\study_list_senior.tsv",sep="\t",row.names=1,header=TRUE)
senior_study_details <- read.csv("notebooks\\figure1\\senior_study_details.csv", row.names = "study_name")
rownames(senior_study_details) <- paste0("senior:",rownames(senior_study_details))
 
#infant_study_details <- read.table("G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\study_list_newborn.tsv",sep="\t",row.names=1,header=TRUE)
infant_study_details <- read.csv("notebooks\\figure1\\infant_study_details.csv", row.names = "study_name")
rownames(infant_study_details) <- paste0("infant:",rownames(infant_study_details))

df_combined_bifido[rownames(senior_study_details),"SequencingType"] <- senior_study_details$seq_type
df_combined_bifido[rownames(adult_study_details),"SequencingType"] <- adult_study_details$seq_type
df_combined_bifido[rownames(infant_study_details),"SequencingType"] <- infant_study_details$seq_type

df_combined_bifido[rownames(senior_study_details),"CohortType"] <- senior_study_details$cohort
df_combined_bifido[rownames(adult_study_details),"CohortType"] <- adult_study_details$cohort
df_combined_bifido[rownames(infant_study_details),"CohortType"] <- infant_study_details$cohort

pcoStudyPrevEuclidean <- dudi.pco(vegdist(df_combined_bifido[,MajorBifsOrdered],method="euclidean"),scannf=FALSE,nf=3)

colorInfant <- rgb(169,209,142,max=255)
colorAdult <- rgb(157,195,230,max=255)
colorSenior <- rgb(255,217,102,max=255)

s.class(pcoStudyPrevEuclidean$li,as.factor(df_combined_bifido$AgeCategory),col=c(colorInfant,colorAdult,colorSenior),plabels.col="black",plabels.cex=1.5)

s.class(pcoStudyPrevEuclidean$li,as.factor(df_combined_bifido$SequencingType),col=c("Pink","Orange"),plabels.col="black",plabels.cex=1.5)

s.class(pcoStudyPrevEuclidean$li[rownames(adult_study_details),],as.factor(df_combined_bifido[rownames(adult_study_details),"CohortType"]),col=c("darkolivegreen1","lightslateblue","hotpink1"),plabels.col="black",plabels.cex=1.2)

#Infant Bifido Plot
temp_mat <- as.matrix(df_infant_bifido[,MajorBifsOrdered])
hmpInfantStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row",col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantImagePrevalence.pdf",height=24,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\InfantImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

#Adult Bifido Plot
temp_mat <- as.matrix(df_adult_bifido[,MajorBifsOrdered])
hmpAdultStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row",col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultImagePrevalence.pdf",height=155.2,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\AdultImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

#Senior Bifido Plot
temp_mat <- as.matrix(df_senior_bifido[,MajorBifsOrdered])
hmpSeniorStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row", col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorImagePrevalence.pdf",height=59.2,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\SeniorImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

# Infant Bifido Prevalence
temp_mat <- as.matrix(df_infant_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\InfantOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\InfantSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantCohortType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\InfantCohortType.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\InfantNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido Prevalence
temp_mat <- as.matrix(df_adult_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\AdultOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\AdultSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultCohortType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\AdultCohortType.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\AdultNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido Prevalence
temp_mat <- as.matrix(df_senior_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\SeniorOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\SeniorSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorCohortType.pdf")
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\SeniorCohortType.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="notebooks\\figure1\\SeniorNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

infant_study_vector <- rev(colnames(hmpInfantStudyBifido$carpet))
adult_study_vector <- rev(colnames(hmpAdultStudyBifido$carpet))
senior_study_vector <- rev(colnames(hmpSeniorStudyBifido$carpet))
save(infant_study_vector,adult_study_vector,senior_study_vector,file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\StudyVector.RData")

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