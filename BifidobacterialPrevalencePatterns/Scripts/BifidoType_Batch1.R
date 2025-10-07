
load("/data/PrevalenceAnalysis_Input1.RData")
load("/data/PrevalenceAnalysis_Input2.RData")
load("/data/metadata.RData")


library(dplyr)
library(pcaPP)
library(psych)
library(vegan)

SpeciesProfile_All <- as.data.frame(rbind(PrevalenceAnalysis_Input1,PrevalenceAnalysis_Input2))

spdata_45809 <- SpeciesProfile_All
metadata_45809 <- metadata

AllBifidoSpecies <- grep("Bifidobacterium",colnames(spdata_45809),value=TRUE)

spdata_Bifs_rel_45809 <- apply(spdata_45809[,AllBifidoSpecies],2,function(x)(ifelse(is.na(x),0,x)))

spdata_Bifs_detect_45809 <- apply(spdata_Bifs_rel_45809,2,function(x)(ifelse(x>0.0001,1,0)))

studies_45809 <- unique(metadata_45809$study_name)

spdata_45809 <- spdata_45809[intersect(rownames(spdata_45809),rownames(metadata_45809)),]

metadata_45809 <- metadata_45809[intersect(rownames(spdata_45809),rownames(metadata_45809)),]

infant <- rownames(metadata_45809[metadata_45809$age_category == "infant",])

adult <- rownames(metadata_45809[metadata_45809$age_category == "adult",])

senior <- rownames(metadata_45809[metadata_45809$age_category == "senior",])

AllSpecies <- colnames(spdata_45809)

spdata_45809$study_name <- metadata_45809$study_name

infant_studies <- unique(metadata_45809[infant,"study_name"])

adult_studies <- unique(metadata_45809[adult,"study_name"])

senior_studies <- unique(metadata_45809[senior,"study_name"])

print("Infants Overall Species Detection")

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

infant_species_detection <- compute_detection(spdata_45809[infant,],AllSpecies,"study_name",infant_studies)
infant_select_species <- names(which(apply(infant_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(infant_species_detection)>=0.33))
print("Adult Overall Species Detection")

adult_species_detection <- compute_detection(spdata_45809[adult,],AllSpecies,"study_name",adult_studies)
adult_select_species <- names(which(apply(adult_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(adult_species_detection)>=0.33))

print("Senior Overall Species Detection")

senior_species_detection <- compute_detection(spdata_45809[senior,],AllSpecies,"study_name",senior_studies)
senior_select_species <- names(which(apply(senior_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(senior_species_detection)>=0.33))


infant_select_non_bif_species <- setdiff(infant_select_species,AllBifidoSpecies)
adult_select_non_bif_species <- setdiff(adult_select_species,AllBifidoSpecies)
senior_select_non_bif_species <- setdiff(senior_select_species,AllBifidoSpecies)

infant_studies_sample_numbers <- table(metadata_45809[infant,"study_name"])
infant_select_studies <- names(which(infant_studies_sample_numbers>=20))

adult_studies_sample_numbers <- table(metadata_45809[adult,"study_name"])
adult_select_studies <- names(which(adult_studies_sample_numbers>=20))

senior_studies_sample_numbers <- table(metadata_45809[senior,"study_name"])
senior_select_studies <- names(which(senior_studies_sample_numbers>=20))

MajorBifidoSpecies <- c("Bifidobacterium_adolescentis","Bifidobacterium_animalis","Bifidobacterium_bifidum","Bifidobacterium_breve","Bifidobacterium_catenulatum","Bifidobacterium_dentium","Bifidobacterium_longum","Bifidobacterium_pseudocatenulatum")

print("Infant Guts")

df_infant_bifido <- as.data.frame(matrix(NA,length(infant_select_studies),13))
rownames(df_infant_bifido) <- infant_select_studies
colnames(df_infant_bifido) <- c(MajorBifidoSpecies,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType")

for(i in 1:length(infant_select_studies))
{
	study_name <- infant_select_studies[i]
	print(study_name)
	study_rows <- intersect(infant,rownames(metadata_45809)[metadata_45809$study_name == study_name])
	df_infant_bifido[study_name,1:8] <- apply(spdata_45809[study_rows,MajorBifidoSpecies],2,function(x)(length(x[(x>0.0001)&(!is.na(x))])))/length(study_rows)
	df_infant_bifido[study_name,"BifPrevalence"] <- length(which(apply(spdata_45809[study_rows,AllBifidoSpecies],1,function(x)(length(x[(x>0.0001)&(!is.na(x))])))>0))/length(study_rows)
	temp_spdata <- spdata_45809[study_rows,infant_select_non_bif_species]
	temp_spdata <- apply(temp_spdata,2,function(x)(ifelse(is.na(x),0,x)))
	temp_spdata <- temp_spdata[rowSums(temp_spdata)>0,colSums(temp_spdata)>0]
	temp_spdata <- temp_spdata/rowSums(temp_spdata)
	BifidoDetection <- apply(spdata_45809[rownames(temp_spdata),AllBifidoSpecies],1,function(x)(length(x[x>=0.0001])))
	print("bray influence")
	tempAdonisBray <- adonis2(vegdist(temp_spdata,method="bray")~BifidoDetection)
	df_infant_bifido[study_name,"BrayAssocBifDetection"] <- tempAdonisBray$Pr[1]
	print("kendall influence")
	tempAdonisKendall <- adonis2(as.dist(1-cor.fk(t(temp_spdata))/2)~BifidoDetection)
	df_infant_bifido[study_name,"KendallAssocBifDetection"] <- tempAdonisKendall$Pr[1]
}



spdata_adult_45809 <- spdata_45809[adult,]
metadata_adult_45809 <- metadata_45809[adult,]

df_adult_bifido <- as.data.frame(matrix(NA,length(adult_select_studies),13))
rownames(df_adult_bifido) <- adult_select_studies
colnames(df_adult_bifido) <- c(MajorBifidoSpecies,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType")

for(i in 1:length(adult_select_studies))
{
  study_name <- adult_select_studies[i]
  print(study_name)
  study_rows <- rownames(metadata_adult_45809)[metadata_adult_45809$study_name == study_name]
  df_adult_bifido[study_name,1:8] <- apply(spdata_adult_45809[study_rows,MajorBifidoSpecies],2,function(x)(length(x[(x>0.0001)&(!is.na(x))])))/length(study_rows)
  df_adult_bifido[study_name,"BifPrevalence"] <- length(which(apply(spdata_adult_45809[study_rows,AllBifidoSpecies],1,function(x)(length(x[(x>0.0001)&(!is.na(x))])))>0))/length(study_rows)
  temp_spdata <- spdata_adult_45809[study_rows,adult_select_non_bif_species]
  temp_spdata <- apply(temp_spdata,2,function(x)(ifelse(is.na(x),0,x)))
  temp_spdata <- temp_spdata[rowSums(temp_spdata)>0,colSums(temp_spdata)>0]
  temp_spdata <- temp_spdata/rowSums(temp_spdata)
  BifidoDetection <- apply(spdata_adult_45809[rownames(temp_spdata),AllBifidoSpecies],1,function(x)(length(x[x>=0.0001])))
  print("bray influence")
  tempAdonisBray <- adonis2(vegdist(temp_spdata,method="bray")~BifidoDetection)
  df_adult_bifido[study_name,"BrayAssocBifDetection"] <- tempAdonisBray$Pr[1]
  print("kendall influence")
  tempAdonisKendall <- adonis2(as.dist(1-cor.fk(t(temp_spdata))/2)~BifidoDetection)
  df_adult_bifido[study_name,"KendallAssocBifDetection"] <- tempAdonisKendall$Pr[1]
}



df_senior_bifido <- as.data.frame(matrix(NA,length(senior_select_studies),13))
rownames(df_senior_bifido) <- senior_select_studies
colnames(df_senior_bifido) <- c(MajorBifidoSpecies,"BifPrevalence","BrayAssocBifDetection","KendallAssocBifDetection","SequencingType","CohortType")

for(i in 1:length(senior_select_studies))
{
  study_name <- senior_select_studies[i]
  print(study_name)
  study_rows <- intersect(senior,rownames(metadata_45809)[metadata_45809$study_name == study_name])
  df_senior_bifido[study_name,1:8] <- apply(spdata_45809[study_rows,MajorBifidoSpecies],2,function(x)(length(x[(x>0.0001)&(!is.na(x))])))/length(study_rows)
  df_senior_bifido[study_name,"BifPrevalence"] <- length(which(apply(spdata_45809[study_rows,AllBifidoSpecies],1,function(x)(length(x[(x>0.0001)&(!is.na(x))])))>0))/length(study_rows)
  temp_spdata <- spdata_45809[study_rows,senior_select_non_bif_species]
  temp_spdata <- apply(temp_spdata,2,function(x)(ifelse(is.na(x),0,x)))
  temp_spdata <- temp_spdata[rowSums(temp_spdata)>0,colSums(temp_spdata)>0]
  temp_spdata <- temp_spdata/rowSums(temp_spdata)
  BifidoDetection <- apply(spdata_45809[rownames(temp_spdata),AllBifidoSpecies],1,function(x)(length(x[x>=0.0001])))
  print("bray influence")
  tempAdonisBray <- adonis2(vegdist(temp_spdata,method="bray")~BifidoDetection)
  df_senior_bifido[study_name,"BrayAssocBifDetection"] <- tempAdonisBray$Pr[1]
  print("kendall influence")
  tempAdonisKendall <- adonis2(as.dist(1-cor.fk(t(temp_spdata))/2)~BifidoDetection)
  df_senior_bifido[study_name,"KendallAssocBifDetection"] <- tempAdonisKendall$Pr[1]
}




save.image("../results/Bifidotypes_Batch1.RData")



