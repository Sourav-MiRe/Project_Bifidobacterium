print("Machine learning derived feature association scores read")
load("G:\\.shortcut-targets-by-id\\1mJ9DmgZE4NfNNbH72zZNyUWw4tVWKsRf\\Bif_Manuscript\\ml_association.RData")

print("Association Data frames with Taxa Names modified read")
load("C:\\Projects\\Bif_Manuscript\\OriginalAssociationScores_Modified.RData")

library(psych)
library(gplots)
library(RColorBrewer)
library(ade4)
library(vegan)
library(ggplot2)
library(ggrepel)

print("adult:wgs")

df_association_adult_detection_wgs$association$scores <- apply(df_association_adult_detection_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_detection_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_longum_wgs$association$scores <- apply(df_association_adult_longum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_longum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_adolescentis_wgs$association$scores <- apply(df_association_adult_adolescentis_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_adolescentis_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_pseudocatenulatum_wgs$association$scores <- apply(df_association_adult_pseudocatenulatum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_pseudocatenulatum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_bifidum_wgs$association$scores <- apply(df_association_adult_bifidum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_bifidum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_breve_wgs$association$scores <- apply(df_association_adult_breve_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_breve_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_dentium_wgs$association$scores <- apply(df_association_adult_dentium_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_dentium_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_catenulatum_wgs$association$scores <- apply(df_association_adult_catenulatum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_catenulatum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_animalis_wgs$association$scores <- apply(df_association_adult_animalis_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_animalis_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

union_species_adult_wgs <- names(which(table(c(rownames(df_association_adult_animalis_wgs$association),rownames(df_association_adult_catenulatum_wgs$association),rownames(df_association_adult_breve_wgs$association),rownames(df_association_adult_bifidum_wgs$association),rownames(df_association_adult_dentium_wgs$association),rownames(df_association_adult_pseudocatenulatum_wgs$association),rownames(df_association_adult_adolescentis_wgs$association),rownames(df_association_adult_longum_wgs$association),rownames(df_association_adult_detection_wgs$association)))==9))

df_association_all_adult_wgs <- data.frame("detection"=df_association_adult_detection_wgs$association[union_species_adult_wgs,"scores"],"longum"=df_association_adult_longum_wgs$association[union_species_adult_wgs,"scores"],"adolescentis"=df_association_adult_adolescentis_wgs$association[union_species_adult_wgs,"scores"],"pseudocatenulatum"=df_association_adult_pseudocatenulatum_wgs$association[union_species_adult_wgs,"scores"],"dentium"=df_association_adult_dentium_wgs$association[union_species_adult_wgs,"scores"],"bifidum"=df_association_adult_bifidum_wgs$association[union_species_adult_wgs,"scores"],"breve"=df_association_adult_breve_wgs$association[union_species_adult_wgs,"scores"],"catenulatum"=df_association_adult_catenulatum_wgs$association[union_species_adult_wgs,"scores"],"animalis"=df_association_adult_animalis_wgs$association[union_species_adult_wgs,"scores"],row.names=union_species_adult_wgs)

select_species_adult_wgs <- names(which(apply(df_association_all_adult_wgs,1,function(x)(length(x[abs(x)>=0.20])))>0))

df_association_all_adult_wgs <- df_association_all_adult_wgs[select_species_adult_wgs,]


print("senior:wgs")

df_association_senior_detection_wgs$association$scores <- apply(df_association_senior_detection_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_detection_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_longum_wgs$association$scores <- apply(df_association_senior_longum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_longum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_adolescentis_wgs$association$scores <- apply(df_association_senior_adolescentis_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_adolescentis_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_pseudocatenulatum_wgs$association$scores <- apply(df_association_senior_pseudocatenulatum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_pseudocatenulatum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_bifidum_wgs$association$scores <- apply(df_association_senior_bifidum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_bifidum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_breve_wgs$association$scores <- apply(df_association_senior_breve_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_breve_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_dentium_wgs$association$scores <- apply(df_association_senior_dentium_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_dentium_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_catenulatum_wgs$association$scores <- apply(df_association_senior_catenulatum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_catenulatum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_animalis_wgs$association$scores <- apply(df_association_senior_animalis_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_animalis_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

union_species_senior_wgs <- names(which(table(c(rownames(df_association_senior_animalis_wgs$association),rownames(df_association_senior_catenulatum_wgs$association),rownames(df_association_senior_breve_wgs$association),rownames(df_association_senior_bifidum_wgs$association),rownames(df_association_senior_dentium_wgs$association),rownames(df_association_senior_pseudocatenulatum_wgs$association),rownames(df_association_senior_adolescentis_wgs$association),rownames(df_association_senior_longum_wgs$association),rownames(df_association_senior_detection_wgs$association)))==9))

df_association_all_senior_wgs <- data.frame("detection"=df_association_senior_detection_wgs$association[union_species_senior_wgs,"scores"],"longum"=df_association_senior_longum_wgs$association[union_species_senior_wgs,"scores"],"adolescentis"=df_association_senior_adolescentis_wgs$association[union_species_senior_wgs,"scores"],"pseudocatenulatum"=df_association_senior_pseudocatenulatum_wgs$association[union_species_senior_wgs,"scores"],"dentium"=df_association_senior_dentium_wgs$association[union_species_senior_wgs,"scores"],"bifidum"=df_association_senior_bifidum_wgs$association[union_species_senior_wgs,"scores"],"breve"=df_association_senior_breve_wgs$association[union_species_senior_wgs,"scores"],"catenulatum"=df_association_senior_catenulatum_wgs$association[union_species_senior_wgs,"scores"],"animalis"=df_association_senior_animalis_wgs$association[union_species_senior_wgs,"scores"],row.names=union_species_senior_wgs)

select_species_senior_wgs <- names(which(apply(df_association_all_senior_wgs,1,function(x)(length(x[abs(x)>=0.20])))>0))

df_association_all_senior_wgs <- df_association_all_senior_wgs[select_species_senior_wgs,]


print("infant:wgs")

df_association_infant_detection_wgs$association$scores <- apply(df_association_infant_detection_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_detection_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_longum_wgs$association$scores <- apply(df_association_infant_longum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_longum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_adolescentis_wgs$association$scores <- apply(df_association_infant_adolescentis_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_adolescentis_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_pseudocatenulatum_wgs$association$scores <- apply(df_association_infant_pseudocatenulatum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_pseudocatenulatum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_bifidum_wgs$association$scores <- apply(df_association_infant_bifidum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_bifidum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_breve_wgs$association$scores <- apply(df_association_infant_breve_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_breve_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_dentium_wgs$association$scores <- apply(df_association_infant_dentium_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_dentium_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_catenulatum_wgs$association$scores <- apply(df_association_infant_catenulatum_wgs$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_catenulatum_wgs$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

union_species_infant_wgs <- names(which(table(c(rownames(df_association_infant_catenulatum_wgs$association),rownames(df_association_infant_breve_wgs$association),rownames(df_association_infant_bifidum_wgs$association),rownames(df_association_infant_dentium_wgs$association),rownames(df_association_infant_pseudocatenulatum_wgs$association),rownames(df_association_infant_adolescentis_wgs$association),rownames(df_association_infant_longum_wgs$association),rownames(df_association_infant_detection_wgs$association)))==8))

df_association_all_infant_wgs <- data.frame("detection"=df_association_infant_detection_wgs$association[union_species_infant_wgs,"scores"],"longum"=df_association_infant_longum_wgs$association[union_species_infant_wgs,"scores"],"adolescentis"=df_association_infant_adolescentis_wgs$association[union_species_infant_wgs,"scores"],"pseudocatenulatum"=df_association_infant_pseudocatenulatum_wgs$association[union_species_infant_wgs,"scores"],"dentium"=df_association_infant_dentium_wgs$association[union_species_infant_wgs,"scores"],"bifidum"=df_association_infant_bifidum_wgs$association[union_species_infant_wgs,"scores"],"breve"=df_association_infant_breve_wgs$association[union_species_infant_wgs,"scores"],"catenulatum"=df_association_infant_catenulatum_wgs$association[union_species_infant_wgs,"scores"],row.names=union_species_infant_wgs)

select_species_infant_wgs <- names(which(apply(df_association_all_infant_wgs,1,function(x)(length(x[abs(x)>=0.20])))>0))

df_association_all_infant_wgs <- df_association_all_infant_wgs[select_species_infant_wgs,]

print("Cohort Metadata read")
load("C:\\Projects\\Bif_Manuscript\\cohort_metadata.RData")

RuralTribal_cohort <- rownames(cohort_metadata)[cohort_metadata$`Cohort-Type`=="RuralTribal"]
IndustrializedUrban_cohort <- rownames(cohort_metadata)[cohort_metadata$`Cohort-Type`=="IndustrializedUrban"]
UrbanRuralMixed_cohort <- rownames(cohort_metadata)[cohort_metadata$`Cohort-Type`=="UrbanRuralMixed"]

print("RuralTribal")
df_detection_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_detection_wgs$dir[intersect(rownames(df_association_adult_detection_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_detection_wgs$dir[intersect(rownames(df_association_adult_detection_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_detection_wgs$dir),RuralTribal_cohort)))

df_detection_adult_wgs_RuralTribal$score <- apply(df_detection_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_detection_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_longum_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_longum_wgs$dir[intersect(rownames(df_association_adult_longum_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_longum_wgs$dir[intersect(rownames(df_association_adult_longum_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_longum_wgs$dir),RuralTribal_cohort)))

df_longum_adult_wgs_RuralTribal$score <- apply(df_longum_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_longum_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_adolescentis_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_adolescentis_wgs$dir[intersect(rownames(df_association_adult_adolescentis_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_adolescentis_wgs$dir[intersect(rownames(df_association_adult_adolescentis_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_adolescentis_wgs$dir),RuralTribal_cohort)))

df_adolescentis_adult_wgs_RuralTribal$score <- apply(df_adolescentis_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_adolescentis_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_pseudocatenulatum_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_pseudocatenulatum_wgs$dir[intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_pseudocatenulatum_wgs$dir[intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),RuralTribal_cohort)))

df_pseudocatenulatum_adult_wgs_RuralTribal$score <- apply(df_pseudocatenulatum_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_pseudocatenulatum_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_dentium_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_dentium_wgs$dir[intersect(rownames(df_association_adult_dentium_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_dentium_wgs$dir[intersect(rownames(df_association_adult_dentium_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_dentium_wgs$dir),RuralTribal_cohort)))

df_dentium_adult_wgs_RuralTribal$score <- apply(df_dentium_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_dentium_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_bifidum_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_bifidum_wgs$dir[intersect(rownames(df_association_adult_bifidum_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_bifidum_wgs$dir[intersect(rownames(df_association_adult_bifidum_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_bifidum_wgs$dir),RuralTribal_cohort)))

df_bifidum_adult_wgs_RuralTribal$score <- apply(df_bifidum_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_bifidum_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_breve_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_breve_wgs$dir[intersect(rownames(df_association_adult_breve_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_breve_wgs$dir[intersect(rownames(df_association_adult_breve_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_breve_wgs$dir),RuralTribal_cohort)))

df_breve_adult_wgs_RuralTribal$score <- apply(df_breve_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_breve_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_catenulatum_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_catenulatum_wgs$dir[intersect(rownames(df_association_adult_catenulatum_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_catenulatum_wgs$dir[intersect(rownames(df_association_adult_catenulatum_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_catenulatum_wgs$dir),RuralTribal_cohort)))

df_catenulatum_adult_wgs_RuralTribal$score <- apply(df_catenulatum_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_catenulatum_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_animalis_adult_wgs_RuralTribal <- data.frame("SP"=apply(df_association_adult_animalis_wgs$dir[intersect(rownames(df_association_adult_animalis_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_animalis_wgs$dir[intersect(rownames(df_association_adult_animalis_wgs$dir),RuralTribal_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_animalis_wgs$dir),RuralTribal_cohort)))

df_animalis_adult_wgs_RuralTribal$score <- apply(df_animalis_adult_wgs_RuralTribal,1,function(x)((x[1]-x[2])/x[3])) * apply(df_animalis_adult_wgs_RuralTribal,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))
df_animalis_adult_wgs_RuralTribal <- as.data.frame(apply(df_animalis_adult_wgs_RuralTribal,2,function(x)(ifelse(is.nan(x),0,x))))

union_species_adult_wgs_RuralTribal <- names(table(c(rownames(df_animalis_adult_wgs_RuralTribal),rownames(df_catenulatum_adult_wgs_RuralTribal),rownames(df_breve_adult_wgs_RuralTribal),rownames(df_bifidum_adult_wgs_RuralTribal),rownames(df_dentium_adult_wgs_RuralTribal),rownames(df_pseudocatenulatum_adult_wgs_RuralTribal),rownames(df_adolescentis_adult_wgs_RuralTribal),rownames(df_longum_adult_wgs_RuralTribal),rownames(df_detection_adult_wgs_RuralTribal))))

df_association_adult_wgs_RuralTribal <- data.frame("detection"=df_detection_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],"longum"=df_longum_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],"adolescentis"=df_adolescentis_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],"pseudocatenulatum"=df_pseudocatenulatum_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],"dentium"=df_dentium_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],"bifidum"=df_bifidum_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],"breve"=df_breve_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],"catenulatum"=df_catenulatum_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],"animalis"=df_animalis_adult_wgs_RuralTribal[union_species_adult_wgs_RuralTribal,"score"],row.names=union_species_adult_wgs_RuralTribal)

select_adult_wgs_RuralTribal <- names(which(apply(df_association_adult_wgs_RuralTribal,1,function(x)(length(x[abs(x)>=0.20])))>1))

df_association_adult_wgs_RuralTribal <- df_association_adult_wgs_RuralTribal[select_adult_wgs_RuralTribal,]

print("UrbanRuralMixed")
df_detection_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_detection_wgs$dir[intersect(rownames(df_association_adult_detection_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_detection_wgs$dir[intersect(rownames(df_association_adult_detection_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_detection_wgs$dir),UrbanRuralMixed_cohort)))

df_detection_adult_wgs_UrbanRuralMixed$score <- apply(df_detection_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_detection_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_longum_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_longum_wgs$dir[intersect(rownames(df_association_adult_longum_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_longum_wgs$dir[intersect(rownames(df_association_adult_longum_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_longum_wgs$dir),UrbanRuralMixed_cohort)))

df_longum_adult_wgs_UrbanRuralMixed$score <- apply(df_longum_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_longum_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_adolescentis_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_adolescentis_wgs$dir[intersect(rownames(df_association_adult_adolescentis_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_adolescentis_wgs$dir[intersect(rownames(df_association_adult_adolescentis_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_adolescentis_wgs$dir),UrbanRuralMixed_cohort)))

df_adolescentis_adult_wgs_UrbanRuralMixed$score <- apply(df_adolescentis_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_adolescentis_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_pseudocatenulatum_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_pseudocatenulatum_wgs$dir[intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_pseudocatenulatum_wgs$dir[intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),UrbanRuralMixed_cohort)))

df_pseudocatenulatum_adult_wgs_UrbanRuralMixed$score <- apply(df_pseudocatenulatum_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_pseudocatenulatum_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_dentium_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_dentium_wgs$dir[intersect(rownames(df_association_adult_dentium_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_dentium_wgs$dir[intersect(rownames(df_association_adult_dentium_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_dentium_wgs$dir),UrbanRuralMixed_cohort)))

df_dentium_adult_wgs_UrbanRuralMixed$score <- apply(df_dentium_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_dentium_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_bifidum_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_bifidum_wgs$dir[intersect(rownames(df_association_adult_bifidum_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_bifidum_wgs$dir[intersect(rownames(df_association_adult_bifidum_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_bifidum_wgs$dir),UrbanRuralMixed_cohort)))

df_bifidum_adult_wgs_UrbanRuralMixed$score <- apply(df_bifidum_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_bifidum_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_breve_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_breve_wgs$dir[intersect(rownames(df_association_adult_breve_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_breve_wgs$dir[intersect(rownames(df_association_adult_breve_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_breve_wgs$dir),UrbanRuralMixed_cohort)))

df_breve_adult_wgs_UrbanRuralMixed$score <- apply(df_breve_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_breve_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_catenulatum_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_catenulatum_wgs$dir[intersect(rownames(df_association_adult_catenulatum_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_catenulatum_wgs$dir[intersect(rownames(df_association_adult_catenulatum_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_catenulatum_wgs$dir),UrbanRuralMixed_cohort)))

df_catenulatum_adult_wgs_UrbanRuralMixed$score <- apply(df_catenulatum_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_catenulatum_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_animalis_adult_wgs_UrbanRuralMixed <- data.frame("SP"=apply(df_association_adult_animalis_wgs$dir[intersect(rownames(df_association_adult_animalis_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_animalis_wgs$dir[intersect(rownames(df_association_adult_animalis_wgs$dir),UrbanRuralMixed_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_animalis_wgs$dir),UrbanRuralMixed_cohort)))

df_animalis_adult_wgs_UrbanRuralMixed$score <- apply(df_animalis_adult_wgs_UrbanRuralMixed,1,function(x)((x[1]-x[2])/x[3])) * apply(df_animalis_adult_wgs_UrbanRuralMixed,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))
df_animalis_adult_wgs_UrbanRuralMixed <- as.data.frame(apply(df_animalis_adult_wgs_UrbanRuralMixed,2,function(x)(ifelse(is.nan(x),0,x))))

union_species_adult_wgs_UrbanRuralMixed <- names(table(c(rownames(df_animalis_adult_wgs_UrbanRuralMixed),rownames(df_catenulatum_adult_wgs_UrbanRuralMixed),rownames(df_breve_adult_wgs_UrbanRuralMixed),rownames(df_bifidum_adult_wgs_UrbanRuralMixed),rownames(df_dentium_adult_wgs_UrbanRuralMixed),rownames(df_pseudocatenulatum_adult_wgs_UrbanRuralMixed),rownames(df_adolescentis_adult_wgs_UrbanRuralMixed),rownames(df_longum_adult_wgs_UrbanRuralMixed),rownames(df_detection_adult_wgs_UrbanRuralMixed))))

df_association_adult_wgs_UrbanRuralMixed <- data.frame("detection"=df_detection_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],"longum"=df_longum_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],"adolescentis"=df_adolescentis_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],"pseudocatenulatum"=df_pseudocatenulatum_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],"dentium"=df_dentium_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],"bifidum"=df_bifidum_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],"breve"=df_breve_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],"catenulatum"=df_catenulatum_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],"animalis"=df_animalis_adult_wgs_UrbanRuralMixed[union_species_adult_wgs_UrbanRuralMixed,"score"],row.names=union_species_adult_wgs_UrbanRuralMixed)

select_adult_wgs_UrbanRuralMixed <- names(which(apply(df_association_adult_wgs_UrbanRuralMixed,1,function(x)(length(x[abs(x)>=0.20])))>1))

df_association_adult_wgs_UrbanRuralMixed <- df_association_adult_wgs_UrbanRuralMixed[select_adult_wgs_UrbanRuralMixed,]

print("IndustrializedUrban")
df_detection_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_detection_wgs$dir[intersect(rownames(df_association_adult_detection_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_detection_wgs$dir[intersect(rownames(df_association_adult_detection_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_detection_wgs$dir),IndustrializedUrban_cohort)))

df_detection_adult_wgs_IndustrializedUrban$score <- apply(df_detection_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_detection_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_longum_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_longum_wgs$dir[intersect(rownames(df_association_adult_longum_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_longum_wgs$dir[intersect(rownames(df_association_adult_longum_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_longum_wgs$dir),IndustrializedUrban_cohort)))

df_longum_adult_wgs_IndustrializedUrban$score <- apply(df_longum_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_longum_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_adolescentis_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_adolescentis_wgs$dir[intersect(rownames(df_association_adult_adolescentis_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_adolescentis_wgs$dir[intersect(rownames(df_association_adult_adolescentis_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_adolescentis_wgs$dir),IndustrializedUrban_cohort)))

df_adolescentis_adult_wgs_IndustrializedUrban$score <- apply(df_adolescentis_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_adolescentis_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_pseudocatenulatum_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_pseudocatenulatum_wgs$dir[intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_pseudocatenulatum_wgs$dir[intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_pseudocatenulatum_wgs$dir),IndustrializedUrban_cohort)))

df_pseudocatenulatum_adult_wgs_IndustrializedUrban$score <- apply(df_pseudocatenulatum_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_pseudocatenulatum_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_dentium_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_dentium_wgs$dir[intersect(rownames(df_association_adult_dentium_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_dentium_wgs$dir[intersect(rownames(df_association_adult_dentium_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_dentium_wgs$dir),IndustrializedUrban_cohort)))

df_dentium_adult_wgs_IndustrializedUrban$score <- apply(df_dentium_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_dentium_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_bifidum_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_bifidum_wgs$dir[intersect(rownames(df_association_adult_bifidum_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_bifidum_wgs$dir[intersect(rownames(df_association_adult_bifidum_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_bifidum_wgs$dir),IndustrializedUrban_cohort)))

df_bifidum_adult_wgs_IndustrializedUrban$score <- apply(df_bifidum_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_bifidum_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_breve_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_breve_wgs$dir[intersect(rownames(df_association_adult_breve_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_breve_wgs$dir[intersect(rownames(df_association_adult_breve_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_breve_wgs$dir),IndustrializedUrban_cohort)))

df_breve_adult_wgs_IndustrializedUrban$score <- apply(df_breve_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_breve_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_catenulatum_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_catenulatum_wgs$dir[intersect(rownames(df_association_adult_catenulatum_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_catenulatum_wgs$dir[intersect(rownames(df_association_adult_catenulatum_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_catenulatum_wgs$dir),IndustrializedUrban_cohort)))

df_catenulatum_adult_wgs_IndustrializedUrban$score <- apply(df_catenulatum_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_catenulatum_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_animalis_adult_wgs_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_animalis_wgs$dir[intersect(rownames(df_association_adult_animalis_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_animalis_wgs$dir[intersect(rownames(df_association_adult_animalis_wgs$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_animalis_wgs$dir),IndustrializedUrban_cohort)))

df_animalis_adult_wgs_IndustrializedUrban$score <- apply(df_animalis_adult_wgs_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_animalis_adult_wgs_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))
df_animalis_adult_wgs_IndustrializedUrban <- as.data.frame(apply(df_animalis_adult_wgs_IndustrializedUrban,2,function(x)(ifelse(is.nan(x),0,x))))

union_species_adult_wgs_IndustrializedUrban <- names(table(c(rownames(df_animalis_adult_wgs_IndustrializedUrban),rownames(df_catenulatum_adult_wgs_IndustrializedUrban),rownames(df_breve_adult_wgs_IndustrializedUrban),rownames(df_bifidum_adult_wgs_IndustrializedUrban),rownames(df_dentium_adult_wgs_IndustrializedUrban),rownames(df_pseudocatenulatum_adult_wgs_IndustrializedUrban),rownames(df_adolescentis_adult_wgs_IndustrializedUrban),rownames(df_longum_adult_wgs_IndustrializedUrban),rownames(df_detection_adult_wgs_IndustrializedUrban))))

df_association_adult_wgs_IndustrializedUrban <- data.frame("detection"=df_detection_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],"longum"=df_longum_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],"adolescentis"=df_adolescentis_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],"pseudocatenulatum"=df_pseudocatenulatum_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],"dentium"=df_dentium_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],"bifidum"=df_bifidum_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],"breve"=df_breve_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],"catenulatum"=df_catenulatum_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],"animalis"=df_animalis_adult_wgs_IndustrializedUrban[union_species_adult_wgs_IndustrializedUrban,"score"],row.names=union_species_adult_wgs_IndustrializedUrban)

select_adult_wgs_IndustrializedUrban <- names(which(apply(df_association_adult_wgs_IndustrializedUrban,1,function(x)(length(x[abs(x)>=0.20])))>1))

df_association_adult_wgs_IndustrializedUrban <- df_association_adult_wgs_IndustrializedUrban[select_adult_wgs_IndustrializedUrban,]

union_all_species_wgs <- names(which(table(c(rownames(df_association_all_adult_wgs),rownames(df_association_all_infant_wgs),rownames(df_association_all_senior_wgs),rownames(df_association_adult_wgs_IndustrializedUrban),rownames(df_association_adult_wgs_UrbanRuralMixed),rownames(df_association_adult_wgs_RuralTribal)))>1))

df_association_all_wgs <- as.data.frame(matrix(0,length(union_all_species_wgs),53))
rownames(df_association_all_wgs) <- union_all_species_wgs
colnames(df_association_all_wgs) <- c("infant_detection","infant_longum","infant_adolescentis","infant_pseudocatenulatum","infant_dentium","infant_bifidum","infant_breve","infant_catenulatum","adult_detection","adult_longum","adult_adolescentis","adult_pseudocatenulatum","adult_dentium","adult_bifidum","adult_breve","adult_catenulatum","adult_animalis","senior_detection","senior_longum","senior_adolescentis","senior_pseudocatenulatum","senior_dentium","senior_bifidum","senior_breve","senior_catenulatum","senior_animalis","IndustrializedUrban_detection","IndustrializedUrban_longum","IndustrializedUrban_adolescentis","IndustrializedUrban_pseudocatenulatum","IndustrializedUrban_dentium","IndustrializedUrban_bifidum","IndustrializedUrban_breve","IndustrializedUrban_catenulatum","IndustrializedUrban_animalis","UrbanRuralMixed_detection","UrbanRuralMixed_longum","UrbanRuralMixed_adolescentis","UrbanRuralMixed_pseudocatenulatum","UrbanRuralMixed_dentium","UrbanRuralMixed_bifidum","UrbanRuralMixed_breve","UrbanRuralMixed_catenulatum","UrbanRuralMixed_animalis","RuralTribal_detection","RuralTribal_longum","RuralTribal_adolescentis","RuralTribal_pseudocatenulatum","RuralTribal_dentium","RuralTribal_bifidum","RuralTribal_breve","RuralTribal_catenulatum","RuralTribal_animalis")

df_association_all_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"infant_detection"] <- df_association_all_infant_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"detection"]

df_association_all_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"infant_longum"] <- df_association_all_infant_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"longum"]

df_association_all_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"infant_adolescentis"] <- df_association_all_infant_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"adolescentis"]

df_association_all_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"infant_pseudocatenulatum"] <- df_association_all_infant_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"pseudocatenulatum"]

df_association_all_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"infant_dentium"] <- df_association_all_infant_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"dentium"]

df_association_all_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"infant_bifidum"] <- df_association_all_infant_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"bifidum"]

df_association_all_wgs[intersect(rownames(df_association_all_wgs),union_all_species_wgs),"infant_breve"] <- df_association_all_infant_wgs[intersect(rownames(df_association_all_wgs),union_all_species_wgs),"breve"]

df_association_all_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"infant_catenulatum"] <- df_association_all_infant_wgs[intersect(rownames(df_association_all_infant_wgs),union_all_species_wgs),"catenulatum"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_detection"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"detection"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_longum"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"longum"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_adolescentis"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adolescentis"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_pseudocatenulatum"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"pseudocatenulatum"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_dentium"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"dentium"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_bifidum"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"bifidum"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_breve"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"breve"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_catenulatum"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"catenulatum"]

df_association_all_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"adult_animalis"] <- df_association_all_adult_wgs[intersect(rownames(df_association_all_adult_wgs),union_all_species_wgs),"animalis"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_detection"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"detection"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_longum"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"longum"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_adolescentis"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"adolescentis"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_pseudocatenulatum"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"pseudocatenulatum"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_dentium"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"dentium"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_bifidum"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"bifidum"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_breve"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"breve"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_catenulatum"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"catenulatum"]

df_association_all_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"senior_animalis"] <- df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_wgs),union_all_species_wgs),"animalis"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_detection"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"detection"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_longum"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"longum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_adolescentis"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"adolescentis"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_pseudocatenulatum"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"pseudocatenulatum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_dentium"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"dentium"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_bifidum"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"bifidum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_breve"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"breve"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_catenulatum"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"catenulatum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"IndustrializedUrban_animalis"] <- df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_wgs_IndustrializedUrban),union_all_species_wgs),"animalis"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_detection"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"detection"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_longum"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"longum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_adolescentis"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"adolescentis"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_pseudocatenulatum"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"pseudocatenulatum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_dentium"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"dentium"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_bifidum"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"bifidum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_breve"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"breve"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_catenulatum"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"catenulatum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"UrbanRuralMixed_animalis"] <- df_association_adult_wgs_UrbanRuralMixed[intersect(rownames(df_association_adult_wgs_UrbanRuralMixed),union_all_species_wgs),"animalis"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_detection"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"detection"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_longum"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"longum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_adolescentis"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"adolescentis"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_pseudocatenulatum"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"pseudocatenulatum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_dentium"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"dentium"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_bifidum"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"bifidum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_breve"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"breve"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_catenulatum"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"catenulatum"]

df_association_all_wgs[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"RuralTribal_animalis"] <- df_association_adult_wgs_RuralTribal[intersect(rownames(df_association_adult_wgs_RuralTribal),union_all_species_wgs),"animalis"]

df_association_all_wgs <- as.data.frame(apply(df_association_all_wgs,2,function(x)(ifelse(is.na(x),0,x))))

save(df_association_all_adult_wgs,df_association_all_senior_wgs,df_association_all_infant_wgs,df_association_adult_wgs_RuralTribal,df_association_adult_wgs_UrbanRuralMixed,df_association_adult_wgs_IndustrializedUrban,df_association_all_wgs,file="C:\\Projects\\Bif_Manuscript\\AssociationScores.RData")

df_association_all_wgs <- as.data.frame(apply(df_association_all_wgs,2,function(x)(ifelse(is.na(x),0,x))))

mat_wgs <- apply(df_association_all_wgs,2,function(x)(ifelse(is.na(x),0,x)))
mat_wgs <- mat_wgs[,rev(colnames(mat_wgs))]
mat_wgs <- t(mat_wgs)
#heatmap.2(mat_wgs,density="none",trace="none",Rowv=FALSE,margins=c(12,10),lhei=c(0.5,5),lwid=c(1,5),col=brewer.pal(8,"PiYG"),srtRow=0,srtCol=90,sepcolor="grey",sepwidth=c(0.05,0.05),rowsep=c(0:nrow(mat_wgs)),colsep=c(0:ncol(mat_wgs)),cellnote=apply(mat_wgs,2,function(x)(ifelse(abs(x)>=0.20,"*",""))),notecol="black")

rank_scale=function(x)
{
	x <- rank(x);
	y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
	y <- ifelse(is.nan(y),0,y)
	return(y);
}

print("adult:16s")

df_association_adult_detection_16s$association$scores <- apply(df_association_adult_detection_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_detection_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_longum_16s$association$scores <- apply(df_association_adult_longum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_longum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_adolescentis_16s$association$scores <- apply(df_association_adult_adolescentis_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_adolescentis_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_pseudocatenulatum_16s$association$scores <- apply(df_association_adult_pseudocatenulatum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_pseudocatenulatum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_bifidum_16s$association$scores <- apply(df_association_adult_bifidum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_bifidum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_breve_16s$association$scores <- apply(df_association_adult_breve_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_breve_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_dentium_16s$association$scores <- apply(df_association_adult_dentium_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_dentium_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_catenulatum_16s$association$scores <- apply(df_association_adult_catenulatum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_catenulatum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_adult_animalis_16s$association$scores <- apply(df_association_adult_animalis_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_adult_animalis_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

union_species_adult_16s <- names(which(table(c(rownames(df_association_adult_animalis_16s$association),rownames(df_association_adult_catenulatum_16s$association),rownames(df_association_adult_breve_16s$association),rownames(df_association_adult_bifidum_16s$association),rownames(df_association_adult_dentium_16s$association),rownames(df_association_adult_pseudocatenulatum_16s$association),rownames(df_association_adult_adolescentis_16s$association),rownames(df_association_adult_longum_16s$association),rownames(df_association_adult_detection_16s$association)))==9))

df_association_all_adult_16s <- data.frame("detection"=df_association_adult_detection_16s$association[union_species_adult_16s,"scores"],"longum"=df_association_adult_longum_16s$association[union_species_adult_16s,"scores"],"adolescentis"=df_association_adult_adolescentis_16s$association[union_species_adult_16s,"scores"],"pseudocatenulatum"=df_association_adult_pseudocatenulatum_16s$association[union_species_adult_16s,"scores"],"dentium"=df_association_adult_dentium_16s$association[union_species_adult_16s,"scores"],"bifidum"=df_association_adult_bifidum_16s$association[union_species_adult_16s,"scores"],"breve"=df_association_adult_breve_16s$association[union_species_adult_16s,"scores"],"catenulatum"=df_association_adult_catenulatum_16s$association[union_species_adult_16s,"scores"],"animalis"=df_association_adult_animalis_16s$association[union_species_adult_16s,"scores"],row.names=union_species_adult_16s)

select_species_adult_16s <- names(which(apply(df_association_all_adult_16s,1,function(x)(length(x[abs(x)>=0.20])))>0))

df_association_all_adult_16s <- df_association_all_adult_16s[select_species_adult_16s,]


print("senior:16s")

df_association_senior_detection_16s$association$scores <- apply(df_association_senior_detection_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_detection_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_longum_16s$association$scores <- apply(df_association_senior_longum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_longum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_adolescentis_16s$association$scores <- apply(df_association_senior_adolescentis_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_adolescentis_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_pseudocatenulatum_16s$association$scores <- apply(df_association_senior_pseudocatenulatum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_pseudocatenulatum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_bifidum_16s$association$scores <- apply(df_association_senior_bifidum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_bifidum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_breve_16s$association$scores <- apply(df_association_senior_breve_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_breve_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_dentium_16s$association$scores <- apply(df_association_senior_dentium_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_dentium_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_catenulatum_16s$association$scores <- apply(df_association_senior_catenulatum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_catenulatum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_senior_animalis_16s$association$scores <- apply(df_association_senior_animalis_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_senior_animalis_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

union_species_senior_16s <- names(which(table(c(rownames(df_association_senior_animalis_16s$association),rownames(df_association_senior_catenulatum_16s$association),rownames(df_association_senior_breve_16s$association),rownames(df_association_senior_bifidum_16s$association),rownames(df_association_senior_dentium_16s$association),rownames(df_association_senior_pseudocatenulatum_16s$association),rownames(df_association_senior_adolescentis_16s$association),rownames(df_association_senior_longum_16s$association),rownames(df_association_senior_detection_16s$association)))==9))

df_association_all_senior_16s <- data.frame("detection"=df_association_senior_detection_16s$association[union_species_senior_16s,"scores"],"longum"=df_association_senior_longum_16s$association[union_species_senior_16s,"scores"],"adolescentis"=df_association_senior_adolescentis_16s$association[union_species_senior_16s,"scores"],"pseudocatenulatum"=df_association_senior_pseudocatenulatum_16s$association[union_species_senior_16s,"scores"],"dentium"=df_association_senior_dentium_16s$association[union_species_senior_16s,"scores"],"bifidum"=df_association_senior_bifidum_16s$association[union_species_senior_16s,"scores"],"breve"=df_association_senior_breve_16s$association[union_species_senior_16s,"scores"],"catenulatum"=df_association_senior_catenulatum_16s$association[union_species_senior_16s,"scores"],"animalis"=df_association_senior_animalis_16s$association[union_species_senior_16s,"scores"],row.names=union_species_senior_16s)

select_species_senior_16s <- names(which(apply(df_association_all_senior_16s,1,function(x)(length(x[abs(x)>=0.20])))>0))

df_association_all_senior_16s <- df_association_all_senior_16s[select_species_senior_16s,]

print("infant:16s")

df_association_infant_detection_16s$association$scores <- apply(df_association_infant_detection_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_detection_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_longum_16s$association$scores <- apply(df_association_infant_longum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_longum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_adolescentis_16s$association$scores <- apply(df_association_infant_adolescentis_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_adolescentis_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_pseudocatenulatum_16s$association$scores <- apply(df_association_infant_pseudocatenulatum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_pseudocatenulatum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_bifidum_16s$association$scores <- apply(df_association_infant_bifidum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_bifidum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_breve_16s$association$scores <- apply(df_association_infant_breve_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_breve_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_dentium_16s$association$scores <- apply(df_association_infant_dentium_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_dentium_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

df_association_infant_catenulatum_16s$association$scores <- apply(df_association_infant_catenulatum_16s$association,1,function(x)((x[1]-x[2])/(x[3])))*apply(df_association_infant_catenulatum_16s$association,1,function(x)((1-(min(x[1:2]+0.0001)/max(x[1:2]+0.0001)))))

union_species_infant_16s <- names(which(table(c(rownames(df_association_infant_catenulatum_16s$association),rownames(df_association_infant_breve_16s$association),rownames(df_association_infant_bifidum_16s$association),rownames(df_association_infant_dentium_16s$association),rownames(df_association_infant_pseudocatenulatum_16s$association),rownames(df_association_infant_adolescentis_16s$association),rownames(df_association_infant_longum_16s$association),rownames(df_association_infant_detection_16s$association)))==8))

df_association_all_infant_16s <- data.frame("detection"=df_association_infant_detection_16s$association[union_species_infant_16s,"scores"],"longum"=df_association_infant_longum_16s$association[union_species_infant_16s,"scores"],"adolescentis"=df_association_infant_adolescentis_16s$association[union_species_infant_16s,"scores"],"pseudocatenulatum"=df_association_infant_pseudocatenulatum_16s$association[union_species_infant_16s,"scores"],"dentium"=df_association_infant_dentium_16s$association[union_species_infant_16s,"scores"],"bifidum"=df_association_infant_bifidum_16s$association[union_species_infant_16s,"scores"],"breve"=df_association_infant_breve_16s$association[union_species_infant_16s,"scores"],"catenulatum"=df_association_infant_catenulatum_16s$association[union_species_infant_16s,"scores"],row.names=union_species_infant_16s)

select_species_infant_16s <- names(which(apply(df_association_all_infant_16s,1,function(x)(length(x[abs(x)>=0.20])))>0))

df_association_all_infant_16s <- df_association_all_infant_16s[select_species_infant_16s,]

df_corr_score <- as.data.frame(matrix(0,9,6))
rownames(df_corr_score) <- c("detection","longum","adolescentis","pseudocatenulatum","dentium","bifidum","breve","catenulatum","animalis")
colnames(df_corr_score) <- c("infant_16s","senior_16s","senior_16s","infant_wgs","senior_wgs","senior_wgs")

common_rows <- intersect(names(rowMeans(apply(df_association_infant_detection_16s$featureImportance,1,rank_scale))),names(apply(df_association_infant_detection_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["detection","infant_16s"] <- corr.test((rowMeans(apply(df_association_infant_detection_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_detection_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_longum_16s$featureImportance,1,rank_scale))),names(apply(df_association_infant_longum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["longum","infant_16s"] <- corr.test((rowMeans(apply(df_association_infant_longum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_longum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_adolescentis_16s$featureImportance,1,rank_scale))),names(apply(df_association_infant_adolescentis_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["adolescentis","infant_16s"] <- corr.test((rowMeans(apply(df_association_infant_adolescentis_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_adolescentis_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_pseudocatenulatum_16s$featureImportance,1,rank_scale))),names(apply(df_association_infant_pseudocatenulatum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["pseudocatenulatum","infant_16s"] <- corr.test((rowMeans(apply(df_association_infant_pseudocatenulatum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_pseudocatenulatum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_dentium_16s$featureImportance,1,rank_scale))),names(apply(df_association_infant_dentium_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["dentium","infant_16s"] <- corr.test((rowMeans(apply(df_association_infant_dentium_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_dentium_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_bifidum_16s$featureImportance,1,rank_scale))),names(apply(df_association_infant_bifidum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["bifidum","infant_16s"] <- corr.test((rowMeans(apply(df_association_infant_bifidum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_bifidum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_breve_16s$featureImportance,1,rank_scale))),names(apply(df_association_infant_breve_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["breve","infant_16s"] <- corr.test((rowMeans(apply(df_association_infant_breve_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_breve_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_catenulatum_16s$featureImportance,1,rank_scale))),names(apply(df_association_infant_catenulatum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["catenulatum","infant_16s"] <- corr.test((rowMeans(apply(df_association_infant_catenulatum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_catenulatum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_detection_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_detection_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["detection","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_detection_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_detection_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_longum_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_longum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["longum","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_longum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_longum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_adolescentis_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_adolescentis_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["adolescentis","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_adolescentis_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_adolescentis_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_pseudocatenulatum_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_pseudocatenulatum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["pseudocatenulatum","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_pseudocatenulatum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_pseudocatenulatum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_dentium_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_dentium_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["dentium","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_dentium_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_dentium_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_bifidum_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_bifidum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["bifidum","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_bifidum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_bifidum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_breve_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_breve_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["breve","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_breve_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_breve_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_catenulatum_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_catenulatum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["catenulatum","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_catenulatum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_catenulatum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_animalis_16s$featureImportance,1,rank_scale))),names(apply(df_association_adult_animalis_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["animalis","adult_16s"] <- corr.test((rowMeans(apply(df_association_adult_animalis_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_animalis_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_detection_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_detection_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["detection","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_detection_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_detection_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_longum_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_longum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["longum","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_longum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_longum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_adolescentis_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_adolescentis_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["adolescentis","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_adolescentis_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_adolescentis_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_pseudocatenulatum_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_pseudocatenulatum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["pseudocatenulatum","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_pseudocatenulatum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_pseudocatenulatum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_dentium_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_dentium_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["dentium","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_dentium_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_dentium_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_bifidum_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_bifidum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["bifidum","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_bifidum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_bifidum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_breve_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_breve_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["breve","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_breve_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_breve_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_catenulatum_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_catenulatum_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["catenulatum","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_catenulatum_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_catenulatum_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_animalis_16s$featureImportance,1,rank_scale))),names(apply(df_association_senior_animalis_16s$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["animalis","senior_16s"] <- corr.test((rowMeans(apply(df_association_senior_animalis_16s$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_animalis_16s$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_detection_wgs$featureImportance,1,rank_scale))),names(apply(df_association_infant_detection_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["detection","infant_wgs"] <- corr.test((rowMeans(apply(df_association_infant_detection_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_detection_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_longum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_infant_longum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["longum","infant_wgs"] <- corr.test((rowMeans(apply(df_association_infant_longum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_longum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_adolescentis_wgs$featureImportance,1,rank_scale))),names(apply(df_association_infant_adolescentis_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["adolescentis","infant_wgs"] <- corr.test((rowMeans(apply(df_association_infant_adolescentis_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_adolescentis_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_pseudocatenulatum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_infant_pseudocatenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["pseudocatenulatum","infant_wgs"] <- corr.test((rowMeans(apply(df_association_infant_pseudocatenulatum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_pseudocatenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_dentium_wgs$featureImportance,1,rank_scale))),names(apply(df_association_infant_dentium_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["dentium","infant_wgs"] <- corr.test((rowMeans(apply(df_association_infant_dentium_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_dentium_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_bifidum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_infant_bifidum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["bifidum","infant_wgs"] <- corr.test((rowMeans(apply(df_association_infant_bifidum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_bifidum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_breve_wgs$featureImportance,1,rank_scale))),names(apply(df_association_infant_breve_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["breve","infant_wgs"] <- corr.test((rowMeans(apply(df_association_infant_breve_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_breve_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_infant_catenulatum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_infant_catenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["catenulatum","infant_wgs"] <- corr.test((rowMeans(apply(df_association_infant_catenulatum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_infant_catenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_detection_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_detection_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["detection","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_detection_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_detection_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_longum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_longum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["longum","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_longum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_longum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_adolescentis_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_adolescentis_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["adolescentis","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_adolescentis_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_adolescentis_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_pseudocatenulatum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_pseudocatenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["pseudocatenulatum","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_pseudocatenulatum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_pseudocatenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_dentium_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_dentium_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["dentium","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_dentium_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_dentium_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_bifidum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_bifidum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["bifidum","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_bifidum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_bifidum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_breve_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_breve_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["breve","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_breve_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_breve_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_catenulatum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_catenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["catenulatum","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_catenulatum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_catenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_adult_animalis_wgs$featureImportance,1,rank_scale))),names(apply(df_association_adult_animalis_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["animalis","adult_wgs"] <- corr.test((rowMeans(apply(df_association_adult_animalis_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_adult_animalis_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_detection_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_detection_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["detection","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_detection_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_detection_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_longum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_longum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["longum","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_longum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_longum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_adolescentis_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_adolescentis_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["adolescentis","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_adolescentis_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_adolescentis_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_pseudocatenulatum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_pseudocatenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["pseudocatenulatum","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_pseudocatenulatum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_pseudocatenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_dentium_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_dentium_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["dentium","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_dentium_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_dentium_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_bifidum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_bifidum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["bifidum","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_bifidum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_bifidum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_breve_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_breve_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["breve","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_breve_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_breve_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_catenulatum_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_catenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["catenulatum","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_catenulatum_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_catenulatum_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

common_rows <- intersect(names(rowMeans(apply(df_association_senior_animalis_wgs$featureImportance,1,rank_scale))),names(apply(df_association_senior_animalis_wgs$association,1,function(x)(max(x[1:2])/x[3]))));df_corr_score["animalis","senior_wgs"] <- corr.test((rowMeans(apply(df_association_senior_animalis_wgs$featureImportance,1,rank_scale)))[common_rows],(apply(df_association_senior_animalis_wgs$association,1,function(x)(max(x[1:2])/x[3])))[common_rows])$r

df_detection_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_detection_16s$dir[intersect(rownames(df_association_adult_detection_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_detection_16s$dir[intersect(rownames(df_association_adult_detection_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_detection_16s$dir),IndustrializedUrban_cohort)))

df_detection_adult_16s_IndustrializedUrban$score <- apply(df_detection_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_detection_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_longum_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_longum_16s$dir[intersect(rownames(df_association_adult_longum_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_longum_16s$dir[intersect(rownames(df_association_adult_longum_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_longum_16s$dir),IndustrializedUrban_cohort)))

df_longum_adult_16s_IndustrializedUrban$score <- apply(df_longum_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_longum_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_adolescentis_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_adolescentis_16s$dir[intersect(rownames(df_association_adult_adolescentis_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_adolescentis_16s$dir[intersect(rownames(df_association_adult_adolescentis_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_adolescentis_16s$dir),IndustrializedUrban_cohort)))

df_adolescentis_adult_16s_IndustrializedUrban$score <- apply(df_adolescentis_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_adolescentis_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_pseudocatenulatum_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_pseudocatenulatum_16s$dir[intersect(rownames(df_association_adult_pseudocatenulatum_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_pseudocatenulatum_16s$dir[intersect(rownames(df_association_adult_pseudocatenulatum_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_pseudocatenulatum_16s$dir),IndustrializedUrban_cohort)))

df_pseudocatenulatum_adult_16s_IndustrializedUrban$score <- apply(df_pseudocatenulatum_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_pseudocatenulatum_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_dentium_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_dentium_16s$dir[intersect(rownames(df_association_adult_dentium_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_dentium_16s$dir[intersect(rownames(df_association_adult_dentium_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_dentium_16s$dir),IndustrializedUrban_cohort)))

df_dentium_adult_16s_IndustrializedUrban$score <- apply(df_dentium_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_dentium_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_bifidum_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_bifidum_16s$dir[intersect(rownames(df_association_adult_bifidum_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_bifidum_16s$dir[intersect(rownames(df_association_adult_bifidum_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_bifidum_16s$dir),IndustrializedUrban_cohort)))

df_bifidum_adult_16s_IndustrializedUrban$score <- apply(df_bifidum_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_bifidum_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_breve_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_breve_16s$dir[intersect(rownames(df_association_adult_breve_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_breve_16s$dir[intersect(rownames(df_association_adult_breve_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_breve_16s$dir),IndustrializedUrban_cohort)))

df_breve_adult_16s_IndustrializedUrban$score <- apply(df_breve_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_breve_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_catenulatum_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_catenulatum_16s$dir[intersect(rownames(df_association_adult_catenulatum_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_catenulatum_16s$dir[intersect(rownames(df_association_adult_catenulatum_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_catenulatum_16s$dir),IndustrializedUrban_cohort)))

df_catenulatum_adult_16s_IndustrializedUrban$score <- apply(df_catenulatum_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_catenulatum_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_animalis_adult_16s_IndustrializedUrban <- data.frame("SP"=apply(df_association_adult_animalis_16s$dir[intersect(rownames(df_association_adult_animalis_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_animalis_16s$dir[intersect(rownames(df_association_adult_animalis_16s$dir),IndustrializedUrban_cohort),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_animalis_16s$dir),IndustrializedUrban_cohort)))

df_animalis_adult_16s_IndustrializedUrban$score <- apply(df_animalis_adult_16s_IndustrializedUrban,1,function(x)((x[1]-x[2])/x[3])) * apply(df_animalis_adult_16s_IndustrializedUrban,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_detection_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_detection_16s$dir[intersect(rownames(df_association_adult_detection_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_detection_16s$dir[intersect(rownames(df_association_adult_detection_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_detection_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_detection_adult_16s_Others$score <- apply(df_detection_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_detection_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_longum_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_longum_16s$dir[intersect(rownames(df_association_adult_longum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_longum_16s$dir[intersect(rownames(df_association_adult_longum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_longum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_longum_adult_16s_Others$score <- apply(df_longum_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_longum_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_adolescentis_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_adolescentis_16s$dir[intersect(rownames(df_association_adult_adolescentis_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_adolescentis_16s$dir[intersect(rownames(df_association_adult_adolescentis_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_adolescentis_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_adolescentis_adult_16s_Others$score <- apply(df_adolescentis_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_adolescentis_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_pseudocatenulatum_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_pseudocatenulatum_16s$dir[intersect(rownames(df_association_adult_pseudocatenulatum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_pseudocatenulatum_16s$dir[intersect(rownames(df_association_adult_pseudocatenulatum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_pseudocatenulatum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_pseudocatenulatum_adult_16s_Others$score <- apply(df_pseudocatenulatum_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_pseudocatenulatum_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_dentium_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_dentium_16s$dir[intersect(rownames(df_association_adult_dentium_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_dentium_16s$dir[intersect(rownames(df_association_adult_dentium_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_dentium_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_dentium_adult_16s_Others$score <- apply(df_dentium_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_dentium_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_bifidum_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_bifidum_16s$dir[intersect(rownames(df_association_adult_bifidum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_bifidum_16s$dir[intersect(rownames(df_association_adult_bifidum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_bifidum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_bifidum_adult_16s_Others$score <- apply(df_bifidum_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_bifidum_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_breve_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_breve_16s$dir[intersect(rownames(df_association_adult_breve_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_breve_16s$dir[intersect(rownames(df_association_adult_breve_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_breve_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_breve_adult_16s_Others$score <- apply(df_breve_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_breve_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_catenulatum_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_catenulatum_16s$dir[intersect(rownames(df_association_adult_catenulatum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_catenulatum_16s$dir[intersect(rownames(df_association_adult_catenulatum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_catenulatum_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_catenulatum_adult_16s_Others$score <- apply(df_catenulatum_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_catenulatum_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

df_animalis_adult_16s_Others <- data.frame("SP"=apply(df_association_adult_animalis_16s$dir[intersect(rownames(df_association_adult_animalis_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(1)]))),"SN"=apply(df_association_adult_animalis_16s$dir[intersect(rownames(df_association_adult_animalis_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort)),],2,function(x)(length(x[x==(-1)]))),"T"=length(intersect(rownames(df_association_adult_animalis_16s$dir),c(RuralTribal_cohort,UrbanRuralMixed_cohort))))

df_animalis_adult_16s_Others$score <- apply(df_animalis_adult_16s_Others,1,function(x)((x[1]-x[2])/x[3])) * apply(df_animalis_adult_16s_Others,1,function(x)(1-(min(x[1:2])+0.00001)/(max(x[1:2])+0.00001)))

union_species_adult_16s_IndustrializedUrban <- names(table(c(rownames(df_animalis_adult_16s_IndustrializedUrban),rownames(df_catenulatum_adult_16s_IndustrializedUrban),rownames(df_breve_adult_16s_IndustrializedUrban),rownames(df_bifidum_adult_16s_IndustrializedUrban),rownames(df_dentium_adult_16s_IndustrializedUrban),rownames(df_pseudocatenulatum_adult_16s_IndustrializedUrban),rownames(df_adolescentis_adult_16s_IndustrializedUrban),rownames(df_longum_adult_16s_IndustrializedUrban),rownames(df_detection_adult_16s_IndustrializedUrban))))

df_association_adult_16s_IndustrializedUrban <- data.frame("detection"=df_detection_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"longum"=df_longum_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"adolescentis"=df_adolescentis_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"pseudocatenulatum"=df_pseudocatenulatum_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"dentium"=df_dentium_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"bifidum"=df_bifidum_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"breve"=df_breve_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"catenulatum"=df_catenulatum_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"animalis"=df_animalis_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],row.names=union_species_adult_16s_IndustrializedUrban)

select_adult_16s_IndustrializedUrban <- names(which(apply(df_association_adult_16s_IndustrializedUrban,1,function(x)(length(x[abs(x)>=0.20])))>1))

df_association_adult_16s_IndustrializedUrban <- df_association_adult_16s_IndustrializedUrban[select_adult_16s_IndustrializedUrban,]

union_species_adult_16s_IndustrializedUrban <- names(table(c(rownames(df_animalis_adult_16s_IndustrializedUrban),rownames(df_catenulatum_adult_16s_IndustrializedUrban),rownames(df_breve_adult_16s_IndustrializedUrban),rownames(df_bifidum_adult_16s_IndustrializedUrban),rownames(df_dentium_adult_16s_IndustrializedUrban),rownames(df_pseudocatenulatum_adult_16s_IndustrializedUrban),rownames(df_adolescentis_adult_16s_IndustrializedUrban),rownames(df_longum_adult_16s_IndustrializedUrban),rownames(df_detection_adult_16s_IndustrializedUrban))))

df_association_adult_16s_IndustrializedUrban <- data.frame("detection"=df_detection_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"longum"=df_longum_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"adolescentis"=df_adolescentis_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"pseudocatenulatum"=df_pseudocatenulatum_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"dentium"=df_dentium_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"bifidum"=df_bifidum_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"breve"=df_breve_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"catenulatum"=df_catenulatum_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],"animalis"=df_animalis_adult_16s_IndustrializedUrban[union_species_adult_16s_IndustrializedUrban,"score"],row.names=union_species_adult_16s_IndustrializedUrban)

select_adult_16s_IndustrializedUrban <- names(which(apply(df_association_adult_16s_IndustrializedUrban,1,function(x)(length(x[abs(x)>=0.30])))>1))

df_association_adult_16s_IndustrializedUrban <- df_association_adult_16s_IndustrializedUrban[select_adult_16s_IndustrializedUrban,]

union_species_adult_16s_Others <- names(table(c(rownames(df_animalis_adult_16s_Others),rownames(df_catenulatum_adult_16s_Others),rownames(df_breve_adult_16s_Others),rownames(df_bifidum_adult_16s_Others),rownames(df_dentium_adult_16s_Others),rownames(df_pseudocatenulatum_adult_16s_Others),rownames(df_adolescentis_adult_16s_Others),rownames(df_longum_adult_16s_Others),rownames(df_detection_adult_16s_Others))))

df_association_adult_16s_Others <- data.frame("detection"=df_detection_adult_16s_Others[union_species_adult_16s_Others,"score"],"longum"=df_longum_adult_16s_Others[union_species_adult_16s_Others,"score"],"adolescentis"=df_adolescentis_adult_16s_Others[union_species_adult_16s_Others,"score"],"pseudocatenulatum"=df_pseudocatenulatum_adult_16s_Others[union_species_adult_16s_Others,"score"],"dentium"=df_dentium_adult_16s_Others[union_species_adult_16s_Others,"score"],"bifidum"=df_bifidum_adult_16s_Others[union_species_adult_16s_Others,"score"],"breve"=df_breve_adult_16s_Others[union_species_adult_16s_Others,"score"],"catenulatum"=df_catenulatum_adult_16s_Others[union_species_adult_16s_Others,"score"],"animalis"=df_animalis_adult_16s_Others[union_species_adult_16s_Others,"score"],row.names=union_species_adult_16s_Others)

select_adult_16s_Others <- names(which(apply(df_association_adult_16s_Others[,1:5],1,function(x)(length(x[abs(x)>=0.30])))>1))

df_association_adult_16s_Others <- df_association_adult_16s_Others[select_adult_16s_Others,1:5]

union_all_species_16s <- names(which(table(c(rownames(df_association_all_senior_16s),rownames(df_association_all_infant_16s),rownames(df_association_all_senior_16s),rownames(df_association_adult_16s_IndustrializedUrban),rownames(df_association_adult_16s_Others)))>1))

df_association_all_16s <- as.data.frame(matrix(0,length(union_all_species_16s),40))
rownames(df_association_all_16s) <- union_all_species_16s
colnames(df_association_all_16s) <- c("infant_detection","infant_longum","infant_adolescentis","infant_pseudocatenulatum","infant_dentium","infant_bifidum","infant_breve","infant_catenulatum","adult_detection","adult_longum","adult_adolescentis","adult_pseudocatenulatum","adult_dentium","adult_bifidum","adult_breve","adult_catenulatum","adult_animalis","senior_detection","senior_longum","senior_adolescentis","senior_pseudocatenulatum","senior_dentium","senior_bifidum","senior_breve","senior_catenulatum","senior_animalis","IndustrializedUrban_detection","IndustrializedUrban_longum","IndustrializedUrban_adolescentis","IndustrializedUrban_pseudocatenulatum","IndustrializedUrban_dentium","IndustrializedUrban_bifidum","IndustrializedUrban_breve","IndustrializedUrban_catenulatum","IndustrializedUrban_animalis","Others_detection","Others_longum","Others_adolescentis","Others_pseudocatenulatum","Others_dentium")

df_association_all_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"infant_detection"] <- df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"detection"]

df_association_all_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"infant_longum"] <- df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"longum"]

df_association_all_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"infant_adolescentis"] <- df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"adolescentis"]

df_association_all_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"infant_pseudocatenulatum"] <- df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"pseudocatenulatum"]

df_association_all_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"infant_dentium"] <- df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"dentium"]

df_association_all_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"infant_bifidum"] <- df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"bifidum"]

df_association_all_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"infant_breve"] <- df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"breve"]

df_association_all_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"infant_catenulatum"] <- df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),union_all_species_16s),"catenulatum"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_detection"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"detection"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_longum"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"longum"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_adolescentis"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adolescentis"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_pseudocatenulatum"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"pseudocatenulatum"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_dentium"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"dentium"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_bifidum"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"bifidum"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_breve"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"breve"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_catenulatum"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"catenulatum"]

df_association_all_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"adult_animalis"] <- df_association_all_adult_16s[intersect(rownames(df_association_all_adult_16s),union_all_species_16s),"animalis"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_detection"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"detection"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_longum"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"longum"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_adolescentis"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"adolescentis"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_pseudocatenulatum"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"pseudocatenulatum"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_dentium"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"dentium"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_bifidum"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"bifidum"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_breve"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"breve"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_catenulatum"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"catenulatum"]

df_association_all_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"senior_animalis"] <- df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),union_all_species_16s),"animalis"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_detection"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"detection"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_longum"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"longum"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_adolescentis"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"adolescentis"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_pseudocatenulatum"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"pseudocatenulatum"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_dentium"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"dentium"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_bifidum"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"bifidum"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_breve"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"breve"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_catenulatum"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"catenulatum"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"IndustrializedUrban_animalis"] <- df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),union_all_species_16s),"animalis"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"Others_detection"] <- df_association_adult_16s_Others[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"detection"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"Others_longum"] <- df_association_adult_16s_Others[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"longum"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"Others_adolescentis"] <- df_association_adult_16s_Others[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"adolescentis"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"Others_pseudocatenulatum"] <- df_association_adult_16s_Others[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"pseudocatenulatum"]

df_association_all_16s[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"Others_dentium"] <- df_association_adult_16s_Others[intersect(rownames(df_association_adult_16s_Others),union_all_species_16s),"dentium"]

df_association_all_16s <- as.data.frame(apply(df_association_all_16s,2,function(x)(ifelse(is.na(x),0,x))))

save(df_association_all_adult_16s,df_association_all_senior_16s,df_association_all_infant_16s,df_association_adult_16s_Others,df_association_adult_16s_IndustrializedUrban,df_association_all_16s,file="C:\\Projects\\Bif_Manuscript\\AssociationScores_16s.RData")

mat_16s <- apply(df_association_all_16s,2,function(x)(ifelse(is.na(x),0,x)))
mat_16s <- t(mat_16s)
mat_16s <- mat_16s[,rev(colnames(mat_16s))]
#heatmap.2(mat_16s,density="none",trace="none",Rowv=FALSE,margins=c(12,10),lhei=c(0.5,5),lwid=c(1,5),col=brewer.pal(8,"PiYG"),srtRow=0,srtCol=90,sepcolor="grey",sepwidth=c(0.05,0.05),rowsep=c(0:nrow(mat_16s)),colsep=c(0:ncol(mat_16s)),cellnote=apply(mat_16s,2,function(x)(ifelse(abs(x)>=0.30,"*",""))),notecol="black")

mean_scores_wgs <- data.frame("mean_score_detection"=apply(df_association_all_wgs[,c("infant_detection","senior_detection","IndustrializedUrban_detection","UrbanRuralMixed_detection","RuralTribal_detection")],1,function(x)(mean(x[abs(x)>0]))),"mean_score_longum"=apply(df_association_all_wgs[,c("infant_longum","senior_longum","IndustrializedUrban_longum","UrbanRuralMixed_longum","RuralTribal_longum")],1,function(x)(mean(x[abs(x)>0]))),"mean_score_adolescentis"=apply(df_association_all_wgs[,c("infant_adolescentis","senior_adolescentis","IndustrializedUrban_adolescentis","UrbanRuralMixed_adolescentis","RuralTribal_adolescentis")],1,function(x)(mean(x[abs(x)>0]))),"mean_score_pseudocatenulatum"=apply(df_association_all_wgs[,c("infant_pseudocatenulatum","senior_pseudocatenulatum","IndustrializedUrban_pseudocatenulatum","UrbanRuralMixed_pseudocatenulatum","RuralTribal_pseudocatenulatum")],1,function(x)(mean(x[abs(x)>0]))),"mean_score_dentium"=apply(df_association_all_wgs[,c("infant_dentium","senior_dentium","IndustrializedUrban_dentium","UrbanRuralMixed_dentium","RuralTribal_dentium")],1,function(x)(mean(x[abs(x)>0]))),"mean_score_bifidum"=apply(df_association_all_wgs[,c("infant_bifidum","senior_bifidum","IndustrializedUrban_bifidum","UrbanRuralMixed_bifidum","RuralTribal_bifidum")],1,function(x)(mean(x[abs(x)>0]))),"mean_score_breve"=apply(df_association_all_wgs[,c("infant_breve","senior_breve","IndustrializedUrban_breve","UrbanRuralMixed_breve","RuralTribal_breve")],1,function(x)(mean(x[abs(x)>0]))),"mean_score_catenulatum"=apply(df_association_all_wgs[,c("infant_catenulatum","senior_catenulatum","IndustrializedUrban_catenulatum","UrbanRuralMixed_catenulatum","RuralTribal_catenulatum")],1,function(x)(mean(x[abs(x)>0]))),"mean_score_animalis"=apply(df_association_all_wgs[,c("senior_animalis","IndustrializedUrban_animalis","UrbanRuralMixed_animalis","RuralTribal_animalis")],1,function(x)(mean(x[abs(x)>0]))))

df_combined_scores_detection_wgs_16s_infants <- data.frame("16s"=df_association_all_infant_16s[intersect(rownames(df_association_all_infant_16s),rownames(df_association_all_infant_wgs)),"detection"],"wgs"=df_association_all_infant_wgs[intersect(rownames(df_association_all_infant_16s),rownames(df_association_all_infant_wgs)),"detection"],row.names=intersect(rownames(df_association_all_infant_16s),rownames(df_association_all_infant_wgs)))

ggplot(df_combined_scores_detection_wgs_16s_infants,aes(x=X16s,y=wgs))+geom_point(col=ifelse((df_combined_scores_detection_wgs_16s_infants$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_infants$wgs>=0.20),"blue",ifelse((df_combined_scores_detection_wgs_16s_infants$X16s<=-0.20)&(df_combined_scores_detection_wgs_16s_infants<=-0.20),"red","grey")))+geom_smooth(method='lm',col="grey50",alpha=0.25)+geom_text_repel(label=ifelse((df_combined_scores_detection_wgs_16s_infants$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_infants$wgs>=0.20),rownames(df_combined_scores_detection_wgs_16s_infants),ifelse((df_combined_scores_detection_wgs_16s_infants$X16s<=-0.20)&(df_combined_scores_detection_wgs_16s_infants$wgs<=-0.20),rownames(df_combined_scores_detection_wgs_16s_infants),"")),col=ifelse((df_combined_scores_detection_wgs_16s_infants$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_infants$wgs>=0.20),"blue4",ifelse((df_combined_scores_detection_wgs_16s_infants$X16s<=-0.20)&(df_combined_scores_detection_wgs_16s_infants$wgs<=-0.20),"firebrick4","")),max.overlaps=30,size=4,box.padding=0.0,point.padding=1e-25,force_pull=10,min.segment.length=4)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_hline(yintercept=0.20,color="blue")+geom_vline(xintercept=0.20,color="blue")+geom_hline(yintercept=-0.20,color="red")+geom_vline(xintercept=-0.20,color="red")

df_combined_scores_detection_wgs_16s_IndustrializedUrban <- data.frame("16s"=df_association_adult_16s_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),rownames(df_association_adult_wgs_IndustrializedUrban)),"detection"],"wgs"=df_association_adult_wgs_IndustrializedUrban[intersect(rownames(df_association_adult_16s_IndustrializedUrban),rownames(df_association_adult_wgs_IndustrializedUrban)),"detection"],row.names=intersect(rownames(df_association_adult_16s_IndustrializedUrban),rownames(df_association_adult_wgs_IndustrializedUrban)))

ggplot(df_combined_scores_detection_wgs_16s_IndustrializedUrban,aes(x=X16s,y=wgs))+geom_point(col=ifelse((df_combined_scores_detection_wgs_16s_IndustrializedUrban$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_IndustrializedUrban$wgs>=0.20),"blue",ifelse((df_combined_scores_detection_wgs_16s_IndustrializedUrban$X16s<=-0.20)&(df_combined_scores_detection_wgs_16s_IndustrializedUrban<=-0.20),"red","grey")))+geom_smooth(method='lm',col="grey50",alpha=0.25)+geom_text_repel(label=ifelse((df_combined_scores_detection_wgs_16s_IndustrializedUrban$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_IndustrializedUrban$wgs>=0.20),rownames(df_combined_scores_detection_wgs_16s_IndustrializedUrban),ifelse((df_combined_scores_detection_wgs_16s_IndustrializedUrban$X16s<=-0.20)|(df_combined_scores_detection_wgs_16s_IndustrializedUrban$wgs<=-0.20),rownames(df_combined_scores_detection_wgs_16s_IndustrializedUrban),"")),col=ifelse((df_combined_scores_detection_wgs_16s_IndustrializedUrban$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_IndustrializedUrban$wgs>=0.20),"blue4",ifelse((df_combined_scores_detection_wgs_16s_IndustrializedUrban$X16s<=-0.20)|(df_combined_scores_detection_wgs_16s_IndustrializedUrban$wgs<=-0.20),"firebrick4","")),max.overlaps=30,size=4,box.padding=0.0,point.padding=1e-25,force_pull=10,min.segment.length=4)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_hline(yintercept=0.20,color="blue")+geom_vline(xintercept=0.20,color="blue")+geom_hline(yintercept=-0.20,color="red")+geom_vline(xintercept=-0.20,color="red")

df_combined_scores_detection_wgs_16s_senior <- data.frame("16s"=df_association_all_senior_16s[intersect(rownames(df_association_all_senior_16s),rownames(df_association_all_senior_wgs)),"detection"],"wgs"=df_association_all_senior_wgs[intersect(rownames(df_association_all_senior_16s),rownames(df_association_all_senior_wgs)),"detection"],row.names=intersect(rownames(df_association_all_senior_16s),rownames(df_association_all_senior_wgs)))

ggplot(df_combined_scores_detection_wgs_16s_senior,aes(x=X16s,y=wgs))+geom_point(col=ifelse((df_combined_scores_detection_wgs_16s_senior$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_senior$wgs>=0.20),"blue",ifelse((df_combined_scores_detection_wgs_16s_senior$X16s<=-0.20)&(df_combined_scores_detection_wgs_16s_senior<=-0.20),"red","grey")))+geom_smooth(method='lm',col="grey50",alpha=0.25)+geom_text_repel(label=ifelse((df_combined_scores_detection_wgs_16s_senior$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_senior$wgs>=0.20),rownames(df_combined_scores_detection_wgs_16s_senior),ifelse((df_combined_scores_detection_wgs_16s_senior$X16s<=-0.20)|(df_combined_scores_detection_wgs_16s_senior$wgs<=-0.20),rownames(df_combined_scores_detection_wgs_16s_senior),"")),col=ifelse((df_combined_scores_detection_wgs_16s_senior$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_senior$wgs>=0.20),"blue4",ifelse((df_combined_scores_detection_wgs_16s_senior$X16s<=-0.20)|(df_combined_scores_detection_wgs_16s_senior$wgs<=-0.20),"firebrick4","")),max.overlaps=30,size=4,box.padding=0.0,point.padding=1e-25,force_pull=10,min.segment.length=4)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_hline(yintercept=0.20,color="blue")+geom_vline(xintercept=0.20,color="blue")+geom_hline(yintercept=-0.20,color="red")+geom_vline(xintercept=-0.20,color="red")

common_rows <- intersect(rownames(df_association_adult_16s_Others),intersect(rownames(df_association_adult_wgs_RuralTribal),rownames(df_association_adult_wgs_UrbanRuralMixed)))

df_combined_scores_detection_wgs_16s_Others <- data.frame("16s"=df_association_adult_16s_Others[common_rows,"detection"],"wgs"=ifelse((df_association_adult_wgs_RuralTribal[common_rows,"detection"]<0)&(df_association_adult_wgs_UrbanRuralMixed[common_rows,"detection"]<0),min(df_association_adult_wgs_RuralTribal[common_rows,"detection"],df_association_adult_wgs_UrbanRuralMixed[common_rows,"detection"]),max(df_association_adult_wgs_RuralTribal[common_rows,"detection"],df_association_adult_wgs_UrbanRuralMixed[common_rows,"detection"])),row.names=common_rows)

ggplot(df_combined_scores_detection_wgs_16s_Others,aes(x=X16s,y=wgs))+geom_point(col=ifelse((df_combined_scores_detection_wgs_16s_Others$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_Others$wgs>=0.20),"blue",ifelse((df_combined_scores_detection_wgs_16s_Others$X16s<=-0.20)&(df_combined_scores_detection_wgs_16s_Others<=-0.20),"red","grey")))+geom_smooth(method='lm',col="grey50",alpha=0.25)+geom_text_repel(label=ifelse((df_combined_scores_detection_wgs_16s_Others$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_Others$wgs>=0.20),rownames(df_combined_scores_detection_wgs_16s_Others),ifelse((df_combined_scores_detection_wgs_16s_Others$X16s<=-0.20)|(df_combined_scores_detection_wgs_16s_Others$wgs<=-0.20),rownames(df_combined_scores_detection_wgs_16s_Others),"")),col=ifelse((df_combined_scores_detection_wgs_16s_Others$X16s>=0.20)&(df_combined_scores_detection_wgs_16s_Others$wgs>=0.20),"blue4",ifelse((df_combined_scores_detection_wgs_16s_Others$X16s<=-0.20)|(df_combined_scores_detection_wgs_16s_Others$wgs<=-0.20),"firebrick4","")),max.overlaps=60,size=4,box.padding=0.0,point.padding=1e-25,force_pull=10,min.segment.length=4)+theme_bw()+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))+geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_hline(yintercept=0.20,color="blue")+geom_vline(xintercept=0.20,color="blue")+geom_hline(yintercept=-0.20,color="red")+geom_vline(xintercept=-0.20,color="red")+ylim(-1,0.6)

df_corr_16s_wgs <- as.data.frame(matrix(0,5,9))
rownames(df_corr_16s_wgs) <- c("infant","adult","senior","IndustrializedUrban","Other")
colnames(df_corr_16s_wgs) <- c("detection","longum","adolescentis","pseudocatenulatum","bifidum","breve","catenulatum","animalis","dentium")

df_p_16s_wgs <- as.data.frame(matrix(1,5,9))
rownames(df_p_16s_wgs) <- c("infant","adult","senior","IndustrializedUrban","Other")
colnames(df_p_16s_wgs) <- c("detection","longum","adolescentis","pseudocatenulatum","bifidum","breve","catenulatum","animalis","dentium")

df_corr_16s_wgs["infant","detection"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_detection"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_detection"])$r
df_p_16s_wgs["infant","detection"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_detection"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_detection"])$p

df_corr_16s_wgs["infant","longum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_longum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_longum"])$r
df_p_16s_wgs["infant","longum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_longum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_longum"])$p

df_corr_16s_wgs["infant","adolescentis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_adolescentis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_adolescentis"])$r
df_p_16s_wgs["infant","adolescentis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_adolescentis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_adolescentis"])$p

df_corr_16s_wgs["infant","pseudocatenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_pseudocatenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_pseudocatenulatum"])$r
df_p_16s_wgs["infant","pseudocatenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_pseudocatenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_pseudocatenulatum"])$p

df_corr_16s_wgs["infant","bifidum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_bifidum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_bifidum"])$r
df_p_16s_wgs["infant","bifidum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_bifidum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_bifidum"])$p

df_corr_16s_wgs["infant","breve"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_breve"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_breve"])$r
df_p_16s_wgs["infant","breve"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_breve"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_breve"])$p

df_corr_16s_wgs["infant","dentium"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_dentium"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_dentium"])$r
df_p_16s_wgs["infant","dentium"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_dentium"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_dentium"])$p

df_corr_16s_wgs["infant","catenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_catenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_catenulatum"])$r
df_p_16s_wgs["infant","catenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_catenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"infant_catenulatum"])$p

df_corr_16s_wgs["adult","detection"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_detection"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_detection"])$r
df_p_16s_wgs["adult","detection"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_detection"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_detection"])$p

df_corr_16s_wgs["adult","longum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_longum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_longum"])$r
df_p_16s_wgs["adult","longum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_longum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_longum"])$p

df_corr_16s_wgs["adult","adolescentis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_adolescentis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_adolescentis"])$r
df_p_16s_wgs["adult","adolescentis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_adolescentis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_adolescentis"])$p

df_corr_16s_wgs["adult","pseudocatenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_pseudocatenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_pseudocatenulatum"])$r
df_p_16s_wgs["adult","pseudocatenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_pseudocatenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_pseudocatenulatum"])$p

df_corr_16s_wgs["adult","bifidum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_bifidum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_bifidum"])$r
df_p_16s_wgs["adult","bifidum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_bifidum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_bifidum"])$p

df_corr_16s_wgs["adult","breve"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_breve"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_breve"])$r
df_p_16s_wgs["adult","breve"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_breve"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_breve"])$p

df_corr_16s_wgs["adult","dentium"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_dentium"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_dentium"])$r
df_p_16s_wgs["adult","dentium"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_dentium"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_dentium"])$p

df_corr_16s_wgs["adult","catenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_catenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_catenulatum"])$r
df_p_16s_wgs["adult","catenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_catenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_catenulatum"])$p

df_corr_16s_wgs["adult","animalis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_animalis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_animalis"])$r
df_p_16s_wgs["adult","animalis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_animalis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"adult_animalis"])$p

df_corr_16s_wgs["senior","detection"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_detection"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_detection"])$r
df_p_16s_wgs["senior","detection"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_detection"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_detection"])$p

df_corr_16s_wgs["senior","longum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_longum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_longum"])$r
df_p_16s_wgs["senior","longum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_longum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_longum"])$p

df_corr_16s_wgs["senior","adolescentis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_adolescentis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_adolescentis"])$r
df_p_16s_wgs["senior","adolescentis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_adolescentis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_adolescentis"])$p

df_corr_16s_wgs["senior","pseudocatenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_pseudocatenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_pseudocatenulatum"])$r
df_p_16s_wgs["senior","pseudocatenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_pseudocatenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_pseudocatenulatum"])$p

df_corr_16s_wgs["senior","bifidum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_bifidum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_bifidum"])$r
df_p_16s_wgs["senior","bifidum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_bifidum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_bifidum"])$p

df_corr_16s_wgs["senior","breve"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_breve"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_breve"])$r
df_p_16s_wgs["senior","breve"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_breve"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_breve"])$p

df_corr_16s_wgs["senior","dentium"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_dentium"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_dentium"])$r
df_p_16s_wgs["senior","dentium"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_dentium"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_dentium"])$p

df_corr_16s_wgs["senior","catenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_catenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_catenulatum"])$r
df_p_16s_wgs["senior","catenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_catenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_catenulatum"])$p

df_corr_16s_wgs["senior","animalis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_animalis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_animalis"])$r
df_p_16s_wgs["senior","animalis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_animalis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"senior_animalis"])$p

df_corr_16s_wgs["IndustrializedUrban","detection"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_detection"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_detection"])$r
df_p_16s_wgs["IndustrializedUrban","detection"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_detection"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_detection"])$p

df_corr_16s_wgs["IndustrializedUrban","longum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_longum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_longum"])$r
df_p_16s_wgs["IndustrializedUrban","longum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_longum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_longum"])$p

df_corr_16s_wgs["IndustrializedUrban","adolescentis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_adolescentis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_adolescentis"])$r
df_p_16s_wgs["IndustrializedUrban","adolescentis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_adolescentis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_adolescentis"])$p

df_corr_16s_wgs["IndustrializedUrban","pseudocatenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_pseudocatenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_pseudocatenulatum"])$r
df_p_16s_wgs["IndustrializedUrban","pseudocatenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_pseudocatenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_pseudocatenulatum"])$p

df_corr_16s_wgs["IndustrializedUrban","bifidum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_bifidum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_bifidum"])$r
df_p_16s_wgs["IndustrializedUrban","bifidum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_bifidum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_bifidum"])$p

df_corr_16s_wgs["IndustrializedUrban","breve"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_breve"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_breve"])$r
df_p_16s_wgs["IndustrializedUrban","breve"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_breve"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_breve"])$p

df_corr_16s_wgs["IndustrializedUrban","dentium"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_dentium"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_dentium"])$r
df_p_16s_wgs["IndustrializedUrban","dentium"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_dentium"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_dentium"])$p

df_corr_16s_wgs["IndustrializedUrban","catenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_catenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_catenulatum"])$r
df_p_16s_wgs["IndustrializedUrban","catenulatum"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_catenulatum"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_catenulatum"])$p

df_corr_16s_wgs["IndustrializedUrban","animalis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_animalis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_animalis"])$r
df_p_16s_wgs["IndustrializedUrban","animalis"] <- corr.test(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_animalis"],df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"IndustrializedUrban_animalis"])$p

df_corr_16s_wgs["Other","detection"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_detection","RuralTribal_detection")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_detection"])$r
df_p_16s_wgs["IndustrializedUrban","detection"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_detection","RuralTribal_detection")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_detection"])$p

df_corr_16s_wgs["Other","detection"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_detection","RuralTribal_detection")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_detection"])$r
df_p_16s_wgs["Other","detection"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_detection","RuralTribal_detection")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_detection"])$p

df_corr_16s_wgs["Other","longum"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_longum","RuralTribal_longum")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_longum"])$r
df_p_16s_wgs["Other","longum"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_longum","RuralTribal_longum")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_longum"])$p

df_corr_16s_wgs["Other","adolescentis"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_adolescentis","RuralTribal_adolescentis")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_adolescentis"])$r
df_p_16s_wgs["Other","adolescentis"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_adolescentis","RuralTribal_adolescentis")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_adolescentis"])$p

df_corr_16s_wgs["Other","pseudocatenulatum"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_pseudocatenulatum","RuralTribal_pseudocatenulatum")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_pseudocatenulatum"])$r
df_p_16s_wgs["Other","pseudocatenulatum"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_pseudocatenulatum","RuralTribal_pseudocatenulatum")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_pseudocatenulatum"])$p

df_corr_16s_wgs["Other","dentium"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_dentium","RuralTribal_dentium")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_dentium"])$r
df_p_16s_wgs["Other","dentium"] <- corr.test(rowMeans(df_association_all_wgs[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),c("UrbanRuralMixed_dentium","RuralTribal_dentium")]),df_association_all_16s[intersect(rownames(df_association_all_wgs),rownames(df_association_all_16s)),"Others_dentium"])$p

save.image("C:\\Projects\\Bif_Manuscript\\bif_ml_analysis_scores.RData")
