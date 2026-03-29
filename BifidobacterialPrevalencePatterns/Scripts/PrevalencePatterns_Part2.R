load("/data/bifido_45809_new_2014_not_normalized_Harnandez_added.RData")
SpeciesProfile_All[is.na(SpeciesProfile_All)] <- 0

names(SpeciesProfile_All)[1265] <- "Intestinibacter_bartlettii"
names(SpeciesProfile_All)[1362] <- "Erysipelatoclostridium_ramosum"

spdata_45809 <- SpeciesProfile_All
metadata_45809 <- metadata_updated

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

load("/data/SpeciesDetection_AgeCategories.RData")

infant_select_species <- names(which(apply(infant_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(infant_species_detection)>=0.33))
adult_select_species <- names(which(apply(adult_species_detection,1,function(x)(length(x[x>=0.05])))/ncol(adult_species_detection)>=0.33))
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

load("/data/df_infant_bifido_new.RData")
load("/data/df_adult_bifido_new.RData")
load("/data/df_senior_bifido_new.RData")


df_infant_bifido_new[,"AgeCategory"] <- "Infant"
df_adult_bifido_new[,"AgeCategory"] <- "Adult"
df_senior_bifido_new[,"AgeCategory"] <- "Senior"

rownames(df_infant_bifido_new) <- paste0("infant:",rownames(df_infant_bifido_new))
rownames(df_adult_bifido_new) <- paste0("adult:",rownames(df_adult_bifido_new))
rownames(df_senior_bifido_new) <- paste0("senior:",rownames(df_senior_bifido_new))

df_combined_bifido <- as.data.frame(rbind(df_infant_bifido_new,df_adult_bifido_new,df_senior_bifido_new))


library(dplyr)
library(ade4)
library(vegan)
library(adegraphics)
library(gplots)
library(RColorBrewer)

MajorBifsOrdered <- colnames(df_infant_bifido_new)[order(apply(rbind(df_infant_bifido_new[,1:8],df_adult_bifido_new[,1:8],df_senior_bifido_new[,1:8]),2,function(x)(length(x[x>=0.3]))))]

summary_df_studies <- metadata_45809 %>%
  group_by(study_name, age_category,seq_type, cohort) %>%
  summarise(samples = n(), .groups = "drop")
summary_df_studies <- as.data.frame(summary_df_studies)

summary_df_studies[102,"cohort"] <- "UrbanRuralMixed"
summary_df_studies[115,"cohort"] <- "RuralTribal"

infant_study_details <- summary_df_studies[summary_df_studies$age_category == 'infant',]

adult_study_details <- summary_df_studies[summary_df_studies$age_category == 'adult',]
adult_study_details <- adult_study_details[-c(62,70),]
senior_study_details <- summary_df_studies[summary_df_studies$age_category == 'senior',]


pcoStudyPrevEuclidean <- dudi.pco(vegdist(df_combined_bifido[,MajorBifsOrdered],method="euclidean"),scannf=FALSE,nf=3)

colorInfant <- rgb(169,209,142,max=255)
colorAdult <- rgb(157,195,230,max=255)
colorSenior <- rgb(255,217,102,max=255)

s.class(pcoStudyPrevEuclidean$li,as.factor(df_combined_bifido$AgeCategory),col=c(colorInfant,colorAdult,colorSenior),plabels.col="black",plabels.cex=1.5)

s.class(pcoStudyPrevEuclidean$li,as.factor(df_combined_bifido$SequencingType),col=c("Pink","Orange"),plabels.col="black",plabels.cex=1.5)

s.class(pcoStudyPrevEuclidean$li[paste0("adult:",adult_study_details$study_name),],as.factor(df_combined_bifido[paste0("adult:",adult_study_details$study_name),"CohortType"]),col=c("darkolivegreen1","lightslateblue","hotpink1"),plabels.col="black",plabels.cex=1.2)

#Infant Bifido Plot
temp_mat <- as.matrix(df_infant_bifido_new[,MajorBifsOrdered])
hmpInfantStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row",col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantImagePrevalence.pdf",height=24,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="InfantImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

#Adult Bifido Plot
temp_mat <- as.matrix(df_adult_bifido_new[,MajorBifsOrdered])
hmpAdultStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row",col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultImagePrevalence.pdf",height=155.2,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="AdultImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

#Senior Bifido Plot
temp_mat <- as.matrix(df_senior_bifido_new[,MajorBifsOrdered])
hmpSeniorStudyBifido <- heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,dendrogram = "row", col=brewer.pal(8,"Greens"))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorImagePrevalence.pdf",height=59.2,width=40)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="SeniorImagePrevalence.svg",bg="transparent", height = svg_height, width=svg_width)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,col=brewer.pal(8,"Greens"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),dendrogram="none",labRow=NA,labCol=NA)
dev.off()

# Infant Bifido Prevalence
temp_mat <- as.matrix(df_infant_bifido_new[rev(colnames(hmpInfantStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="InfantOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="InfantSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantCohortType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="InfantCohortType.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Infant Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpInfantStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\InfantNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="InfantNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido Prevalence
temp_mat <- as.matrix(df_adult_bifido_new[rev(colnames(hmpAdultStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="AdultOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="AdultSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultCohortType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="AdultCohortType.svg",bg="transparent", height = svg_height, width=svg_width*8)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Adult Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpAdultStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\AdultNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="AdultNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido Prevalence
temp_mat <- as.matrix(df_senior_bifido_new[rev(colnames(hmpSeniorStudyBifido$carpet)),c("BifPrevalence","BifPrevalence")])
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorOverallPrevalence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="SeniorOverallPrevalence.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=brewer.pal(8,"Blues"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido SequencingType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("SequencingType","SequencingType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="wgs",1,2)))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorSequencingType.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="SeniorSequencingType.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("cornflowerblue","darkgoldenrod"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido CohortType
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("CohortType","CohortType")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x=="IndustrializedUrban",1,ifelse(x=="UrbanRuralMixed",2,3))))
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorCohortType.pdf")
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="SeniorCohortType.svg",bg="transparent", height = svg_height, width=svg_width*4)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("darkolivegreen1","hotpink1","lightslateblue"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

# Senior Bifido Association BifDetection
temp_mat <- as.matrix(df_combined_bifido[rev(colnames(hmpSeniorStudyBifido$carpet)),c("BrayAssocBifDetection","KendallAssocBifDetection")])
temp_mat <- apply(temp_mat,2,function(x)(ifelse(x<=0.05,1,0)))
colnames(temp_mat) <- c("Bray","Kendall")
#pdf(file="G:\\My Drive\\Lab\\Projects\\ProbioArchaea\\Bif_Manuscript\\AnalysisData\\SeniorNonBifInfluence.pdf",height=24,width=10)
svg_height = nrow(temp_mat)
svg_width = ncol(temp_mat)
svg(file="SeniorNonBifInfluence.svg",bg="transparent", height = svg_height, width=svg_width*2)
heatmap.2(temp_mat,density="none",trace="none",Colv=FALSE,Rowv=FALSE,col=c("white","mediumpurple3"),key=FALSE,margins=c(5,10),lhei=c(1,5),lwid=c(1,5),notecex=6,sepcolor="grey40",sepwidth=c(0.02,0.04),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)),labRow=NA,labCol=NA)
dev.off()

infant_study_vector <- rev(colnames(hmpInfantStudyBifido$carpet))
adult_study_vector <- rev(colnames(hmpAdultStudyBifido$carpet))
senior_study_vector <- rev(colnames(hmpSeniorStudyBifido$carpet))

save(infant_study_vector,adult_study_vector,senior_study_vector,file="StudyVector.RData")

######################################################################

# #Figure S1A: World Map 
# 
# For re-generating the MAP, please run the code below  ##
# if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
# if (!requireNamespace("rnaturalearth", quietly = TRUE)) install.packages("rnaturalearth")
# if (!requireNamespace("rnaturalearthdata", quietly = TRUE)) install.packages("rnaturalearthdata")
# if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
# 
# library(ggplot2)
# library(dplyr)
# library(sf)
# library(rnaturalearth)
# library(ggrepel)
# library(RColorBrewer)
# library(openxlsx)
# 
# df <- read.xlsx("/data/country_counts.xlsx")
# colnames(df) <- c("ISO3", "value")
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# world_data <- merge(
#   world,
#   df,
#   by.x = "iso_a3_eh",
#   by.y = "ISO3",
#   all.x = TRUE
# )
# 
# world_min <- world_data %>%
#   select(iso_a3_eh, name, value, geometry) %>%
#   mutate(has_data = !is.na(value), fill_country = ifelse(is.na(value), "No data", name))
# 
# 
# # Countries with data
# countries_with_data <- unique(world_min$fill_country)
# countries_with_data <- countries_with_data[countries_with_data != "No data"]
# 
# # Generate qualitative palette
# pal <- colorRampPalette(
#   brewer.pal(12, "Set3")
# )(length(countries_with_data))
# 
# names(pal) <- countries_with_data
# 
# # Add grey for no data
# pal <- c(pal, "No data" = "grey90")
# 
# points_sf <- st_point_on_surface(world_min)
# 
# label_df <- cbind(
#   as.data.frame(st_coordinates(points_sf)),
#   country = world_min$name,
#   iso_a3_eh  = world_min$iso_a3_eh,
#   value   = world_min$value,
#   has_data = world_min$has_data
# )
# 
# label_df <- label_df[!is.na(label_df$value), ]
# 
# european_iso3 <- c(
#   "FRA","DEU","ITA","ESP","CHE","AUT","NLD","BEL","LUX","DNK",
#   "SWE","FIN","NOR","ISL","IRL","GBR","POL","SVK","SVN",
#   "HUN","EST","LTU","LVA","CZE","GRC"
# )
# 
# set.seed(123)
# label_df <- label_df %>%
#   mutate(
#     is_europe = iso_a3_eh %in% european_iso3,
#     X_label = ifelse(is_europe, X + 20, X),
#     Y_label = ifelse(is_europe, Y + runif(n(), 8, 25), Y)
#   )
# 
# label_nudge_tbl <- data.frame(
#   iso_a3_eh = c("GBR","FRA","DEU","ITA","NLD","POL"),
#   nudge_x = c(-12, 8, 25, 18, 15, 18),
#   nudge_y = c(18, 20, 23, 10, 15, 13)
# )
# 
# label_df <- label_df %>%
#   left_join(label_nudge_tbl, by = "iso_a3_eh") %>%
#   mutate(
#     nudge_x = ifelse(is.na(nudge_x), 0, nudge_x),
#     nudge_y = ifelse(is.na(nudge_y), 0, nudge_y),
#     X_label = X_label + nudge_x,
#     Y_label = Y_label + nudge_y
#   )
# 
# p <- ggplot(world_min) +
#   geom_sf(
#     aes(fill = fill_country),
#     color = "grey50",
#     size = 0.4
#   ) +
#   scale_fill_manual(
#     values = pal,
#     guide = "none"
#   ) +
#   
#   # anchor points
#   geom_point(
#     data = label_df,
#     aes(x = X, y = Y),
#     size = 2.2,
#     color = "black"
#   ) +
#   
#   # arrows
#   geom_segment(
#     data = label_df,
#     aes(x = X, y = Y, xend = X_label, yend = Y_label),
#     arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
#     size = 0.6,
#     color = "grey10"
#   ) +
#   
#   # labels
#   geom_text(
#     data = label_df,
#     aes(x = X_label, y = Y_label,
#         label = paste0(country, " (", value, ")")),
#     size = 4.2
#   ) +
#   
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid = element_blank(),
#     plot.title = element_text(face = "bold", hjust = 0.5),
#     plot.margin = margin(12, 12, 12, 12)
#   ) +
#   labs(title = "Sample Count by Country")
# 
# ggsave(filename = "Sample_Count_by_Country_WorldMap.pdf", plot = p, device = pdf, width = 16, height = 9, units = "in")

####################################################################

#Figure S1B: Bifidobacterium Prevalence barplot

library(ggplot2)

percent_df <- read.csv("/data/percentage_bif_species.csv")

#----------- bar plot ------------------#

make_percentage_barplot <- function(df, filename) {
  df_long <- df %>%
    mutate(label = sprintf("%.2f", percent_studies))
  
  p <- ggplot(df_long, aes(x = reorder(Bifidobacterium, -percent_studies),
                           y = percent_studies)) +
    geom_col(fill = "slateblue") +
    geom_text(aes(label = label),
              vjust = -0.3,
              color = "black", size = 4) +
    labs(y = "Percentage of studies (%)",
         x = "Bifidobacterium species") +
    scale_y_continuous(
      breaks = seq(0, 100, 20),     # y-axis ticks every 20 up to 100
      limits = c(0, 100),           # force range to 0–100
      expand = expansion(mult = c(0, 0.05))  # remove bottom space, keep small top margin
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 12),
      legend.position = "none"
    )
  
  # Save as PDF
  ggsave(filename, plot = p, width = 18, height = 7)
}

# Run the function on percent_df
make_percentage_barplot(percent_df, "bif_percentage_barplot_new.pdf")

####################################################################

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

heatmap.2(temp_mat,density="none",trace="none",Rowv=F,Colv=F,col=c("hotpink","white","cornflowerblue"),margins=c(15,10),sepcolor="grey40",sepwidth=c(0.02,0.02),rowsep=c(0:nrow(temp_mat)),colsep=c(0:ncol(temp_mat)))

colorInfant <- rgb(169,209,142,max=255)
colorAdult <- rgb(157,195,230,max=255)
colorSenior <- rgb(255,217,102,max=255)

s.class(pcoStudyPrevEuclidean$li,as.factor(df_combined_bifido$AgeCategory),col=c(colorInfant,colorAdult,colorSenior),plabels.col="black",plabels.cex=0.5)

colorWGS <- rgb(100,149,237,max=255)
color16S <- rgb(184,134,11,max=255)

s.class(pcoStudyPrevEuclidean$li,factor(df_combined_bifido$SequencingType,levels=c("WGS","16s")),col=c(colorWGS,color16S),plabels.col="black",plabels.cex=1.5)





