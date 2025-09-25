# Project_Bifidobacterium
Identifying the non-Bifidobacterial hallmarks of a Bifidobacterium-receptive gut microbiome

This repository contains the R scripts, input data, and supporting files used in the analyses described in the manuscript.  
The workflow is organized into five sections corresponding to different analyses.

# Section 1: Bifidobacterial Prevalence Patterns
- Scripts:  
  - Bifido_Type_Batch1.R  
  - Bifido_Type_Batch2.R  
  - Bifido_Type_Batch3.R  
- Input: RData files from "data/"
- Purpose: Computes prevalence and distribution patterns of _Bifidobacterium_ types.

# Section 2: Overall Association Patterns
- Scripts:  
  - bif_association_patterns.R  
  - Figures_pcoa.R  
  - rf_association_heatmaps.R  
- Input: RData files from "data/"
- Purpose: Identifies global association patterns, generates PCoA plots and heatmaps.

# Section 3: Disease Association Analysis
- Script:  
  - DiseaseAssociationAnalysis.R 
- Input: RData files from "data/"
- Purpose: Computation of Association-Scores for different Disease-groups.

# Section 4: Intervention Analysis
- Scripts:  
  - Intervention_analysis_Part1.R 
  - Intervention_analysis_Part2.R  
- Input: RData files from "data/"  
- Purpose: Examines intervention trial datasets for receptive microbiome shifts.

# Section 5: Functional Analysis
- Scripts:  
  - FunctionalAnalysis_final.R 
- Input: RData files from "data/" 
  - Protein FASTA files for functional groups are provided in 5 zip folders
- Purpose: Characterizes non-Bifidobacterial protein functional groups.  
