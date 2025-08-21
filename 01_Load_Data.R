#Project: scRNA-seq Analysis of T Cell Subsets in Murine Colitis
#Author: PincakM33
# Date: 18/08/25
# Description: This script performs a single-cell RNA-seq analysis of a mouse
#              colitis dataset (GSE264408) to characterize T cell subsets and
#              their expression of CD38.

#          Section 1: Data Loading and Preparation
# This section loads the raw data and metadata, merges the samples,
# and prepares the Seurat object for downstream analysis.

# Load necessary libraries
library(Seurat)
library(tidyverse)

#      Load Metadata
# Read the metadata from the 'data' folder
metadata <- read.csv("GSE264408_metadata.csv")
View(metadata)

#Creating Seurat objects 
GSM8217725.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217725_Normal1")
GSM8217725 <- Seurat::CreateSeuratObject(counts = GSM8217725.data, project = "GSM8217725")

GSM8217726.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217726_Normal2")
GSM8217726 <- Seurat::CreateSeuratObject(counts = GSM8217726.data, project = "GSM8217726")

GSM8217727.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217727_Normal3")
GSM8217727 <- Seurat::CreateSeuratObject(counts = GSM8217727.data, project = "GSM8217727")

GSM8217728.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217728_AcuteColitis1")
GSM8217728 <- Seurat::CreateSeuratObject(counts = GSM8217728.data, project = "GSM8217728")

GSM8217729.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217729_AcuteColitis2")
GSM8217729 <- Seurat::CreateSeuratObject(counts = GSM8217729.data, project = "GSM8217729")

GSM8217730.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217730_AcuteColitis3")
GSM8217730 <- Seurat::CreateSeuratObject(counts = GSM8217730.data, project = "GSM8217730")

GSM8217731.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217731_ChronicColitis1")
GSM8217731 <- Seurat::CreateSeuratObject(counts = GSM8217731.data, project = "GSM8217731")

GSM8217732.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217732_ChronicColitis2")
GSM8217732 <- Seurat::CreateSeuratObject(counts = GSM8217732.data, project = "GSM8217732")

GSM8217733.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217733_ChronicColitis3")
GSM8217733 <- Seurat::CreateSeuratObject(counts = GSM8217733.data, project = "GSM8217733")

GSM8217734.data <- Seurat::Read10X(data.dir = "raw_data/GSM8217734_ChronicColitis4")
GSM8217734 <- Seurat::CreateSeuratObject(counts = GSM8217734.data, project = "GSM8217734")

seurat_list <- list(GSM8217725, GSM8217726, GSM8217727, GSM8217728, GSM8217729, GSM8217730 ,
  GSM8217731, GSM8217732, GSM8217733, GSM8217734)

merged_seurat <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)])

#Attempt to add metadata
correct_lookup <- data.frame(
  short_gsm = c("GSM8217725", "GSM8217726", "GSM8217727", "GSM8217728", "GSM8217729", "GSM8217730", "GSM8217731", "GSM8217732", "GSM8217733", "GSM8217734"),
  metadata_name = c("HC1", "HC2", "HC3", "AC1", "AC2", "AC3", "CC1", "CC2", "CC3", "CC4")
)
merged_seurat$metadata_name <- correct_lookup$metadata_name[match(merged_seurat$orig.ident, correct_lookup$short_gsm)]
table(merged_seurat$metadata_name)
merged_seurat$condition <- metadata$condition[match(merged_seurat$metadata_name, metadata$sample)]
table(merged_seurat$condition)
#Save
save.image(file = "project_data_ready_for_QC.RData")
#To load in for next session
load("project_data_ready_for_QC.RData")




