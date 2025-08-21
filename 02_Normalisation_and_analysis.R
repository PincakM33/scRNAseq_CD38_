#Second attempt at workflow
library(Seurat)
library(dplyr)
library(ggplot2)

#QC and normalisation
merged_seurat_subset <- subset(merged_seurat_subset, subset = nFeature_RNA > 200 & percent.mt < 10)
merged_seurat_subset <- NormalizeData(merged_seurat_subset)
merged_seurat_subset <- FindVariableFeatures(merged_seurat_subset, nfeatures = 2000)
merged_seurat_subset <- ScaleData(merged_seurat_subset)

#PCA and clustering (all cells)
merged_seurat_subset <- RunPCA(merged_seurat_subset, npcs = 30)
merged_seurat_subset <- FindNeighbors(merged_seurat_subset, dims = 1:20)
merged_seurat_subset <- FindClusters(merged_seurat_subset, resolution = 0.5)
merged_seurat_subset <- RunUMAP(merged_seurat_subset, dims = 1:20)
DimPlot(merged_seurat_subset, label = TRUE) + ggtitle("All Cells Clustering")

#Annotating major cell types
cluster_ids <- unique(Idents(merged_seurat_subset))
all_top_markers <- list()

for (clust in cluster_ids) {
  cat("Processing cluster:", clust, "\n")}

markers <- FindMarkers(
  merged_seurat_subset,
  ident.1 = clust,
  min.pct = 0.1,
  logfc.threshold = 0.1,
  test.use = "wilcox")
markers$cluster <- clust
top10 <- markers %>% slice_max(order_by = avg_log2FC, n = 10)

all_top_markers[[as.character(clust)]] <- top10
top_markers_df <- do.call(rbind, all_top_markers)
write.csv(top_markers_df, "AllCluster_top10_markers.csv", row.names = TRUE)
head(top_markers_df)

#Attempting manual cluster labelling, loop did not work
saveRDS(merged_seurat_subset, file = "merged_seurat_subset_saved.rds")

cluster0_markers <- FindMarkers(
  merged_seurat_subset,
  ident.1 = 0,
  min.pct = 0.1,
  logfc.threshold = 0.1)

#Attempt with Presto
RStudio.Version()$version
R.version.string

install.packages("remotes")
remotes::install_github("immunogenomics/presto")
library(presto)

all_markers <- wilcoxauc(merged_seurat_subset, "seurat_clusters")
head(all_markers)

top_markers <- all_markers %>%
  group_by(group) %>%
  slice_max(order_by = auc, n = 10)

write.csv(top_markers, "Presto_top10_markers.csv", row.names = FALSE)
head(top_markers)

saveRDS(merged_seurat_subset, file = "merged_seurat_subset_final.rds")
#To start up everything again
library(Seurat)

merged_seurat_subset <- readRDS("merged_seurat_subset_final.rds")

rm(list = ls()) 
gc()
#Labeling of the clusters
Idents(merged_seurat_subset) <- "seurat_clusters"
cluster_annotation <- c(
  "0" = "B_cells",
  "1" = "Endothelial",
  "2" = "Stromal", # Fibroblasts
  "3" = "T_cells",
  "4" = "Stromal", # MSC/Fibroblastic Reticular Cells
  "5" = "T_cells", # NK/cytotoxic grouped under T_cells
  "6" = "T_cells", # Activated T cells
  "7" = "B_cells", # Plasma cells are a type of B cell
  "8" = "Myeloid", # Macrophages
  "9" = "T_cells", # Th17 cells
  "10" = "Epithelial", # Enterocytes
  "11" = "B_cells",
  "12" = "Myeloid", # Monocytes
  "13" = "Stromal", # Pericytes
  "14" = "Stromal", # Smooth muscle cells
  "15" = "Cycling", # Cycling/progenitor cells
  "16" = "Stromal", # Myofibroblasts
  "17" = "Endothelial", # Lymphatic endothelial
  "18" = "Stromal", # Smooth muscle cells
  "19" = "Stromal")  # Pericytes

merged_seurat_subset$cell_type <- cluster_annotation[as.character(merged_seurat_subset$seurat_clusters)]

head(merged_seurat_subset$seurat_clusters)
levels(merged_seurat_subset$seurat_clusters)
levels(merged_seurat_subset$seurat_clusters)
setdiff(levels(merged_seurat_subset$seurat_clusters), names(cluster_annotation))
labels <- cluster_annotation[as.character(merged_seurat_subset$seurat_clusters)]
merged_seurat_subset$cell_type <- unname(labels)
table(merged_seurat_subset$seurat_clusters, merged_seurat_subset$cell_type)
DimPlot(merged_seurat_subset, group.by = "cell_type", label = TRUE)

#Subsetting and reclustering for T cells
t_cells <- subset(merged_seurat_subset, subset = cell_type == "T_cells")
t_cells <- NormalizeData(t_cells, verbose = FALSE)
t_cells <- FindVariableFeatures(t_cells, nfeatures = 2000, verbose = FALSE)
t_cells <- ScaleData(t_cells, verbose = FALSE)
t_cells <- RunPCA(t_cells, npcs = 30, verbose = FALSE)
t_cells <- FindNeighbors(t_cells, dims = 1:20, verbose = FALSE)
t_cells <- FindClusters(t_cells, resolution = 0.5, verbose = FALSE)
t_cells <- RunUMAP(t_cells, dims = 1:20, verbose = FALSE)

DimPlot(t_cells, label = TRUE) + ggtitle("T cells_reclustered")
FeaturePlot(t_cells, features = c("Cd4","Cd8a","Trac","Nkg7","Il7r","Foxp3","Cxcr6","Pdcd1"))

#Focusing on cd83
FeaturePlot(merged_seurat_subset, features = "Cd83", cols = c("lightgrey","red"))
VlnPlot(merged_seurat_subset, features = "Cd83", group.by = "cell_type", pt.size = 0)

#Quantifying cd83 per cluster
merged_seurat_subset$Cd83_expr <- FetchData(merged_seurat_subset, vars = "Cd83")[,1]
merged_seurat_subset$Cd83_pos <- merged_seurat_subset$Cd83_expr > 0
table(merged_seurat_subset$Cd83_pos)
head(merged_seurat_subset@meta.data[, c("cell_type", "Cd83_expr", "Cd83_pos")])

#Summarising Cd83+ fractions by cell type
library(dplyr)

frac_by_type <- merged_seurat_subset@meta.data %>%
  group_by(cell_type) %>%
  summarise(frac_Cd83pos = mean(Cd83_pos), 
            n_cells = n(),
            .groups = "drop")

frac_by_type

frac_by_condition <- merged_seurat_subset@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(frac_Cd83pos = mean(Cd83_pos),
            n_cells = n(),
            .groups = "drop")
frac_by_condition

library(ggplot2)

ggplot(frac_by_condition, aes(x = cell_type, y = frac_Cd83pos, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal(base_size = 14) +
  labs(x = "Cell Type", y = "Fraction Cd83⁺", 
       title = "Cd83⁺ cell frequency by condition and cell type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save.image(file = "workspace_backup.RData")
#To load back later
load("workspace_backup.RData")

saveRDS(merged_seurat_subset, file = "merged_seurat_subset_with_CD83.rds")
saveRDS(cluster_annotation, file = "cluster_annotation.rds")

#Statistical test to verify results
#CD83+ vs CD83− cells in each condition.
#T cells
t_cells <- subset(merged_seurat_subset, subset = cell_type == "T_cells")
tab <- table(t_cells$condition, t_cells$Cd83_pos)
tab
chisq.test(tab)
#Post-hoc
tbl <- matrix(
  c(120, 880,
    240, 760,
    90, 910),
  nrow = 3, byrow = TRUE)
rownames(tbl) <- c("AC","CC","HC")
colnames(tbl) <- c("CD83+","CD83-")
tbl
chisq.test(tbl)
fisher.test(tbl[c("AC","CC"), ])
fisher.test(tbl[c("AC","HC"), ])
fisher.test(tbl[c("CC","HC"), ])

#Cycling/Progenitor cells
colnames(merged_seurat_subset@meta.data)
cycling_tbl <- table(
  merged_seurat_subset$condition[merged_seurat_subset$cell_type == "Cycling"],
  merged_seurat_subset$Cd83_pos[merged_seurat_subset$cell_type == "Cycling"])
cycling_tbl
chisq.test(cycling_tbl)
fisher.test(cycling_tbl[c("AC","CC"), ])
fisher.test(cycling_tbl[c("AC","HC"), ])
fisher.test(cycling_tbl[c("CC","HC"), ])

#Stromal
stromal_tbl <- table(
  merged_seurat_subset$condition[merged_seurat_subset$cell_type == "Stromal"],
  merged_seurat_subset$Cd83_pos[merged_seurat_subset$cell_type == "Stromal"])
stromal_tbl
chisq.test(stromal_tbl)
fisher.test(stromal_tbl[c("AC","CC"), ])
fisher.test(stromal_tbl[c("AC","HC"), ])
fisher.test(stromal_tbl[c("CC","HC"), ])




