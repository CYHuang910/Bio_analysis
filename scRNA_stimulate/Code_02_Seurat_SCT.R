#!/usr/bin/env Rscript
#Author : chien-yu
#Date : 20220503

#rm(list=ls(all=TRUE))
#gc()
library(optparse)

cat("\nSTART!! ",date(), "\n")
start.time <- Sys.time()
option_list <- list(
   make_option(c("-r", "--rds"), type = "character", default = FALSE, help = "Simulate RNA RDS data"),
   make_option(c("-m", "--metadata"), type = "character", default = FALSE, help = "Simulate RNA metadata"),
   make_option(c("-d", "--date"), type = "character", default = FALSE, help = "output file date" ))
#   make_option(c("-n", "--topn"), type = "integer", default = 2000, help = "select top n varalbe genes" ),

arguments <- parse_args(OptionParser(option_list=option_list))
args.r <- arguments$rds
args.m <- arguments$metadata
args.d <- arguments$date
#args.n <- arguments$topn
plot_o <- '/staging/biology/u8324386/HCC_scRNA/2.UMAP_for_scRNA/'
rds_o <- '/staging/biology/u8324386/HCC_scRNA/3.harmonization/'

library(Seurat)
library(patchwork)
library(aricode)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
# Read in simulate bulk RNA-seq normalized expression and metadata
exprs_norm = readRDS(args.r)
metadata = read.csv(args.m)

#exprs_norm <- LogNormalize(exprs_norm, scale.factor = 10000)
row.names(metadata) <- metadata[[1]]
cat("\nsimulate bulk RNA-seq :", dim(exprs_norm))
cat("\nsimulate bulk RNA-seq metadata :", dim(metadata),"\n\n")

ref_exp_full <- exprs_norm
ref_metadata <- metadata
colnames(ref_metadata)[2] <- "donor"
rm(exprs_norm, metadata)

scrna.obj <- CreateSeuratObject(counts = ref_exp_full)
donor <- ref_metadata$donor
tissue <- ref_metadata$Tissue
celltype <- ref_metadata$celltype_sub
names(donor) <- colnames(scrna.obj)
names(tissue) <- colnames(scrna.obj)
names(celltype) <- colnames(scrna.obj)

scrna.obj <- AddMetaData(object = scrna.obj, metadata = donor, col.name ="Donor")
scrna.obj <- AddMetaData(object = scrna.obj, metadata = tissue, col.name ="Tissue")
scrna.obj <- AddMetaData(object = scrna.obj, metadata = celltype, col.name ="Celltype")

## Runnning Lognormalize on a UMI matrix
scrna_log <- scrna.obj
scrna_log <- NormalizeData(scrna_log, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
scrna_log <- FindVariableFeatures(scrna_log, verbose = FALSE)
scrna_log <- ScaleData(scrna_log, verbose = FALSE)
scrna_log <- RunPCA(scrna_log, npcs = 50, verbose = FALSE)
scrna_log <- RunUMAP(scrna_log, reduction = "pca", dims = 1:50, verbose = FALSE)

## Runnning sctransform on a UMI matrix
#Split data by Donor
obj.list <- SplitObject(scrna.obj, split.by = "Donor")

scrna_1 <- obj.list[["D20171109"]]
scrna_2 <- obj.list[["D20171215"]]
scrna_3 <- obj.list[["D20180108"]]
scrna_4 <- obj.list[["D20180110"]]
scrna_5 <- obj.list[["D20180116"]]

# normalize and run dimensionality reduction on control dataset
scrna_1 <- SCTransform(scrna_1, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:50, verbose = FALSE) %>%
    FindClusters(resolution = 0.7, verbose = FALSE)

p1 <- DimPlot(scrna_1, label = T, repel = T) +
  ggtitle("Unsupervised clustering") +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=10))
p2 <- DimPlot(scrna_1, label = T, repel = T, group.by = "Celltype") +
  ggtitle("Annotated celltypes") +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=7))
p3 <- ggarrange(p1, p2,
          labels = c("A", "B"), widths = c(1, 1.5),
          ncol = 2, nrow = 1)
ggsave(paste(plot_o, args.d, "_scRNA_D20171109_after_SCT_umap.pdf", sep=''),
       p3, width=13, height=6, units="in", scale=1.5)
rm(p1,p2,p3)
gc()

scrna_2 <- SCTransform(scrna_2, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
scrna_3 <- SCTransform(scrna_3, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
scrna_4 <- SCTransform(scrna_4, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
scrna_5 <- SCTransform(scrna_5, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)
gc()
cat("\n======================================================\n",
    "   Finish SCTransform: ", date(), "\n",
    "   Now combine data. ",
    "\n======================================================\n\n")

obj.list <- list(scrna_1=scrna_1, scrna_2=scrna_2, scrna_3=scrna_3, scrna_4=scrna_4, scrna_5=scrna_5)
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

#To integrate the five datasets
immune.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
#Perform an integrated analysis
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:50, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:50)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)
cat("\nFinish combine SCT data.\n")

#Use AMI to calculate batch effect
umap_log <- scrna_log@reductions[['umap']]
umap_log <- umap_log@cell.embeddings
#kmeans
log.donor <- kmeans(umap_log, centers=5)
log.type <- kmeans(umap_log, centers=25)
ami_log.d <- AMI(log.donor$cluster, scrna_log@meta.data[['Donor']])
ami_log.c <- AMI(log.type$cluster, scrna_log@meta.data[['Celltype']])

umap_sct <- immune.combined.sct@reductions[['umap']]
umap_sct <- umap_sct@cell.embeddings
#kmeans
sct.donor <- kmeans(umap_sct, centers=5)
sct.type <- kmeans(umap_sct, centers=25)
ami_sct.d <- AMI(sct.donor$cluster, immune.combined.sct@meta.data[['Donor']])
ami_sct.c <- AMI(sct.type$cluster, immune.combined.sct@meta.data[['Celltype']])

table.ami <- matrix(NA, nrow = 2, ncol = 2)
colnames(table.ami) <- c("Donor", "Cell types")
row.names(table.ami) <- c("Log-normalize", "Sctransform")
table.ami[1,] <- c(ami_log.d,ami_log.c)
table.ami[2,] <- c(ami_sct.d,ami_sct.c)
write.csv(table.ami, paste(plot_o, args.d, "_compare_Log_SCT_AMI_table.csv", sep=''))
rm(umap_log, log.donor, log.type, ami_log.d, ami_log.c, umap_sct, sct.donor, sct.type,
   ami_sct.d, ami_sct.c, table_ami)
gc()

cat("\n========================================================\n",
    "   Now draw SCT's UMAP : ", date(), "\n",
    "\n========================================================\n")
#To visualize the 5 batch
p <- DimPlot(immune.combined.sct, reduction = "umap", split.by = "Donor", ncol=2)
ggsave(paste(plot_o, args.d, "_scRNA_each_donor_after_SCT_umap.pdf", sep=''),
       p, width=12, height=15, units="in", scale=1.5)

p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Donor")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T) +
  theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Celltype", label = T, repel = T) +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=15))
p4 <- ggarrange(p3,
                ggarrange(p1, p2, ncol = 2, labels = c("B", "C")),
                labels = "A", nrow = 2)
ggsave(paste(plot_o, args.d, "_scRNA_after_SCT_umap.pdf", sep=''),
       p4, width=12, height=12, units="in", scale=1.5)
rm(p,p1,p2,p3,p4)
gc()

cat("\n========================================================\n",
    "   Now draw Log's UMAP : ", date(), "\n",
    "\n========================================================\n")

#Visualize the clustering results on the sctransform and log-normalized
scrna_log$clusterID <- Idents(immune.combined.sct)

p1 <- DimPlot(scrna_log, reduction = "umap", group.by = "Donor") +
  ggtitle("Log-normalize (Donor)") +  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Donor") +
  ggtitle("Sctransform (Donor)") +
  theme(plot.title = element_text(hjust = 0.5))
p_d <- ggarrange(p1, p2,
                labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave(paste(plot_o, args.d, "_compare_log_SCT_donor_umap.pdf", sep=''),
       p_d, width=12, height=6, units="in", scale=1.5)

p3 <- DimPlot(scrna_log, reduction = "umap", group.by = "Celltype", label = T, repel = T) +
  ggtitle("Log-normalize (Celltype)") +
  theme(plot.title = element_text(hjust = 0.5),legend.text=element_text(size=7))
p4 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Celltype", label = T, repel = T) +
  ggtitle("Sctransform (Celltype)") +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=7))
p_c <- ggarrange(p3, p4,
                 labels = c("A", "B"), ncol = 1, nrow = 2)
ggsave(paste(plot_o, args.d, "_compare_log_SCT_celltype_umap.pdf", sep=''),
       p_c, width=9, height=10, units="in", scale=1.5)

p5 <- DimPlot(scrna_log, reduction = "umap", group.by = "clusterID", label = T, repel = T) +
  ggtitle("Log-normalize (seurat_cluster)")
p6 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = T, repel = T) +
  ggtitle("Sctransform (seurat_cluster)") +
  theme(plot.title = element_text(hjust = 0.5))
p_sc <- ggarrange(p5, p6,
                  labels = c("A", "B"), ncol = 1, nrow = 2)
ggsave(paste(plot_o, args.d, "_compare_log_SCT_seurat_cluster_umap.pdf", sep=''),
       p_sc, width=6, height=10, units="in", scale=1.5)

#write.csv(sct_exp, paste("/staging/biology/u8324386/HCC_scRNA/3.harmonization/", file_name, ".csv", sep=""))
saveRDS(immune.combined.sct, paste(plot_o, args.d, "_scRNA_after_SCT.rds", sep=''))
saveRDS(scrna_log, paste(plot_o, args.d, "_scRNA_after_Lognormalize.rds", sep=''))

end.time <- Sys.time()
time.taken <- end.time - start.time

cat("\n======================================================\n",
    "  All code excute: ", time.taken, "hr",
    "\n======================================================\n\n")

