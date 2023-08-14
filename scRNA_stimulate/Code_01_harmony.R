#!/usr/bin/env Rscript
#Author : chien-yu
#Date : 20220503

library(optparse)
start.time <- Sys.time()
option_list <- list(
   make_option(c("-r", "--rds"), type = "character", default = FALSE, help = "Simulate RNA RDS data"),
   make_option(c("-m", "--metadata"), type = "character", default = FALSE, help = "Simulate RNA metadata"),
   make_option(c("-o", "--output"), type = "character", default = FALSE, help = "output file direction" ),
   make_option(c("-n", "--topn"), type = "integer", default = 2000, help = "select top n varalbe genes" ),
   make_option(c("-p", "--PC"), type = "integer", default = 20, help = "how many pc be used"),
   make_option(c("-t", "--theta"), type = "integer", default = 2, help = "cluster diversity enforcement"),
   make_option(c("-k", "--cluster"), type = "integer", default = 100, help = "number of clusters in Harmony model"))

arguments <- parse_args(OptionParser(option_list=option_list))
args.r <- arguments$rds
args.m <- arguments$metadata
args.o <- arguments$output
args.n <- arguments$topn
args.p <- arguments$PC
args.t <- arguments$theta
args.k <- arguments$cluster


suppressPackageStartupMessages({
    source('/staging/biology/u8324386/HCC_scRNA/libs.R') # imports
#    source('/staging/biology/u8324386/HCC_scRNA/utils.R') # color definitions and plotting functions
})

# Read in simulate bulk RNA-seq normalized expression and metadata
exprs_norm = readRDS(args.r)
metadata = read.csv(args.m)

#exprs_norm <- LogNormalize(exprs_norm, scale.factor = 10000)
row.names(metadata) <- metadata[[1]]
cat("\nsimulate bulk RNA-seq :", dim(exprs_norm), "\n")
cat("simulate bulk RNA-seq metadata :", dim(metadata), "\n")

ref_exp_full <- exprs_norm
ref_metadata <- metadata
colnames(ref_metadata)[2] <- "donor"
rm(exprs_norm, metadata)
gc()

## Function to find variable genes using variance stabilizing transform (vst) method
var_genes = vargenes_vst(ref_exp_full, groups = as.character(ref_metadata[['donor']]), topn = args.n)
ref_exp = ref_exp_full[var_genes, ]
cat("variable gene :", dim(ref_exp), "\n")
num_gene <- dim(ref_exp)[1]

#vst_exp <- t(ref_exp)
#umap_vst <- umap(vst_exp)
#df <- data.frame(x = umap_vst$layout[,1],
#                 y = umap_vst$layout[,2])
#df_pick <- ref_metadata[,c(2,4)]
#df <- cbind(df, df_pick)

vst_shape <- c(seq(0,14),33,35,36,37,38,60,61,62,63,64)
#p = ggplot(df, aes(x, y)) +
#   geom_point(alpha=0.5, aes(shape = celltype_sub, colour = donor)) +
#   scale_shape_manual(values = vst_shape) +
#   scale_colour_manual(values = c("#F8766D","#A3A500","#00BF7D","#A6CEE3","#E76BF3")) +
#   ggtitle(paste(" scRNA expression \n (variable genes : ", num_gene, ", before harmony)", sep="")) +
#   theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=7))
#ggsave(paste("/staging/biology/u8324386/HCC_scRNA/2.UMAP_for_scRNA/0602_scRNA_expression_pick_variable_genes_umap.pdf", sep="" ), p, width=4, height=4, units="in", scale=2)
#rm(vst_exp, umap_vst, df, df_pick, p)
#gc()

####
vargenes_means_sds = tibble::tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
vargenes_means_sds$stddev = symphony::rowSDs(ref_exp, vargenes_means_sds$mean)
head(vargenes_means_sds)

#Scale data using calculated gene means and standard deviations
ref_exp_scaled = symphony::scaleDataWithStats(ref_exp, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
#ref_exp_scaled[1:5, 1:3]

####
set.seed(0)
s = irlba::irlba(ref_exp_scaled, nv = args.p)
Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
loadings = s$u

set.seed(0)
ref_harmObj = harmony::HarmonyMatrix(
  data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
  meta_data = ref_metadata, ## dataframe with cell labels
  theta = c(args.t),             ## cluster diversity enforcement
  vars_use = c('donor'),    ## variable to integrate out
  nclust = args.k,             ## number of clusters in Harmony model
  max.iter.harmony = 20,
  return_object = TRUE,     ## return the full Harmony model object
  do_pca = FALSE            ## don't recompute PCs
)

reference = symphony::buildReferenceFromHarmonyObj(
  ref_harmObj,            # output object from HarmonyMatrix()
  ref_metadata,           # reference cell metadata
  vargenes_means_sds,     # gene names, means, and std devs for scaling
  loadings,               # genes x PCs matrix
  verbose = TRUE,         # verbose output
  do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
  save_uwot_path = paste(args.o, "uwot_model_1", sep=""))

reference$normalization_method = 'log(CP10k+1)'
#x <- reference$loadings %*% reference$Z_corr
#gene <- reference$vargenes$symbol
#l <- which(reference$meta_data$cell_subtype=="X")

file_name <- paste("0901_scRNA_25type_top", args.n, "_PC", args.p, "_t", args.t, "_K", args.k, sep="")
saveRDS(reference, paste(args.o, file_name, "_reference.rds", sep=""))


##---- UMAP ----##
df <- data.frame(x = reference$umap$embedding[,1],
                 y = reference$umap$embedding[,2])
df_pick <- ref_metadata[,c(2,4)]
df <- cbind(df, df_pick)
dim(df)
type_num <- seq(1,length(unique(ref_metadata$celltype_sub)))

#kmeans
kmeans.h.donor <- kmeans(df[,c(1,2)], centers=5)
kmeans.h.type <- kmeans(df[,c(1,2)], centers=25)
a<-AMI(kmeans.h.donor$cluster, df$donor)
b<-AMI(kmeans.h.type$cluster, df$celltype_sub)

table.kmeans <- matrix(NA, nrow = 1, ncol = 2)
colnames(table.kmeans) <- c("Donor", "Cell types")
row.names(table.kmeans) <- file_name
table.kmeans[1,] <- c(a,b)
write.csv(table.kmeans, paste("/staging/biology/u8324386/HCC_scRNA/4.kmeans/", file_name, "_AMI.csv", sep=""))
write.csv(kmeans.h.donor$cluster, paste("/staging/biology/u8324386/HCC_scRNA/4.kmeans/", file_name, "_donor_cluster.csv", sep=""))
write.csv(kmeans.h.type$cluster, paste("/staging/biology/u8324386/HCC_scRNA/4.kmeans/", file_name, "_type_cluster.csv", sep=""))

##--- export the plot ---##
p = ggplot(df, aes(x, y)) +
  geom_point(alpha=0.5, aes(shape = celltype_sub, colour = donor)) +
  scale_shape_manual(values = vst_shape) +
  scale_colour_manual(values = c("#F8766D","#A3A500","#00BF7D","#A6CEE3","#E76BF3")) +
  ggtitle(paste(file_name, "\n (", num_gene, " genes, after harmony)", sep="")) +
  labs(tag = paste("Donor's AMI: ", a,"\n\n Cell type's AMI: ", b, sep = "")) +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=8),
        plot.tag.position = c(0.8, 0.9))

ggsave(paste("/staging/biology/u8324386/HCC_scRNA/2.UMAP_for_scRNA/", file_name, "_umap.pdf", sep="" ),
       p, width=6, height=5, units="in", scale=2)

end.time <- Sys.time()
time.taken <- end.time - start.time
cat("\nfile name :", file_name, "\n")
cat(time.taken, "\n")

rm(df, df_pick, type_num, p)
rm(s, reference, ref_harmObj, Z_pca_ref)
gc()
