#!/usr/bin/env Rscript
#Author : chien-yu
#Date : 20220503

library(optparse)
start.time <- Sys.time()
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = FALSE, help = "Simulate RNA RDS data"),
  make_option(c("-m", "--metadata"), type = "character", default = FALSE, help = "Simulate RNA metadata"),
  make_option(c("-o", "--output"), type = "character", default = FALSE, help = "output for bulk RNA-seq direction"))

arguments <- parse_args(OptionParser(option_list=option_list))
args.r <- arguments$rds
args.m <- arguments$metadata
args.o <- arguments$output

suppressPackageStartupMessages({
  source('/staging/biology/u8324386/HCC_scRNA/libs.R') # imports
  source('/staging/biology/u8324386/HCC_scRNA/utils.R') # color definitions and plotting functions
})

name <- gsub("^.*/", "", args.r)
name <- str_split(name, "\\.")[[1]][1]
#name <- gsub("\\_reference", "", name)
name

##--------  先找出scRNA 種類含量 ------
ref_exp_full = readRDS('/staging/biology/u8324386/HCC_scRNA/GSE140228_log_droplet_data.rds')
ref_metadata = read.csv('/staging/biology/u8324386/HCC_scRNA/GSE140228_droplet_25_cellinfo.csv')
row.names(ref_metadata) <- ref_metadata[[1]]

unique_table <- as.data.frame(table(ref_metadata$celltype_sub))
colnames(unique_table) <- c("celltype", "All_donor")
donor_table <- as.data.frame(table(ref_metadata$donor))
colnames(donor_table) <- c("donor", "sample_count")

summary_type_table = unique_table
for (i in c(1:nrow(donor_table))) {
  a <- as.data.frame(table(ref_metadata[which(ref_metadata$donor == donor_table$donor[i]),]$celltype_sub))
  colnames(a) <- c("celltype", paste(donor_table$donor[i]))
  summary_type_table <- merge(summary_type_table, a, by = "celltype", all = T)
}

df_type = summary_type_table[,-2]
for (i in 2:ncol(df_type)) {
  df_type[,i] <- df_type[,i]/sum(df_type[,i])
}
df_type[,2:6] <- round(df_type[,2:6], digits = 2)
df_type[df_type == 0.0] <- 0.01
rm(a, i, ref_exp_full, ref_metadata, unique_table)
gc()

##-------- 進入模擬樣本階段 ---------
##--- 25 cell type ---
##--- scRNA SCTransform ---
ref_exp_edit <- readRDS(args.r)
ref_metadata = read.csv('/staging/biology/u8324386/HCC_scRNA/GSE140228_droplet_25_cellinfo.csv')
 
row.names(ref_metadata) <- ref_metadata[[1]]
colnames(ref_metadata)[2] <- "donor"
cat("\nref_sct_exp: ",dim(ref_exp_edit), "\n")

cat("\nref_sct_exp: ",dim(ref_exp_edit), "\n")

#move the value
ref_exp_edit[1:3, 1:2]
cat("\ncheck matrix negative value (before): ", sum(apply( ref_exp_edit , 2 , function(x) sum ( x < 0 ) )), "\n")
#留下和 CGL2 病人交集的基因
cgl2_sct_gene <- read.table('/staging/biology/u8324386/HCC_scRNA/cgl2_sct_gene.txt')
cgl2_sct_gene <- unlist(cgl2_sct_gene, use.names = FALSE)
ref_exp_edit <- ref_exp_edit[cgl2_sct_gene, ]
cat("\nref_exp_edit remain cgl2 gene\n", dim(ref_exp_edit))
ref_exp_edit[1:3,1:2]
gc()

cat("\n=============================================\n",
    "    START to simulate bulk RNA-seq: \n    ",
    date(),
    "\n=============================================\n")

unique_type <- sort(unique(ref_metadata$celltype_sub))
#cat("\ncell type: ", unique_type, "\n")
celltype_num <- length(unique(ref_metadata$celltype_sub))
one_donor_sample_num <- 8000
over <- c(2,4,5,6,7,8,9,11,12,14,15,20,21,22,23,25)

###--------- 根據原始樣本細胞種類含量 --------
wall_table <- matrix(NA, nrow = 1, ncol = celltype_num)
for (i in c(1:nrow(donor_table))) {
  percent_donor <- df_type[, i+1]
  for (n in c(1:one_donor_sample_num)) {
    sample_wall = c()
    for (c in c(1:celltype_num)) {
      if (c %in% over) {
        #找出有超過 5% 的 16 種細胞
        #樣本數做本身x(1+-10%)
        tmp_w <- percent_donor[c]
        tmp_w <- tmp_w*(1 + sample(seq(-10, 10, 1)/100, size = 1))
        sample_wall[c] <- tmp_w
      }else{
        tmp_w <- percent_donor[c]
        tmp_w <- tmp_w*(1 + sample(seq(-10, 0, 1)/100, size = 1))
        sample_wall[c] <- tmp_w
      }
    }
    wall_table <- rbind(wall_table,sample_wall)
  }
  print(paste("simulate", i, "donor"))
}
wall_table <- wall_table[-1, ]

##調整含量已符合，每個樣本數都為 400 個細胞
n_cell = 400
new_wall <- wall_table*n_cell
new_wall <- round(new_wall, digits = 0)
for (i in c(1:nrow(new_wall))) {
  row.names(new_wall)[i] <- paste("Simulate_RNA_after_scRNA_SCT_", i, sep = "")
}
colnames(new_wall) <- df_type$celltype

for (i in 1:nrow(new_wall)) {
  check <- sum(new_wall[i,])
  if (check > n_cell) {
    c_max = check- n_cell
    if (c_max > 25) {
      new_wall[i,] <- new_wall[i,] -1
      tmp <- sample(seq(1, 25, 1), size = c_max - 25, replace = F)
      new_wall[i,tmp] <- new_wall[i,tmp] -1
    }else{
      tmp <- sample(seq(1, 25, 1), size = c_max, replace = F)
      new_wall[i, tmp] <- new_wall[i, tmp] -1
    }
  }
  if(check < n_cell) {
    c_min = n_cell - check
    if (c_min > 25) {
      new_wall[i,] <- new_wall[i,] +1
      tmp <- sample(seq(1, 25, 1), size = c_min-25, replace = F)
      new_wall[i,tmp] <- new_wall[i,tmp] +1
    }else{
      tmp <- sample(seq(1, 25, 1), size = c_min, replace = F)
      new_wall[i, tmp] <- new_wall[i, tmp] + 1
    }
  }
}
df_sum <- as.data.frame(rowSums(new_wall))
a<- which(df_sum[,1] == n_cell)
cat("\new_wall 400 cells: ", length(a), "\n")
cat("\ncheck matrix negative value: ", sum(apply(new_wall , 2 , function(x) sum ( x < 0 ) )), "\n")

wall_table <- new_wall
rm(a,c,i,n,tmp,tmp_w, check, new_wall, df_sum, c_min,c_max)
gc()
#write.table(wall, paste(args.o, "0720_SCT_filter_wall_list.txt", sep=""))
#wall[1:25]


###--------- 模擬 Bulk RNA-seq --------
#必須同一 batch 才可以混合
#表現量資料為去除Ｘ類別的樣本數: 55072
simulate_rna_with_batch <- matrix(NA, nrow = nrow(ref_exp_edit), ncol = 1)
pick_cell_id <- matrix(NA, nrow = nrow(wall_table), ncol = celltype_num)
meta_table <- matrix(NA, nrow = nrow(wall_table), ncol = 2)

for (i in c(1:nrow(donor_table))) {
  diff_donor <- ref_metadata[which(ref_metadata$donor == donor_table$donor[i]),]
  
  for (n1 in c(1:one_donor_sample_num)) {
    proportion = c()
    for (n in c(1:celltype_num)) {
      #Find each types of locations
      location <- which(diff_donor$celltype_sub == unique_type[n])
      check <- wall_table[one_donor_sample_num*(i-1)+ n1, n]
      #random pick sample
      #若是要取的數目 > 實際的樣本數 --> 則重複抽樣
      if (check > length(location)) {
        tmp_id <- sample(location, size = wall_table[one_donor_sample_num*(i-1)+ n1, n], replace = T)
      }
      else{
        tmp_id <- sample(location, size = wall_table[one_donor_sample_num*(i-1)+ n1, n], replace = F)
      }
      proportion <- append(proportion, tmp_id)
      pick_cell_id[one_donor_sample_num*(i-1) + n1, n] <- str_c(tmp_id, collapse = ",")
      meta_table[one_donor_sample_num*(i-1) + n1, 2] <- donor_table$donor[i]
    }
    tmp <- ref_exp_edit[ ,proportion]
    data <- as.data.frame(rowSums(tmp))
    simulate_rna_with_batch <- cbind(simulate_rna_with_batch, data)
    if (n1%%2000 == 0) {
      print(paste("simulate", i, "donor", n1, "sample"))
    }
  }
}
simulate_rna_with_batch <- simulate_rna_with_batch[ ,-1]
for (i in c(1:ncol(simulate_rna_with_batch))) {
  colnames(simulate_rna_with_batch)[i] <- paste("Simulate_RNA_after_scRNA_SCT_", i, sep = "")
}
row.names(simulate_rna_with_batch) <- row.names(ref_exp_edit)
cat("\nsimulate_rna_with_batch: ", dim(simulate_rna_with_batch), "\n")
simulate_rna_with_batch[1:3,1:2]

meta_table[,1] <- colnames(simulate_rna_with_batch)
ccc = c()
for (i in 1:5) {
  cc <- c(replicate(one_donor_sample_num, as.character(donor_table$donor[i])))
  ccc <- append(ccc,cc)
}
meta_table[,2] <- ccc

colnames(meta_table) <- c("barcode", "donor")
row.names(pick_cell_id) <- row.names(wall_table) <- colnames(simulate_rna_with_batch)
colnames(pick_cell_id) <- colnames(wall_table) <- unique_type

#Save SCT data
write.csv(wall_table, paste(args.o, "-1.Simulate_RNA_after_SCT_filter_wall_table.csv", sep=""), row.names = T)
write.csv(pick_cell_id, paste(args.o, "-2.Simulate_RNA_after_SCT_filter_cell_id_table.csv", sep=""), row.names = T)
write.csv(meta_table, paste(args.o, "-3.Simulate_RNA_after_SCT_filter_metadata.csv", sep=""),row.names = F)
write.csv(simulate_rna_with_batch, paste(args.o, "_Simulate_RNA_after_SCT_filter.csv", sep=""), row.names = T)

rm(i,n,n1,tmp_id,w,pick_cell_id,wall_table,cc,ccc)
rm(tmp,data,check,celltype_num)
rm(location,proportion,unique_type,one_donor_sample_num)
gc()

##---- UMAP ----##
umap_wd <- "/staging/biology/u8324386/HCC_scRNA/2.UMAP_for_scRNA/0901"
cat("\nStart to draw umap: ", date())
num_gene <- dim(ref_exp_edit)[1]
simulate_exp <- t(simulate_rna_with_batch)
umap_s <- umap(simulate_exp)
df <- data.frame(x = umap_s$layout[,1],
                 y = umap_s$layout[,2])
df_pick <- meta_table[,2]
df <- cbind(df, df_pick)
colnames(df)[3] <- "donor"

kmeans.donor <- kmeans(df[,c(1,2)], centers=5)
a<-AMI(kmeans.donor$cluster, df$donor)
#u_name <- gsub("0529_scRNA_25type_","",name)

cat("\nExport pdf: ", dim(df), "\n")
p = ggplot(df, aes(x, y)) +
  geom_point(alpha=0.6, size=3, aes(colour = donor)) +
  scale_colour_manual(values = c("#F8766D","#A3A500","#00BF7D","#A6CEE3","#E76BF3")) +
  ggtitle(paste(" Simulate bulk RNA \n",
                "after scRNA SCT \n", name, "\n",
                "(", num_gene, " x 8000)", sep="")) +
  labs(tag = paste("Donor's AMI: \n", a, sep = "")) +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=10),
        plot.tag.position = c(0.9, 0.97))
ggsave(paste(umap_wd, "_Simulate_RNA_after_SCT_filter_umap.pdf", sep=""),
       p, width=6, height=6, units="in", scale=1.7)

rm(umap_s,df,df_pick,kmeans.donor,a,p)
gc()

x <- simulate_rna_with_batch
row.names(x) <- NULL
x <- Matrix(as.matrix(x), sparse = TRUE)
row.names(x) <- row.names(simulate_rna_with_batch)
colnames(x) <- colnames(x)
obj <- CreateSeuratObject(x, project = "Simulate_RNA_from_HCC_droplet")
#obj[["RNA"]]@data[1:4,1:3]
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
#obj[["RNA"]]@data[1:4,1:3]

cat("\nDraw umap after lognormalize: ", date())
num_gene <- dim(ref_exp_edit)[1]
x_2 <- t(as.matrix(obj[["RNA"]]@data))
x_2[is.nan(x_2)] = 0
cat("\nx_2: ",dim(x_2))
umap_x <- umap(x_2)
dfx <- data.frame(x = umap_x$layout[,1],
                  y = umap_x$layout[,2])
dfx_pick <- meta_table[,2]
dfx <- cbind(dfx, dfx_pick)
colnames(dfx)[3] <- "donor"

kmeans.donor_x <- kmeans(dfx[,c(1,2)], centers=5)
a_x<-AMI(kmeans.donor_x$cluster, dfx$donor)

cat("\nExport pdf after lognormalize: ", dim(df))
p_x = ggplot(dfx, aes(x, y)) +
  geom_point(alpha=0.6, size=3, aes(colour = donor)) +
  #scale_shape_manual(value= 1) +
  scale_colour_manual(values = c("#F8766D","#A3A500","#00BF7D","#A6CEE3","#E76BF3")) +
  ggtitle(paste(" Simulate bulk RNA \n",
                "after scRNA SCT and lognormalize\n", name, "\n",
                "(", num_gene, " x 8000)", sep="")) +
  labs(tag = paste("Donor's AMI: \n", a_x, sep = "")) +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=10),
        plot.tag.position = c(0.9, 0.97))
ggsave(paste(umap_wd, "_Simulate_RNA_lognor_after_SCT_filter_umap.pdf", sep=""),
       p_x, width=6, height=6, units="in", scale=1.7)

end.time <- Sys.time()
time.taken <- end.time - start.time

cat("\n======================================================\n",
    "  Finish the bulk RNA umap: ", date(), "\n",
    "  All code excute: ", time.taken, "min",
    "\n======================================================\n")



