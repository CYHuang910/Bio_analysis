
setwd("/Users/qianyu/Documents/ymwork/3.project/1.mRNA/20201007_Cryocord/edge_R/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

# Load package
library("limma")
library("edgeR")

# Load read_count.txt data
read_count <- read.delim("read_count.txt", header = T)
colnames(read_count) <- c("", "N2", "N3", "N1", "O2", "O3", "O1")
#remove column 1-3
read_matrix <- read_count[-c(1, 2, 3), ]
read_matrix <- read_matrix[, c(1, 4, 2, 3, 7, 5, 6)]
rownames(read_matrix) = seq(length=nrow(read_matrix))
as.matrix(read_matrix)

#DGEList
y <- DGEList(counts = read_matrix[,2:7], genes = read_matrix[,1])
#TMM normalization
y <- calcNormFactors(y, method = "TMM")
y$samples

#define our design matrix based on the experimental design.
condition <- factor(c(rep("N", 3), rep("O", 3)))
ID <- factor(c(seq(1, 3), seq(1, 3)))
data.frame(sample = colnames(read_matrix[,-c(1)]), condition, ID)

#############--- pair t test ---#################
design_matrix <- model.matrix(~ID+condition)
design_condition <- model.matrix(~0 + condition)
rownames(design_matrix) <- colnames(y)
design_matrix

#plot MDS
pdf('0.sample_plot_MDS.pdf')
target <- data.frame(sample = colnames(read_matrix[,-c(1)]), condition, ID)
group <- factor(target$sample)
points <- c(0,1,2,15,16,17)
colors <- rep(c("blue", "darkgreen", "red"), 2)
plotMDS(y, col=colors[group],pch = points[group])
legend("topright", legend=levels(group),pch=points, col=colors, ncol=2)
dev.off()

#estimate the NB dispersion for the dataset
y <- estimateDisp(y, design_matrix, robust = TRUE)
y$common.dispersion

#Differential expression
fit <- glmFit(y, design_matrix)
lrt <- glmLRT(fit)
#Conduct likelihood ratio tests for two group differences and show the top genes:
topTags(lrt)

#export DEG data
ouput_DEG <- topTags(lrt, n = nrow(read_matrix))
DEG_table <- ouput_DEG$table

#filter the DEG data
DEG_filter_fdr <- DEG_table[DEG_table$FDR < 0.05, ]
write.csv(DEG_filter_fdr, "1.pair_t_DEG_filter_0.5.csv")


#############--- biomaRt ---#################
getwd()
library(biomaRt)

dt_mart <- useMart("ensembl")
#list database 
list_database <- listMarts(dt_mart)
# list the available datasets in this Mart
list_dataset <- listDatasets(dt_mart)
#list the version of ensembl
list_version <- listEnsemblArchives()

mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "https://may2021.archive.ensembl.org")
#select output content
list_ouput_col <- listAttributes(mart)
list_filter <- listFilters(mart)

#Load DEG data filter < 0.5
deg_data <- read.csv("1.pair_t_DEG_filter_0.5.csv")
deg_data <- deg_data[,-c(1)]
ensembl_id <- deg_data[[1]]
head(ensembl_id)

#serch the mart , if the ensembl gene id is same
mart_result <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                filters = "ensembl_gene_id",
                values = ensembl_id,
                mart = mart)

nrow(mart_result)
#check biomart con't find gene 
nrow(mart_result[which(mart_result[,2] == ""),])

#merge the DEG matrix and biomaRt result by ensembl id
#when merge, column's name must be same 
colnames(deg_data)[1] <- "ensembl_gene_id"

# all=TRUE 兩個資料取聯集
combine_table <- merge(x = deg_data, y = mart_result , by = "ensembl_gene_id", all = TRUE)

write.csv(combine_table, "2.combine_pair_t_DEG_and_biomart.csv")


#############--- Compare N,O ---#################
design_condition <- model.matrix(~0 + condition)
rownames(design_condition) <- colnames(y)
design_condition

#estimate the NB dispersion for the dataset
y_con <- estimateDisp(y, design_condition, robust = TRUE)
y_con$common.dispersion

#Differential expression
fit_con <- glmFit(y_con, design_condition)
lrt_con <- glmLRT(fit_con)
#Conduct likelihood ratio tests for two group differences and show the top genes:
topTags(lrt_con)

#export DEG data
ouput_DEG_con <- topTags(lrt_con, n = nrow(read_matrix))
DEG_table_con <- ouput_DEG_con$table

#filter the DEG data
DEG_filter_fdr_con <- DEG_table_con[DEG_table_con$FDR < 0.05, ]
write.csv(DEG_filter_fdr_con, "3.compare_N-O_DEG_filter_0.5.csv")

#############--- biomaRt ---#################
deg_con_data <- read.csv("3.compare_N-O_DEG_filter_0.5.csv")
deg_con_data <- deg_con_data[,-c(1)]
ensembl_id_con <- deg_con_data[[1]]
head(ensembl_id)

#serch the mart , if the ensembl gene id is same
mart_result_con <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                     filters = "ensembl_gene_id",
                     values = ensembl_id_con,
                     mart = mart)

nrow(mart_result_con)
#check biomart con't find gene 
nrow(mart_result_con[which(mart_result_con[,2] == ""),])

#merge the DEG matrix and biomaRt result by ensembl id
#when merge, column's name must be same 
colnames(deg_con_data)[1] <- "ensembl_gene_id"

# all=TRUE 兩個資料取聯集
combine_table_con <- merge(x = deg_con_data, y = mart_result_con , by = "ensembl_gene_id", all = TRUE)

write.csv(combine_table_con, "4.combine_N-O_DEG_and_biomart.csv")








