#!/usr/bin/env Rscript
#Author : chien-yu
#Date : 20220503

suppressPackageStartupMessages({
  source('/staging/biology/u8324386/HCC_scRNA/libs.R') # imports
  #source('utils.R') # color definitions and plotting functions
})
start.time <- Sys.time()
option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = FALSE, help = "Simulate RNA RDS data"),
  make_option(c("-o", "--output"), type = "character", default = FALSE, help = "output for bulk RNA-seq direction"))
arguments <- parse_args(OptionParser(option_list=option_list))
args.r <- arguments$rds
args.o <- arguments$output

plot_o <- '/staging/biology/u8324386/HCC_scRNA/2.UMAP_for_scRNA/'

minmax <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

name <- gsub("^.*/", "", args.r)
name <- str_split(name, "\\.")[[1]][1]
name

#將 harmony,SCT 表現量直接轉成 dataframe
##----Harmony 有平移，以防負數----
ref_scrna = readRDS(args.r)
ref_exp_edit <- ref_scrna$loadings %*% ref_scrna$Z_corr
row.names(ref_exp_edit) <- ref_scrna$vargenes[['symbol']]
dim(ref_exp_edit)
exp_min <- which(ref_exp_edit == min(ref_exp_edit), arr.ind = TRUE)
exp_min <- abs(ref_exp_edit[exp_min])
ref_exp_edit[1:5, 1:3]
ref_exp_edit <- ref_exp_edit + exp_min
cat("\ncheck matrix negative value: ",
    sum(apply( ref_exp_edit, 2 , function(x) sum ( x < 0 ) )), "\n")

##--------
dfs1 <- as.data.frame(ref_exp_edit[,1:5000])
dfs2 <- as.data.frame(ref_exp_edit[,5001:10000])
dfs <- cbind(dfs1, dfs2)
row.names(dfs) <- row.names(ref_exp_edit)

dfs1 <- as.data.frame(ref_exp_edit[,10001:15000])
dfs <- cbind(dfs, dfs1)
dfs1 <- as.data.frame(ref_exp_edit[,15001:20000])
dfs <- cbind(dfs, dfs1)

dfs1 <- as.data.frame(ref_exp_edit[,20001:25000])
dfs <- cbind(dfs, dfs1)
dfs1 <- as.data.frame(ref_exp_edit[,25001:30000])
dfs <- cbind(dfs, dfs1)

dfs1 <- as.data.frame(ref_exp_edit[,30001:35000])
dfs <- cbind(dfs, dfs1)
dfs1 <- as.data.frame(ref_exp_edit[,35001:40000])
dfs <- cbind(dfs, dfs1)

dfs1 <- as.data.frame(ref_exp_edit[,40001:45000])
dfs <- cbind(dfs, dfs1)
dfs1 <- as.data.frame(ref_exp_edit[,45001:50000])
dfs <- cbind(dfs, dfs1)

dfs1 <- as.data.frame(ref_exp_edit[,50001:55000])
dfs <- cbind(dfs, dfs1)
dfs1 <- as.data.frame(ref_exp_edit[,55001:60000])
dfs <- cbind(dfs, dfs1)

dfs1 <- as.data.frame(ref_exp_edit[,60001:65000])
dfs <- cbind(dfs, dfs1)
dfs1 <- as.data.frame(ref_exp_edit[,65001:66187])
dfs <- cbind(dfs, dfs1)

rm(dfs1,dfs2)
gc()

#加總基因在每個樣本的表現量
df_sum1 <- as.data.frame(rowSums(dfs[,1:20000]))
df_sum2 <- as.data.frame(rowSums(dfs[,20001:40000]))
df_sum <- cbind(df_sum1,df_sum2)
df_sum1 <- as.data.frame(rowSums(dfs[,40001:60000]))
df_sum <- cbind(df_sum,df_sum1)
df_sum1 <- as.data.frame(rowSums(dfs[,60001:66187]))
df_sum <- cbind(df_sum,df_sum1)
df_sum <- as.data.frame(rowSums(df_sum))

rm(df_sum1,df_sum2)
gc()

#去除表現量太高或太低的基因
q.001 <- quantile(df_sum[,1], 0.01)
q.099 <- quantile(df_sum[,1], 0.99)

a <- which(df_sum[,1] > q.001)
b <- which(df_sum[,1] < q.099)
#聯及
q_filter <- intersect(a, b)

#去除cv太高或太低的基因
coeff_var <- function(x, na.rm=TRUE) {
  return(sd(x)/mean(x))
}

#使用 harmony 用這邊
a <- as.data.frame(t(dfs[1:2000, ]))
co1 <- sapply(a, coeff_var )
a <- as.data.frame(t(dfs[2001:3692, ]))
co2 <- sapply(a, coeff_var )
co <- append(co1,co2)

#做整合
t <- as.array(co)
df_co <- data.frame(t)
co_na <- which(is.na(df_co$CV))
sum(is.na(df_co[,1])) #6574

rm(a,co,co1,co2,t)
gc()


#去除表現量太高或太低的基因
q_symbol <- row.names(dfs)[q_filter]
df_co <- as.data.frame(df_co[q_symbol, ])
row.names(df_co) <- q_symbol
colnames(df_co) <- "CV"
sum(is.na(df_co$CV)) #0

#利用剩下的15602 genes 過濾
# 去除cv太高或太低的基因
cv.001 <- quantile(df_co$CV, 0.01)
cv.099 <- quantile(df_co$CV, 0.99)

#Harmony 過濾
a <- which(df_co[,1] > cv.001)
b <- which(df_co[,1] < cv.099)

cv_filter <- intersect(a, b)
cv_symbol <- row.names(df_co)[cv_filter]
exp_filter <- ref_exp_edit[cv_symbol, ]
dim(exp_filter)

## Harmony儲存
saveRDS(exp_filter, paste(args.o, name, '_filter.rds', sep=''))
cv_symbol <- as.data.frame(cv_symbol)
write.csv(cv_symbol, paste(args.o, name, '_filter_gene.csv', sep=''))
which(row.names(ref_exp_edit) == cv_symbol[1])


####----bar plot for CV-----
library(tidyverse)
max(df_co$CV)
min(df_co$CV)
quantile(df_co$CV, c(0.25,0.5,0.75))

#harmony
ll <- c(0,0.05,0.1,0.15,0.2,0.25)

##計算每區段的CV，基因數量
a = c()
cv_level = c()

for (t in 2:length(ll)) {
  tmp <- df_co %>%
    filter(CV >= ll[t-1], CV < ll[t])
  a <- append(a,nrow(tmp))  
  tmp2 <- paste(ll[t-1], '-', ll[t], sep = '')
  cv_level <- append(cv_level, tmp2)
}
a
cv_table <- as.data.frame(a)

#畫圖
# Basic barplot
cv_table$num <- cv_level

library(ggplot2)
## Hamrony
p<-ggplot(data=cv_table, aes(x=factor(num, level = cv_level), y=a)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = a), vjust = -0.3) +
  labs(title = "Harmony coefficient of variation \n(3618 genes)",
       x = "interval", y="counts") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 30, vjust = 0.5))
p
ggsave(paste(plot_o, "Harmony_CV_bar_plot.pdf", sep=''),
       p, width=6, height=5, units="in", scale=1.5)
