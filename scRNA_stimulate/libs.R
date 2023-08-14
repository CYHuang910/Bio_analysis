# Analysis
library(optparse)
library(symphony)
library(harmony)
library(singlecellmethods)
library(irlba)
library(tidyverse)
library(data.table)
library(matrixStats)
library(Matrix)
library(plyr)
library(dplyr)
library(Seurat)

#clustering
library(factoextra)
library(cluster)
library(aricode)

# Plotting
library(umap)
library(ggplot2)
library(ggthemes)
library(ggrastr)
library(RColorBrewer)
library(patchwork)
library(ggpubr)

fig.size <- function (height, width) {
  options(repr.plot.height = height, repr.plot.width = width)
}
