
getwd()
setwd("/Users/me5547lon/Documents/ymwork/4.microarray/20200701")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("geneplotter")
BiocManager::install("gplots")
BiocManager::install("limma")
BiocManager::install("simpleaffy")
#Plot package for microarray data
#對於microarray資料畫圖
library(geneplotter)

#Fill the label for the figuer and table
#填入圖表上各項資訊的標記
library(gplots)

#Calculate the liner microarray data gene experssion difference
#計算線性microarray data基因表現量的不同
library(limma)

#Analyze package for affymatrix
#單一分析affymatrix的函式庫
library(simpleaffy)

#Rename the file name automatically
#Define a functoin that substitute the space with underline
rename.files <- function(x){
  file.rename(x, gsub(" ","_",x))
}
#collect the file name automatically
cel.files <- list.files(pattern="CEL$") #find the working directory with CEL file.
sapply(cel.files,rename.files) #apply the function to the list.
#Get Sample_id
tmp_sample.names <- {}
split_file_name <- {}
for (i in c(1:length(cel.files))){
  split_file_name <- strsplit(cel.files[i], split = "_")
  if (length(split_file_name[[1]])<4){
    tmp_sample.names <- c(tmp_sample.names, split_file_name[[1]][1])
  }else{
    tmp_sample.names <- c(tmp_sample.names, paste(split_file_name[[1]][1],
                                                  split_file_name[[1]][2],
                                                  sep = "_"))
  }
}
tmp_sample.names <- gsub("_L-T", "LT", tmp_sample.names)
tmp_sample.names <- gsub("_L-N", "LN", tmp_sample.names)
tmp_sample.names <- gsub("m_Liv", "mLiv", tmp_sample.names)
tmp_sample.names <- gsub("m_LT", "mLT", tmp_sample.names)
tmp_sample.names <- gsub("m_LN", "mLN", tmp_sample.names)
sample.names <- {}
split_file_name <- {}
for (i in c(1:length(tmp_sample.names))){
  split_file_name <- strsplit(tmp_sample.names[i], split = "-")
  if (length(split_file_name[[1]])<4){
    sample.names <- c(sample.names, split_file_name[[1]][3])
  }else{
    sample.names <- c(sample.names, paste(split_file_name[[1]][3],
                                          split_file_name[[1]][4],
                                          sep = "_"))
  }
}

#---**---8/27

##Automatically give B6, 122+/-, HBV_Tg, hybrid
mice_table <- read.csv("mice_group.csv", stringsAsFactors = F)
tmp_cel.group <- mice_table[,5][which(mice_table[,2] %in% sample.names)]
#replace HBV Tg -> nonhybrid, 122+/- -> miR122, B6 -> contorl
cel.group <-gsub("HBV Tg","nonhybrid",gsub("122\\+\\/\\-","miR122",gsub("B6","control",tmp_cel.group)))
cel.batch <- vector(mode="character",length=length(cel.files))
rm(i, split_file_name, tmp_sample.names, tmp_cel.group, mice_table)
## bind the cel.file, and three vector above by row.
cel.des.tab <- cbind(cel.files,sample.names,cel.group,cel.batch)
## give the matrix the colname.
colnames(cel.des.tab) <- c("FileName","SampleName","Group","Batch")
## output csv file
write.csv(cel.des.tab,"CEL.des.tab.csv",row.names=F,quote=F)

#Show the package version 
package_name <- c("geneplotter", "limma", "simpleaffy", "gplots")
each_packgae_version <- matrix(NA, nrow = length(package_name), ncol = 2 )
colnames(each_packgae_version) <- c("Package name", "Version")
for (i in 1:length(package_name)){
  each_packgae_version[i, 1] <- package_name[i]
  each_packgae_version[i, 2] <- package.version(package_name[i])
}
write.csv(each_packgae_version, "Package_version.csv")
rm (cel.des.tab, each_packgae_version, cel.batch, cel.files, cel.group, sample.names)

#Import CEL description table for working directory.
#And put each column to individual variable
cel.des.tab <- read.csv("CEL.des.tab.csv")
cel.files <- as.character(cel.des.tab[,1])
sample.names <- as.character(cel.des.tab[,2])
cel.group <- as.character(cel.des.tab[,3])
cel.batch <- as.numeric(cel.des.tab[,4])
#Set the color
#Assign different colors to samples of different groups
#unique:remove the same value, sort:sort by value
group.names <- sort(unique(cel.group))
#number of group
ng <- length(group.names)
#Give the color to each group
#Extract color from color name which R knows about, but we need to remove the similar color and the color 
#which is hard to identify from white background.
clist <- colors()
clist <- clist[-c(grep("gray", clist), 
                  grep("grey", clist),
                  grep("white", clist),
                  grep("black", clist),
                  grep("1", clist),
                  grep("2", clist), 
                  grep("3", clist),
                  grep("4", clist),
                  grep("^light", clist),
                  grep("^dark", clist),
                  grep("^medium", clist),
                  grep("^corn", clist),
                  grep("^pale", clist),
                  grep("mintcream", clist))]
group.color.palettes <- clist[sample(5:length(clist), ng)]
#Set label's color->label was taken from cel.group
sample.col <- as.character(factor(cel.group,labels=group.color.palettes))
#Assign different colors to samples of different batches
batch.names <- sort(unique(cel.batch))
nb <- length(batch.names)
if(nb >1 ){
  batch.color.palettes <- clist[sample(length(clist), (ng-2))]
  batch.col <- as.character(factor(cel.batch,labels = batch.color.palettes))
}
#Generate the required file 'covdesc'
#Comebind the cel.file and cel.group by column.
covdesc <- cbind(cel.files,cel.group)
#Change covdesc's second column name to "Treatment
colnames(covdesc) <- c("","Treatment")
#write.table(x,filename)
#sep=set separator string. quote.character or factor columns will be surrounded by double quotes.(logic value)
write.table(covdesc,"covdesc",row.names=F,col.names=T,sep=" ",quote=F)



#---**----



#Import data:use read.affy to read the CEL file from covdesc file list, and read in a special class.
raw.data <- read.affy("covdesc") # affy can read CEL file
#QC plots
#Normalized the affy data by mas5 algorithm
x.mas <- call.exprs(raw.data,"mas5")
#Do the affymatrix QC
#QC the normalized data(x.mas) by qc function. qc(unnormalised, ...)
qcs <- qc(raw.data,x.mas)
'''
more info about qcs

data(qcs)
ratios(qcs)
avbg(qcs)
maxbg(qcs)
minbg(qcs)
spikeInProbes(qcs)
qcProbes(qcs)
percent.present(qcs)
plot(qcs)
sfs(qcs)
target(qcs)
ratios(qcs)
'''
#Generate an output PDF file
#Filename, width 11 inch, height 7inch.If the data more than 10, the plot margin will be too large.
#Generate the pdf file bit without figure
pdf("1.Affymetrix_qc.pdf",width=11,height=(7+max(0,(length(cel.files)-10))*7/10))# height setting 
#plot
plot(qcs)
#dev.off write the figure to the file and save
dev.off()

#Unnromalized box plot
#Unnromalized data between RMA and MAS5 should be the same. So use the function in RMA to get unnormalized data
raw.exprs <- exprs(justRMA(filenames=cel.files,normalize=FALSE,background=FALSE))
#Extract the probe name
probe.names <- rownames(raw.exprs)
#Use the sample name to replace the column name
colnames(raw.exprs) <- sample.names
#Generate a boxplot (before MAS Normalization)
#Set the plot width
plot.width <- 7+(max(0,length(cel.files)-20)*7/20)
pdf("2.Boxplot Plot before MAS5 Normalization.pdf",width=plot.width)
#par=set the figure's parameter mar=c(bottom, left, top, right) 
#xpd=Show the icon or not
par(mar=c(10, 4, 4, 10)+0.1,xpd=T)
#frame:dd the outer line on boxplot
#las:Set the title direction, 0=parallel with the figure, 1=vertical with the figure
boxplot(raw.exprs,las=2,col=sample.col,main="Boxplot (before MAS Normalization)",frame=F)
#lengend->add on icon ncol?Gnumber of mas.exprs max?Gmaxium value of mas.exprs
#(ncol,max)=(x,y), pch=set the icon use in the figure
legend(ncol(raw.exprs)+1,max(raw.exprs),group.names,pch=15,col=group.color.palettes)
#if the batch over 2 then go this for loop to process
if(nb>1){
  par(mar=c(10, 4, 4, 10)+0.1,xpd=T)
  boxplot(raw.exprs,las=2,col=batch.col,main="Boxplot (before MAS Normalization)",frame=F)
  legend(ncol(raw.exprs)*1.2,max(raw.exprs),batch.names,pch=15,col=batch.color.palettes)
}
dev.off()

#Extract the exprs from mas5 normalized data
mas.exprs <- exprs(x.mas)
#Extract the probe name
probe.names <- rownames(mas.exprs)
#Use the sample name to replace the 
colnames(mas.exprs) <- sample.names
#Generate a boxplot (after MAS Normalization)
#Set the plot width
plot.width <- 7+(max(0,length(cel.files)-20)*7/20)
pdf("3.Boxplot Plot after MAS5 Normalization.pdf",width=plot.width)
#par=set the figure's parameter mar=c(bottom, left, top, right) 
#xpd=Show the icon or not
par(mar=c(10, 4, 4, 10)+0.1,xpd=T)
#frame->Add the outer line on boxplot
#las->Set the title direction, 0:parallel with the figure, 1:vertical with the figure
boxplot(mas.exprs,las=2,col=sample.col,main="Boxplot (after MAS Normalization)",frame=F)
#lengend?Gadd on icon ncol=number of mas.exprs max=maxium value of mas.exprs
#(ncol,max)=(x,y), pch=Gset the icon use in the figure
legend(ncol(mas.exprs)+1,max(mas.exprs),group.names,pch=15,col=group.color.palettes)
##if the batch over 2 then go this for loop to process
if(nb>1){
  par(mar=c(10, 4, 4, 10)+0.1,xpd=T)
  boxplot(mas.exprs,las=2,col=batch.col,main="Boxplot (after MAS Normalization)",frame=F)
  legend(ncol(mas.exprs)*1.2,max(mas.exprs),batch.names,pch=15,col=batch.color.palettes)
}
dev.off()

#Check Affymetrix spike-in control probes
#Find the probe is start with AFF, use re to find the row number.
#Here use the method2 data to plot the heatmap
affy.spike.in.ii<-grep("^AFF",probe.names)
#Extract the spike-in control row data from raw.exprs
affy.spike.in.mat.raw<- raw.exprs[affy.spike.in.ii,]
#Generate a heatmap of Affymetrix spike-in control probes (before MAS Normalization)
pdf("4.Heatmap of Affymetrix Spike-in Control Probes (before MAS5 Normalization).pdf",width=14.5,height=14.5)
#col：set the color,  dendrogram：select for 'none', 'row', 'column' or 'both'
#Rowv：do dendorgram for row Colv?Gdo dendorgram for row.
#ColSideColors：set the column side color in the front code had set the color list.
#lhei：heatmap's height & width
#trace:character string indicating whether a solid "trace" line should be drawn across 'row's or down 'column's, 'both' or 'none'.
#Keysize：numeric value indicating the size of the key
#margins：numeric vector of length 2 containing the margins (see par(mar= *)) for column and row names, respectively.
## cexRow, cexCo：used as cex.axis in for the row or column axis labeling. 
heatmap.2(affy.spike.in.mat.raw,col=dChip.colors(19),
          dendrogram="both",Rowv=T,Colv=T,ColSideColors=sample.col,lhei=c(5,15),
          trace="none",keysize=1.2,margins=c(8,12),cexRow=0.7,cexCol=1,
          main="Heatmap of Affymetrix Spike-in Control Probes")
legend(0.85,0.99, legend = group.names,fill = group.color.palettes ,cex = 0.8)
dev.off()

#Check Affymetrix spike-in control probes
## 輸出probe.name內開頭為AFF的位置取出來
affy.spike.in.ii<-grep("^AFF",probe.names)
#Extract the spike-in control row data from mas.exprs
affy.spike.in.mat<- mas.exprs[affy.spike.in.ii,]
#Generate a heatmap of Affymetrix spike-in control probes (after MAS Normalization)
pdf("5.Heatmap of Affymetrix Spike-in Control Probes (after MAS5 Normalization).pdf",width=14.5,height=14.5)
#col：set the color,  dendrogram：select for 'none', 'row', 'column' or 'both'
#Rowv：do dendorgram for row Colv?Gdo dendorgram for row.
#ColSideColors：set the column side color in the front code had set the color list.
##https://www.plob.org/article/10045.html
heatmap.2(affy.spike.in.mat,col=dChip.colors(19),
          dendrogram="both",Rowv=T,Colv=T,ColSideColors=sample.col,lhei=c(5,15),
          trace="none",keysize=1.2,margins=c(8,12),cexRow=0.7,cexCol=1,
          main="Heatmap of Affymetrix Spike-in Control Probes")
legend(0.85,1, legend = group.names,fill = group.color.palettes ,cex = 0.7)
dev.off()

## Remove Affymetrix QC related variables
rm(raw.data, x.mas, qcs, mas.exprs, raw)
gc()


# MDS without Normalization
unnor.mds <- justRMA(filenames=cel.files, normalize=FALSE, background=FALSE)
#Extract the value 
mds.exprs <- exprs(unnor.mds)
#Use the sample name to replace the colname
colnames(mds.exprs)<-sample.names
#Calculate the distance between the matrix
#t Transpose, rma.exprs. Because R calculation need by row
#distance method ? euclidean?
d<-dist(t(mds.exprs))
#Use MDS to calculate the distance between sample and sample
#eigenvalue, add ?? MDS plot? Su
dS<-cmdscale(d)
t(mds.exprs)
#Generate a MDS plot
pdf("6.MDS Plot before RMA Normalization.pdf")
par(mar=c(5, 4, 4, 6)+0.1, xpd=T)
plot(dS,xlab="Dimension 1",ylab="Dimension 2",
     pch=20,col=sample.col,main="MDS Plot (before RMA Normalization)",
     xlim = c(-200, 200), ylim = c(-200,200))
#axis(1,seq(-100, 100, 50), las = 1)
#axis(2,seq(-100, 100, 50), las = 1)
legend("topright",group.names,pch=20,col=group.color.palettes,cex =1)

if(nb>1){
#if batch more then 1
  par(mar=c(5, 4, 4, 6)+0.1, xpd=T)
  plot(dS,xlab="Dimension 1",ylab="Dimension 2",
       pch=20,col=sample.col,main="MDS Plot (before RMA Normalization)",
       xlim = c(-200, 200), ylim = c(-200, 200))
  legend(max(dS[,1])*1.1,max(dS[,2]),batch.names,pch=20,col=batch.color.palettes)
}
dev.off()

#plot the unnormalized boxplot  (RMA)
#Generate a boxplot
plot.width<-7+(max(0,length(cel.files)-20)*7/20)
pdf("7.Boxplot Plot before RMA Normalization.pdf",width=plot.width)
par(mar=c(10, 4, 4, 10)+0.1,xpd=T)
boxplot(mds.exprs,las=2,col=sample.col,main="Boxplot (before RMA Normalization)",frame=F)
legend(ncol(mds.exprs)+1,max(mds.exprs),group.names,pch=15,col=group.color.palettes)
#if batch more then 1
if(nb>1){
  par(mar=c(10, 4, 4, 10)+0.1,xpd=T)
  boxplot(mds.exprs,las=2,col=batch.col,main="Boxplot (before RMA Normalization)",frame=F)
  legend(ncol(mds.exprs)+1,max(mds.exprs),batch.names,pch=15,col=batch.color.palettes)
}
dev.off()

#RMA EXPRS Export
#Background defult? or set ?
rma<-justRMA(filenames=cel.files,normalize=TRUE)
#Extract the value 
rma.exprs<-exprs(rma)
#Extract the probe name
probe.names<-rownames(rma.exprs)
#Use the sample name to replace the colname
colnames(rma.exprs)<-sample.names
#Combine the probe name and the colmane(had replaced by sample )
out.mat<-cbind(probe.names,rma.exprs)
#Write the name to the figure and make the figure format
write.table(out.mat,"RMA_normalized_expression_matrix.txt",row.names=F,col.names=T,sep="\t",quote=F)
rm(out.mat)

#Calculate the distance between the matrix
d<-dist(t(rma.exprs))
#Use MDS to calculate the distance between sample and sample
dS<-cmdscale(d)
#Generate a Zoom in MDS plot
pdf("8.MDS Plot after RMA Normalization.pdf")
par(mar=c(5, 4, 4, 6)+0.1, xpd=F)
plot(dS,xlab="Dimension 1",ylab="Dimension 2",
     pch=20,col=sample.col,main="MDS Plot (after RMA Normalization)",
     xlim = c(-100, 100), ylim = c(-100, 100), xaxt = "n", yaxt = "n")
axis(1,seq(-100, 100, 50), seq(-100, 100, 50), las = 1)
axis(2,seq(-100, 100, 50), seq(-100, 100, 50), las = 1)
legend("topright",group.names,pch=20,col=group.color.palettes,cex =1)
#if batch more then 1
if(nb>1){
  par(mar=c(5, 4, 4, 6)+0.1, xpd=F)
  plot(dS,xlab="Dimension 1",ylab="Dimension 2",pch=20,
       col=batch.col,main="MDS Plot (after RMA Normalization)")
  legend(max(dS[,1])*1.1,max(dS[,2]),batch.names,pch=20,col=batch.color.palettes)
}
dev.off()

#畫3D圖
#library(rgl)
#dS3D <- cmdscale(d, k = 3)
#plot3d(dS3D,xlab="Dimension 1",ylab="Dimension 2",zlab = "Dimension 3",
#              type = "s",col = sample.col, size = 1,
#              main="MDS Plot (after RMA Normalization)")
#legend(group.names,pch=20,col=group.color.palettes)
#dev.off()

#boxplot
#Generate a boxplot
plot.width<-7+(max(0,length(cel.files)-20)*7/20)
pdf("9.Boxplot Plot after RMA Normalization.pdf",width=plot.width)
par(mar=c(10, 4, 4, 10)+0.1,xpd=T)
boxplot(rma.exprs,las=2,col=sample.col,main="Boxplot (after RMA Normalization)",frame=F)
legend(ncol(rma.exprs)+1,max(rma.exprs),group.names,pch=15,col=group.color.palettes)
#if batch more then 1
if(nb>1){
  par(mar=c(10, 4, 4, 10)+0.1,xpd=T)
  boxplot(rma.exprs,las=2,col=batch.col,main="Boxplot (after RMA Normalization)",frame=F)
  legend(ncol(rma.exprs)+1,max(rma.exprs),batch.names,pch=15,col=batch.color.palettes)
}
dev.off()

#Heatmap before RMA(Quality control) 
## Check Affymetrix spike-in control probes
## 輸出probe.name內開頭為AFF的位置取出來
unnor.rma<-justRMA(filenames=cel.files,normalize=FALSE,background=FALSE)
#Extract the value 
unnor.rma.exprs<-exprs(unnor.rma)
#Use the sample name to replace the colname
colnames(unnor.rma.exprs)<-sample.names
affy.spike.in.ii<-grep("^AFF",probe.names)
## 將mas.exprs與affy.spike.in.ii相同的列取出，並列出當中的row
affy.spike.in.mat.unnor<-unnor.rma.exprs[affy.spike.in.ii,]
## Generate a heatmap of Affymetrix spike-in control probes
pdf("10.Heatmap of Affymetrix Spike-in Control Probes (before RMA Normalization).pdf",width=14.5,height=14.5)
heatmap.2(affy.spike.in.mat.unnor,col=dChip.colors(19),
          dendrogram="both",Rowv=T,Colv=T,ColSideColors=sample.col,lhei=c(5,15),
          trace="none",keysize=1.2,margins=c(8,12),cexRow=0.7,cexCol=1,
          main="Heatmap of Affymetrix Spike-in Control Probes(before RMA Normalization)")
legend(0.85,1, legend = group.names,fill = group.color.palettes ,cex = 0.7)
dev.off()

#Heatmap after RMA(Quality control) 
## Check Affymetrix spike-in control probes
## 輸出probe.name內開頭為AFF的位置取出來
affy.spike.in.ii<-grep("^AFF",probe.names)
## 將mas.exprs與affy.spike.in.ii相同的列取出，並列出當中的row
affy.spike.in.mat<-rma.exprs[affy.spike.in.ii,]
## Generate a heatmap of Affymetrix spike-in control probes
pdf("11.Heatmap of Affymetrix Spike-in Control Probes (after RMA Normalization).pdf",width=14.5,height=14.5)
heatmap.2(affy.spike.in.mat,col=dChip.colors(19),
          dendrogram="both",Rowv=T,Colv=T,ColSideColors=sample.col,lhei=c(5,15),
          trace="none",keysize=1.2,margins=c(8,12),cexRow=0.7,cexCol=1,
          main="Heatmap of Affymetrix Spike-in Control Probes(after RMA Normalization)")
legend(0.85,1, legend = group.names,fill = group.color.palettes ,cex = 0.7)
dev.off()

#3.3.2	Limma
#3.3.4.1 Compute
#### Identify differentially expressed probes by limma ####

# Define the design matrix
cel.group
## 將group取出，會將相同的分類用標示
treatment<-as.factor(cel.group)
#establish a unique design matrix
#Why use = not the <-?
#Design matrix?
design = model.matrix(~0 + treatment)
#change column name 
colnames(design) <- gsub("treatment","",colnames(design))
#fit the linear model
#use linear model to calculate
fit = lmFit(rma, design)
#改成自動化
#To compare between the conditions
ng<-length(group.names)
## ng的內容分割成兩個row
combn.mat<-combn(ng,2)

#make contrasts matrix
## 
contrasts.maker<-function(idx,group.names){
  return(paste(group.names[as.numeric(idx)],collapse="-"))
}
## 
#my.contrasts <- apply(combn.mat,2,contrasts.maker,group.names)
my.contrasts <- c("hybrid-control")
#levels :the parameter names as column names.
contrast.matrix<-makeContrasts(contrasts=my.contrasts,levels=design)
# Estimate those contrasts using contrasts.fit from limma
#Contrasts fit :Given a linear model fit to microarray data, compute estimated coefficients and standard errors 
      #for a given set of contrasts.
fit2 = contrasts.fit(fit, contrast.matrix)

# Obtain the moderated statistics using eBayes from limma
fit2 = eBayes(fit2)


z.transformation <- function(x){
  return((x-mean(x))/sd(x))
}
#3.3.4.2 P-values (0.01, 0.05, 0.1, 0.3, 0.5)
## 設定要找出的p-value
adj.p.threshold <- c(0.01, 0.05, 0.1)
if(ncol(contrast.matrix)<=5){
  pdf("12.Result of limma - Venn Diagram.pdf")
}
for(i in 1:length(adj.p.threshold)){
  DE = decideTests(fit2,method="separate",adjust.method="BH", p.value = adj.p.threshold[i])    
  print(adj.p.threshold[i])
  ##  vennDiagram圖
  if(ncol(contrast.matrix)<=5){
    vennDiagram(DE,main=paste("Adjust p-value Threshold =", adj.p.threshold[i]))    
  }
  #Extract a table of the top-ranked genes from a linear model fit.
  tT<-topTable(fit2,adjust="BH",sort.by="B",p.value=adj.p.threshold[i],number=45101)
  
  if(nrow(tT)>0){
    # Write out the summary table
    my.probe.names <- rownames(tT)
    np <- length(my.probe.names)
    my.row.cex<-max(0.3,0.7-round(max(0,(np-50)/50),0)*0.07)
    tT <- cbind(my.probe.names,tT)
    colnames(tT)[1] <- c("ProbeNames")
    write.csv(tT,paste("13.Result of limma - TopTable - Adjust p-value Threshold = ",
                       adj.p.threshold[i],".csv",sep=""),row.names=F)	
    
    # Generate a heatmap with these DE genes
    my.pii <- which(probe.names %in% my.probe.names)
    exprs.mat <- rma.exprs[my.pii,]
    z.exprs.mat <- t(apply(exprs.mat,1,z.transformation))	
    
    pdf(paste("14.Result of limma - Heatmap of the DE Probes (Adjust p-value Threshold = ",adj.p.threshold[i],").pdf",sep=""),
        width=14,height=14)
    # original scale
    heatmap.2(exprs.mat,col=dChip.colors(19),
              dendrogram="both",Rowv=T,Colv=T,ColSideColors=sample.col,lhei=c(3,15),
              trace="none",keysize=1.2,margins=c(8,12),cexRow=my.row.cex,cexCol=1,
              main=paste("Heatmap of the DE Probes (Adjust p-value Threshold = ",adj.p.threshold[i],")\nNo. of Probes = ",
                         length(my.probe.names),sep=""))
    legend(0.85,1, legend = group.names,fill = group.color.palettes ,cex = 0.7)
    # with z-score transformation
    heatmap.2(z.exprs.mat,col=dChip.colors(19),
              dendrogram="both",Rowv=T,Colv=T,ColSideColors=sample.col,lhei=c(3,15),
              trace="none",keysize=1.2,margins=c(8,12),cexRow=my.row.cex,cexCol=1,
              main=paste("Heatmap of the DE Probes (Adjust p-value Threshold = ",
                         adj.p.threshold[i],")\nz-score transformation by probe sets\nNo. of Probes = "
                         ,length(my.probe.names),sep=""))
    legend(0.85,1, legend = group.names,fill = group.color.palettes ,cex = 0.7)
    dev.off()
    
    rm(my.probe.names,my.pii,exprs.mat,z.exprs.mat,np,my.row.cex)
  }
  rm(DE,tT)
}

if(ncol(contrast.matrix)<=5){
  dev.off()
}


#---------------------
## 當teatment > 2才會執行以下程式，也是p-values的圖形呈現
if(ncol(contrast.matrix) >=2){
  adj.p.threshold <- c(0.01, 0.05, 0.1)
  for(i in 1:length(adj.p.threshold)){
    tT.matrix<-NULL
    for(j in 1:ncol(contrast.matrix)){  
      tT<-topTable(fit2,j,adjust="BH",sort.by="B",p.value=adj.p.threshold[i],number=45101)
      if(nrow(tT)>0){
        my.probe.names<-rownames(tT)
        my.contrast<-rep(colnames(contrast.matrix)[j],length(my.probe.names))
        tT<-cbind(my.probe.names,my.contrast,tT)
        colnames(tT)[1:2]<-c("ProbeNames","Contrast")
        tT.matrix<-rbind(tT.matrix,as.matrix(tT))
        rm(my.contrast,tT)
      }
    }
    write.csv(tT.matrix,paste("15.Result of limma - TopTable of different contrasts - Adjust p-value Threshold = ",adj.p.threshold[i],".csv",sep=""), row.names=F)	
  
    my.probe.names<-rownames(tT.matrix)
    np<-length(my.probe.names)
    my.row.cex<-max(0.3,0.7-round(max(0,(np-50)/50),0)*0.07)
  
    # Generate a heatmap with these DE genes
    my.pii<-which(probe.names %in% my.probe.names)
    exprs.mat<-rma.exprs[my.pii,]
    z.exprs.mat<-t(apply(exprs.mat,1,z.transformation))	
  
    pdf(paste("16.Result of limma - Heatmap of the DE Probes of different contrasts (Adjust p-value Threshold = ",adj.p.threshold[i],").pdf",sep=""),width=14,height=21)
    # original scale
    heatmap.2(exprs.mat,col=dChip.colors(19),
              dendrogram="both",Rowv=T,Colv=T,ColSideColors=sample.col,lhei=c(3,15),
              trace="none",keysize=1.2,margins=c(8,12),cexRow=my.row.cex,cexCol=1,
              main=paste("Heatmap of the DE Probes (Adjust p-value Threshold = ",adj.p.threshold[i],")\nNo. of Probes = ",length(my.probe.names),sep=""))
    legend(0.85,1, legend = group.names,fill = group.color.palettes ,cex = 0.7)
    # with z-score transformation
    heatmap.2(z.exprs.mat,col=dChip.colors(19),
              dendrogram="both",Rowv=T,Colv=T,ColSideColors=sample.col,lhei=c(3,15),
              trace="none",keysize=1.2,margins=c(8,12),cexRow=my.row.cex,cexCol=1,
              main=paste("Heatmap of the DE Probes (Adjust p-value Threshold = ",adj.p.threshold[i],")\nz-score transformation by probe sets\nNo. of Probes = ",length(my.probe.names),sep=""))
    legend(0.85,1, legend = group.names,fill = group.color.palettes ,cex = 0.7)
    dev.off()
  rm(my.probe.names,my.pii,exprs.mat,z.exprs.mat,np,my.row.cex,tT.matrix)
  }
}


tT.matrix<-NULL
for(i in 1:ncol(contrast.matrix)){        
  tT<-topTable(fit2,i,adjust="BH",sort.by="B",p.value=1,number=45101)
  if(nrow(tT)>0){
    my.probe.names<-rownames(tT)
    my.contrast<-rep(colnames(contrast.matrix)[i],length(my.probe.names))
    tT<-cbind(my.probe.names,my.contrast,tT)
    colnames(tT)[1:2]<-c("ProbeNames","Contrast")
    tT.matrix<-rbind(tT.matrix,as.matrix(tT))
    rm(my.contrast,tT)
  }
}
write.csv(tT.matrix,paste("15.Result of limma - TopTable of different contrasts - Without Adjust p-value Threshold",".csv",sep=""), row.names=F)	

## 將整支程式執行的時間輸出
end.time<-Sys.time()
runTime<-difftime(end.time,start.time,units="secs")
print(runTime)
# Record the processing time in seconds
write.table(runTime,"Time.log.txt",row.names=F,col.names=F,sep="\t",quote=F)

#semiscatterplot
file_num = round((nrow(rma.exprs)/1000))
for (i in (1:file_num)){
  if (file_num > i ){
    pdf(paste(15+i,".Result of limma - Scatter plot of the DE Gene ", (i-1)*1000+1 ,"-", i*1000, ".pdf",sep=""))
    a = (i-1)*1000+1
    b = i*1000
    for (j in c(a:b)){
      no.row <- which(probe.names %in% dt_table[j,1])
      p = dt_table[j, 6]
      plot(c(1:ncol(rma.exprs)), (rma.exprs[no.row, ])[order(treatment)]
           , main =paste(dt_table[j,1], "\nAdjust-p.value: ", p,"\nLogFC: ", dt_table[j,2]),
           xlab = "Sample",
           ylab = "Expression",
           pch = 20,col = sample.col[order(treatment)])
      legend("topright", legend = c("NT","TU"),fill = sample.col ,cex = 1)
    }
    dev.off()
  }
  else{
    pdf(paste(15+i,".Result of limma - Scatter plot of the DE Gene ", (i-1)*1000+1 ,"-", i*1000, ".pdf",sep=""))
    a = (i-1)*1000+1
    for (j in c(a:b)){
      no.row <- which(probe.names %in% dt_table[j,1])
      p = dt_table[j, 6]
      plot(c(1:ncol(rma.exprs)), (rma.exprs[no.row, ])[order(treatment)]
           , main =paste(dt_table[j,1], "\nAdjust-p.value: ", p,"\nLogFC: ", dt_table[j,2]),
           xlab = "Sample",
           ylab = "Expression",
           pch = 20,col = sample.col[order(treatment)])
      legend("topright" , legend = c("NT","TU"),fill = sample.col ,cex = 1)
    }
    dev.off()
  }
}
