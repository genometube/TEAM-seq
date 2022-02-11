########## Load libraries ########## 
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(reshape2)
setwd("C:/Users/fangyitong/Desktop/2021Research/LIANTI/tables")

########## 12.15 human_mouse_collision ########## 
df <- readxl::read_excel("scLIANTI_hg_mm.xlsx")
colnames(df) <- c("barcode","demultiplexed_reads","hg_or_mm","unique_hg","unique_mm","unique_mapping")
ggplot(df,mapping=aes(x=unique_mm/1000000,y=unique_hg/1000000,group=hg_or_mm))+
  geom_point(aes(shape=hg_or_mm,color=hg_or_mm,size=hg_or_mm)) +
  scale_size_manual(values=c(1.5,1.5,1.5)) +
  xlab("Mouse unique aligned reads (in million)") +
  ylab("Human unique aligned reads (in million)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_blank(),
        axis.line=element_line(colour="black"))

########## 12.16 SF2.Unique reads percentage ########## 
df <- readxl::read_excel("scLIANTI_hg_mm.xlsx")
colnames(df) <- c("barcode","demultiplexed_reads","hg_or_mm","unique_hg","unique_mm","unique_mapping")
tmp <- data.frame(barcode=df$barcode,unique_mapping=df$unique_mapping,unique_reads=df$unique_hg+df$unique_mm)
ggplot(tmp,mapping=aes(x=unique_mapping,y=log10(unique_reads)))+
  geom_point() +
  geom_density_2d() +
  scale_size_manual(values=c(1.5,1.5,1.5)) +
  xlab("Percentage of uniquely mapped reads (%)") +
  ylab("Log10 unique reads") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_blank(),
        axis.line=element_line(colour="black"))

########## 12.16 SF3.Single cell discrimination by unique read count ########## 
library(mixtools)
df <- readxl::read_excel("scLIANTI_hg_mm.xlsx")
colnames(df) <- c("barcode","demultiplexed_reads","hg_or_mm","unique_hg","unique_mm","unique_mapping")
#Human cells kmean clustering
tmp <- subset(df,hg_or_mm=="human")[,c(1,3,4)]
km_cluster <- kmeans(x=log10(tmp$unique_hg),centers=3)
tmp$cluster <- as.character(km_cluster$cluster)
ggplot(tmp,mapping=aes(x=log10(unique_hg))) +
  geom_histogram(aes(y=..density..,color=cluster,fill=cluster), alpha=0.5, position="identity") +
  geom_density(alpha=0.2) +
  xlab("Log10 unique reads") +
  ylab("Density") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_text("cluster"),
        axis.line=element_line(colour="black"))

#Human mouse mixed cells kmean clustering
tmp <- df[,c(1,3,4,5)]
tmp$unique_read <- tmp$unique_hg+tmp$unique_mm
km_cluster <- kmeans(x=log10(tmp$unique_read),centers=3)
tmp$cluster <- as.character(km_cluster$cluster)
ggplot(tmp,mapping=aes(x=log10(unique_read))) +
  geom_histogram(aes(y=..density..,color=cluster,fill=cluster), alpha=0.5, position="identity") +
  geom_density(alpha=0.2) +
  xlab("Log10 unique reads") +
  ylab("Density") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_text("cluster"),
        axis.line=element_line(colour="black"))

########## Test code calculate weighted mC ratio of each ENCODE
#/zfssz4/BC_RD_P1/PROJECT/fangyitong/LIANTI/scLIANTI/code/weighted_mCratio.R
df <- fread("Build_ENCFF092FNE.bed.gz",header=F)
build_name <- c()
build_wa <- c()
for (build in unique(dd$V6)) {
  build_name <- c(build_name,build)
  tmp <- subset(dd,V6==build)
  x <- c(x,length(tmp$V6))
  wa <- sum(tmp$V4*tmp$V5/100)/sum(tmp$V4)
  build_wa <- c(build_wa,wa)
}
out <- data.frame(build=build_name,mean=build_wa)

########## 12.23 Variance vs mC ratio ##########
#(positions on chr1 and covered in >80% cells)
df <- read.csv("chr1.csv",header=F,sep="\t")
colnames(df) <- c("chr","start","end","ctype","mc","cover","cells")
df$mcratio <- df$mc/df$cover
df$var <- (df$mc*(1-df$mcratio)^2+(df$cover-df$mc)*df$mcratio^2)/(df$cover-1)

ggplot(df,mapping=aes(x=mcratio*100,y=var))+
  geom_point() +
  scale_size_manual(values=c(1.5,1.5,1.5)) +
  xlab("mC ratio (%)") +
  ylab("Variance across cells") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_blank(),
        axis.line=element_line(colour="black"))

#(CG positions covered in >80% cells)
df <- read.csv("CG.csv",header=F,sep="\t")
colnames(df) <- c("chr","start","end","ctype","mc","cover","cells")
df$mcratio <- df$mc/df$cover
df$var <- (df$mc*(1-df$mcratio)^2+(df$cover-df$mc)*df$mcratio^2)/(df$cover-1)

ggplot(df,mapping=aes(x=mcratio*100,y=var))+
  geom_point() +
  xlab("mC ratio (%)") +
  ylab("Variance across cells") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_blank(),
        axis.line=element_line(colour="black"))

dd <- subset(df,var>=0.22570)
write.csv(dd,file="CGvar75pct.csv",quote=F,row.names=F)
########## 12.30 Correlation of Regulatory build of 30 cells ##########
df <- read.csv("chr1build_30cells.csv",header=F,sep="\t") %>% rename(cell=V1,build=V2,mean=V3) %>% spread(build,mean)
tmp <- data.frame(t(df)) %>% drop_na()
colnames(tmp) <- tmp[1,]
tmp <- tmp[-1,] 
tmp <- tmp %>% mutate_if(is.character, as.numeric)
pearson <- as.data.frame(cor(tmp,method="pearson"))
get_upper_tri <- function(pearson){
  pearson[lower.tri(pearson)]<-NA
  return(pearson)
}
dd <- melt(as.matrix(get_upper_tri(pearson)),na.rm=T)
ggplot(data=dd,aes(Var2,Var1,fill=value))+
  geom_tile(color="white")+
  scale_fill_gradient2(high="red",mid="yellow",limit = c(0,1),space="Lab", 
                       name="Pearson Correlation")+
  theme_minimal()+ 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=90,vjust=1,size=8,hjust=1),
        axis.text.y=element_text(size=8))+
  coord_fixed()

write.csv(pearson,file="TEAMseq_merged_pearson.csv",quote=F)

########## 01.06 Correlation of chr1 cover>20 merged TEAMseq and WGBS hela########## 
library(dplyr)
df <- read.csv("cor_pearson_final.csv",header=F,sep="\t",quote="")
df$V4[59] <- as.numeric(df$V3[59])
df$V3 <- paste(df$V1,df$V2,sep="_")
ggplot(df,aes(x=V3,y=V4)) +
  geom_bar(stat="identity",fill="#999999",color="white") +
  ylab("Pearson correlation") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_blank(),
        axis.line=element_line(colour="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))


########## 01.11 SF2.Unique reads percentage ########## 
df <- readxl::read_excel("scLIANTI_Hela293T.xlsx")
ggplot(df,mapping=aes(x=mapping,y=log10(unique_reads)))+
  geom_point() +
  geom_density_2d() +
  scale_size_manual(values=c(1.5,1.5,1.5)) +
  xlab("Percentage of uniquely mapped reads (%)") +
  ylab("Log10 unique reads") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_blank(),
        axis.line=element_line(colour="black"))

########## 01.11 SF3.Single cell discrimination by unique read count ########## 
library(mixtools)
df <- readxl::read_excel("scLIANTI_Hela293T.xlsx")
#Hela&HEK293T mixed cells kmean clustering
km_cluster <- kmeans(x=log10(df$unique_reads),centers=3)
df$cluster <- as.character(km_cluster$cluster)
ggplot(df,mapping=aes(x=log10(unique_reads))) +
  geom_histogram(aes(y=..density..,color=cluster,fill=cluster), alpha=0.5, position="identity") +
  geom_density(alpha=0.2) +
  xlab("Log10 unique reads") +
  ylab("Density") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_text("cluster"),
        axis.line=element_line(colour="black"))

########## 01.12 Basic statistics ##########
df <- readxl::read_excel("scLIANTI_Hela293T.xlsx")
#Mapping rate
ggplot(data=df,aes(x=cell_type,y=mapping,color=cell_type)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylab("Mapping rate (%)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="none",
        axis.line=element_line(colour="black"))
#Duplication
ggplot(data=df,aes(x=cell_type,y=duplication,color=cell_type)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylab("Duplication (%)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="none",
        axis.line=element_line(colour="black"))
#Coverage
ggplot(data=df,aes(x=cell_type,y=coverage,color=cell_type)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylab("Coverage (%)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="none",
        axis.line=element_line(colour="black"))
#Depth
ggplot(data=df,aes(x=cell_type,y=depth,color=cell_type)) +
  geom_violin() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  ylab("Depth (X)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="none",
        axis.line=element_line(colour="black"))

tmp <- subset(df,cell_type=="HEK293T")
ggplot(data=tmp,aes(x=coverage,y=depth)) +#,color=cell_type)) +
  geom_point() +
  xlab("Coverage (%)") +
  ylab("Depth (X)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.title=element_blank(),
        legend.position=c(0.15,0.9),
        axis.line=element_line(colour="black"))

ggplot(data=tmp,aes(x=log10(unique_reads),y=duplication)) +#,color=cell_type)) +
  geom_point() +
  xlab("Log10 Unique Reads") +
  ylab("Duplication (%)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.title=element_blank(),
        legend.position=c(0.15,0.9),
        axis.line=element_line(colour="black"))

########## 01.12 Barcode edit distance ##########
library(DNABarcodes)
library(pheatmap)
df <- readxl::read_excel("barcodes.xlsx")
analyse.barcodes(as.character(as.array(df$barcode)))
bar_mat<-barcode.set.distances(as.character(as.array(df$barcode)),metric="hamming")
tmp <- subset(as.data.frame(bar_mat@x),`bar_mat@x`>0)
ggplot(tmp,aes(x=`bar_mat@x`)) + 
  geom_density(alpha=0.6, bw=0.5) +
  xlab("Edit distances") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(colour="black"))

########## 01.13 TEST: mean mCratio of each build in H4 ##########
df <- fread("H4.bed_build.bed",header=F,sep="\t") %>% 
  select(V1,V2,V3,V5,V6,V7,V12) %>%
  rename(chr=V1,start=V2,end=V3,Ctype=V5,mC=V6,cover=V7,build=V12)
tmp <- df %>% group_by(build) %>%
  summarise(mCsum=sum(mC),coversum=sum(cover)) %>%
  mutate(mCratio=(mCsum+1)/(coversum+2))

########## 01.13 Correlation of Regulatory build of 32 HEK293T cells ##########
df <- fread("293T_build.csv",header=F,sep="\t") %>% 
  select(V5,V1,V4) %>% rename(cell=V5,build=V1,mean=V4) %>% spread(build,mean)
tmp <- data.frame(t(df)) %>% drop_na()
colnames(tmp) <- tmp[1,]
tmp <- tmp[-1,] 
x <- tmp %>% 
  #slice(grep("CTCF", row.names(.))) %>% 
  #slice(grep("enhancer:", row.names(.))) %>% 
  #slice(grep("promoter:", row.names(.))) %>% 
  #slice(grep("promoter_flanking_region:", row.names(.))) %>%
  mutate_if(is.character, as.numeric)

pearson <- as.data.frame(cor(x,method="pearson"))
pearsonCTCF <- as.data.frame(cor(x,method="pearson"))
pearsonenhancer <- as.data.frame(cor(x,method="pearson"))
pearsonpromoter <- as.data.frame(cor(x,method="pearson"))
pearsonPFR <- as.data.frame(cor(x,method="pearson"))

write.table(pearsonPFR,"Rout/pearson32cellsCG_PromoterFlankingRegion.csv",quote=F,sep="\t")
pheatmap(pearson)
pheatmap(pearsonCTCF)
pheatmap(pearsonenhancer)
pheatmap(pearsonpromoter)
pheatmap(pearsonPFR)
get_upper_tri <- function(pearson){
  pearson[lower.tri(pearson)]<-NA
  return(pearson)
}
dd <- melt(as.matrix(get_upper_tri(pearson)),na.rm=T)
ggplot(data=dd,aes(Var2,Var1,fill=value))+
  geom_tile(color="white")+
  scale_fill_gradient2(high="red",mid="yellow",limit = c(0,1),space="Lab", 
                       name="Pearson Correlation")+
  theme_minimal()+ 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=90,vjust=1,size=8,hjust=1),
        axis.text.y=element_text(size=8))+
  coord_fixed()

########## 01.17 NMF and clustering ########## 
##### TEST: generate NMF table ##### 
#A8 <- fread("A8",header=F,sep="\t") %>% select(V1,V4) %>% rename(Build=V1,A10=V4)
#A9 <- fread("A9",header=F,sep="\t") %>% select(V1,V4) %>% rename(Build=V1,A11=V4)
#alist <- list(A8,A9)
#full <- Reduce(
#  function(x, y, ...) full_join(x, y, by = "Build"),
#  alist
#)
df293T <- fread("NMFtable293T.csv",header=T,sep="\t") #%>% select(-D9)
dfhela <- fread("NMFtableHela.csv",header=T,sep="\t") #%>% select(-H4)
##### Plot NMF coverage of 293T #####
df <- nmf %>% mutate(cover = rowSums(!is.na(select(.,-Build)))) %>% 
  count(cover) %>% arrange(desc(cover)) %>% mutate(csum=cumsum(n))
ggplot(data=df,aes(x=cover,y=csum/1000000)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=(seq(1,32,by=3)),"Number of cells") +
  scale_y_continuous(breaks=(seq(0,0.5,by=0.1)),"Number of Ensembl Regulatory Builds (in Million)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_text("cluster"),
        axis.line=element_line(colour="black"))
rm(df)
##### NMF #####
library(NMF)
df <- merge(df293T,dfhela,by="Build")
nmfin <- df[complete.cases(df),] %>% select(-Build)
rm(df293T,dfhela,df)
#Look for the best k value k=2 
#refer to https://www.pnas.org/content/101/12/4164
bestk <- nmf(nmfin,2:6,nrun=10,seed=123)
plot(2:6,bestk$measures$cophenetic,type="b",col="purple")
rm(bestk)
#NMF decomposition
nmfout <- nmf(nmfin,2,nrun=10,seed=1234)
basismap(nmfout)
coefmap(nmfout)

##### kmeans #####
library(factoextra)
set.seed(1234)
kmeansout <- kmeans(t(nmfin),centers=2)
fviz_cluster(kmeansout,data=t(nmfin),labelsize=8,ellipse=F,show.clust.cent=F)

kmeansout <- kmeans(t(nmfin[1:(211*0.6),]),centers=2)
fviz_cluster(kmeansout,data=t(nmfin[1:(211*0.6),]),labelsize=8,ellipse=F,show.clust.cent=F)

kmeansout <- kmeans(t(nmfin[1:(211*0.3),]),centers=2)
fviz_cluster(kmeansout,data=t(nmfin[1:(211*0.3),]),labelsize=8,ellipse=F,show.clust.cent=F)

kmeansout <- kmeans(t(nmfin[1:(211*0.1),]),centers=2)
fviz_cluster(kmeansout,data=t(nmfin[1:(211*0.1),]),labelsize=8,ellipse=F,show.clust.cent=F)

########## 01.09 Variance vs mC ratio in HEK293T ##########
#(C positions on chr7 and covered in >80%, ie. >=26 cells)
df <- fread("chr7_26cells.csv",header=F,sep="\t") %>%
  rename(chr=V1, start=V2, end=V3, strand=V4, ctype=V5, mc=V6,cover=V7,cells=V8) %>%
  mutate(mcratio=mc/cover)
df$var <- (df$mc*(1-df$mcratio)^2+(df$cover-df$mc)*df$mcratio^2)/(df$cover-1)

ggplot(df,mapping=aes(x=mcratio*100,y=var))+
  geom_point() +
  scale_size_manual(values=c(1.5,1.5,1.5)) +
  xlab("mC ratio (%)") +
  ylab("Variance across cells") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_blank(),
        axis.line=element_line(colour="black"))
dd <- subset(df,var>0.030592)
write.csv(dd,file="Cvar_large.csv",quote=F,row.names=F)
dd <- subset(df,var<=0.030592)
write.csv(dd,file="Cvar_small.csv",quote=F,row.names=F)

#(CG positions covered in >80% cells)
df <- fread("chr7_26cells.csv",header=F,sep="\t") %>%
  rename(chr=V1, start=V2, end=V3, strand=V4, ctype=V5, mc=V6,cover=V7,cells=V8) %>%
  filter(ctype=="CG") %>%
  mutate(mcratio=mc/cover)
df$var <- (df$mc*(1-df$mcratio)^2+(df$cover-df$mc)*df$mcratio^2)/(df$cover-1)

ggplot(df,mapping=aes(x=mcratio*100,y=var))+
  geom_point() +
  xlab("mC ratio (%)") +
  ylab("Variance across cells") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8,0.8),
        legend.title=element_blank(),
        axis.line=element_line(colour="black"))

dd <- subset(df,var>=0.22570)
write.csv(dd,file="CGvar75pct.csv",quote=F,row.names=F)