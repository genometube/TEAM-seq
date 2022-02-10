library(ggplot2)
library(dplyr)
library(data.table)
library(rio)
setwd("D:/项目/nanopore+DNA+methylation/结果/cpG/")

########## Prepare files ########## 
A1 <- fread("allc_LIANTI_A1_new_type_anno1.csv.gz",header=T,sep="\t")
B1 <- fread("allc_LIANTI_B1_new_type.gz",header=F,sep="\t")
B1 <- fread("all_HUMyepUWGBS_count.csv.gz",header = F,sep = "\t")
C1 <-fread("SRR1769036_1_bismark_count.bed.gz",header=F,sep="\t")


names(B1) <-c("chr","pos","Ctype","mc","nmc")
              #"5-UTR","CDS","CGI","intron","mRNA","ncRNA","PseudoGene"
              #,"tandem_repeats","transposons")#"methylated"
#names(A1) <-c("chr","pos","Ctype","mc","allc","3-UTR",
              #"5-UTR","CDS","CGI","intron","mRNA","ncRNA","PseudoGene"
              #,"tandem_repeats","transposons")
names(C1) <-c("chr","pos","Ctype","mc","nmc")#,"3-UTR",
              #"5-UTR","CDS","CGI","intron","mRNA","ncRNA","PseudoGene"
              #,"tandem_repeats","transposons")
#####注释Ctype ---> CG/CHG/CHH/CNN########################
B1$Ctype[grepl(".CG.",B1$Ctype)] <-'CG'
B1$Ctype[grepl( ".C[ACT]G",B1$Ctype)] <-'CHG'
B1$Ctype[grepl(".C[ACT][ACT]",B1$Ctype)] <-'CHH'
B1$Ctype[grepl(".CN.|.C.N",B1$Ctype)] <-"CNN"
fwrite(B1,file = "allc_LIANTI_B1_new_type_anno1.csv.gz",sep="\t",compress="gzip")

A1$Ctype[grepl(".CG.",A1$Ctype)] <-"CG"
A1$Ctype[grepl(".C[ACT]G",A1$Ctype)] <-"CHG"
A1$Ctype[grepl(".C[ACT][ACT]",A1$Ctype)] <-"CHH"
A1$Ctype[grepl(".CN.|.C.N",A1$Ctype)] <-"CNN"
fwrite(A1,file = "allc_LIANTI_A1_new_type_anno2.csv.gz",sep="\t",compress="gzip")

########## Coverage vs Modification plot ########## 
#LIANTI测到每个CG位点的深度和他们的平均甲基化率
p_A1 <- data.table(Ctype=A1$Ctype,depth=A1$allc,ratio=A1$mc/A1$allc) %>% 
  filter(Ctype=="CG") %>% group_by(depth) %>% summarise(mean_ratio=mean(ratio*100))

p_B1 <- data.table(Ctype=B1$Ctype,depth=B1$allc,ratio=B1$mc/B1$allc) %>% 
  filter(Ctype=="CG") %>% group_by(depth) %>% summarise(mean_ratio=mean(ratio*100))

write.table(p_A1,"D:/项目/nanopore+DNA+methylation/结果/cpG/coverage_bias/TEAM-seq100pg_CGsite_covVSdep.txt",row.names = F,sep = '\t',quote=F)
write.table(p_B1,"D:/项目/nanopore+DNA+methylation/结果/cpG/coverage_bias/TEAM-seq20pg_CGsite_covVSdep.txt",row.names = F,sep = '\t',quote=F)


#其他技术测到每个CG位点的深度和他们的平均甲基化率
p_B1 <- data.table(Ctype=B1$Ctype,depth=B1$mc+B1$nmc,ratio=B1$mc/(B1$mc+B1$nmc)) %>% 
  filter(Ctype=="CG") %>% group_by(depth) %>% summarise(mean_ratio=mean(ratio*100))

write.table(p_B1,"D:/项目/nanopore+DNA+methylation/结果/cpG/coverage_bias/WGBS_CGsite_covVSdep.txt",row.names = F,sep = '\t',quote=F)

p_C1 <- data.table(Ctype=C1$Ctype,depth=C1$mc+C1$nmc,ratio=C1$mc/(C1$mc+C1$nmc)) %>% 
  filter(Ctype=="CG") %>% group_by(depth) %>% summarise(mean_ratio=mean(ratio*100))
write.table(p_C1,"D:/项目/nanopore+DNA+methylation/结果/cpG/coverage_bias/scWGBS_CGsite_covVSdep.txt",row.names = F,sep = '\t',quote=F)


########## Chromosome vs modification plot ########## 
bin_num <- 249250621/2000000 #根据chr1的长度计算应该划多少个bin
#bin_num <- 243199373/100000 #根据chr2的长度计算应该划多少个bin
#bin_num <- 198022430/100000 #根据chr3的长度计算应该划多少个bin
bin_num <- 191154276/2000000 #根据chr4的长度计算应该划多少个bin
#bin_num <- 180915260/100000 #根据chr5的长度计算应该划多少个bin
#100pg测到的CG/CH/CHH位点，划bin，计算每个bin的平均甲基化率
#df_4a <- data.table(chr=A1$chr,pos=A1$pos,Ctype=A1$Ctype,coverage=A1$mc+A1$nmc,ratio=A1$mc/(A1$mc+A1$nmc)) %>%
  #filter(chr=="chr4")  %>% filter(Ctype=="CHG") %>% filter(coverage>3)
df4_a1 <- data.table(chr=df4_a$chr,pos=df4_a$pos,Ctype=df4_a$Ctype,coverage=df4_a$coverage,ratio=df4_a$mc/df4_a$coverage) %>%
  filter(coverage>=3) %>% filter(Ctype=="CHH")
  
p_A1 <- df4_a1 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:114))#mutate(chr4=c(0:floor(bin_num)))

#df_4b <- data.table(chr=B1$chr,pos=B1$pos,Ctype=B1$Ctype,coverage=B1$mc+B1$nmc,ratio=B1$mc/(B1$mc+B1$nmc)) %>% 
 # filter(chr=="chr4") %>% filter(coverage>3) %>% filter(Ctype=="CHG")
df4_b1<- data.table(chr=df4$chr,pos=df4$pos,Ctype=df4$Ctype,coverage=df4$coverage,ratio=df4$mc/df4$coverage) %>%
  filter(coverage>=3) #%>% filter(Ctype=="CHH")#chr4

#df4_b1<- data.table(chr=df1$chr,pos=df1$pos,Ctype=df1$Ctype,coverage=df1$coverage,ratio=df1$mc/df1$coverage) %>%
  filter(coverage>=3) %>% filter(Ctype=="CHH")#chr1
p_B1 <- df4_b1 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:114)) #mutate(chr4=c(0:floor(bin_num)),na.rm=FALSE)

#############其他技术
C1=fread("3_chr_coverage/scWGBS_chr1_count.csv.gz",header = T,sep = "\t")
df <- data.table(chr=C1$chr,pos=C1$pos,Ctype=C1$Ctype,coverage=C1$coverage,ratio=C1$ratio) %>% 
   filter(coverage>=3) %>% filter(Ctype=="CG")
p_C1 <- df %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:114)) #mutate(chr4=c(0:floor(bin_num)),na.rm=FALSE)


########## Chromosome vs coverage plot ########## 
bin_num <- 249250621/2000000 #根据chr1的长度计算应该划多少个bin
#bin_num <- 243199373/100000 #根据chr2的长度计算应该划多少个bin
#bin_num <- 198022430/100000 #根据chr3的长度计算应该划多少个bin
bin_num <- 191154276/2000000 #根据chr4的长度计算应该划多少个bin
#bin_num <- 180915260/100000 #根据chr5的长度计算应该划多少个bin
#BS测到的CG/CH/CHH位点，划bin，计算每个bin的平均覆盖深度
df4 <- data.table(chr=A1$chr,pos=A1$pos,Ctype=A1$Ctype,coverage=A1$allc,mc=A1$mc) %>%
  filter(chr=="chr4")# %>% filter(Ctype=="CHH")
df1 <- data.table(chr=A1$chr,pos=A1$pos,Ctype=A1$Ctype,coverage=A1$allc,mc=A1$mc) %>%
  filter(chr=="chr1")# %>% filter(Ctype=="CHH")

df4 <- data.table(chr=B1$chr,pos=B1$pos,Ctype=B1$Ctype,coverage=B1$allc,mc=B1$mc) %>%
  filter(chr=="chr4")# %>% filter(Ctype=="CHH")
df1 <- data.table(chr=B1$chr,pos=B1$pos,Ctype=B1$Ctype,coverage=B1$allc,mc=B1$mc) %>%
  filter(chr=="chr1")# %>% filter(Ctype=="CHH")

fwrite(df1,file = "D:/项目/nanopore+DNA+methylation/结果/cpG/3_chr_coverage/TEAM-seq20pg_chr1_count.txt",sep="\t",compress="gzip")
fwrite(df4,file = "D:/项目/nanopore+DNA+methylation/结果/cpG/3_chr_coverage/TEAM-seq20pg_chr4_count.txt",sep="\t",compress="gzip")


p_A1 <- df4_a %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:114))#mutate(chr1=c(0:floor(bin_num-1)))
###画CG CHG CHH
df4_a1 <- data.table(chr=df4_a$chr,pos=df4_a$pos,Ctype=df4_a$Ctype,coverage=df4_a$coverage,mc=df4_a$mc) %>%
   filter(Ctype=="CHH")
p_A1 <- df4_a1 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:114))#mutate(chr1=c(0:floor(bin_num-1)))

#20pg 测到的CG/CH/CHH位点，划bin，计算每个bin的平均覆盖深度
#df <- data.table(chr=B1$chr,Ctype=B1$Ctype,pos=B1$pos,coverage=B1$mc+B1$nmc,ratio=B1$mc/(B1$mc+B1$nmc)) %>% 
 # filter(chr=="chr4") %>% filter(Ctype=="CHH")
df4_b1 <- data.table(chr=df4$chr,Ctype=df4$Ctype,pos=df4$pos,coverage=df4$coverage,mc=df4$mc) #%>% 
  #filter(Ctype=="CHG")

df4_b1 <- data.table(chr=df1$chr,Ctype=df1$Ctype,pos=df1$pos,coverage=df1$coverage,mc=df1$mc) %>%
  filter(Ctype=="CHH")

p_B1 <- df4_b1 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:114))#mutate(chr4=c(0:floor(bin_num)))

############其他技术
C1=fread("scWGBS_chr4_count.csv.gz",header = T,sep = "\t")
df <- data.table(chr=C1$chr,Ctype=C1$Ctype,pos=C1$pos,coverage=C1$coverage,
      ratio=C1$ratio)# %>%filter(Ctype=="CHH") # filter(chr=="chr4")
p_C1 <- df %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:94))#mutate(chr4=c(0:floor(bin_num)))


############B1,B2,B3之间的coverage###############
B1 = read.table("2_chr_coverage/TEAM-seq20pg_chr1_count.txt",header = T,sep = "\t")
B2 = fread("B2_chr1_count.csv.gz",header = T,sep = "\t")
B3 = fread("B3_chr1_count.csv.gz",header = T,sep = "\t")

df_B1 <- B1 %>% filter(Ctype=="CHH")
df_B2 <- B2 %>% filter(Ctype=="CHH")
df_B3 <- B3 %>% filter(Ctype=="CHH")

p_B1 <- df_B1 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:94))
p_B2 <- df_B2 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:94))
p_B3 <- df_B3 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:94))

##methylation ratio
B1_1=data.table(chr=B1$chr,pos=B1$pos,Ctype=B1$Ctype,mc=B1$mc,coverage=B1$coverage,ratio=B1$mc/B1$coverage)
df_B1 <- B1_1 %>% filter(coverage>=3) %>%filter(Ctype=="CG")
df_B2 <- B2 %>% filter(coverage>=3) %>%filter(Ctype=="CG")
df_B3 <- B3 %>% filter(coverage>=3) %>%filter(Ctype=="CG")

p_B1 <- df_B1 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:94))
p_B2 <- df_B2 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:94))
p_B3 <- df_B3 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:94))


########## CpG coverage vs Number of CpGs plot ########## 
#LIANTI在每个深度下测到C位点的个数
p_A1 <- data.table(Ctype=A1$Ctype,depth=A1$allc) %>% 
  #filter(Ctype=="CG") %>% #只看CG/CHG/CHH位点,注释掉本行看all C位点
  group_by(depth)
p_B1 <- data.table(Ctype=B1$Ctype,depth=B1$allc) %>% 
  #filter(Ctype=="CG") %>% #只看CG/CHG/CHH位点,注释掉本行看all C位点
  group_by(depth)


#其他技术在每个深度下测到C位点的个数
p_B1 <- data.table(Ctype=B1$Ctype,depth=B1$mc+B1$nmc) %>% 
  #filter(Ctype=="CHH") %>% #只看CG/CHG/CHH位点,注释掉本行看all C位点
  group_by(depth)
p_C1 <- data.table(Ctype=C1$Ctype,depth=C1$mc+C1$nmc) %>% 
  filter(Ctype=="CHG") %>% #只看CG/CHG/CHH位点,注释掉本行看all C位点
  group_by(depth)


aa=table(p_A1$depth)
nam=rownames(aa)
num=as.numeric(aa)
final=as.data.frame(cbind(nam,num),stringAsFactor=F)
p_A1=final[c(1:30),]
p_A1$num=as.numeric(p_A1$num)

######### Barplot mC ratio of functional elements ########## 
A1 <- read.csv("5_genome_element/allc_LIANTI_A1_genome_element.txt",sep="\t",header=F)
B1 <- read.csv("5_genome_element/allc_LIANTI_B1_genome_element.txt",sep="\t",header=F)
C1 <- read.csv("5_genome_element/all_HUMyepUWGBS_genome_element.txt",sep="\t",header=F)
D1 <- read.csv("5_genome_element/scWGBS_genome_element.txt",sep="\t",header=F)

##head(A1)如下:##
#3-UTR 	 allc 	 0.04073796 	 CG 	 0.6330561 	 CHG 	 0.01221225 	 CHH 	 0.01973421 
#5-UTR 	 allc 	 0.0403175 	 CG 	 0.6043971 	 CHG 	 0.01201075 	 CHH 	 0.02017632 

A1 <- A1[,c(1,3,5,7,9)]
names(A1) <- c("type","allc","CG","CHG","CHH")
B1 <- B1[,c(1,3,5,7,9)]
names(B1) <- c("type","allc","CG","CHG","CHH")
C1 <- C1[,c(1,3,5,7,9)]
names(C1) <- c("type","allc","CG","CHG","CHH")
D1 <- D1[,c(1,3,5,7,9)]
names(D1) <- c("type","allc","CG","CHG","CHH")

df <- data.frame(seq=c(rep("TEAM-seq_100pg",10),rep("TEAM-seq_20pg",10)),
                 type=c(A1$type,B1$type),
                 ratio=c(A1$CG,B1$CG))
df <- data.frame(seq=c(rep("WGBS",10),rep("scWGBS",10)),
                 type=c(C1$type,D1$type),
                 ratio=c(C1$CG,D1$CG))
##all method
df <- data.frame(seq=c(rep("TEAM-seq_100pg",10),rep("TEAM-seq_20pg",10),rep("WGBS",10),rep("scWGBS",10)),
                 type=c(A1$type,B1$type,C1$type,D1$type),
                 ratio=c(A1$allc,B1$allc,C1$allc,D1$allc))

####### chromosone C information #############################
P_A1 <- data.table(chr=A1$chr,Ctype=A1$Ctype) #%>% filter(Ctype=="CG")####ctype来计算total C/CG/CHG/CHH
P_A11 <- data.table(chr=A11$chr,Ctype=A11$Ctype)#%>% filter(Ctype=="CG")
table(P_A1$chr)
P_A1_1<- data.table(chr=P_A1$chr,Ctype=P_A1$Ctype) %>% filter(Ctype=="CHH")
P_A11_1<- data.table(chr=P_A11$chr,Ctype=P_A11$Ctype)%>% filter(Ctype=="CHG")
table(P_A1_1$chr)
table(P_A11_1$chr)

P_B1 <- data.table(chr=B1$chr,Ctype=B1$Ctype) %>% filter(Ctype=="CHH") 
table(P_B1$chr)

######################### gene element###################
p_A1 <- data.table(element=A1$`3-UTR`,Ctype=A1$Ctype) %>%
  filter(element > 0)
table(p_A1$Ctype)


p_B1 <- data.table(element=B1$`3-UTR`,
                   Ctype=B1$Ctype) %>%
  filter(element > 0)  ##### coverage

p_C1 <- data.table(element=C1$`3-UTR`,Ctype=C1$Ctype) %>%
  filter(element > 0)

p_B1 <- data.table(element=B1$CDS,mc=B1$mc,allc=B1$mc+B1$nmc,ratio=B1$mc/(B1$mc+B1$nmc),Ctype=B1$Ctype) %>%
  filter(element > 0) %>%filter(Ctype=="CG") %>% filter(allc >3)#####methylation ratio


####################################sc_LIANTI##############################
#depth
setwd("D:/项目/nanopore+DNA+methylation/单细胞/cell64/")
a=read.table("depth_coverage.txt",header = T,sep = '\t')
a$unique.reads=as.numeric(a$unique.reads)/1000000
a$CG_cover_num=a$CG_cover_num/1000000
plot(x=a$unique.reads/10000000,y = a$CG_cover_num/1000000)
aa<- drm(CG_cover_num~unique.reads,data = a,fct=L.3())
plot(aa,col ="DarkCyan",ylab=expression("CGs covered("~"x10"^"6"~")"),
     xlab=expression("Align unique reads ("~"x10"^"7"~")"),pch =20,lwd=1)
#x=a$unique.reads,y=a$CG_cover_num

##cells_c_count 
#log10(a$V3) Log10 C counts cell_coverCG_counts.txt xlab="Number of shared Cells",ylab="chromosomal CG coverage(%)",col='SkyBlue1')
a=read.table("coveredC_cell.txt",header = T,sep = '\t')
plot(x=a$merge_cell,y=a$coverage,pch=19,type="o",lwd=2,xlab="Number of combine cells",ylab="Covered C(%)",col='SteelBlue3')

#######  NMF ###########
library(NMF)
library(tsne)
library(doMPI)
a=read.csv("NMFtable.csv",header = T,sep = "\t")

rownames(a) <- a[,1]
a=a[,-1]
a=a[complete.cases(a),]
a=as.matrix(a)
res=nmf(a,2,seed=78,nrun=10)
basismap(res)
coefmap(res)
consensusmap(res)
#k-means
library(ggfortify)
library(ggplot2)
aa=t(a)
re=kmeans(aa,2)
autoplot(kmeans(aa, 2,nstart = 10), data = aa)
#test
data("esGolub")
esGolub <- esGolub[1:200,]
esGolub$Sample <- NULL
estim.r <- nmf(esGolub, 4, nrun=10, seed=123456)
plot(estim.r)
plot(2:6,estim.r$measures$cophenetic, type="b", col="purple")
#test finish
