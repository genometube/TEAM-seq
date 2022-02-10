library(ggplot2)
library(dplyr)
library(data.table)
library(rio)
setwd("D:/项目/nanopore+DNA+methylation/结果/cpG/")
a=read.table("3_CGI_coverage/A1B1_scWGBS_WGBS_windepth.txt",header = F,sep = "\t")##/3.9
#CGI aver coverage A1_scWGBS_total.txt
a=read.table("8_20pg_total_analysis/all_20pg_aver_win_depth.txt",header = F,sep = "\t")
ggplot(data=a,mapping = aes(x=V1,y=V9,colour=V5))+geom_line(size=1,alpha=0.8)+
  labs(x="CGI",y="Average coverage depth")+scale_color_manual(values = c("#868686FF","blue","#0073C2FF","orange"))+
  geom_vline(xintercept = c(50,70),linetype="dashed")+
  theme_bw() + theme(legend.title=element_blank(),axis.text.x =element_blank(),axis.ticks.x=element_blank(),legend.position = "top",panel.border=element_blank(),
                     panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                     axis.line = element_line(colour = "black"))
a=read.table("2_CGI_coverage/all_lianti_scwgbs.txt",header = F,sep = "\t") ########"Ara_BS_EM.txt" name="sample",
a=read.table("8_20pg_total_analysis/all20pg_avermethy.txt",header = F,sep = "\t") ########"Ara_BS_EM.txt" name="sample",

#CGI aver methylation
ggplot(data=a,mapping = aes(x=V1,y=V4*100,colour=V5))+geom_line(size=0.8,alpha=0.8)+  ##"dodgerblue", "goldenrod1","ForestGreen"  ##"#868686FF","#0073C2FF","Chartreuse","orange"
  labs(x="CGI",y="Average modification(%)")+scale_color_manual(values = c("dodgerblue", "goldenrod1","ForestGreen"))+
  theme_bw() + theme(axis.text.x =element_blank(),axis.ticks.x=element_blank(),legend.position = "top",panel.border=element_blank(),
                     panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                     axis.line = element_line(colour = "black"))
#sc_LANITI CGI methy rate
a=read.table("../../单细胞/cell64/scWGBS_HEK293T.enhancer.final.txt",header=F,sep='\t')
ggplot(data=a,mapping = aes(x=V5,y=V4*100,group=V6))+geom_line(size=1,alpha=0.8,colour="black")+  ##"dodgerblue", "goldenrod1","ForestGreen"  ##"#868686FF","#0073C2FF","Chartreuse","orange"
  labs(x="Enhancer",y="Average modification(%)")+#scale_color_manual(values = c("dodgerblue",'DeepSkyBlue3', "goldenrod1",alpha=0.8,"ForestGreen"，DeepSkyBlue2))+
  geom_vline(xintercept = c(50,70),linetype="dashed")+ylim(0,100)+
  theme_bw() + theme(axis.text.x =element_blank(),axis.ticks.x=element_blank(),legend.position = "top",panel.border=element_blank(),
                     panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                     axis.line = element_line(colour = "black"))
#deviation
a=read.table("../../单细胞/all_CGI_CGfinal.bed",header=F,sep='\t')
aa=data.table(window=a$V5,methyratio=a$V4) %>% group_by(window) %>% summarise(sd_win=sd(methyratio))
export(aa,"clipboard")
#plot(aa$window,aa$sd_win)

#sc_LANITI CGI aver coverage
a=read.table("../../单细胞/cell64/239T_CGI_averdepth_final.txt",header=F,sep='\t')
a$V5=as.numeric(a$V5)
ggplot(data=a,mapping = aes(x=V5,y=V4,group=V6))+geom_line(size=1,alpha=0.6,colour="#0073C2FF")+
  labs(x="CGI",y="Average coverage")+scale_color_manual(values = c("blue","blue","#0073C2FF","orange"))+
  geom_vline(xintercept = c(50,70),linetype="dashed")+ylim(0,2.5)+
  theme_bw() + theme(axis.text.x =element_blank(),axis.ticks.x=element_blank(),legend.position = "top",panel.border=element_blank(),
                     panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                     axis.line = element_line(colour = "black"))


########combine
a=read.table("new_coverage.txt",header = F,sep = "\t")
a$V3 <-factor(a$V3,ordered=TRUE,levels=c("B1","B1+B2","B1+B2+B3","A1"))
ggplot(a,aes(x=V1,y=V2,fill=V3)) + xlab("X_depth")+ylab("coverage(%)")+
  scale_fill_manual(values = c("#E889BD","#FC8A61","#67C2A3","#8EA0C9"))+
  geom_bar(stat="identity",position=position_dodge(),alpha=0.9)+
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.direction = "horizontal",
        legend.position = c(0.28,1),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"))

a=read.table("B1_conbine.txt",header = T,sep = "\t")
a$V3 <-factor(a$V3,ordered=TRUE,levels=c("B1","B1+B2","B1+B2+B3","A1"))
ggplot(a,aes(x=average.depth,y=coverage)) + xlab("X_depth")+ylab("coverage(%)")+
  scale_color_manual(values = c("#E889BD","#FC8A61","#67C2A3","#8EA0C9"))+
  geom_line(linetype="twodash",color="#67C2A3",cex=1)+geom_point(color="#67C2A3",alpha=0.9,)+
  #stat_summary(geom="line",) + 
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.direction = "horizontal",
        legend.position = c(0.28,1),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"))
plot(x=a$average.depth,y=a$coverage,pch=19,type="o",lwd=2,xlab="X_depth",ylab="Coverage(%)",col='#67C2A3')

########## Prepare files ########## 
A1 <- fread("allc_LIANTI_A1_new_type_anno1.csv.gz",header=T,sep="\t")
B1 <- fread("allc_LIANTI_B1_new_type.gz",header=F,sep="\t")
B1 <- fread("all_HUMyepUWGBS_count.csv.gz",header = F,sep = "\t")
C1 <-fread("SRR1769036_1_bismark_count.bed.gz",header=F,sep="\t")
#A11 <- fread("allc_LIANTI_A1_new_type_anno2.csv.gz",header=T,sep="\t")


names(A1) <- c("Ctype","mc","allc")
names(B1) <-c("chr","pos","Ctype","mc","allc")# B1
B1=subset(B1,B1$chr != "Control")



names(B1) <-c("chr","pos","Ctype","mc","nmc")#,"3-UTR",#######WGBS
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

#画图
p_A1=read.table("A1_CGsite_covVSdep.txt",header = T,sep = "\t")
ggplot(p_A1,mapping=aes(x=depth,y=mean_ratio),) + #画BS
  geom_point(aes(colour="TEAM-seq_100pg"),alpha=0.7,shape=17) +#geom_smooth(method=lm)+ #画BS scale_y_continuous(breaks = seq(0,100,20))+
  geom_point(p_B1,mapping=aes(x=depth,y=mean_ratio,colour="TEAM-seq_20pg"),alpha=0.6,shape=19)+ #geom_smooth(method=lm)+#画EM
  #geom_point(p_C1,mapping=aes(x=depth,y=mean_ratio,colour="scWGBS_100pg"),alpha=0.6,shape=18)+ #c("darkgrey", "orange","#0073C2FF")  c("darkred","#0073C2FF")
  xlim(0,100) + ylim(0,100)+
  xlab("Coverage") +
  scale_color_manual(values = c("IndianRed3","SteelBlue3"))+
  ylab("Methylation ratio (%)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"))

########## Chromosome vs modification plot ########## 
bin_num <- 249250621/2000000 #根据chr1的长度计算应该划多少个bin
#bin_num <- 243199373/100000 #根据chr2的长度计算应该划多少个bin
#bin_num <- 198022430/100000 #根据chr3的长度计算应该划多少个bin
bin_num <- 191154276/2000000 #根据chr4的长度计算应该划多少个bin
#bin_num <- 180915260/100000 #根据chr5的长度计算应该划多少个bin
#100pg测到的CG/CH/CHH位点，划bin，计算每个bin的平均甲基化率
#df <- data.table(chr=A1$chr,pos=A1$pos,Ctype=A1$Ctype,coverage=A1$mc+A1$nmc,ratio=A1$mc/(A1$mc+A1$nmc)) %>%
  #filter(chr=="chr4")  %>% filter(Ctype=="CHG") %>% filter(coverage>3)
df4_a1 <- data.table(chr=df4_a$chr,pos=df4_a$pos,Ctype=df4_a$Ctype,coverage=df4_a$coverage,ratio=df4_a$mc/df4_a$coverage) %>%
  filter(coverage>=3) %>% filter(Ctype=="CHH")
  
p_A1 <- df4_a1 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:114))#mutate(chr4=c(0:floor(bin_num)))
#EM测到的CG/CH/CHH位点，划bin，计算每个bin的平均甲基化率
#df <- data.table(chr=B1$chr,pos=B1$pos,Ctype=B1$Ctype,coverage=B1$mc+B1$nmc,ratio=B1$mc/(B1$mc+B1$nmc)) %>% 
 # filter(chr=="chr4") %>% filter(coverage>3) %>% filter(Ctype=="CHG")
df4_b1<- data.table(chr=df4$chr,pos=df4$pos,Ctype=df4$Ctype,coverage=df4$coverage,ratio=df4$mc/df4$coverage) %>%
  filter(coverage>=3) #%>% filter(Ctype=="CHH")#chr4

df4_b1<- data.table(chr=df1$chr,pos=df1$pos,Ctype=df1$Ctype,coverage=df1$coverage,ratio=df1$mc/df1$coverage) %>%
  filter(coverage>=3) %>% filter(Ctype=="CHH")#chr1
p_B1 <- df4_b1 %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:114)) #mutate(chr4=c(0:floor(bin_num)),na.rm=FALSE)

#############其他技术
C1=fread("3_chr_coverage/scWGBS_chr1_count.csv.gz",header = T,sep = "\t")
df <- data.table(chr=C1$chr,pos=C1$pos,Ctype=C1$Ctype,coverage=C1$coverage,ratio=C1$ratio) %>% 
   filter(coverage>=3) %>% filter(Ctype=="CG")
p_C1 <- df %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:114)) #mutate(chr4=c(0:floor(bin_num)),na.rm=FALSE)


#画图
ggplot(p_C1,mapping=aes(x=chr4,y=mean_ratio,)) + #画BS
  geom_line(aes(colour="scWGBS"),alpha=0.7) + #画BS "#0073C2FF"
  #geom_line(p_B1,mapping=aes(x=chr4,y=mean_ratio,colour="TEAM-seq_20pg"),alpha=0.7) + #画EM
  xlab("Chr1") +scale_color_manual(values = c( "orange","SteelBlue3"))+
  ylab("Methylation ratio (%)") +labs(title = "CG")+ ylim(0,100)+
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
##Repeat for CG,CHG,CHH
##Repeat for Chr1,2,3,4,5


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


p_B1 <- p_B1[-95,]
p_A1 <- p_A1[-95,]

############其他技术
C1=fread("scWGBS_chr4_count.csv.gz",header = T,sep = "\t")
df <- data.table(chr=C1$chr,Ctype=C1$Ctype,pos=C1$pos,coverage=C1$coverage,
      ratio=C1$ratio)# %>%filter(Ctype=="CHH") # filter(chr=="chr4")
p_C1 <- df %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*2000000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:94))#mutate(chr4=c(0:floor(bin_num)))

p_C1 <- p_C1[-95,]

#画图
ggplot(p_C1,mapping=aes(x=chr4,y=mean_coverage)) + #画BS"darkred", "#0073C2FF"
  geom_line(aes(colour="scWGBS")) + #画WGBS TEAM-seq_100pg
  #geom_line(p_B1,mapping=aes(x=chr4,y=mean_coverage,colour="TEAM-seq_20pg")) + #画scWGBS  "orange","#0073C2FF"
  xlab("Chr1") +labs(title = "All C sites coverage")+
  ylab("Coverage") +scale_color_manual(values = c("orange"))+ #"IndianRed3","SteelBlue3"
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
##Repeat for CG,CHG,CHH
##Repeat for Chr1,2,3,4,5

############B1,B2,B3之间的coverage###############3
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
p_B1 <- p_B1[-95,]
p_B2 <- p_B2[-95,]
p_B3 <- p_B3[-95,]



ggplot(p_B1,mapping=aes(x=chr4,y=mean_coverage)) + #画BS"darkred", "#0073C2FF"
  geom_line(aes(colour="TEAM-seq_20pg_1")) + #画WGBS TEAM-seq_100pg
  geom_line(p_B2,mapping=aes(x=chr4,y=mean_coverage,colour="TEAM-seq_20pg_2")) + #画scWGBS  "orange","#0073C2FF"
  geom_line(p_B3,mapping=aes(x=chr4,y=mean_coverage,colour="TEAM-seq_20pg_3")) +
  xlab("Chr4") +labs(title = "CHH")+ ylim(0,10)+
  ylab("Coverage") +scale_color_manual(values = c("#0073C2FF","ForestGreen","DarkGoldenrod3"))+ #"IndianRed3","SteelBlue3"
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
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

ggplot(p_B1,mapping=aes(x=chr4,y=mean_ratio)) + #画BS
  geom_line(aes(colour="TEAM-seq_20pg_1"),alpha=0.8) + #画BS "#0073C2FF"
  geom_line(p_B2,mapping=aes(x=chr4,y=mean_ratio,colour="TEAM-seq_20pg_2"),alpha=0.8) + #画EM
  geom_line(p_B3,mapping=aes(x=chr4,y=mean_ratio,colour="TEAM-seq_20pg_3"),alpha=0.8) + #画EM
  xlab("Chr4") +scale_color_manual(values = c( "#0073C2FF","ForestGreen","DarkGoldenrod3"))+
  ylab("Methylation ratio (%)") +labs(title = "CpG")+ ylim(0,100)+
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
####################截取染色体部分位置###########

df4_a=read.table("3_chr_coverage/TEAM-seq100pg_chr4_count.txt",header = T,sep = "\t")
C1=fread("WGBS_chr4_count.csv.gz",header = T,sep = "\t")
D1=fread("scWGBS_chr4_count.csv.gz",header = T,sep = "\t")

p_A1 <- data.table(chr=df4_a$chr,pos=df4_a$pos,Ctype=df4_a$Ctype,coverage=df4_a$coverage,mc=df4_a$mc) %>%
  filter(6500000< pos & pos< 7000000)# %>% filter(Ctype=="CHH")

p_C1 <- data.table(chr=C1$chr,Ctype=C1$Ctype,pos=C1$pos,coverage=C1$coverage,
                 ratio=C1$ratio) %>%filter(6500000< pos & pos< 7000000) # filter(chr=="chr4")
p_D1 <- data.table(chr=D1$chr,Ctype=D1$Ctype,pos=D1$pos,coverage=D1$coverage,
                   ratio=D1$ratio) %>%filter(6500000< pos & pos< 7000000) # filter(chr=="chr4")
ggplot(p_C1,aes(x=pos,y=coverage)) + #画BS"darkred", "#0073C2FF"
  geom_point(aes(colour="WGBS"),alpha=0.8,size=0.6) + #画WGBS TEAM-seq_100pg
  #geom_line(p_B1,mapping=aes(x=chr4,y=mean_coverage,colour="TEAM-seq_20pg")) + #画scWGBS  "orange","#0073C2FF"
  xlab("Chr4") +labs(title = "All C sites coverage")+  #ylim(0,16)+
  ylab("Coverage") +scale_color_manual(values = c("SteelBlue3"))+ #"IndianRed3","SteelBlue3" 
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black")
        )



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

aa=table(p_B1$depth)
nam=rownames(aa)
num=as.numeric(aa)
final=as.data.frame(cbind(nam,num),stringAsFactor=F)
p_B1=final[c(1:30),]
p_B1$num=as.numeric(p_B1$num)

ggplot(p_A1,aes(x=order(as.numeric(nam)),y=num/10000000,fill="TEAM-seq_100pg")) + #酌情除以10的5次方或者6次方
  geom_bar(stat="identity",position="stack",alpha=0.5) + #画BS
  geom_bar(p_B1,mapping=aes(x=order(as.numeric(nam)),y=num/10000000,fill="TEAM-seq_20pg"),alpha=0.5,stat="identity",position="stack") + #酌情除以10的5次方或者6次方
  xlim(0,20) +scale_fill_manual(values = c("IndianRed3","SteelBlue3"))+#c("orange", "#0073C2FF"))+
  #scale_x_continuous(expand=c(0,0),limits=c(0,30)) +
  #scale_y_continuous(expand=c(0,0),limits=c(0,NA)) + #limits值可以根据实际情况调整
  annotate("text",x=0,y=32,label="1e7",size=3.5) + #1e5和1e6的位置，y值可以根据实际情况调整
  #xlab("Coverage per CG position") +
  xlab("Coverage per C") +
  #ylab("Number of CG positions") + 
  ylab("Number of Cs") + 
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8, 0.9),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"))
##Repeat for allC,CG,CHG,CHH

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


##画图
ggplot(df,aes(x=type,y=ratio*100,fill=seq)) +
  geom_bar(stat="identity",color="black",position=position_dodge(),alpha=0.6) +
  #ylim(0,50)+
  labs(title = "All C sites")+
  ylab("Mean methylation ratio (%)") +scale_fill_manual(values = c("#868686FF","#0073C2FF","ForestGreen","DarkGoldenrod3"))+ #DarkGoldenrod1
  #scale_y_continuous(expand=c(0,0),limits=c(0,NA)) +
  theme(panel.background=element_blank(),
        #legend.position=c(0.2, 0.9),
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 0.65, vjust = 0.8),
        axis.line=element_line(colour="black"),
        axis.title.x=element_blank()) 

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

###################cover overlap###########################
B1 <- B1 %>% filter(Ctype=="CG")
B2 <- B2 %>% filter(Ctype=="CG")
B3 <- B3 %>% filter(Ctype=="CG")

B1$anno= paste(B1$chr,sep = "_",B1$pos)
B2$anno= paste(B2$chr,sep = "_",B2$pos)
B3$anno= paste(B3$chr,sep = "_",B3$pos)
test=venn.diagram(list(df_B1=B1$anno,df_B2=B2$anno,df_B3=B3$anno),
                  col = "transparent",cat.dist = c(-0.4, -0.4,-0.04),
                  cat.pos = c(170,-160,170) ,cat.default.pos = "outer",category = c( "B1","B2","B3"),cat.cex = 1,filename = NULL,
                  #print.mode="percent",
                  fill = c("dodgerblue", "goldenrod1","ForestGreen"))
grid.newpage()
grid.draw(test)
draw.triple.venn(area1=1282918,area2=1340749,area3=1263547,
                 n12 = 509897,n23=472494,n13=483049,n123=216399,
                 #n12=9291513,n23=8575875,n13=8936118,n123=3845985,  b2=25470852
                 category=c("B1","B2","B3"),
                 print.mode="percent",
                 fill = c("dodgerblue", "goldenrod1","ForestGreen"))
library(eulerr)
s=list(df_B1=B1$anno,df_B2=B2$anno,df_B3=B3$anno)
plot(euler(s,shape = "ellipse"), quantities = TRUE)
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
####sc_chr_coverage############
a=read.table("chr7_cover_CpG_final.txt",header = T,sep = '\t')
ggplot(a,mapping=aes(x=chr4,y=mean_coverage,group=cells)) + #画BS"darkred", "#0073C2FF"
  geom_line(aes(colour="scTEAM-seq"),alpha=0.5,size=0.3) + #画WGBS TEAM-seq_100pg
  #geom_line(p_B2,mapping=aes(x=chr4,y=mean_coverage,colour="TEAM-seq_20pg_2")) + #画scWGBS  "orange","#0073C2FF"
  #geom_line(p_B3,mapping=aes(x=chr4,y=mean_coverage,colour="TEAM-seq_20pg_3")) +
  xlab("Chr7") +labs(title = "CpG")+#ylim(0,5)+#xlim(0,93)+
  ylab("Coverage") +scale_color_manual(values = c("SteelBlue3","DarkCyan","SteelBlue3"))+ #"IndianRed3","SteelBlue3"
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#### New  bam_chr_coverage################
p_C1=read.table("2_chr_coverage/bam_chr_coverage/scWGBS_chr4_cover_allc.txt",header = T,sep = "\t")
p_B1=read.table("2_chr_coverage/bam_chr_coverage/B1_cover_CpG.txt",header = T,sep = "\t")
p_C1=p_C1[-95,]
p_B1=p_B1[-95,]
ggplot(p_C1,mapping=aes(x=chr1,y=mean_coverage)) + #画BS"darkred", "#0073C2FF"
  geom_line(aes(colour="TEAM-seq_100pg")) + #画WGBS TEAM-seq_100pg
  #geom_line(p_B1,mapping=aes(x=chr1,y=mean_coverage,colour="TEAM-seq_20pg")) + #画scWGBS  "orange","#0073C2FF"
  xlab("Chr1") +labs(title = "All CpG sites coverage")+#ylim(0,15)+
  ylab("Coverage") +scale_color_manual(values = c("IndianRed3","SteelBlue3"))+ #""orange"","SteelBlue3"
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
####### C tpye sites coverage violin#########
a=read.table("CG_ocver.txt",header = F,sep = "\t")
ggplot(a,aes(x=V2,y=V1,col =V2 ,group=V2))+geom_violin() +
  geom_jitter(position=position_jitter(0.2))+
  xlab("CpG CHG CHH")+ylab("Coverage(%)")+ylim(0,25)+scale_color_manual(values = c( "#0073C2FF","ForestGreen","DarkGoldenrod3"))+
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="none",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#"dodgerblue", "goldenrod1","ForestGreen"

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