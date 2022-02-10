########## Load libraries ########## 
library(ggplot2)
library(dplyr)
library(data.table)
setwd("C:/Users/fangyitong/Desktop/2021Research/LIANTI/tables")

########## Prepare files ########## 
BS <- fread("Ara_BS.csv.gz",header=TRUE,sep="\t")
EM <- fread("Ara_EM.csv.gz",header=TRUE,sep="\t")
#把功能元件改成0（如果该位点不是某功能元件）和1
elements <- names(BS)[6:11]
BS[,(elements) := lapply(.SD, function(x) ifelse(x>0,1,0)), .SDcols=(elements)]
EM[,(elements) := lapply(.SD, function(x) ifelse(x>0,1,0)), .SDcols=(elements)]
rm(elements)

########## Coverage vs Modification plot ########## 
#BS测到每个CG位点的深度和他们的平均甲基化率
p_BS <- data.table(Ctype=BS$Ctype,depth=BS$mc+BS$nmc,ratio=BS$mc/(BS$mc+BS$nmc)) %>% 
  filter(Ctype=="CG") %>% group_by(depth) %>% summarise(mean_ratio=mean(ratio*100))
#EM测到每个CG位点的深度和他们的平均甲基化率
p_EM <- data.table(Ctype=EM$Ctype,depth=EM$mc+EM$nmc,ratio=EM$mc/(EM$mc+EM$nmc)) %>% 
  filter(Ctype=="CG") %>% group_by(depth) %>% summarise(mean_ratio=mean(ratio*100))
#画图
ggplot(p_BS,mapping=aes(x=depth,y=mean_ratio)) + #画BS
  geom_point(aes(colour="BS"),size=0.9) + #画BS 
  geom_point(p_EM,mapping=aes(x=depth,y=mean_ratio,colour="EM"),size=0.9) + #画EM
  xlim(0,80) + 
  xlab("Coverage (CG)") +
  ylab("Methylation ratio (%)") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"))
#下载图片的时候长=500宽=400

########## Chromosome vs modification plot ########## 
#bin_num <- 30427671/100000 #根据chr1的长度计算应该划多少个bin
#bin_num <- 19698289/100000 #根据chr2的长度计算应该划多少个bin
#bin_num <- 23459830/100000 #根据chr3的长度计算应该划多少个bin
bin_num <- 18585056/100000 #根据chr4的长度计算应该划多少个bin
#bin_num <- 26975502/100000 #根据chr5的长度计算应该划多少个bin
#BS测到的CG/CH/CHH位点，划bin，计算每个bin的平均甲基化率
df <- data.table(chr=BS$chr,pos=BS$pos,Ctype=BS$Ctype,coverage=BS$mc+BS$nmc,ratio=BS$mc/(BS$mc+BS$nmc)) %>%
  filter(chr=="Chr4") %>% filter(coverage>3) %>% filter(Ctype=="CHG")
p_BS <- df %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*100000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:floor(bin_num)))
#EM测到的CG/CH/CHH位点，划bin，计算每个bin的平均甲基化率
df <- data.table(chr=EM$chr,pos=EM$pos,Ctype=EM$Ctype,coverage=EM$mc+EM$nmc,ratio=EM$mc/(EM$mc+EM$nmc)) %>% 
  filter(chr=="Chr4") %>% filter(coverage>3) %>% filter(Ctype=="CHG")
p_EM <- df %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*100000)) %>%
  group_by(bin) %>% summarise(mean_ratio=mean(ratio*100)) %>% mutate(chr4=c(0:floor(bin_num)))
#画图
ggplot(p_BS,mapping=aes(x=chr4,y=mean_ratio)) + #画BS
  geom_line(aes(colour="BS")) + #画BS 
  geom_line(p_EM,mapping=aes(x=chr4,y=mean_ratio,colour="EM")) + #画EM
  xlab("Chr4") +
  ylab("Methylation ratio (%)") +
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
#下载图片的时候长=500宽=400

########## Chromosome vs coverage plot ########## 
#bin_num <- 30427671/100000 #根据chr1的长度计算应该划多少个bin
#bin_num <- 19698289/100000 #根据chr2的长度计算应该划多少个bin
#bin_num <- 23459830/100000 #根据chr3的长度计算应该划多少个bin
bin_num <- 18585056/100000 #根据chr4的长度计算应该划多少个bin
#bin_num <- 26975502/100000 #根据chr5的长度计算应该划多少个bin
#BS测到的CG/CH/CHH位点，划bin，计算每个bin的平均覆盖深度
df <- data.table(chr=BS$chr,pos=BS$pos,Ctype=BS$Ctype,coverage=BS$mc+BS$nmc) %>%
  filter(chr=="Chr4") %>% filter(Ctype=="CHG")
p_BS <- df %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*100000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:floor(bin_num)))
#EM测到的CG/CH/CHH位点，划bin，计算每个bin的平均覆盖深度
df <- data.table(chr=EM$chr,Ctype=EM$Ctype,pos=EM$pos,coverage=EM$mc+EM$nmc,ratio=EM$mc/(EM$mc+EM$nmc)) %>% 
  filter(chr=="Chr4") %>% filter(Ctype=="CHG")
p_EM <- df %>% mutate(bin=cut(pos,breaks=c(0:bin_num)*100000)) %>%
  group_by(bin) %>% summarise(mean_coverage=mean(coverage)) %>% mutate(chr4=c(0:floor(bin_num)))
#画图
ggplot(p_BS,mapping=aes(x=chr4,y=mean_coverage)) + #画BS
  geom_line(aes(colour="BS")) + #画BS 
  geom_line(p_EM,mapping=aes(x=chr4,y=mean_coverage,colour="EM")) + #画EM
  xlab("Chr4") +
  ylab("Coverage") +
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
#下载图片的时候长=500宽=400


########## CpG coverage vs Number of CpGs plot ########## 
#BS在每个深度下测到C位点的个数
p_BS <- data.table(Ctype=BS$Ctype,depth=BS$mc+BS$nmc) %>% 
  filter(Ctype=="CHH") %>% #只看CG/CHG/CHH位点,注释掉本行看all C位点
  group_by(depth)
#EM在每个深度下测到C位点的个数
p_EM <- data.table(Ctype=EM$Ctype,depth=EM$mc+EM$nmc) %>% 
  filter(Ctype=="CHH") %>% #只看CG/CHG/CHH位点,注释掉本行看all C位点
  group_by(depth)

ggplot(p_BS,mapping=aes(x=depth,y=stat(count/1000000))) + #酌情除以10的5次方或者6次方
  geom_histogram(aes(fill="BS"),alpha=0.5) + #画BS
  geom_histogram(p_EM,mapping=aes(x=depth,y=stat(count/1000000),fill="EM"),alpha=0.5) + #酌情除以10的5次方或者6次方
  xlim(0,100) + 
  scale_x_continuous(expand=c(0,0),limits=c(0,101)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,NA)) + #limits值可以根据实际情况调整
  annotate("text",x=3,y=4,label="1e6",size=3.5) + #1e5和1e6的位置，y值可以根据实际情况调整
  xlab("Coverage per CHH position") +
  #xlab("Coverage per C") +
  ylab("Number of CHH positions") + 
  #ylab("Number of Cs") + 
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.8, 0.9),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"))
##Repeat for allC,CG,CHG,CHH



########## Barplot mC ratio of functional elements ########## 
BS <- read.csv("Ara_BS_Element_mCratio.csv",sep=",",header=T)
EM <- read.csv("Ara_EM_Element_mCratio.csv",sep=",",header=T)
##head(BS)如下:##
#BS      C     CG    CHG    CHH#
#5-UTR  6.286  4.862  3.043  7.423#
#CDS  6.773 19.178  3.966  4.243#

df <- data.frame(seq=c(rep("BS",7),rep("EM",7)),
                 type=c(BS$BS,EM$EM),
                 ratio=c(BS$CHH,EM$CHH))
##head(df)如下:##
#seq   type  ratio#
#BS  5-UTR  7.423#
#BS    CDS  4.243#

##画图
ggplot(df,aes(x=type,y=ratio,fill=seq)) +
  geom_bar(stat="identity",color="black",position=position_dodge(),alpha=0.5) +
  ylab("Mean methylation ratio (%)") +
  scale_y_continuous(expand=c(0,0),limits=c(0,NA)) +
  theme(panel.background=element_blank(),
        legend.position=c(0.1, 0.9),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.title.x=element_blank()) 