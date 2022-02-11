########## Load libraries ########## 
library(ggplot2)
library(dplyr)
library(data.table)
library(VennDiagram)
library(eulerr)
library(vecsets)
setwd("C:/Users/fangyitong/Desktop/2021Research/LIANTI/tables")

########## Prepare files ########## 
A1 <- fread("A1chr1.csv.gz",header=FALSE,sep="\t")
colnames(A1) <- c("chr","pos","Ctype","mc","cover")
WGBS <- fread("WGBSchr1.csv.gz",header=FALSE,sep="\t")
colnames(WGBS) <- c("chr","pos","Ctype","mc","nmc")

########## mC ratio density plot ########## 
p_B3 <- B3 %>% filter(cover>=3) %>% mutate(ratio=mc/cover) %>% select(ratio)
ggplot(p_B1) +
  geom_density(mapping=aes(x=ratio,colour="TEAM20pg_1"),alpha=0.5) +
  geom_density(p_B2,mapping=aes(x=ratio,color="TEAM20pg_2"),alpha=0.5) +
  geom_density(p_B3,mapping=aes(x=ratio,color="TEAM20pg_3"),alpha=0.5) +
  ylim(0,8) +
  xlab("Methylation ratio") +
  ylab("Density") +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.position=c(0.9,7),
        axis.line=element_line(colour="black"))

########## binned mC ratio vs frequency plot ########## 
allC <- B1 %>% filter(cover>=3) %>% mutate(ratio=mc/cover) %>% select(ratio) %>%
  mutate(bin=cut(ratio,breaks=seq(0,1,by=0.1),include.lowest=TRUE,labels=seq(10,100,by=10))) %>% 
  count(bin) %>% mutate(freq=n/sum(n)*100)
CG <- B1 %>% filter(cover>=3) %>% filter(grepl(".CG.",seq)) %>% mutate(ratio=mc/cover) %>% 
  select(ratio) %>% mutate(bin=cut(ratio,breaks=seq(0,1,by=0.1),include.lowest=TRUE,labels=seq(10,100,by=10))) %>% 
  count(bin) %>% mutate(freq=n/sum(n)*100)
CHG <- B1 %>% filter(cover>=3) %>% filter(grepl(".CAG|.CTG|.CCG",seq)) %>% mutate(ratio=mc/cover) %>% 
  select(ratio) %>% mutate(bin=cut(ratio,breaks=seq(0,1,by=0.1),include.lowest=TRUE,labels=seq(10,100,by=10))) %>% 
  count(bin) %>% mutate(freq=n/sum(n)*100)
CHH <- B1 %>% filter(cover>=3) %>% filter(grepl(".CAA|.CAT|.CAC|.CTA|.CTT|.CTC|.CCA|.CCT|.CCC",seq)) %>% mutate(ratio=mc/cover) %>% 
  select(ratio) %>% mutate(bin=cut(ratio,breaks=seq(0,1,by=0.1),include.lowest=TRUE,labels=seq(10,100,by=10))) %>% 
  count(bin) %>% mutate(freq=n/sum(n)*100)

ggplot() +
  geom_point(allC,mapping=aes(x=bin,y=freq,group=1,color="C")) +
  geom_line(allC,mapping=aes(x=bin,y=freq,group=1,color="C")) +
  geom_point(CG,mapping=aes(x=bin,y=freq,group=1,color="CG")) +
  geom_line(CG,mapping=aes(x=bin,y=freq,group=1,color="CG")) +
  geom_point(CHG,mapping=aes(x=bin,y=freq,group=1,color="CHG")) +
  geom_line(CHG,mapping=aes(x=bin,y=freq,group=1,color="CHG")) +
  geom_point(CHH,mapping=aes(x=bin,y=freq,group=1,color="CHH")) +
  geom_line(CHH,mapping=aes(x=bin,y=freq,group=1,color="CHH")) +
  xlab("Methylation level (%)") +
  ylab("Percentage (%) in corresponding mC type") +
  theme_bw() +
  theme(panel.background=element_blank(),
        panel.border=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        legend.position = c(0.9, 0.85))

########## CONSISTENCY B1_B2_B3 ##########
B1 <- fread("B1chr1.csv.gz",header=FALSE,sep="\t")
colnames(B1) <- c("chr","pos","Ctype","mc","cover")
B2 <- fread("B2chr1.csv.gz",header=FALSE,sep="\t")
colnames(B2) <- c("chr","pos","Ctype","mc","cover")
B3 <- fread("B3chr1.csv.gz",header=FALSE,sep="\t")
colnames(B3) <- c("chr","pos","Ctype","mc","cover")

B123 <- merge(merge(B1,B2,by="pos"),B3,by="pos") %>%
  mutate(B1ratio=mc.x/cover.x,B2ratio=mc.y/cover.y,B3ratio=mc/cover) %>%
  select(pos,B1ratio,B2ratio,B3ratio) %>%
  mutate(B1=ifelse(B1ratio>0.5,1,0),B2=ifelse(B2ratio>0.5,1,0),B3=ifelse(B3ratio>0.5,1,0))

common123=dim(subset(subset(B123,B1==B2),B1==B3))[1]
common12=dim(subset(B123,B1==B2))[1]
common13=dim(subset(B123,B1==B3))[1]
common23=dim(subset(B123,B2==B3))[1]
v <- draw.triple.venn(#area1=3845985,area2=3845985,area3=3845985,
                      area1=100,area2=100,area3=100,
                      #n12=common12,n23=common23,n13=common13,n123=common123,
                      n12=95.14,n13=95.03,n23=95.3,n123=92.73,
                      category=c("B1","B2","B3"),
                      fill = c("dodgerblue","goldenrod1","ForestGreen"),
                      lty="blank", print.mode="percent",
                      cat.pos=c(315,45,180),cat.dist=c(0.03,0.03,0.03))
common123=round(dim(subset(subset(B123,B1==B2),B1==B3))[1]/dim(B123)[1]*100,2)
common12=round(dim(subset(B123,B1==B2))[1]/dim(B123)[1]*100,2)
common13=round(dim(subset(B123,B1==B3))[1]/dim(B123)[1]*100,2)
common23=round(dim(subset(B123,B2==B3))[1]/dim(B123)[1]*100,2)
