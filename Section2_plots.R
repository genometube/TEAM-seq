library(ggplot2)
library(tidyr)
setwd("C:/Users/fangyitong/Desktop/2021Research/LIANTI/tables")

########## Barplot mC ratio of functional elements ########## 
all <- readxl::read_xlsx("AAA 10ng_1ng_100pg_10pg_NEB.xlsx")

df <- data.frame(sample=c(rep("NEB",9),rep("100pg",9),rep("10pg",9)),
                 element=rep(all$Element,9*3),
                 ratio=c(all$nebCHH,all$`100pgCHH`,all$`10pgCHH`))

ggplot(df,aes(x=element,y=ratio,fill=sample)) +
  geom_bar(stat="identity",color="black",position=position_dodge(),alpha=0.8) +
  ylab("Mean methylation ratio (%)") +
  scale_fill_manual(values=c("#D55E00","#E69F00","#999999")) + #"#56B4E9","#0072B2"
  scale_y_continuous(expand=c(0,0),limits=c(0,NA)) +
  theme(panel.background=element_blank(),
        #legend.position=c(0.1, 0.9),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.line=element_line(colour="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=60,hjust=1))