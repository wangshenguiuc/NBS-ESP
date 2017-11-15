library(survival)
#install.packages('ggfortify', dependencies=)
library(ggfortify)
library(ggplot2)
library(plyr)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

setwd("E:\\swang141\\project\\SupervisedSubgraph\\Sheng\\src")

args = commandArgs(trailingOnly=TRUE)

sur_data = args[1]
outpout_file = args[2]
nclst = as.numeric(args[3])
#KIPAN_6
#sur_data = '../output/subtype/KIPAN/KIPAN_6_0.5_3_4_2_50_1ME_Network_based.txt';
#outpout_file = '../paper/Nat.Genetic/figure/km_plot//KIPAN_6.pdf';
#nclst = 3


OV <- read.table(file=sur_data,header=TRUE,sep="\t",row.names=NULL)

tim <- 120
exp<-survdiff(Surv(OV$time, OV$event) ~ factor(OV$cluster))
nclst <- length(exp$n)
pv<-1 - pchisq(exp$chisq, nclst - 1)
pv_new<-prettyNum(pv, digits=3, width=4, format="fg")

OV <- read.table(file=sur_data,header=TRUE,sep="\t",row.names=NULL)
zselect60 <- which(OV$time > tim)
OV$time[zselect60] <- tim
OV$event[zselect60] <- 0
#file_name <- paste(outpout_file,tim,'_our.pdf',sep="")
#jpeg(file_name, width = 4, height = 4, units = 'in', res = 300)

#print(file_name)

fitMeta <- survfit(Surv(OV$time, OV$event) ~ factor(OV$cluster))
autoplot(fitMeta,conf.int = FALSE,surv.size=1.5,censor.size=3)+
  scale_y_continuous(labels = scales::percent,breaks=seq(0, 1, by=0.2),limits = c(0,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="top", legend.direction="horizontal") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.text=element_text(size=15))+
  theme(axis.text.y = element_text(size=15))+ 
  theme(axis.text.x = element_text(size=15))+ 
  theme(axis.title.x = element_text(size=20))+ 
  theme(axis.title.y = element_text(size=20))+ 
  theme(plot.title = element_text(size =20))+
  scale_x_continuous(breaks=seq(0, 60, by=10),limits = c(0,60))+
  #expand_limits(x = 0, y = 0)+
  #scale_y_continuous(labels = 'percent',breaks=seq(0, 100, by=10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title=element_blank())+
  scale_color_hue(labels = paste(levels(as.factor(OV$cluster))," (",table(OV$cluster),')',sep="")) +
  labs(title = paste('p=',pv_new,sep=""), x ="Time (months)", y = "Survival probability")+
  theme(plot.margin=unit(c(0,1,0,0),"cm"))
ggsave(outpout_file)
#plot(fitMeta,col= rainbow(nclst),mark.time = T,mark.color='black',xlab='Time (months)',ylab="Survival probablity",main = paste('Log-rank P =',pv_new,sep=""),cex.lab=1.3, cex.axis=0.8, cex.main=2, cex.sub=1.5)
#llgened.label <- paste(levels(as.factor(OV$our_cluster))," (",table(OV$our_cluster),')',sep="")
#legend(x=55,y=0.99,llgened.label,col=rainbow(nclst),lty=1,lwd=3,cex=1.0, pch=1, pt.cex = 1)
#xx<-dev.off()


