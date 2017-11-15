library(survival)

args = commandArgs(trailingOnly=TRUE)

sur_data = args[1]

outpout_file = args[2]
print('outpout_file')
write_figure = as.numeric(args[3])
print('write_figure')


ts1 = as.numeric(args[4])
ts2 = as.numeric(args[5])
ts3 = as.numeric(args[6])
ts4 = as.numeric(args[7])
nclst = as.numeric(args[8])
OV <- read.table(file=sur_data,header=TRUE,sep="\t",row.names=NULL)

time_threshold <- c(ts1,ts2,ts3,ts4)
pv_our_ret <- c(0,0,0,0)
chi_our_ret <- c(0,0,0,0)
pv_nbs_ret <- c(0,0,0,0)
chi_nbs_ret <- c(0,0,0,0)
for (i in 4:length(time_threshold))
{
  tim <- time_threshold[i]
  OV <- read.table(file=sur_data,header=TRUE,sep="\t",row.names=NULL)
  zselect60 <- which(OV$time > tim)
  OV$time[zselect60] <- tim
  OV$event[zselect60] <- 0
  
  exp<-survdiff(Surv(OV$time, OV$event) ~ factor(OV$our_cluster))
  nclst <- length(exp$n)
  pv<-1 - pchisq(exp$chisq, nclst - 1)
  pv_new<-prettyNum(pv, digits=3, width=4, format="fg")
  print(as.numeric(pv_new))
  print(tim)
  pv_our_ret[i] <- pv_new
  chi_our_ret[i] <- exp$chisq
  if ((write_figure==1))
  {
    file_name <- paste(outpout_file,tim,'_our.jpg',sep="")
    print(file_name)
    jpeg(file_name, width = 4, height = 4, units = 'in', res = 300)
    OV <- read.table(file=sur_data,header=TRUE,sep="\t",row.names=NULL)
    zselect60 <- which(OV$time > tim)
    OV$time[zselect60] <- tim
    OV$event[zselect60] <- 0
    fitMeta <- survfit(Surv(OV$time, OV$event) ~ factor(OV$our_cluster))
    plot(fitMeta,col= rainbow(nclst),xlab='Time (months)',ylab="Probablity of survival",main = paste('Log-rank P =',pv_new,sep=""),cex.lab=1.3, cex.axis=0.8, cex.main=2, cex.sub=1.5)
    
    llgened.label <- paste('Subtype ',levels(as.factor(OV$our_cluster))," (",table(OV$our_cluster),')',sep="")
    legend(x=10,y=0.2,llgened.label,col=rainbow(nclst),lty=1,lwd=3,cex=0.5, pch=1, pt.cex = 0.5)
    xx<-dev.off()
  }
  
}
for (i in 4:length(time_threshold))
{
  tim <- time_threshold[i]
  OV <- read.table(file=sur_data,header=TRUE,sep="\t",row.names=NULL)
  zselect60 <- which(OV$time > tim)
  OV$time[zselect60] <- tim
  OV$event[zselect60] <- 0
  
  exp<-survdiff(Surv(OV$time, OV$event) ~ factor(OV$nbs_cluster))
  nclst <- length(exp$n)
  pv<-1 - pchisq(exp$chisq, nclst - 1)
  pv_new<-prettyNum(pv, digits=3, width=4, format="fg")
  print(as.numeric(pv_new))
  pv_nbs_ret[i] <- pv_new
  chi_nbs_ret[i] <- exp$chisq
  if ((write_figure==1) & (pv_new<0.01) & FALSE)
  {
    OV <- read.table(file=sur_data,header=TRUE,sep="\t",row.names=NULL)
    zselect60 <- which(OV$time > tim)
    OV$time[zselect60] <- tim
    OV$event[zselect60] <- 0
    file_name <- paste(outpout_file,tim,'_nbs.jpg',sep="")
    jpeg(file_name)
    
    fitMeta <- survfit(Surv(OV$time, OV$event) ~ factor(OV$nbs_cluster))
    plot(fitMeta,col= rainbow(nclst),xlab='Time(months)',ylab="Probablity of survival")
    title(main = paste('nbs p=',pv_new,sep=""))
    llgened.label <- paste('Subtype ',levels(as.factor(OV$our_cluster))," (",table(OV$our_cluster),')',sep="")
    legend(x=5,y=0.3,llgened.label,col=rainbow(nclst),lty=1,lwd=3)
    xx<-dev.off()
  }
  
}
res<-cbind(pv_our_ret,chi_our_ret,pv_nbs_ret,chi_nbs_ret)
options(scipen=10)
write(res, file = outpout_file,append = FALSE, sep = " ")
options(scipen=0)

