rootdir <- '/Volumes/arraymapMirror/arraymap/hg19'
series <- list.files(rootdir)
series <- series[grep('GSE',series)]
quality <- data.frame()
for (s in series) {
  f <- file.path(rootdir,s,'CNARA_output.txt')
  if (file.exists(f)){
    stable <- read.table(f,header = T,sep = '\t')
    quality <- rbind(quality, stable[,c(1,9,10)])
  }
}
quality <- quality[!duplicated(quality$Sample),]

setwd('/Users/pgweb/arraydata/aroma/EvaluationPack')
source('stepFilter/helper.R')
ev_summ <- read.table('output/sample_evaluation_summary.tsv',header = T,stringsAsFactors = F)
ev_summ <- ev_summ[ev_summ$Array %in% quality$Sample,]
ev_summ$quality <- quality$Good_Poor[match(ev_summ$Array,quality$Sample)]
ev_summ$Qclass <- quality$Case_Diagnosis[match(ev_summ$Array,quality$Sample)]
levels(ev_summ$Qclass) <- sapply(levels(ev_summ$Qclass),function(x)substring(x,1,1))
ev_summ$logCNsegments <- log(ev_summ$CNsegments)
ev_summ$logmax1cn <- log(ev_summ$max1_cn)
ev_summ$logmax2cn <- log(ev_summ$max2_cn)
library(reshape2)
library(ggplot2)
ncol(ev_summ)
mev <- melt(ev_summ[,c(2,3,4,9,21:22,25)])
p <- ggplot(mev,aes(x=quality,y=value))
p + geom_violin() + facet_grid(variable~platform,scales = 'free') + ggtitle( paste(levels(quality$Case_Diagnosis),collapse = '\n'))
ggsave('output/Rplot2.pdf',height = 10,width = 15)

p <- ggplot(mev,aes(x=Qclass,y=value))
p + geom_violin() + facet_grid(variable~platform,scales = 'free') + ggtitle( paste(levels(quality$Case_Diagnosis),collapse = '\n'))
ggsave('output/Rplot3.pdf',height = 10,width = 15)

p <- ggplot(mev,aes(value))
p + geom_histogram() + facet_grid(quality~variable,scales = 'free') 
p + geom_histogram() + facet_wrap(~platform,scales = 'free') 

plot(table(ev_summ[ev_summ$CNsegments > 1000,3 ]))

hist(ev_summ$CNsegments[ev_summ$platform=="GPL16131"])
boxplot(CNsegments~as.numeric(quality),data=ev_summ)
plot(as.numeric(ev_summ$quality),(ev_summ$CNsegments))
table(ev_summ[,ncol(ev_summ)])
plot(as.numeric(ev_summ$quality))

######################################
#analyze the lambdas for all samples#
######################################

fyles<-list.files('output/')
fyles<- fyles[grep('minLambda',fyles)]
minlambda <- data.frame()
for (f in fyles) {
  f <- file.path('output',f)
  print(f)
  stable <- read.table(f,sep = '\t',stringsAsFactors = F)
  minlambda <- rbind(minlambda, stable)
}
sum(minlambda$V2%in%quality$Sample) #35867 #39831
minlambda$quality <- quality$Case_Diagnosis[match(minlambda$V2,quality$Sample)]
smallminlambda <- minlambda[minlambda$V3<=0.5,]
bigminlambda <- minlambda[minlambda$V3>0.5,]

pdf('output/MinLambdaAnalysis/big_samples.pdf')
plot(bigminlambda$V3,log(bigminlambda$V4),col=bigminlambda$quality)
legend(2,9.3,unique(bigminlambda$quality),col=unique(bigminlambda$quality),pch=1)
dev.off()

pdf('output/MinLambdaAnalysis/overall_samples.pdf')
plot(log(minlambda$V3),log(minlambda$V4),col=minlambda$quality,ylab = 'log_noseg',xlab='minLambda')
legend(1,6.3,unique(minlambda$quality),col=unique(minlambda$quality),pch=1)
dev.off()

plot(smallminlambda$V3,log(smallminlambda$V4),col=smallminlambda$quality,ylab = 'log_noseg',xlab='minLambda')

colnames(smallminlambda)[3] <- 'min.lambda'
smallminlambda$quality <- swr(smallminlambda$quality)
ggplot(smallminlambda,aes(x=log(min.lambda)))+geom_histogram()+facet_wrap(~quality)
dir.create('../output/MinLambdaAnalysis')
ggsave('output/MinLambdaAnalysis/smalllambdasamples.pdf')

colnames(bigminlambda)[3] <- 'min.lambda'
bigminlambda$quality <- swr(bigminlambda$quality)
ggplot(bigminlambda,aes(x=log(min.lambda)))+geom_histogram()+facet_wrap(~quality)
ggsave('output/MinLambdaAnalysis/biglambdasamples.pdf')

boxplot(log(minlambda$V4),log(minlambda$V10))
lmd0.3 <- apply(minlambda,1,function(x){
  if (as.numeric(x[10]) != 0){
    c(x[c(1,2,4)],0.3,x[10])
  }else{
    tmp <- 0.3/as.numeric(x[3])
    idx <- which.min(abs(as.numeric(tmp)-c(1,2,10,100,1000))) ## idx that is closest to tmp
    c(x[c(1,2,4)],as.numeric(x[3])*(c(1,2,10,100,1000)[idx]) ,x[idx+4])
  }
})
lmd0.3 <- as.data.frame(t(lmd0.3))
lmd0.3$group <- cut(as.numeric(lmd0.3[,3]),c(0,200,500,1000,2000,5000,10000,100000))
table(lmd0.3$group)
levels(lmd0.3$group) <- paste(levels(lmd0.3$group), paste(table(lmd0.3$group),'samples'),sep='\n')
cols <- c(6)
lmd0.3[,-cols] <- apply(lmd0.3[,-cols],2,function(x) as.character(x)) 
colnames(lmd0.3)[-cols] <- c('series','array','original','lambda','new')
lmd0.3$quality <- quality$Good_Poor[match(lmd0.3$array,quality$Sample)]
mlmd <- melt(lmd0.3,id.vars = c('series','array','lambda','group','quality'),value.name = 'No.segments')
mlmd$No.segments <- as.numeric(mlmd$No.segments)
ggplot(mlmd,aes(x=variable,y=log(No.segments))) +geom_violin() +facet_grid(quality~group)+theme(axis.text.x = element_text(angle=45),axis.title.x = element_blank())
ggsave('output/MinLambdaAnalysis/allsamplesrange.pdf')
ggplot(mlmd[mlmd$lambda=='0.3',],aes(x=variable,y=log(No.segments))) +geom_violin() +facet_grid(quality~group)+theme(axis.text.x = element_text(angle=45),axis.title.x = element_blank())
ggsave('output/MinLambdaAnalysis/smalllambda_samplesrange.pdf')
