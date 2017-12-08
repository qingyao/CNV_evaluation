setwd('/Users/pgweb/arraydata/aroma/EvaluationPack/')
library(R.utils)
library(reshape2)
library(ggplot2)
lmds <- c(0.3,0.5,1)
gps <- c('1e4','5e4','1e5')
conds <- length(cid)*length(lmds)*length(gps)
segCompare <- data.frame()
for (cid in cids){
  for (lmd in lmds){
    for (gp in gps){
      before_lasso <- 0
      before_lasso <- countLines(sprintf('test_data/combine_series/%s/segments,cn,5_sdundo_1.tsv',cid))[1]
      after_lasso <- 0
      after_lasso <- countLines(sprintf('test_data/combine_series/%s/segments,cn,5_sdundo_1,lasso%s_%s.tsv',cid,lmd,gp))
      segCompare <- rbind(segCompare,data.frame(cid,lmd,gp,before_lasso,after_lasso))
    }
  }
}
gp_dict=c('1e4'='10kb','5e4'='50kb','1e5'='100kb')
segCompare$gp<- gp_dict[segCompare$gp]
segCompare$lmd_gp <- paste(segCompare$lmd,segCompare$gp,sep=',')

plf_search <- read.table('select_3_sample_platform.txt',header=T,stringsAsFactors = F)
platforms <- sapply(cids,function(x) plf_search[plf_search$sample==x,3])
segCompare$platform <- platforms[segCompare$cid]
pdata <- segCompare[,c(1,4:7)]
pdata <- melt(pdata)
pdata$lmd_gp<-factor(pdata$lmd_gp,levels=unique(pdata$lmd_gp))
ggplot(pdata, aes(x=variable, y=value, group=cid)) +
  geom_point(aes(colour=platform), size=3, position=position_dodge(width=0.1)) +
  geom_line(size=0.5, alpha=0.5, position=position_dodge(width=0.1)) +
  xlab('lambda, gapSize') +
  ylab('Number of segments') +
  scale_colour_manual(values=c("#009E73", "#D55E00")) + 
  scale_y_log10(limits=c(10,10000),breaks=c(10,31,100,316,1000,3162,1e4)) +
  theme_bw() + 
  facet_wrap(~lmd_gp)
ggsave('test_data/1dlasso_seg_reduction.pdf')
