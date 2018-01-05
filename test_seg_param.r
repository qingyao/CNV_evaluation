##########################################
#test for the parameters for segmentation#
##########################################

##First test: 
# for all 9 platforms
# try min.width 2,5 (no change)
# try smooth.sd , no change
# try undo.prune, too slow
# run undo.sd 1 and 2, smooth.region 5 and 10, in loop
# generate log.txt with segment number (followed by significant ones p<0.1) per chromosome and process time
# for GPL16131, GPL6801, too many segments. Other platforms are fine. Mostly smooth.region 5 and undo.sd 1 is the best option.

##Second test:
# test for GPL6801, GPL16131
# try undo.sd 3 and 4, 
# they removed important segments, not desirable
# so will merge closeby segments

library(DNAcopy)
library(reshape2)
seriesName='GSE15264'
seriesName='GSE17359'
seriesName='combine_series'
setwd('/Users/pgweb/arraydata/aroma/EvaluationPack/')
workingdir = 'test_data/'
cids=list.files(file.path(workingdir,seriesName))
cids=cids[-grep("_",cids)]
cid = 'GSM381297'
cid = 'GSM433918'
chipType='GenomeWideSNP_6'
workingdir='/Users/pgweb/arraydata/aroma/hg19/processed'
sourcedir='/Users/pgweb/arraydata/aroma/AromaPack/'
filedir = '/Volumes/arraymapMirror/arraymap/hg19/'
chrs = 1:23
# chrname = chrs
# minw <- 2
# smrg <- 5
# smsd <- 1
# undosplit <- "sdundo"
# undopr <- NA
# undosd  <- 4
cnsegPerArray <- function(workingdir,seriesName, cid, chrs, minw, smrg, smsd, undosplit, undopr, undosd){
  # dir.create(file.path(workingdir,seriesName,cid), showWarnings = FALSE)
  fn <- file.path(workingdir,seriesName,cid, sprintf('segments,cn,%s_%s_%s.tsv',smrg, undosplit,undosd))
  fn <- file.path(workingdir,seriesName,cid, sprintf('segments,cn,sd,%s_%s_%s.tsv',smrg, undosplit,undosd))
  cat("sample_id","chromosome","start", "end", "value", "probes\n",sep="\t",file=fn,append = F)
  cat("sample_id","chromosome","start", "end", "probes", "seg.mean","seg.sd\n",sep="\t",file=fn,append = F)
  t <- 0 #
  logf <- file.path(workingdir,seriesName,cid, 'log.txt')
  cat(minw,smrg, smsd, undosplit, undopr,undosd,'',sep = '\t',file = logf,append = TRUE)
  for (chrname in 1:23){
    data<- read.table(sprintf('%s/%s/%s/probes,cn,chr%d.tsv',workingdir,seriesName,cid,chrname),header=T)
    cna1 <- CNA(genomdat=data$VALUE,
                chrom=rep(chrname, length(data$VALUE)),
                maploc=data$BASEPOS,
                data.type="logratio",
                sampleid=cid)
    t <- t + system.time(smoo1 <- smooth.CNA(cna1, smooth.region=smrg,smooth.SD.scale=smsd))[3]#
    message(paste("Processed CN segmentation for sample:",cid,"Chr:",chrname))
                        
    t <- t + system.time(seg1 <- segment(smoo1, min.width= minw, verbose=0, undo.splits = undosplit, undo.prune = undopr,undo.SD = undosd))[3]#
    ss1 <- segments.summary(seg1)[c(1:4,6,5)] #standard output format, without sd
    # ss1 <- segments.summary(seg1)[c(1:7)] # for test 3
    write.table(ss1, file=fn, sep="\t", quote=FALSE,
                append = T, row.names=FALSE, col.names=F)
    
    cat(nrow(ss1),sum(segments.p(seg1)$pval < 0.05,na.rm = T),'',sep = '\t',file = logf,append = TRUE)
  }
  cat(t, '\n',file = logf,append = TRUE,sep = '')
}
for (cid in cids){
  for (undosplit in c("none")){ #"sdundo",
    for (minw in c( 5)){
      # for (smrg in c(5,10)){
      for (smrg in c(5)){  ## test for GPL16131, GPL6801
        for (smsd in c(2) ){ 
          
          if (undosplit == "prune"){
            for (undopr in c(0.05, 0.1)){
              undosd <- NA
              cat(minw, smrg, smsd, undosplit, undopr, undosd)
              cnsegPerArray(workingdir,seriesName, cid, chrs, minw, smrg, smsd, undosplit, undopr, undosd)
              }
          }else if(undosplit=="none") {
            undosd <- NA
            undopr <- NA
            cat(minw, smrg, smsd, undosplit, undopr, undosd)
            cnsegPerArray(workingdir,seriesName, cid, chrs, minw, smrg, smsd, undosplit, undopr, undosd)
            
            }else{
          for (undosd in c(1,2)){
          # for (undosd in c(3,4)){ ## test for GPL16131, GPL6801
            undopr <- NA
            cat(minw, smrg, smsd, undosplit, undopr, undosd)
            cnsegPerArray(workingdir,seriesName, cid, chrs, minw, smrg, smsd, undosplit, undopr, undosd)
          }
          
          }
        }
      }
    }
  }
}

seg_df <- data.frame()
sig_seg_df <- data.frame()
for (cid in cids){
  logf <- file.path(workingdir,seriesName,cid, 'log.txt')
  output <- read.table(logf,stringsAsFactors = F)
  nr <- nrow(output) 
  cols <- 6+(1:23)*2-1
  seg_num <- apply(output[c(1:4,nr),cols], 1, sum) ##some tests in some samples in the middle rows
  # sig_seg_num <- apply(output[,cols+1], 1, sum)+23
  sig_seg_num <-seg_num-apply(output[c(1:4,nr),cols], 1, max)
  seg_df <- rbind(seg_df,seg_num)
  sig_seg_df <- rbind(sig_seg_df, sig_seg_num)
}

plf_search <- read.table('select_3_sample_platform.txt',header=T,stringsAsFactors = F)
colnames(seg_df) = c('5,1','5,2','10,1','10,2','none')
colnames(sig_seg_df) = colnames(seg_df)

platforms <- sapply(cids,function(x) plf_search[plf_search$sample==x,3])
seg_df$platform <- factor(platforms)
sig_seg_df$platform <- factor(platforms)

seg_df <- melt(seg_df,id.vars = 6)
sig_seg_df <- melt(sig_seg_df,id.vars = 6)
colnames(seg_df) <- c('platform','smoothRegion,undoSD','No.segments')
colnames(sig_seg_df) = colnames(seg_df)
seg_df$`smoothRegion,undoSD` <- factor(seg_df$`smoothRegion,undoSD`,levels = c('none','5,1','5,2','10,1','10,2'))
sig_seg_df$`smoothRegion,undoSD` <- factor(sig_seg_df$`smoothRegion,undoSD`,levels = c('none','5,1','5,2','10,1','10,2'))
p = ggplot(seg_df,aes(x=`smoothRegion,undoSD`,y=No.segments))
p = ggplot(sig_seg_df,aes(x=`smoothRegion,undoSD`,y=No.segments))
p+geom_violin()+ geom_jitter(height = 0, width = 0.1)+
  facet_wrap(~platform, scales = 'free') + ggtitle('total CN segments')
p+geom_violin()+ geom_jitter(height = 0, width = 0.1)+
  facet_wrap(~platform, scales = 'free') + ggtitle('total CN segments-maxchr')
ggsave('test_data/segment_test_plot2.pdf',width = 12, height = 7)
ggsave('test_data/segment-maxchr_test_plot2.pdf',width = 10, height = 6)

##Third test:
# merging
# function which input size limit for gaps between real segment, try 100kb, output real segment with sd,mean,accumluative length
stepwiseDifference <- function(gapSize,segmentFile){
  seg <- read.table(segmentFile,stringsAsFactors = F,header = T)
  print(nrow(seg))
  accumulatePos <- 0
  newseg <-data.frame()
  for (row in 1:nrow(seg)){
    segLen <- seg[row,4]-seg[row,3]
    if ( segLen> gapSize) {
      accumulatePos <- accumulatePos + segLen
       newseg<-rbind(newseg,cbind(accumulatePos,seg[row,c(2:7)]))
    }
  }
  return(newseg)
} 
sampleName <- 'GSM1146803'
segStep <- stepwiseDifference(1e+5,sprintf('test_data/combine_series/%s/segments,cn,sd,10_sdundo_2.tsv',sampleName))
p=ggplot()
p+geom_step(data=segStep,mapping=aes(x=accumulatePos, y=seg.mean),color='black')+
  geom_step(data=segStep,mapping=aes(x=accumulatePos, y=seg.mean+seg.sd),color='red') +
  geom_step(data=segStep,mapping=aes(x=accumulatePos, y=seg.mean-seg.sd),color='red')
ggsave(sprintf('test_data/combine_series/%s/segment_step_plot.pdf',sampleName))

statsSmallSeg <- function(gapSizes,segmentFile){
  seg <- read.table(segmentFile,stringsAsFactors = F,header = T)
  smallSeg <- lapply(gapSizes, function(gapSize){
    newseg <-data.frame()
    for (row in 1:nrow(seg)){
      segLen <- seg[row,4]-seg[row,3]
      if ( segLen< gapSize) {
        newseg<-rbind(newseg,data.frame(length=seg[row,4]-seg[row,3],probes=seg[row,5]))
      }
    }
    return (newseg)
  })
  names(smallSeg) <- gapSizes
  return(smallSeg)
} 

sampleName <- 'GSM1146803'
a=statsSmallSeg(c(1e4,5e4,1e5),sprintf('test_data/combine_series/%s/segments,cn,sd,10_sdundo_2.tsv',sampleName))
lengt<-lapply(a,function(x)x$length)
l.melt <- melt(lengt,na.rm = T) #melt list
colnames(l.melt) <- c('segLength','gapSize')
l.melt$gapSize<-factor(l.melt$gapSize,levels=unique(l.melt$gapSize))
ggplot(l.melt,aes(x=gapSize,y=segLength)) + geom_violin() + 
  # geom_jitter(shape=16,position=position_jitter(0.2)) +
  scale_y_log10() + 
  labs(title="size threshold for small: number of small segments",subtitle=paste(c("10kb", "50kb", "100kb"),summary(l.melt$gapSize),sep = ": ", collapse = '\n'))

ggsave(sprintf('test_data/%s_small_segment_len.pdf',sampleName),width = 4,height = 4)

probenumber<-lapply(a,function(x)x$probe)
p.melt <- melt(probenumber,na.rm = T) #melt list
colnames(p.melt) <- c('probeNumber','gapSize')
p.melt$gapSize<-factor(l.melt$gapSize,levels=unique(p.melt$gapSize))
ggplot(p.melt,aes(x=gapSize,y=probeNumber)) + geom_violin() + 
  # geom_jitter(shape=16,position=position_jitter(0.2)) + 
  scale_y_log10() + 
  labs(title="size threshold for small: number of small segments",subtitle=paste(c("10kb", "50kb", "100kb"),summary(l.melt$gapSize),sep = ": ", collapse = '\n'))

ggsave(sprintf('test_data/%s_small_segment_probeNr.pdf',sampleName),width = 4,height = 4)

library(ggplot2)
library(readr)
Summary <- read_delim("/Users/pgweb/arraydata/aroma/EvaluationPack/sample_evalutaion_summary.tsv", 
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)
S <- Summary[Summary$Array %in% cids,]
colnames(Summary)
p = ggplot(S)
p + geom_point(aes(x=CNsegments,y=fracbsegments)) + facet_wrap(~platform) + scale_x_log10() + scale_y_log10()

p + geom_point(aes(x=CNsegments,y=kurtosis_prob)) + facet_wrap(~platform) + scale_x_log10() + scale_y_log10()

p + geom_point(aes(x=CNsegments,y=kurtosis_seg)) + facet_wrap(~platform) + scale_x_log10() + scale_y_log10()

p + geom_point(aes(x=kurtosis_seg,y=kurtosis_prob)) + facet_wrap(~platform) + scale_x_log10() + scale_y_log10()

