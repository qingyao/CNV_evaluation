library(ggplot2)
library(gridExtra)
# set.seed (123)
# xvar <- c(rnorm (100, 50, 30), rnorm (100, 40, 10), rnorm (100, 70, 10))
# yvar <-   xvar + rnorm (length (xvar), 0, 20)
# myd <- data.frame (xvar, yvar)
# 
# pMain <- ggplot(myd,aes(x=xvar,y=yvar))+
#   stat_density2d(aes(fill=..level..), geom="polygon") +
#   coord_cartesian(c(0, 150), c(0, 150)) +
#   theme(legend.position = "none")
# 
# pTop <- ggplot(myd, aes(x = xvar)) + stat_density() +
#   coord_cartesian(c(0, 150))
# pRight <- ggplot(myd, aes(x = yvar)) + stat_density() + 
#   coord_flip(c(0, 150))
# pEmpty <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#   geom_blank() +
#   theme(axis.text = element_blank(),
#         axis.title = element_blank(),
#         line = element_blank(),
#         panel.background = element_blank())
# 
# grid.arrange(pTop, pEmpty, pMain, pRight,
#              ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))

library(readr)
library(plyr)
rootdir <- '/Users/pgweb/arraydata/aroma/hg19/processed/GSE22305/'
samples <- list.files(rootdir)
samples <- samples[grep('GSM',samples)]
for (sample in samples){
  pdf(file.path(rootdir,sample,sprintf('%s_densityplot.pdf',sample)))
  data <- data.frame()
  for (typeData in c("probe-level","segment-level")){
    if (typeData =="probe-level"){
      
      probef <- read_delim(file.path(rootdir,sample,'probes,fracb.tsv'),delim = '\t')
      probecn <- read_delim(file.path(rootdir,sample,'probes,cn.tsv'),delim = '\t')
      probecn <- probecn[which(probecn$PROBEID %in% probef$ID),]
      
      data <- data.frame(pos = 1:nrow(probef),fracb = probef[[4]], cn = probecn[[4]]) #probe
    } else if(typeData=="segment-level"){
      
      segf<- read_delim(file.path(rootdir,sample,'segments,fracb.tsv'),delim = '\t')
      segcn <- read_delim(file.path(rootdir,sample,'segments,cn.tsv'),delim = '\t')
      ls_newsegment <- getOverlapSegments(segcn,segf)
      segcn <- ls_newsegment[[1]]
      segf<-ls_newsegment[[2]]
      segcn$chrpos <- paste(segcn[[2]],round(as.numeric(segcn[[3]])/100000,1),round((as.numeric(segcn[[4]])-as.numeric(segcn[[3]]))/100000),sep=',')
      segf$chrpos <- paste(segf[[2]],round(as.numeric(segf[[3]])/100000,1),round((as.numeric(segf[[4]])-as.numeric(segf[[3]]))/100000),sep=',')
      
      segdata<- join(segcn,segf,by='chrpos')
      segdata<- segdata[!sapply(segdata[,9],is.null),]
      
      data <- data.frame(pos=1:nrow(segdata),fracb = as.numeric(segdata$fracB),cn=as.numeric(segdata[[5]])) #segment
    }
    
    pMain <- ggplot(data,aes(x=fracb,y=cn))+
      stat_density2d(aes(fill=..level..), geom="polygon") +
      coord_cartesian(c(0, 1), c(-2,2)) +
      theme(legend.position = "none")
    pTop <- ggplot(data, aes(x = fracb)) + stat_density() +
      coord_cartesian(c(0, 1))
    pRight <- ggplot(data, aes(x = cn)) + stat_density() + 
      coord_flip(c(-2, 2))
    pEmpty <- ggplot(data, aes(x = fracb, y = cn)) +
      geom_blank() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            line = element_blank(),
            panel.background = element_blank())
    
    # g <- arrangeGrob(pTop, pEmpty, pMain, pRight,
    # ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
    # ggsave(file = file.path(rootdir,sample,sprintf('%s_segment_densityplot.pdf',sample)),g)
    grid.arrange(pTop, pEmpty, pMain, pRight, top = typeData, ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
    
  }
  dev.off()
}

getOverlapSegments <-function(segment1,segment2){
  newsegment1 <- data.frame()
  newsegment2 <- data.frame()
  for (chr in 1:23){
    allsegpos <- sort(unique(c(segment1[which(segment1[,2]==chr),][[3]],segment2[which(segment2[,2]==chr),][[3]])))
    for (line in 1:sum(segment1[,2]==chr)){
      startpos <- as.numeric(segment1[line,3])
      endpos <- as.numeric(segment1[line,4])
      idx <- which(allsegpos > startpos & allsegpos < endpos) 
      if (length(idx) > 0){
        tmpseg <- cbind(matrix(rep(segment1[line,],length(idx)+1),byrow =T, ncol = ncol(segment1)))
        tmppos <- c(startpos, allsegpos[idx],endpos)
        for (tmpline in 1:nrow(tmpseg)){
          tmpseg [tmpline,3] <- tmppos[tmpline]
          tmpseg [tmpline,4] <- tmppos[tmpline+1]
        }
        tmpseg <-as.data.frame(tmpseg)
        colnames(tmpseg) <- colnames(segment1)
        newsegment1 <- rbind(newsegment1,tmpseg)
        
      }else{
        newsegment1 <- rbind(newsegment1,segment1[line,])
      }
    }
    for (line in 1:sum(segment2[,2]==chr)){
      startpos <- as.numeric(segment2[line,3])
      endpos <- as.numeric(segment2[line,4])
      idx <- which(allsegpos > startpos & allsegpos < endpos) 
      if (length(idx) > 0){
        tmpseg <- cbind(matrix(rep(segment2[line,],length(idx)+1),byrow =T, ncol = ncol(segment2)))
        tmppos <- c(startpos, allsegpos[idx],endpos)
        for (tmpline in 1:nrow(tmpseg)){
          tmpseg [tmpline,3] <- tmppos[tmpline]
          tmpseg [tmpline,4] <- tmppos[tmpline+1]
        }
        tmpseg <-as.data.frame(tmpseg)
        colnames(tmpseg) <- colnames(segment2)
        newsegment2 <- rbind(newsegment2,tmpseg)
        
      }else{
        newsegment2<- rbind(newsegment2,segment2[line,])
      }
    }
  }
  return(list(newsegment1,newsegment2))
}
  
