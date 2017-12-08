stepwiseDifference <- function(gapSize,segmentFile){
  seg <- read.table(segmentFile,stringsAsFactors = F,header = T)
  print(nrow(seg))
  accumulatePos <- 0
  newseg <-data.frame()
  for (row in 1:nrow(seg)){
    segLen <- seg[row,4]-seg[row,3]
    if ( segLen> gapSize) {
      accumulatePos <- accumulatePos + segLen
      newseg<-rbind(newseg,cbind(accumulatePos,seg[row,c(1:6)]))
    }
  }
  return(newseg)
}

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