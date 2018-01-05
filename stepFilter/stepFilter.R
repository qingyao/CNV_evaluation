setwd('~/arraydata/aroma/EvaluationPack/')
source('seg_stats.r')
library(genlasso)

cids= c("GSM1704973","GSM1704975","GSM1704979","GSM1146803","GSM337641", "GSM381297", "GSM433918")
options("scipen"=100, "digits"=4)
for (cid in cids){
  for (gp in c('1e4','5e4','1e5')){ 
    segStep = stepwiseDifference(as.numeric(gp),sprintf('test_data/combine_series/%s/segments,cn,5_sdundo_1.tsv',cid))
    segStep <- segStep[,c(2:ncol(segStep),1)]
    print (nrow(segStep))
    
    x = segStep$value
    for (lmd in c(0.3)){#0.5, 1
      ### 1dlasso
      
      idx <- calculateLasso(x,lmd)[[2]]
      wm <- calculateWeightedMean(idx,segStep)
      newseg <- data.frame()
      idx <- c(idx,nrow(segStep)) ##each idx is the end of segment and add the last row of all segments
      for (i in 1:length(idx)){
        if (i==1){
          tmpseg <- vector()
          tmpseg <-segStep[idx[i],]
          tmpseg[,3] <- segStep[1,3]
          tmpseg[,5] <- wm[i]
          tmpseg[,6] <- sum(segStep[c(1:idx[i]),6])
          newseg<- rbind(newseg, tmpseg)
        }else if (segStep[idx[i],2] == segStep[idx[i-1]+1,2]){
          tmpseg <- vector()
          tmpseg<-segStep[idx[i],]
          tmpseg[,3] <- segStep[idx[i-1]+1,3] #start of last row +1's start
          tmpseg[,5] <- wm[i] # weighted mean 
          tmpseg[,6] <- sum(segStep[c((idx[i-1]+1):idx[i]),6]) # sum of probe numbers since last row +1
          newseg<- rbind(newseg, tmpseg)
        }else{
          if (idx[i-1]+1 != idx[i]) {
            previousChr <- which(segStep[,2] == segStep[idx[i-1],2])
            endChridx <-previousChr[length(previousChr)]
            tmpseg <- vector()
            tmpseg<-segStep[endChridx,]
            tmpseg[,3] <- segStep[idx[i-1]+1,3] 
            tmpseg[,5] <- wm[i]
            tmpseg[,6] <- sum(segStep[c((idx[i-1]+1):endChridx),6])
            newseg<- rbind(newseg, tmpseg)}
          
          tmpseg <- vector()
          tmpseg <- segStep[idx[i],] 
          tmpseg[,3] <- segStep[endChridx+1,3]
          tmpseg[,5] <- wm[i]
          tmpseg[,6] <- sum(segStep[c((endChridx+1):idx[i]),6])
          newseg <- rbind(newseg, tmpseg)
        }
        
      }
      write.table(newseg, sprintf('test_data/combine_series/%s/segments,cn,5_sdundo_1,lasso%s_%s.tsv',cid,lmd,gp), sep="\t", quote=FALSE,row.names=FALSE)
    }
    }
  }

  

for (cid in cids){
  for (gp in c('1e4','1e5')){
    for (lmd in c(0.5, 1)){
      segStep = stepwiseDifference(as.numeric(gp),sprintf('test_data/combine_series/%s/segments,cn,5_sdundo_1.tsv',cid))
      print (nrow(segStep))
      x = segStep$value
      beta1<- calculateLasso(x,lmd)[[1]]
      idx <- calculateLasso(x,lmd)[[2]]
      wm <- calculateWeightedMean(idx,segStep)
      plotLasso(cid,lmd,gp,idx,segStep,wm,beta1)
    }
  }
}

calculateLasso <- function(segmentData, lmd) {
  set.seed(1)
  segmentData<-segStep[segStep$chromosome==23,]
  segmentData<-segStep
  x <- segmentData$value
  pos <- segmentData$accumulatePos
  out = fusedlasso1d(x,pos = pos)
  out = fusedlasso1d(x)
  beta1 = coef(out, lambda=lmd)$beta
  # beta2 = softthresh(out, lambda=1.5, gamma=0.1)
  # plot(out, lambda=0.2)
  # coef(out)
  # sum(diff(beta1) > 1e-4)
  # segStep[which(diff(beta1) > 1e-4),]
  # segmentValue[which(diff(beta1) > 1e-4)]
  # beta1[which(diff(beta1) > 1e-4)]
  
  ## calculate weighted mean in each segment
  beta1 <- round(beta1,4)
  idx = which(abs(diff(beta1)) > 1e-4)
  return(list(beta1,idx))
}

calculateWeightedMean <- function(idx, segmentData) {
  lasti=1 
  wm=rep(0,length(idx)+1)
  for(i in 1:(length(idx)+1)){
    wm[i] <- ifelse(i!=length(idx)+1,
                    weighted.mean(segmentData$value[lasti:idx[i]],w=segmentData$end[lasti:idx[i]]-segmentData$start[lasti:idx[i]]),
                    weighted.mean(segmentData$value[lasti:nrow(segmentData)],w=segmentData$end[lasti:nrow(segmentData)]-segmentData$start[lasti:nrow(segmentData)]))
    
    lasti <- idx[i]+1}
  return (wm)
}


plotLasso <- function(cid, lmd, gp, idx, segmentData, wm, beta1) {
    pdf(sprintf('test_data/combine_series/%s_%s_%s.pdf',cid,lmd,gp))
    x <- segmentData$value
    plot(x)
    lines(beta1,col='blue')
    x1 <- do.call(c,apply(cbind(diff(c(0,idx,nrow(segStep))),wm),1, function(x) rep(x[2],x[1]))) ## create vector where all data points within segment are weighted mean
    lines(x1,col='red')
    abline(0.15,0,lty=2)
    abline(-0.15,0,lty=2)
    title(sprintf("%s: %s segments",cid,length(wm)))
    legend('topright',legend=c('point-wise beta1','weighted by len'),col = c('blue','red'),lty=c(1,1),cex = 0.6)
    dev.off()
    
    
    sprintf('test_data/combine_series/%s/segments,cn,5_sdundo_1.tsv',cid)
}
