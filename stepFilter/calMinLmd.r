suppressWarnings(suppressMessages(library(genlasso)))
source('seg_stats.r')
options(digits =4)

getNoSeg <- function(lassoOut, lmd) {
  beta1 = coef(lassoOut, lambda=lmd)$beta
  beta1 <- round(beta1,4)
  idx = which(abs(diff(beta1)) > 1e-4)
  return(length(idx))
}
  
CalculateMinLambda <- function(workingdir,seriesName,arrayName,gp){
  segStep <- stepwiseDifference(as.numeric(gp),file.path(workingdir,seriesName,arrayName,'segments,cn.tsv'))
  segStep <- segStep[,c(2:ncol(segStep))]

  x = segStep[,5]
  noseg <- length(x)
  set.seed(1)
  out <- fusedlasso1d(x)
  minl <- min(out$lambda)
  minl <- signif(minl,digits = 4)
  i0 <- getNoSeg(out,minl)
  i1 <- getNoSeg(out,minl*2)
  i2 <- getNoSeg(out,minl*10)
  i3 <- getNoSeg(out,minl*100)
  i4 <- getNoSeg(out,minl*1000)
  if (minl * 1000 < 0.3) {
    i5 <- getNoSeg(out,0.3)
  }else{
    i5 <- 0
  }
  
  cat(paste(c(seriesName,arrayName,minl,noseg,i0,i1,i2,i3,i4,i5), collapse ='\t'),'\n',sep = '',append = T, file = 'output/minLambda.tsv')
}

CalibrateLambda <- function(workingdir,seriesName,arrayName,gp){
  segStep <- stepwiseDifference(as.numeric(gp),file.path(workingdir,seriesName,arrayName,'segments,cn.tsv'))
  segStep <- segStep[,c(2:ncol(segStep))]
  
  x = segStep[,5]
  noseg <- length(x)
  set.seed(1)
  out <- fusedlasso1d(x)
  minl <- min(out$lambda)
  minl <- signif(minl,digits = 4)
  if (noseg < 200) {
    l <- 0
  } else if(noseg <500){
    l <- max(0.3,minl)
  } else if(noseg <1000){
    l <- max(0.5,minl)
  } else {
    l <- max(1, minl)
  }
  if (l == 0) {
    newnoseg <- noseg
    } else{
      newnoseg <- getNoSeg(out,l)
    }
  cat(paste(c(seriesName,arrayName,l,noseg,newnoseg,paste(summary(out$lambda),collapse = ','),paste(summary(segStep$value[segStep$value>0]),collapse = ','),paste(summary(segStep$value[segStep$value<0]),collapse = ',')), collapse ='\t'),'\n',sep = '',append = T, file = 'output/caliLambda.tsv')
}


args <- commandArgs(trailingOnly = T)
workingdir <- args[1]
seriesName <- args[2]
arrayName <- args[3]
gp <- args[4]

CalibrateLambda(workingdir,seriesName,arrayName,gp)
