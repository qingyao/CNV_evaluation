library('robfilter')
library('DNAcopy')
library('kernlab')
library('e1071')

workPath <- file.path(getwd(),'CNARAsrc')
options("scipen"=100,"digits"=4)
source(file.path(workPath,'CNARA.R'))

args <- commandArgs(trailingOnly = T)
series <- args[1]
samplePath <- args[2] #'/Volumes/arraymapMirror/arraymap/hg19/'
update <- args[3]
try(if (!dir.exists(samplePath)) stop("check server directory"))
probeFileName <- "probes,cn.tsv"
segFileName <- "segments,cn.tsv"

trainingFile <- file.path(workPath, "trainingSet.txt")

arrayDir <- file.path(samplePath, series)
arrayList <- list.dirs(arrayDir)[-1]
#arrayList <- unlist(regmatches(arrayList,gregexpr('GSM[0-9]+',arrayList,perl=T)))
arrayList <- sapply(arrayList,function(x) {
    allfolders <- strsplit(x,"/",fixed=T)[[1]]
    return(allfolders[length(allfolders)])
})
names(arrayList) = NULL
noArray <- length(arrayList)

summaryFile <- file.path(getwd(),'CNARA_summary.tsv')
outputFile <- file.path(arrayDir, "CNARA_output.txt")
tableHeader <- paste("Sample", "Series", "Speak", "Breakpoint_Step", "Breakpoint_CBS", "Spread", "Classification_Label", "Decision_Value", "Good_Poor", "Case_Diagnosis", sep = "\t")
if (!file.exists(outputFile) | update == 1){
    write.table(tableHeader, file = outputFile, append = F, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
if (!file.exists(summaryFile) | update == 1){
    write.table(tableHeader, file = summaryFile, append = F, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
classifier <- trainSVM(trainingFile)

for (arr in 1:noArray) {

  if (file.exists(file.path(arrayDir,arrayList[arr],'sample_evaluation.tsv')) & update != 1) next
  print(paste0("processing ", arrayList[arr], " (", arr, " of ", noArray, ")"))
  probeFile <- file.path(arrayDir, arrayList[arr], probeFileName)
  segFile <- file.path(arrayDir, arrayList[arr], segFileName)

  newCNProbe <- readProbe(probeFile=probeFile, sampleID=arrayList[arr])
  newSpeakCNAno <- calSpeakCNAno(newCNProbe)

  #plot S graph
  plot(quality(newSpeakCNAno), xlab = "number of iterations", ylab = "S", main=arrayList[arr])

  segNumberCBS <- calCBSBreakpoints(newCNProbe)
  segSpread <- calSpread(newCNProbe, segFile=segFile)
  segSpread <- round(segSpread,4)
  #segSpread <- calSpread(newCNProbe) # if you don't have segmentation file

  CNAno <- numberOfCNA(newSpeakCNAno)
  Speak <- speak(newSpeakCNAno)
  Speak <- round(Speak,4)
  #create a new object "Metrics" for the copy number profile for quality assessment
  CNProfileMetrics <- createMetrics(sampleID=arrayList[arr], speak=Speak, numberOfCNA=CNAno, cbsBreakpoints=segNumberCBS, spread=segSpread)

  #quality assessment for the copy number profile
  assessment <- assessQuality(CNProfileMetrics, svmClassifier=classifier)

  tmp <- paste(arrayList[arr], series, Speak, CNAno, segNumberCBS, segSpread, assessment$label, round(assessment$decision.values,4), assessment$flag, assessment$caseDiag, sep = "\t")
  write.table(tmp, file = outputFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(tmp, file = summaryFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  ### write QA log file for each sample
  sampleLog <- file.path(arrayDir, arrayList[arr], "sample_evaluation.tsv")
  for (entry in 1:length(strsplit(tableHeader, "\t")[[1]])){
    write.table(paste(strsplit(tableHeader, "\t")[[1]][entry], strsplit(tmp, "\t")[[1]][entry], sep="\t"), file = sampleLog, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

    #plot a best fit and counter-fit for a copy number profile. Call calSpeakCNAno again, can be omitted if you just want the assessment without the plot.
    # plotFCF(calSpeakCNAno(newCNProbe, iterations=CNAno))
  }
}
