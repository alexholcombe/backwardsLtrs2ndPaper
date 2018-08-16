library(R.matlab)
library(ggplot2)
library(reshape2)
rm(list=ls())

setwd('~/Google Drive/Backwards paper/secondPaper/E2/Data/')

datafiles <- list.files(pattern = '^[A-Z][A-Z]')

chars <- 'ABCDEFGJKLPQRTUVY'
chars <- unlist(strsplit(chars,split=''))

IDs <- c()

nRows <- 200* length(datafiles)*2

masterDF <- expand.grid(ID = character(nRows), condition = factor(x=character(1), levels = c(1,2), labels = c('Normal', 'Rotated')), responseOrientation = numeric(1),trial = numeric(1), stream = numeric(1), responsePosRelative = numeric(1), stringsAsFactors = F)

start <- 1
for(file in datafiles){
  ID <- paste(unlist(strsplit(file, split =''))[1:2], collapse = '')
  IDs <- c(IDs, ID)
  
  
  matFile <- readMat(file) 
  
  
  responsePosRelativeMatrix <- matrix(NA, nrow = 200, ncol = 5)
  responsePosMatrix <- matrix(NA, nrow = 200, ncol = 5)
  
  for(stream in 1:2){
    for(trial in 1:200){
      thisResponse <- matFile$allResponses[trial, stream]
      thisStream <- matFile$allLetterOrder[trial, stream,]
      responsePos <- which(thisStream == thisResponse)
      responsePosRelative <- responsePos - matFile$allTargets[trial, stream]
      responsePosRelativeMatrix[trial,stream] <- responsePosRelative
      responsePosMatrix[trial,stream] <- responsePos
    }
  }
  responsePosRelativeMatrix[,3] <- matFile$allConditions[1,] 
  responsePosRelativeMatrix[,4] <- matFile$allResponseOrientations[,1]
  responsePosRelativeMatrix[,5] <- 1:200
  #SPE left, SPE right, Target orientation, response orientation, trial
  
  responsePosMatrix[,3] <- matFile$allConditions[1,]
  responsePosMatrix[,4] <- 1:200
  
  responsePosRelativeDF <- data.frame(responsePosRelativeMatrix)
  responsePosRelativeDF <- cbind(rep(ID, times=200), responsePosRelativeDF)
  
  #responsePosRelativeDF cols, ID, SPE left, SPE right, Target orientation, response orientation, trial
  
  colnames(responsePosRelativeDF) <- c('ID','responsePosRelative0', 'responsePosRelative1', 'condition', 'responseOrientation', 'trial')
  responsePosRelativeDF <- melt(responsePosRelativeDF,measure.vars = c('responsePosRelative0','responsePosRelative1'), stringsAsFactors=F)
  responsePosRelativeDF$ID <- as.character(responsePosRelativeDF$ID) #Default behaviour for melt is stringsAsFactors=False
  responsePosRelativeDF$condition <- factor(responsePosRelativeDF$condition, levels = c(1,2), labels  = c('Normal', 'Rotated'))
  colnames(responsePosRelativeDF)[5:6] <- c('stream', 'responsePosRelative')
  end <- start + nrow(responsePosRelativeDF) -1
  print(start)
  print(end)
  masterDF[start:end,] <- responsePosRelativeDF
  start <- end+1
}


masterDF$correctOrientation <- ifelse(as.numeric(masterDF$condition) == masterDF$responseOrientation+1, TRUE, FALSE)

#masterDF$ID <- rep(c('CB','CL','PP'), each = 400)

ggplot(masterDF[masterDF$correctOrientation,], aes(x=responsePosRelative))+
  geom_histogram(binwidth = 1)+
  facet_wrap(~ID+condition+stream,labeller = 'label_both')+
  labs(title = 'Histograms for responses with correct orientation')


ggplot(masterDF[!masterDF$correctOrientation,], aes(x=responsePosRelative))+
  geom_histogram(binwidth = 1)+
  facet_wrap(~ID+condition+stream,labeller = 'label_both')+
  labs(title = 'Histograms for responses with incorrect orientation')


ggplot(masterDF, aes(x=responsePosRelative))+
  geom_histogram(binwidth = 1, aes(fill = factor(correctOrientation, levels = c('TRUE', 'FALSE'))), alpha= .5, position = 'identity')+
  facet_wrap(~ID+condition+stream,labeller = 'label_both')+
  labs(fill = 'correctOrientation')+
  scale_fill_manual(values = c("TRUE" = '#0D7526','FALSE'= '#F8AD3F'))