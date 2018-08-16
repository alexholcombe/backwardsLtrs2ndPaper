
#Import raw data
dataPath<- file.path("E2/Data")
#Experiment was administered by MATLAB
#.mat file been preprocessed into melted long dataframe by importE1data.Rmd
file.exists( file.path(dataPath, "backwards2E2_rawDataFromMAT.rda") )
getwd()
data<- readRDS( file.path(dataPath, "backwards2E2_rawDataFromMAT.rda") ) 

#tidy data
library(dplyr)
df<- data
#to work with dplyr, can't have array field like letterSeq
df$letterSeq<- NULL

numItemsInStream<- length( data$letterSeq[1,] )  
minSPE<- -17; maxSPE<- 17

dg<- df %>% filter(subject < "AC")


library(mixRSVP)
analyzeOneCondition(df %>% filter(subject=="AA"),numItemsInStream,parameterBounds(),nReplicates=5)
