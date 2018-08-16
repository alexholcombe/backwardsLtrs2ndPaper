#Loads raw data from MATLAB files for second backwards paper first experiment
require(R.matlab)

#define function to be used below
turnMATintoMeltedDataframe<- function(fromMAT) {
  #It was imported as a list of lists
  mydata<-fromMAT
  numTrials<- length( unlist( mydata$allConditions ) )
  #Some of the list members are one-off fields pertaining to entire experiment. Remove them, put in miscInfo
  miscInfo <- list()
  membersToDelete<- list()
  #loop through MATLAB .mat file data structure
  for (i in 1:length(mydata)) {
    thisListMember <- mydata[i]
    name<- names(thisListMember)
    len<- length(unlist(thisListMember))
    #cat(name, len,'\n')           
    if (len < numTrials) { #one-off information, not whole experiment
      miscInfo[name] <- thisListMember
      membersToDelete<-c(membersToDelete,name)
    } 
  }
  for (m in membersToDelete) {
    mydata[m]<-NULL #remove from list to leave only every-trial fields
  }
  
  #Remaining fields are matrices, one dimension of which is the trial number. Other dimension is whichTarget in dual-stream case
  #To melt, need to unpack the matrix, going through the second dimension and taking it out, putting at end to create flat dataframe
  # with one row for each trial.
  meltedDf <- data.frame()
  for (target in 1:2) {
    dfThis <- data.frame(trial=seq(1,numTrials))
    dfThis$target<- target
    for ( i in 1:length(mydata) ) {  #Go through each field, turn into column
      thisMember <- mydata[i]
      name <- names( thisMember )
      thisMember <- thisMember[[name]] #for some reason it's always a list of 1 (with the name name)
      #if only one dimension then it is a condition name applying to both streams
      thisMember <- drop(thisMember) #remove any singleton dimensions
      ndim <- length( dim( thisMember ) ) #number of dimensions
      #if (name=="allLetterOrientations" | name=="allTargets") {
      #  print(paste0("name=",name," thisMember=",thisMember," ndim=",ndim))
      #}
      if ((ndim==0)) {  #it is a condition name applying to both streams
        dfThis[[name]] <- thisMember
      } else if (ndim==2) { #else it has separate values for each target
        if (dim(thisMember)[2] > 2) { #allLetterOrientations has 2 dimensions but 16 different values for second dimension instead of 2
          dfThis[[name]] <- thisMember #entire array
        } else {
          dfThis[[name]] <- thisMember[,target]
        }
      } else if (ndim==3) { #e.g. allLetterOrder has 3 dimensions
        dfThis[[name]] <-  thisMember[,target,]
      } else { print("unexpected dims")}
      #The above should work for any variables, don't need to know their name or specific formats
      
      #The below is based on kmnowledge of the variables in this case
      dfThis$subject<- miscInfo$participantID[[1]]
      
      #calculate the serial position of the response from the allLetterOrder
      thisStreamLetterOrder<- mydata$allLetterOrder[,target,]
      allResponses<- mydata$allResponses[,target]
      
      #first dimension is trial. second dimension is serial position
      #for (trial in 1:numTrials) {
      #  #for each trial in allResponses, do the match
      #  respSP<- match(allResponses[trial], thisStreamLetterOrder[trial,])
      #}
      #mapply doesn't work because thisStreamLetterOrder is not a list of things I want to operate on
      #mapply(match, 1, list(c(3,2,1))) #this works. So, need to turn thisStreamLetterOrder into many lists
      eachTrialLetterOrderList<- split(thisStreamLetterOrder, row(thisStreamLetterOrder))
      responsePositions<- mapply(match, allResponses, eachTrialLetterOrderList )
      dfThis$respSP<- responsePositions
    }
    meltedDf<-rbind(meltedDf,dfThis)
  }
  
  E<-meltedDf
  library(purrr) #for map_if
  #data.frames and tibbles must have all atomic or list entries, not arrays, therefore
  #turn any remaining arrays (like letterSequence) into list so tidyverse operations can be done
  E<- map_if(E,is.matrix,~split(.,seq(nrow(.)))) %>% as_tibble #https://stackoverflow.com/questions/51622017/change-all-array-columns-of-a-data-frame-into-lists
  
  names(E)[names(E) == 'allTargets'] <- 'targetSP'
  # allLetterOrder         : num [1:216, 1:2, 1:16] 13 6 9 1 10 16 11 1 12 8 ...
  names(E)[names(E) == 'allLetterOrder'] <- 'letterSeq'
  #will have to set letterSeq to NULL later to do dataframe operations because not allowed to have a list
  
  #The only other variable that is a list even after splitting by target is allLetterOrientations
  # allLetterOrientations  : num [1:216, 1:16] 1 0 1 0 0 0 1 0 1 0 ...
  # Because doesn't have second dimension, unlike allLetterOrder, seems only the target stream was saved?
  # We will have to set it to NULL to avoid list, but want to preserve the orientations at
  # least for a few items around the target so can assess if takes time to build up.
  # Need to ask Chris to verify it's the target stream it's recording.
  names(E)[names(E) == 'allLetterOrientations'] <- 'targetStreamOrientations'
  names(E)[names(E) == 'allResponses'] <- 'resp'
  names(E)[names(E) == 'allConditions'] <- 'condition'
  names(E)[names(E) == 'allTargetOrientations'] <- 'targetOrient'
  names(E)[names(E) == 'allResponseOrientations'] <- 'respOrient'
  names(E)[names(E) == 'allResponseOrder'] <- 'respOrder'
  return(E)
}
