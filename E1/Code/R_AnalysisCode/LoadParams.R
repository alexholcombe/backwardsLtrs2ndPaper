# Backwards Letters - Study 2 - Experiment 1
# This script loads the parameters (efficacy, latency, and precision) 
# from the matlab data file--output from mixture modelling code. 
# Creates a seperate data frame for each parameter in both wide and long format. 

library(R.matlab)

E1.data <- readMat('/Users/admin/Documents/MATLAB/BackwardsLetters/Study_2/Exp1_Random_0and180degs/Data/Exp1_InvertedLtrs_ParamsForAnalysis.mat')
Subject <- c('AA','AB', 'AC', 'AD', 'AH',  'AI', 'AJ', 'AK', 'AL', 'AO', 'AU', 'AV', 'AW', 'AX', 'AZ', 'BA', 'BC', 'BD', 'BE', 'BF', 'BG', 'BI', 'BJ', 'BK', 'BN', 'BO', 'BP', 'BQ', 'BR', 'BS')


# Create dataFrame for efficacy
LeftCanon <- E1.data$E1.allEstMinusExcluded[1,1,,1]*100 # [Condition, Stream, Participant, Parameter]
RightCanon <- E1.data$E1.allEstMinusExcluded[1,2,,1]*100 # Muliply by 100 to convert to percentage
LeftInvert <- E1.data$E1.allEstMinusExcluded[2,2,,1]*100 # Left and right stream data are flipped in the
RightInvert <- E1.data$E1.allEstMinusExcluded[2,1,,1]*100 # inverted condition so we unflip them here
efficacy.wide <- data.frame(Subject, LeftCanon, RightCanon, LeftInvert, RightInvert)
remove(LeftCanon, RightCanon, LeftInvert, RightInvert)


# Create dataFrame for latency
LeftCanon <- E1.data$E1.allEstMinusExcluded[1,1,,2]*(1000/10) # [Condition, Stream, Participant, Parameter]
RightCanon <- E1.data$E1.allEstMinusExcluded[1,2,,2]*(1000/10) # Mutiply by 1000msec divided by item rate (10)
LeftInvert <- E1.data$E1.allEstMinusExcluded[2,2,,2]*(1000/10) # Left and right stream data are flipped in the
RightInvert <- E1.data$E1.allEstMinusExcluded[2,1,,2]*(1000/10) # inverted condition so we unflip them here
latency.wide <- data.frame(Subject, LeftCanon, RightCanon, LeftInvert, RightInvert)
remove(LeftCanon, RightCanon, LeftInvert, RightInvert)


# Create dataFrame for precision
LeftCanon <- E1.data$E1.allEstMinusExcluded[1,1,,3]*(1000/10) # [Condition, Stream, Participant, Parameter]
RightCanon <- E1.data$E1.allEstMinusExcluded[1,2,,3]*(1000/10) # Mutiply by 1000msec divided by item rate  (10)
LeftInvert <- E1.data$E1.allEstMinusExcluded[2,2,,3]*(1000/10) # Left and right stream data are flipped in the
RightInvert <- E1.data$E1.allEstMinusExcluded[2,1,,3]*(1000/10) # inverted condition so we unflip them here
precision.wide <- data.frame(Subject, LeftCanon, RightCanon, LeftInvert, RightInvert)
remove(LeftCanon, RightCanon, LeftInvert, RightInvert)


# Convert efficacy to long format
library(reshape2)
efficacy.long <- melt(data=efficacy.wide, id.var="Subject",
                  measure.vars=c("LeftCanon", "RightCanon", "LeftInvert", "RightInvert"),
                  variable.name="Condition")
names(efficacy.long)[names(efficacy.long)=="value"] <- "Efficacy"

# Split Condition column into Stream and Orientation
efficacy.long$Stream <- NA
efficacy.long$Stream[grepl("^Left",  efficacy.long$Condition)] <- "Left"
efficacy.long$Stream[grepl("^Right", efficacy.long$Condition)] <- "Right"
efficacy.long$Stream <- factor(efficacy.long$Stream)

efficacy.long$Orientation <- NA
efficacy.long$Orientation[grepl("Canon$",  efficacy.long$Condition)] <- "Canonical"
efficacy.long$Orientation[grepl("Invert$", efficacy.long$Condition)] <- "Inverted"
efficacy.long$Orientation <- factor(efficacy.long$Orientation, levels=c("Canonical","Inverted"))

# Remove the Condition column now
efficacy.long$Condition <- NULL

# Convert latency to long format
library(reshape2)
latency.long <- melt(data=latency.wide, id.var="Subject",
                    measure.vars=c("LeftCanon", "RightCanon", "LeftInvert", "RightInvert"),
                    variable.name="Condition")
names(latency.long)[names(latency.long)=="value"] <- "Latency"

# Split Condition column into Stream and Orientation
latency.long$Stream <- NA
latency.long$Stream[grepl("^Left",  latency.long$Condition)] <- "Left"
latency.long$Stream[grepl("^Right", latency.long$Condition)] <- "Right"
latency.long$Stream <- factor(latency.long$Stream)

latency.long$Orientation <- NA
latency.long$Orientation[grepl("Canon$",  latency.long$Condition)] <- "Canonical"
latency.long$Orientation[grepl("Invert$", latency.long$Condition)] <- "Inverted"
latency.long$Orientation <- factor(latency.long$Orientation, levels=c("Canonical","Inverted"))

# Remove the Condition column now
latency.long$Condition <- NULL


# Convert precision to long format
library(reshape2)
precision.long <- melt(data=precision.wide, id.var="Subject",
                      measure.vars=c("LeftCanon", "RightCanon", "LeftInvert", "RightInvert"),
                      variable.name="Condition")
names(precision.long)[names(precision.long)=="value"] <- "Precision"

# Split Condition column into Stream and Orientation
precision.long$Stream <- NA
precision.long$Stream[grepl("^Left",  precision.long$Condition)] <- "Left"
precision.long$Stream[grepl("^Right", precision.long$Condition)] <- "Right"
precision.long$Stream <- factor(precision.long$Stream)

precision.long$Orientation <- NA
precision.long$Orientation[grepl("Canon$",  precision.long$Condition)] <- "Canonical"
precision.long$Orientation[grepl("Invert$", precision.long$Condition)] <- "Inverted"
precision.long$Orientation <- factor(precision.long$Orientation, levels=c("Canonical","Inverted"))

# Remove the Condition column now
precision.long$Condition <- NULL





