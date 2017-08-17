# Backwards Letters - Study 2 - Experiment 1

setwd("/Users/admin/Desktop/WorkInProgress/BackwardsLetterStudy-2/E1/code/Rcode/")
source("StatsFunctions.R")

# Make efficacy graph
efficacySummary <- summarySEwithin(efficacy.long, measurevar="Efficacy", withinvars=c("Stream","Orientation"), idvar="Subject")

library(ggplot2)
ggplot(efficacySummary, aes(x=Orientation, y=Efficacy, fill=Stream)) +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=Efficacy-ci, ymax=Efficacy+ci)) +
  coord_cartesian(ylim=c(0,100)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  theme_bw() + 
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  theme(legend.key = element_rect(colour = "black"))


# Make latency graph
latencySummary <- summarySEwithin(latency.long, measurevar="Latency", withinvars=c("Stream","Orientation"), idvar="Subject")

library(ggplot2)
ggplot(latencySummary, aes(x=Orientation, y=Latency, fill=Stream)) +
  ylab("Latency (ms)") +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=Latency-ci, ymax=Latency+ci)) +
  coord_cartesian(ylim=c(0,150)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_y_continuous(breaks = seq(0, 150, by = 25)) +
  theme_bw() + 
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  theme(legend.key = element_rect(colour = "black"))


# Make precsion graph
precisionSummary <- summarySEwithin(precision.long, measurevar="Precision", withinvars=c("Stream","Orientation"), idvar="Subject")

library(ggplot2)
ggplot(precisionSummary, aes(x=Orientation, y=Precision, fill=Stream)) +
  ylab("Precision (ms)") +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=Precision-ci, ymax=Precision+ci)) +
  coord_cartesian(ylim=c(0,150)) +
  scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  scale_y_continuous(breaks = seq(0, 150, by = 25)) +
  theme_bw() + 
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  theme(legend.key = element_rect(colour = "black"))

