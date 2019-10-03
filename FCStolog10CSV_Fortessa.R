rm(list=ls())
library(flowCore)
library(flowStats)
library(flowUtils)
library(flowViz)
library(flowWorkspace)

#install.packages("tools") #currently installed; need it for file_path_sans_ext
library (tools)

parent.dir <- ""
setwd(parent.dir)

filenames <- list.files(".", pattern=".fcs", recursive = TRUE, include.dirs = TRUE)
filename <- NULL;

i <- 0;

for (file in filenames)    
{
  data <- read.FCS(file, transformation = FALSE)  
  dataframe <- exprs(data)
  dataframe[dataframe==0] <- NA
  dataframe <- dataframe[complete.cases(dataframe[,c("FSC-A", "V_450/50-A", "B_530/30-A")]),] #FORTESSA
  logdataframe <- log10(dataframe)
  
  i<-i+1;
  filename[i] <- filenames[i]
  name <- paste(file_path_sans_ext(basename(filename[i])),".csv") #get filename without path
  name <- gsub(" ","", name) #delete spaces
  write.csv(logdataframe, file=name, na="", row.names=FALSE)
}
