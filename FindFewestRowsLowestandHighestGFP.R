#This code does:
#1) Find the minimum number of rows in all the files
#2) Finds the minimum and maximum GFP values

rm(list=ls()) #clears the environment/workspace (especially important when using R studio)
#install.packages("tools") #currently installed; need it for file_path_sans_ext
library (tools)

parent.dir <- "" #insert pathname here
setwd(parent.dir)

filenames <- list.files(".", pattern=".csv", recursive = TRUE, include.dirs = TRUE) #get all filenames ending in csv
filename <- NULL;
fewestrows=10000000;
fewestrowsfilename <- NULL;
lowestGFP <- 100;
highestGFP <- 0;
i <- 0;

for (file in filenames)    
{

  i<-i+1;
  filename[i] <- filenames[i]
  basename <- file_path_sans_ext(basename(filename[i])) #get filename without path
  
  data <- read.csv(file, header=T)  #read the csv file
  
  cat(paste("\nCurrent lowest # rows: ", fewestrows))
  cat(paste("\nCurrent lowestGFP: ", lowestGFP))
  cat(paste("\nCurrent highestGFP: ", highestGFP))
  
  #if current number of rows are less than fewestrows, change fewestrows to current number of rows and store filename
  if (nrow(data)<fewestrows)
  {
    fewestrows = nrow(data)
    fewestrowsfilename = basename
  }
  
  #if current minimum GFP value is less than lowestGFP, change lowestGFP to current minimum GFP value
  currentmingfp <- min(data$'B_530_30.A')
  if(currentmingfp < lowestGFP)
  {
    lowestGFP <- min(data$'B_530_30.A')
  }
  
  #if current maxium GFP value is highest than highestGFP, change highestGFP to current highest GFP value
  currentmaxgfp <- max(data$'B_530_30.A')
  if(currentmaxgfp > highestGFP)
  {
    highestGFP <- max(data$'B_530_30.A')
  }
}

cat (paste("\n~Fewest Rows: ",fewestrows, "\n~Name: ", fewestrowsfilename))
cat (paste("\n~Lowest GFP: ",lowestGFP))
cat (paste("\n~Higest GFP: ",highestGFP))