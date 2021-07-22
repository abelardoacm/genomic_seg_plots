suppressMessages(library("tidyverse"))

PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
#PosArgs <- "Alpha"
dataset_name = PosArgs[1]
infile <- paste("../results/LCR_counts/",dataset_name,"_LCR_matrix.csv",sep="")

LCR_count <- read.table(infile, sep = ',', header = TRUE, na.strings = "", fill=TRUE,quote = "")

LCR_count$date <- gsub("^.*202","202",LCR_count$virus)
LCR_count$year <- gsub("_.*","",LCR_count$date)
LCR_count$month <- gsub("0","",gsub("_.*","",sub(".*?_","",LCR_count$date)))
LCR_count$day <-  as.numeric(gsub("_","",sub("^.*?_","",sub("^.*?_","",LCR_count$date))))
LCR_count[is.na(LCR_count)] <- 0

LCR_count <- LCR_count[order( LCR_count$year, LCR_count$month, LCR_count$day),]

plot(LCR_count$delta_kappa)
