#install and load libraries/packages.
library(seqinr)

##set working directory##
getwd()
setwd()

##Read files containing covid sequences and IDs##
sequences <- read.fasta(file="1_msa.fasta", set.attributes = FALSE, as.string=TRUE)

duplicate_counter(sequences,"table.NB_of_duplicates_funcitontest.txt")

#---------function------------------------------------------------------------------------------------------------
##Produces a tab delimited file with the sequence and the frequency of each sequence##
#x is a list containing DNA sequences and their strain codes/IDs
#path is the file path to save the data.
duplicate_counter <- function(x,path){
  
table.sequences <- data.frame(matrix(unlist(x),nrow=length(x),byrow=TRUE))#create a data frame containing just the sequences
colnames(table.sequences) <- c("sequences")#change column title to sequences
table.NB_of_duplicates <- table(as.data.frame(table.sequences))#generate frequency of each sequence
table.NB_of_duplicates <- data.frame(table.NB_of_duplicates)#change list to data frame
colnames(table.NB_of_duplicates) <- c("sequences","number_of_identical_sequences")

#save dataframe as excelsheet
temp <- write.table(table.NB_of_duplicates, path, sep = "\t", row.names = FALSE)
}
