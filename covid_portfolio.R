#-----------------housekeeping--------------------------------------------------------------------------------------

##set working directory##
setwd()

##load libraries##

library(seqinr)
library(dplyr)
library(tidyr)
library(ampir)
library("Biostrings") #this is a specialized package containing functions for analysing and performing sequence alignments.

##Read file contain covid sequences and IDs##
assay1 <- read.fasta(file="1_msa.fasta", set.attributes = FALSE, as.string=TRUE)
assay2 <- read.fasta(file="2_msa.fa", set.attributes = FALSE, as.string=TRUE)
assay3 <- read.fasta(file="3_msa.fa", set.attributes = FALSE, as.string=TRUE)

#----------Table1-------------------------------------------------------------------------------------------------------------
##Produces a table that contains the strain code and the associated sequences.

#table.1 is the original table where column names is the strain code and the rows are the sequence.
table.1 <- data.frame(assay2) #create data frame from assay data.
table.1 <- tibble::rownames_to_column(table.1, "strain_code")#add column with column name "strain code".
table.1$strain_code="sequence"#add "sequence" value so it can become the column name.
table.strain_code_and_sequence <- gather(table.1, key = "strain_code",value = "sequences")#reshape the data so that the column names become the row names witht he column name "strain_code" and the dna sequences become a column with the column name "sequences".

#save dataframe as tab delimited file
temp <- write.table(table.strain_code_and_sequence,"table.strain_code_and_sequence_assay2.txt", sep = "\t", row.names = FALSE)

#---------Table2------------------------------------------------------------------------------------------------
##Produces a table with the sequence and the frequency of each sequence##

table.sequences.assay1 <- data.frame(matrix(unlist(assay1),nrow=length(assay1),byrow=TRUE))#create a data frame containing just the sequences
colnames(table.sequences.assay1) <- c("sequences")#change column title to sequences

table.sequences.assay2 <- data.frame(matrix(unlist(assay2),nrow=length(assay2),byrow=TRUE))#create a data frame containing just the sequences
colnames(table.sequences.assay2) <- c("sequences")#change column title to sequences


table.sequences.assay3 <- data.frame(matrix(unlist(assay3),nrow=length(assay3),byrow=TRUE))#create a data frame containing just the sequences
colnames(table.sequences.assay3) <- c("sequences")#change column title to sequences

table.NB_of_duplicates <- table(as.data.frame(table.sequences))#generate frequency of each sequence
table.NB_of_duplicates <- data.frame(table.NB_of_duplicates)#change list to data frame
colnames(table.NB_of_duplicates) <- c("sequences","Number_of_identical_sequences")

#save dataframe as excelsheet
temp2 <- write.table(table.NB_of_duplicates,"table.NB_of_duplicatess_assay2.txt", sep = "\t", row.names = FALSE)
temp5 <- write.table(table.sequences,"table.sequences_assay3.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

#-----------Table3--------------------------------------------------------------------------------------------------------------------
##produces table with the sequence and the strain code with this particular sequence.

table.unique <- unique.data.frame(table.sequences)#create dataframe containing only the unique sequences
list.unique <- as.vector(table.unique$sequences)#convert the data frame into a list containing all 925 unique sequences

#create a list of lists where for each sequence and list of the indices where table.strain_code_and_sequences$sequences is equal to sequences in list.unique.
Index <- lapply(list.unique, function(x) which(table.strain_code_and_sequence$sequences == x))
table.sequences_and_strain_code <- mutate(table.unique, strains_with_identical_sequences=lapply(Index, function(x) table.strain_code_and_sequence$strain_code[x]))
table.sequences_and_strain_code$strains_with_identical_sequences <-vapply(table.sequences_and_strain_code$strains_with_identical_sequences, paste, collapse = ",", character(1L))

temp3 <- write.table(table.sequences_and_strain_code,"table.sequences_and_strain_code_assay2.txt", sep = "\t", row.names = FALSE)

#----testing---------
test2$strains_with_identical_sequences <-vapply(test2$strains_with_identical_sequences, paste, collapse = ", ", character(1L))
test2 <- table.sequences_and_strain_code[10:11,]
test2 <- as.character(test2)
class(test2$strains_with_identical_sequences)
test1 <- table.sequences_and_strain_code$strains_with_identical_sequences[1:1000]

#test1<-as.data.frame(table.strain_code_and_sequences[1:10000,])#select a samll sample of the data to speed up testing of code.
#colnames(test1)<-c("sequences")#change column name

#head(test2)
#--------Primer extraction--------------------------------------------------------------------------

#primers assay 1
HARVARD_As1e_F3 <- "cggtggacaaattgtcac"
HARVARD_As1e_B3 <- "aagtgtgttaaatccagagaag"
HARVARD_As1e_FIPpart_1 <- "taaatttttggctttgtgtgctga"
HARVARD_As1e_FIPpart_2 <- "ctgtgcaaaggaaattaaggag"
HARVARD_As1e_BIPpart_1 <- "tattggtggagctaaacttaaagcctt"
HARVARD_As1e_BIPpart_2 <- "cactcaaagggattgtacagaa"
HARVARD_As1e_LF <- "agtgttcagacattctttaagcttgtaa"
HARVARD_As1e_LB <- "ttgaatttaggtgaaacatttgtcacg"

#primer assay 2
COLOR_N_F3 <- "aacacaagctttcggcag"
Color_N_B3 <- "ggatgacaaagatccaaatttc"
Color_N_FIPpart1 <- "ctgattacaaacattggccgca"
Color_N_FIPpart2 <- "ccaaggaaattttggggac"
Color_N_BIPpart1 <- "cgcattggcatggaagtcac"
Color_N_BIPpart2 <- "ctacacaggtgccatcaaa"
COLOR_N_LF <- "gaactaatcagacaaggaa"
COLOR_N_LB <- "accttcgggaacgtggtt"

#primer assay 3

NEB_Gene_E1_F3 <- "tgagtacgaacttatgtactcat"
NEB_Gene_E1_B3 <- "actctcgtgttaaaaatctgaa"
NEB_Gene_E1_FIPpart1 <- "acttctttttcttgctttcgtggt"
NEB_Gene_E1_FIPpart2 <- "tcgtttcggaagagacag"
NEB_Gene_E1_BIPpart1 <- "ttgctagttacactagccatcctta"
NEB_Gene_E1_BIPpart2 <- "acgtgagtcttgtaaaacct"
NEB_Gene_E1_LF <- "cgttaatagttaatagcg"
NEB_Gene_E1_LB <- "gcgcttcgattgtgtgcgt"



test1 <- slice_head(table.sequences,n=100)


#find primer indicies
primer_indicies(seq1,HARVARD_As1e_B3)

#code from primer indicies#
seq1<-"cggtggacaaattgtcacctgtgcaaaggaaattaaggagagtgttcagacattctttaagcttgtaaataaatttttggctttgtgtgctgactctatcattattggtggagctaaacttaaagccttgaatttaggtgaaacatttgtcacgcactcaaagggattgtacagaaagtgtgttaaatccagagaag"

start <- regexpr(HARVARD_As1e_FIPpart_2,seq1)#function finds the exact match of primer in the dna sequence
end <- start + attr(start,"match.length") - 1#need to subtract as range include the first base, would include one extra base otherwise.

#testing/quality assurance
nchar(NEB_Gene_E1_LB)
start
end

#funciton to subset all sequences at where the primer matches.
#convert all sequences to DNAstrings
#table.sequences.asDNAString<-apply(table.sequences, 1, function(x) DNAString(x))

subset_sequences <- lapply(table.sequences.assay1, function(x) subseq(x, start, end))#subset DNA sequence according to where primer matches (numbers got from local alignment pairwiseAlignmnet())

#save data
subset_sequences <- lapply(subset_sequences, function(x) as.character(x))#dna are in format DNAString() so must be converted before saving.
subset_sequences <- data.frame(matrix(unlist(subset_sequences),nrow=length(subset_sequences),byrow=TRUE))#turn list into table for easy saving.
subset_sequences <- data.frame(matrix(unlist(subset_sequences),nrow=length(subset_sequences),byrow=TRUE))#turn list into table for easy saving.


temp <- write.table(subset_sequences, file = "subset_sequences_HARVARD_As1e_FIPpart_2.txt", sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)


length(subset_sequences$matrix.unlist.subset_sequences...nrow...length.subset_sequences...)


testing <- subset_sequences$matrix.unlist.subset_sequences...nrow...length.subset_sequences...[subset_sequences$matrix.unlist.subset_sequences...nrow...length.subset_sequences...=="CGGTGGACAAATTGTCACC"]
