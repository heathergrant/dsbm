
###################################################################################################################################################
#
# n525 clean up scueal output 
###################################################################################################################################################

library(seqinr)
library(tidyverse)
library(ape)
library(dplyr)
require(ggplot2)

cleanUpScuealRawFile <- function(scueal.output.fn) {
  SCUEAL <- read.table(scueal.output.fn, header=T, as.is=T, row.names=1, sep='\t')
  SCUEAL <- SCUEAL[, c(1,2,3,8)]
  SCUEAL$Name <- trimws(as.character(SCUEAL$Name), which="both")
  SCUEAL$Simplified.Subtype <- gsub(" intra-subtype recombinant", "", SCUEAL$Simplified.Subtype, fixed=T)
  SCUEAL$Simplified.Subtype <- gsub(" inter-subtype recombinant", "", SCUEAL$Simplified.Subtype, fixed=T)
  SCUEAL$Simplified.Subtype <- trimws(SCUEAL$Simplified.Subtype, which="both")
  SCUEAL$Simplified.Subtype <- gsub(" recombinant","", SCUEAL$Simplified.Subtype, fixed=T)
  
  SCUEAL$Subtype <- gsub(" intra-subtype recombinant", "", SCUEAL$Subtype, fixed=T)
  SCUEAL$Subtype <- gsub(" inter-subtype recombinant", "", SCUEAL$Subtype, fixed=T)
  SCUEAL$Subtype <- trimws(SCUEAL$Subtype, which="both") #trimws = trim white space - leading or trailing (in our case, 'both')
  SCUEAL$Subtype <- gsub(" recombinant","", SCUEAL$Subtype, fixed=T)
  
  SCUEAL$Breakpoints <- gsub("[(][0-9]+[-][0-9]+[)]", "", SCUEAL$Breakpoints)

  return(SCUEAL)
}

n525<-cleanUpScuealRawFile("../DYSBM/n525/scueal_genome/525.fa.out")
n525$PDF<-rownames(n525)
head(n525)
###################################################################################################################################################
#
# Add length of sequence 
seqs <- read.fasta("../DYSBM/n525/scueal_genome/525.fa")
subty<-n525
numBP <- c()
for (i in 1:nrow(subty)){
  nam <- subty$Name[i]
  seq1 <- seqs[[which(attr(seqs,"name")==nam)]]
  numBP <- c(numBP, length(which(seq1!="-")))
  numBP<-as.numeric(numBP)
}
n525<-cbind(n525, numBP) 

###################################################################################################################################################
#compare to scueal.csv 
scueal_orig<-read.csv("../DYSBM/n525/scueal.csv")
#this is the scueal for all 3900 though   
head(scueal_orig)
sum(n525$Name %in% scueal_orig$Name)
scueal_orig525<-scueal_orig[scueal_orig$Name %in% n525$Name, ]
###Loop which edits the subtype
##remove intra-subtype breakpoints 
###And find the new breakpoints adjusted for gaps. 


subty2 <- subty
simpSub<-c()

for (i in 1:nrow(subty)){
  nam <- subty$Name[i]
  subt <- strsplit(subty$Subtype[i],",")[[1]]
  brk <- as.numeric(strsplit(subty$Breakpoints[i],";")[[1]])

  #In sequences with inter-SCUEAL.Subtype recomb,
  #get rid of intra-SCUEAL.Subtype recomb (compress SCUEAL.Breakpoints so only beteween SCUEAL.Subtype change, not within)
  #those with only inter-SCUEAL.Subtype recomb will retain their SCUEAL.Breakpoints
  if(length(subt)>1){
      keep <- 1
      for(j in 2:length(subt)){
        if(subt[j]!=subt[j-1]){#if the subtype is not the same as the one before it 
          keep <- c(keep, j) #add the subtype to the keep list 
        }
      } #keep is the NUMBER of subtypes to keep.e.g. h D A1 C D makes 1 2, 3, 4, 5
  keepBr <- keep-1 # this makes 0,1,2,3,4 (sequentially)
  subt <- subt[keep]#index of the subtypes to keep
  brk <- brk[keepBr]
  #brk <- brk[keepBr]#index of breakpoints to keep
    } else {
  #it's a pure subtype - but want to get rid of "(1 breakpoints)" etc
  subt <- gsub("[[:blank:]][(][0-9]+[[:blank:]][[:alpha:]]+[)]", "", subt)
  #intra<-which(grep("breakpoints",subt,fixed=T))
  brk<-"NA"#get rid of breakpoints - intra-subtype breakpoints unreliable 
  }
  # 
  cleanSub <- unique(subt)
  cleanSub <- cleanSub[order(cleanSub)]
  simpSub<-c(simpSub, paste(cleanSub, collapse=','))
  #totSubs<-c() #added?
  #totSubs <- c(totSubs, paste(cleanSub,collapse=","))
  #numBrkPt <- c(numBrkPt, length(brk))
  
  subty2[i,"Subtype"] <- paste(subt,collapse=",")
  subty2[i,"Breakpoints"] <- paste(brk,collapse=";")
  #newBrks <- c(newBrks, paste(newBrk, collapse=";"))
  
}

#new dataframe 
Genome_breakpoints <- as.data.frame(cbind(subty$Name, subty$Subtype, subty$Breakpoints, subty2$Subtype, subty2$Breakpoints,      simpSub,  subty$PDF))
colnames(Genome_breakpoints) <-            c("Name", "OriginalSub", "OriginalBreakpoints", "EditedSubtype", "EditedBreakpoints","CleanSubtype","PDF")

head(Genome_breakpoints)
hist(as.numeric(Genome_breakpoints$EditedBreakpoints))
unique(Genome_breakpoints$CleanSub)


head(scueal_orig)
scueal_orig<-scueal_orig[, c(1,2)]
colnames(scueal_orig)[2]<-"pol_scueal"
compare<-merge(Genome_breakpoints, scueal_orig)
write.csv(compare, "Scueal_525_full_length_versoin_edited.csv")
