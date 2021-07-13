#make windows AND TN93 distances 
## R version of sliding_window.py PLUS tn93.R in one go - heather's version 
## making sliding windows 
library(ips) #for running mafft in R 
library(ape)
library(seqinr)
library(parallel)
setwd("~/Desktop/MAPPED_DRIVE_COPY/DYSBM/UGANDA/D_intra")
input.fasta <- read.dna("D_genomes_aligned_to_win81_pre-00s_no=ref.fasta", format = "fasta", as.character=T, as.matrix=T)
windowsize <- 500
shift<-100
numSeq <- nrow(input.fasta)
window <- 1
numWindows <- seq(1,(length(input.fasta[1,])%/%shift)*shift-windowsize+1,by=shift)
dir.create("refs_41_pairwiseF")
setwd("refs_41_pairwiseF") #create this directory and set 
dir.create("mafft")
dir.create("tn93")
for (i in numWindows) {
  print(paste(Sys.time(), " Processing: ", as.character(window), "/" , as.character(length(numWindows)), sep=""))
  subseq <- as.DNAbin(input.fasta[,seq(i,i+windowsize-1)]) #cut window 
  #### Maybe remove mafft step because I have already very carefully aligned the windows? 
  subseqM<-mafft(subseq, thread=10, op=3, ep=0.123, maxiterate=1000) #align
  ####
  #dir.create(paste0("win", as.character(window)))     #make directory 
  write.dna(subseqM, paste0("mafft/win_", as.character(window),  ".fasta"), format='fasta') #write file
  tn93 <- dist.dna(subseqM, model = "TN93", variance = FALSE,
                   gamma = FALSE, pairwise.deletion = TRUE, #i think pairwise.del = FAL was causing NA 
                   base.freq = NULL, as.matrix = FALSE) #make tn93 distance 
  write.table(data.frame(as.matrix(tn93)),  
              file= paste0("tn93/win_", as.character(window), ".tn93.csv"), sep = ",", col.names = NA)
  window <- window + 1
  
}
 

