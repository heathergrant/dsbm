#make windows AND TN93 distances 
## R version of sliding_window.py PLUS tn93.R in one go - heather's version 
## making sliding windows 
library(ips) #for running mafft in R 
library(ape)
library(seqinr)
library(parallel)
setwd("~/Desktop/MAPPED_DRIVE_COPY/DYSBM")
input.fasta <- read.dna("../aln/allA1Dgenomes.fasta", format = "fasta", as.character=T, as.matrix=T)
input.fasta <- read.dna("aln/recombs164_plus_AxAyAz.fasta", format = "fasta", as.character=T, as.matrix=T)
input.fasta <- input.fasta[c(10,20,45,100, 134, 160),] #just do with 6 random sequences 
input.fasta <- input.fasta[100:151,] #just do with 51 random sequences 

windowsize <- 500
shift<-250
numSeq <- nrow(input.fasta)
window <- 1
numWindows <- seq(1,(length(input.fasta[1,])%/%shift)*shift-windowsize+1,by=shift)
dir.create("recombs_164")
setwd("recombs_164") #create this directory and set 
for (i in numWindows) {
  print(paste(Sys.time(), " Processing: ", as.character(window), "/" , as.character(length(numWindows)), sep=""))
  subseq <- as.DNAbin(input.fasta[,seq(i,i+windowsize-1)]) #cut window 
  subseqM<-mafft(subseq, thread=10, op=3, ep=0.123, maxiterate=1000) #align
  dir.create(paste0("win", as.character(window)))     #make directory 
  write.dna(subseqM, paste0("win", as.character(window), "/", "mafft_", as.character(window),  ".fasta"), format='fasta') #write file
  tn93 <- dist.dna(subseqM, model = "TN93", variance = FALSE,
                   gamma = FALSE, pairwise.deletion = FALSE,
                   base.freq = NULL, as.matrix = FALSE) #make tn93 distance 
  write.table(data.frame(as.matrix(tn93)),  
              file= paste0("win", as.character(window), "/", "tn93_", as.character(window),  ".csv"), sep = ",", col.names = NA)
  window <- window + 1
  
}
 



