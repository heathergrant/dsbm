#make frankenstein genome
#heather grant 15 june 2021 
setwd("~/Dropbox/My Mac (sce-bio-c04853)/Desktop/MAPPED_DRIVE_COPY/dsbm")

setwd("../DYSBM/walk_through/tn93/")
files <- Sys.glob('*.tn93.csv')
slices <- as.integer(sapply(files, function(x) gsub(".*sw([0-9]+).*", "\\1", x)))
files <- files[order(slices)]
infileG<-"../json/gene_dynamics.tsv"
infileJ<-"../json/bestoutput.json"

setwd("../mafft/")
files.dna <- Sys.glob('*.mafft')
slices.dna <- as.integer(sapply(files.dna, function(x) gsub(".*sw([0-9]+).*", "\\1", x)))
files.dna <- files.dna[order(slices.dna)]


plot_heatmap("../json/bestoutput.json", "../json/gene_dynamics.tsv", "../heatmap_new_colours.pdf")

tt <- read.csv(infileG, header = T, sep = "\t")
json <- read_json(infileJ)
mb <- json$membership
df <- t(as.data.frame(matrix(unlist(mb), nrow=length(unlist(mb[1])))))
rownames(df) <- tt$X
#rownames(df)<-gsub("Ref.", "", rownames(df))
df<-df[order(rownames(df)),]   
m <- json$trans
n<-length(m)
df[1,] #so this is the first sequence and it's membership along the genome 
df[,1] #this is the first window membership 

length(unique(df[,1]))#this is the number of k along the genome - 0 means no membership - so white in figure
k_along<-data.frame(win=1:ncol(df), k=rep(0,ncol(df)))
for(i in 1:ncol(df)){
  k_along[,i]<-length(unique(df[,i]))
}
k_along
library(ape)

setwd("../Frankenstein/")

#make dna object to start off with. 
starter<-read.dna("../ref_genomes.afa", format='fasta')
starter<-starter[1,]#jsut the first one
library(seqinr)
write.fasta(rep('-', dim(starter)[2]), 'blank', 'blank-dash.fasta', open = "w", nbchar = 60, as.string = FALSE)
blank<-read.dna("blank-dash.fasta", format='fasta')
numWindows <- seq(1,((length(starter[1,])+500)%/%100)*100-500+1,by=100)-1
#starter[,1]<-#list(rep("n", 3))# dim(starter)[2]
#loop through windows - let's just make community 1 to begin with 
#Frankenstein<-starter
k<-1 #first community 
for(k in 1:n){
  rm(Frankenstein)
  Frankenstein<-starter
  for(w in 1:ncol(df)){
    
    kw<-df[,w][df[,1]==k] #only 5 sequences 
    #which genome of this lot has the most '1's elsewhere in the genome? 
    mem1<-df[names(kw),] #subset to just ones of interest 
    counts<-c()
    for(c in 1:length(rownames(mem1))){
      counts<-c(counts, sum(mem1[c,]==k))
    }
    top<-which(counts==max(counts))
    top<-top[sample(length(top), 1)] #pick one if there's more than one (will have to revisit this, quick fix)
    frank<-rownames(mem1)[top] #4th sequence is best representative because it has other 31 windows elsewhere in its genome that are also this community 
    dna.1<- read.dna(files.dna[w], format='fasta')
    frank.dna<-dna.1[frank,]
    rownames(frank.dna)<-paste0("community_", k)
    #some windows are 503 base pairs long! Presumably because of re-maffting them. 
    frank.dna<-frank.dna[,1:500]
    start<- numWindows[w]    #       or #100*(w+1)-200#start 
    end<- numWindows[w+5]+1   #      or100*(w+1)+300 #end
    blank1<-blank[,0:start]
    rownames(blank1)<-paste0("community_", k)
    blank2<-blank[,end:dim(starter)[2]]
    rownames(blank2)<-paste0("community_", k)
    rownames(frank.dna)<-paste0("community_", k)
    #stitch blank + window + blank sequence to make an alignment of equal length Frankenstein to make cons
    record.dna<-cbind.DNAbin(blank1,frank.dna)
    record.dna<-cbind.DNAbin(record.dna, blank2)

    Frankenstein<-rbind.DNAbin(Frankenstein, record.dna)
  }
  Frankenstein<-Frankenstein[-1,]#remove starter seq 
  write.dna(Frankenstein, paste0("Community_", k, ".fasta"), format='fasta')

}


#blast approach to frankenstein references 
library(ips) #for running mafft in R 
library(ape)
library(seqinr)
library(parallel)
#Q<-read.dna("../../hist_plus_ref/aln/hist_109_plus_ref_mafft_eye_leave_tough_bits_tackle.fasta", format='fasta')
#Q<-Q[1,] #query sequence 
Q<-read.dna("../Frankenstein/Query.fasta", format='fasta')#aligned to 27 reference genomes so should fit 
w<-1
query_community<-data.frame(win=1:ncol(df), k=rep(0,ncol(df)), ref=rep('U', ncol(df)))
for(w in 1:ncol(df)){ #in each window 
      
    win.dna<- read.dna(files.dna[w], format='fasta') 
    win.Q<-Q[,numWindows[w]:numWindows[w+5]-500+ncol(win.dna)]
    together<-rbind(win.dna, win.Q)
    align<-mafft(together, thread=10, op=3, ep=0.123, maxiterate=1000) #align
    tn93 <- dist.dna(align, model = "TN93", variance = FALSE,
                     gamma = FALSE, pairwise.deletion = FALSE,
                     base.freq = NULL, as.matrix = TRUE) #make tn93 distance 
    d<-as.data.frame(tn93) #distances in a data_frame for easier manipulation
    d$Seq1<-rownames(d) # create new column with seq1 name
    table<-gather(d, Seq2, distance, 1:38) # turns distance matrix into a table of pairwise comparisons
    table2<- dplyr::filter(table, Seq1 != Seq2) #remove ones where seqs being compared with themselves
    top<-table2[order(table2$distance),] #order by distance
    top<-top[1,]$Seq2 #this is the closest seuqence in the reference alignment. Find out what community it belongs 
    df[top,1]#widow 1, best match is in community 6 
    k<-df[top,1]#widow 1, best match is in community 6 
    query_community[w,2]<-k
    query_community[w,3]<-top
}
    
####

####
    
    kw<-df[,w][df[,1]==k] #only 5 sequences 
    #which genome of this lot has the most '1's elsewhere in the genome? 
    mem1<-df[names(kw),] #subset to just ones of interest 
    counts<-c()
    for(c in 1:length(rownames(mem1))){
      counts<-c(counts, sum(mem1[c,]==k))
    }
    top<-which(counts==max(counts))
    top<-top[sample(length(top), 1)] #pick one if there's more than one (will have to revisit this, quick fix)
    frank<-rownames(mem1)[top] #4th sequence is best representative because it has other 31 windows elsewhere in its genome that are also this community 
    dna.1<- read.dna(files.dna[w], format='fasta')
    frank.dna<-dna.1[frank,]
    rownames(frank.dna)<-paste0("community_", k)
    #some windows are 503 base pairs long! Presumably because of re-maffting them. 
    frank.dna<-frank.dna[,1:500]
    start<- numWindows[w]    #       or #100*(w+1)-200#start 
    end<- numWindows[w+5]+1   #      or100*(w+1)+300 #end
    blank1<-blank[,0:start]
    rownames(blank1)<-paste0("community_", k)
    blank2<-blank[,end:dim(starter)[2]]
    rownames(blank2)<-paste0("community_", k)
    rownames(frank.dna)<-paste0("community_", k)
    #stitch blank + window + blank sequence to make an alignment of equal length Frankenstein to make cons
    record.dna<-cbind.DNAbin(blank1,frank.dna)
    record.dna<-cbind.DNAbin(record.dna, blank2)
    
    Frankenstein<-rbind.DNAbin(Frankenstein, record.dna)
  }
  Frankenstein<-Frankenstein[-1,]#remove starter seq 
  write.dna(Frankenstein, paste0("Community_", k, ".fasta"), format='fasta')
  
}






#probs best to check by eye anyway? 
#image(Frankenstein)
#Frankenstein<-Frankenstein[-1,]#remove starter seq 
#Frankenstein_condensed<-data.frame(1:dim(Frankenstein)[2])
#df <- data.frame(ID=labels(Frankenstein),
#                 seq=sapply(Frankenstein, paste, collapse=","))

#for(f in 1:dim(Frankenstein)[2]){
#  nt<-unique(Frankenstein[,f])
#  for(n in 1:length(nt)){
#  Frankenstein_condensed[n,f]<-as.character(nt[n])
#  }
#}
#BiocManager::install("DECIPHER")
#library(DECIPHER)
#dnastring<-DNAStringSet(matailn)
#ConsensusSequence(matailn,
#
#Frankenstein<-as.matrix.alignment(Frankenstein)
#write.dna(Frankenstein, "Community1.fasta", format='fasta')
#matailn<-read.alignment("Community1.fasta", format='fasta')
#com1<-consensus(matailn, method =  "threshold", warn.non.IUPAC = FALSE, type = "DNA", threshold = 1)
#com1
#probably should doubel check by eye anyway. Make consensus in geneious... can't get consensus in R without gaps very easily. 



