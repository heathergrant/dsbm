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

#okay, let's look at cluster 1 in window 1 and get a central node? graph concept of centrality 
k1w1<-df[,1][df[,1]=='1'] #only 5 sequences 
?igraph_closeness
df[,1]=='1'
names(k1w1)
#files[1]#
#x <- read.table(files[1], sep=',', header=T, na.strings=c("NA", "NaN"), row.names=1)
#x1<-x[names(k1w1),names(k1w1)]
#library(igraph)
#y<-graph_from_data_frame(x1, directed = TRUE, vertices = NULL)
#plot(y)
#x1
#closeness(y, vids = V(y), mode = c("out", "in", "all",
#                                           "total"), weights = NULL, normalized = FALSE)
#not sure this works very well

#which genome of this lot has the most '1's elsewhere in the genome? 
mem1<-df[names(k1w1),]
counts<-c()
for(c in 1:length(rownames(mem1))){
  counts<-c(counts, sum(mem1[c,]=='1'))
}
frank<-rownames(mem1)[which(counts==max(counts))] #4th sequence is best representative because it has other windows (31) that are also this community 
#"Ref.G.PT.x.PT2695.AY612637"
library(ape)
dna.1<- read.dna(files.dna[1], format='fasta')
frank.dna<-dna.1[frank,]

