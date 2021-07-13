#window exploration 
#heathergrant 10/06/21 
library(tidyr)
library(ggplot2)
library(dplyr)
library(forcats) 
setwd("~/Dropbox/My Mac (sce-bio-c04853)/Desktop/MAPPED_DRIVE_COPY/dsbm")
#setwd
setwd("../DYSBM/UGANDA/D_intra/tn93/")
setwd("../DYSBM/walk_through/tn93/")
setwd("../DYSBM/NEW_ALIGN/500_new/tn93/")
setwd("../44_refs2/tn93/")
setwd('tn93')
files <- Sys.glob('*.tn93.csv')
slices <- as.integer(sapply(files, function(x) gsub(".+_([0-9]+).*$", "\\1", x)))
files <- files[order(slices)]

# Read TN93 files and get distances 
list.dist <- lapply(files, function(f) {
  x <- read.table(f, sep=',', header=T, na.strings=c("NA", "NaN"), row.names=1)
  as.list(as.numeric(unlist(x[lower.tri(x)])))
})

df<-do.call(rbind,list.dist)
dim(df)
df.dist<-t(df)
colnames(df.dist)<-paste0("win_", 1:length(files))
df.dist2<-as.data.frame(df.dist)
df3<-df.dist2 %>% gather(key="Window", value='dist') 
df4<-df3 %>% mutate(Window = fct_relevel(Window, paste0("win_", 1:length(files))))
df4<-df4 %>% na.omit() #na in here - but why would there be na? gaps in seq?
ggplot(df4, aes(x=fct_inorder(Window), y=as.numeric(dist)), fill='Window') +
  geom_violin(fill='coral')+
  theme(axis.text.x = element_text(angle=90))+
  labs(title='diversity tn93 in windows of D intra genomes pre 2000', x="", y="tn93")
  
##
## rectangle annotations 
hmin=0.5
plot.new()
locator(1)
rect(0, hmin+0.1, 396, hmin+0.15, lwd=1) #, col="dodgerblue1") #gag 
rect(397, hmin+0.1, 1304, hmin+0.15, lwd=1) #, col="dodgerblue2") #gag 
rect(1305, hmin+0.1, 1500, hmin+0.15, lwd=1) #, col="dodgerblue3") #gag 

#rect(1302, hmin+0.05, 4317, hmin+0.1, lwd=1) #, col='darkorange') #pol 
rect(1302, hmin, 4317, hmin+0.05, lwd=1) #, col='darkorange') #pol 
rect(1472, hmin, 1769, hmin+0.05, lwd=1) #, col='darkorange1') #pol 
rect(1770, hmin, 3450, hmin+0.05, lwd=1) #, col='darkorange2') #pol 
rect(3451, hmin, 4317, hmin+0.05, lwd=1) #, col='darkorange3') #pol 

rect(4265, hmin+0.1, 4842, hmin+0.15, lwd=1) #, col='gray') #vif 
rect(4785, hmin, 5072, hmin+0.05, lwd=1) #, col = "dodgerblue1") #vpr 
rect(5057, hmin+0.05, 5273, hmin+0.1, lwd=1) #, col="magenta4")#tat1
rect(7264, hmin+0.1, 7305, hmin+0.15, lwd=1) #, col="magenta4")#tat2 
rect(5290, hmin+0.05, 5534, hmin+0.1, lwd=1) #, col= "brown")# vpu
rect(5198, hmin+0.0, 5273, hmin+0.05, lwd=1) #, col= "chartreuse4")# rev1
rect(7264, hmin+0.05, 7538, hmin+0.1, lwd=1) #, col= "chartreuse4")# rev2

rect(5455, hmin, 6644, hmin+0.05, lwd=1) #, col="coral1") #gp120 
rect(6645, hmin, 7680, hmin+0.05, lwd=1) #, col="coral2") #gp41
rect(7579, hmin+0.1, 8292, hmin+0.15, lwd=1) #, col="cyan4") #nef


text(200, (hmin+0.125),"p17", col="black", cex=1)
text(850, (hmin+0.125),"p24", col="black", cex=1)

text(1600, (hmin+0.025),"prot", col="black", cex=1)
text(2500, (hmin+0.025),"RT", col="black", cex=1)
text(3900, (hmin+0.025),"int", col="black", cex=1)
text(4500, (hmin+0.125),"vif", col="black", cex=1)
text(4900, (hmin+0.025),"vpr", col="black", cex=1)

#text(5150, (hmin+0.11), "tat", col='magenta4', cex=1)
#text(7150, (hmin+0.125), "tat", col='magenta4', cex=1)
text(5150, (hmin+0.075), "tat", col='black', cex=1)
text(7400, (hmin+0.125), "tat", col='black', cex=1)
#curveGrob(5150, hmin+0.11, 7200, hmin+0.125)

#text(5400, (hmin+0.108), "vpu", col='brown', cex=1)

#text(5250, (hmin-0.0125), "rev", col='chartreuse4', cex=1)
#text(7150, (hmin+0.09), "rev", col='chartreuse4', cex=1)
text(5400, (hmin+0.075), "vpu", col='black', cex=1)

text(5250, (hmin-0.0125), "rev", col='black', cex=1)
text(7400, (hmin+0.075), "rev", col='black', cex=1)
text(7900, (hmin+0.125),"nef", col="black", cex=1)

text(6100, (hmin+0.025),"g120", col="black", cex=1)
text(7100, (hmin+0.025),"gp41", col="black", cex=1)



####end figure 
#extra text 

###



win17<-read.table(files[17], sep=',', header=T, na.strings=c("NA", "NaN"), row.names=1)
win17<-as.matrix(win17)
heatmap(win17)

heatmap(win17[450:500,450:500])

heatmap(as.matrix(df.dist2), Rowv=T, Colv=NA, hclustfun = hclust,
        labCol= T, scale='none', col=c('white', pal[idx]),
        cexCol = 0.3, cexRow = 0.45,
        xlab="", ylab="",
        cex.lab = 1.3) -> h
#any(is.na(df4$dist))
length(y[is.na(y)]) #1437 ? why - big gap in maj of seqs 
##number of K 
#par(mfrow = c(5, 1))
getwd()
infileG<-("gene_dynamics.tsv")
infileJ<-("bestoutput.json")
plot_heatmap(infileJ, infileG, "../refs41_pairwiseF_new_colours.pdf")
rm(df)
json <- read_json(infileJ)
tt <- read.csv(infileG, header = T, sep = "\t")
mb <- json$membership
df <- t(as.data.frame(matrix(unlist(mb), nrow=length(unlist(mb[1])))))
rownames(df) <- tt$X
#rownames(df)<-gsub("Ref.", "", rownames(df))
df<-df[order(rownames(df)),]   

df[1,] #so this is the first sequence and it's membership along the genome 
df[,1] #this is the first window membership 

length(unique(df[,1]))#this is the number of k along the genome - 0 means no membership - so white in figure
k_along<-data.frame(win=1:ncol(df))#, k=rep(0,ncol(df)))
k_along$k<-""
for(i in 1:ncol(df)){
  k_along[i,2]<-length(unique(df[,i]))
}
k_along
pdf('../K_along.pdf', width=5, height=6)
plot(k_along)
dev.off()
dim(k_along)



