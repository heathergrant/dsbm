#window exploration 
#heathergrant 10/06/21 
setwd("~/Dropbox/My Mac (sce-bio-c04853)/Desktop/MAPPED_DRIVE_COPY/dsbm")
#setwd
setwd("../DYSBM/UGANDA/genomes_200/tn93/")

files <- Sys.glob('*.tn93.csv')
#slices <- as.numeric(sapply(files, function(x) gsub(".+_([0-9]+).+", "\\1", x)))
slices <- as.integer(sapply(files, function(x) gsub(".*sw([0-9]+).*", "\\1", x)))
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
library(tidyr)
library(ggplot2)
library(dplyr)
df.dist2<-as.data.frame(df.dist)

df3<-df.dist2 %>% gather(key="Window", value='dist') 
library(forcats)  
df4<-df3 %>% mutate(Window = fct_relevel(Window, paste0("win_", 1:length(files))))

ggplot(df4, aes(x=fct_inorder(Window), y=as.numeric(dist)), fill='Window') +
  geom_violin()+
  theme(axis.text.x = element_text(angle=90))
  
#any(is.na(df4$dist))

