#library(ape)
#library(gplots)
require(jsonlite)
#library(RColorBrewer)
#library(dplyr)

getwd()
setwd("../../genomes_500_bioinf7/json")
plot_heatmap<-function(infileJ, infileG, outfile){

tt <- read.csv(infileG, header = T, sep = "\t")
#nam <- lapply(strsplit(as.character(tt$V1), "_"), "[", 2)
#nam <- gsub("Ref_", "", as.character(tt$V1))

json <- read_json(infileJ)

#json <- read_json("../../recombs_50_164/bestoutput.json")
mb <- json$membership
df <- t(as.data.frame(matrix(unlist(mb), nrow=length(unlist(mb[1])))))
rownames(df) <- tt$X
rownames(df)<-gsub("Ref.", "", rownames(df))
df<-df[order(rownames(df)),]   
json$self.loop
m <- json$trans
n<-length(m)
m <- matrix(unlist(m), ncol=n, byrow=T)
m <- m + t(m)  # convert transition rates into a symmetric distance matrix
m <- m/2
hc1 <- hclust(as.dist(log(m), diag=FALSE), method='ward.D2')

# brings together clusters linked by higher transition rates
idx <- order(hc1$order)
idx
# draw heatmap
pdf(outfile, width=5, height=6)
par(mar=c(1,1,1,1))
pal <- hcl.colors(n, 'Sunset')
heatmap(as.matrix(df), Rowv=T, Colv=NA, hclustfun = hclust,
        labCol=T, scale='none', col=c('white', pal[idx]),
        cexCol = 0.3, cexRow = 0.45,
        xlab="", ylab="",
        cex.lab = 1.3) -> h
dev.off()
print(paste0("value of K is" , n))
}