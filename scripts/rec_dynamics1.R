library(igraph)
library(dynsbm)
#library(RJSONIO)
library(jsonlite)

#####################################################
# params to change to tn93 folder 
setwd("")
#number of cores to use 
cores<-1
#range of K to try out (1 to 10)
Qmin<-1
Qmax<-10
#####################################################

#start time
start.time <- Sys.time()

# load TN93 files
files <- Sys.glob('*.tn93.csv')
#slices <- as.numeric(sapply(files, function(x) gsub(".+_([0-9]+).+", "\\1", x)))
slices <- as.integer(sapply(files, function(x) gsub(".*sw([0-9]+).*", "\\1", x)))
files <- files[order(slices)]

# Read TN93 files and get the maximum distance value
list.m <- lapply(files, function(f) {
  x <- read.table(f, sep=',', header=T, na.strings=c("NA", "NaN"), row.names=1)
  quantile(as.matrix(x), probs = 0.4, na.rm = TRUE)
})

list.g <- vector("list", length = length(list.m))
for (v in 1:length(list.m)) {
  f <- read.table(files[v], sep=',', header=T, na.strings="NA", row.names=1)
  co <- (as.numeric(list.m[[v]]))
  list.g[[v]] <- graph_from_adjacency_matrix(f <= co, weighted=TRUE, diag=FALSE, mode='undirected')
}  

# Print graphs
#for (g in 1:length(list.g)) {
# plot(list.g[[g]], vertex.size = 3, vertex.label = NA)
#}


# converting to R array
all.vertices.name <- unique(sort(unlist(lapply(list.g, function(g) V(g)$name))))
T <- length(list.g)
N <- length(all.vertices.name)
Y <- array(dim=c(T,N,N))

for (t in 1:T){
  g.tmp <- graph.empty(n=N, directed=F)
  V(g.tmp)$name <- all.vertices.name
  Y[t,,] <- as.matrix(get.adjacency(union(g.tmp,list.g[[t]])))
}


compute.icl <- function(dynsbm){    
  T <- ncol(dynsbm$membership)
  Q <- nrow(dynsbm$trans)
  N <- nrow(dynsbm$membership)
  pen <- 0.5*Q*log(N*(N-1)*T/2) + 0.25*Q*(Q-1)*T*log(N*(N-1)/2) # binary case
  if ("sigma" %in% names(dynsbm)) pen <- 2*pen # continuous case
  return(dynsbm$loglikelihood - ifelse(T>1,0.5*Q*(Q-1)*log(N*(T-1)),0) - pen)    
}

print("Got here ........................................................ list.dynsbm")

# dynSBM analysis
list1.dynsbm <- select.dynsbm(Y, Qmin=1, Qmax=31, edge.type="binary", directed=FALSE, self.loop=FALSE,
                              nb.cores=12, iter.max=20, nstart=25, perturbation.rate=0.2, fixed.param=FALSE,
                              bipartition=NULL, plot=TRUE)
#save(list1.dynsbm, file = "dynsbm.Rdata")

# manual ICL and loglikelihood plot
Qmin <- ncol(list1.dynsbm[[1]]$trans)
Qmax <- ncol(list1.dynsbm[[length(list1.dynsbm)]]$trans)
dir.create("../json")
setwd("../json")
#create json files for each group
for (Q in Qmin:Qmax){
  json <- toJSON(list1.dynsbm[[Q]], pretty = TRUE)
  write(json, file= paste0(Q,"_dynsbm.json"))
}

#analyze the group with the highest ICL
dynsbm <- list1.dynsbm[[which.max(sapply(list1.dynsbm, compute.icl))]]
# == which.max(sapply(list1.dynsbm, compute.icl)) == 6
#useful for later to manually map out seq names to data 
colnames(dynsbm$membership) <- paste0("window_", 1:ncol(dynsbm$membership))
write.table(data.frame(dynsbm$membership), file = "gene_dynamics.tsv", row.names = all.vertices.name, quote = FALSE, sep = "\t", col.names = NA)
#json <- toJSON(list1.dynsbm[[6]], pretty = TRUE)
#write(dynsbm, file = "bestoutput.json")
write(toJSON(dynsbm, pretty=TRUE), file = "bestoutput.json")

# stop and calculate time
end.time <- Sys.time()
time.taken <- end.time - start.time
xx <- time.taken
print(xx)

write(paste(xx, attr(xx, "units")), file = paste0("time.txt"))
