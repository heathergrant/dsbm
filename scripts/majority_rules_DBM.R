#majority rules DBM ? 
#heather grant 13 July 2021 
setwd("~/Dropbox/My Mac (sce-bio-c04853)/Desktop/MAPPED_DRIVE_COPY/dsbm")

setwd("../DYSBM/UGANDA/D_intra/json_40pc/") #df1 
infileG<-"gene_dynamics.tsv"
infileJ<-"bestoutput.json"
tt <- read.csv(infileG, header = T, sep = "\t")
json <- read_json(infileJ)
mb <- json$membership
df <- t(as.data.frame(matrix(unlist(mb), nrow=length(unlist(mb[1])))))
rownames(df) <- tt$X
df<-df[order(rownames(df)),]   
m <- json$trans
n<-length(m)
df[1,] #so this is the first sequence and it's membership along the genome 
df[,1] #this is the first window membership 
df1<-df
length(unique(df[,1]))#this is the number of k along the genome - 0 means no membership - so white in figure
k_along<-data.frame(win=1:ncol(df), k=rep(0,ncol(df)))
for(i in 1:ncol(df)){
  k_along[,i]<-length(unique(df[,i]))
}
k_along

####
greekletters<-c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma")
numbers1<-c(1:n, 0) #add case where k is unassinged = 0 
newlabels<-c(greekletters[1:n], "blank")
reps<-c()
for(k in 1:n){
  counts<-c()
  for(c in 1:length(rownames(df))){
    counts<-c(counts, sum(df[c,]==k))
  }
  top<-which(counts==max(counts))
  top<-top[sample(length(top), 1)] #pick one at random if there's more than one (will have to revisit this, quick fix)
  rep<-rownames(df)[top] #this is the best representative genome of group k of n 
  
  reps<-c(reps, rep)
}
reps<-c(reps, "blank")
bestRefs<-data.frame(df1_label=numbers1,newlabels=newlabels, BestRep=reps)

#then repeat - read in another one #so let's try another one where n=7
setwd("../DYSBM/UGANDA/D_intra/json_30pc/") #df2
infileG<-"gene_dynamics.tsv"
infileJ<-"bestoutput.json"
tt <- read.csv(infileG, header = T, sep = "\t")
json <- read_json(infileJ)
mb <- json$membership
df <- t(as.data.frame(matrix(unlist(mb), nrow=length(unlist(mb[1])))))
rownames(df) <- tt$X
df<-df[order(rownames(df)),]   
m <- json$trans
n<-length(m)
df2<-df

#now. translate into alpha etc  
df2[rownames(df2)==bestRefs[1,3],] #this is genome best looking like alpha community 
library(DescTools)
#mode<-Mode(df2[rownames(df2)==bestRefs[1,2],])[1] #alpha translates as 1 
modes<-c() #this is the blank one 
for(k in 1:n){
  df2[rownames(df2)==bestRefs[k,3],] #this is genome best looking like alpha community 
  mode<-Mode(df2[rownames(df2)==bestRefs[k,3],])[1] #alpha translates as 4 
  modes<-c(modes, mode)
}
modes<-c(modes, 0)
bestRefs<-cbind(bestRefs, modes)
colnames(bestRefs)[4]<-"df2_label"
#okay so now komunity 1 and 4 best match both look like df1 number 2... 



#okay so how to merge then ?
df1g<-df1
#convert into alpha beta gamma first? 
for(i in 1:nrow(df1g)){
  for(j in 1:ncol(df1g)){
    for(k in 1:n){
    if(df1g[i,j]==as.character(k)){
      df1g[i,j]<-bestRefs[bestRefs$newlabels==as.character(k),2]
      }
    }
  }
}
df1g
df2g<-df2
#now translate backwards - three in this context is zeta 
for(i in 1:nrow(df2g)){
  for(j in 1:ncol(df2g)){
    for(k in 1:n){
      if(df2g[i,j]==as.character(k)){
        df2g[i,j]<-bestRefs[bestRefs$df2_label==as.character(k),2]
      }
    }
  }
}

df2
df1
#how to merge? 
df3<-df #copy of original 1 with numbers not greek letters 
for(i in 1:nrow(df3)){
  for(j in 1:ncol(df3)){
    for(k in 1:n){
      if(df1[i,j]==df2[i,j]){
        df2[i,j]<-greekletters[k]
      }
    }
  }
}


