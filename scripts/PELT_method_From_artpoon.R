

#library(randomcoloR)
library(jsonlite)
set.seed(1)
#setwd("~/Desktop/git/skycross/paper/data_scripts/")
setwd("~/git/DYSBM/walk_through/json_1/")

gene <- read.table('gene_dynamics.tsv', sep='\t', header=T)
cahr <- read.csv('cahr_updated.csv')

require(ape)
fasta <- read.FASTA("../../data/rstep_0.fa")

# map subtypes to gene data frame, which should map 1-to-1 to dynsbm JSON
cahr$accno <- gsub("^(.+)_.+$", "\\1", cahr$names)


gene$subtype <- cahr$subtype[match(gene$pos, cahr$accno)]
gene$is.recomb <- cahr$is.recomb[match(gene$pos, cahr$accno)]
gene$year <- cahr$year[match(gene$pos, cahr$accno)]

# 3 missing entries
cahr$names[match(gene$pos, cahr$accno)][490:492]
gene$pos[490:492]


# read json file containing community memberships at the highest ICL
# extract community data and convert to dataframe
json <- read_json(paste0("original_data/dynsbm_og/25_dynsbm.json"))
mb <- json$membership
df <- t(as.data.frame(matrix(unlist(mb), nrow=length(unlist(mb[1])))))



# convert transition rates into a symmetric distance matrix
#ds <- fromJSON('original_data/dynsbm_og/25_dynsbm.json')
m <- matrix(unlist(json$trans), ncol=25, byrow=T)
#m2 <- matrix(0, nrow=nrow(m), ncol=ncol(m))
m <- m + t(m)
m <- m/2
#m[upper.tri(m)] <-m[upper.tri(m)]/2
#m[lower.tri(m)] <- 0
#m[lower.tri(m)] <- t(m)[lower.tri(m)]
#m <- 1/m  # mean waiting time = 1/rate
#diag(m) <- 0

# hierarchical clustering on distance matrix
hc1 <- hclust(as.dist(log(m), diag=FALSE), method='ward.D2')

# brings together clusters linked by higher transition rates
idx <- order(hc1$order)

plot(hclust(as.dist(log10(m)), method='ward.D2'))

barplot(table(unlist(df[gene$subtype=='A1',])), las=2)
barplot(table(unlist(df[gene$subtype=='B',])), las=2)
barplot(table(unlist(df[gene$subtype=='C',])), las=2)
barplot(table(unlist(df[gene$subtype=='F',])), las=2)

#plot(NA, xlim=c(0,26), ylim=c(1,9))
subs <- c('A1', 'AE', 'B', 'C', 'D', 'F1', 'G', 'H', 'J', 'K', 'Complex')
mx <- matrix(NA, nrow=length(subs), ncol=26)
for (sb in 1:length(subs)) {
  y <- table(factor(df[gene$subtype==subs[sb]], levels=c(0:25)))
  y.denom <- 82*sum(gene$subtype==subs[sb], na.rm=T)
  for (j in 0:25) {
    z <- 1-y[[as.character(j)]]/y.denom
    mx[sb,(j+1)] <- z
    #rect(xleft=j, xright=j+1, ybottom=sb, ytop=sb+1, 
    #     col=rgb(z,z,z), border=NA)
  }
}
rownames(mx) <- subs
colnames(mx) <- as.character(0:25)
heatmap(mx)


sub.A1 <- t(apply(df[gene$subtype=='A1',], 1, function(x) {
  sapply(x, function(y) ifelse(is.na(y), NA, idx[y]))
}))
plot(NA, xlim=c(1, 82), ylim=c(1, 25))
for (i in 1:nrow(sub.A1)) {
  lines(sub.A1[i,], col=rgb(1,0,0,0.1))
}

# define number of communities (k)
k <- 25

#palette <- distinctColorPalette(k, runTsne = T)
#palette1 <- sort(palette)
#palette2 <- palette1[hc1$order]
#typ = sapply(strsplit(as.character(df$name),"_"), `[`, 2)

#pdf(file='heatmap_og1.pdf', width=5, height=5)
#par(mar=c(2,2,2,2))
#heatmap(as.matrix(df), Rowv=T, Colv=NA, hclustfun = hclust,
#reorderfun = function(d, w) reorder(d, w),
#labRow = NA, labCol=T, scale='column', col=palette, cexCol = 0.4)
#dev.off()

pal <- c("#6867C6", "#80A390", "#8AFBA3", "#DBCEE5", "#69E4BC",
         "#D087CF", "#CFC386", "#E5FB86", "#C0B3EC", "#BAE5F7", 
         "#DB4FD7", "#D1DBAF", "#FDFA60", "#DF68B6", "#89F2F2",
         "#ECB75B", "#E5D8CE", "#0885BA", "#AB68FD", "#B9EEA7",
         "#F07383", "#87F164", "#F4A588", "#6EE1FE", "#C78DAA")

#cola <- ifelse(df == 0,'black', pal)

#pdf(file='heatmap_og1.pdf', width=5, height=5)
d <- dist(as.matrix(df))

# generates a linear gradient by interpolating colours
#pal <- colorRampPalette(c('firebrick', 'goldenrod', 'forestgreen', 'dodgerblue', 'black'))(25)
pal <- hcl.colors(n=25, palette='viridis')
require(dichromat)
pal <- dichromat(pal, type='deutan')
require(ggfree)
#pal <- gg.rainbow(30, c=100, l=100)


pdf(file="~/git/skycross/paper/heatmap-accessible.pdf", width=5, height=5)
heatmap(as.matrix(df), Rowv=T, Colv=NA, hclustfun = hclust,
        #reorderfun = function(d, w) reorder(d, w),
        labRow = NA, labCol=T, scale='none', col=c('white', pal[idx]), cexCol = 0.5)
legend(x=0.9, y=1.1, legend=order(idx), fill=pal, bty='n', xpd=NA, 
       cex=0.6, border='white', x.intersp=0.1)
par(xpd=NA)
draw.hiv(x0=0.048, x1=0.91, y0=1.12, y1=1.22, label=T, cex.lab=0.6)
par(xpd=F)
dev.off()


traces <- apply(as.matrix(df), 1:2, function(i) ifelse(i==0, NA, idx[i]))

pdf(file='~/git/skycross/paper/traces-small.pdf', width=6, height=6)
par(mar=c(1,1,1,1), mfrow=c(4,4))
for (i in 1:16) {
  tr <- traces[i,]
  plot(tr, type='n', ylim=c(1, 25), 
       main=i, adj=0, cex.main=1.2, font.main=1,
       bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
  res <- cpt.mean(tr[!is.na(tr)], minseglen=10)
  na.count <- cumsum(is.na(tr))
  x <- c(0, res@cpts + na.count[res@cpts])
  y <- res@param.est$mean
  for (j in 1:length(y)) {
    segments(x0=x[j], x1=x[j+1], y0=y[j], col='salmon', lwd=2)  
  }
  lines(tr, col='grey20', lwd=1)
}
dev.off()


# deal with missing data by dropping those windows
# limit number of missing values to 10 (no more than)
plot(names(x), cumsum(x)/525, type='b')
text(as.integer(names(x))-1, cumsum(x)/525, label=names(x), cex=0.5)
abline(h=0.9)


require(changepoint)
breakpoints <- function(traces, max.missing=10, min.seg.len=10) {
  amoc <- rep(NA, nrow(traces))
  pelt <- rep(NA, nrow(traces))
  pelt.loc <- list()
  
  for (i in 1:nrow(traces)) {
    tr <- traces[i,]
    if (sum(is.na(tr)) > max.missing) {
      next
    }
    tra <- tr[!is.na(tr)]  # simply drop missing values
    res <- cpt.mean(tra, method='PELT', minseglen = min.seg.len)
    pelt[i] <- length(res@cpts)-1  # number of breakpoints
    
    # map back to original coords (if any NAs)
    loc <- res@cpts[-length(res@cpts)]
    na.count <- cumsum(is.na(tr))
    pelt.loc[[i]] <- loc + na.count[loc]
    
    res <- cpt.mean(tra, minseglen=min.seg.len)  # default method, AMOC
    amoc[i] <- length(res@cpts)-1
  }
  list('amoc'=amoc, 'pelt'=pelt, 'pelt.loc'=pelt.loc)
}


res.2 <- breakpoints(traces, min.seg.len=2)
res.3 <- breakpoints(traces, min.seg.len=3)
res.5 <- breakpoints(traces, min.seg.len=5)
res.10 <- breakpoints(traces, min.seg.len=10)

write.csv(data.frame(
  nbpt3=res.3[['pelt']], 
  nbpt5=res.5[['pelt']], 
  nbpt10=res.10[['pelt']]
),
file='breakpoints.csv', quote=F
)


boxplot(split(res.3[['pelt']], gene$is.recomb))
boxplot(split(res.5[['pelt']], gene$is.recomb), add=T)
boxplot(split(res.10[['pelt']], gene$is.recomb))
wilcox.test(res.3[['pelt']][gene$is.recomb], res.3[['pelt']][!gene$is.recomb])
wilcox.test(res.5[['pelt']][gene$is.recomb], res.5[['pelt']][!gene$is.recomb])
wilcox.test(res.10[['pelt']][gene$is.recomb], res.10[['pelt']][!gene$is.recomb])

plot(jitter(gene$year), jitter(res.3[['pelt']]))

require(MASS)
temp <- data.frame(year=gene$year, nbp=res.5[['pelt']])
tcc <- temp[complete.cases(temp), ]
k <- kde2d(tcc$year, tcc$nbp)

pdf(file='~/git/skycross/paper/time.pdf', width=5, height=5)
par(mar=c(5,5,1,1))
image(k, xlab='Year of sample collection', ylab='Number of breakpoints (DSBM)')
points(jitter(tcc$year), jitter(tcc$nbp), pch=19, cex=0.25)
dev.off()

summary(lm(nbp~year, data=tcc))



require(ggfree)


# draw HIV genome
draw.hiv <- function(x0, x1, y0, y1, label=F, cex.lab=1, ...) {
  width <- x1-x0
  height <- y1-y0
  map.x <- function(nt) {
    x0 + width * (nt-790)/(9719-790)  # scaled to (0,1)
  }
  map.y <- function(line) {
    y0 + height * line/3
  }
  # gag 
  rect(xleft=map.x(790), xright=map.x(2292), ybottom=map.y(2), ytop=map.y(3), ...)
  if (label) text(map.x(mean(c(790, 2292))), y=map.y(2.5), label='gag', cex=cex.lab)
  # pol
  rect(xleft=map.x(2085), xright=map.x(5096), ybottom=map.y(0), ytop=map.y(1), ...)
  if (label) text(map.x(mean(c(2085, 5096))), y=map.y(0.5), label='pol', cex=cex.lab)
  # vif
  rect(xleft=map.x(5041), xright=map.x(5619), ybottom=map.y(2), ytop=map.y(3), ...)
  # vpr
  rect(xleft=map.x(5559), xright=map.x(5850), ybottom=map.y(0), ytop=map.y(1), ...)
  # tat
  rect(xleft=map.x(5831), xright=map.x(6045), ybottom=map.y(1), ytop=map.y(2), ...)
  rect(xleft=map.x(8379), xright=map.x(8469), ybottom=map.y(2), ytop=map.y(3), ...)
  # rev
  rect(xleft=map.x(5970), xright=map.x(6045), ybottom=map.y(0), ytop=map.y(1), ...)
  rect(xleft=map.x(8379), xright=map.x(8653), ybottom=map.y(1), ytop=map.y(2), ...)
  # vpu
  rect(xleft=map.x(6062), xright=map.x(6310), ybottom=map.y(1), ytop=map.y(2), ...)
  # env
  rect(xleft=map.x(6225), xright=map.x(8795), ybottom=map.y(0), ytop=map.y(1), ...)
  if (label) text(map.x(mean(c(6225, 8795))), y=map.y(0.5), label='env', cex=cex.lab)
  # nef
  rect(xleft=map.x(8797), xright=map.x(9417), ybottom=map.y(2), ytop=map.y(3), ...)
}



pdf(file='~/git/skycross/paper/breakpoints.pdf', width=10, height=5)

par(mfrow=c(1,2), mar=c(5,5,1,1))

hist(res.5[['pelt']], right=F, border='white', xlim=c(0, 20), ylim=c(0, 0.35),
     col=rgb(0,0,0,0.5), breaks=0:13, main=NA,
     xlab='Number of breakpoints', ylab='Frequency', freq=F)
hist(res.10[['pelt']],  border='white', col=rgb(0,0,1,0.3), add=T, 
     breaks=0:7, right=F, freq=F)
hist(res.3[['pelt']], right=F, border='white', col=rgb(1,0,0,0.3), add=T, 
     breaks=0:20, freq=F)

legend(x=9, y=0.3, legend=c('10 windows', '5 windows', '3 windows'), 
       bty='n', border=NA, title='Minimum segment length',
       fill=c(rgb(0,0,1,0.3), rgb(0,0,0,0.75), rgb(1,0,0,0.3)))

denom <- apply(traces, 2, function(x) sum(!is.na(x)))

locs <- unlist(res.10[['pelt.loc']])
x1 <- sort(unique(locs))
y1 <- as.numeric(table(locs) / denom[min(locs):max(locs)])
plot(x=x1, y=y1, type='s', main=NA, xlim=c(1, 82),
     col=rgb(0,0,1,0.25), ylim=c(0, 0.18), 
     bty='n', xlab='Breakpoint location (window)', ylab='Frequency')
lines(smooth.spline(x1, y1), col='royalblue', lwd=2)

locs <- unlist(res.5[['pelt.loc']])
x2 <- sort(unique(locs))
y2 <- as.numeric(table(locs) / denom[min(locs):max(locs)])
lines(x=x2, y=y2, type='s', col=rgb(0,0,0,0.25))
lines(smooth.spline(x2, y2), col='grey', lwd=2)
#lines(density(locs, bw=4, kernel='rectangular'))

locs <- unlist(res.3[['pelt.loc']])
x3 <- sort(unique(locs))
y3 <- as.numeric(table(locs) / denom[min(locs):max(locs)])
lines(x=x3, y=y3, type='s', col=rgb(1,0,0,0.25))
lines(smooth.spline(x3, y3), col='firebrick3', lwd=2)
#lines(density(locs, bw=2))

# windows go from HXB2 790 to 9465
draw.hiv(x0=0, x1=82 * (9719-790)/(9465-790), y0=0, y1=0.025, label=T, cex.lab=0.8, col='grey')

dev.off()



plot(y1, y2[match(x1, x2)]) 
cor.test(y1, y2[match(x1, x2)], method='spearman') 

plot(y1, y3[match(x1, x3)])
cor.test(y1, y3[match(x1, x3)], method='spearman')

plot(y2, y3[match(x2, x3)])
cor.test(y2, y3[match(x2, x3)], method='spearman')


#####################################################################

# Hierachical Clustering
df1 <- read.csv(file = 'original_data/gene_dynamics_updated.tsv',
                head = T, sep=",", check.names=F, row.names=1,
                na.strings='0')

hc <- hclust(dist(df1, method="euclidean"))
labels <- gsub("^.+_(.+)$", "\\1", hc$labels)
plot(hc, labels=labels, cex=0.5)

#####################################################
# PCA
d <- dist(df1, method="euclidean")
pca <- princomp(d)
#biplot(pca, var.axes=F)
subtypes <- gsub("^.+_(.+)$", "\\1", row.names(pca$scores))

par(mar=c(5,5,1,1))
plot(pca$scores[,1], pca$scores[,2], type='n', xlab = "PC1", ylab = "PC2")
text(pca$scores[,1], pca$scores[,2],
     label=subtypes, cex=0.5)

subC <- which(subtypes=='C')
hpts <- chull(pca$scores[subC,1], pca$scores[subC,2])
hpts <- subC[c(hpts, hpts[1])]
polygon(pca$scores[hpts,1], pca$scores[hpts,2], border=NA, 
        col=rgb(0,1,0,0.1))

subB <- which(subtypes=='B')
hpts <- chull(pca$scores[subB,1], pca$scores[subB,2])
hpts <- subC[c(hpts, hpts[1])]
polygon(pca$scores[hpts,1], pca$scores[hpts,2], border=NA, 
        col=rgb(1,0,0,0.1))

subA <- which(subtypes=='A1')
hpts <- chull(pca$scores[subA,1], pca$scores[subA,2])
hpts <- subC[c(hpts, hpts[1])]
polygon(pca$scores[hpts,1], pca$scores[hpts,2], border=NA, 
        col=rgb(0,0,1,0.1))


##################################################
# Transition Plot
require(jsonlite)
library(RColorBrewer)

ds <- fromJSON('original_data/dynsbm_og/25_dynsbm.json')

pdf(file='~/git/skycross/paper/transition.pdf', width=5, height=5)
heatmap(ds$trans, cexRow = 0.8, cexCol = 0.8, #symm=T, scale="none",
        col=colorRampPalette(brewer.pal(8,'Blues'))(30))
dev.off()
#mtext(text = "Community", side = 1, line = 1.5, cex = 1.2)
#mtext(text = "Clusters", side = 2, line = 2.2, cex = 1.2)



# barplots.pdf
pal <- c("#6867C6", "#80A390", "#8AFBA3", "#DBCEE5", "#69E4BC",
         "#D087CF", "#CFC386", "#E5FB86", "#C0B3EC", "#BAE5F7", 
         "#DB4FD7", "#D1DBAF", "#FDFA60", "#DF68B6", "#89F2F2",
         "#ECB75B", "#E5D8CE", "#0885BA", "#AB68FD", "#B9EEA7",
         "#F07383", "#87F164", "#F4A588", "#6EE1FE", "#C78DAA")

pdf(file='~/git/skycross/paper/distribution.pdf', width=5, height=5)
par(mar=c(5,5,1,2))
barplot(apply(ds$membership, 2, table), col=c('white', pal[idx]), border=NA, 
        ylab = "", cex.axis = 1., names.arg=1:ncol(foo), cex.names=0.8,
        xlab="Genome window", las=1)
legend(x=100, y=500, legend=25:1, fill=rev(pal[idx]), bty='n', xpd=NA, 
       cex=0.6, border='white', x.intersp=0.1)
#barplot(foo[order(hm$rowInd), ], col=rainbow(25,s=0.5,v=.9)[hm$rowInd], border=NA)
#add text 
#mtext(text = "Genome window", side = 1, line = 2.5, cex = 1.1)
mtext(text = "Cluster memberships", side = 2, line = 3, cex = 1.1)
dev.off()




# Are there significant trends?
foo <- apply(ds$membership, 2, table)

plot(NA, xlim=c(1,ncol(foo)), ylim=c(0, 55))
for (i in c(24)) { lines(foo[i,], col=pal[i]) }  # cluttered

require(MASS)
fits <- lapply (1:nrow(foo), function(i) {
  temp <- data.frame(y=foo[i,], x=1:ncol(foo))
  glm.nb(y~x, data=temp)
})


# use centered log-ratio transform

# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
geo.mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

clr <- apply(foo, 1, function(x) log(x)-log(geo.mean(x)))
dclr <- apply(clr, 2, diff)
corrplot(cor(clr), diag=F, type='upper', sig.level=0.05/26, 
         tl.pos='n', method='ellipse', cl.pos='n')
res1 <- cor.mtest(clr, conf.level = .99)
summary(res1$p[upper.tri(res1$p)])
hist(log10(res1$p[upper.tri(res1$p)]))
alpha <- 0.05/(26*25/2)  # 0.0001538462

# 8-23 (clusters 7 and 22) and 8-24 are significant
cor.test(clr[,8], clr[,23])
plot(clr[,8], clr[,24])

fits <- lapply(1:ncol(clr), function(i) {
  temp <- data.frame(y=clr[,i], x=1:nrow(clr))
  glm(y~x, data=temp)
})




# subtypes.pdf - use `traces`
pdf(file='~/git/skycross/paper/subtype.pdf', width=6, height=8)
par(mfrow=c(3,1), mar=c(3,5,2,1))
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

# A1 = 4, 13, 24
for (s in c('A1', 'B', 'C')) {
  plot(NA, xlim=c(1,82), ylim=c(0,25), yaxt='n', ylab='Cluster', 
       xlab='Window', cex.lab=1.2)
  #  text(x=1.3, y=24, label=s, cex=2)
  title(main=paste('Subtype', s), adj=0, font=2, cex.main=1.5, line=0.5, xpd=NA)
  
  for (i in which(gene$subtype==s)) {
    lines(x=1:82, y=jitter(traces[i,]), 
          col=rgb(0,0,0,0.1), lwd=2)
  }
  axis(side=2, at=1:25, label=order(idx), las=2, cex.axis=0.6)
}
dev.off()


###############################################################################################

library(dplyr)
library("ggpubr")

# get frequency of breakpoints for each sequence
df1 <- read.csv(file = 'original_data/gene_dynamics.tsv',
                head = T, sep="\t", check.names=F, row.names=1,
                na.strings='0')

# breakpoint category data
bpts <- read.csv("original_data/breakpoints.csv")
rownames(bpts) <- rownames(df1)

# link frequency with sequence subtype classification
sub <- read.csv("original_data/cahr_updated.csv", sep = ',', header = T, quote = "", na.strings = T)
nam <- sub$accno[match(rownames(bpts), sub$accno)]
subc <- sub$class[match(rownames(bpts), sub$accno)]
subt <- sub$subtype[match(rownames(bpts), sub$accno)]
ndf <- data.frame(subc, bpts$nbpt3,  bpts$nbpt5,  bpts$nbpt10)

# store in pdf
pdf(file='original_data/pureVScrf_new.pdf', width=5, height=5)
par(mar=c(6,5,1,1), mfrow=c(1,1), xpd=F)

# color
pal <- brewer.pal(3, 'Set2')
df.m <- melt(ndf, id.var = "subc")

# draw boxplot
boxplot(data=df.m, value ~ subc + variable, ylab="Number of Breakpoints",
        xlab="", cex.lab = 1.2, boxwex=0.4,
        names = rep(c('Pure', 'CRF'), 3), las = 2,
        cex.axis=0.8, col = pal)

# add additional information
par(xpd=T)

segments(x0=0.9, x1=2.2, y0=-4.4, lwd=2)
text(x=1.5, y=-5, adj=0.5, label="1 breakpoint", cex=0.8)
segments(x0=2.8, x1=4.2, y0=-4.4, lwd=2)
text(x=3.5, y=-5, adj=0.5, label="2 breakpoints", cex=0.8)
segments(x0=4.8, x1=6.2, y0=-4.4, lwd=2)
text(x=5.5, y=-5, adj=0.5, label="3 breakpoints", cex=0.8)

par(xpd=F)
dev.off()

# wilcox test
test <- wilcox.test(ndf$bpts~factor(ndf$subc), data=ndf)

# median bpts
pure <- median(na.omit(ndf$bpts[ndf$subc == "Pure"]))
crf <- median(na.omit(ndf$bpts[ndf$subc == "Recombinant"]))
raing <- range(bpts)

################################################################################

# analyse pure and recobinant bps in the context of years
year <- sub$year[match(dt$names.freq., sub$accno)]
dat <- data.frame(nam, year, bpts)

ggscatter(dat, x = "year", y = "bpts", ylim = c(0, 17), xlim = c(1980, 2020),
          add = "reg.line", conf.int = TRUE, conf.int.level = 0.95,
          font.label = 18, cor.coef = TRUE, cor.method = "spearman", size = 1.2,
          xlab = "Year", ylab = "Number of Breakpoints") +  font("xlab", size = 14) +
  font("ylab", size = 14) + font("xy.text", size = 14)


#cor.test(ndf$bpts, dat$year, method="pearson")
#dev.copy2pdf2(file = paste0('time'), height = 5, width = 5)


###################

# are there any genomes where all windows are assigned to the same cluster?
baz <- apply(df, 1, function(x) {
  length(unique(x[x>0]))
})

# Can we define new references?
foo <- apply(df, 1, function(x) {
  tab <- table(x)
  c(max(tab), names(tab)[which.max(tab)])
})
foof <- as.data.frame(t(foo))
names(foof) <- c('max.count', 'max.type')
foof$max.count <- as.integer(foof$max.count)
foof$max.type <- as.integer(foof$max.type)

foof$names <- gene$pos


# what is the highest max.count such that all clusters are represented?
x <- seq(min(foof$max.count), max(foof$max.count))
y <- sapply(x, function(xx) length(unique(foof$max.type[foof$max.count >= xx])))
plot(x, y, type='s')

require(changepoint)
breakpoints <- function(traces, max.missing=10, min.seg.len=10) {
  amoc <- rep(NA, nrow(traces))
  pelt <- rep(NA, nrow(traces))
  pelt.loc <- list()
  for (i in 1:nrow(traces)) {
    tr <- traces[i,]
    if (sum(is.na(tr)) > max.missing) {
      next
    }
    tra <- tr[!is.na(tr)]  # simply drop missing values
    res <- cpt.mean(tra, method='PELT', minseglen = min.seg.len)
    pelt[i] <- length(res@cpts)-1  # number of breakpoints
    # map back to original coords (if any NAs)
    loc <- res@cpts[-length(res@cpts)]
    na.count <- cumsum(is.na(tr))
    pelt.loc[[i]] <- loc + na.count[loc]
    res <- cpt.mean(tra, minseglen=min.seg.len)  # default method, AMOC
    amoc[i] <- length(res@cpts)-1
  }
  list('amoc'=amoc, 'pelt'=pelt, 'pelt.loc'=pelt.loc)
}


