# Transition Plot
require(jsonlite)
library(RColorBrewer)
ds <- fromJSON("bestoutput.json")
pdf(file='transition.pdf', width=5, height=5)
heatmap(ds$trans, cexRow = 0.8, cexCol = 0.8, #symm=T, scale="none",
        col=colorRampPalette(brewer.pal(8,'Blues'))(30))
dev.off()