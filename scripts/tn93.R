library(seqinr)
library(ape)
#library(strataG)
library(parallel)

setwd('')

files <- list.files(paste0("mafft/"))

mclapply(files, function(f) {
  my_align <- read.alignment(file = paste0("mafft/", f), format='fasta')
  y <- as.DNAbin(my_align)
  tn93 <- dist.dna(y, model = "TN93", variance = FALSE,
                   gamma = FALSE, pairwise.deletion = FALSE,
                   base.freq = NULL, as.matrix = FALSE)
  
  write.table(data.frame(as.matrix(tn93)), 
              file = paste0("tn93/", f, ".tn93.csv"), 
              append = F, 
              quote = FALSE, 
              sep = ",", 
              col.names = NA)
},
mc.cores = 8
)
