library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg19)

args<-commandArgs(TRUE)

bamdir <- args[1]
output <- args[2]

Aneufinder(inputfolder = bamdir, outputfolder = output, numCPU = 1, method = c("dnacopy"), binsizes = 500 * 1000, chromosomes = c(1:22,'X','Y'), correction.method='GC', GC.BSgenome = BSgenome.Hsapiens.UCSC.hg19, cluster.plots=FALSE)