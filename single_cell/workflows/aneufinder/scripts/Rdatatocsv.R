library(GenomicRanges)

args <- commandArgs(TRUE)
segments_filename <- args[1]
segments_outfile <- args[2]
load(segments_filename)
write.csv(model$segments, file=segments_outfile)

reads_outfile <- args[3]
write.csv(model$bins, file=reads_outfile)