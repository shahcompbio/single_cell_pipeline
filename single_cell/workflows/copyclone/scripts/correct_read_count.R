library(HMMcopy)

args <- commandArgs(TRUE)

infile <- args[1]
gc <- args[2]
map <- args[3]
outfile <- args[4]
map_cutoff <- as.numeric(args[5])
sample <- args[6]

chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

format_read_count_table <- function(samp.corrected, chromosomes) {
        hmmcopy.table <- as.data.frame(samp.corrected)
        colnames(hmmcopy.table)[1] <- "chr"
        hmmcopy.table$chr <- factor(hmmcopy.table$chr, levels=chromosomes, ordered=T)
        hmmcopy.table <- hmmcopy.table[order(hmmcopy.table$chr),]
}

# read alignment data
samp.uncorrected <- wigsToRangedData(infile, gc, map, verbose=F)

# correct and segment data
samp.corrected <- try(correctReadcount(samp.uncorrected, verbose=F), silent=T)
if (inherits(samp.corrected, "try-error") || length((which(samp.corrected$cor.map == Inf))) > 0) {
        warning(paste("Low coverage sample results in loess regression failure, unable to correct and segment ", infile, sep=""))

        uncorrected.table <- format_read_count_table(samp.uncorrected, chromosomes)
        uncorrected.table["cor_gc"] <- NA
        uncorrected.table["cor_map"] <- NA
        uncorrected.table["ideal"] <- NA
        uncorrected.table["valid"] <- NA
        uncorrected.table["state"] <- NA
        uncorrected.table["copy"] <- NA
        uncorrected.table["integer_copy_number"] <- NA
        uncorrected.table["integer_copy_scale"] <- NA
        uncorrected.table["modal_curve"] <- NA
        uncorrected.table["modal_quantile"] <- NA
        uncorrected.table["cell_id"] <- sample

        write.table(format(uncorrected.table, scientific=F, trim=T), file=outfile, quote=F, sep=",", col.names=T, row.names=F)

} else{

samp.corrected <- format_read_count_table(samp.corrected)
colnames(samp.corrected) <- c("chr","start","end","width","reads","gc","map","valid","ideal","cor_gc","cor_map","copy")
samp.corrected["modal_curve"] <- NA
samp.corrected["modal_quantile"] <- NA
samp.corrected["cell_id"] <- sample

samp.corrected["state"] <- NA
samp.corrected["integer_copy_number"] <- NA
samp.corrected["integer_copy_scale"] <- NA

samp.corrected["copy"] <- samp.corrected["cor_gc"]

# if mappability cutoff given, remove bins with mappability below cutoff value
if (!is.null(map_cutoff)) {
        samp.corrected$copy[samp.corrected$map < map_cutoff] <- NA
}


write.table(samp.corrected, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep =",")

}
