# Args:
#   [1] : FASTQ - FASTQ file
#   [2] : OUT_HTML - File where .html report will be written
#   [3] : OUT_ZIP - File where .zip report files will be written

BAM=$1

METRICS=$2

HIST=$3

picard CollectInsertSizeMetrics INPUT=$BAM OUTPUT=$METRICS HISTOGRAM_FILE=$HIST ASSUME_SORTED=True VALIDATION_STRINGENCY=LENIENT

