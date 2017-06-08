# Args:
#   [1] : FASTQ - FASTQ file
#   [2] : OUT_HTML - File where .html report will be written
#   [3] : OUT_ZIP - File where .zip report files will be written

JAVA=$1

PICARD=$2

BAM=$3

METRICS=$4

HIST=$5

R=$6

export PATH=$(dirname $R):$PATH



$JAVA -Xmx4g -jar $PICARD CollectInsertSizeMetrics INPUT=$BAM OUTPUT=$METRICS HISTOGRAM_FILE=$HIST ASSUME_SORTED=True VALIDATION_STRINGENCY=LENIENT

