# Args:
#   [1] : FASTQ - FASTQ file
#   [2] : OUT_HTML - File where .html report will be written
#   [3] : OUT_ZIP - File where .zip report files will be written

BAM=$1

METRICS=$2

REFERENCE=$3

SUMM=$4

CHART=$5

picard CollectGcBiasMetrics INPUT=$BAM OUTPUT=$METRICS REFERENCE_SEQUENCE=$REFERENCE S=$SUMM CHART_OUTPUT=$CHART VALIDATION_STRINGENCY=LENIENT

