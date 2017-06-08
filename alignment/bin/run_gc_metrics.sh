# Args:
#   [1] : FASTQ - FASTQ file
#   [2] : OUT_HTML - File where .html report will be written
#   [3] : OUT_ZIP - File where .zip report files will be written

JAVA=$1

PICARD=$2

BAM=$3

METRICS=$4

REFERENCE=$5

SUMM=$6

CHART=$7

R=$8

export PATH=$(dirname $R):$PATH



$JAVA -Xmx4g -jar $PICARD CollectGcBiasMetrics INPUT=$BAM OUTPUT=$METRICS REFERENCE_SEQUENCE=$REFERENCE S=$SUMM CHART_OUTPUT=$CHART VALIDATION_STRINGENCY=LENIENT

