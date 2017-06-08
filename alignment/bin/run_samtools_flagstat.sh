# Args:
#   [1] : BAM_FILE - Path to BAM file
#   [2] : OUT_FILE - Path where output will be written

BAM_FILE=$1

OUT_FILE=$2

SAMTOOLS=$3

$SAMTOOLS flagstat $BAM_FILE > $OUT_FILE || { echo 'samtools flagstat command failed'; exit 1; }
