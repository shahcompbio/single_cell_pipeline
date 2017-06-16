# Args:
#   [1] : READ_COUNTER - Path to readCounter executable shipped with HMMCopy
#   [2] : BAM_FILE - Path to BAM file to analyse
#   [3] : OUT_FILE - Path where results will be written
#   [4] : WINDOW - Size of non-overlapping windows
#   [5] : MIN_MQUAL - Minimum mapping quality for an alignment to be counted
#   [6] : CHROMS - Comma separated list of chromosomes to build WIG for.

BAM_FILE=$1

OUT_FILE=$2

WINDOW=$3

MIN_MQUAL=$4

CHROMS=$5

readCounter -w $WINDOW -q $MIN_MQUAL -c $CHROMS $BAM_FILE > $OUT_FILE
