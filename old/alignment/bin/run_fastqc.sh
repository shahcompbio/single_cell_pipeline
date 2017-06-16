# Args:
#   [1] : FASTQ - FASTQ file
#   [2] : OUT_HTML - File where .html report will be written
#   [3] : OUT_ZIP - File where .zip report files will be written

FASTQ=$1

OUT_HTML=$2

OUT_ZIP=$3

FASTQC=$4

JAVA=$5

export PATH=$(dirname $JAVA):$PATH
# run fastqc



OUT_DIR=$(dirname $OUT_HTML)

FASTQ_TMP=$(echo ${OUT_DIR}/$(basename $FASTQ) | sed 's/\.fastq\.gz/.fastq/')

echo $FASTQ

gunzip -c $FASTQ > $FASTQ_TMP || { echo 'gunzip command failed'; exit 1; }

$FASTQC --outdir=${OUT_DIR} $FASTQ_TMP || { echo 'fastqc command failed'; exit 1; }

if [ -f $FASTQ_TMP ]; then
	rm $FASTQ_TMP
fi

# change file names

HTML_FILE=$(echo ${OUT_DIR}/$(basename $FASTQ) | sed 's/\.fastq\.gz/_fastqc\.html/')
ZIP_FILE=$(echo ${OUT_DIR}/$(basename $FASTQ) | sed 's/\.fastq\.gz/_fastqc\.zip/')

if [ -f $HTML_FILE ]; then
	mv $HTML_FILE $OUT_HTML
fi

if [ -f $ZIP_FILE ]; then
	mv $ZIP_FILE $OUT_ZIP
fi
