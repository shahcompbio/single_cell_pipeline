# Args:
#   [1] : REF_GENOME - Path to reference genome BAM file was aligned to
#   [2] : FASTQ_1 - First FASTQ with paired end data
#   [3] : FASTQ_2 - Second FASTQ with paired end data
#   [4] : OUT_FILE - Path where output will be written in BAM format
#   [5] : READ_GROUP - Optional tab-separated read group information line to add to the BAM header
#         "@RG\tID:someid\tPL:illumina\tPU:someplatformunit\tLB:somelibraryid\tSM:somesample\tCN:somecentreid"

REF_GENOME=$1

FASTQ_1=$2

FASTQ_2=$3

OUT_FILE=$4

OUT_SAI_1=$(echo "$OUT_FILE" | sed 's/\.bam/_R1.sai/')

OUT_SAI_2=$(echo "$OUT_FILE" | sed 's/\.bam/_R2.sai/')

# align with BWA

bwa aln $REF_GENOME $FASTQ_1 > $OUT_SAI_1 || { echo 'bwa aln command failed for read 1'; exit 1; }

bwa aln $REF_GENOME $FASTQ_2 > $OUT_SAI_2 || { echo 'bwa aln command failed for read 2'; exit 1; }

bwa sampe ${5:+-r "${5}"} $REF_GENOME $OUT_SAI_1 $OUT_SAI_2 $FASTQ_1 $FASTQ_2 | samtools view -bSh - > $OUT_FILE || { echo 'bwa sampe command failed'; exit 1; }

# remove temporary files

if [ -f $OUT_SAI_1 ]; then
	rm $OUT_SAI_1
fi

if [ -f $OUT_SAI_2 ]; then
	rm $OUT_SAI_2
fi
