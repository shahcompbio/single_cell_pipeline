#!/bin/bash

INPUTVCF=$1
OUTPUTMAF=$2
FASTA=$3
VEPDATA=$4
BUFFERSIZE=$5

vcf2maf.pl --input-vcf $1 --output-maf $2 --ref-fasta $3 --vep-data $4  --vep-path $(dirname `which vep`) --buffer-size $5
