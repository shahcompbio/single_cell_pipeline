#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

NUMCORES=`nproc --all`

mkdir -p ALIGN/ref_test_data

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/alignment ALIGN/ref_test_data --recursive --quiet


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata --rm \
  $1/single_cell_pipeline_alignment:$TAG \
  single_cell alignment --input_yaml single_cell/tests/codebuild/align/inputs.yaml \
  --library_id A97318A --maxjobs $NUMCORES --nocleanup --sentinel_only  \
  --submit local --loglevel DEBUG \
  --tmpdir ALIGN/temp \
  --pipelinedir ALIGN/pipeline \
  --submit local \
  --out_dir ALIGN/output \
  --bams_dir ALIGN/bams \
  --sequencing_center TEST

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata --rm \
  $1/single_cell_pipeline_alignment:$TAG \
  python single_cell/tests/codebuild/align/test_alignment.py ALIGN/output A97318A  ALIGN/ref_test_data/refdata/bwa-mem


docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline_alignment:$TAG rm -rf ALIGN
