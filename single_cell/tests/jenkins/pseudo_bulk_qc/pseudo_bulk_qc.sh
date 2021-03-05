#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`
NUMCORES=`nproc --all`

mkdir -p PSEUDO_BULK_QC/ref_test_data

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/pseudo_bulk_qc PSEUDO_BULK_QC/ref_test_data --recursive

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline_qc:$TAG \
  single_cell pseudo_bulk_qc --input_yaml single_cell/tests/jenkins/pseudo_bulk_qc/inputs.yaml \
  --maxjobs $NUMCORES --nocleanup --sentinel_only  \
  --submit local --loglevel DEBUG \
  --tmpdir PSEUDO_BULK_QC/temp \
  --pipelinedir PSEUDO_BULK_QC/pipeline \
  --submit local \
  --out_dir PSEUDO_BULK_QC/output \
  --config_override '{"annotation": {"chromosomes": ["6", "8", "17"]}, "version": '\"$TAG\"'}' \

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline:$TAG rm -rf PSEUDO_BULK_QC
