#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`
NUMCORES=`nproc --all`

mkdir -p COHORT_QC/testdata

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/cohort_qc COHORT_QC/testdata --recursive

ls -l COHORT_QC/testdata
echo $ONCOKB_KEY


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata \
  --rm \
  $1/single_cell_pipeline_qc:$TAG \
  single_cell cohort_qc --input_yaml single_cell/tests/jenkins/cohort_qc/inputs.yaml \
  --maxjobs $NUMCORES --nocleanup --sentinel_only  \
  --submit local --loglevel DEBUG \
  --tmpdir COHORT_QC/temp \
  --pipelinedir COHORT_QC/pipeline \
  --submit local \
  --out_dir COHORT_QC/output \
  --config_override '{"refdir":"refdata"}'
  --api_key $ONCOKB_KEY

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline_qc:$TAG rm -rf COHORT_QC
