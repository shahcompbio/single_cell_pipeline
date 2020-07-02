#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
DOCKER=`which docker`

mkdir -p INFER_HAPS/ref_test_data

docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v $PWD:$PWD -w $PWD $1/awscli:v0.0.1 \
  aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/infer-haps INFER_HAPS/ref_test_data/ --recursive

docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline:$TAG \
  single_cell infer_haps --input_yaml single_cell/tests/jenkins/infer_haps/inputs.yaml \
  --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config single_cell/tests/jenkins/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir INFER_HAPS/temp \
  --pipelinedir INFER_HAPS/pipeline \
  --submit local \
  --out_dir INFER_HAPS/output \
  --config_override '{"infer_haps":{"chromosomes":["15"], "extract_seqdata": {"genome_fai_template": "/refdata/human/infer_haps/GRCh37-lite.fa.fai", "genome_fasta_template": "/refdata/human/infer_haps/GRCh37-lite.fa"}, "ref_data_dir": "/refdata/human/infer_haps/"}, "version": '\"$TAG\"'}' \


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v $DOCKER:$DOCKER --rm \
  $1/single_cell_pipeline:$TAG \
  python single_cell/tests/jenkins/infer_haps/test_infer_haps.py INFER_HAPS/output INFER_HAPS/ref_test_data

docker run -w $PWD -v $PWD:$PWD --rm $1/single_cell_pipeline:$TAG rm -rf INFER_HAPS
