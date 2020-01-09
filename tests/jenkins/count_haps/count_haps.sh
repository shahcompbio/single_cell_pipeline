#!/bin/bash

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`


mkdir -p COUNT_HAPS/data

docker run -v $PWD:$PWD -w $PWD singlecellpipeline/azurecli:v0.0.1 \
  az storage blob download-batch -s count-haps  -d COUNT_HAPS/data/ --account-name $1 --account-key $2


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  singlecellpipeline/single_cell_pipeline:$TAG \
  single_cell count_haps --input_yaml tests/jenkins/count_haps/inputs.yaml \
  --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config tests/jenkins/count_haps/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir COUNT_HAPS/temp \
  --config_override '{"count_haps":{"chromosomes":["15"], "extract_seqdata": {"genome_fai_template": "/refdata/human/infer_haps/GRCh37-lite.fa.fai", "genome_fasta_template": "/refdata/human/infer_haps/GRCh37-lite.fa"}, "ref_data_dir": "/refdata/human/infer_haps/"}}' \
  --pipelinedir COUNT_HAPS/pipeline --submit local --out_dir COUNT_HAPS/output
