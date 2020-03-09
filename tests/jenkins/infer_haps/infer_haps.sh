#!/bin/bash
set -e
set -o pipefail

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`


mkdir -p INFER_HAPS/bam

docker run -v $PWD:$PWD -w $PWD singlecellpipeline/azurecli:v0.0.1 \
  az storage blob download-batch -s infer-haps  -d INFER_HAPS/bam/ --account-name $1 --account-key $2


docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  singlecellpipeline/single_cell_pipeline:$TAG \
  single_cell infer_haps --input_yaml tests/jenkins/infer_haps/inputs.yaml \
  --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config tests/jenkins/context_config.yaml \
  --submit local --loglevel DEBUG \
  --tmpdir INFER_HAPS/temp \
  --pipelinedir INFER_HAPS/pipeline --submit local --out_dir INFER_HAPS/output \
  --config_override '{"infer_haps":{"chromosomes":["15"], "extract_seqdata": {"genome_fai_template": "/refdata/human/infer_haps/GRCh37-lite.fa.fai", "genome_fasta_template": "/refdata/human/infer_haps/GRCh37-lite.fa"}, "ref_data_dir": "/refdata/human/infer_haps/"}}' \
   --rerun
