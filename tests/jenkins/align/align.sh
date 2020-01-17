!/bin/bash

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

#make test data directory
mkdir -p ALIGN/ref_test_data

#pull test data and ref data
docker run -v $PWD:$PWD -w $PWD singlecellpipeline/azurecli:v0.0.1 \
  az storage blob download-batch -s alignment -d ALIGN/ref_test_data --account-name $1 --account-key $2

#run alignment
docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
  -v /usr/bin/docker:/usr/bin/docker --rm \
  singlecellpipeline/single_cell_pipeline:$TAG \
  single_cell alignment --input_yaml tests/jenkins/align/inputs.yaml \
  --library_id A97318A --maxjobs 4 --nocleanup --sentinel_only  \
  --context_config tests/jenkins/context_config.yaml \
  --submit local --loglevel DEBUG \
  --config_override '{"pypeliner_storage_account": "scdnadev"}' --tmpdir ALIGN/temp \
  --pipelinedir ALIGN/pipeline --submit local --out_dir ALIGN/output --bams_dir ALIGN/bams


# docker run -w $PWD -v $PWD:$PWD -v /refdata:/refdata -v /var/run/docker.sock:/var/run/docker.sock \
#   -v /usr/bin/docker:/usr/bin/docker --rm \
#   singlecellpipeline/single_cell_pipeline:$TAG \
python tests/jenkins/align/test_alignment.py ALIGN/output A97318A  ALIGN/ref_test_data/refdata

# rm -r ALIGN
