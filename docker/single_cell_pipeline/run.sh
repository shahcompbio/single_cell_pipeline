#!/bin/bash

REGISTRY=$1
TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

cat docker/single_cell_pipeline/dockerfile_template \
 | sed "s/{container_registry}/$REGISTRY/g" \
 | sed "s/{git_commit}/jenkins/g" \
 > docker/single_cell_pipeline/dockerfile

docker build -t single_cell_pipeline -f docker/single_cell_pipeline/dockerfile . --no-cache

docker tag single_cell_pipeline $REGISTRY/singlecellpipeline/single_cell_pipeline:$TAG

docker push $REGISTRY/singlecellpipeline/single_cell_pipeline:$TAG
