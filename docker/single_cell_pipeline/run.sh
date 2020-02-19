#!/bin/bash

REGISTRY=$1
ORG=$2
TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

cat docker/single_cell_pipeline/dockerfile_template \
 | sed "s/{git_commit}/$TAG/g" \
 > docker/single_cell_pipeline/dockerfile

docker build -t $REGISTRY/$ORG/single_cell_pipeline:$TAG -f docker/single_cell_pipeline/dockerfile . --no-cache

docker push $REGISTRY/$ORG/single_cell_pipeline:$TAG
