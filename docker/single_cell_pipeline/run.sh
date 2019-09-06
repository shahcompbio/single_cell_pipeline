#!/bin/bash

REGISTRY=$1

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

docker build -t single_cell_pipeline -f docker/single_cell_pipeline/dockerfile .

docker tag single_cell_pipeline $REGISTRY/singlecellpipeline/single_cell_pipeline:$TAG

docker push $REGISTRY/singlecellpipeline/single_cell_pipeline:$TAG
