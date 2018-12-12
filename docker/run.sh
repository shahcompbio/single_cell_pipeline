#!/bin/bash

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

docker build -t single_cell_pipeline -f docker/dockerfile .

docker tag single_cell_pipeline shahlab.azurecr.io/scp/single_cell_pipeline:$TAG

docker push shahlab.azurecr.io/scp/single_cell_pipeline:$TAG
