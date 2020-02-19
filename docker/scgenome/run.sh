#!/bin/bash

REGISTRY=$1
ORG=$2
sed -i "s/{container_registry}/$1/g" docker/scgenome/dockerfile

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

docker build -t scgenome -f docker/scgenome/dockerfile .

docker tag scgenome $REGISTRY/$ORG/scgenome:$TAG

docker push $REGISTRY/$ORG/scgenome:$TAG
