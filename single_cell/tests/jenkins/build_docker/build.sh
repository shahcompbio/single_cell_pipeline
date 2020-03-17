# lets assume that the commit included the tag for the latest version
# build the container, tag with version
# push it to azure container registry
# usage: bash build_docker.sh dockerhub_username dockerhub_password
#!/usr/bin/env bash
echo "\n LOGIN \n"
docker login -u $1 --password $2

REGISTRY=$3
ORG=$4
TAG=`git describe --tags $(git rev-list --tags --max-count=1)`
TAG="${TAG}.beta"

COMMIT=`git rev-parse HEAD`

cat tests/jenkins/build_docker/dockerfile_template \
 | sed "s/{git_commit}/$COMMIT/g" \
 > tests/jenkins/build_docker/dockerfile

docker build -t $REGISTRY/$ORG/single_cell_pipeline:$TAG -f tests/jenkins/build_docker/dockerfile . --no-cache

docker push $REGISTRY/$ORG/single_cell_pipeline:$TAG

