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

COMMIT=`git rev-parse HEAD`

cat single_cell/tests/jenkins/build_docker_staging/dockerfile_template \
 | sed "s/{git_commit}/$COMMIT/g" \
 > single_cell/tests/jenkins/build_docker_staging/dockerfile

docker build -t $REGISTRY/$ORG/single_cell_pipeline:$TAG -f single_cell/tests/jenkins/build_docker_staging/dockerfile . --no-cache

docker push $REGISTRY/$ORG/single_cell_pipeline:$TAG

