# lets assume that the commit included the tag for the latest version
# build the container, tag with version
# push it to azure container registry
# usage: bash build_docker.sh dockerhub_username dockerhub_password
#!/usr/bin/env bash
echo "\n LOGIN \n"
docker login -u $1 --password $2

echo "\n RUN DOCKER \n"
sh docker/single_cell_pipeline/run.sh docker.io singlecellpipelinetest

sh docker/scgenome/run.sh docker.io singlecellpipelinetest
