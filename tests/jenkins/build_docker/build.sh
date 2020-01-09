# lets assume that the commit included the tag for the latest version
# build the container, tag with version
# push it to azure container registry
# usage: bash build_docker.sh dockerhub_username dockerhub_password azure_client_id azure_secret_key azure_registry_url
#!/usr/bin/env bash
echo "\n LOGIN \n"
docker login scdnastagingcr.azurecr.io -u $3 --password $4
docker login -u $1 --password $2

echo "\n RUN DOCKER \n"
sh docker/single_cell_pipeline/run.sh $5
sh docker/single_cell_pipeline/run.sh docker.io

sh docker/scgenome/run.sh $5
sh docker/scgenome/run.sh docker.io
