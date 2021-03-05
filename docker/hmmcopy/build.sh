REGISTRY=$1
ORG=$2

echo "\n LOGIN \n"
docker login $REGISTRY -u $3 --password $4

TAG=`git describe --tags $(git rev-list --tags --max-count=1)`

COMMIT=`git rev-parse HEAD`

cat dockerfile_template \
 | sed "s/{git_commit}/$COMMIT/g" \
 > dockerfile

docker build -t $REGISTRY/$ORG/single_cell_pipeline_hmmcopy:$TAG . --no-cache

docker push $REGISTRY/$ORG/single_cell_pipeline_hmmcopy:$TAG

