#!/bin/bash

if [ ! -d "/refdata" ]; then
    docker run -v /refdata:/refdata singlecellpipeline/azurecli:v0.0.1 az storage blob download-batch -s refdata -d /refdata --account-name $1 --account-key $2
fi
