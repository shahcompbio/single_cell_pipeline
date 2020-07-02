#!/bin/bash

if [ ! -d "/refdata" ]; then
    docker run -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_DEFAULT_REGION -v /refdata:/refdata singlecellpipeline/awscli:v0.0.1 aws s3 cp s3://singlecelltestsets/TESTDATA_CODEBUILD/refdata /refdata --recursive
fi

