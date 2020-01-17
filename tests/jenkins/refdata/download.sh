#!/bin/bash


if [ ! -d "/refdata" ]; then
    docker run -v /refdata:/refdata singlecellpipeline/azurecli:v0.0.1 az storage blob download-batch -s refdata -d /refdata --account-name $1 --account-key $2
fi

#better solution gets permission error
# docker run singlecellpipeline/azurecli:v0.0.1 az storage blob list -c refdata --account-name $1 --account-key $2 > files.json

# python test/jenkins/refdata/parse_filelist.py files.json > files.txt && rm files.json

# REFDIR="/refdata"
# nfiles_exist=0
# nfiles_total=0
# while IFS= read -r f; do
#     if [ -e "$REFDIR/$f" ]; then
#         $((nfiles_total=$nfiles_total+1))
#         $((nfiles_exists=$nfiles_exists+1))
#     else
#         $((nfiles_total=$nfiles_total+1))
#     fi
# done < "files.txt"


#  if [  $nfiles_exist -lt $nfiles_total ]; then
#      docker run -v /refdata:/refdata singlecellpipeline/azurecli:v0.0.1 az storage blob download-batch -s refdata -d /refdata --account-name $1 --account-key $2
#  fi
