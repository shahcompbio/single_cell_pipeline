# Test data download and setup on and off the cluster

#### Requirements:

1. installation of wget or a similar tool

### Download the test data

All test data is stored on Azure in public storage containers. To download test
data sets for different parts of the pipeline, simply use wget + the appropriate URL
for the test data of your choice. 

## Single Cell Pipeline test data

All test data and example input yamls are available on the public container 
```
https://singlecelltestsets.blob.core.windows.net/public-testdata
```
in tarball form. See below for sub-pipeline-specific download paths.

*alignment:*  
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/alignment.tar.gz 
```
*annotation:*
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/annotation.tar.gz 
```
*hmmcopy:*
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/hmmcopy.tar.gz 
```
*variant calling:*
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/variant-calling.tar.gz 
```
*breakpoint calling:*
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/breakpoint-calling.tar.gz 
```
*infer haplotypes:*
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/infer-haps.tar.gz 
```
*count haplotypes:*
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/count-haps.tar.gz 
```
*merge bams:*
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/merge-bams.tar.gz 
```
*split bams:*
``` 
https://singlecelltestsets.blob.core.windows.net/public-testdata/split-bams.tar.gz 
```

Each tarball contains test data and a premade input.yaml that can be fed into the pipeline.  


