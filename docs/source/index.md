# Single Cell Pipeline


Welcome to the home page for the single cell pipelines documentation.


## Quick Setup
### Reference Data
Before you can run a subpipeline you must acquire the necessary reference data. If you are working from `juno`, 
reference data can be found at `/juno/work/shah/reference/singlecellpipeline`. If you are not working from juno, 
you may mount the path locally, or download the reference data from the singlecelltestsets azure storage account. 

Run the following command to download reference data from the account:
```
docker run -v /refdata:/refdata singlecellpipeline/azurecli:v0.0.1 az storage blob download-batch -s refdata -d /refdata --account-name singlecelltestsets --account-key {$ACCOUNT KEY}
```
Ask us for the storage account password and replace `{$ACCOUNT KEY}` with it.

By default, the above command downloads the reference data to a directory named `/refdata`. 
If you choose to change it, you will have to manually update the pipeline's internal configuration with your chosen directory, 
found in a file named `config.yaml`, which is generated automatically when the pipeline is run.

### Test Data
Test data for each of the subpipelines can be downloaded from azure storage. 
Just execute the below commands for the subpipeline you which to run.
#### Alignment
Description: Aligns cell bam files, performs QC
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/alignment.tar.gz && tar -xvf alignment.tar.gz 
rm alignment.tar.gz && cd alignment
```
#### Hmmcopy
Description: generate read count wig files, perform GC correction, predict copynumber states
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/hmmcopy.tar.gz && tar -xvf hmmcopy.tar.gz 
rm hmmcopy.tar.gz && cd hmmcopy
```
#### Annotation
Description: Assign cell-specific quality scores and cell states, consolidate alignment and hmmcopy metrics and qc data 
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/annotation.tar.gz && tar -xvf annotation.tar.gz 
rm annotation.tar.gz && cd annotation
```
#### Merge Bams
Description: Merge cell bams into region-specific bams for downstream analysis
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/merge-bams.tar.gz && tar -xvf merge-bams.tar.gz 
rm merge-bams.tar.gz && cd merge-bams
```
#### Split Bam
Description: Split bulk wgs bam into region-specific bams for downstream analysis
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/split-bam.tar.gz && tar -xvf split-bam.tar.gz 
rm split-bam.tar.gz && cd split-bam
```
#### Variant Calling
Description: Generate consensus variant calls from strelka and museq, perform vcf annotation
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/variant-calling.tar.gz && tar -xvf variant-calling.tar.gz 
rm variant-calling.tar.gz && cd variant-calling
```
#### Breakpoint Calling
Description: Generate consensus breakpoint calls using destruct and lumpy
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/breakpoint-calling.tar.gz && tar -xvf breakpoint-calling.tar.gz 
rm breakpoint-calling.tar.gz && cd breakpoint-calling
```
#### Infer Haps
Description: Infer haplotypes.
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/infer-haps.tar.gz && tar -xvf infer-haps.tar.gz 
rm infer-haps.tar.gz && cd infer-haps
```
#### Count Haps
Description: Call haplotypes.
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/count-haps.tar.gz && tar -xvf count-haps.tar.gz 
rm count-haps.tar.gz && cd count-haps
```

At this point, you should be in a new working directory with all the testdata you need to run the pipeline and the necessary reference data ready to go.
In addition to the test data, you will have two `yaml` files: `input.yaml`, and `context.yaml`. 
<br/><br/>`input.yaml`: specifies the testdata files you are using.
<br/><br/>`context.yaml`: specifies information used by the pipeline to run locally in docker.

#### Running the pipeline 
In addition to all the other files you will have had downloaded, there will be two `.sh` files.
The first `final_run.sh` runs the pipeline within the pipeline docker container. The latter, `run_pipeline.sh`,
specifies all the pipeline inputs.

Simply run 
```
sh final_run.sh
```
To launch the pipeline.
