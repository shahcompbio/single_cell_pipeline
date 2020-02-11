# Single Cell Pipeline

Welcome to the home page for the single cell pipelines documentation.

## Quick Setup
The single cell pipeline is composed of 10 subpipelines that can be run individually. Here you will find everything you need to run any of the subpipelines with small testdata and expected outputs.

### Reference Data
Before you can run a subpipeline you must acquire the necessary reference data. If you are working from the `juno` cluster, 
reference data can be found at `/juno/work/shah/reference/singlecellpipeline`. If you are not working from juno, 
you must download the data locally from microsoft azure. 

Run the following command to download reference data from the account to a local path, replacing `{YOUR_REF_DATA_PATH}` with your preferred destination.
```
docker run -v /refdata:/refdata singlecellpipeline/azurecli:v0.0.1 az storage blob download-batch -s refdata -d {YOUR_REF_DATA_PATH} --account-name singlecelltestsets --account-key {$ACCOUNT KEY}
```
Ask us for the storage account password and replace `{$ACCOUNT KEY}` with it.


### Test Data
Test data for each of the subpipelines can also be downloaded from azure storage. 
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
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/merge_bams.tar.gz && tar -xvf merge-bams.tar.gz 
rm merge-bams.tar.gz && cd merge-bams
```
#### Split Bam
Description: Split bulk wgs bam into region-specific bams for downstream analysis
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/split_bam.tar.gz && tar -xvf split-bam.tar.gz 
rm split-bam.tar.gz && cd split-bam
```
#### Variant Calling
Description: Generate consensus variant calls from strelka and museq, perform vcf annotation
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/variant_calling.tar.gz && tar -xvf variant-calling.tar.gz 
rm variant-calling.tar.gz && cd variant-calling
```
#### Breakpoint Calling
Description: Generate consensus breakpoint calls using destruct and lumpy
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/breakpoint_calling.tar.gz && tar -xvf breakpoint-calling.tar.gz 
rm breakpoint-calling.tar.gz && cd breakpoint-calling
```
#### Infer Haps
Description: Infer haplotypes.
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/infer_haps.tar.gz && tar -xvf infer-haps.tar.gz 
rm infer-haps.tar.gz && cd infer-haps
```
#### Count Haps
Description: Call haplotypes.
```
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/count_haps.tar.gz && tar -xvf count-haps.tar.gz 
rm count-haps.tar.gz && cd count-haps
```

#### SV Genotyping
```buildoutcfg
wget https://singlecelltestsets.blob.core.windows.net/public-testdata/sv_genotyping.tar && tar -xvf sv_genotyping.tar 
rm sv_genotyping.tar && cd sv_genotyping
```

You should now be in a named directory with everything you need to run the pipeline.

#### Customizing with your chosen paths

Before you can run the pipeline, you need to specify your chosen ref data path in the input files. This happens in two files in particular which are needed to run all subpipelines: in `context_config.yaml` and in `final_run.sh`. In both, replace `REF_DATA` with your chosen path.
Finally, replace `WORKING_DIR` with your working directory in `context_config.yaml`

#### Executing a pipeline

Simply run 
```
sh final_run.sh
```
To launch the pipeline.

#### Verifying your setup

Each downloadable `tar` files contains expected output from each pipeline in a directory labeled `expected_outputs`. You can use this to verify your local setup by comparing it with the result of your current pipeline run. Results from the current run can be found in a directory labeled `output`. 


