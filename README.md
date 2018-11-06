# Single Cell Pipeline 

For a detailed guide see [INSTALL](INSTALL/README.md)

For azure documentation see [azure](azure/README.md)

## Setup and installation

Set up conda with the required packages.

First ensure you have the correct channels:

```
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels https://conda.anaconda.org/aroth85
conda config --add channels 'bioconda'
conda config --add channels 'r'
conda config --add channels 'conda-forge'

```

### From Source

Then create an environment with the required packages:

```
conda create --name singlecellpipeline --file conda_packages.txt
```

Activate the environment:

```
source activate singlecellpipeline
```

Add the single cell pipeline into the current site packages:

```
python setup.py develop
```

Use develop mode to allow for modifying of the code during the debugging phase.
Later versions will be released through conda.

### From anaconda.org

To install the latest production version from anaconda.org:

```
conda install single_cell_pipeline
```

## 1. Align

![align](readme_data/align.png)

The single cell analysis pipeline runs on a list of pairs of fastq file (paired end) and performs the following steps:

* Run fastqc on the fastq files
* Align the fastq pairs with bwa (supports both mem and aln)
* Merge all lanes together
* create targets and realign around those targets with GATK
* sort bams files, mark duplicate reads and index final bams
* collect metrics from samtools flagstat, picardtools CollectInsertSizeMetrics, picardtools CollectWgsMetrics and picardtools CollectGcBiasMetrics
* generate QC plot

### Input
The pipeline accepts a yaml file as input. The yaml file contains the input paths and metadata for each cell, the format for each cell is as follows:
```
SA12345-A12345-R01-C01:
  bam: /path/to/aligned/SA12345-A12345-R01-C01.bam
  column: 01
  condition: A
  fastqs:
    LANE_ID_1:
      fastq_1: /path/to/fastqfile/ACTACT-AGTAGT_1.fastq.gz
      fastq_2: /path/to/fastqfile/ACTACT-AGTAGT_1.fastq.gz
      sequencing_center: CENTERID
      sequencing_instrument: INSTRUMENT_TYPE
    LANE_ID_2:
      fastq_1: /path/to/fastqfile/ATTATT-ACTACT_1.fastq.gz
      fastq_2: /path/to/fastqfile/ATTATT-ACTACT_1.fastq.gz
      sequencing_center: CENTERID
      sequencing_instrument: INSTRUMENT_TYPE
  img_col: 10
  index_i5: i5-INDEX
  index_i7: i7-INDEX
  pick_met: CELLCALL
  primer_i5: ACTACTATT
  primer_i7: AGTAGTACT
  row: 01
  sample_type: C
```

### Run 

```
single_cell align \
--input_yaml inputs/SC-1234/bams.yaml \
--tmpdir temp/SC-1234/tmp \
--pipelinedir pipeline/SC-1234  \
--out_dir results/SC-1234/results \
 --library_id A123123 \
 ...
```


## 2. Hmmcopy

![align](readme_data/hmmcopy.png)

Hmmcopy runs the hmmcopy package on the aligned data and calls copy number:

1. generate read count wig files from the bam files
2. perform  GC correction
3. Run Hmmcopy
4. generate segment and bias plots, kernel density plot and heatmaps 

### Input
The pipeline accepts a yaml file as input. The yaml file contains the input paths and metadata for each cell, the format for each cell is as follows:
```
SA12345-A12345-R01-C01:
  bam: /path/to/aligned/SA12345-A12345-R01-C01.bam
  column: 01
  condition: A
  img_col: 10
  index_i5: i5-INDEX
  index_i7: i7-INDEX
  pick_met: CELLCALL
  primer_i5: ACTACTATT
  primer_i7: AGTAGTACT
  row: 01
  sample_type: C
```
The yaml file for the align step should work as well


### Run 

```
single_cell hmmcopy \
--input_yaml inputs/SC-1234/bams.yaml \
--tmpdir temp/SC-1234/tmp \
--pipelinedir pipeline/SC-1234  \
--out_dir results/SC-1234/results \
 --library_id A123123 \
 ...
```



## 3. Copyclone
### Input
The pipeline accepts a yaml file as input. The yaml file contains the input paths and metadata for each cell, the format for each cell is as follows:
```
SA12345-A12345-R01-C01:
  bam: /path/to/aligned/SA12345-A12345-R01-C01.bam
  column: 01
  condition: A
  img_col: 10
  index_i5: i5-INDEX
  index_i7: i7-INDEX
  pick_met: CELLCALL
  primer_i5: ACTACTATT
  primer_i7: AGTAGTACT
  row: 01
  sample_type: C
```


### Run 

```
single_cell copyclone \
--input_yaml inputs/SC-1234/bams.yaml \
--tmpdir temp/SC-1234/tmp \
--pipelinedir pipeline/SC-1234  \
--out_dir results/SC-1234/results \
 --library_id A123123 \
 ...
```


## 4. Aneufinder
### Input
The pipeline accepts a yaml file as input. The yaml file contains the input paths and metadata for each cell, the format for each cell is as follows:
```
SA12345-A12345-R01-C01:
  bam: /path/to/aligned/SA12345-A12345-R01-C01.bam
  column: 01
  condition: A
  img_col: 10
  index_i5: i5-INDEX
  index_i7: i7-INDEX
  pick_met: CELLCALL
  primer_i5: ACTACTATT
  primer_i7: AGTAGTACT
  row: 01
  sample_type: C
```
The yaml file for the align step should work as well


### Run 

```
single_cell aneufinder \
--input_yaml inputs/SC-1234/bams.yaml \
--tmpdir temp/SC-1234/tmp \
--pipelinedir pipeline/SC-1234  \
--out_dir results/SC-1234/results \
 --library_id A123123 \
 ...
```



## 5. merge bams
The tumour needs to be simultaneously merged across cells and split by region. The input for this step is the per cell bam yaml and the template for the merged bams by region.

```
SA501X5XB00877-A95670A-R04-C03:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C03.bam
SA501X5XB00877-A95670A-R04-C05:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C05.bam
SA501X5XB00877-A95670A-R04-C07:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C07.bam
SA501X5XB00877-A95670A-R04-C09:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C09.bam
SA501X5XB00877-A95670A-R04-C10:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C10.bam
```


```
single_cell merge_bams \
--input_yaml inputs/SC-1234/bams.yaml \
--merged_bam_template temp/SC-1234/pseudo_wgs/TUMOUR_{region}.bam \
--tmpdir temp/SC-1234/tmp \
--pipelinedir pipeline/SC-1234  \
--out_dir results/SC-1234/results \
 ...
```


## 6. split bams

The normal also needs to be split by region from an input data path to an output per region template.
```
single_cell split_bam \
 --wgs_bam data/NORMAL.bam \
 --split_bam_template temp/SC-1234/pseudo_wgs/NORMAL_{region}.bam \
 --tmpdir temp/SC-1234/tmp \
 --pipelinedir pipeline/SC-1234 \
 --out_dir results/SC-1234/results \
```


## 7. Variant Calling

```
SA501X5XB00877-A95670A-R04-C03:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C03.bam
SA501X5XB00877-A95670A-R04-C05:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C05.bam
SA501X5XB00877-A95670A-R04-C07:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C07.bam
SA501X5XB00877-A95670A-R04-C09:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C09.bam
SA501X5XB00877-A95670A-R04-C10:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C10.bam
```

The variant calling takes in both the per cell bam yaml, using the per cell bams for variant allele counting, and the tumour and normal region templates for calling snvs in parallel by region.
```
single_cell variant_calling \
 --input_yaml inputs/SC-1234/bams.yaml \
 --tumour_template temp/SC-1234/pseudo_wgs/TUMOUR_{region}.bam \
 --normal_template temp/SC-1234/pseudo_wgs/NORMAL_{region}.bam \
 --tmpdir temp/SC-1234/tmp \
 --out_dir results/SC-1099/results \
 --pipelinedir pipeline/SC-1234 \
 ...
```


## 8. Copy Number Calling

The copy number analysis takes in per cell bam yaml, with the ability to provide per cell bams for normal and tumour. If a WGS bam is used for normal that can be a single entry in the normal bam yaml.


per cell bam yaml format:
```
SA501X5XB00877-A95670A-R04-C03:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C03.bam
SA501X5XB00877-A95670A-R04-C05:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C05.bam
SA501X5XB00877-A95670A-R04-C07:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C07.bam
SA501X5XB00877-A95670A-R04-C09:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C09.bam
SA501X5XB00877-A95670A-R04-C10:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C10.bam
```

single wgs bam format:
```
SA501:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C03.bam
```

Run:
```
single_cell copy_number_calling \
 --tumour_yaml inputs/SC-1234/bams.yaml \
 --normal_yaml inputs/SC-1234/normal.yaml \
 --clone_id 1 \
 --tmpdir temp/SC-1234/tmp \
 --pipelinedir pipeline/SC-1234 \
 --out_dir results/SC-1234/results
 ...
```

## 9. Germline Calling

input yaml:
```
SA501X5XB00877-A95670A-R04-C03:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C03.bam
SA501X5XB00877-A95670A-R04-C05:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C05.bam
SA501X5XB00877-A95670A-R04-C07:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C07.bam
SA501X5XB00877-A95670A-R04-C09:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C09.bam
SA501X5XB00877-A95670A-R04-C10:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C10.bam
```
Run:
```
single_cell germline_calling \
 --input_yaml inputs/SC-1234/bams.yaml \
 --input_template temp/SC-1234/pseudo_wgs/INPUT_{region}.bam \
 --tmpdir temp/SC-1234/tmp \
 --pipelinedir pipeline/SC-1234 \
 --out_dir results/SC-1234/results
 ...
```


## 10. Destruct

```
SA501X5XB00877-A95670A-R04-C03:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C03.bam
SA501X5XB00877-A95670A-R04-C05:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C05.bam
SA501X5XB00877-A95670A-R04-C07:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C07.bam
SA501X5XB00877-A95670A-R04-C09:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C09.bam
SA501X5XB00877-A95670A-R04-C10:
  bam: data/single_cell_indexing/bam/A95670A/grch37/bwa-aln/SA501X5XB00877-A95670A-R04-C10.bam
```

The breakpoint analysis takes in per cell bam yaml in addition to the unsplit matched normal bam filename.
```
single_cell breakpoint_calling \
 --input_yaml inputs/SC-1234/bams.yaml \
 --matched_normal data/NORMAL.bam \
 --tmpdir temp/SC-1234/tmp \
 --pipelinedir pipeline/SC-1234 \
 --out_dir results/SC-1234/results
 ...
```


## 11. Generate Config 

The pipeline auto generates a config file with the default parameters before every run. Some of the values in the config file can be updated by using the ``--config_override`` option.  ```generate_config``` option allows users to generate the config files. These configs can then be specified as input to the pipeline after making the required changes.
```
single_cell generate_config --pipeline_config pipeline_config.yaml --batch_config batch_config.yaml
```
the pipeline config file contains all pipeline defaults and the batch config specifies all azure batch specific settings.

the pipeline config can be specified manually when running the pipeline with ```--config_file``` option and the batch config with ```--submit_config``` option.



## 12. Clean Sentinels
the pipeline will skip any successful tasks from previous runs when run again. The ``--rerun`` flag force run all tasks including the successful tasks from the previous runs while the ```clean_sentinels``` option provides a more fine grained control.

```
single_cell clean_sentinels --mode list --pattern "*" --tmpdir temp/SC-1234/tmp/hmmcopy
```
will list all successfully completed hmmcopy tasks

```
single_cell clean_sentinels --mode delete --pattern "*" --tmpdir temp/SC-1234/tmp/hmmcopy
```
is the same as ```--rerun``` flag

running
```
single_cell clean_sentinels --mode delete --pattern "*plot_heatmap*" --tmpdir temp/SC-1234/tmp/hmmcopy
```
before launching the hmmcopy will rerun the heatmap plotting  and any tasks that haven't completed yet.


### Common options

* add ``` --storage azureblob``` to read the inputs and write the outputs to Microsoft Azure Blob Storage
* add ``` --submit local``` to run the pipeline on the current node
* add ``` --submit asyncqsub ``` to run the pipeline on a SGE cluster
* add ``` --submit azurebatch``` to run the tasks on Microsoft Azure Batch
* add ``` --nocleanup``` to disable temporary file deletion
* add ``` --loglevel DEBUG``` for verbose logging
* ```--submit_config <path to yaml config>``` for custom azure batch submission configuration file. see 10. Generate config for more details
* ```--config_file <path to yaml config>``` to specify pipeline config file.

## 13. Demultiplex BAM

![align](readme_data/align.png)

Demultiplexes a BAM file by cell identifier and converts the demultiplexed bam files to paired end fastq files.

* Demultiplex bam files based on the optional CB (cell identifier) tag.
* User shall specify list of cell ids/barcodes on which to demultiplex.
* Aligned reads associated with barcodes not contained in user specified list must be filtered out.
* Aligned reads that don't have CB tags must be assigned the special tag="undetermined".
* Reconstruct the paired reads for each originating cell and output paired reads to _1 and _2 FASTQ files.

### Input
All input is via command line arguments. 

### Run 

```
single_cell demultiplex_bam \
--bam tests/data/demultiplex_bam/bj_mkn45_10pct_possorted_bam_10k_snippet.bam \
--out_dir tests/data/demultiplex_bam/output/ \
--barcode_csv tests/data/demultiplex_bam/bj_mkn45_10pct_per_cell_summary_metrics.csv \
--submit local \
--tmpdir temp \
--pipelinedir pipeline \
--maxjobs 1 \
--config_file config.yaml \
--loglevel DEBUG
```

 