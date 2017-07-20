# Single Cell Pipeline

## Setup and installation

Set up conda with the required packages.

First ensure you have the correct channels:

```
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels 'bioconda'
conda config --add channels 'r'
```

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

## Run the pipeline

### Inputs

TODO: cell id or sample id???

#### Configuration file

The configuration file contains global options including:

* reference genome
* chromosomes to operate on
* hmm parameters
* thresholds
* adapters

See config_shahlab.yaml and config_cloud.yaml for examples.

#### Input sample info table

The input sample info table provides meta data per cell.  The columns should be:

- cell_id
- cell_call
- experimental_condition
- sample_type
- sample_well
- sample_plate
- i5_barcode
- i7_barcode

#### Input fastq table

The fastq file paths should be provided in a CSV file, and the path to this input
file given on the command line.  The following columns are required:

- sample_id
- lane_id
- fastq_1
- fastq_2

#### Library id

The library id should be provided as an argument on the command line and will be
added to the header of the bam file.

### Command line interface

The pipeline takes 4 positional arguments:

1. input sample info table
2. input fastq filename table
3. library id
4. output directory
5. config filename

The remaining arguments are for controlling execution using [pypeliner](http://pypeliner.readthedocs.org/).

For example, run the pipeline as follows:

```
single_cell \
  /shahlab/amcpherson/single_cell_test/A90696ABC_sample_info.csv \
  /shahlab/amcpherson/single_cell_test/A90696ABC_fastqs.csv \
  A90696ABC \
  /shahlab/amcpherson/single_cell_test/results/ \
  /shahlab/dgrewal/test_andrew_sc/data/config_shahlab.yaml \
  --loglevel DEBUG \
  --submit asyncqsub \
  --maxjobs 1000 \
  --nocleanup
```

## Setup in the cloud

Run the following in an instance:

```
sudo mkdir /mnt/software
sudo chown shahlab:shahlab /mnt/software
cd /mnt/software/
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p /mnt/software/miniconda2

source ~/.bashrc

# conda environment
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels 'bioconda'
conda config --add channels 'r'

source activate singlecellpipeline

# Install GATK
wget -O GATK.tar.bz2 https://software.broadinstitute.org/gatk/download/auth?package=GATK
tar -jxvf GATK.tar.bz2
gatk-register GenomeAnalysisTK.jar

git clone https://amcpherson@svn.bcgsc.ca/bitbucket/scm/sc/single_cell_pipeline.git
cd single_cell_pipeline/
conda create --name singlecellpipeline --file conda_packages.txt --yes
python setup.py develop

# Download reference genome
sudo mkdir /mnt/refdata
sudo chown shahlab:shahlab /mnt/refdata
cd /mnt/refdata/
wget ftp://ftp.bcgsc.ca/public/shahlab/singlecellpipeline/*
picard CreateSequenceDictionary R= /mnt/refdata/GRCh37-lite.fa O= /mnt/refdata/GRCh37-lite.dict

# Create analysis space
sudo mkdir /mnt/analysis
sudo chown shahlab:shahlab /mnt/analysis
```

Copy the fastq files.  Run the following on thost:

```
scp -r /genesis/shahlab/dgrewal/test_andrew_sc/data/hiseq/PX0577/ sccompute:/mnt/analysis/
scp -r /genesis/shahlab/amcpherson/single_cell_test/*.csv sccompute:/mnt/analysis/
```

Reformat the fastq list:

```
cd /mnt/analysis
sed 's#/genesis/shahlab/dgrewal/test_andrew_sc/data/hiseq#/mnt/analysis#g' < A90696ABC_fastqs.csv > A90696ABC_fastqs_local.csv
```

Subset:

```
head -1 A90696ABC_fastqs_local.csv > A90696ABC_fastqs_local_subset.csv
grep '\(R03-C63\|R08-C63\|R04-C63\|R11-C63\|R13-C63\|R15-C63\|R20-C63\|R45-C63\)' A90696ABC_fastqs_local.csv >> A90696ABC_fastqs_local_subset.csv
head -1 A90696ABC_sample_info.csv > A90696ABC_sample_info_subset.csv
grep '\(R03-C63\|R08-C63\|R04-C63\|R11-C63\|R13-C63\|R15-C63\|R20-C63\|R45-C63\)' A90696ABC_sample_info.csv >> A90696ABC_sample_info_subset.csv
```

Run the pipeline:

```
single_cell_nextseq \
  /mnt/analysis/A90696ABC_sample_info_subset.csv \
  /mnt/analysis/A90696ABC_fastqs_local_subset.csv \
  A90696ABC \
  /mnt/analysis/results/ \
  /mnt/software/single_cell_pipeline/config_cloud.yaml \
  --loglevel DEBUG \
  --submit local \
  --maxjobs 32 \
  --nocleanup \
  --realign
```




