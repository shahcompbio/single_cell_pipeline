# Single Cell Pipeline


![build](https://www.bcgsc.ca/bamboo/plugins/servlet/wittified/build-status/SC-SCP)
## Setup and installation

Set up conda with the required packages.

First ensure you have the correct channels:

```
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels https://conda.anaconda.org/aroth85
conda config --add channels 'bioconda'
conda config --add channels 'r'
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
- sample_type # No longer necessary
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
sudo mkdir /datadrive/software
sudo chown shahlab:shahlab /datadrive/software
cd /datadrive/software/
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p /datadrive/software/miniconda2
echo "export PATH=/datadrive/software/miniconda2/bin:\$PATH" >> ~/.bashrc
echo "export GIT_SSL_NO_VERIFY=1" >> ~/.bashrc
source ~/.bashrc

# Conda setup
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels 'bioconda'
conda config --add channels 'r'

# Setup environment
cd /datadrive/software/
git clone https://amcpherson@svn.bcgsc.ca/bitbucket/scm/sc/single_cell_pipeline.git
cd single_cell_pipeline/
conda create --name singlecellpipeline --file conda_packages.txt --yes
source activate singlecellpipeline
python setup.py develop

# Install GATK
wget -O GATK.tar.bz2 https://software.broadinstitute.org/gatk/download/auth?package=GATK
tar -jxvf GATK.tar.bz2
gatk-register GenomeAnalysisTK.jar

# Download reference genome
sudo mkdir /datadrive/refdata
sudo chown shahlab:shahlab /datadrive/refdata
cd /datadrive/refdata/
wget ftp://ftp.bcgsc.ca/public/shahlab/singlecellpipeline/*
picard CreateSequenceDictionary R= /datadrive/refdata/GRCh37-lite.fa O= /datadrive/refdata/GRCh37-lite.dict

# Create analysis space
sudo mkdir /datadrive/analysis
sudo chown shahlab:shahlab /datadrive/analysis
```

Copy the fastq files.  Run the following on thost:

```
scp -r /genesis/shahlab/dgrewal/test_andrew_sc/data/hiseq/PX0577/ sccompute1:/datadrive/analysis/
scp -r /genesis/shahlab/amcpherson/single_cell_test/*.csv sccompute1:/datadrive/analysis/
```

Reformat the fastq list:

```
cd /datadrive/analysis
sed 's#/genesis/shahlab/dgrewal/test_andrew_sc/data/hiseq#/datadrive/analysis#g' < A90696ABC_fastqs.csv > A90696ABC_fastqs_local.csv
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
  /datadrive/analysis/A90696ABC_sample_info_subset.csv \
  /datadrive/analysis/A90696ABC_fastqs_local_subset.csv \
  A90696ABC \
  /datadrive/analysis/results/ \
  /datadrive/software/single_cell_pipeline/config_cloud.yaml \
  --loglevel DEBUG \
  --submit local \
  --maxjobs 32 \
  --nocleanup \
  --realign
```




