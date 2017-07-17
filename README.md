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

Add the single cell nextseq pipeline into the current site packages:

```
python setup.py develop
```

Use develop mode to allow for modifying of the code during the debugging phase.
Later versions will be released through conda.

## Run the pipeline

### Inputs

#### Configuration file

The configuration file contains global options including:

* reference genome
* chromosomes to operate on
* hmm parameters
* thresholds
* adapters

See config_shahlab.yaml and config_cloud.yaml for examples.

#### Input Fastq Directory

The fastqs are assumed to be in a single directory with a specific structure.

```
/.../{sequencing_library_id}/{lane_id}/{sequencing_library_id}_{index_1}-{index_2}/{lane_id}_{index_1}-{index_2}_{read_end}.fastq.gz
```

The directory provided should be `/.../{sequencing_library_id}`, and `sequencing_library_id` will be inferred from the directory name.

#### Lanes to analyze

The set of lanes to be analyzed should be provided as a list of `lane_id`.

#### Sample Sheet

The sample sheet provides information about the indices (`index_1` and `index_2`) of each cell, in
addition to per cell metadata.  The library id is also extracted from the header and added to bam
header.

### Command line interface

The pipeline takes 4 positional arguments:

1. input fastq directory
2. sample sheet filename
3. output directory
4. config filename

In addition a list of lanes should be provided using `--lanes`.

The remaining arguments are for controlling execution using [pypeliner](http://pypeliner.readthedocs.org/).

For example, run the pipeline as follows:

```
single_cell_nextseq \
  /genesis/shahlab/dgrewal/test_andrew_sc/data/hiseq/PX0577/ \
  /genesis/shahlab/dgrewal/test_andrew_sc/data/hiseq/SampleSheet_CB643ANXX_3.csv \
  /genesis/shahlab/dgrewal/test_andrew_sc/test_output \
  /shahlab/dgrewal/test_andrew_sc/data/config_shahlab.yaml \
  --lanes CB643ANXX_3 CB8HDANXX_4 CB8HDANXX_5 \
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

conda create --name singlecellpipeline python=2.7 \
  hmmcopy_utils \
  samtools \
  bwa \
  fastqc \
  picard \
  cutadapt \
  trim-galore \
  pypeliner \
  bioconductor-hmmcopy \
  r-plyr \
  r-getopt \
  pyyaml \
  pandas \
  statsmodels
  

source activate singlecellpipeline

git clone https://dranew@bitbucket.org/dranew/single_cell_nextseq.git
cd single_cell_nextseq/
python setup.py develop

# Download reference genome
sudo mkdir /mnt/refdata
sudo chown shahlab:shahlab /mnt/refdata
cd /mnt/refdata/
wget ftp://ftp.bcgsc.ca/public/shahlab/singlecellpipeline/*

# Create analysis space
sudo mkdir /mnt/analysis
sudo chown shahlab:shahlab /mnt/analysis
```

Copy the fastq files.  Run the following on thost:

```
scp -r /shahlab/amcpherson/single_cell_nextseq1/test_output_new/fastq sccompute:/mnt/analysis/
scp -r /shahlab/archive/single_cell_indexing/NextSeq/161230_NS500668_0153_AHHHWJAFXX/SampleSheet.csv sccompute:/mnt/analysis/
```




