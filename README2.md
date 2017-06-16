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
  r-getopt
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

The pipeline takes 3 arguments.  A config file (example for genesis provided),
the nextseq input directory, and the output directory.

Run the pipeline as follows:

```
single_cell_nextseq \
  /shahlab/archive/single_cell_indexing/NextSeq/161230_NS500668_0153_AHHHWJAFXX/ \
  /shahlab/amcpherson/single_cell_nextseq1/test_output_new/ \
  /shahlab/amcpherson/single_cell_nextseq1/config_shahlab_new.yaml \
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




