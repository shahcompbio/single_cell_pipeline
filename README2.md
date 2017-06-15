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
  pypeliner
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
