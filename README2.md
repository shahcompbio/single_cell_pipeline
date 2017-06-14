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

## Run the pipeline

The pipeline takes 3 arguments.  A config file (example for genesis provided),
the nextseq input directory, and the output directory.

Run the pipeline as follows:

```
```
