## Build and Upload Conda Packages 
This directory contains recipes for building conda packages. Here are the steps of building, uploading and installing the conda packages.

#### 1. Launch a docker shell:
cd into the current directory (```single_cell_pipeline/conda```) and launch a docker shell by running: 
```docker run -it -v $PWD:$PWD -w $PWD singlecellpipeline/conda_build:v0.0.1 bash```

#### 2. To build: 
##### 2.1 For alignement workflow: 

```conda build single_cell_pipeline_align -c bioconda -c shahcompbio```

##### 2.2 For hmmcopy workflow: 

```conda update -n base -c defaults conda -y && conda build single_cell_pipeline_hmmcopy -c bioconda -c shahcompbio``` 

##### 2.3 For annotation workflow: 

```conda update -n base -c defaults conda -y && conda build single_cell_pipeline_annotation -c shahcompbio``` 

#### 3. To upload: 
Copy and run the command printed in the console. It looks like:
```anaconda upload /opt/conda/conda-bld/linux-64/single_cell_pipeline_annotation-0.6.8-py37_1.tar.bz2 --force```

#### (Optional) 4. To install the conda package: 
Create a conda environment and run 

4.1 For alignment: ```conda update -n base -c defaults conda -y && conda install -c shahcompbio -c conda-forge single_cell_pipeline_align ```

4.2 For hmmcopy: ```conda update -n base -c defaults conda -y && conda install -c bioconda -c shahcompbio -c dranew -c conda-forge single_cell_pipeline_hmmcopy```

4.3 For annotation: ```conda update -n base -c defaults conda -y && conda install -c conda-forge -c bioconda -c shahcompbio single_cell_pipeline_annotation -y ```

