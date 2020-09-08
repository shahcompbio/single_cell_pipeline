## Build:


```
 docker run -it -v $PWD:$PWD -w $PWD singlecellpipeline/conda_build:v0.0.1 bash
```

```
conda update -n base -c defaults conda -y
conda build -c shahcompbio -c bioconda single_cell_pseudo_bulk_qc/
```



## Install:

```
conda install -c bioconda single_cell_pseudo_bulk_qc
pip install -r requirements.txt
```
