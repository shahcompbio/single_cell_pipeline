set -e

docker run -e ANACONDA_API_TOKEN  -v $PWD:$PWD -w $PWD singlecellpipeline/conda_build:v0.0.1 bash -c "conda update -n base -c defaults conda -y &&  conda config --set anaconda_upload yes && conda build -c bioconda -c shahcompbio conda/single_cell_pipeline_align"

docker run -e ANACONDA_API_TOKEN  -v $PWD:$PWD -w $PWD singlecellpipeline/conda_build:v0.0.1 bash -c "conda update -n base -c defaults conda -y &&  conda config --set anaconda_upload yes && conda build conda/single_cell_pipeline_hmmcopy -c bioconda -c shahcompbio -c conda-forge"

docker run -e ANACONDA_API_TOKEN  -v $PWD:$PWD -w $PWD singlecellpipeline/conda_build:v0.0.1 bash -c "conda update -n base -c defaults conda -y &&  conda config --set anaconda_upload yes && conda build conda/single_cell_pipeline_annotation  -c shahcompbio"

docker run -e ANACONDA_API_TOKEN  -v $PWD:$PWD -w $PWD singlecellpipeline/conda_build:v0.0.1 bash -c "conda update -n base -c defaults conda -y &&  conda config --set anaconda_upload yes && conda build conda/single_cell_pipeline_pseudo_bulk_qc -c bioconda -c shahcompbio"
