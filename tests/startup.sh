git clone https://dgrewal@bitbucket.org/dgrewal/single_cell_pipeline.git


sudo apt-get update
sudo apt-get -y install bzip2 wget libkeyutils-dev ssh ttf-dejavu fontconfig vim make build-essential libpng-dev zlib1g-dev

#instal azure cli
AZ_REPO=$(lsb_release -cs)
echo "deb [arch=amd64] https://packages.microsoft.com/repos/azure-cli/ $AZ_REPO main" | sudo tee /etc/apt/sources.list.d/azure-cli.list
curl -L https://packages.microsoft.com/keys/microsoft.asc | sudo apt-key add -
sudo apt-key adv --keyserver packages.microsoft.com --recv-keys 52E16F86FEE04B979B07E28DB02C46DF417A0893
sudo apt-get install apt-transport-https
sudo apt-get update && sudo apt-get install azure-cli


sudo rm -rf /usr/local/miniconda2
wget -q https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
sudo bash Miniconda2-latest-Linux-x86_64.sh -b -p /usr/local/miniconda2
echo "export PATH=/usr/local/miniconda2/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
sudo chmod -R 777 /usr/local/miniconda2
export PATH="/usr/local/miniconda2/bin":$PATH
conda upgrade conda -y
conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels 'bioconda'
conda config --add channels 'r'
conda config --add channels 'conda-forge'
conda config --add channels https://conda.anaconda.org/aroth85


#conda install --file single_cell_pipeline/tests/conda_packages.txt -y -q

gatk3-register vmresources/GenomeAnalysisTK.jar
pip install azure-storage azure-batch futures azure-mgmt

conda install networkx pyyaml dill -y -q 


