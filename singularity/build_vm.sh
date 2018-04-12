echo "GIT_SSL_NO_VERIFY=1" >> ~/.bashrc

sudo apt-get update
sudo apt-get -y install bzip2 wget libkeyutils-dev ssh ttf-dejavu fontconfig vim make build-essential libpng-dev zlib1g-dev

#instal azure cli
AZ_REPO=$(lsb_release -cs)
echo "deb [arch=amd64] https://packages.microsoft.com/repos/azure-cli/ $AZ_REPO main" | sudo tee /etc/apt/sources.list.d/azure-cli.list
curl -L https://packages.microsoft.com/keys/microsoft.asc | sudo apt-key add -
sudo apt-key adv --keyserver packages.microsoft.com --recv-keys 52E16F86FEE04B979B07E28DB02C46DF417A0893
sudo apt-get install apt-transport-https
sudo apt-get update && sudo apt-get install azure-cli

#download data
cd ~
az storage blob download -c vmresources -f vmresources.tar.gz -n vmresources.tar.gz
tar -xvf vmresources.tar.gz


cd ~/vmresources
git clone -b master https://dgrewal@svn.bcgsc.ca/bitbucket/scm/sc/single_cell_pipeline.git
git clone -b azure_production https://dgrewal@bitbucket.org/dranew/pypeliner.git
cd ~


#conda install
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
sudo bash Miniconda2-latest-Linux-x86_64.sh -b -p /usr/local/miniconda2
echo "export PATH=/usr/local/miniconda2/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
sudo chmod -R 777 /usr/local/miniconda2
conda upgrade conda -y
 conda config --add channels https://conda.anaconda.org/dranew
conda config --add channels 'bioconda'
conda config --add channels 'r'
conda config --add channels 'conda-forge'
conda config --add channels https://conda.anaconda.org/aroth85

conda install --file ~/vmresources/conda_packages.txt -y

gatk-register vmresources/GenomeAnalysisTK.jar
pip install azure-storage azure-batch futures
cd ~/vmresources/pypeliner; python setup.py install
cd ~/vmresources/single_cell_pipeline; python setup.py install
sudo    mv ~/vmresources/museq /usr/local/museq
cd /usr/local/museq/
make clean
make PYTHON=python BOOSTPATH=~/vmresources/boost_1_57_0/
cd ~

pip install scikit-learn==0.13.1
conda install strelka==1.0.14 -y



rm -rf vmresources vmresources.tar.gz Miniconda2-latest-Linux-x86_64.sh
