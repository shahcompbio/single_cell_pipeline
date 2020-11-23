conda create -y -n bwa_0_7_17 -c bioconda bwa=0.7.17
conda create -y -n picard -c bioconda picard
conda create -y -n trimgalore -c bioconda trim-galore
conda create -y -n scp_versionx -c bioconda -c shahcompbio  python drmaa openssl=1.0 matplotlib pandas=0.25.3 pypdf2 pypeliner pysam pyyaml seaborn statsmodels scikit-learn

