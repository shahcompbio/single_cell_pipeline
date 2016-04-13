# NextSeq single-cell alignment pipeline

The NextSeq alignment pipeline is designed for multiplexed low-depth single-cell data in BCL format. It processes and aligns raw sequencing data and outputs both a table of sequencing metrics and some QC plots.

## Inputs

The only required input to the pipeline is a `config.yaml` file, which contains information necessary for the run.

### Config

The `config.yaml` file specifies a NextSeq directory (which must include a `SampleSheet.csv` file) and an output directory, as well as several other paths and parameters which do not need to be modified for every run. 

Note that the `SampleSheet.csv` file must be formatted according to specific requirements that go beyond Illumina's standards. Example `config.yaml` and `SampleSheet.csv` files are found in the `/alignment/templates` directory.

## Requirements

### Python

Python 2.7 is required and should be on the path. The [Anaconda](https://docs.continuum.io/anaconda) distribution is recommended.

* Python packages:
	* [pipelines](https://bitbucket.org/aroth85/pipelines)
	* [ruffus](http://www.ruffus.org.uk)
	* [drmaa](https://pypi.python.org/pypi/drmaa) (required for running pipeline in cluster mode)
	* [pandas](http://pandas.pydata.org)
	* [matplotlib](http://matplotlib.org) (included in Anaconda distribution, required for generating plots)
	* [seaborn](https://stanford.edu/~mwaskom/software/seaborn) (required for generating plots)

### Illumina's bcl2fastq2

[bcl2fastq](http://support.illumina.com/downloads/bcl2fastq-conversion-software-v217.html) is required and should be on the path.

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc) is required and should be on the path.

### BWA

[BWA](http://bio-bwa.sourceforge.net) is required and should be on the path.

### Samtools

[Samtools](http://www.htslib.org) is required and should be on the path.

### Picard

[Picard](http://broadinstitute.github.io/picard) is required. A path must be specified in the config.yaml file.

## Running the pipeline on the Shah lab infrastructure

After installation, please follow the steps below for each pipeline run:

* Download the NextSeq run from the UBC BRC web server to lustre archive

		cd /share/lustre/archive/single_cell_indexing/NextSeq/bcl
		wget -m -np -nH --cut-dirs=1 --reject="index.html*" <webserver_address>/sequencing/<run_id>/

* Create a `SampleSheet.csv` file that conforms to specifications, and save it in the NextSeq run directory on lustre archive
* Use rsync to transfer the NextSeq run directory to genesis scratch space

		ssh xhost
		ssh thost05
		rsync --dry-run -avzL <your_user_id>@<beast_ip_address>:/share/lustre/archive/single_cell_indexing/NextSeq/bcl/<run_id> /genesis/extscratch/shahlab/<your_user_id>/<some_directory>

* Create a run output directory on genesis scratch space and copy the `config.yaml` template and `run_pipeline.sh` template to it

		mkdir /genesis/extscratch/shahlab/<your_user_id>/<run_analysis_directory>
		cp /path/to/pipeline/installation/single_cell_nextseq/alignment/templates/config.yaml /genesis/extscratch/shahlab/<your_user_id>/<analysis_directory>
		cp /path/to/pipeline/installation/single_cell_nextseq/alignment/templates/run_pipeline.sh /genesis/extscratch/shahlab/<your_user_id>/<analysis_directory>

* Specify the paths to the NextSeq run directory and output directory in the `config.yaml` file
* Specify the path to the `config.yaml` file in the `run_pipeline.sh` file
* Open a screen on genesis

		screen -S <name_of_screen>

* Source your bashrc file and ensure the required programs are on the path

		source ~/.bashrc
		which python

* Launch the pipeline from your run directory

		cd /genesis/extscratch/shahlab/<your_user_id>/<run_analysis_directory>
		bash run_pipeline.sh

* You can check on your jobs with

		qstat -u <your_user_id>

* You can check for pipeline completion status with

		cat run_pipeline.error.txt

* After successful completion, rsync the output files to your directory on lustre projects

		rsync --dry-run -avzL --exclude='tmp/' /genesis/extscratch/shahlab/<your_user_id>/<analysis_directory> <your_user_id>@<beast_ip_address>:/share/lustre/projects/<some_directory>

* Note that you should NOT rsync the NextSeq directory back, we will only store the original BCL files and `SampleSheet.csv` on lustre archive
* Open the metric table and figures on lustre and check that they look sensible
* Clean up your temporary genesis directories
* Close your screen

		exit

## Notes

* If the index sequences are specified incorrectly in `SampleSheet.csv`, the pipeline will not be aware of this and will proceed to the alignment steps. If you suspect that demultiplexing was not successful, you can check the file sizes with

	du -sh /genesis/extscratch/shahlab/<your_user_id>/<some_directory>/<run_id>/Data/Intensities/BaseCalls/*.fastq.gz

* File sizes should be in Mb for all non-NTC samples, and in Kb for NTC samples. It is normal for the Undetermined file to be larger than any given single-cell file
* If you find that demultiplexing was indeed unsuccessful, you MUST do the following:
	* Delete the NextSeq run directory on genesis
	* Fix the SampleSheet.csv file on lustre archive
	* Re-transfer the NextSeq directory to genesis and launch the pipeline again
* It is very important that the correct `SampleSheet.csv` is saved on archive for future reference