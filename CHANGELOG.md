# Change Log

### v0.8.14
1. bugfix in read attrition metrics

### v0.8.13
1. dtype for is_control in annotation

### v0.8.12
1. dtype for state should be float to NA values from failed cells

### v0.8.11
1. add subsampling prior to museq to remove read sinks

### v0.8.10
1. adding track type header to wig
2. hmmcopy: updating state to NA, allow failed cells in heatmap


### v0.8.9
1. adding is_control and read attrition metrics

### v0.8.8

### Changes:
1. annotation pipeline for variant calling can now be configured in config file

### v0.8.7

### Changes:
1. fastqscreen tags limited to [0,1]
2. fastqscreen now can filter out custom tag combinations
3.  supports list of references per organism 

### v0.8.6 

### Changes:
1. added ref and alt to infer haps output

### v0.8.5

###  changes:
1. merged development QC pipeline changes

### v0.8.4

### changes:
1. updated test datasets
2. added softclipped read filter to merge bams


### v0.8.3

### Changes:
1. updated annotation to handle configurable genomes

### v0.8.2

### Changes:
1. added filtering per genome in fastqscreen
2. more cutomizable genomes in fastqscreen
3. update pypeliner to 0.6.2

### v0.8.1

#### Changes:
1. moved hmmcopy R into repo

### v0.8.0

#### changes:
1. picard: quiet mode
2. pypeliner update to v0.6.1
3. docker: added azure libs
4. updated input file for snv genotyping
5. snv genotyping in testing now

#### bugs:
1. fastqscreen supports non gzipped fastqs now

### v0.7.6

#### Changes:
1. updated pypeliner to v0.6.0

### v0.7.5

#### Changes:
1. updated pypeliner

### v0.7.4

#### Changes:
1. added cohort_qc to testing
2. changes from andrew

### v0.7.3

#### Changes:
1. update destruct to v0.4.19

### v0.7.2

#### Changes:
1. networkx version for nreakpoint docker

### v0.7.1

#### Changes:
1. added pseudobulk QC to codebuild

### v0.7.0

 - Outputs do not change with this release
#### Changes:
1. deprecated conda package
2. deprecated docker in docker 
3. new docker containers (one per pipeline) and new org in quay.io


### v0.6.46:

####  Changes:
1. Added sample id and library id to alignment and hmm metrics

### v0.6.45:

#### Changes:
1. remove unused code in hmmcopy
2. review of QC codebase
3. pseudobulk QC is in codebuild

### v0.6.44

#### Changes:
1. raise exception when reference type is unknown
2. set all dtypes as str when reading maf 

#### Bug:
1. errors in cohort QC

### v0.6.43

#### changes: 
1. update to latest pypeliner (v0.6.27)

### v0.6.42

#### Changes
1. produces all outputs necessary to load cbioportal with cna + an oncoplot from maftools
2. updated HMMcopy to R=4

### v0.6.41

#### Changes
1. trim and sequencing center is now a pipeline level flag

#### Bug
1. alignment and hmmcopy tarballs were broken

### v0.6.40

#### Bug:
1. bugfix: trim field in seqinfo incorrect (#119)

#### Changes:
1. remove mem_retry_factor overrides in pipeline
2. deprecation: sequencing instrument in input.yaml is deprecated, replaced with trim (boolean)

### v0.6.39

Bugfix:
1. Error in fastqscreen when removing contaminated reads

### v0.6.38

Changes:
1. destruct updated to remove secondary reads
2. updated pseudobulk QC on jun o (singularity)
3. updated pseudobulk plots
4. updated pseudobulk documentation



### v0.6.37

Changes:
1. fixed snpeff call in  germline calling
2. updated docs


### v0.6.36

Changes:
1. updated hmmcopy

### v0.6.35

Changes:
1. updated pypeliner to v0.5.23 to support latest azure python sdk

### v0.6.34

Changes: 
1. update pypeliner to v0.5.22 which lowers lsf query volume

### v0.6.33

### Bugs:
1. fixed issue with missing contamination table in annotation html

### v0.6.32

### Bugs:
1. csvutils couldnt handle tsv

### v0.6.31

### Bugs:
1. hmmcopy plots were hardcoded to human genomes

### v0.6.30

### Bugs:
1. update cell cycle classifier

### v0.6.29

### Bugs:
1. csvutils annotate_csv: annotation was in incorrect order


### v0.6.28

### Bugs:
1. added Trim to dtypes

### v0.6.27

### Changes:
1. pseudo_bulk_qc: allele data loading is done in chunks to improve memory usage
2. pseudo_bulk_qc: added rearrangement types for lumpy
3. pseudo_bulk_qc: additional qcplots
4. pseudo_bulk_qc: better plot organization 
5. pseudo_bulk_qc: resized report output
6. pseudo_bulk_qc: code can now handle incomplete data
7. pseudo_bulk_qc: placeholder text for missing plots


### v0.6.26

### Changes:
1. remove bwa aln

### v0.6.25

### Changes:
1. update pypeliner
2. deleted redundant dockerfiles

### v0.6.24

### Changes:
1. clean up pseudobulk QC
2. add pseudoculk QC to master build

### BUG:
1. unused ete3 code

### v0.6.23

### Changes:
1. conda packages: added drmaa, simpler recipes

### v0.6.22

### Changes:
1. pseudo bulk QC pipeline
2. remove ete3 dependency
3. updated conda recipe. removed ete3 added drmaa

### v0.6.21

### Bug:
1. KDE plot can't handle bad data

### v0.6.20

### Changes:
1. conda packages are now built on top on py 3.8

### Bug:
1. minor heatmap axis labeling bug (assumes chrom1 is starting point)


### v0.6.19

#### Changes:
1. produces a MT bam file
2. setup doc also includes information to switch pipeline to production datasets


### v0.6.18

#### Bugs:
1. missing data error in species classifier
2. build: pip cannot install from commits in a PR coming from a fork

#### Changes:
1. conda: annotation pipeline requires jinja2 dependency
2. conda: align requires ete3 and R
2. added a quickstart guide based on conda

#### v0.6.17

#### Bugs:
1. missing fastqscreen training data for mouse in config
2. data type bug in fastqscreen summary

### v0.6.16

#### Changes:
1. new build process (AWS codebuild, on PR)


### v0.6.15

#### Bug:
1. there were no counts for cells with no data in fastqscreen summary

### v0.6.14

#### Changes:
1. removed joblib files from repo
2. added conda recipes

### v0.6.13

#### Changes:
1. updated remixt


### v0.6.12

#### Bug:
1. .tmp in input.yaml file path in  metadata yaml file

#### Changes:
1. test data set for count haps
2. updated infer haps testing
3. removed hardcoded genomes in HTML report
4. infer and count haps use proper chromosomes

### v0.6.11

#### Changes:
1. pypeliner version update to v0.5.19



### v0.6.10

#### Bug:
1. missing cell id in count haps

### v0.6.9

#### Bug:
1. fixed snv annotation issues

#### Changes:
1. biowrappers upgraded to v0.2.8

### v0.6.8

#### Bug:
1. incorrect rebase on is_contaminated flag code
2. fixed issues with tag name in docker build

#### Changes:
1. updated QC generation code, added classifier
2. updated docs: read the docs compatible

### v0.6.7

#### Bug:
1. issue with dtypes for lumpy

#### Changes:
1. updated biowrappers to v0.2.7

### v0.6.6

#### Bug:
1. missing yaml extension for count haps input


### v0.6.5

#### Bug:
1. issues with missing yaml files

#### Changes:
1. add yaml validators
2. added docs


### v0.6.4:

#### Bug:
1. Issue with fastqscreen counts
2. remove corrupt tree files from metadata
3. dtype for chr is str not int in lumpy

### v0.6.3:

#### Bug:
1. issue with yaml missing in extensions for snv file

### v0.6.2:

#### Bug:
1. collect metrics introduced nan into int columns

### v0.6.1

### Changes:
1. per sample genotyping results

### v0.6.0

### Changes:
1. jenkins refactor
2. csvutils refactor

### V0.5.21

#### Changes:
1. upgraded remixt to v0.5.11

### v0.5.20

#### Changes:
1. upgraded destruct and remixt to latest pypeliner.

### v0.5.19

####Changes:
1. upgraded remixt and cell cycle classifier

### v0.5.18

#### Changes:
1. upgraded pypeliner to v0.5.18

### v0.5.17

#### Changes:
1. updated the calculation for is_contaminated
2. moved contaminated flagging to annotation
3. pypeliner updated to v0.5.17
4. removed org from scp docker container names in config

### v0.5.16

####Changes:
1. updated pypeliner to v0.5.16
2. conda build fixed


### v0.5.15

#### bug:
1. metadata regions are none sometimes

### v0.5.14

#### bug:
1. hdf to csv was only writing header for first output file


### v0.5.13

#### bug:
1. splitting heatmap into 1000 cells per page doesnt account for cases where the next page only has 1 cell which cant be plotted.

### v0.5.12

#### bug:
1. destruct read indexing was broken

#### Changes:
1. flags to selectively run destruct or lumpy

### v0.5.11:

#### Changes:
1. tempdir is now separated by subcommands
2. file names updated for breakpoint calling


### v0.5.10

#### Bug:
1. conversion to csv from hdf was overwriting not appending

#### Changes:
1. destruct version updated
2. docker container is built from tag, not from commit id

### v0.5.9

#### Bug:
1. excessive memory usage in h5 to csv conversion
2. fixed config issue for snpeff

### v0.5.8

#### Changes:
1. updated biowrappers to v0.2.4
2. updated pypeliner version to v0.5.15
3. optimized QC VM size for cost savings
4. snpeff now uses data from the refdata dir instead of downloading
5. added test datasets


### v0.5.7

#### Changes:
1. updated biowrappers version to v0.2.2
2. updated destruct and remixt container versions
3. added sv genotyping (experimental)

#### bugs:
1. missing last line in lumpy parsed output
2. issues with fastqscreen tags 


### v0.5.6:

#### bugs:
1. fixed an issue with double headers in snv calling output
2. issues in plotting due to gc data dtypes set to string

#### Changes:
1. destruct container does not need single cell pipeline
2. remixt container does not need single cell pipeline 
3. some updates to documentation
4. updated to biowrappers v0.2.1
5. split infer and count haps

### v0.5.5:

#### changes:
1. updated pypeliner to v0.5.13
2. updated documentation

### v0.5.4:

#### bugs:
1. mismatching csv type error in empty fastq files

#### changes:
1. refactor errors with empty fastq screen files when fastq is empty


### v0.5.3:

#### bugs:
1. fixed base64 encode issue for images in qc reports
2. issue deleting non-existent bam key
3. missing .gz.csv.yaml extensions for lumpy files
4. fixed container for lumpy
5. pd.concat on empty set of dataframes fixed
6. plotting bugs in hmmcopy
7. pandas loading issues related to dtypes specified but not names
8. use new hmmcopy script, better params format

### v0.5.2:

#### Changes:
1. alignment for all lanes is run in a single job
2. optimized destruct fastq reindexing

#### bugs:
1. removed redundant replace ? by 0 in plotting
2. annotation works when sample_info is none
3. fixed paths for annotation low mad yaml in tests
4. updated dtypes in integrity tests

#### bugs:
1. issues with NA handling in csvutils
2. all median cols are float now

### v0.5.1:

#### bugs:
1. issues with NA handling in csvutils
2. all median cols are float now

### v0.5.0: 
CLI refactor

#### added:
1. germline calling mode

#### Changes:
1. refactored the CLI
2. removed: 
    a) QC and
    b) multi_sample_pseudo_bulk
3. added in commands for
    a) alignment,
    b) hmmcopy,
    c) annotation,
    d) merge_cell_bams,
    e) split_normal,
    f) variant_calling,
    g) variant_counting (multi-sample)
    h) germline_calling
    i) infer_haps
    j) breakpoint_calling
4. add predefined dtypes to all workflows
5. metadata yaml files are generated within the pipeline
6. added a sentinel file with some provenance information as a teardown job

#### removed:
1. aneufinder
2. copynumber_calling

### v0.4.2
#### added
1. Enforced data dtypes in csvutils
2. added input yaml to output and metadata.yaml
3. jenkins: added output integrity check

#### Changes
1. removed autoscale
2. fixes to travis builds, travis now builds with python3
3. merged alignment tasks
4. removed local indel realignment

### v0.4.1
#### Changes:
* variant calling: switch h5 output to csv
* containers: picard and hmmcopy containers now use base R container

#### bugs:
* fixed error in QC generation on low quality datasets.
* destruct needs more disk
* metadata: code didn't upload when storage is set to local.

### v0.4.0

#### Changes:
* pipeline now runs on python3. python2 is not supported anymore
* pandas call for converting to categoricals
* refactor config generation
* added docs for LSF and singularity
* updated type from alignment to align in metadata
* library level snv counting removed from variant calling

#### bugs: 
* check fastq screen output directory for files from older runs and delete them
* fixed error raised when uploading meta yaml where storage is not specified
* readcounter: can handle non tagged bams now
* fastqscreen: can handle fastq files with multiple periods in name

### v0.3.1

#### added
* added a column indicating if a cell is contaminated
* added a column indicating if a segment is low mappability
* filtered contaminated cells from heatmap
* added extensions to metadata yaml files
#### Changes:
* seqdata files from haplotype calling are now temporary files
* fastqscreen counts column names begin with fastqscreen_


### v0.3.0

##### added
* added fastq screen
	* runs fastqscreen with `--tag`
	* all downstream analysis is run on the tagged data.
	* bam headers contain required information for parsing fastq screen tag
	* by default, pipeline removes all reads that belong to another organism
	* generates a detailed table and adds summary metrics to alignment table 
	* more details at [organism filter](docs/description/organism_filter.md)
* added salmon reference to images.
* added conda package for corrupt tree. updated docker container to use the conda package
* added newick support to heatmap
* added cell order based on corrupt tree to output
* added this changelog
* added metadata yaml files to output directories
* added flag to disable corrupt tree
##### changes:
* hmmcopy segments plots have a global max for ylim per run (library)
* standardized page size for corrupt tree output, annotated each page.
* replaced yaml.load with yaml.safe_load
* replaced nan values in QC html with 0
* removed biobloom
* destruct can now handle empty/small fastq files.
* fixed strelka filename issue (missing _)
* refactor alignment workflow
* added a tarball output with all hmmcopy outputs except autoploidy (multipliers 1-6)
* merged all picard based metrics into a single tarball
* reorganized reference data
* now uses miniconda docker image to delete files in batch
* arguments changes for QC 
	* removed --out_dir
	* add --alignment_output
	* add --hmmcopy_output
	* add --annotation_output
* argument changes for pseudo bulk
    * removed --out_dir
    * added --variants_output
    * added --haps_output
    * added --destruct_output
    * added --lumpy_output
##### bugs:
* fixed missing header issue with destruct outputs
* pipeline can now handle tsv files.
* fixed issues with missing cell cycle data in outputs

### v0.2.25
##### added
* Added Corrupt Tree
##### changes
* Reorganized QC pipeline outputs
* updated to newest biobloom container (v0.0.2). biobloom container now runs as root user.
* QC html doesn't require reference GC curve data
* give more memory to biobloom
* load input yaml with safe_load
##### bugs
* bugfix: fixed a merge issue with trim galore running script.

### v0.2.24
##### added
* * Added Cell Cycle Classifier
##### changes
* disable biobloom by default

### v0.2.23
##### bugs
* bugfix: Destruct was not tagging reads with cell ids
##### changes
* removed reference fasta index from github

### v0.2.22
##### changes
* lumpy can now handle empty bams
### v0.2.21
##### added
* added Html QC output
* added biobloom
* added a single 'QC' command to run alignment and hmmcopy
##### changes
* merged the alignment and metrics workflows.
* remove hmmcopy multipliers, only use autoploidy downstream
* removed option to specify multiple hmmcopy parameter sets
##### bugs
* bugfix: issue with automatic dtype detection in csv yaml files.

### v0.2.20
##### changes
* updated input yaml format for pseudowgs. The normal section now follows same schema as tumour (with sample and library id).

### v0.2.19
##### bugs
* bug: missing header in allele_counts file.

### v0.2.18
##### added
* added travis build.
* added smarter dtype merging for csv files.
##### changes
* updated conda recipe
* updated destruct output format from h5 to csv
* fixed destruct to generate counts from filtered output to remove normal reads
* optimized breakpoint calling, normal preprocessing runs only once per run.
* cleaned up raw_dir in output folder
* updated to conda based hmmcopy and mutationseq containers
* updated to latest version of lumpy with correct bed output
* now supports multiple libraries per normal in pseudowgs

### v0.2.17
##### changes
* updated lumpy bed file parsing.
* changed lumpy output file format from h5 to csv.

### v0.2.16
##### changes
* added mutationseq parameters to config. now users can override default settings.

### v0.2.15
##### added
* added parallelization over libraries in pseudowgs
##### bugs
* fixed an issue with read tagging that caused int overflow in bowtie
* some pickling issues due to python compatibility updates in pypeliner
* fixes in csv and yaml generation code

### v0.2.14
##### changes
* refactored lumpy workflow, only run normal preprocessing once per run
* merges in destruct require more disk space
* destruct: read indexes are now unique int
* destuct: reindex both reads in a single job to reduce number of jobs
* destruct: prepocess normal once per run
* faster csv file concatenation
* updated batch config to match pypeliner v0.5.6. now pool selection also accounts for disk usage.
	* Each pool will have available disk space. jobs will be scheduled in a pool based on requirements. production will have smaller disk in standard pool to save on costs.
* classifier now supports csv inputs
##### bugs
* bwa couldnt parse readgroup when not running in docker

### v0.2.13
##### changes
* split and merge bams only when running snv calling
* refactored to make main workflow calling functions standalone subworkflows
* revamped destruct workflow for pseudobulk

### v0.2.12
##### changes
##### changes
* switched to gzipped csv from H5 due to compatibility issues
* order IGV segs file by the clustering order, filter on quality

### v0.2.11
##### added
* added flags to only run parts of pseudowgs workflow
##### bugs
* fixed issue with infer haps where some parameters werent specified correctly.

### v0.2.10
##### changes
* destruct and remixt containers now use the same versioning as single cell pipeline
* lumpy accepts normal cells
* refactor: haplotype calling workflows
* all psudo wgs commands use the same input format as multi sample pseudo bulk.
* bam merge now supports merging larger number of files.

### v0.2.9
##### added
* destruct now supports list of cells as normal
* separate pools based on disk sizes.

##### changes
* separate docker container for destruct
* haplotype calling supports list of cells as normal
* parallel runs support more cells now.
##### bugs
* bug: fixed an issue that caused low mappability mask to disappear in the heatmap

### v0.2.8
##### changes
* replaced python based multiprocessing with gnu parallel
##### bugs
* bug: remixt path fixed in config, mkdir doesnt cause failures anymore in batch vm startup

### v0.2.7
##### added
* added trim galore container
* added a flag to switch disk to 1TB for all batch nodes 
* added a flag to specify whether to trim the fastqs. The flag overrides the sequencer based trimming logic.
* row, column, cell_call and experimental condition can now be null
* switched to gnu parallel for parallel on node runs
* feature: pools are now chosen automatically
##### changes
* updated readgroup string.
* renamed total_mapped_reads column in hmmcopy to total_mapped_reads_hmmcopy to avoid clashed with column of same name in alignment metrics
* snv calling: allow overlaps in vcf files
* remove meta yaml file
* moved autoploidy segment plot to top of page
* added option to launch the pipeline with docker by just adding `--run_with_docker`
* h5 dtype casting uses less memory now
* alignment metrics plot: now faster, plots atmost 1000 cells per page. extra cells overflow onto to the next page.
* support for pypeliner auto detect batch pool
* updated docker 
* updated vcfutils to use pypeliner to handle vcf index files.
* subworkflow resolution now runs in a docker container on compute nodes.
* switched from warnings to logging. the logs from compute now gets reported in main pypeliner log file.
* pipeline now uses OS disk in azure to store temporary files
* updated docker configuration changes in pypeliner. the container doesnt require the prefix anymore.
* added info.yaml file with some metadata per run
* merged andrew's pseudobulk changes.
##### bugs
* fix: strelka uses chromosome size instead of genome size

### v0.2.6
##### changes
* now supports plain text fastq files

### v0.2.5
##### added
* added: heatmap and boxplots for cell quality score
##### changes
* switching to smaller 256GB disks

### v0.2.4
##### changes
* use `table` format in h5 files
##### bugs
* fix: non unique categories error fixed by specifying categories at initialization

### v0.2.3
##### changes
* properly deletes file on batch node after task completes
* cell_id is now a categorical in output
* supports `.fq` and `.fq.gz` fastq file extensions
* heatmap is generated even if all cells are `nan`
* default for non-azure environments is not singularity anymore.

### v0.2.2
##### added
* added docker container info to info yaml files
* added info yaml files 
##### changes
* interprets `?` in picard tools output as 0
* casts all columns in h5 to their correct dtypes.

### v0.2.1
##### added
* added multisample pseudobulk
##### changes
* divided alignment into 2 separate workfloes
* replaced segments and bias pdf with per cell plots with a tarball of png files.

### v0.2.0
##### added
* added singularity support
* added docker support for whole genome
* added docker support for alignment and hmmcopy
* added test data set
##### changes
* project wide refactor
* switched to image with a single disk. removed startup mount commands.
* pseudowgs: added infoer Haplotype code
* clip copynumber in hmmcopy plots to 40

### v0.1.5
##### added
* added LTM
* added haploid poison to hmmcopy
##### changes
* yaml files use block style format
* added cell quality classifier
* bwa-mem is now the default aligner. bwa-aln is also supported
##### bugs
* fix: handles empty segments in hmmcopy segment plot

### v0.1.4
##### added
* pseudowgs: added option to merge bams 
* pseudowgs: added option to split bam by reads (pairs next to each other)
##### changes
* smaller segments plots file size, consistent colormap across plots
* added mask for low mappability regions in heatmap
* one pdf file for segments and bias plots per row of cells 

### v0.1.3
##### added
* added VM image URI and SKU to batch yaml file
* hmmcopy: added autoploidy 
* added titan to pseudowgs
* hmmcopy can now run independently from alignment
##### changes
* rename pick_met to cell_call and condition to experimental_condition
* chooses VM image based on the pipeline version
* cleanly exit hmmcopy script if data is not enough/missing
* aneufinder,alignment, hmmcopy output is now h5
##### bugs
* fixed reouding issue in autoscale formula (there is no round method)

### v0.1.2
##### changes
* streka now runs over split bam files
* can now save split bams using a template specified at run time

### v0.1.1
##### changes
* switch to 1-based state in hmmcopy

### v0.1.0
##### added
* added filtering to copynumber heatmap
* added classifier
* each lane now has a sequencing centre 
##### changes
* metrics heatmaps are not restricted to 72*72
* exposed all hmmcopy params to config file
* modal correction can run on empty datasets

