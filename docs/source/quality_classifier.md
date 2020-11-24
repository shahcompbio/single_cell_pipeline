# Cell Quality Classifier


|Feature name|Source|Description |
| ----| ----|----|
|percent_duplicate_reads|picard|percentage of reads marked as PCR duplicate by MarkDuplicates|
|total_mapped_reads | samtools|number of reads mapped by the bwa mem alignment algorithm |
|total_duplicate_reads | samtools|number of reads marked as PCR duplicate by MarkDuplicates |
|standard_deviation_insert_size| picard| read insert size standard deviation |
|MSRSI_non_integerness| hmmcopy| median of segment residuals from segment integer copy number states|
|MBRSI_dispersion_non_integerness| hmmcopy| median of bin residuals from segment integer copy number states|
|MBRSM_dispersion| hmmcopy | median of bin residuals from segment median copy number values|
|autocorrelation_hmmcopy| | autocorrelation of CNV results|
|cv_hmmcopy| hmmcopy| coefficient of variation of CNV results|
|mad_hmmcopy| hmmcopy| mean absolute deviation of CNV results|
|total_halfiness|hmmcopy | halfiness score but without copy number state scaling|
|scaled_halfiness| hmmcopy| a scaled metric to assess integer goodness of fit, described in text|
|mean_state_mads| hmmcopy| the mean across all MADs of each copy number state|
|mean_state_vars| hmmcopy| the mean across all variances of each copy number state|
|breakpoints| hmmcopy| number of intrachromosomal breakpoints|
|mean_copy| hmmcopy| mean copy number of all genomic bin segments|
|state_mode| hmmcopy| the most commonly occuring copy nubmer state|
|log_likelihood| hmmcopy| log-likelihood of HMMcopy CNV fit|


## Percent Duplicate Reads

 Calculated from the output of Mark Duplicates from picard tools. Please see [mark duplicates](http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics)


Formula:

UNPAIRED_READ_DUPLICATES + ((READ_PAIR_DUPLICATES + READ_PAIR_OPTICAL_DUPLICATES)*2) / (UNPAIRED_READS_EXAMINED  + (READ_PAIRS_EXAMINED * 2)) 

