# QC pipeline metrics

|Column|Description|
|------|-----------|
|multiplier|during parameter searching, the set [1..6] that was chosen|
|MSRSI_non_integerness|median of segment residuals from segment integer copy number states|
|MBRSI_dispersion_non_integerness|median of bin residuals from segment integer copy number states|
|MBRSM_dispersion|median of bin residuals from segment median copy number values|
|autocorrelation_hmmcopy|hmmcopy copy autocorrelation|
|cv_hmmcopy||
|empty_bins_hmmcopy|number of empty bins in hmmcopy|
|mad_hmmcopy|median absolute deviation of hmmcopy copy|
|mean_hmmcopy_reads_per_bin|mean reads per hmmcopy bin|
|median_hmmcopy_reads_per_bin|median reads per hmmcopy bin|
|std_hmmcopy_reads_per_bin|standard deviation value of reads in hmmcopy bins|
|total_mapped_reads_hmmcopy|total mapped reads in all hmmcopy bins|
|total_halfiness|summed halfiness penality score of the cell|
|scaled_halfiness|summed scaled halfiness penalty score of the cell|
|mean_state_mads|mean value for all median absolute deviation scores for each state|
|mean_state_vars|variance value for all median absolute deviation scores for each state|
|mad_neutral_state|median absolute deviation score of the neutral 2 copy state|
|breakpoints|number of breakpoints, as indicated by state changes not at the ends of chromosomes|
|mean_copy|mean hmmcopy copy value|
|state_mode|the most commonly occuring state|
|log_likelihood|hmmcopy log likelihood for the cell|
|true_multiplier|the exact decimal value used to scale the copy number for segmentation|
|cell_id|label of the cell|
|order|order of the cell in the hierarchical clustering tree|
|index_sequence|index sequence of the adaptor sequence|
|column|column of the cell on the nanowell chip|
|img_col|column of the cell from the perspective of the microscope|
|index_i5|id of the i5 index adapter sequence|
|sample_type|type of the sample|
|primer_i7|id of the i5 index primer sequence|
|experimental_condition|experimental treatment of the cell, includes controls|
|index_i7|id of the i7 index adapter sequence|
|cell_call|living/dead classification of the cell based on staining usually, C1 == living, C2 == dead|
|sample_id|name of the sample|
|primer_i5|id of the i5 index primer sequence|
|row|row of the cell on the nanowell chip|
|quality|random forest classifier proability score that cell is good|
|estimated_library_size|scaled total number of mapped reads|
|total_mapped_reads|total number of mapped reads|
|nohit|number of reads with no organism match|
|salmon_multihit|number of reads that were classified as salmon and something else|
|total_duplicate_reads|number of duplicate reads|
|percent_duplicate_reads|percentage of duplicate reads|
|total_properly_paired|number of properly paired reads|
|mean_insert_size|mean insert size between paired reads|
|coverage_breadth|percentage of genome covered by some read|
|grch37|number of reads that were classified as human|
|unpaired_duplicate_reads|number of unpaired duplicated reads|
|unpaired_mapped_reads|number of unpaired mapped reads|
|unmapped_reads|number of unmapped reads|
|coverage_depth|average reads per nucleotide position in the genome|
|median_insert_size|median insert size between paired reads|
|salmon|number of reads that were classified as salmon|
|grch37_multihit|number of reads that were classified as human and something else|
|mm10|number of reads that were classified as mouse|
|total_reads|total number of reads, regardless of mapping status|
|standard_deviation_insert_size|standard deviation of the insert size between paired reads|
|paired_mapped_reads|number of mapped reads that were properly paired|
|mm10_multihit|number of reads classified as mouse and something else|
|paired_duplicate_reads|number of paired reads that were also marked as duplicate|
|is_contaminated|boolean, set to True if most reads belong to a different genome|

# HMMCopy Reads

|Column|Description|
|------|-----------|
|chr|chromosome|
|start|start position|
|end|end position|
|width|width of genomie segment that comprises the bin|
|reads|number of reads that start in the bin|
|gc|average GC content of all bases in the bin, -1 if N is present|
|map|average mappability value of bin|
|cor_gc|gc-corrected copy number value|
|copy|final output copy number value|
|valid|TRUE if reads > 0 & gc > 0, else FALSE|
|ideal|TRUE if bin is VALID with good mappability and non-outlier gc and read values|
|modal_curve|value of the gc-correction modal curve given the bin's gc|
|modal_quantile||
|cor_map|mappability-corrected gc-corrected copy number value|
|multiplier|hmmcopy parameter set used [1..6]|
|state|the copy number state of the bin|
|cell_id|label of the cell|
|is_low_mappability|bool, set to True if the segment has a low mappability score|

# HMMCopy Segments
|Column|Description|
|------|-----------|
|chr|chromosome|
|start|start position|
|end|end position|
|state|copy number state|
|median|median copy number value of segment|
|multiplier|hmmcopy parameter set used [1..6]|
|cell_id|label of the cell|
