# hmmcopy metrics

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
|is_contaminated|boolean, set to True if most reads belong to a different genome|
