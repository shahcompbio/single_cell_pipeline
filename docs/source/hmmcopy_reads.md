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
