# QC pipeline metrics

|Column|Description|
|------|-----------|
|cell_id|label of the cell|
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
