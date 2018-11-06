A reduced dataset for testing demultiplex_bam functionality. It should get processed in less than a minute.  

bj_mkn45_10pct_possorted_bam_10k_snippet.bam/bai:

Head 10,000 lines of bj_mkn45_10pct_possorted_bam.bam from 10X Genomics - https://support.10xgenomics.com/single-cell-dna/datasets/1.0.0/bj_mkn45_10pct

bj_mkn45_10pct_per_cell_summary_metrics.csv:

A csv file containing barcodes (cell identifiers) in the first column. The bam file is demultiplexed using these barcodes. This file is also from the 10x Genomics dataset and is truncated after the first 10 barcodes of the original file.  