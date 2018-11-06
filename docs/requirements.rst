============
Requirements
============

Demultiplex BAM files:
-----------------------
**Functional**:

1. Demultiplex bam files based on the optional CB (cell identifier) tag.

2. User shall specify list of cell ids/barcodes on which to demultiplex.

3. Aligned reads associated with barcodes not contained in user specified list must be filtered out.

4. Aligned reads that don't have CB tags must be assigned the special tag="undetermined".

5. Reconstruct the paired reads for each originating cell and output paired reads to _1 and _2 FASTQ files.


**Non-functional**:

1. Use parallelization where possible.

2. TODO: Application level logging.

3. TODO: Unit tests.

4. Cleanup any temp files in post processing step.

5. TODO: Output stats obtained from processing to log and console out.
