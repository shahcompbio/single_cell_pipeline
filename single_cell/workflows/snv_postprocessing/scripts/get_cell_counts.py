'''
Created on Jul 14, 2017

@author: dgrewal
'''
from collections import Counter

import csv
import pysam
import gzip



class GetCounts(object):
    def __init__(self, bam, positions, output, sample_id, **kwargs):
        self.sample_id = sample_id
        self.bam = bam
        self.positions = positions
        self.output = output

        self.max_pileup_depth = kwargs.get('max_pileup_depth')
        self.min_bqual = kwargs.get('min_bqual')
        self.min_mqual = kwargs.get('min_mqual')
        self.count_duplicates = kwargs.get('count_duplicates')
        self.count_qc_failures = kwargs.get('count_qc_failures')
        self.strand_counts = kwargs.get('strand_counts')
        self.min_counts = kwargs.get('min_counts')

        self.max_pileup_depth = self.max_pileup_depth if self.max_pileup_depth else int(1e7)
        self.min_bqual = self.min_bqual if self.min_bqual else 0
        self.min_mqual = self.min_mqual if self.min_mqual else 0
        self.min_counts = self.min_counts if self.min_counts else 0

    
    def main(self):
        sampid = self.sample_id
        
        bam_file = pysam.Samfile(self.bam, 'rb')
    
        bases, fields = self.get_bases_and_fields()
    
        writer = csv.DictWriter(open(self.output, 'w'), fields, delimiter=',')
    
        writer.writeheader()
    
        positions, target_bases = self.load_positions(self.positions)
    
        pileup_iterator = self.positions_iterator(
            bam_file, positions, self.max_pileup_depth)
    
        for i,pileup_column in enumerate(pileup_iterator):

            if not pileup_column:
                continue

            counts = self.get_counts(pileup_column,
                                self.min_bqual,
                                self.min_mqual,
                                self.count_duplicates,
                                self.count_qc_failures,
                                self.strand_counts)

            total_counts = sum([counts[b] for b in bases])

            if total_counts < self.min_counts:
                continue

            out_row = {}
    
            out_row['sample_id'] = self.sample_id
            out_row['chrom'] = bam_file.getrname(pileup_column.tid)
    
            # One based coordinate
            out_row['coord'] = pileup_column.pos + 1
    
            pos = (out_row['chrom'], out_row['coord'])
    
            ref_base = target_bases[pos]['ref_base']
    
            var_base = target_bases[pos]['var_base']
    
            out_row['ref_base'] = ref_base
    
            out_row['var_base'] = var_base
    
            if self.strand_counts:
                out_row['ref_counts_forward'] = counts[ref_base.upper()]
    
                out_row['ref_counts_reverse'] = counts[ref_base.lower()]
    
                out_row['var_counts_forward'] = counts[var_base.upper()]
    
                out_row['var_counts_reverse'] = counts[var_base.lower()]
    
            else:
                out_row['ref_counts'] = counts[ref_base]
    
                out_row['var_counts'] = counts[var_base]
    
            writer.writerow(out_row)
    
    
    def get_bases_and_fields(self):
        if self.strand_counts:
            bases = ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't']
        else:
            bases = ['A', 'C', 'G', 'T']
    
    
        fields = ['chrom', 'coord', 'ref_base', 'var_base']
    
        if self.strand_counts:
            fields += ["sample_id", 'ref_counts_forward', 'ref_counts_reverse',
                       'var_counts_forward', 'var_counts_reverse']
        else:
            fields += ["sample_id", 'ref_counts', 'var_counts']
    
        return bases, fields
    
    
    def get_counts(self, pileup_column, min_bqual, min_mqual, count_duplicates, count_qc_failures, strand_counts):
        bases = []
    
        for pileup_read in pileup_column.pileups:
            # Skip flagged duplicates if we don't want them
            if pileup_read.alignment.is_duplicate and not count_duplicates:
                continue
    
            # Skip QC failures if we don't want them
            if pileup_read.alignment.is_qcfail and not count_qc_failures:
                continue
    
            if pileup_read.is_del:
                continue
    
            mqual = pileup_read.alignment.mapq
    
            # Skip positions with low mapping quality
            if mqual < min_mqual:
                continue
    
            # Nucleotide handling
            else:
                bqual = ord(pileup_read.alignment.qual[pileup_read.query_position]) - 33
    
                if bqual < min_bqual:
                    continue
    
                base = pileup_read.alignment.seq[pileup_read.query_position]
    
                if strand_counts:
                    if pileup_read.alignment.is_reverse:
                        base = base.lower()
    
                    else:
                        base = base.upper()
    
                else:
                    base = base.upper()
    
                bases.append(base)
    
        return Counter(bases)
    
    
    def is_gzip(self, filename):
        """
        Uses the file contents to check if the file is gzip or not.
        The magic number for gzip is 1f 8b
        See KRONOS-8 for details
        """
        with open(filename) as f:
            file_start = f.read(4)
        
            if file_start.startswith("\x1f\x8b\x08"):
                return True
            return False
    
    def get_var_base(self, counts, ref_base):
        # Get rid of strand information
        counts = Counter([x.upper() for x in counts.elements()])
    
        # Remove reference base
        del counts[ref_base]
    
        if len(counts) == 0:
            var_base = 'N'
    
        else:
            # Get most common non-reference base
            var_base, var_counts = counts.most_common()[0]
    
        return var_base
    
    
    def load_positions(self, file_name):
    
        positions = []
        target_bases = {}
    
    
        if self.is_gzip(file_name):
            merge_file = gzip.open(file_name, 'rb')
        else:
            merge_file = open(file_name)
    
        header = merge_file.readline().strip().split(",")
        chr_idx = header.index('chromosome')
        pos_idx = header.index('start')
        ref_idx = header.index('ref')
        alt_idx = header.index('alt')
    
    
    
        for line in merge_file:
            columns = line.split(",")
            chrom = columns[chr_idx]
            pos = int(columns[pos_idx])
            ref = columns[ref_idx]
            alt  = columns[alt_idx]
    
            pos = (chrom, pos)
            target_bases[pos] = {'ref_base': ref, 'var_base': alt}
            
            positions.append(pos)
        
        merge_file.close()
    
        return positions, target_bases
    
    
    
    def positions_iterator(self, bam_file, positions, max_pileup_depth):
        for chrom, coord in positions:
            if bam_file.count(chrom, coord - 1, coord) == 0:
                pos = coord - 1
    
                tid = bam_file.gettid(chrom)
    
                pileup_column = DummyPileupProxy(pos, tid)
    
            else:
                pileup_iterator = bam_file.pileup(chrom,
                                                  coord - 1,
                                                  coord,
                                                  truncate=True,
                                                  max_depth=max_pileup_depth,
                                                  mask=False,
                                                  stepper='all')

                pileup_column = None
                for pileup_column in pileup_iterator:
                    break
    
            yield pileup_column


class DummyPileupProxy(object):

    def __init__(self, pos, tid):
        self.pileups = []

        self.pos = pos

        self.tid = tid

