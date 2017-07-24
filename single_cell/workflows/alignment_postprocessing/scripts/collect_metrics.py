'''
Extract metrics table.
'''

from __future__ import division
import pandas as pd
import os

class CollectMetrics(object):
    def __init__(self, wgs_metrics, insert_metrics, flagstat_metrics, markdups_metrics, output, samplesheet, sample_id):
        self.wgs_metrics = wgs_metrics
        self.flagstat_metrics = flagstat_metrics
        self.insert_metrics = insert_metrics
        self.markdups_metrics = markdups_metrics
        self.output = output
        self.samplesheet = samplesheet
        self.sample_id = sample_id
    

    def extract_wgs_metrics(self):
        """
        get the coverage_depth (mean_coverage column)
        get the coverage_breadth (count/genome_territory)
        """
    
        mfile = open(self.wgs_metrics)
        
        metrics = []
        hist = {}
        
        addmetrics = False
        addhist = False
       
       
        for line in mfile:
            if line.strip() == '':
                continue
            if line.startswith('## METRICS CLASS'):
                addmetrics = True
                addhist = False
                continue
    
            if line.startswith('## HISTOGRAM'):
                addhist = True
                addmetrics = False
                continue
            
            if addmetrics:
                metrics.append(line.strip().split('\t'))
            if addhist:
                line = line.strip().split('\t')
                if line[0] == 'coverage':
                    continue
                hist[int(line[0])] = int(line[1])
    
        mfile.close()
        header, data = metrics
    
        header = [v.lower() for v in header]
        header = {v:i for i,v in enumerate(header)}
    
        gen_territory = int(data[header['genome_territory']])
        cov_depth  = float(data[header['mean_coverage']])
        count = int(hist[0])
        cov_breadth = (gen_territory - count) / gen_territory 
    
        return cov_breadth, cov_depth
    
    
    def extract_flagstat_metrics(self):
        """
        extract from flagstat
        """
    
        df = pd.read_table(self.flagstat_metrics,
                           sep=r'\s\+\s0\s',
                           header=None,
                           names=['value', 'type'],
                           engine='python')
    
        tot_reads = df[df['type']=='in total (QC-passed reads + QC-failed reads)']['value']
        tot_mpd_reads = df[(df['type'].str.contains('mapped') == True ) & ( df['type'].str.contains('mate mapped') == False)]
        tot_dup_reads = df[df['type']=='duplicates']['value']
        tot_prop_paired = df[df['type'].str.contains('properly paired') ]
    
        assert len(tot_reads) == 1
        assert len(tot_mpd_reads) == 1
        assert len(tot_dup_reads) == 1
        assert len(tot_prop_paired) == 1
    
        tot_reads = tot_reads.iloc[0]
        tot_mpd_reads = tot_mpd_reads['value'].iloc[0]
        tot_dup_reads = tot_dup_reads.iloc[0]
        tot_prop_paired = tot_prop_paired['value'].iloc[0]
    
    
        return tot_reads, tot_mpd_reads, tot_dup_reads, tot_prop_paired
    
    def extract_duplication_metrics(self):
        """
        extract from markdups
        """
    
        mfile = open(self.markdups_metrics)
    
        targetlines = []
    
        line = mfile.readline()
    
        while line != '':
            if line.startswith('## METRICS CLASS'):
                targetlines.append(mfile.readline().strip('\n').split('\t'))
                targetlines.append(mfile.readline().strip('\n').split('\t'))
                break
            line = mfile.readline()
    
        mfile.close()
    
        header, data = targetlines
    
        header = [v.lower() for v in header]
        header = {v:i for i,v in enumerate(header)}
    
        unprd_mpd_rds = int(data[header['unpaired_reads_examined']])
        prd_mpd_rds = int(data[header['read_pairs_examined']])
        unprd_dup_rds = int(data[header['unpaired_read_duplicates']])
        prd_dup_rds = int(data[header['read_pair_duplicates']])
        unmpd_rds = data[header['unmapped_reads']]
        est_lib_size = data[header['estimated_library_size']]
    
        rd_pair_opt_dup = int(data[header['read_pair_optical_duplicates']])
    
        try:
            perc_dup_reads = (unprd_dup_rds + ((prd_dup_rds + rd_pair_opt_dup) * 2)) / (unprd_mpd_rds + (prd_mpd_rds * 2))
        except ZeroDivisionError:
            perc_dup_reads = 0
    
        outdata = (unprd_mpd_rds, prd_mpd_rds, unprd_dup_rds, prd_dup_rds,
                   unmpd_rds, perc_dup_reads, est_lib_size)
    
        outdata = tuple(['nan' if val == '' else val for val in outdata])
        return outdata
    
    def extract_insert_metrics(self):
        ''' Extract median and mean insert size '''
    
        # picardtools insertmetrics completes with code 0 and doesn't generate metrics file
        # if inputs don't have sufficient read count
        if not os.path.isfile(self.insert_metrics):
            return 0, 0, 0
    
        mfile = open(self.insert_metrics)
        
        targetlines = []
        
        line = mfile.readline()
        
        while line != '':
            if line.startswith('## METRICS CLASS'):
                targetlines.append(mfile.readline().strip().split('\t'))
                targetlines.append(mfile.readline().strip().split('\t'))
                break
            line = mfile.readline()
    
        mfile.close()
    
        header, data = targetlines
    
        header = [v.lower() for v in header]
        header = {v:i for i,v in enumerate(header)}
    
        median_ins_size = data[header['median_insert_size']]
        mean_ins_size = data[header['mean_insert_size']]
        std_dev_ins_size = data[header['standard_deviation']]
    
        return median_ins_size, mean_ins_size, std_dev_ins_size
    
    
    def write_data(self, header, data):
        """
        write to the output
        """
        assert len(header) == len(data)
        #replace empty vals with NA
        data = [v if v != '' else 'NA' for v in data]
    
        writer = open(self.output, 'w')
        writer.write(','.join(header) + '\n')
        writer.write(','.join([str(v) for v in data]))
        writer.close()
    
    
    def extract_sample_info(self):
        """
        get info
        """
        header = True
    
        with open(self.samplesheet) as sample_info:
            for line in sample_info:
                line = line.strip().split(',')
        
                if header:
                    if line[0] == "[Data]":
                        header = False
                else:
                    if line[0] in ['Sample_ID', 'Sample-ID']:
                        continue
                    sampid = line[0]
                    
                    if sampid.replace('-','_') != self.sample_id.replace('-','_'):
                        continue
                    plate = line[2]
                    well = line[3]
                    desc = line[9]
                    i7 = line[6]
                    i5 = line[8]
    
    
        if ';' in desc:
            desc = desc.split(';')
            desc = [val.split('=') for val in desc]
    
            cell_call = [val[1] for val in desc if val[0] == 'CC']
            cell_call = cell_call[0] if cell_call else 'NA'
    
            exp_cond = [val[1] for val in desc if val[0] == 'EC']
            exp_cond = exp_cond[0] if exp_cond else 'NA'
    
    
            samp_typ = [val[1] for val in desc if val[0] == 'ST']
            samp_typ = samp_typ[0] if samp_typ else 'NA'
        elif desc == '':
            cell_call = 'C1'
            exp_cond = 'A'
            samp_typ = '1'
        else:
            cell_call = desc
            exp_cond = 'NA'
            samp_typ = 'NA'
    
        well = well if well != '' else 'R1_C1'
        plate = plate if plate !='' else 'R1-C1'
        i5 = i5 if i5!='' else 'i5-1'
        i7 = i7  if i7 != '' else 'i7-1'
    
    
        return cell_call, exp_cond, samp_typ, well, plate, i5, i7
    
    #=========================================================================
    # Run script
    #=========================================================================
    
    def main(self):
    
        duplication_metrics = self.extract_duplication_metrics()
        flagstat_metrics = self.extract_flagstat_metrics()
        wgs_metrics = self.extract_wgs_metrics()
    
        header = ['cell_id', 'cell_call', 'experimental_condition', 'sample_type', 'sample_well',
                  'sample_plate', 'i5_barcode', 'i7_barcode', 'unpaired_mapped_reads',
                  'paired_mapped_reads', 'unpaired_duplicate_reads',
                  'paired_duplicate_reads', 'unmapped_reads', 'percent_duplicate_reads',
                  'estimated_library_size', 'total_reads', 'total_mapped_reads',
                  'total_duplicate_reads', 'total_properly_paired',
                  'coverage_breadth', 'coverage_depth']
    
        if self.sample_id:
            sample_id = self.sample_id
            cell_call, exp_cond, samp_typ, samp_well, samp_plate, i5, i7 = self.extract_sample_info()
        else:
            sample_id = cell_call = exp_cond = samp_typ = samp_well = samp_plate = i5 = i7 = 'NA'
    
        output = (sample_id, cell_call, exp_cond, samp_typ, samp_well, samp_plate, i5, i7) + duplication_metrics + flagstat_metrics + wgs_metrics
    
        if self.insert_metrics:
            insert_metrics = self.extract_insert_metrics()
            output += insert_metrics
            header += ['median_insert_size',
                       'mean_insert_size',
                       'standard_deviation_insert_size']
    
        self.write_data(header, output)

