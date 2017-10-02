'''
Created on Sep 29, 2017

@author: svatrt
'''
import csv

class ConvertCSVToSEG(object):
	def __init__(self, filtered_segs, filtered_reads, output_seg):
		self.filtered_segs = open(filtered_segs, 'r')
		self.filtered_reads = open(filtered_reads, 'r')
		self.output_seg = open(output_seg, 'w')

	def main(self):
		# Get bin_width
		reads = csv.reader(self.filtered_reads)
		bin_width = [row for row_num, row in enumerate(reads) if row_num == 1][0][-1]
		self.filtered_reads.close()
		
		# Write the first line of the output file
		self.output_seg.write('\'ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n')

		# Skip the first line of the segs file
		segs = csv.reader(self.filtered_segs)
		lines = enumerate(segs)
		next(lines)
		
		# Read the segs file and write to the output file
		for row_num, row in lines:
			segment_length = int(row[2]) - int(row[1]) + 1
			num_mark = segment_length / int(bin_width)
			self.output_seg.write(row[-1] + '\t' + row[0] + '\t' + row[1] + '\t' + row[2] + '\t' + str(num_mark) + '\t' + row[3] + '\n')
		
		self.filtered_segs.close()
		self.output_seg.close() 

# Testing
# test = ConvertCSVToSEG('filtered_segs.csv', 'filtered_reads.csv', 'output.seg')
# test.main()