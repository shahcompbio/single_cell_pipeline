'''
Created on Sep 29, 2017

@author: svatrt
'''
import csv
import os

class ConvertCSVToSEG(object):
	def __init__(self, filtered_segs, filtered_reads, output_seg):
		self.filtered_segs = filtered_segs
		self.filtered_reads = filtered_reads
		self.output_seg = open(output_seg, 'w')

	def main(self):

		# Check if filtered_reads file is empty
		if os.stat(self.filtered_reads).st_size != 0:
			self.filtered_reads = open(self.filtered_reads, 'r')
			reads = csv.reader(self.filtered_reads)
			
			# Get bin_width
			bin_width = [row for row_num, row in enumerate(reads) if row_num == 1][0][-1]
			self.filtered_reads.close()
		
			# Write the first line of the output file
			self.output_seg.write('\'ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n')
		
		# Check if filtered_segs file is empty
		if os.stat(self.filtered_segs).st_size != 0:
			self.filtered_segs = open(self.filtered_segs, 'r')
			segs = csv.reader(self.filtered_segs)
			
			# Skip the first line of the segs file
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
