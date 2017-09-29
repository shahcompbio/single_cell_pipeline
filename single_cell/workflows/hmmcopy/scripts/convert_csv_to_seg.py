'''
Created on Sep 29, 2017

@author: svatrt
'''

class ConvertCSVToSEG(object):
	def __init__(self, filtered_segs, filtered_reads, output_seg):
		self.filtered_segs = filtered_segs
		self.filtered_reads = filtered_reads
		self.output_seg = output_seg

	def main(self):
		pass

x = ConvertCSVToSEG(1, 2, 3)
print x.filtered_reads