'''
Created on Feb 20, 2018

@author: dgrewal
'''


from PyPDF2 import PdfFileMerger
import os
        
def merge_pdfs(infiles, outfile):

    merger = PdfFileMerger()

    for infile in infiles:
        merger.append(open(infile, 'rb'))

    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))

    with open(outfile, 'wb') as fout:
        merger.write(fout)
