'''
Created on Feb 20, 2018

@author: dgrewal
'''


from PyPDF2 import PdfFileMerger
import helpers
import os

def merge_pdfs(infiles, outfile):

    if isinstance(infiles, dict):
        infiles = infiles.values()

    merger = PdfFileMerger()

    for infile in infiles:
        #add it to list if not empty. skip empty files to avoid errors later
        if os.path.getsize(infile):
            merger.append(open(infile, 'rb'))

    helpers.makedirs(outfile, isfile=True)

    with open(outfile, 'wb') as fout:
        merger.write(fout)
