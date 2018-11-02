'''
Created on November 01, 2018

@author: pashaa
'''
from scripts import DemultiplexBam
import pysam
import logging
from logging.config import fileConfig


def demultiplex(ibam, obams, barcode_csv):
    fileConfig("logging_config.ini")
    logging.info('Demultiplexing bams...')
    with DemultiplexBam(ibam, obams, barcode_csv) as demux:
        demux.demultiplex_bam()

def collate(ifile, ofile):
    """
    Collate bam file

    :param ifile: input bam
    :param ofile: output collated bam
    """

    logging.info('collating bam file ' + str(ifile))

    pysam.collate("-o", ofile, ifile)

def fastq(ifile, ofile1, ofile2):
    """
    :param ifile: input bam file
    :param ofile1: first file of fastq pair
    :param ofile2: second file of fastq pair
    """

    logging.info('generating fastq file for ' + str(ifile))

    pysam.fastq('-1', ofile1,
                '-2', ofile2,
                '-n',
                ifile)
