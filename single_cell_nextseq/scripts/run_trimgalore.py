'''
Created on Jan 4, 2017

@author: dgrewal
'''

import argparse
from subprocess import Popen
import os
import warnings

class RunTrimGalore(object):
    """
    Trim fastq files with trimgalore
    """

    def __init__(self, args):
        self.args = args

        self.check_inputs()

        if not os.path.exists(self.args.tempdir):
            os.makedirs(self.args.tempdir)


    @classmethod
    def run_cmd(cls, cmd):
        """
        run a command with subprocess and return outputs
        """
        import sys

        proc = Popen(cmd, stdout=sys.stdout, stderr=sys.stderr)

        proc.communicate()

    @staticmethod
    def check_file(path):
        """
        check
        * if they exist or not,
        * are they too small
        * are they gzipped
        """
        if not os.path.exists(path):
            raise IOError("Couldn't find file: %s" % path)

        if os.path.getsize(path) < 1000:
            warnings.warn("extremely small file detected: %s" % path)

        ext = os.path.splitext(path)[1][1:]
        # get the last extension
        ext = ext[-1]

        return ext

    def check_inputs(self):
        """
        check input files
        * if they exist or not,
        * are they too small
        * are they gzipped
        """

        ext1 = self.check_file(self.args.seq1)
        ext2 = self.check_file(self.args.seq2)

        assert ext1 == ext2, "input fastqfiles should have the same extension"

    def run_trimgalore(self):
        """
        launch trimgalore
        """
        cmd = [self.args.trimgalore_path,
               '--fastqc',
               '--paired',
               '--path_to_cutadapt', self.args.cutadapt_path,
               '--output_dir', self.args.tempdir + '/',
               ]

        if self.args.adapter:
            cmd.extend(['--adapter', self.args.adapter])

        if self.args.adapter2:
            cmd.extend(['--adapter2', self.args.adapter2])

        cmd.extend([self.args.seq1, self.args.seq2, ])

        self.run_cmd(cmd)

    def move_files(self, fname, outpath):
        """
        move files from from temp dir to the expected path
        """
        path = os.path.join(self.args.tempdir, fname)
        os.rename(path, outpath)
        assert os.path.isfile(outpath)

        print outpath
        print os.stat(outpath)


    def get_file(self, r1_out, r2_out, ext, trimreport=False):
        """
        find the file in outdir and rename to r1_out and r2_out
        """
        outfiles = os.listdir(self.args.tempdir)

        reps = [v for v in outfiles if ext in v]
        assert reps != [], "Couldn't move %s files" % ext

        seq1 = os.path.basename(self.args.seq1)
        seq2 = os.path.basename(self.args.seq2)

        if trimreport:
            for rep in reps:
                if seq1 in rep:
                    self.move_files(rep, r1_out)
                elif seq2 in rep:
                    self.move_files(rep, r2_out)
                else:
                    raise Exception("Couldn't move %s files" % ext)
            return

        for rep in reps:
            if "val_1" in rep:
                self.move_files(rep, r1_out)
            elif "val_2" in rep:
                self.move_files(rep, r2_out)
            else:
                raise Exception("Couldn't move %s files" % ext)

    def gather_outputs(self):
        """
        rename the output metrics files.
        """
        # trimming reports
        self.get_file(
            self.args.report_r1,
            self.args.report_r2,
            "_trimming_report.txt",
            trimreport=True)

        # fastqc html files
        self.get_file(
            self.args.fastqc_report_r1,
            self.args.fastqc_report_r2,
            ".html")

        # fastqc zip files
        self.get_file(self.args.fastqc_zip_r1, self.args.fastqc_zip_r2, ".zip")

        # trimmed fastq files
        self.get_file(self.args.fastq_r1, self.args.fastq_r2, ".fq.gz")


def parse_args():
    """
    parse cmd line params
    """
    #=========================================================================
    # make a UI
    #=========================================================================
    parser = argparse.ArgumentParser(prog='run_trim_galore',
                                     description="""This script runs trim_galore with options --paired --nextera""")

    parser.add_argument('seq1',
                        help='FASTQ file for read 1')

    parser.add_argument('seq2',
                        help='FASTQ file for read 2')

    parser.add_argument('fastq_r1',
                        help="Path to output file (R1)")

    parser.add_argument('fastq_r2',
                        help="Path to output file (R2)")

    parser.add_argument('fastqc_report_r1',
                        help="Path to fastqc_report (R1)")

    parser.add_argument('fastqc_report_r2',
                        help="Path to fastqc_report (R2)")

    parser.add_argument('fastqc_zip_r1',
                        help="Path to fastqc_zip (R1)")

    parser.add_argument('fastqc_zip_r2',
                        help="Path to fastqc_zip (R2)")

    parser.add_argument('report_r1',
                        help="Path to report (R1)")

    parser.add_argument('report_r2',
                        help="Path to report (R2)")

    parser.add_argument('tempdir',
                        help="Path to output directory (for metrics)")

    parser.add_argument('--adapter',
                        help="specify adapter for R1")

    parser.add_argument('--adapter2',
                        help="specify adapter for R2")

    parser.add_argument('--trimgalore_path',
                        default='trimgalore',
                        help="specify path to trimgalore")

    parser.add_argument('--cutadapt_path',
                        default='cutadapt',
                        help="specify path to cutadapt")

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = parse_args()
    run_tg = RunTrimGalore(args)
    run_tg.run_trimgalore()
    run_tg.gather_outputs()
