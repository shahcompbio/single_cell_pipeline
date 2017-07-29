'''
Created on Jan 4, 2017

@author: dgrewal
'''

from subprocess import Popen
import os
import warnings

class RunTrimGalore(object):
    """
    Trim fastq files with trimgalore
    """

    def __init__(self, seq1, seq2, fq_r1, fq_r2, trimgalore, cutadapt, tempdir,
                adapter, adapter2, report_r1, report_r2, qc_report_r1,
                qc_report_r2, qc_zip_r1, qc_zip_r2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.trimgalore_path = trimgalore
        self.cutadapt_path = cutadapt
        self.tempdir = tempdir
        self.adapter = adapter
        self.adapter2 = adapter2
        self.fastqc_report_r1 = qc_report_r1
        self.fastqc_report_r2 = qc_report_r2
        self.fastqc_zip_r1 = qc_zip_r1
        self.fastqc_zip_r2 = qc_zip_r2
        self.fastq_r1 = fq_r1
        self.fastq_r2 = fq_r2
        self.report_r1 = report_r1
        self.report_r2 = report_r2

        self.check_inputs()

        if not os.path.exists(self.tempdir):
            os.makedirs(self.tempdir)

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

        ext1 = self.check_file(self.seq1)
        ext2 = self.check_file(self.seq2)

        assert ext1 == ext2, "input fastqfiles should have the same extension"

    def run_trimgalore(self):
        """
        launch trimgalore
        """
        cmd = [self.trimgalore_path,
               '--fastqc',
               '--paired',
               '--path_to_cutadapt', self.cutadapt_path,
               '--output_dir', self.tempdir + '/',
               ]

        if self.adapter:
            cmd.extend(['--adapter', self.adapter])

        if self.adapter2:
            cmd.extend(['--adapter2', self.adapter2])

        cmd.extend([self.seq1, self.seq2, ])

        self.run_cmd(cmd)

    def move_files(self, fname, outpath):
        """
        move files from from temp dir to the expected path
        """
        path = os.path.join(self.tempdir, fname)
        os.rename(path, outpath)
        assert os.path.isfile(outpath)

        print outpath
        print os.stat(outpath)


    def get_file(self, r1_out, r2_out, ext, trimreport=False):
        """
        find the file in outdir and rename to r1_out and r2_out
        """
        outfiles = os.listdir(self.tempdir)

        reps = [v for v in outfiles if ext in v]
        assert reps != [], "Couldn't move %s files" % ext

        seq1 = os.path.basename(self.seq1)
        seq2 = os.path.basename(self.seq2)

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
            self.report_r1,
            self.report_r2,
            "_trimming_report.txt",
            trimreport=True)

        # fastqc html files
        self.get_file(
            self.fastqc_report_r1,
            self.fastqc_report_r2,
            ".html")

        # fastqc zip files
        self.get_file(self.fastqc_zip_r1, self.fastqc_zip_r2, ".zip")

        # trimmed fastq files
        self.get_file(self.fastq_r1, self.fastq_r2, ".fq.gz")
