import os
import pandas as pd
import pypeliner
import shutil
from single_cell.workflows.hmmcopy.scripts.plot_hmmcopy import GenHmmPlots
from PyPDF2 import PdfFileMerger

scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
run_aneufinder_rscript = os.path.join(scripts_directory, 'Aneufinder.R')
rdata_to_csv_rscript = os.path.join(scripts_directory, 'Rdatatocsv.R')


def symlink(actual_file, symlink):
    shutil.copy(actual_file, symlink)
#    if not os.path.exists(symlink):
#        os.symlink(actual_file, symlink)


def makedir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)


def convert_segments_to_hmmcopy_format(
    csv_name,
    cell_id):

    rdata_df = pd.read_csv(csv_name)

    df = rdata_df[['seqnames', 'start', 'end', 'mean.counts', 'state']]
    df = df.rename(columns={'seqnames': 'chr'})
    df['state'] = df['state'].apply(lambda x: int(x[0]))
    df['cell_id'] = cell_id

    df.to_csv(csv_name, index=False)


def convert_reads_to_hmmcopy_format(
    csv_name,
    cell_id):

    rdata_df = pd.read_csv(csv_name)

    df = rdata_df[['seqnames', 'start', 'end', 'counts', 'GC', 'state']]
    df = df.rename(columns={'seqnames': 'chr', 'counts': 'reads', 'GC': 'gc'})
    df['state'] = df['state'].apply(lambda x: int(x[0]))
    df['cell_id'] = cell_id

    df.to_csv(csv_name, index=False)


def run_aneufinder(
    bam_file,
    working_dir,
    sample_id,
    aneufinder_output,
    segments,
    reads,
    dnacopy_plot):

    # Create an output folder for temp storage
    makedir(working_dir)

    bam_name = os.path.basename(bam_file)
    # Symlink to the original bam file
    fake_bam_file = os.path.join(working_dir, bam_name)
    symlink(bam_file, fake_bam_file)

    bai_file = bam_file + '.bai'
    if os.path.exists(bai_file):
        bai_name = bam_name + '.bai'
        fake_bai_file = os.path.join(working_dir, bai_name)
        symlink(bai_file, fake_bai_file)

    # Give aneufinder a place to output
    temp_output = os.path.join(working_dir, 'output')
    makedir(temp_output)

    # Run Aneufinder
    cmd = ['Rscript', run_aneufinder_rscript, working_dir, temp_output]

    try:
        pypeliner.commandline.execute(*cmd)
    except:
        print('Aneufinder failed on {}'.format(bam_file))
        open(segments, 'w').close()
        open(reads, 'w').close()
        file(dnacopy_plot, 'w').close()

    # Copy out the plot files. Grabs the first bin size that it sees
    all_plots = os.listdir(os.path.join(temp_output, 'PLOTS', 'method-dnacopy'))
    cell_plot = [x for x in all_plots if x.startswith('profiles')][0]
    cell_plot = os.path.join(temp_output, 'PLOTS', 'method-dnacopy', cell_plot)

    shutil.move(cell_plot, dnacopy_plot)

    # Get the segment data into a csv
    segments_rdata = os.listdir(os.path.join(temp_output, 'MODELS', 'method-dnacopy'))
    if len(segments_rdata) != 1:
        raise Exception("Wrong number of segment files")
    segments_rdata = segments_rdata[0]
    segments_rdata = os.path.join(temp_output, 'MODELS', 'method-dnacopy', segments_rdata)

    cmd = ['Rscript', rdata_to_csv_rscript, segments_rdata, segments, reads]

    pypeliner.commandline.execute(*cmd)

    convert_segments_to_hmmcopy_format(segments, sample_id)
    convert_reads_to_hmmcopy_format(reads, sample_id)


def merge_pdf(in_filenames, out_filename):

    for in_files, out_file in zip(in_filenames, out_filename):

        outdir = os.path.dirname(out_file)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        merger = PdfFileMerger()

        for samp, infile in in_files.iteritems():
            if os.stat(infile).st_size > 0:
                merger.append(open(infile, 'rb'))
        if not os.path.exists(os.path.dirname(out_file)):
            os.makedirs(os.path.dirname(out_file))
        with open(out_file, 'wb') as fout:
            merger.write(fout)
