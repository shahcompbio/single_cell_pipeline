import os
import shutil

import pandas as pd
from single_cell.utils import csvutils
from single_cell.utils import helpers, hdfutils
from single_cell.utils import pdfutils

import pypeliner

scripts_directory = os.path.join(
    os.path.realpath(
        os.path.dirname(__file__)),
    'scripts')
run_aneufinder_rscript = os.path.join(scripts_directory, 'Aneufinder.R')
rdata_to_csv_rscript = os.path.join(scripts_directory, 'Rdatatocsv.R')


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
        cell_id,
        segments,
        reads,
        dnacopy_plot,
        docker_image=None):
    # Create an output folder for temp storage
    helpers.makedirs(working_dir)

    bam_name = os.path.basename(bam_file)
    # Symlink to the original bam file
    fake_bam_file = os.path.join(working_dir, bam_name)
    helpers.copy_file(bam_file, fake_bam_file)

    bai_file = bam_file + '.bai'
    if os.path.exists(bai_file):
        bai_name = bam_name + '.bai'
        fake_bai_file = os.path.join(working_dir, bai_name)
        helpers.copy_file(bai_file, fake_bai_file)

    # Give aneufinder a place to output
    temp_output = os.path.join(working_dir, 'output')
    helpers.makedirs(temp_output)

    # Run Aneufinder
    cmd = ['Rscript', run_aneufinder_rscript, working_dir, temp_output]

    try:
        pypeliner.commandline.execute(*cmd, docker_image=docker_image)
    except:
        print('Aneufinder failed on {}'.format(bam_file))

        with open(segments, 'w') as segfile:
            header = ['chr', 'start', 'end', 'mean.counts', 'state', 'cell_id']
            header = ','.join(header)
            segfile.write(header)
        with open(reads, 'w') as readsfile:
            header = ['chr', 'start', 'end', 'reads', 'GC', 'state']
            header = ','.join(header)
            readsfile.write(header)
        file(dnacopy_plot, 'w').close()
        return

    # Copy out the plot files. Grabs the first bin size that it sees
    all_plots = os.listdir(
        os.path.join(
            temp_output,
            'PLOTS',
            'method-dnacopy'))

    cell_plot = [x for x in all_plots if x.startswith('profiles')][0]
    cell_plot = os.path.join(temp_output, 'PLOTS', 'method-dnacopy', cell_plot)

    shutil.move(cell_plot, dnacopy_plot)

    # Get the segment data into a csv
    segments_rdata = os.listdir(
        os.path.join(
            temp_output,
            'MODELS',
            'method-dnacopy'))
    if len(segments_rdata) != 1:
        raise Exception("Wrong number of segment files")
    segments_rdata = segments_rdata[0]
    segments_rdata = os.path.join(
        temp_output,
        'MODELS',
        'method-dnacopy',
        segments_rdata)

    cmd = ['Rscript', rdata_to_csv_rscript, segments_rdata, segments, reads]

    pypeliner.commandline.execute(*cmd, docker_image=docker_image)

    convert_segments_to_hmmcopy_format(segments, cell_id)
    convert_reads_to_hmmcopy_format(reads, cell_id)


def merge_pdf(in_filenames, out_filename):
    for in_files, out_file in zip(in_filenames, out_filename):
        pdfutils.merge_pdfs(in_files, out_file)


def merge_outputs_to_hdf(
        reads_files, segs_files, outfile, tempdir):
    helpers.makedirs(tempdir)

    reads_csv_merged = os.path.join(tempdir, "merged_reads.csv")

    csvutils.concatenate_csv_lowmem(reads_files, reads_csv_merged)

    segs_csv_merged = os.path.join(tempdir, "merged_segs.csv")

    csvutils.concatenate_csv_lowmem(segs_files, segs_csv_merged)

    hdfutils.concat_csvs_to_hdf([reads_csv_merged,
                                 segs_csv_merged],
                                outfile,
                                ['/aneufinder/reads/',
                                 '/aneufinder/segments/'])
