import pypeliner
import pypeliner.managed as mgd
from single_cell.workflows.extract_allele_readcounts.dtypes import dtypes


def extract_allele_readcounts(
        haplotypes_filename,
        cell_bams,
        allele_counts_filename,
        config,
):
    baseimage = {'docker_image': config['docker']['single_cell_pipeline']}

    remixt_image = config['docker']['remixt']

    remixt_config = config.get('extract_seqdata', {})
    remixt_ref_data_dir = config['ref_data_dir']

    workflow = pypeliner.workflow.Workflow(ctx=baseimage)

    workflow.set_filenames('cell.bam', 'cell_id', fnames=cell_bams)

    workflow.setobj(
        obj=mgd.OutputChunks('cell_id'),
        value=list(cell_bams.keys()),
    )

    workflow.subworkflow(
        name='create_seqdata_readcounts',
        axes=('cell_id',),
        func='remixt.workflow.create_extract_seqdata_workflow',
        ctx={'docker_image': remixt_image},
        args=(
            mgd.InputFile('cell.bam', 'cell_id', extensions=['.bai']),
            mgd.TempOutputFile('seqdata.h5', 'cell_id', axes_origin=[]),
            remixt_config,
            remixt_ref_data_dir,
        ),
        kwargs={'no_parallelism': True}
    )

    # TODO Segments with bin width from single cell
    workflow.transform(
        name='create_segments',
        func='remixt.analysis.segment.create_segments',
        ctx={'mem': 16, 'docker_image': remixt_image},
        args=(
            mgd.TempOutputFile('segments.tsv'),
            remixt_config,
            remixt_ref_data_dir,
        ),
    )

    workflow.transform(
        name='generate_haplotypes_tsv',
        func='single_cell.workflows.extract_allele_readcounts.tasks.convert_csv_to_tsv',
        args=(
            mgd.InputFile(haplotypes_filename),
            mgd.TempOutputFile('haplotypes.tsv')
        )
    )

    workflow.transform(
        name='haplotype_allele_readcount',
        axes=('cell_id',),
        ctx={'mem': 16, 'docker_image': remixt_image},
        func='remixt.analysis.readcount.haplotype_allele_readcount',
        args=(
            mgd.TempOutputFile('allele_counts.tsv', 'cell_id', axes_origin=[]),
            mgd.TempInputFile('segments.tsv'),
            mgd.TempInputFile('seqdata.h5', 'cell_id'),
            mgd.TempInputFile('haplotypes.tsv'),
            remixt_config,
        ),
    )

    workflow.transform(
        name='prep_readcount_csv',
        axes=('cell_id',),
        func='single_cell.utils.csvutils.rewrite_csv_file',
        args=(
            mgd.TempInputFile('allele_counts.tsv', 'cell_id'),
            mgd.TempOutputFile('allele_counts.csv.gz', 'cell_id', extensions=['.yaml']),
        ),
        kwargs={
            'write_header': True,
            'dtypes': dtypes()['readcount']
        },
    )

    workflow.transform(
        name='merge_allele_readcount',
        ctx={'mem': 16},
        func='single_cell.utils.csvutils.concatenate_csv',
        args=(
            mgd.TempInputFile('allele_counts.csv.gz', 'cell_id', extensions=['.yaml']),
            mgd.OutputFile(allele_counts_filename, extensions=['.yaml']),
        ),
        kwargs={
            'write_header': True
        },
    )

    return workflow
