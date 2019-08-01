import os

import pypeliner
import pypeliner.managed as mgd
from single_cell.utils import helpers

from single_cell import infer_haps
from single_cell import variant_calling


def load_cell_data(yamldata, key):
    for sample_id, sample_data in yamldata[key].items():
        for library_id, library_data in sample_data.items():
            for cell_id, cell_data in library_data.items():
                yield sample_id, library_id, cell_id, cell_data['bam']


def load_tumour_data(yamldata):
    sample_data = {}

    for sample_id, library_id, cell_id, cell_bam in load_cell_data(yamldata, 'tumour_cells'):
        sample_data[(sample_id, library_id, cell_id)] = cell_bam

    return sample_data


def load_normal_data(yamldata):
    libraries = set()
    if 'normal_wgs' in yamldata:
        assert len(yamldata['normal_wgs'].keys()) == 1
        sample_id = list(yamldata['normal_wgs'].keys())[0]
        assert len(yamldata['normal_wgs'][sample_id].keys()) == 1
        library_id = list(yamldata['normal_wgs'][sample_id].keys())[0]
        libraries.add(library_id)
        cell_bams = yamldata['normal_wgs'][sample_id][library_id]['bam']
    else:
        cell_bams = {}

        if not len(yamldata['normal_cells'].keys()) == 1:
            raise Exception("Pipeline does not support multiple normal samples")

        for sample_id, library_id, cell_id, cell_bam in load_cell_data(yamldata, 'normal_cells'):
            libraries.add(library_id)
            if cell_id in cell_bams:
                raise Exception("non unique cell id {} encountered".format(cell_id))
            cell_bams[cell_id] = cell_bam

    return sample_id, sorted(libraries), cell_bams


def resolve_template(filepath, ids):
    filepaths = []
    for id_keys in ids:
        try:
            filepaths.append(filepath.format(**id_keys))
        except KeyError:
            filepaths.append(filepath)
    return filepaths


def generate_meta_files(normal_sample_id, normal_libraries, tumour_cell_bams, args):
    destruct_dir = args['destruct_output']
    lumpy_dir = args['lumpy_output']
    haps_dir = args['haps_output']
    varcall_dir = args["variants_output"]

    tumour_ids = sorted(set([(val[0], val[1]) for val in tumour_cell_bams]))
    tumour_ids = [
        {'sample_id': val[0], 'library_id': val[1]} for val in tumour_ids
    ]

    metadata = {
        'normal_sample': {
            'sample_id': normal_sample_id,
            'library_id': normal_libraries,
        },
        'tumour_samples': tumour_ids
    }

    if destruct_dir:
        destruct_files = get_file_paths(destruct_dir, 'destruct')
        filepaths = [resolve_template(v, tumour_ids) for v in destruct_files.values()]
        filepaths = [x for v in filepaths for x in v]
        metadata['type'] = 'destruct'
        helpers.generate_and_upload_metadata(
            args,
            destruct_dir,
            filepaths,
            metadata,
        )

    if lumpy_dir:
        lumpy_files = get_file_paths(lumpy_dir, 'lumpy')
        filepaths = [resolve_template(v, tumour_ids) for v in lumpy_files.values()]
        filepaths = [x for v in filepaths for x in v]
        metadata['type'] = 'lumpy'
        helpers.generate_and_upload_metadata(
            args,
            lumpy_dir,
            filepaths,
            metadata,
        )

    if haps_dir:
        haps_files = get_file_paths(haps_dir, 'haps_calling')
        filepaths = [resolve_template(v, tumour_ids) for v in haps_files.values()]
        filepaths = [x for v in filepaths for x in v]
        metadata['type'] = 'haplotype_calling'
        helpers.generate_and_upload_metadata(
            args,
            haps_dir,
            filepaths,
            metadata,
        )

    if varcall_dir:
        var_files = get_file_paths(varcall_dir, 'variant_calling')
        metadata['type'] = 'variant_calling'
        filepaths = [resolve_template(v, tumour_ids) for v in var_files.values()]
        filepaths = [x for v in filepaths for x in v]
        helpers.generate_and_upload_metadata(
            args,
            varcall_dir,
            filepaths,
            metadata,
        )


def multi_sample_pipeline(args):
    data = helpers.load_yaml(args['input_yaml'])

    tumour_cell_bams = load_tumour_data(data)
    normal_sample_id, normal_libraries, normal_bams = load_normal_data(data)

    pyp = pypeliner.app.Pypeline(config=args)

    workflow = create_multi_sample_workflow(
        normal_bams,
        tumour_cell_bams,
        helpers.load_config(args),
        destruct_dir=args['destruct_output'],
        lumpy_dir=args['lumpy_output'],
        haps_dir=args['haps_output'],
        varcall_dir=args["variants_output"],
        normal_sample_id=normal_sample_id
    )

    pyp.run(workflow)

    generate_meta_files(normal_sample_id, normal_libraries, tumour_cell_bams, args)


def get_file_paths(root_dir, pipeline_type):
    if pipeline_type == 'variant_calling':
        data = {
            'museq_vcfs': os.path.join(root_dir, '{sample_id}_{library_id}_museq.vcf.gz'),
            'strelka_snvs': os.path.join(root_dir, '{sample_id}_{library_id}_strelka_snv.vcf.gz'),
            'strelka_indels': os.path.join(root_dir, '{sample_id}_{library_id}_strelka_indel.vcf.gz'),
            'snv_annotations': os.path.join(root_dir, '{sample_id}_{library_id}_snv_annotations.h5'),
            'snv_counts': os.path.join(root_dir, '{sample_id}_{library_id}_snv_counts.h5'),
        }
    elif pipeline_type == 'haps_calling':
        data = {
            'haplotypes_file': os.path.join(root_dir, 'haplotypes.tsv'),
            'allele_counts_template': os.path.join(
                root_dir, '{sample_id}_{library_id}_allele_counts.csv'
            ),
        }
    elif pipeline_type == 'destruct':
        data = {
            'breakpoints_template': os.path.join(
                root_dir, '{sample_id}_{library_id}_destruct.csv.gz'
            ),
            'breakpoints_library_template': os.path.join(
                root_dir, '{sample_id}_{library_id}_destruct_library.csv.gz'
            ),
            'cell_counts_template': os.path.join(
                root_dir, '{sample_id}_{library_id}_cell_counts_destruct.csv.gz'
            )
        }
    elif pipeline_type == 'lumpy':
        data = {
            'lumpy_breakpoints_bed': os.path.join(
                root_dir, '{sample_id}_{library_id}_lumpy_breakpoints.bed'
            ),
            'lumpy_breakpoints_csv': os.path.join(
                root_dir, '{sample_id}_{library_id}_lumpy_breakpoints.csv.gz'
            ),
            'lumpy_breakpoints_evidence': os.path.join(
                root_dir, '{sample_id}_{library_id}_lumpy_breakpoints_evidence.csv.gz'),
        }
    else:
        raise Exception('unknown pipeline type: {}'.format(pipeline_type))

    return data


def create_multi_sample_workflow(
        normal_wgs_bam,
        tumour_cell_bams,
        config,
        destruct_dir=None,
        lumpy_dir=None,
        haps_dir=None,
        varcall_dir=None,
        normal_sample_id='normal',
):
    """ Multiple sample pseudobulk workflow. """

    baseimage = config['multi_sample']['docker']['single_cell_pipeline']
    ctx = {'mem_retry_increment': 2, 'disk_retry_increment': 50, 'ncpus': 1, 'docker_image': baseimage}

    workflow = pypeliner.workflow.Workflow(default_ctx=ctx)

    if isinstance(normal_wgs_bam, dict):
        workflow.setobj(
            obj=mgd.OutputChunks('normal_cell_id'),
            value=list(normal_wgs_bam.keys()),
        )
        workflow.set_filenames('normal_cells.bam', 'normal_cell_id', fnames=normal_wgs_bam)
        normal_bam = mgd.InputFile('normal_cells.bam', 'normal_cell_id', extensions=['.bai'])
    else:
        normal_bam = mgd.InputFile(normal_wgs_bam, extensions=['.bai'])

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id', 'library_id', 'cell_id'),
        value=list(tumour_cell_bams.keys()),
    )

    if varcall_dir:
        varcall_files = get_file_paths(varcall_dir, 'variant_calling')

        workflow.subworkflow(
            name='variant_calling_multi_sample',
            func=variant_calling.variant_calling_multi_sample_workflow,
            args=(
                config,
                normal_wgs_bam,
                tumour_cell_bams,
                varcall_dir,
                mgd.OutputFile('museq.vcf', 'sample_id', 'library_id', axes_origin=[],
                               template=varcall_files['museq_vcfs']),
                mgd.OutputFile('strelka_snv.vcf', 'sample_id', 'library_id', axes_origin=[],
                               template=varcall_files['strelka_snvs']),
                mgd.OutputFile('strelka_indel.vcf', 'sample_id', 'library_id', axes_origin=[],
                               template=varcall_files['strelka_indels']),
                mgd.OutputFile('snv_annotations.h5', 'sample_id', 'library_id', axes_origin=[],
                               template=varcall_files['snv_annotations']),
                mgd.OutputFile('snv_counts.h5', 'sample_id', 'library_id', axes_origin=[],
                               template=varcall_files['snv_counts']),
            )
        )

    if haps_dir:
        haps_files = get_file_paths(haps_dir, 'haps_calling')
        workflow.set_filenames('allele_counts.csv', 'sample_id', 'library_id',
                               template=haps_files['allele_counts_template'])

        workflow.subworkflow(
            name='infer_haps_from_cells_normal',
            func=infer_haps.infer_haps,
            args=(
                normal_bam,
                mgd.OutputFile(haps_files['haplotypes_file']),
                mgd.TempOutputFile("allele_counts.csv"),
                config['infer_haps'],
            ),
            kwargs={'normal': True}
        )

        workflow.subworkflow(
            name='extract_allele_readcounts',
            func=infer_haps.extract_allele_readcounts,
            axes=('sample_id', 'library_id',),
            args=(
                mgd.InputFile(haps_files['haplotypes_file']),
                mgd.InputFile('tumour_cells.bam', 'sample_id', 'library_id', 'cell_id', extensions=['.bai'],
                              fnames=tumour_cell_bams),
                mgd.OutputFile('allele_counts.csv', 'sample_id', 'library_id',
                               template=haps_files['allele_counts_template']),
                config['infer_haps'],
            ),
        )

    if destruct_dir:
        destruct_files = get_file_paths(destruct_dir, 'destruct')
        workflow.set_filenames(
            'breakpoints.h5', 'sample_id', 'library_id', template=destruct_files['breakpoints_template']
        )
        workflow.set_filenames('breakpoints_library.h5', 'sample_id', 'library_id',
                               template=destruct_files['breakpoints_library_template']
                               )
        workflow.set_filenames('cell_counts.h5', 'sample_id', 'library_id',
                               template=destruct_files['cell_counts_template'])

        destruct_pipeline_config = config['breakpoint_calling']
        destruct_config = destruct_pipeline_config.get('destruct_config', {})
        destruct_ref_data_dir = destruct_pipeline_config['ref_data_directory']
        destruct_raw_data_template = os.path.join(destruct_dir, 'rawdata')

        workflow.subworkflow(
            name='run_destruct_multi_sample',
            func='single_cell.workflows.destruct_singlecell.destruct_multi_sample_workflow',
            ctx={'docker_image': destruct_pipeline_config['docker']['destruct']},
            args=(
                normal_bam,
                mgd.InputFile('tumour_cells.bam', 'sample_id', 'library_id', 'cell_id', extensions=['.bai'],
                              axes_origin=[], fnames=tumour_cell_bams),
                destruct_config,
                destruct_pipeline_config,
                destruct_ref_data_dir,
                mgd.OutputFile('breakpoints.h5', 'sample_id', 'library_id', axes_origin=[]),
                mgd.OutputFile('breakpoints_library.h5', 'sample_id', 'library_id', axes_origin=[]),
                mgd.OutputFile('cell_counts.h5', 'sample_id', 'library_id', axes_origin=[]),
                destruct_raw_data_template,
            ),
            kwargs={'normal_sample_id': normal_sample_id}
        )

    if lumpy_dir:
        lumpy_files = get_file_paths(lumpy_dir, 'lumpy')

        workflow.set_filenames(
            'lumpy_breakpoints.csv.gz', 'sample_id', 'library_id',
            template=lumpy_files['lumpy_breakpoints_csv']
        )
        workflow.set_filenames('lumpy_breakpoints_evidence.csv.gz', 'sample_id', 'library_id',
                               template=lumpy_files['lumpy_breakpoints_evidence'])
        workflow.set_filenames('lumpy_breakpoints.bed', 'sample_id', 'library_id',
                               template=lumpy_files['lumpy_breakpoints_bed'])

        lumpy_pipeline_config = config['breakpoint_calling']

        workflow.subworkflow(
            name='run_lumpy_multi_sample',
            func='single_cell.workflows.lumpy.lumpy_multi_sample_workflow',
            ctx={'docker_image': lumpy_pipeline_config['docker']['single_cell_pipeline']},
            args=(
                lumpy_pipeline_config,
                normal_bam,
                mgd.InputFile('tumour_cells.bam', 'sample_id', 'library_id', 'cell_id', extensions=['.bai'],
                              axes_origin=[], fnames=tumour_cell_bams),
                mgd.OutputFile('lumpy_breakpoints.csv.gz', 'sample_id', 'library_id', axes_origin=[]),
                mgd.OutputFile('lumpy_breakpoints_evidence.csv.gz', 'sample_id', 'library_id', axes_origin=[]),
                mgd.OutputFile('lumpy_breakpoints.bed', 'sample_id', 'library_id', axes_origin=[]),
            )
        )

    return workflow
