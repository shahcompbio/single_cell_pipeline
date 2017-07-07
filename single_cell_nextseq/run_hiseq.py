import os
import argparse
import utils
import pypeliner
import pypeliner.managed as mgd
from workflows import fastqc_workflow



scripts_directory = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'scripts')
plot_hmmcopy_script = os.path.join(scripts_directory, 'plot_hmmcopy.py')



def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    pypeliner.app.add_arguments(parser)

    parser.add_argument('hiseq_dir',
                        help='''Path to input hiseq directory.''')

    parser.add_argument('samplesheet',
                        help='''Path to the samplesheet.''')

    parser.add_argument('out_dir',
                        help='''Path to output files.''')

    parser.add_argument('config_file',
                        help='''Path to yaml config file.''')

    args = vars(parser.parse_args())
    args['tmpdir'] = os.path.join(args['out_dir'], 'tmp')

    return args


def main():

    args = parse_args()

    pyp = pypeliner.app.Pypeline(config=args)

    config = utils.load_config(args)

    fastq_directory = os.path.join(args['out_dir'], 'fastq')


    _, _, sample_ids, fastq_1_filenames, fastq_2_filenames = utils.read_samplesheet_hiseq(args,fastq_directory)

    trimgalore_results_template_r1 = os.path.join(args['out_dir'], 'trim', '{sample_id}_R1.fastq.gz')
    trimgalore_results_template_r2 = os.path.join(args['out_dir'], 'trim', '{sample_id}_R2.fastq.gz')


    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=sample_ids,
    )

    workflow.subworkflow(
        name='fastqc_workflow',
        axes=('sample_id',),
        func=fastqc_workflow.create_fastqc_workflow,
        args=(mgd.InputFile('fastq_1', 'sample_id', fnames=fastq_1_filenames),
              mgd.InputFile('fastq_2', 'sample_id', fnames=fastq_2_filenames),
              mgd.OutputFile('fastq_trim_1', 'sample_id', template=trimgalore_results_template_r1),
              mgd.OutputFile('fastq_trim_2', 'sample_id', template=trimgalore_results_template_r2),
              config,
              mgd.InputInstance('sample_id'),
              args
            ),
        )

    pyp.run(workflow)

if __name__ == '__main__':
    main()
