from single_cell.utils import helpers
import os

def annotate_ref_alt(haps_csv, refdir, output_csv):
    thousand_genomes = os.path.join(refdir, 'thousand_genomes_snps.tsv')

    annotation_data = {}

    with helpers.getFileHandle(thousand_genomes, 'rt') as db:
        for line in db:
            line = line.strip().split('\t')

            chrom, pos, ref, alt = line

            annotation_data[(chrom, pos)] = (ref, alt)

    with helpers.getFileHandle(haps_csv, 'rt') as reader, helpers.getFileHandle(output_csv, 'wt') as writer:

        header = reader.readline()

        header += ',ref,alt\n'
        writer.write(header)

        for line in reader:
            line = line.strip()
            l_split = line.split('\t')

            chrom = l_split[0]
            pos = l_split[1]

            if (chrom, pos) in annotation_data:
                ref, alt = annotation_data[(chrom, pos)]
            else:
                ref = 'NA'
                alt = 'NA'

            line += ',{},{}\n'.format(ref, alt)

            writer.write(line)
