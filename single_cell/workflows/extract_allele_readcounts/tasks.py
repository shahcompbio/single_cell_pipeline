from single_cell.utils import csvutils


def convert_csv_to_tsv(csv_infile, tsv_outfile):
    csvinput = csvutils.CsvInput(csv_infile)

    csvdata = csvinput.read_csv()

    csvdata.to_csv(tsv_outfile, sep='\t', index=False)
