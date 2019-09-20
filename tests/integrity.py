
import logging
from scgenome.loaders.qc import load_qc_data
from single_cell.utils import storageutils
from single_cell.utils import helpers

LOGGING_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


dtypes_check = {
    'align_metrics': {
        'cell_id': 'category',
        'total_mapped_reads': 'int64',
        'library_id': 'category',
    },
    'hmmcopy_reads': {
        'chr': 'category',
        'start': 'int64',
        'end': 'int64',
        'reads': 'int64',
        'gc': 'float64',
        'copy': 'float64',
        'state': 'int64',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
    },
    'hmmcopy_segs': {
        'chr': 'category',
        'start': 'int64',
        'end': 'int64',
        'state': 'int64',
        'median': 'float64',
        'multiplier': 'int64',
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
    },
    'hmmcopy_metrics': {
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'order': 'float64',
        'total_mapped_reads_hmmcopy': 'int64',
        'experimental_condition': 'object',
        'mean_copy': 'float64',
        'state_mode': 'int64',
    },
    'annotation_metrics': {
        'cell_id': 'category',
        'sample_id': 'category',
        'library_id': 'category',
        'total_mapped_reads_hmmcopy': 'int64',
        'mean_copy': 'float64',
        'state_mode': 'int64',
        'order': 'float64',
        'experimental_condition': 'object',
        'quality': 'float64',
        'is_s_phase': 'boolean',
    },
    'gc_metrics': {
    },
}


def test_qc_data(results_tables):
    for table_name, table_data in results_tables.items():
        logging.info(f'table {table_name} has size {len(table_data)}')
        for column_name, dtype_name in dtypes_check[table_name].items():
            column_dtype = str(results_tables[table_name][column_name].dtype)
            if not column_dtype == dtype_name:
                raise Exception(f'{column_name} has dtype {column_dtype} not {dtype_name}')


def test_load_local_qc_data(results_dir):
    results_tables = load_qc_data(results_dir)
    test_qc_data(results_tables)


def integrity_check(blob_prefix, tempdir):
    helpers.makedirs(tempdir)

    storageutils.download_blobs(blob_prefix, tempdir, storage='azureblob')

    test_load_local_qc_data(tempdir)
