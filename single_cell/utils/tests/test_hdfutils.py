import pandas as pd
import numpy as np

from single_cell.utils.hdfutils import convert_hdf_to_csv


def test_convert_hdf_to_csv(tmp_path):
    temp_h5_path = (tmp_path / 'test.h5').as_posix()

    data = pd.DataFrame(np.random.randint(0,100,size=(100, 4)), columns=list('ABCD'))

    with pd.HDFStore(temp_h5_path, 'w') as store:
        store.put('test', data, format='table')

    temp_csv_path = (tmp_path / 'test.csv').as_posix()
    convert_hdf_to_csv(temp_h5_path, {'test': temp_csv_path}, chunksize=2)

    data_result = pd.read_csv(temp_csv_path)

    assert data.equals(data_result)
