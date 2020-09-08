cd $SRC_DIR/single_cell_pipeline
$PYTHON setup.py install

cd $SRC_DIR/scgenome
$PYTHON setup.py install

cd $SRC_DIR/wgs_analysis
$PYTHON setup.py install

cd $SRC_DIR/pypeliner
$PYTHON setup.py install
