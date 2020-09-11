#!/bin/bash

# Build and install python package
cd $SRC_DIR/single_cell_pipeline
$PYTHON setup.py install

cd $SRC_DIR/cell_cycle_classifier
$PYTHON setup.py install
