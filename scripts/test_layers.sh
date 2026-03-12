#!/usr/bin/env bash

set -e
set -x

threads=8
ranks=1

$INSTALL_DIR/smesh/bin/cube TET4 100 100 100 -1 -1 -1 1 1 1 vol
$INSTALL_DIR/smesh/bin/cube TET4 10 10 10 -0.5 -0.5 -0.5 0.5 0.5 0.5 vol_temp
$INSTALL_DIR/smesh/bin/skin vol_temp surf

mesh=vol
surface=surf

# mesh=/Users/patrickzulian/Desktop/code/smesh/build_debug/pump
# surface=/Users/patrickzulian/Desktop/code/smesh/build_debug/elastic

# mesh=/Users/patrickzulian/Desktop/code/smesh/build_release/pump
# $CODE_DIR/smesh/python/smesh/raw_to_db.py $surface "$surface".vtk

OMP_NUM_THREADS=$threads OMP_PROC_BIND=true  mpiexec -np $ranks ./bin/mesh_layers $mesh $surface "$mesh"_edf
$CODE_DIR/smesh/python/smesh/raw_to_db.py $mesh "$mesh".vtk -p "$mesh"_edf/layer.float32