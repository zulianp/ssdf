#!/usr/bin/env bash

set -e
set -x

threads=1
ranks=8

# threads=8
# ranks=1

# rm -rf vol vol_temp surf vol_edf
# $INSTALL_DIR/smesh/bin/cube TET4 10 10 300 -2 -2 -2       2   2   58 vol
# $INSTALL_DIR/smesh/bin/cube TET4 10 10 10 -0.5 -0.5 -0.5 0.5 0.5 0.5 vol_temp
# $INSTALL_DIR/smesh/bin/skin vol_temp surf

mesh=vol
surface=surf


$INSTALL_DIR/smesh/bin/skin /Users/patrickzulian/Desktop/code/smesh/build_release/elastic /Users/patrickzulian/Desktop/code/smesh/build_release/elastic_skin
mesh=/Users/patrickzulian/Desktop/code/smesh/build_release/pump
surface=/Users/patrickzulian/Desktop/code/smesh/build_release/elastic_skin

# mesh=/Users/patrickzulian/Desktop/code/smesh/build_release/pump
# $CODE_DIR/smesh/python/smesh/raw_to_db.py $surface "$surface".vtk

OMP_NUM_THREADS=$threads OMP_PROC_BIND=true  mpiexec -np $ranks ./bin/mesh_layers $mesh $surface "$mesh"_edf

layer="$mesh"_edf/layer.float32
sdf="$mesh"_edf/sdf.float32
$CODE_DIR/smesh/python/smesh/raw_to_db.py $mesh "$mesh".vtk -p $layer,$sdf