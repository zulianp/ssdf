#!/usr/bin/env bash

set -e
set -x

threads=1
ranks=8


# mesh=/Users/patrickzulian/Desktop/code/smesh/build_debug/pump
mesh=/Users/patrickzulian/Desktop/code/smesh/build_release/pump

OMP_NUM_THREADS=$threads OMP_PROC_BIND=true  mpiexec -np $ranks ./bin/mesh_distance_to_boundary $mesh "$mesh"_edf
$CODE_DIR/smesh/python/smesh/raw_to_db.py $mesh "$mesh".vtk -p "$mesh"_edf/edf.float32