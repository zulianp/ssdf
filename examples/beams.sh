#!/usr/bin/env bash

set -e

./create_beam_mesh.py


X_PROP=100
Y_PROP=100
Z_PROP=338
RES=6

set -x

# Usage: ../build_openmp/bin/beam_to_edf <beam_mesh_folder> <nx> <ny> <nz> <output_folder>
SSDF_MARGIN=0.3 ../build_openmp/bin/beam_to_edf . $(( RES *  X_PROP)) $(( RES *  Y_PROP)) $(( RES *  Z_PROP)) beam_edf
