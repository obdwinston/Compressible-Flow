#!/bin/bash

# requirements
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
mkdir -p bin
mkdir -p data

# pre-processing
python -u mesh/body.py
python -u mesh/geo.py
python -u mesh/su2.py
gfortran -Jmods \
    mods/mod_mesh.f90 \
    mesh/mesh.f90 \
    -o bin/mesh && \
    ./bin/mesh

# simulation
gfortran -Jmods \
    mods/mod_mesh.f90 \
    mods/mod_config.f90 \
    mods/mod_utils.f90 \
    mods/mod_flux.f90 \
    mods/mod_solve.f90 \
    main.f90 \
    -o bin/main && \
    ./bin/main

# post-processing
python -u read.py
