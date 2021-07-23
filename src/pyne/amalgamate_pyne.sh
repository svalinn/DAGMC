#!/bin/bash

# take 2 arguments:
# $1 Python executable
# $2 location of the pyne source folder
# $3 destination of the amalgamated files

set -e

# Update amalgamated pyne

cd $2/pyne/src
$1 atomicgen.py
cd ..
$1 amalgamate.py -f license.txt src/utils.* src/extra_types.h src/h5wrap.h \
    src/state_map.cpp src/nucname.* src/rxname.* src/particle.* src/data.* \
    src/json-forwards.h src/json.h src/jsoncpp.cpp src/jsoncustomwriter.* \
    src/material.* src/material_library.* src/tally.* src/atomic_data.* src/measure.* \
    src/source_sampling.*
cp pyne.cpp pyne.h $3

githash=`git rev-parse HEAD`
cd $3
$1 $2/remove_unsupported.py
mv -fv pyne.cpp.new pyne.cpp

# Update source.F90
cp -v $2/pyne/share/source.F90       $3/../mcnp/mcnp5/pyne_source/source.F90
cp -v $2/pyne/share/source_mcnp6.F90 $3/../mcnp/mcnp6/pyne_source/source.F90

