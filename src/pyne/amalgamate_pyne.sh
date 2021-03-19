#!/bin/bash

set -e

# Update amalgamated pyne
ls $1/src
cd $1/src
# ls
# python atomicgen.py
# cd ..
# python amalgamate.py -f license.txt src/utils.* src/extra_types.h src/h5wrap.h \
#     src/state_map.cpp src/nucname.* src/rxname.* src/particle.* src/data.* \
#     src/json-forwards.h src/json.h src/jsoncpp.cpp src/jsoncustomwriter.* \
#     src/material.* src/material_library.* src/tally.* src/atomic_data.* src/measure.* \
#     src/source_sampling.*
# cp pyne.cpp pyne.h $2
# githash=`git rev-parse HEAD`
# cd ..
# astyle --options=../../astyle_google.ini --suffix=none --verbose --formatted \
#        "pyne.cpp" "pyne.h"
# python remove_unsupported.py
# mv -fv pyne.cpp.new pyne.cpp

# # Update source.F90
# cp -v pyne/share/source.F90       ../mcnp/mcnp5/pyne_source/source.F90
# cp -v pyne/share/source_mcnp6.F90 ../mcnp/mcnp6/pyne_source/source.F90

# # Delete pyne repo
# rm -rf pyne

