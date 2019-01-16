#!/bin/bash

set -e

git clone https://github.com/pyne/pyne
cd pyne/src
python atomicgen.py
cd ..
python amalgamate.py -f license.txt src/utils.* src/extra_types.h src/h5wrap.h \
    src/state_map.cpp src/nucname.* src/rxname.* src/particle.* src/data.* \
    src/json-forwards.h src/json.h src/jsoncpp.cpp src/jsoncustomwriter.* \
    src/material.* src/tally.* src/atomic_data.*
cp pyne.cpp pyne.h ..
cd ..
rm -rf pyne
astyle --options=../../astyle_google.ini --suffix=none --verbose --formatted \
       "pyne.cpp" "pyne.h"
python remove_unsupported.py
mv -f pyne.cpp.new pyne.cpp
