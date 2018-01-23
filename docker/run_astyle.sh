#!/bin/bash

set -e

# Run astyle to see if there are any differences
astyle_deb=astyle_3.0.1-1ubuntu1_amd64.deb
wget http://archive.ubuntu.com/ubuntu/pool/universe/a/astyle/${astyle_deb}
dpkg -i ${astyle_deb}
rm -f ${astyle_deb}
astyle --options=astyle_google.ini \
       --exclude=gtest \
       --exclude=mcnp/mcnp5/Source \
       --exclude=mcnp/mcnp6/Source \
       --ignore-exclude-errors \
       --recursive \
       --verbose \
       --formatted \
       "*.cc" "*.cpp" "*.h" "*.hh" "*.hpp"
# Exit if astyle found diffs
git diff --exit-code