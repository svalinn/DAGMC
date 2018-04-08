#!/bin/bash

set -e

# Run astyle to see if there are any differences
astyle_deb=astyle_3.1-1ubuntu2_amd64.deb
wget http://archive.ubuntu.com/ubuntu/pool/universe/a/astyle/${astyle_deb}
dpkg -i ${astyle_deb}
rm -f ${astyle_deb}
astyle --options=astyle_google.ini \
       --exclude=gtest \
       --exclude=mcnp/mcnp5/Source \
       --exclude=mcnp/mcnp6/Source \
       --ignore-exclude-errors \
       --suffix=none \
       --recursive \
       --verbose \
       --formatted \
       "*.cc" "*.cpp" "*.h" "*.hh" "*.hpp"

# Exit if astyle found diffs
diffs=`git status --porcelain`
if [ -z "${diffs}" ]; then
  echo "Style guide checker passed!"
else
  echo "ERROR: Style guide checker failed. Please run astyle."
  git diff
  exit 1
fi
