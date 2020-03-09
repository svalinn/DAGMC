#!/bin/bash

set -e

source ${docker_env}

rm -rf ${double_down_build_dir}/bld ${double_down_install_dir}
mkdir -p ${double_down_build_dir}/bld
cd ${double_down_build_dir}
git clone https://github.com/pshriwise/double_down
ln -snf double_down src
/bin/bash ${double_down_build_dir}/src/ci/embree-install.sh
cd bld
LD_LIBRARY_PATH=${moab_install_dir}/lib:$HOME/EMBREE/lib/:$LD_LIBRARY_PATH \
cmake ../src -DCMAKE_INSTALL_PREFIX=${double_down_install_dir}
make -j${jobs}
make test
make install
cd
