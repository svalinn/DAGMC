#!/bin/bash

# $1: compiler (gcc-4.8, gcc-5, gcc-6, clang-3.8)
# $2: moab version (4.9.2, master)

set -e

source /root/etc/$1.env
moab_version=$2

export PATH=${install_dir}/hdf5-${hdf5_version}/bin:${PATH}
export PATH=${install_dir}/moab-${moab_version}/bin:${PATH}
export PATH=${install_dir}/geant4-${geant4_version}/bin:${PATH}
export LD_LIBRARY_PATH=${install_dir}/hdf5-${hdf5_version}/lib
export LD_LIBRARY_PATH=${install_dir}/moab-${moab_version}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${install_dir}/geant4-${geant4_version}/lib:${LD_LIBRARY_PATH}

echo $PATH
echo $LD_LIBRARY_PATH

ls -l ${install_dir}/moab-${moab_version}/bin
ls -l ${install_dir}/moab-${moab_version}/lib

mkdir -p ${build_dir}/DAGMC-moab-${moab_version}/bld
rm -rf ${install_dir}/DAGMC-moab-${moab_version}
cd ${build_dir}/DAGMC-moab-${moab_version}
#git clone https://github.com/svalinn/DAGMC -b develop
ln -snf DAGMC src
cd bld
cmake ../src -DBUILD_TALLY=ON \
             -DBUILD_CI_TESTS=ON \
             -DBUILD_GEANT4=ON \
             -DGEANT4_DIR=${install_dir}/geant4-${geant4_version} \
             -DCMAKE_INSTALL_PREFIX=${install_dir}/DAGMC-moab-${moab_version}
make -j`grep -c processor /proc/cpuinfo`
make install
#rm -rf ${build_dir}/DAGMC-moab-${moab_version}
