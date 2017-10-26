#!/bin/bash

# $1: compiler (gcc-4.8, gcc-5, gcc-6, gcc-7, clang-4.0, clang-5.0)
# $2: build static (OFF, ON)
# $3: moab version (5.0, master)
# $4: build Dag-Geant4 (OFF, ON)

set -e

source /root/etc/$1.env
build_static=$2
moab_version=$3
build_daggeant4=$4
hdf5_version=1.8.13
geant4_version=10.04

export PATH=${install_dir}/hdf5-${hdf5_version}/bin:${PATH}
export PATH=${install_dir}/moab-${moab_version}/bin:${PATH}
export PATH=${install_dir}/geant4-${geant4_version}/bin:${PATH}
export LD_LIBRARY_PATH=${install_dir}/hdf5-${hdf5_version}/lib
export LD_LIBRARY_PATH=${install_dir}/moab-${moab_version}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${install_dir}/geant4-${geant4_version}/lib:${LD_LIBRARY_PATH}

mkdir -p ${build_dir}/DAGMC-moab-${moab_version}/bld
rm -rf ${install_dir}/DAGMC-moab-${moab_version}
cd ${build_dir}/DAGMC-moab-${moab_version}
#git clone https://github.com/svalinn/DAGMC -b develop
ln -snf DAGMC src
cd bld
cmake ../src -DMOAB_DIR=${install_dir}/moab-${moab_version} \
             -DBUILD_GEANT4=${build_daggeant4} \
             -DGEANT4_DIR=${install_dir}/geant4-${geant4_version} \
             -DBUILD_CI_TESTS=ON \
             -DBUILD_STATIC_EXE=${build_static} \
             -DCMAKE_C_COMPILER=${CC} \
             -DCMAKE_CXX_COMPILER=${CXX} \
             -DCMAKE_Fortran_COMPILER=${FC} \
             -DCMAKE_INSTALL_PREFIX=${install_dir}/DAGMC-moab-${moab_version}
make -j`grep -c processor /proc/cpuinfo`
make install
#rm -rf ${build_dir}/DAGMC-moab-${moab_version}
rm -rf ${build_dir}/DAGMC-moab-${moab_version}/bld
