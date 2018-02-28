#!/bin/bash

# $1: compiler (gcc-4.8, gcc-5, gcc-6, gcc-7, clang-4.0, clang-5.0)
# $2: moab version (5.0, master)
# $3: build Dag-Geant4 (OFF, ON)

set -e

source /root/etc/$1.env
moab_version=$2
build_daggeant4=$3
hdf5_version=1.8.13
geant4_version=10.04

DAGMC_dir=${install_dir}/DAGMC-moab-${moab_version}


# clean out config test directory for next build
git clean -dxf . 

export PATH=${install_dir}/hdf5-${hdf5_version}/bin:${PATH}
export PATH=${install_dir}/moab-${moab_version}/bin:${PATH}
export PATH=${DAGMC_dir}/bin:${PATH}
export LD_LIBRARY_PATH=${install_dir}/hdf5-${hdf5_version}/lib
export LD_LIBRARY_PATH=${install_dir}/moab-${moab_version}/lib:${LD_LIBRARY_PATH}

# test DAGMC CMake configuration file using LD_LIBRARY_PATH to find DAGMC
cd ${build_dir}/DAGMC-moab-${moab_version}/src/cmake/test_config
cmake . -DDAGMC_ROOT=${DAGMC_dir}
make all test

# add DAGMC installation dir to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${DAGMC_dir}/lib:${LD_LIBRARY_PATH}

# test DAGMC CMake configuration file using LD_LIBRARY_PATH to find DAGMC
cd ${build_dir}/DAGMC-moab-${moab_version}/src/cmake/test_config
cmake .
make all test

# Run the tests
cd ${DAGMC_dir}/tests
./dagmc_unit_tests
./dagmc_pointinvol_test
./dagmc_rayfire_test
./dagmc_simple_test
if [ "$build_daggeant4" == "ON" ]; then
  PATH=${install_dir}/geant4-${geant4_version}/bin:${PATH}
  LD_LIBRARY_PATH=${install_dir}/geant4-${geant4_version}/lib:${LD_LIBRARY_PATH}
  ./dagsolid_unit_tests
fi
#./fludag_unit_tests  # no fludag yet
./test_CellTally
./test_KDEKernel
./test_KDEMeshTally
./test_KDENeighborhood
./test_PolynomialKernel
./test_Quadrature
./test_Tally
./test_TallyData
./test_TallyEvent
./test_TrackLengthMeshTally
./uwuw_unit_tests
./uwuw_unit_tests_tally
./uwuw_unit_tests_preprocessor
./make_watertight_no_curve_sphere_tests
./make_watertight_sphere_n_box_test
./make_watertight_cylinder_tests
./make_watertight_cone_tests
./dagmcnp_unit_tests

# If this is not a pull request, run regression tests
if [ ! -z $TRAVIS_PULL_REQUEST ] && [ $TRAVIS_PULL_REQUEST == "false" ] ; then
  wget $MW_REG_TEST_MODELS_URL -O mw_reg_test_files.tar.gz -o wget.out
  tar xzvf mw_reg_test_files.tar.gz
  ./make_watertight_regression_tests
  rm -rf *.h5m mw_reg_test_files.tar.gz
  if [ $? != 0 ]; then
    exit 1
  fi
fi
