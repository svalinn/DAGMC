#!/bin/bash

# $1: compiler (gcc-4.8, gcc-4.9, gcc-5, gcc-6, clang-3.8, clang-3.9, clang-4.0)
# $2: moab version (4.9.2, master)
# $3: geant4 version (10.01.p03, 10.02.p03, 10.03.p01)

set -e

source /root/etc/$1.env
hdf5_version=1.8.13
moab_version=$2
geant4_version=$3

DAGMC_dir=${install_dir}/DAGMC-moab-${moab_version}

export PATH=${install_dir}/hdf5-${hdf5_version}/bin:${PATH}
export PATH=${install_dir}/moab-${moab_version}/bin:${PATH}
export PATH=${install_dir}/geant4-${geant4_version}/bin:${PATH}
export PATH=${DAGMC_dir}/bin:${PATH}
export LD_LIBRARY_PATH=${install_dir}/hdf5-${hdf5_version}/lib
export LD_LIBRARY_PATH=${install_dir}/moab-${moab_version}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${install_dir}/geant4-${geant4_version}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${DAGMC_dir}/lib:${LD_LIBRARY_PATH}

# Run astyle to see if there are any differences
${DAGMC_dir}/tools/astyle --options=${DAGMC_dir}/tools/google.ini \
                          --exclude=gtest \
                          --exclude=tools/astyle \
                          --exclude=mcnp/mcnp5/Source \
                          --exclude=mcnp/mcnp6/Source \
                          --ignore-exclude-errors \
                          --recursive \
                          --verbose \
                          --formatted \
                          "*.cc" "*.cpp" "*.h" "*.hh" "*.hpp"
git diff --exit-code  # Exit if astyle found diffs

# Run the tests
cd ${DAGMC_dir}/tests
./dagmc_unit_tests
./dagmc_pointinvol_test
./dagmc_rayfire_test
./dagmc_simple_test
./dagsolid_unit_tests
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
  if [ $? != 0 ]; then
    exit 1
  fi
fi

# Cleanup
cd
rm -rf ${build_dir}/DAGMC-moab-${moab_version}
