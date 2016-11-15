#!/bin/bash

set -e

export LD_LIBRARY_PATH=/root/geant4.10.00.p02/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/root/moab/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/root/dagmc/lib

# move to test dir
cd tests
# run each test
./dagsolid_unit_tests
# no fludag yet
# ./fludag_unit_tests
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
./make_watertight_cylinder_tests
./make_watertight_cone_tests
# if this is not a pull request, run regression tests
if [ ! -z $TRAVIS_PULL_REQUEST ] && [ $TRAVIS_PULL_REQUEST == "false" ] ; then
    wget $MW_REG_TEST_MODELS_URL -O mw_reg_test_files.tar.gz -o wget.out
    tar xzvf mw_reg_test_files.tar.gz
    ./make_watertight_regression_tests
fi
# move to the base directory
cd ..
# remove the bld dir
rm -rf bld
# run astyle to see if there are any differences
tools/astyle_google --style=linux --indent=spaces=2 --exclude=gtest -r *.cpp *.h *.hpp *.cc *.hh	
# checks for C++ style guide adherence
git diff --exit-code
