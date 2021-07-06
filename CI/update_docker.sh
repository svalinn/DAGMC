#!/bin/bash

set -ex

PUSH=true

# check for dry-run flag
for arg in "$@"
do
  case $arg in
      -n|--no-push)
      PUSH=false
      shift
      ;;
  esac
done

ubuntu_versions="16.04 18.04"
compilers="gcc clang"
hdf5_versions="1.10.4"
moab_versions="9c96d17 develop master"
for ubuntu_version in ${ubuntu_versions}; do
  image_name="svalinn/test_dagmc-ci-ubuntu-${ubuntu_version}"
  docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} \
               -f CI/Dockerfile_0_base .
  if [ "$PUSH" = true ]; then docker push ${image_name}; fi
  for compiler in $compilers; do
    image_name="svalinn/dagmc-ci-ubuntu-${ubuntu_version}-${compiler}-ext"
    docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} \
                                  --build-arg COMPILER=${compiler} \
                 -f CI/Dockerfile_1_external_deps .
    if [ "$PUSH" = true ]; then docker push ${image_name}; fi
    for hdf5 in $hdf5_versions; do
      image_name="svalinn/dagmc-ci-ubuntu-${ubuntu_version}-${compiler}-ext-hdf5_${hdf5}"
      docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} \
                                    --build-arg COMPILER=${compiler} \
                                    --build-arg HDF5=${hdf5} \
                   -f CI/Dockerfile_2_hdf5 .
      if [ "$PUSH" = true ]; then docker push ${image_name}; fi
      for moab in $moab_versions; do
        image_name="svalinn/dagmc-ci-ubuntu-${ubuntu_version}-${compiler}-ext-hdf5_${hdf5}-moab_${moab}"
        docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} \
                                      --build-arg COMPILER=${compiler} \
                                      --build-arg HDF5=${hdf5} \
                                      --build-arg MOAB=${moab} \
                     -f CI/Dockerfile_3_moab .
        if [ "$PUSH" = true ]; then docker push ${image_name}; fi
      done
    done
  done
done

# Building image for houe keeping
docker build -t svalinn/test_dagmc-ci-ubuntu-18.04-housekeeping -f CI/Dockerfile_1_housekeeping .
if [ "$PUSH" = true ]; then docker push svalinn/test_dagmc-ci-ubuntu-18.04-housekeeping; fi
