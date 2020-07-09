#!/bin/bash

set -ex

ubuntu_versions="16.04" # 18.04 20.04"
compilers="gcc" # clang"
hdf5_versions="1.10.4"
moab_versions="5.1.0"
for ubuntu_version in ${ubuntu_versions}; do
  image_name="svalinn/test_dagmc-ci-ubuntu-${ubuntu_version}"
  docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} \
               -f CI/Dockerfile_0_base .
  docker push ${image_name}
  for compiler in $compilers; do
    image_name="svalinn/test_dagmc-ci-ubuntu-${ubuntu_version}-${compiler}-g4"
    docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} \
                                  --build-arg COMPILER=${compiler} \
                 -f CI/Dockerfile_1_g4 . --no-cache
    docker push ${image_name}
    for hdf5 in $hdf5_versions; do
      image_name="svalinn/test_dagmc-ci-ubuntu-${ubuntu_version}-${compiler}-g4-hdf5_${hdf5}"
      docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} \
                                    --build-arg COMPILER=${compiler} \
                                    --build-arg HDF5=${hdf5} \
                   -f CI/Dockerfile_2_hdf5 . --no-cache
      docker push ${image_name}
      for moab in $moab_versions; do
        image_name="svalinn/test_dagmc-ci-ubuntu-${ubuntu_version}-${compiler}-g4-hdf5_${hdf5}-moab_${moab}"
        docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} \
                                      --build-arg COMPILER=${compiler} \
                                      --build-arg HDF5=${hdf5} \
                                      --build-arg MOAB=${moab} \
                     -f CI/Dockerfile_2_hdf5 . --no-cache
        docker push ${image_name}
      done
    done
  done
done
