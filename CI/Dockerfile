# Svalinn DAGMC Dockerfile

ARG UBUNTU_VERSION=18.04

ARG ci_jobs=8
ARG MOAB_BRANCH=5.4.0
ARG EMBREE_BRANCH='v3.6.1'
ARG geant4_version=10.5.1

ARG HDF5_VERSION=1.10.4
# currently major version must be manually set to match HDF5 version
ARG HDF5_VERSION_major=1.10

ARG build_dir=/root/build_dir
ARG install_dir=/root/opt

# clang and gcc are options for the compiler but to use them CC and CXX must be set
# if clang then change cc to clang and cxx to clang++
# if gcc then change cc to gcc and cxx to g++
ARG CC=gcc
ARG CXX=g++

FROM ubuntu:${UBUNTU_VERSION} AS base

# Use bash as the default shell
SHELL ["/bin/bash", "-c"]
# Ubuntu Setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Update core packages
RUN apt-get -y update; \
    apt-get -y install autoconf \
                       clang \
                       cmake \
                       g++ \
                       gcc \
                       gfortran \
                       libhdf5-dev \
                       libtool \
                       libeigen3-dev\
                       python3-numpy \
                       python3 \
                       python3-pip \
                       python3-setuptools \
                       python3-dev \
                       libpython3-dev \
                       wget \
                       software-properties-common; \
    add-apt-repository ppa:git-core/ppa; \
    apt-get update; \
    apt-get install -y git; \
    update-alternatives --install /usr/bin/python python /usr/bin/python3 10; \
    update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 10; \
    pip install cython;


FROM base as external_deps

ARG build_dir
ARG install_dir 
ARG ci_jobs
ARG geant4_version
ARG CC
ARG CXX
ARG EMBREE_BRANCH

ENV geant4_basename=geant4-v${geant4_version}
ENV geant4_tarball=${geant4_basename}.tar.gz
ENV geant4_shasum=2397eb859dc4de095ff66059d8bda9f060fdc42e10469dd7890946293eeb0e39
ENV geant4_build_dir=${build_dir}/geant4
ENV geant4_install_dir=${install_dir}/geant4

RUN mkdir -p ${geant4_build_dir}/bld
RUN cd ${geant4_build_dir} && \
    wget https://gitlab.cern.ch/geant4/geant4/-/archive/v${geant4_version}/${geant4_tarball} --no-check-certificate && \
    tar -xzf ${geant4_tarball} && \
    cd bld && \
    cmake ../${geant4_basename} -DBUILD_STATIC_LIBS=ON \
                -DGEANT4_USE_SYSTEM_EXPAT=OFF \
                -DCMAKE_BUILD_TYPE=Release \
                -DCMAKE_C_COMPILER=${CC} \
                -DCMAKE_CXX_COMPILER=${CXX} \
                -DCMAKE_INSTALL_RPATH=${geant4_install_dir}/lib \
                -DCMAKE_INSTALL_PREFIX=${geant4_install_dir} && \
    make -j${ci_jobs} && \
    make install && \
    cd && \
    rm -rf ${geant4_build_dir}

# these two double down ENVs don't appear to be used yet
# However they were in orginal dockerfile, so they have been kept
# TODO add double down to the dockerfile
ENV double_down_build_dir=${build_dir}/double-down/
ENV double_down_install_dir=${install_dir}/double-down/

ENV embree_install_dir=${install_dir}/embree

RUN mkdir EMBREE && \
    cd EMBREE && \
    git clone -b ${EMBREE_BRANCH} https://github.com/embree/embree && \
    mkdir build && \
    cd build && \
    cmake ../embree -DCMAKE_INSTALL_PREFIX=${embree_install_dir} \
        -DEMBREE_ISPC_SUPPORT=OFF \
        -DEMBREE_TASKING_SYSTEM=INTERNAL \
        -DEMBREE_TUTORIALS=OFF \
        -DEMBREE_TBB_ROOT=/usr && \
    make -j${ci_jobs} && make -j${ci_jobs} install && \
    rm -rf $HOME/EMBREE/embree


FROM external_deps AS hdf5

ARG build_dir
ARG install_dir
ARG ci_jobs
ARG HDF5_VERSION
ARG HDF5_VERSION_major
ARG CC
ARG CXX

ENV hdf5_build_dir=${build_dir}/hdf5
ENV hdf5_install_dir=${install_dir}/hdf5

RUN mkdir -p ${hdf5_build_dir}/bld && \
    cd ${hdf5_build_dir} && \
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION_major}/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz && \
    tar -xzf hdf5-${HDF5_VERSION}.tar.gz && \
    cd bld && \
    ../hdf5-${HDF5_VERSION}/configure --enable-shared \
        --prefix=${hdf5_install_dir} \
        CC=${CC} \
        CXX=${CXX} && \
    make -j${ci_jobs} && \
    make install && \
    cd && \
    rm -rf ${hdf5_build_dir}


FROM hdf5 AS moab

ARG ci_jobs
ARG MOAB_BRANCH
ARG CC
ARG CXX
ARG build_dir
ARG install_dir

# Set MOAB env variable
ENV moab_build_dir=${build_dir}/moab
ENV moab_install_dir=${install_dir}/moab

RUN mkdir -p ${moab_build_dir}/bld
RUN cd ${moab_build_dir} && \
    git clone -b ${MOAB_BRANCH} --depth 1 https://bitbucket.org/fathomteam/moab && \
    cd bld && \
    cmake ../moab -DENABLE_HDF5=ON -DHDF5_ROOT=${hdf5_install_dir} \
                -DENABLE_BLASLAPACK=OFF \
                -DENABLE_FORTRAN=OFF \
                -DCMAKE_INSTALL_PREFIX=${moab_install_dir} \
                -DCMAKE_C_COMPILER=${CC} \
                -DCMAKE_CXX_COMPILER=${CXX} \
                -DBUILD_SHARED_LIBS=OFF && \
    make -j${ci_jobs} && \
    make install && \
    rm -rf * && \
    cmake ../moab -DENABLE_HDF5=ON -DHDF5_ROOT=${hdf5_install_dir} \
                -DENABLE_PYMOAB=ON \
                -DENABLE_BLASLAPACK=OFF \
                -DENABLE_FORTRAN=OFF \
                -DCMAKE_INSTALL_PREFIX=${moab_install_dir} \
                -DCMAKE_C_COMPILER=${CC} \
                -DCMAKE_CXX_COMPILER=${CXX} \
                -DBUILD_SHARED_LIBS=ON \
                -DCMAKE_INSTALL_RPATH=${hdf5_install_dir}/lib:${moab_install_dir}/lib && \
    # ensure that the installation directory for python packages exist
    for d in $(echo $PYTHONPATH | tr : '\n'); do mkdir -p $d; done && \
    make -j${ci_jobs} && \
    make install && \
    cd && \
    rm -rf ${moab_build_dir}

ENV dagmc_build_dir=${build_dir}/DAGMC
ENV dagmc_install_dir=${install_dir}/DAGMC

ENV dagmc_build_dir_shared=${dagmc_build_dir}-shared
ENV dagmc_install_dir_shared=${dagmc_install_dir}-shared

ENV dagmc_build_dir_static=${dagmc_build_dir}-static
ENV dagmc_install_dir_static=${dagmc_install_dir}-static

