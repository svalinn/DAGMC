ARG UBUNTU_VERSION=18.04
FROM ubuntu:${UBUNTU_VERSION}

# Use bash as the default shell
SHELL ["/bin/bash", "-c"]

# Update core packages
ARG UBUNTU_VERSION=18.04
RUN apt-get -y update; \
    apt-get -y install autoconf clang cmake g++ gcc gfortran git libblas-dev \
                       libhdf5-dev liblapack-dev libpython2.7-dev libtool libeigen3-dev\
                       python-numpy python-pip python-setuptools wget; \
    if [ "${UBUNTU_VERSION}" = "18.04" ]; then \
        pip install sphinx; \
        apt-get -y install astyle; \
    fi; \
    pip install cython

# Copy scripts to docker image
RUN mkdir -p /root/etc/
COPY CI/ /root/etc/CI
ENV docker_env=/root/etc/CI/env.sh

# Build HDF5 then Geant4 then MOAB 5.1.0
RUN COMPILER=gcc   HDF5_VERSION=1.10.4 /root/etc/CI/docker/build_hdf5.sh; \
    COMPILER=clang HDF5_VERSION=1.10.4 /root/etc/CI/docker/build_hdf5.sh; \
    COMPILER=gcc   /root/etc/CI/docker/build_geant4.sh; \
    COMPILER=clang /root/etc/CI/docker/build_geant4.sh; \
    COMPILER=gcc   HDF5_VERSION=system MOAB_VERSION=5.1.0 /root/etc/CI/docker/build_moab.sh; \
    COMPILER=clang HDF5_VERSION=system MOAB_VERSION=5.1.0 /root/etc/CI/docker/build_moab.sh; \
    COMPILER=gcc   HDF5_VERSION=1.10.4 MOAB_VERSION=5.1.0 /root/etc/CI/docker/build_moab.sh; \
    COMPILER=clang HDF5_VERSION=1.10.4 MOAB_VERSION=5.1.0 /root/etc/CI/docker/build_moab.sh
