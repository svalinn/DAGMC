ARG UBUNTU_VERSION=18.04
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


# Copy scripts to docker image
RUN mkdir -p /root/etc/
COPY CI/ /root/etc/CI
#TODO move sh file contents into this Dockerfile
ENV docker_env=/root/etc/CI/env.sh

ENV build_dir=/root/build_dir
ENV install_dir=/root/opt


FROM base as external_deps

#setting the COMPILER variable
ARG COMPILER=gcc
ENV COMPILER=${COMPILER}

# Set Geant4 env variable
ENV geant4_build_dir=${build_dir}/geant4
ENV geant4_install_dir=${install_dir}/geant4

# Build Geant4
#TODO move sh file contents into this Dockerfile
RUN /root/etc/CI/docker/build_geant4.sh

ENV double_down_build_dir=${build_dir}/double-down/
ENV double_down_install_dir=${install_dir}/double-down/

# Build Embree
#TODO move sh file contents into this Dockerfile
RUN /root/etc/CI/docker/build_embree.sh

FROM external_deps AS hdf5

# Set HDF5 env variable
ENV hdf5_build_dir=${build_dir}/hdf5
ENV hdf5_install_dir=${install_dir}/hdf5

# Build HDF5
# HDF5 argument possible value: 1.10.4 or system
ARG HDF5=1.10.4
ENV HDF5_VERSION=${HDF5}
#TODO move sh file contents into this Dockerfile
RUN /root/etc/CI/docker/build_hdf5.sh

FROM hdf5 AS moab

# Set MOAB env variable
ENV moab_build_dir=${build_dir}/moab
ENV moab_install_dir=${install_dir}/moab

ARG MOAB=5.3.0
ENV MOAB_VERSION ${MOAB}
#TODO move sh file contents into this Dockerfile
RUN if [ "${MOAB_VERSION}" != "master" ] && [ "${MOAB_VERSION}" != "develop" ]; then \
        /root/etc/CI/docker/build_moab.sh; \
    fi;

