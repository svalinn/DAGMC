ARG UBUNTU_VERSION=18.04
FROM ubuntu:${UBUNTU_VERSION}

# Use bash as the default shell
SHELL ["/bin/bash", "-c"]

# Update core packages
RUN apt-get -y update

# Get some apt packages
RUN apt-get -y install autoconf clang cmake g++ gcc gfortran git libblas-dev \
                       libhdf5-dev liblapack-dev libpython2.7-dev libtool \
                       python-numpy python-pip python-setuptools wget

# Get astyle
RUN if [ "${UBUNTU_VERSION}" == "18.04" ]; then \
      apt-get -y install astyle; \
    else \
      astyle_deb=astyle_3.1-1ubuntu2_amd64.deb && \
      wget https://launchpad.net/ubuntu/+source/astyle/3.1-1ubuntu2/+build/14532685/+files/${astyle_deb} && \
      dpkg -i ${astyle_deb} && \
      rm -f ${astyle_deb}; \
    fi

# Get some pip packages
RUN pip install cython sphinx

# Copy scripts to docker image
RUN mkdir -p /root/etc/
COPY build_hdf5.sh build_geant4.sh build_moab.sh env.sh /root/etc/
ENV docker_env=/root/etc/env.sh

# Build HDF5
RUN COMPILER=gcc   HDF5_VERSION=1.10.4 /root/etc/build_hdf5.sh
RUN COMPILER=clang HDF5_VERSION=1.10.4 /root/etc/build_hdf5.sh

# Build Geant4
RUN COMPILER=gcc   /root/etc/build_geant4.sh
RUN COMPILER=clang /root/etc/build_geant4.sh

# Build MOAB 5.1.0
RUN COMPILER=gcc   HDF5_VERSION=system MOAB_VERSION=5.1.0 /root/etc/build_moab.sh
RUN COMPILER=clang HDF5_VERSION=system MOAB_VERSION=5.1.0 /root/etc/build_moab.sh
RUN COMPILER=gcc   HDF5_VERSION=1.10.4 MOAB_VERSION=5.1.0 /root/etc/build_moab.sh
RUN COMPILER=clang HDF5_VERSION=1.10.4 MOAB_VERSION=5.1.0 /root/etc/build_moab.sh
