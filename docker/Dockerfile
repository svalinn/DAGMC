FROM ubuntu:16.04

# Use bash as the default shell
RUN ln -snf /bin/bash /bin/sh

RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get -y dist-upgrade

# Get some common utils
RUN apt-get -y install build-essential software-properties-common

# Add the GCC repo
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get -y update

# Get the GNU compilers
RUN apt-get -y install gcc-4.8 g++-4.8 gfortran-4.8
RUN apt-get -y install gcc-5   g++-5   gfortran-5
RUN apt-get -y install gcc-6   g++-6   gfortran-6
RUN apt-get -y install gcc-7   g++-7   gfortran-7

# Get the Clang compilers
RUN apt-get -y install clang-4.0
RUN apt-get -y install clang-5.0

# Get some more specific utils
RUN apt-get -y install autoconf cmake gfortran git libtool unzip wget
RUN apt-get -y install libblas-dev libexpat-dev libhdf5-dev liblapack-dev libz-dev

# Copy the environment variable and build instruction files
RUN mkdir -p /root/etc/
COPY *.env /root/etc/

# Build HDF5
COPY build_hdf5.sh /root/etc/
RUN bash /root/etc/build_hdf5.sh gcc-4.8
RUN bash /root/etc/build_hdf5.sh gcc-5
RUN bash /root/etc/build_hdf5.sh gcc-6
RUN bash /root/etc/build_hdf5.sh gcc-7
RUN bash /root/etc/build_hdf5.sh clang-4.0
RUN bash /root/etc/build_hdf5.sh clang-5.0

# Build Geant4
COPY build_geant4.sh /root/etc/
RUN bash /root/etc/build_geant4.sh gcc-7
RUN bash /root/etc/build_geant4.sh clang-5.0

CMD ["bash"]
