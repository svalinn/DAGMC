FROM ubuntu:16.04

# Use bash as the default shell
RUN ln -snf /bin/bash /bin/sh

# Add the Ubuntu 14.04 sources
RUN echo ""                                                                   >> /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu/ trusty universe"              >> /etc/apt/sources.list && \
    echo "deb-src http://archive.ubuntu.com/ubuntu/ trusty universe"          >> /etc/apt/sources.list

# Add the Clang 3.9 and 4.0 sources
RUN echo ""                                                                   >> /etc/apt/sources.list && \
    echo "deb http://apt.llvm.org/xenial/ llvm-toolchain-xenial-3.9 main"     >> /etc/apt/sources.list && \
    echo "deb-src http://apt.llvm.org/xenial/ llvm-toolchain-xenial-3.9 main" >> /etc/apt/sources.list && \
    echo "deb http://apt.llvm.org/xenial/ llvm-toolchain-xenial-4.0 main"     >> /etc/apt/sources.list && \
    echo "deb-src http://apt.llvm.org/xenial/ llvm-toolchain-xenial-4.0 main" >> /etc/apt/sources.list

RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get -y dist-upgrade

# Get some common utils
RUN apt-get -y install apt-utils build-essential software-properties-common

# Add the GCC repo
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get -y update

# Get the GNU compilers
RUN apt-get -y install gcc-4.4 g++-4.4 gfortran-4.4 \
                       gcc-4.6 g++-4.6 gfortran-4.6 \
                       gcc-4.7 g++-4.7 gfortran-4.7 \
                       gcc-4.8 g++-4.8 gfortran-4.8 \
                       gcc-4.9 g++-4.9 gfortran-4.9 \
                       gcc-5   g++-5   gfortran-5 \
                       gcc-6   g++-6   gfortran-6

# Get the Clang compilers
RUN apt-get -y install clang-3.5 clang-3.6 clang-3.7 clang-3.8
RUN apt-get -y install clang-3.9 clang-4.0 --allow-unauthenticated

# Get some more specific utils
RUN apt-get -y install autoconf cmake curl cython emacs expat git htop \
                       libtool man unzip vim wget
RUN apt-get -y install libblas-dev libexpat1-dev libhdf5-dev liblapack-dev
RUN apt-get -y install python-jinja2 python-matplotlib python-nose \
                       python-numpy python-setuptools python-scipy python-tables

# Copy the environment variable and build instruction files
RUN mkdir -p /root/etc/
COPY *.env /root/etc/

# Build HDF5
COPY build_hdf5.sh /root/etc/
RUN bash /root/etc/build_hdf5.sh gcc-4.8
RUN bash /root/etc/build_hdf5.sh gcc-4.9
RUN bash /root/etc/build_hdf5.sh gcc-5
RUN bash /root/etc/build_hdf5.sh gcc-6
RUN bash /root/etc/build_hdf5.sh clang-3.8
RUN bash /root/etc/build_hdf5.sh clang-3.9
RUN bash /root/etc/build_hdf5.sh clang-4.0

# Build Geant4
COPY build_geant4.sh /root/etc/
RUN bash /root/etc/build_geant4.sh gcc-4.8   10.03.p01
RUN bash /root/etc/build_geant4.sh gcc-4.9   10.03.p01
RUN bash /root/etc/build_geant4.sh gcc-5     10.03.p01
RUN bash /root/etc/build_geant4.sh gcc-6     10.03.p01
RUN bash /root/etc/build_geant4.sh clang-3.8 10.03.p01
RUN bash /root/etc/build_geant4.sh clang-3.9 10.03.p01
RUN bash /root/etc/build_geant4.sh clang-4.0 10.03.p01

CMD ["bash"]
