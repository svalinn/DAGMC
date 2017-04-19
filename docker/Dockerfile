FROM ubuntu:16.04

# Use bash as the default shell
RUN ln -snf /bin/bash /bin/sh

# Add the Ubuntu 14.04 sources
RUN echo ""                                                          >> /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu/ trusty universe"     >> /etc/apt/sources.list && \
    echo "deb-src http://archive.ubuntu.com/ubuntu/ trusty universe" >> /etc/apt/sources.list

RUN apt-get -y update
RUN apt-get -y upgrade
RUN apt-get -y dist-upgrade

# Get some common utils
RUN apt-get -y install apt-utils build-essential software-properties-common

# Add the GCC repo
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get -y update

# Get all the GNU compilers
RUN apt-get -y install gcc-4.4 g++-4.4 gfortran-4.4 \
                       gcc-4.6 g++-4.6 gfortran-4.6 \
                       gcc-4.7 g++-4.7 gfortran-4.7 \
                       gcc-4.8 g++-4.8 gfortran-4.8 \
                       gcc-4.9 g++-4.9 gfortran-4.9 \
                       gcc-5   g++-5   gfortran-5 \
                       gcc-6   g++-6   gfortran-6

# Get all the Clang compilers
RUN apt-get -y install clang-3.5 clang-3.6 clang-3.7 clang-3.8

# Setup easy switching between different GNU compiler versions
# To use 4.4:    echo 1 | update-alternatives --config gcc
# To use 4.6:    echo 2 | update-alternatives --config gcc
# To use 4.7:    echo 3 | update-alternatives --config gcc
# To use 4.8:    echo 4 | update-alternatives --config gcc
# To use 4.9:    echo 5 | update-alternatives --config gcc
# To use 5:      echo 6 | update-alternatives --config gcc
# To use 6:      echo 7 | update-alternatives --config gcc
# To use latest: echo 0 | update-alternatives --config gcc
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.4 440 --slave /usr/bin/g++ g++ /usr/bin/g++-4.4 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-4.4 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.6 460 --slave /usr/bin/g++ g++ /usr/bin/g++-4.6 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-4.6 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.7 470 --slave /usr/bin/g++ g++ /usr/bin/g++-4.7 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-4.7 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 480 --slave /usr/bin/g++ g++ /usr/bin/g++-4.8 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-4.8 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 490 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9 --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-4.9 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5   500 --slave /usr/bin/g++ g++ /usr/bin/g++-5   --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-5 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6   600 --slave /usr/bin/g++ g++ /usr/bin/g++-6   --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-6

# Setup easy switching between different Clang compiler versions
# To use 3.5:    echo 1 | update-alternatives --config clang
# To use 3.6:    echo 2 | update-alternatives --config clang
# To use 3.7:    echo 3 | update-alternatives --config clang
# To use 3.8:    echo 4 | update-alternatives --config clang
# To use latest: echo 0 | update-alternatives --config clang
RUN update-alternatives --install /usr/bin/clang clang /usr/bin/clang-3.5 350 --slave /usr/bin/clang++ clang++ /usr/bin/clang++-3.5 && \
    update-alternatives --install /usr/bin/clang clang /usr/bin/clang-3.6 360 --slave /usr/bin/clang++ clang++ /usr/bin/clang++-3.6 && \
    update-alternatives --install /usr/bin/clang clang /usr/bin/clang-3.7 370 --slave /usr/bin/clang++ clang++ /usr/bin/clang++-3.7 && \
    update-alternatives --install /usr/bin/clang clang /usr/bin/clang-3.8 380 --slave /usr/bin/clang++ clang++ /usr/bin/clang++-3.8

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
RUN bash /root/etc/build_hdf5.sh gcc-5
RUN bash /root/etc/build_hdf5.sh gcc-6
RUN bash /root/etc/build_hdf5.sh clang-3.8

# Build MOAB
COPY build_moab.sh /root/etc/
RUN bash /root/etc/build_moab.sh gcc-4.8   4.9.2
RUN bash /root/etc/build_moab.sh gcc-5     4.9.2
RUN bash /root/etc/build_moab.sh gcc-6     4.9.2
RUN bash /root/etc/build_moab.sh clang-3.8 4.9.2

# Build Geant4
COPY build_geant4.sh /root/etc/
RUN bash /root/etc/build_geant4.sh gcc-4.8
RUN bash /root/etc/build_geant4.sh gcc-5
RUN bash /root/etc/build_geant4.sh gcc-6
RUN bash /root/etc/build_geant4.sh clang-3.8

CMD ["bash"]
