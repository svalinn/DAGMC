#
# Ubuntu Dockerfile
#
# https://github.com/dockerfile/ubuntu
#

# pull the base image
FROM ubuntu:14.04

# install
RUN \
    apt-get update && \ 
    apt-get -y upgrade && \
    apt-get install -y build-essential && \
    apt-get install -y software-properties-common && \
    apt-get install -y cmake curl git htop man unzip vim wget && \
    apt-get install -y libhdf5-dev gfortran g++ libblas-dev liblapack-dev && \
    apt-get install -y expat libexpat1-dev
# install geant4
RUN \ 
    wget http://geant4.cern.ch/support/source/geant4.10.00.p02.tar.gz && \
    tar -zxf geant4.10.00.p02.tar.gz && \
    cd geant4.10.00.p02 && \
    mkdir bld && \
    cd bld && \
    cmake ../. -DCMAKE_INSTALL_PREFIX=/root/geant4.10.00.p02 && \
    make -j && \
    make install 

# Add files.
#ADD root/.bashrc /root/.bashrc
#A DD root/.gitconfig /root/.gitconfig
#ADD root/.scripts /root/.scripts

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]
    