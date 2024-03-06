# Svalinn DAGMC Dockerfile

# build from repository root directory with command
# docker build -t dagmc -f CI/Dockerfile .

# Global arguments set before the first build stage are accessable by all build stages
ARG UBUNTU_VERSION=22.04
ARG HDF5_VERSION=hdf5_1_14_3
ARG MOAB_VERSION=5.5.1
# Setting GEANT4_VERSION=off will disable GEANT4 build
ARG GEANT4_VERSION=v11.2.1
# Setting DOUBLE_DOWN_VERSION=off will disable DOUBLE_DOWN build
ARG DOUBLE_DOWN_VERSION=v1.1.0
# Setting DOUBLE_DOWN_VERSION=off will disable EMBREE build 
ARG EMBREE_VERSION=v4.0.1

ARG BUILD_DIR=/root/build
ARG INSTALL_DIR=/opt

ARG COMPILER=gcc
ARG CI_JOBS=4

FROM ubuntu:${UBUNTU_VERSION} AS base

# Use bash as the default shell
SHELL ["/bin/bash", "-c"]

# timezone setup
ENV TZ=America/Chicago
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install Build Dependencies
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends \
        git \
        ca-certificates \
        make \
        cmake \
        g++ \
        clang \
        gfortran \
        python3 \
        libpython3-dev \
        zlib1g-dev \
        libeigen3-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*;

ARG BUILD_DIR
ARG INSTALL_DIR

ENV BUILD_DIR=${BUILD_DIR}
ENV INSTALL_DIR=${INSTALL_DIR}


# set default compiler to gcc
FROM base AS compiler-gcc

ENV CC=gcc
ENV CXX=g++


# set default compiler to clang
FROM base AS compiler-clang

ENV CC=clang
ENV CXX=clang++


# Install External Dependencies
FROM compiler-${COMPILER} AS external_deps

# accessing gloabl ARGs in build stage
ARG INSTALL_DIR 
ARG BUILD_DIR
ARG CI_JOBS
ARG CXX
ARG CC

# Build Embree
ARG DOUBLE_DOWN_VERSION
ARG EMBREE_VERSION

ENV EMBREE_INSTALL_DIR=${INSTALL_DIR}/Embree
ENV EMBREE_BUILD_DIR=${BUILD_DIR}/Embree

# Clone and install Embree
RUN if [ "${DOUBLE_DOWN_VERSION}" != "off" ]; then \
        mkdir -p ${EMBREE_BUILD_DIR}/build && \
        cd ${EMBREE_BUILD_DIR} && \
        git clone -b ${EMBREE_VERSION} --depth 1 https://github.com/embree/embree.git && \
        cd build && \
        cmake ../embree \
            -DCMAKE_INSTALL_PREFIX=${EMBREE_INSTALL_DIR} \
            -DEMBREE_TASKING_SYSTEM=INTERNAL \
            -DEMBREE_ISPC_SUPPORT=OFF \
            -DEMBREE_TUTORIALS=OFF \
            -DEMBREE_TBB_ROOT=/usr && \
        make -j${CI_JOBS} && \
        make -j${CI_JOBS} install && \
        cd && \
        rm -rf ${EMBREE_BUILD_DIR}; \
    fi


# Build Geant4
ARG GEANT4_VERSION

ENV GEANT4_BUILD_DIR=${BUILD_DIR}/geant4
ENV GEANT4_INSTALL_DIR=${INSTALL_DIR}/geant4

RUN if [ "${GEANT4_VERSION}" != "off" ]; then \
        mkdir -p ${GEANT4_BUILD_DIR}/build && \
        cd ${GEANT4_BUILD_DIR} && \
        git clone -b ${GEANT4_VERSION} --depth 1 https://gitlab.cern.ch/geant4/geant4.git && \
        cd build && \
        cmake ../geant4 \
            -DCMAKE_INSTALL_RPATH=${GEANT4_INSTALL_DIR}/lib \
            -DCMAKE_INSTALL_PREFIX=${GEANT4_INSTALL_DIR} \
            -DGEANT4_USE_SYSTEM_EXPAT=OFF \
            -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_CXX_COMPILER=${CXX} \
            -DCMAKE_C_COMPILER=${CC} \
            -DBUILD_STATIC_LIBS=ON && \
        make -j${CI_JOBS} && \
        make install && \
        cd && \
        rm -rf ${GEANT4_BUILD_DIR}; \
    fi


# Build HDF5
ARG HDF5_VERSION

ENV HDF5_BUILD_DIR=${BUILD_DIR}/hdf5
ENV HDF5_INSTALL_DIR=${INSTALL_DIR}/hdf5

RUN mkdir -p ${HDF5_BUILD_DIR}/build && \
    cd ${HDF5_BUILD_DIR} && \
    git clone -b ${HDF5_VERSION} --depth 1 https://github.com/HDFGroup/hdf5.git && \
    cd build && \
    cmake ../hdf5 \
        -DCMAKE_INSTALL_PREFIX=${HDF5_INSTALL_DIR} \
        -DBUILD_SHARED_LIBS=ON && \
    make -j${CI_JOBS} && \
    make install && \
    cd && \
    rm -rf ${HDF5_BUILD_DIR}


# Build MOAB
ARG MOAB_VERSION

ENV MOAB_BUILD_DIR=${BUILD_DIR}/moab
ENV MOAB_INSTALL_DIR=${INSTALL_DIR}/moab

RUN mkdir -p ${MOAB_BUILD_DIR}/build && \
    cd ${MOAB_BUILD_DIR} && \
    git clone -b ${MOAB_VERSION} --depth 1 https://bitbucket.org/fathomteam/moab.git && \
    cd build && \
    cmake ../moab \
        -DCMAKE_INSTALL_RPATH=${HDF5_INSTALL_DIR}/lib:${MOAB_INSTALL_DIR}/lib \
        -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_HDF5=ON \
        -DHDF5_ROOT=${HDF5_INSTALL_DIR} \
        -DCMAKE_INSTALL_PREFIX=${MOAB_INSTALL_DIR} \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_C_COMPILER=${CC} \
        -DENABLE_BLASLAPACK=OFF \
        -DBUILD_SHARED_LIBS=ON \
        -DENABLE_FORTRAN=OFF \
        -DENABLE_PYMOAB=OFF && \
    make -j${CI_JOBS} && \
    make install && \
    cd && \
    rm -rf ${MOAB_BUILD_DIR}


# Build Double Down
ARG DOUBLE_DOWN_VERSION

ENV DOUBLE_DOWN_BUILD_DIR=${BUILD_DIR}/double-down
ENV DOUBLE_DOWN_INSTALL_DIR=${INSTALL_DIR}/double-down

# Clone and install Double-Down
# Performed after moab and embree install as double-down requires moab and embree
RUN if [ "${DOUBLE_DOWN_VERSION}" != "off" ]; then \
        mkdir -p ${DOUBLE_DOWN_BUILD_DIR}/build && \
        cd ${DOUBLE_DOWN_BUILD_DIR} && \
        git clone -b ${DOUBLE_DOWN_VERSION} --depth 1 https://github.com/pshriwise/double-down.git && \
        cd build && \
        cmake ../double-down \
            -DMOAB_DIR=${MOAB_INSTALL_DIR} \
            -DCMAKE_INSTALL_PREFIX=${DOUBLE_DOWN_INSTALL_DIR} \
            -DEMBREE_DIR=${EMBREE_INSTALL_DIR} && \
        make -j${CI_JOBS} && \
        make -j${CI_JOBS} install && \
        cd && \
        rm -rf ${DOUBLE_DOWN_BUILD_DIR}; \
    fi


# Build DAGMC
FROM external_deps AS dagmc

# accessing gloabl ARGs in build stage
ARG INSTALL_DIR
ARG BUILD_DIR
ARG CXX
ARG CC
ARG GEANT4_VERSION
ARG DOUBLE_DOWN_VERSION

ENV DAGMC_BUILD_DIR=${BUILD_DIR}/dagmc
ENV DAGMC_INSTALL_DIR=${INSTALL_DIR}/dagmc

# Copy the enitre git repo into the dockerfile to ensure the submodules can also be found
COPY . ${DAGMC_BUILD_DIR}

RUN mkdir -p ${DAGMC_BUILD_DIR}/build && \
    cd ${DAGMC_BUILD_DIR} && \
    git submodule update --init && \
    cd build && \
    cmake ../ \
        -DMOAB_DIR=${MOAB_INSTALL_DIR} \
        -DBUILD_GEANT4=$([ "${GEANT4_VERSION}" != "off" ] && echo "ON" || echo "OFF") \
        -DGEANT4_DIR=${GEANT4_INSTALL_DIR} \
        -DBUILD_CI_TESTS=ON \
        # reg tests disabled see issue https://github.com/svalinn/DAGMC/issues/858
        -DBUILD_MW_REG_TESTS=OFF \
        -DBUILD_STATIC_EXE=OFF \
        -DBUILD_STATIC_LIBS=OFF \
        -DCMAKE_C_COMPILER=${CC} \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_Fortran_COMPILER=gfortran \
        -DCMAKE_INSTALL_PREFIX=${DAGMC_INSTALL_DIR} \
        -DDOUBLE_DOWN=$([ "${DOUBLE_DOWN_VERSION}" != "off" ] && echo "ON" || echo "OFF") \
        -Ddd_ROOT=${DOUBLE_DOWN_INSTALL_DIR} && \
    make -j${CI_JOBS} && \
    make install


# Test DAGMC
FROM dagmc AS dagmc_test

ENV DAGMC_BUILD_DIR=${BUILD_DIR}/dagmc
ENV DAGMC_INSTALL_DIR=${INSTALL_DIR}/dagmc

# Test DAGMC executables
RUN cd ${DAGMC_BUILD_DIR}/build && \
    PATH=${DAGMC_INSTALL_DIR}/bin:$PATH CTEST_OUTPUT_ON_FAILURE=1 make test

# clean out config test directory for next build
RUN cd ${DAGMC_BUILD_DIR} && \
    git clean -dxf .

# Test DAGMC CMake configuration file
RUN cd ${DAGMC_BUILD_DIR}/cmake/test_config && \
    cmake . -DDAGMC_ROOT=${DAGMC_INSTALL_DIR} && \
    CTEST_OUTPUT_ON_FAILURE=1 make all test


# Release image
FROM ubuntu:${UBUNTU_VERSION} AS final

ARG INSTALL_DIR

COPY --from=dagmc_test ${INSTALL_DIR} ${INSTALL_DIR}

ENV PATH=${DAGMC_INSTALL_DIR}/bin:$PATH
ENV PATH=${MOAB_INSTALL_DIR}/bin:$PATH
ENV PATH=${EMBREE_INSTALL_DIR}/bin:$PATH
ENV PATH=${DOUBLE_DOWN_INSTALL_DIR}/bin:$PATH
ENV PATH=${HDF5_INSTALL_DIR}/bin:$PATH
ENV PATH=${GEANT4_INSTALL_DIR}/bin:$PATH

WORKDIR ${INSTALL_DIR}

ENTRYPOINT ["/bin/bash"]