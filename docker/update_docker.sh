#!/bin/bash

set -e

ubuntu_versions="16.04 18.04"
for ubuntu_version in ${ubuntu_versions}; do
  docker build -f Dockerfile-Ubuntu-${ubuntu_version} .
  docker tag `docker images -q | head -1` ljacobson64/dagmc-ci-ubuntu-${ubuntu_version}:latest
  docker push ljacobson64/dagmc-ci-ubuntu-${ubuntu_version}:latest
done
