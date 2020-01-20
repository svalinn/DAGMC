#!/bin/bash

set -e

ubuntu_versions="16.04 18.04"
for ubuntu_version in ${ubuntu_versions}; do
  image_name="svalinn/dagmc-ci-ubuntu-${ubuntu_version}:latest"
  docker build -t ${image_name} --build-arg UBUNTU_VERSION=${ubuntu_version} .
  docker push ${image_name}
done
