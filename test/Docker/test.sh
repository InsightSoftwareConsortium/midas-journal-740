#!/bin/bash

# This is a script to build the modules and run the test suite in the base
# Docker container.

set -x
set -e

cd /usr/src/ITKCuberille-build

cmake \
  -G Ninja \
  -DITK_DIR:PATH=/usr/src/ITK-build \
  -DITK_USE_SYSTEM_SWIG:BOOL=ON \
  -DSWIG_EXECUTABLE:FILEPATH=/usr/bin/swig3.0 \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DBUILDNAME:STRING=External-Cuberille \
    /usr/src/ITKCuberille
ctest -VV -D Experimental
