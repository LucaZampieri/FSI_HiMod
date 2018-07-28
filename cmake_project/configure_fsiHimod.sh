#!/bin/bash
#
# Andrea Bartezzaghi, 14-Nov-2017 modified by luca zampieri
#

# check arguments
CONFIG_FILE="fsiHimod_config.sh"

# TO LOAD AND DO LIBRAIRIES :D
#LIBRARY_CONFIG_FILE=libs/config.sh
# load config
#source $LIBRARY_CONFIG_FILE

# stop on errors
set -e

# provide default values
POSTFIX=


# load parameters from custom config file
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: could not find $CONFIG_FILE with custom configuration!"
    exit 1
fi
source $CONFIG_FILE

# directories
SCRIPT_PATH=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo "Script path is: "
echo ${SCRIPT_PATH}
SOURCE_DIR=${SCRIPT_PATH}/${FSIHIMOD}
BUILD_DIR=${SOURCE_DIR}${POSTFIX}-build
INSTALL_DIR=${SOURCE_DIR}${POSTFIX}-install

#CMAKE_MODULE_PATH =


# enter build directory
mkdir -pv ${BUILD_DIR} && cd ${BUILD_DIR}

#export CMAKE_BIN=${CMAKE_INSTALL_DIR}/bin/cmake #/usr/local/bin/cmake
export CMAKE_BIN=cmake
#CMAKE_PREFIX_PATH
echo "prefix path is:" ${CMAKE_PREFIX_PATH}
# configure
#cmake -D CMAKE_PREFIX_PATH=/
cmake -G "Unix Makefiles" ..
      #-Wno-deprecated
      #-Wno-dev
      #${SOURCE_DIR}

# build
#make -j${NUM_PROC}
