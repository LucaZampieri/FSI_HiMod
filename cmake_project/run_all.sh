#!/bin/bash

./configure_fsiHimod.sh
make -C src-build -j 4
rm -r src-build/

