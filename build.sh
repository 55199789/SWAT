#!/bin/bash
set -e

# rm -rf bin build
cmake -DRECORD_LENGTH_IN_BYTE=4096 -S . -B build 

START=$(date +%s)
cmake --build build -j
END=$(date +%s)
echo "Total build time (real) = $(( $END - $START )) seconds"