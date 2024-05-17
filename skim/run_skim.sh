#!/bin/bash

# Check if the number of arguments is as expected
if [ "$#" -ne 4 ]; then
  echo "Usage: ${0} <INPUT_FILE> <RUN_NUMBER> <OUTPUT_FOLDER> <isTB>"
  exit 1
fi

INPUT_FILE=$1
RUN_NUMBER=$2
OUTPUT_FOLDER=$3
isTB=$4


# Define an array of commands
root -l -b -x <<EOF
.L libs/SciFiPlaneView.cpp
.L libs/USPlaneView.cpp
.L skim.cpp
skim("$INPUT_FILE", $RUN_NUMBER, "$OUTPUT_FOLDER", $isTB)
.q
EOF
