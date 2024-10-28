#!/bin/bash

HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/


# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

# I3BASE=/data/user/enpaudel/icecube_software/icetray_main
I3BASE=/data/user/enpaudel/icecube_software/icetray_filter/
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh

GCD=/cvmfs/icecube.opensciencegrid.org/data/i3-test-data-svn/trunk/2023data/Level2_IC86.2023_data_Run00138821_84_781_GCD_withSLCcal2.i3.zst

FILTER_ICETOP_PY=$I3SRC/offline_filterscripts/resources/scripts/Filter_IceTop.py
# INPUT=/home/enpaudel/dataExp/run2023/IceTopTrig/PFFilt_PhysicsFiltering_Run00138615_Subrun00000000_00000000.i3.gz
INPUT=$1

filepath="/path/to/your/file.txt"
filename=$(basename "$INPUT")
echo $filename
OUTPUT=/home/enpaudel/dataExp/run2023/new_filter_test/filtered/$filename
FLAGS="-i $INPUT -o $OUTPUT -g $GCD -u --it-inclined-split"


$ICETRAY_ENV ${FILTER_ICETOP_PY} $FLAGS