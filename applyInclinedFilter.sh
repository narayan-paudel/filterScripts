#!/bin/bash

HERE=$(dirname $(realpath -s $0))
BASEDIR=$HERE/

echo $BASEDIR


# eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

# I3BASE=/data/user/enpaudel/icecube_software/icetray_main
I3BASE=/data/user/enpaudel/icecube_software/icetray_filter/
I3SRC=$I3BASE/src
I3BUILD=$I3BASE/build
ICETRAY_ENV=$I3BUILD/env-shell.sh

GCD=/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz

FILT_INCL_PY=$BASEDIR/applyInclinedFilter.py
# INPUT=/home/enpaudel/dataExp/run2023/IceTopTrig/PFFilt_PhysicsFiltering_Run00138615_Subrun00000000_00000000.i3.gz
INPUT=$1

filename=$(basename "$INPUT")
echo $filename
OUTPUT=/data/sim/IceTop/2023/generated/untriggered/filterStudy/filtered/inclinedSplit/$filename

FLAGS="-i $INPUT -o $OUTPUT -g $GCD"

echo $INPUT $OUTPUT $GCD 
echo $FLAGS


$ICETRAY_ENV ${FILT_INCL_PY} $FLAGS