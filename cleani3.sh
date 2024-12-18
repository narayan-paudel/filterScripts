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

CLEAN_I3_PY=$BASEDIR/cleani3.py

# INPUT_FOLDER=/data/user/enpaudel/triggerStudy/simFiles/dataSetClean1_6/
INPUT_LIST="$@"
OUTPUT_FOLDER=/data/sim/IceTop/2023/generated/untriggered/IceTop7HGSim/
for ifile in $INPUT_LIST;do
	# echo $(basename $ifile)
	outfile=$OUTPUT_FOLDER$(basename $ifile)
	$ICETRAY_ENV ${CLEAN_I3_PY} -i $ifile -o $outfile
	echo $outfile
	# echo $ifile
done
 # to run; ./cleani3.sh /data/user/enpaudel/triggerStudy/simFiles/dataSetClean1_6/*.i3*
 # to run; ./cleani3.sh /data/user/enpaudel/triggerStudy/simFiles/dataSetClean1_6/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz
