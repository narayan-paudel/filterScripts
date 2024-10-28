


#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube import topeventcleaning


from icecube.offline_filterscripts.filter_segments.inclined_filter import (
    inclined_filter,)  # For InclinedFilter

import numpy as np
import glob

from icecube.phys_services.which_split import which_split

KEEP_FILTERED_ONLY = True


#usage python selectInclinedEvents.py --vertEvents


fileDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck_afterSplit/"
outputDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck_afterFilter/"

GCD = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"


fileList = sorted(glob.glob(fileDir+"Fe*.i3.*"))[:1]
name = ""

def keep_filtered(frame):
  pass_filter = False
  if frame.Stop == icetray.I3Frame.Physics:
    if abs(frame["IceTopInclined_24"].value - 1) <= 0.00001:
      # print(frame["IceTopInclined_24"].value-1)
      pass_filter = True
  return pass_filter


for ifile in fileList:
  fileName = ifile.split("/")[-1]
  tray = I3Tray()
  tray.AddModule("I3Reader","reader",
               # filenameList = args.input,
               filenameList = [GCD]+[ifile],
              )
  InclinedSplit = "IceTopSplitIncl"
  tray.Add(inclined_filter, name+"_inclined_boolean_only",
                 InputTriggerHierarchy="TriggerHierarchy",
                 Input="CleanedTankPulses",
                 zenith_threshold=45*I3Units.degree,
                 Output="IceTopInclined_24",
                 If=which_split(split_name=InclinedSplit))
  if KEEP_FILTERED_ONLY == True:
    tray.AddModule(keep_filtered,
      streams=[icetray.I3Frame.Physics])
    tray.Add("I3OrphanQDropper")

  tray.AddModule("I3Writer","i3writer",
                filename=str(outputDir)+fileName,
                streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
                # DropOrphanStreams=[icetray.I3Frame.DAQ]
                )

  tray.Execute()
  tray.Finish()