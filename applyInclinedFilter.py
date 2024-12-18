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


import argparse
parser = argparse.ArgumentParser(
    prog="applySplittingInclined",
    description="apply inclined splitting")
parser.add_argument("-i", "--input", action="store", default=None, nargs="+",
                    dest="INPUT", help="Input i3 file(s) to process, separated by spaces", required=True)
parser.add_argument("-o", "--output", action="store", default=None,
                    dest="OUTPUT", help="Output i3 file", required=True)
parser.add_argument("-g", "--gcd", action="store", default="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz",
                    dest="GCD", help="GCD file for input i3 file", required=True)
args = parser.parse_args()



name = ""

def keep_filtered(frame):
  pass_filter = False
  if frame.Stop == icetray.I3Frame.Physics:
    if abs(frame["IceTopInclined_24"].value - 1) <= 0.00001:
      # print(frame["IceTopInclined_24"].value-1)
      pass_filter = True
  return pass_filter



tray = I3Tray()
tray.AddModule("I3Reader","reader",
             # filenameList = args.input,
             filenameList = [args.GCD]+args.INPUT,
            )
InclinedSplit = "IceTopSplitIncl"
fullEventSplit = "fullevent"
tray.Add(inclined_filter, name+"_inclined_boolean_only",
               InputTriggerHierarchy="QTriggerHierarchy",
               Input="CleanedTankPulses",
               zenith_threshold=45*I3Units.degree,
               Output="IceTopInclined_24",
               Split = InclinedSplit,
               # If=which_split(split_name=InclinedSplit)
               )
tray.Add(inclined_filter, name+"_full_event",
               InputTriggerHierarchy="QTriggerHierarchy",
               Input="IceTopTankPulses",
               zenith_threshold=45*I3Units.degree,
               Output="IceTopInclined_24",
               Split = fullEventSplit,
               # If=which_split(split_name=fullEventSplit)
               )
# if KEEP_FILTERED_ONLY == True:
#   tray.AddModule(keep_filtered,
#     streams=[icetray.I3Frame.Physics],
#     If=which_split(split_name=InclinedSplit))
#   tray.Add("I3OrphanQDropper")

tray.AddModule("I3Writer","i3writer",
              filename=args.OUTPUT,
              streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
              # DropOrphanStreams=[icetray.I3Frame.DAQ]
              )

tray.Execute()
tray.Finish()