#!/usr/bin/env python3

"""
This script produces p frames from q frame using null split and cluster cleaning.
"""

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube import topeventcleaning

import numpy as np
import glob

from icecube.phys_services.which_split import which_split



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
parser.add_argument("-f","--full_event", dest="FULL_EVENT",default=False, action="store_true", required=False,
                        help="create a p frame from Q frame")

args = parser.parse_args()





name = ""

tray = I3Tray()
tray.AddModule("I3Reader","reader",
             # filenameList = args.input,
             filenameList = [args.GCD]+args.INPUT,
            )
if args.FULL_EVENT == True:
  tray.AddModule("I3NullSplitter", "fullevent")

tray.AddModule("I3TopHLCClusterCleaning", name + "_Inclined",
               SubEventStreamName="IceTopSplitIncl",
               InputPulses="IceTopTankPulses",
               OutputPulses="CleanedTankPulses",
               BadTankList="TankPulseMergerExcludedTanksIncl",
               ExcludedTanks="ClusterCleaningExcludedTanksIncl",
               InterStationTimeTolerance=200.0 * I3Units.ns,  # Default
               IntraStationTimeTolerance=200.0 * I3Units.ns,  # Default
               # If=lambda frame: "IceTopTankPulses" in frame and len(frame["IceTopTankPulses"])>0,
               If=lambda frame: "IceTopTankPulses" in frame,
               )
class PhysicsCopyTriggers(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox("OutBox")

    def Configure(self):
        pass

    def Physics(self, frame):
        if frame.Has("QTriggerHierarchy"):
            myth = dataclasses.I3TriggerHierarchy.from_frame(frame,
                                                             "QTriggerHierarchy")
            if frame.Has("TriggerHierarchy"):
                icetray.logging.log_error("ERROR: PhysicsCopyTriggers: triggers in frame")
            else:
                frame.Put("TriggerHierarchy", myth)
        else:
            icetray.logging.log_error("Error: PhysicsCopyTriggers: Missing QTriggerHierarchy")
        self.PushFrame(frame)

# Null split and IceTop split need the triggerHierarchy put back in P frame

tray.AddModule(PhysicsCopyTriggers, name + "_IceTopTrigCopy",
               If=which_split(split_name="IceTopSplitIncl"))

tray.Add("I3OrphanQDropper")


tray.AddModule("I3Writer","i3writer",
              filename=args.OUTPUT,
              streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
              # DropOrphanStreams=[icetray.I3Frame.DAQ]
              )

tray.Execute()
tray.Finish()