#!/usr/bin/env python3

'''
This example shows how to include IceTop DOMSets to the GCD file in addition to default InIce DOMSets
'''
from icecube.icetray import I3Tray, I3Units
from icecube import dataclasses,icetray,dataio
from icecube.trigger_sim import GetDefaultDOMSets
from icecube.trigger_sim.InjectDefaultDOMSets import InjectDefaultDOMSets
# from icecube.filterscripts.offlineL2 import LaputopStandard
# from icecube.toprec import LaputopStandard
# from icecube.toprec import LaputopSmallShower
from icecube.icetop_Level3_scripts.segments.level3_IceTop import level3_IceTop
from icecube.recclasses import I3LaputopParams, LaputopParameter as Par

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import numpy as np


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input',"-i",nargs="+",type=str,default="",help="input simulation GCD for IceTop")
parser.add_argument('--output',"-o",type=str,default="../testInput.i3.gz",help='output simulation GCD for IceTop')
args = parser.parse_args()

# GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"
GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HGDomsets.i3.gz"

plotFolder = "/home/enpaudel/icecube/filterReco/plots/"
output = args.output.split("Gen")[0]+"IceTop7HG.i3.gz"
print("output",output)

keep_list = ["CleanIceTopRawData","CleanedHLCTankPulses","FilterMask","H4aWeight","tank7_3000","QTriggerHierarchy","QFilterMask"
"HLC6_5000","I3EventHeader","I3TriggerHierarchy","MCPrimary","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge",
"OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopHLCTankPulsesCleanTimeCleanCharge","InIcePulses","IceTopPulses","InIceDSTPulses","IceTopDSTPulses",
"I3SuperDST",
]

@icetray.traysegment
def clean_rename(tray,name=""):
  tray.AddModule("Keep","keep unusual time",
               keys = keep_list,
               # If = haveUnusualTime,
               )
  tray.AddModule("Rename","renameFiles",
               keys = ["OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopSLCVEMPulses",
               "OfflineIceTopSLCTankPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulses",
               "OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopHLCVEMPulses",
               "OfflineIceTopHLCTankPulsesCleanTimeCleanCharge","OfflineIceTopHLCTankPulses"],
               # If = haveUnusualTime,
               )




tray = I3Tray()
tray.AddModule("I3Reader","reader",
             filenameList=args.input,
            )
tray.Add(clean_rename)

tray.AddModule("I3Writer","i3writer",
            filename=output,
            streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            )

tray.Execute()
tray.Finish()

