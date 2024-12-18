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
from cleani3 import clean_rename

import glob


file_path = "/data/user/enpaudel/triggerStudy/simFiles/dataSetClean1_6/"
input_files = sorted(glob.glob(file_path+"*.i3.*"))
# GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"
GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HGDomsets.i3.gz"

# primary = args.input[0].split("/")[-1][0]
# example python reconstructIceTopS125.py -i 
# /home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique1_6/HeDAT004171GenDetFiltProcUnique.i3.bz2
# -o /home/enpaudel/icecube/triggerStudy/simFiles/dataSetReco/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz

plotFolder = "/home/enpaudel/icecube/filterReco/plots/"

keep_list = ["CleanIceTopRawData","CleanedHLCTankPulses","FilterMask","H4aWeight","tank7_3000","QTriggerHierarchy","QFilterMask"
"HLC6_5000","I3EventHeader","I3TriggerHierarchy","MCPrimary","OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge",
"OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopHLCTankPulsesCleanTimeCleanCharge","InIcePulses","IceTopPulses","InIceDSTPulses","IceTopDSTPulses",
"I3SuperDST",
]

for ifile in file_list:
  ioutput = ifile.split("/")[-1][0]  
  # tray = I3Tray()
  # tray.AddModule("I3Reader","reader",
  #              filenameList=ifile,
  #             )
  # tray.Add(clean_rename)

  # tray.AddModule("I3Writer","i3writer",
  #             filename=args.output,
  #             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
  #             )

  # tray.Execute()
  # tray.Finish()