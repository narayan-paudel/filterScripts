#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio

import numpy as np
import glob


#usage python selectInclinedEvents.py --vertEvents

inclinationCut = 45 #degree
energyCut = 10**16 #eV
fileDir = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/"
outputDir = "/data/sim/IceTop/2023/generated/untriggered/filterStudy/dataSet/"

fileList = sorted(glob.glob(fileDir+"*.i3.*"))[:]
for ifile in fileList:
  fileName = ifile.split("/")[-1].replace("GenDetFiltProcUniqueCleanVEMEvts","IceTop7HG")
  tray = I3Tray()
  tray.AddModule("I3Reader","reader",
               # filenameList = args.input,
               filename = ifile,
              )
  tray.AddModule("Delete", "pulse removal before rename",
                 keys = ["OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses",
                 "OfflineIceTopHLCVEMPulses","OfflineIceTopSLCVEMPulses"]
                 )


  tray.AddModule("Rename","renameFiles",
               keys = ["OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopSLCVEMPulses",
               "OfflineIceTopSLCTankPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulses",
               "OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopHLCVEMPulses",
               "OfflineIceTopHLCTankPulsesCleanTimeCleanCharge","OfflineIceTopHLCTankPulses"],
               # If = haveUnusualTime,
               )
  def Unify(frame, Keys, Output):
    """
    Simple utility to merge RecoPulseSerieses into a single Union.
    """
    extants = [k for k in Keys if k in frame]
    union = dataclasses.I3RecoPulseSeriesMapUnion(frame, extants)
    frame[Output] = union

  tray.Add(Unify,"UnionHLCSLC",
    Keys=["OfflineIceTopHLCTankPulses","OfflineIceTopSLCTankPulses"],
    Output='IceTopTankPulses',
    streams=[icetray.I3Frame.DAQ]
    )



  relevant_keeps = ["tank7_3000","SMT102","SMT273","QTriggerHierarchy","QFilterMask",
  "OfflineIceTopSLCVEMPulses","OfflineIceTopSLCTankPulses",
  "OfflineIceTopHLCVEMPulses","OfflineIceTopHLCTankPulses",
  "MCPrimary","IceTopPulses","I3Triggers","I3EventHeader","H4aWeight","HLC6_5000","IceTopTankPulses",
  "IceTopDSTPulses"]

  tray.AddModule("Keep", "keep_relevant",
                 keys = relevant_keeps,
                 If=lambda x:True
                 )


  tray.AddModule("I3Writer","i3writer",
                filename=str(outputDir)+fileName,
                streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ],
                )

  tray.Execute()
  tray.Finish()