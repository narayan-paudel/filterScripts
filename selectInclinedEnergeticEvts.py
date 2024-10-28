#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio

import numpy as np
import glob


#usage python selectInclinedEvents.py --vertEvents

inclinationCut = 45 #degree
energyCut = 10**16 #eV
# fileDir = "/data/sim/IceTop/2023/generated/untriggered/dataSetClean/"
fileDir = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/"

outputDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck/"


def testZenith(frame,zcut):
  # print(np.rad2deg(frame["MCPrimary"].dir.zenith))
  return np.rad2deg(frame["MCPrimary"].dir.zenith) >= zcut

def testEnergy(frame,Ecut):
  # print(frame["MCPrimary"].energy)
  return frame["MCPrimary"].energy * 10**9 >= Ecut

def testTrigger(frame,trigFlag):
  # print(frame[trigFlag].value)
  return abs(frame[trigFlag].value -1.0) < 0.001

def test_trigger_configID(frame,trigger_hierarchy_name,configID):
  if not frame.Has(trigger_hierarchy_name):
    icetray.logging.log_error(f"No {trigger_hierarchy_name} in frame")
  triggerHierarchy = frame[trigger_hierarchy_name]
  triggers = [t for t in triggerHierarchy if (t.key.config_id == configID and t.fired)]
  return len(triggers)>0

def selectEvents(frame,zcut,Ecut):
  # return testZenith(frame,zcut) and testEnergy(frame,Ecut) and testTrigger(frame,"tank7_3000")
  return testZenith(frame,zcut) and testEnergy(frame,Ecut) and test_trigger_configID(frame,"QTriggerHierarchy",30043)

def selectVertEvents(frame,zcut,Emin,Emax):
  return np.rad2deg(frame["MCPrimary"].dir.zenith) <= zcut and frame["MCPrimary"].energy * 10**9 >= Emin and frame["MCPrimary"].energy * 10**9 < Emax

fileList = sorted(glob.glob(fileDir+"Fe*.i3.*"))[:1]
for ifile in fileList:
  fileName = ifile.split("/")[-1]
  tray = I3Tray()
  tray.AddModule("I3Reader","reader",
               # filenameList = args.input,
               filename = ifile,
              )
  tray.AddModule(selectEvents,"selEvts",
                Ecut = 10**16,
                zcut = 45,
                streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics])

  def Unify(frame, Keys, Output):
    """
    Simple utility to merge RecoPulseSerieses into a single Union.
    """
    extants = [k for k in Keys if k in frame]
    union = dataclasses.I3RecoPulseSeriesMapUnion(frame, extants)
    frame[Output] = union

  tray.Add(Unify,"UnionHLCSLC",
    Keys=["OfflineIceTopHLCTankPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge"],
    Output='IceTopTankPulses',
    streams=[icetray.I3Frame.DAQ]
    )

  relevant_keeps = ["tank7_3000","SMT102","SMT273","QTriggerHierarchy","QFilterMask",
  "OfflineIceTopSLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopSLCTankPulsesCleanTimeCleanCharge",
  "OfflineIceTopHLCVEMPulsesCleanTimeCleanCharge","OfflineIceTopHLCTankPulsesCleanTimeCleanCharge",
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