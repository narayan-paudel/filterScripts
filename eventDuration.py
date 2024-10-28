#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube import topeventcleaning
from icecube.phys_services.which_split import which_split

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob



plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = "/home/enpaudel/icecube/filterReco/plots/"

from customColors import qualitative_colors

colorsCustom = qualitative_colors(5)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']



#usage python selectInclinedEvents.py --vertEvents

# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('--input',"-i", type=str, 
#   default="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305.i3.gz",
#  help="input simulation GCD")
# parser.add_argument('--output',"-o", type=str,
#  default="/data/user/enpaudel/triggerStudy/simFiles/modified_GCD/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HG.i3.gz",
#   help='output simulation GCD')
# args = parser.parse_args()

inclinationCut = 45 #degree
energyCut = 10**16 #eV
sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
energyBins = 10**np.linspace(6, 8.0, 21)
# fileDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck/"
fileDir = "/data/sim/IceTop/2023/generated/untriggered/dataSetClean/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz"

outputDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck_afterSplit/"

GCD = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"


# fileList = sorted(glob.glob(fileDir+"Fe*.i3.*"))[:1]
fileList = ["/data/sim/IceTop/2023/generated/untriggered/dataSetClean/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz"]

def ifHighGain(domKey):
  return (domKey.om == 61+(domKey.string in [26,39,74])) or (domKey.om == 63+(domKey.string == 67))

class durEvent:
  def __init__(self,eventID,duration,zenith,azimuth,energy):
    self.eventID = eventID
    self.duration = duration
    self.zenith = zenith
    self.azimuth = azimuth
    self.energy = energy

class EventDuration(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox("OutBox")

    def Configure(self):
        self.event_list = []

    def DAQ(self, frame):
      event_id = frame["I3EventHeader"].event_id
      zenith = frame["MCPrimary"].dir.zenith
      azimuth = frame["MCPrimary"].dir.azimuth
      energy = frame["MCPrimary"].energy
      if frame.Has("IceTopTankPulses"):
        psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,"IceTopTankPulses")
      else:
        icetray.logging.log_error("Error: EventDuration: Missing IceTopTankPulses")
      times_list = []
      hit_stations = []
      hit_omkeys = []
      for om,pulses in psm:
        if ifHighGain(om):
          for pulse in pulses:
            times_list.append(pulse.time)
      if len(times_list) > 3:
        dur = max(times_list) - min(times_list)
        self.event_list.append(durEvent(event_id,dur,zenith,azimuth,energy))

      # print("pulse times",times_list)
      self.PushFrame(frame)
    def Finish(self):
      print("len",len(self.event_list))

      fig = plt.figure(figsize=(8,5))
      gs = gridspec.GridSpec(nrows=1,ncols=1)
      ax = fig.add_subplot(gs[0])
      durList = [ievt.duration for ievt in self.event_list]
      # bins = np.linspace(0,int(max(durList)+1))
      bins = np.arange(0,max(durList)+1)
      ax.hist(durList,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom[0],label='',alpha=0.8)
      ax.tick_params(axis='both',which='both', direction='in', labelsize=20, pad=8)
      ax.set_ylabel("count", fontsize=20)
      ax.set_xlabel("time [ns]", fontsize=20)
      # ax.set_xticks([1,3],["2023 December","2024 January"])
      ax.legend(loc="lower right",fontsize=20,ncol=2)
      # ax.set_ylim(0,35)
      ax.set_yscale("log")
      # ax.set_xlim(0,4)
      ax.grid(True,alpha=0.5)
      plt.savefig(plotFolder+"/durationHist.pdf",transparent=False,bbox_inches='tight')
      plt.savefig(plotFolder+"/durationHist.png",transparent=False,bbox_inches='tight')




def testTrigger(frame,trigFlag):
  # print(frame[trigFlag].value)
  return abs(frame[trigFlag].value -1.0) < 0.001



name = ""
tray = I3Tray()
tray.AddModule("I3Reader","reader",
             # filenameList = args.input,
             # filenameList = [GCD]+[ifile],
             filenameList = [GCD]+fileList,
            )

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

tray.AddModule(EventDuration, name + "_IceTopTrigCopy",
               # If=which_split(split_name="IceTopSplitIncl")
               )


# tray.AddModule("I3Writer","i3writer",
#               filename=str(outputDir)+fileName,
#               streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#               )

tray.Execute()
tray.Finish()