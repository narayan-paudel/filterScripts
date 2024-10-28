#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube import topeventcleaning
from icecube.phys_services.which_split import which_split

icetray.set_log_level(icetray.logging.I3LogLevel.LOG_INFO)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob



plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = "/home/enpaudel/icecube/filterReco/plots/"

from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
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
fileDir = "/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/"

outputDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck_afterSplit/"

GCD = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"


# fileList = sorted(glob.glob(fileDir+"FeDAT000*.i3.*"))[:1]
fileList = sorted(glob.glob(fileDir+"FeDAT000*.i3.*"))
# fileList = ["/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz"]

def ifHighGain(domKey):
  return (domKey.om == 61+(domKey.string in [26,39,74])) or (domKey.om == 63+(domKey.string == 67))

# def domLocation()

class durEvent:
  def __init__(self,eventID,time_offset,zenith,azimuth,energy):
    self.eventID = eventID
    self.time_offset = time_offset
    self.zenith = zenith
    self.azimuth = azimuth
    self.energy = energy

class PulseLocTime:
  def __init__(self,om,time,pos):
    self.om = om
    self.time = time
    self.pos = pos


class hitDuration(icetray.I3ConditionalModule):
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
        icetray.logging.log_error("Error: hitDuration: Missing IceTopTankPulses")
      times_list = []
      hit_stations = []
      hit_omkeys = []
      for om,pulses in psm:
        if ifHighGain(om):
          for pulse in pulses:
            times_list.append(pulse.time)
      if len(times_list) > 3:
        dur_diff = np.diff(times_list)
        # dur = max(times_list) - min(times_list)
        self.event_list.append(durEvent(event_id,max(dur_diff),dur_diff,zenith,azimuth,energy))

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
      plt.savefig(plotFolder+"/deltaHitHist.pdf",transparent=False,bbox_inches='tight')
      plt.savefig(plotFolder+"/deltaHitHist.png",transparent=False,bbox_inches='tight')


class hitOffset(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox("OutBox")
        self.AddParameter("input_pulse","pulse","IceTopTankPulses")

    def Configure(self):
        self.event_list = []
        self.input_pulse = self.GetParameter("input_pulse")

    def get_position(self,omgeo,omkey):
      return omgeo[omkey].position

    def get_distance(self,pos1,pos2):
      return np.sqrt((pos1.x-pos2.x)**2+(pos1.y-pos2.y)**2+(pos1.z-pos2.z)**2)

    def contained(self,frame,radius):
      return frame["MCPrimary"].pos.y**2 + frame["MCPrimary"].pos.x**2 <= radius**2

    def energySelect(self,frame,energy):
      return frame["MCPrimary"].energy > energy

    def DAQ(self, frame):
      event_id = frame["I3EventHeader"].event_id
      zenith = frame["MCPrimary"].dir.zenith
      azimuth = frame["MCPrimary"].dir.azimuth
      energy = frame["MCPrimary"].energy
      triggerHierarchy = frame["QTriggerHierarchy"]
      HG7Triggers = [t for t in triggerHierarchy if (t.key.config_id == 30043 and t.fired)]
      if len(HG7Triggers)>0 and self.contained(frame,410) and self.energySelect(frame,10**7):
        if not frame.Has("I3Geometry"):
          icetray.logging.log_info("No geometry")      
          # icetray.logging.log_debug("No geometry")
          # icetray.logging.log_warn("No geometry")
        omgeo = frame["I3Geometry"].omgeo

        # print("omgeo",omgeo[icetray.OMKey(81,61,0)].position.x)
        geom = frame["I3Geometry"]
        if frame.Has(self.input_pulse):
          psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,self.input_pulse)
          # psm = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,"OfflineIceTopHLCTankPulsesCleanTimeCleanCharge")
        else:
          icetray.logging.log_error(f"Error: hitDuration: Missing {self.input_pulse}")
        if "HLC" in self.input_pulse:
          self.name_suffix = "HLC"
        elif "SLC" in self.input_pulse:
          self.name_suffix = "SLC"
        else:
          self.name_suffix = self.input_pulse
        times_list = []
        hit_stations = []
        hit_omkeys = []
        timeLoc = []
        time_offset_station = []
        time_offset_tank = []
        for om,pulses in psm:
          if ifHighGain(om):
            for pulse in pulses:
              timeLoc.append(PulseLocTime(om,pulse.time,self.get_position(omgeo,om)))
        if len(timeLoc) > 3:
          timeLoc.sort(key=lambda x:x.time,reverse=False)
          for i in range(len(timeLoc)-1):
            # print("time",timeLoc[i+1].time,timeLoc[i+1].pos,timeLoc[i+1].time - timeLoc[i].time,self.get_distance(timeLoc[i+1].pos,timeLoc[i].pos),self.get_distance(timeLoc[i+1].pos,timeLoc[i].pos)/dataclasses.I3Constants.c)
            # print("time",timeLoc[i].om,timeLoc[i].pos,timeLoc[i+1].om,timeLoc[i+1].pos,self.get_distance(timeLoc[i+1].pos,timeLoc[i].pos),self.get_distance(timeLoc[i+1].pos,timeLoc[i].pos)/dataclasses.I3Constants.c)
            # timeOff = timeLoc[i+1].time - (timeLoc[i].time + self.get_distance(timeLoc[i+1].pos,timeLoc[i].pos)/dataclasses.I3Constants.c)
            timeOff = abs(timeLoc[i+1].time - timeLoc[i].time) - (self.get_distance(timeLoc[i+1].pos,timeLoc[i].pos)/dataclasses.I3Constants.c)
            # timeOff = abs(timeLoc[i+1].time - timeLoc[i].time)
            # timeOff = timeLoc[i+1].time - (timeLoc[i].time)
            # if timeOff > 7000:
            #   print(event_id,timeLoc[i].om,timeLoc[i].time,timeLoc[i+1].om,timeLoc[i+1].time,timeLoc[i+1].time - (timeLoc[i].time + self.get_distance(timeLoc[i+1].pos,timeLoc[i].pos)/dataclasses.I3Constants.c))
            # print("timeoff",timeOff)
            if timeLoc[i+1].om.string == timeLoc[i].om.string:
              time_offset_tank.append(timeOff)
            else:
              time_offset_station.append(timeOff)

          # dur = max(times_list) - min(times_list)
          # print(time_offset_tank)
          # print(f"time offset {max(time_offset):.1f} {min(time_offset):.1f}")
          self.event_list.append(durEvent(event_id,[time_offset_tank,time_offset_station],np.rad2deg(zenith),np.rad2deg(azimuth),energy))

      # print("pulse times",times_list)
      self.PushFrame(frame)

    def plot_offset_zen(self):
      whatToPlot = "minmax" #Mean,Median 
      # whatToPlot = "mean" #mean,median 
      fig = plt.figure(figsize=(8,5))
      gs = gridspec.GridSpec(nrows=1,ncols=1)
      ax = fig.add_subplot(gs[0])
      # durList = [ievt.duration for ievt in self.event_list]
      # bins = np.linspace(0,int(max(durList)+1))
      # bins = np.arange(0,max(durList)+1)
      # ax.hist(durList,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom[0],label='',alpha=0.8)
      if whatToPlot == "minmax":
        ax.plot([ielt.zenith for ielt in self.event_list],[min(ielt.time_offset[1]) for ielt in self.event_list],"o",color=colorsCustom[0],label='min inter',alpha=0.8)
        ax.plot([ielt.zenith for ielt in self.event_list],[max(ielt.time_offset[1]) for ielt in self.event_list],"o",color=colorsCustom[1],label='max inter',alpha=0.8)
        ax.plot([ielt.zenith for ielt in self.event_list if len(ielt.time_offset[0])>0],[min(ielt.time_offset[0]) for ielt in self.event_list if len(ielt.time_offset[0])>0],"o",color=colorsCustom[2],label='min intra',alpha=0.8)
        ax.plot([ielt.zenith for ielt in self.event_list if len(ielt.time_offset[0])>0],[max(ielt.time_offset[0]) for ielt in self.event_list if len(ielt.time_offset[0])>0],"o",color=colorsCustom[3],label='max intra',alpha=0.8)
      elif whatToPlot == "mean":
        ax.plot([ielt.zenith for ielt in self.event_list],[np.median(ielt.time_offset[0]) for ielt in self.event_list],"o",color=colorsCustom[0],label='tank offset',alpha=0.8)
        ax.plot([ielt.zenith for ielt in self.event_list],[np.median(ielt.time_offset[1]) for ielt in self.event_list],"o",color=colorsCustom[1],label='station offset',alpha=0.8)
      elif whatToPlot == "median":
        ax.plot([ielt.zenith for ielt in self.event_list],[np.mean(ielt.time_offset[0]) for ielt in self.event_list],"o",color=colorsCustom[0],label='tank offset',alpha=0.8)
        ax.plot([ielt.zenith for ielt in self.event_list],[np.mean(ielt.time_offset[1]) for ielt in self.event_list],"o",color=colorsCustom[1],label='station offset',alpha=0.8)

      ax.tick_params(axis='both',which='both', direction='in', labelsize=20, pad=8)
      ax.set_ylabel("time offset [ns]", fontsize=20)
      ax.set_xlabel(r"$\theta^{\circ}$", fontsize=20)
      # ax.set_xticks([1,3],["2023 December","2024 January"])
      ax.legend(loc="lower right",fontsize=20,ncol=1)
      # ax.set_ylim(0,5000)
      ax.set_yscale("log")
      # ax.set_xlim(0,4)
      ax.grid(True,alpha=0.5)
      plt.savefig(plotFolder+f"/offsetHitTime{whatToPlot}.pdf",transparent=False,bbox_inches='tight')
      plt.savefig(plotFolder+f"/offsetHitTime{whatToPlot}.png",transparent=False,bbox_inches='tight')
    def plot_offset_hist_station(self):
      whatToPlot = "all_zenith"
      whatToPlot = "zenith_binned"
      fig = plt.figure(figsize=(8,5))
      gs = gridspec.GridSpec(nrows=1,ncols=1)
      ax = fig.add_subplot(gs[0])
      station_offset = []
      tank_offset = []
      for ievt in self.event_list:
        station_offset += ievt.time_offset[1]
        tank_offset += ievt.time_offset[0]

      # durList = [ievt.duration for ievt in self.event_list]
      # bins = np.linspace(0,int(max(durList)+1))
      # bins = np.linspace(0,int(max(station_offset)+1))
      print("offsetCheck",min(tank_offset),max(tank_offset),min(station_offset),max(station_offset))
      bins = np.linspace(int(min(station_offset)-1),int(max(station_offset)+1),(int(max(station_offset)+1)-int(min(station_offset)+1))+1)
      # bins = np.linspace(-3500,1000,4501)
      # bins = np.linspace(-500,500,1001)
      # bins = np.linspace(int(min(station_offset)-1),int(max(station_offset)+1),1200)
      # bins = np.linspace(0,1000,1001)
      if whatToPlot == "all_zenith":
        ax.hist(station_offset,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom[0],label='station',alpha=0.8)
        ax.hist(tank_offset,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom[1],label='tank',alpha=0.8)
      else:
        for nbin, binStart in enumerate(sin2ZenBins[:-1]):
          lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
          highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
          evtZenBin = [ievt for ievt in self.event_list if lowEdge <= np.deg2rad(ievt.zenith) < highEdge]
          istation_offset = []
          itank_offset = []
          for ievt in evtZenBin:
            istation_offset += ievt.time_offset[1]
            itank_offset += ievt.time_offset[0]

          ax.hist(istation_offset,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom2[nbin],label=\
            r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi)\
            ,alpha=0.4)
          # ax.hist(itank_offset,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom2[nbin],label=\
          #   r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi)\
          #   ,alpha=0.4)
          # ax.hist(itank_offset,bins=bins,histtype="step",ls="--",lw = 2.5,color=colorsCustom2[nbin],label='',alpha=0.4)

      ax.tick_params(axis='both',which='both', direction='in', labelsize=20, pad=8)
      ax.set_ylabel("count", fontsize=20)
      ax.set_xlabel(r"offset [ns]", fontsize=20)
      # ax.set_xticks([1,3],["2023 December","2024 January"])
      # ax.legend(loc="lower right",fontsize=20,ncol=1)
      ax.legend(fontsize=12,ncol=1)
      ax.set_ylim(10**(-0.1),10**3)
      ax.set_yscale("log")
      ax.set_xlim(-300,300)
      ax.grid(True,alpha=0.5)
      plt.savefig(plotFolder+f"/offsetHitTimeHistStation{whatToPlot}_{self.name_suffix}.pdf",transparent=False,bbox_inches='tight')
      plt.savefig(plotFolder+f"/offsetHitTimeHistStation{whatToPlot}_{self.name_suffix}.png",transparent=False,bbox_inches='tight')

    def plot_offset_hist_tank(self):
      whatToPlot = "all_zenith"
      whatToPlot = "zenith_binned"
      fig = plt.figure(figsize=(8,5))
      gs = gridspec.GridSpec(nrows=1,ncols=1)
      ax = fig.add_subplot(gs[0])
      station_offset = []
      tank_offset = []
      for ievt in self.event_list:
        station_offset += ievt.time_offset[1]
        tank_offset += ievt.time_offset[0]

      # durList = [ievt.duration for ievt in self.event_list]
      # bins = np.linspace(0,int(max(durList)+1))
      # bins = np.linspace(0,int(max(station_offset)+1))
      print("offsetCheck",min(tank_offset),max(tank_offset),min(station_offset),max(station_offset))
      bins = np.linspace(int(min(station_offset)-1),int(max(station_offset)+1),(int(max(station_offset)+1)-int(min(station_offset)+1))+1)
      # bins = np.linspace(-3500,1000,4501)
      # bins = np.linspace(-500,500,1001)
      # bins = np.linspace(int(min(station_offset)-1),int(max(station_offset)+1),1200)
      # bins = np.linspace(0,1000,1001)
      if whatToPlot == "all_zenith":
        ax.hist(station_offset,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom[0],label='station',alpha=0.8)
        ax.hist(tank_offset,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom[1],label='tank',alpha=0.8)
      else:
        for nbin, binStart in enumerate(sin2ZenBins[:-1]):
          lowEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin]))
          highEdge = np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))
          evtZenBin = [ievt for ievt in self.event_list if lowEdge <= np.deg2rad(ievt.zenith) < highEdge]
          istation_offset = []
          itank_offset = []
          for ievt in evtZenBin:
            istation_offset += ievt.time_offset[1]
            itank_offset += ievt.time_offset[0]

          # ax.hist(istation_offset,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom2[nbin],label=\
          #   r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi)\
          #   ,alpha=0.4)
          ax.hist(itank_offset,bins=bins,histtype="step",ls="-",lw = 2.5,color=colorsCustom2[nbin],label=\
            r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(np.arcsin(np.sqrt(sin2ZenBins[nbin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[nbin+1]))*180.0/np.pi)\
            ,alpha=0.4)
          # ax.hist(itank_offset,bins=bins,histtype="step",ls="--",lw = 2.5,color=colorsCustom2[nbin],label='',alpha=0.4)

      ax.tick_params(axis='both',which='both', direction='in', labelsize=20, pad=8)
      ax.set_ylabel("count", fontsize=20)
      ax.set_xlabel(r"offset [ns]", fontsize=20)
      # ax.set_xticks([1,3],["2023 December","2024 January"])
      # ax.legend(loc="lower right",fontsize=20,ncol=1)
      ax.legend(fontsize=12,ncol=1)
      ax.set_ylim(10**(-0.1),10**3)
      ax.set_yscale("log")
      ax.set_xlim(-300,300)
      ax.grid(True,alpha=0.5)
      plt.savefig(plotFolder+f"/offsetHitTimeHistTank{whatToPlot}_{self.name_suffix}.pdf",transparent=False,bbox_inches='tight')
      plt.savefig(plotFolder+f"/offsetHitTimeHistTank{whatToPlot}_{self.name_suffix}.png",transparent=False,bbox_inches='tight')    





    def Finish(self):
      print("len",len(self.event_list))
      # self.plot_offset_zen()
      self.plot_offset_hist_station()
      self.plot_offset_hist_tank()










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

# tray.AddModule(hitDuration, name + "_hitDur",
#                # If=which_split(split_name="IceTopSplitIncl")
#                )

tray.AddModule(hitOffset, name + "_hitDur",
                input_pulse = "IceTopTankPulses"
               # If=which_split(split_name="IceTopSplitIncl")
               )
# tray.AddModule(hitOffset, name + "_hitDurHLC",
#                 input_pulse = "OfflineIceTopHLCTankPulsesCleanTimeCleanCharge"
#                # If=which_split(split_name="IceTopSplitIncl")
#                )
# tray.AddModule(hitOffset, name + "_hitDurSLC",
#                 input_pulse = "OfflineIceTopSLCTankPulsesCleanTimeCleanCharge"
#                # If=which_split(split_name="IceTopSplitIncl")
#                )


# tray.AddModule("I3Writer","i3writer",
#               filename=str(outputDir)+fileName,
#               streams=[icetray.I3Frame.Simulation,icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#               )

tray.Execute()
tray.Finish()