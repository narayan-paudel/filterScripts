#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube import recclasses

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


import glob

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--vertEvents", dest="vertEvents",
                    default=False, action="store_true", required=False,
                    help="plot vertical events")
args = parser.parse_args()

from customColors import qualitative_colors
plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)


icetray.set_log_level(icetray.logging.I3LogLevel.LOG_INFO)

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = "/home/enpaudel/icecube/filterReco/plots/"

from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


inclinationCut = 45 #degree
energyCut = 10**16 #eV
sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
energyBins = 10**np.linspace(6, 8.0, 21)

fileDir = "/data/sim/IceTop/2023/generated/untriggered/filterStudy/filtered/inclinedSplit/"

# outputDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck_afterSplit/"

GCD = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"

plotSuffix = ""
# fileList = sorted(glob.glob(fileDir+"FeDAT000*.i3.*"))[:1]
fileList = sorted(glob.glob(fileDir+"Fe*DAT0*.i3.*"))[:5]
# fileList = ["/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz"]

def openingAngle(theta1,phi1,theta2,phi2):
  return np.arccos(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2))

def test7HG(frame):
  triggerHierarchy = frame["QTriggerHierarchy"]
  icetop7HG = [t for t in triggerHierarchy if (t.key.config_id == 30043 and t.fired)]
  return len(icetop7HG)>0

def test_filter(frame,filter_bool):
  return abs(frame[filter_bool]-1)< 0.001


def excludeITSMT(frame):
  triggerHierarchy = frame["QTriggerHierarchy"]
  itSMT = [t for t in triggerHierarchy if (t.key.config_id == 102 and t.fired)]
  return len(itSMT)<1

class RecoEvent(object):
  """docstring for RecoEvent"""
  def __init__(self,event_id,x_true,y_true,x_reco,y_reco,energy,zenith_true,azimuth_true,zenith_reco,azimuth_reco,opening_angle,chi2,
    Qtot_HLC,Ntanks_HLC,Qtot_SLC,Ntanks_SLC,Qtot_both,Ntanks_both):
    super(RecoEvent, self).__init__()
    self.event_id = event_id
    self.x_true = x_true
    self.y_true = y_true
    self.x_reco = x_reco
    self.y_reco = y_reco
    self.r_diff = np.sqrt((self.x_true-self.x_reco)**2+(self.y_true-self.y_reco)**2)/I3Units.m
    self.energy = energy
    self.zenith_true = zenith_true
    self.azimuth_true = azimuth_true
    self.zenith_reco = zenith_reco
    self.azimuth_reco = azimuth_reco
    self.opening_angle = opening_angle
    self.chi2 = chi2
    self.Qtot_HLC = Qtot_HLC
    self.Ntanks_HLC = Ntanks_HLC
    self.Qtot_SLC = Qtot_SLC
    self.Ntanks_SLC = Ntanks_SLC
    self.Qtot_both = Qtot_both
    self.Ntanks_both = Ntanks_both

def containedEvents(evtList,radius):
  """only selects the events with core with in given radius of IceTop center"""
  return [ievt for ievt in evtList if ((ievt.x_true**2 + ievt.y_true**2) <= radius**2)]
    

class zenithCheck(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
    self.AddParameter("split","which split to use","IceTopSplitIncl")
    self.AddParameter("containment","use showers core within IceTop",True)
  def Configure(self):
    self.split = self.GetParameter("split")
    self.containment = self.GetParameter("containment")
    self.openingAngleList = []
    self.r_diffList = []
    self.chi2_list = []
    self.RecoEventList = []

  def GetTotalCharge(self,pulse_series):
    return sum([sum([c.charge for c in pulse_series[i]]) for i in pulse_series.keys()])

  def GetNTank(self,pulse_series):
    Nch = 0
    for om,pulses in pulse_series:
      for pulse in pulses:
        Nch +=1
        break
    return Nch

  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
    # print(frame["ShowerPlaneParams"].chi2)
    # print(frame["I3EventHeader"].sub_event_stream,self.split,frame["I3EventHeader"].sub_event_stream == self.split)    
    if frame["I3EventHeader"].sub_event_stream == self.split:
      event_id = frame["I3EventHeader"].event_id
      energy = frame["MCPrimary"].energy
      xcore = frame["MCPrimary"].pos.x
      ycore = frame["MCPrimary"].pos.y
      xcore_reco = frame["ShowerCOG"].pos.x
      ycore_reco = frame["ShowerCOG"].pos.y
      r_diff = np.sqrt((xcore-xcore_reco)**2+(ycore-ycore_reco)**2)/I3Units.m
      # print("cores",(xcore,ycore),(xcore_reco,ycore_reco),(xcore-xcore_reco,ycore-ycore_reco),r_diff)
      zenith_true = frame["MCPrimary"].dir.zenith
      azimuth_true = frame["MCPrimary"].dir.azimuth
      zenith_reco = frame["ShowerPlane"].dir.zenith
      azimuth_reco = frame["ShowerPlane"].dir.azimuth
      chi2_reco = frame["ShowerPlaneParams"].chi2
      openAngle = openingAngle(zenith_true,azimuth_true,zenith_reco,azimuth_reco)
      HLC_pulse = frame["OfflineIceTopHLCVEMPulses"]
      Qtot_HLC = self.GetTotalCharge(HLC_pulse)
      Ntanks_HLC = self.GetNTank(HLC_pulse)
      SLC_pulse = frame["OfflineIceTopSLCVEMPulses"]
      Qtot_SLC  = self.GetTotalCharge(SLC_pulse)
      Ntanks_SLC = self.GetNTank(SLC_pulse)
      if frame.Has("CleanedTankPulses"):
        cleaned_tank_pulse = frame["CleanedTankPulses"]
      else:
        print("No cleaned Tank pulses")
      cleaned_tank_pulse = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, "CleanedTankPulses")
      Qtot_both = self.GetTotalCharge(cleaned_tank_pulse)
      Ntanks_both = self.GetNTank(cleaned_tank_pulse)

      # print("zenith True",zenith_reco,zenith_true,np.arcsin(np.sqrt(self.zenithBin[0])),np.arcsin(np.sqrt(self.zenithBin[1])))
      # if np.arcsin(np.sqrt(self.zenithBin[0])) <= zenith_true < np.arcsin(np.sqrt(self.zenithBin[1])) and not np.isnan(openAngle):
      # if not np.isnan(openAngle) and not np.isnan(r_diff) and chi2_reco < 600:
      if not np.isnan(openAngle) and not np.isnan(r_diff):
        self.openingAngleList.append(openAngle*180.0/np.pi)
        self.r_diffList.append(r_diff)
        self.chi2_list.append(chi2_reco)
        self.RecoEventList.append(RecoEvent(event_id,xcore,ycore,xcore_reco,ycore_reco,energy,zenith_true,azimuth_true,zenith_reco,azimuth_reco,openAngle,chi2_reco,
          Qtot_HLC,Ntanks_HLC,Qtot_SLC,Ntanks_SLC,Qtot_both,Ntanks_both))


  def plotChi2Zenith(self):
    """ chi2 plot binned in zenith"""
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    p68 = np.percentile(self.r_diffList,68)
    print(p68)
    bins = np.linspace(-1,2000,2002)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(openingAngleList,chi2_list,"o",c=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"opening angle $\psi^{\circ}$", fontsize=22)
    self.ax.set_ylabel(r"$\chi^{2}$", fontsize=22)
    self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    self.ax.legend()
    plt.savefig(plotFolder+"/chi2OpenAngle"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/chi2OpenAngle"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')

  def plotChi2Energy(self):
    """ chi2 plot binned in zenith"""
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    p68 = np.percentile(self.r_diffList,68)
    print(p68)
    bins = np.linspace(-1,2000,2002)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(openingAngleList,chi2_list,"o",c=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"opening angle $\psi^{\circ}$", fontsize=22)
    self.ax.set_ylabel(r"$\chi^{2}$", fontsize=22)
    self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    self.ax.legend()
    plt.savefig(plotFolder+"/chi2OpenAngle"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/chi2OpenAngle"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')


  def plotOpenAngleEnergy(self):
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,52,102)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      self.ax.hist(openingAngleList,bins=bins,histtype="step",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}$^{{\circ}}$".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"opening angle $\psi^{\circ}$", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_xlim(0,100)
    self.ax.set_ylim(0.9,3*10**2)
    self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/openAngleHG7Only"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/openAngleHG7Only"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')
    plt.close()


  def plotOpenAngleZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      self.ax.hist(openingAngleList,bins=bins,histtype="step",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}$^{{\circ}}$".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"opening angle $\psi^{\circ}$", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_xlim(0,100)
    self.ax.set_ylim(0.9,3*10**2)
    self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/openAngleHG7Only"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/openAngleHG7Only"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')
    plt.close()


  def plotCoreRecoX(self):
    """ chi2 plot binned in zenith"""
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    p68 = np.percentile(self.r_diffList,68)
    print(p68)
    bins = np.linspace(-1,2000,2002)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      x_true_list = [ievt.x_true for ievt in evtEBin]
      y_true_list = [ievt.y_true for ievt in evtEBin]
      x_reco_list = [ievt.x_reco for ievt in evtEBin]
      y_reco_list = [ievt.y_reco for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(x_true_list,x_reco_list,"o",c=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"x$_{\mathrm{true}}$", fontsize=22)
    self.ax.set_ylabel(r"x$_{\mathrm{reco}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreReco_x"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreReco_x"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')
  def plotCoreRecoY(self):
    """ chi2 plot binned in zenith"""
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    p68 = np.percentile(self.r_diffList,68)
    print(p68)
    bins = np.linspace(-1,2000,2002)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      x_true_list = [ievt.x_true for ievt in evtEBin]
      y_true_list = [ievt.y_true for ievt in evtEBin]
      x_reco_list = [ievt.x_reco for ievt in evtEBin]
      y_reco_list = [ievt.y_reco for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(y_true_list,y_reco_list,"o",c=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"y$_{\mathrm{true}}$", fontsize=22)
    self.ax.set_ylabel(r"y$_{\mathrm{reco}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreReco_y"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreReco_y"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')


  def plotCoreRecoXZenith(self):
    """ chi2 plot binned in zenith"""
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      x_true_list = [ievt.x_true for ievt in evtZBin]
      y_true_list = [ievt.y_true for ievt in evtZBin]
      x_reco_list = [ievt.x_reco for ievt in evtZBin]
      y_reco_list = [ievt.y_reco for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(x_true_list,x_reco_list,"o",c=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"x$_{\mathrm{true}}$", fontsize=22)
    self.ax.set_ylabel(r"x$_{\mathrm{reco}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreReco_x"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreReco_x"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')
  def plotCoreRecoYZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      x_true_list = [ievt.x_true for ievt in evtZBin]
      y_true_list = [ievt.y_true for ievt in evtZBin]
      x_reco_list = [ievt.x_reco for ievt in evtZBin]
      y_reco_list = [ievt.y_reco for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(y_true_list,y_reco_list,"o",c=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"y$_{\mathrm{true}}$", fontsize=22)
    self.ax.set_ylabel(r"y$_{\mathrm{reco}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreReco_y"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreReco_y"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')


  def plotCoreRecoRZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.hist(r_list,bins=bins,histtype="step",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"r$_{\mathrm{diff}}$", fontsize=22)
    self.ax.set_ylabel(r"count", fontsize=22)
    # self.ax.set_yscale("log")
    self.ax.set_xlim(0,600)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreReco_r"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreReco_r"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')

  def plotCoreRecoREnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.hist(r_list,bins=bins,histtype="step",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"r$_{\mathrm{diff}}$", fontsize=22)
    self.ax.set_ylabel(r"count", fontsize=22)
    # self.ax.set_yscale("log")
    self.ax.set_xlim(0,600)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreReco_r"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreReco_r"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')

  def plotPlaneRecoZenZenith(self):
    """ chi2 plot binned in zenith"""
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      zenith_true_list = [ievt.zenith_true*180.0/np.pi for ievt in evtZBin]
      azimuth_true_list = [ievt.azimuth_true*180.0/np.pi for ievt in evtZBin]
      zenith_reco_list = [ievt.zenith_reco*180.0/np.pi for ievt in evtZBin]
      azimuth_reco_list = [ievt.azimuth_reco*180.0/np.pi for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(zenith_true_list,zenith_reco_list,"o",c=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\theta$$_{\mathrm{true}}$", fontsize=22)
    self.ax.set_ylabel(r"$\theta$$_{\mathrm{reco}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/planeReco_zen"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/planeReco_zen"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')

  def plotPlaneRecoAziZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      zenith_true_list = [ievt.zenith_true*180.0/np.pi for ievt in evtZBin]
      azimuth_true_list = [ievt.azimuth_true*180.0/np.pi for ievt in evtZBin]
      zenith_reco_list = [ievt.zenith_reco*180.0/np.pi for ievt in evtZBin]
      azimuth_reco_list = [ievt.azimuth_reco*180.0/np.pi for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(azimuth_true_list,azimuth_reco_list,"o",c=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\phi$$_{\mathrm{true}}$", fontsize=22)
    self.ax.set_ylabel(r"$\phi$$_{\mathrm{reco}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/planeReco_azi"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/planeReco_azi"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')

  def plotPlaneRecoZenEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenith_true_list = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      azimuth_true_list = [ievt.azimuth_true*180.0/np.pi for ievt in evtEBin]
      zenith_reco_list = [ievt.zenith_reco*180.0/np.pi for ievt in evtEBin]
      azimuth_reco_list = [ievt.azimuth_reco*180.0/np.pi for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(zenith_true_list,zenith_reco_list,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\theta$$_{\mathrm{true}}$", fontsize=22)
    self.ax.set_ylabel(r"$\theta$$_{\mathrm{reco}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/planeReco_zen"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/planeReco_zen"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')
    
  def plotPlaneRecoAziEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenith_true_list = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      azimuth_true_list = [ievt.azimuth_true*180.0/np.pi for ievt in evtEBin]
      zenith_reco_list = [ievt.zenith_reco*180.0/np.pi for ievt in evtEBin]
      azimuth_reco_list = [ievt.azimuth_reco*180.0/np.pi for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(azimuth_true_list,azimuth_reco_list,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\phi$$_{\mathrm{true}}$", fontsize=22)
    self.ax.set_ylabel(r"$\phi$$_{\mathrm{reco}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/planeReco_azi"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/planeReco_azi"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')

  def plotCoreOpenAngleEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(openingAngleList,r_list,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_ylabel(r"r$_{\mathrm{diff}}$", fontsize=22)
    self.ax.set_xlabel(r"opening angle $\psi^{\circ}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreOpenAngle"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreOpenAngle"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')
  def plotCoreOpenAngleZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(openingAngleList,r_list,"o",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_ylabel(r"r$_{\mathrm{diff}}$", fontsize=22)
    self.ax.set_xlabel(r"opening angle $\psi^{\circ}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreOpenAngle"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreOpenAngle"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')

  def plotCoreChi2Energy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(chi2_list,r_list,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_ylabel(r"r$_{\mathrm{diff}}$", fontsize=22)
    self.ax.set_xlabel(r"$\chi^{2}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreChi2"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreChi2"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')

  def plotCoreChi2Zenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(chi2_list,r_list,"o",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_ylabel(r"r$_{\mathrm{diff}}$", fontsize=22)
    self.ax.set_xlabel(r"$\chi^{2}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/coreChi2"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreChi2"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')



  def plotZenCoreEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(zenithList,r_list,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_ylabel(r"r$_{\mathrm{diff}}$", fontsize=22)
    self.ax.set_xlabel(r"$\theta^{\circ}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/zenCore"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/zenCore"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')

  def plotZenOpenAngleEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(zenithList,openingAngleList,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_ylabel(r"opening angle $\psi^{\circ}$", fontsize=22)
    self.ax.set_xlabel(r"$\theta^{\circ}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/zenOpenAngle"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/zenOpenAngle"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')

  def plotZenChi2Energy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(zenithList,chi2_list,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_ylabel(r"$\chi^{2}$", fontsize=22)
    self.ax.set_xlabel(r"$\theta^{\circ}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/zenChi2"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/zenChi2"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')

  def plotnTankSLCEnergyZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      Ntanks_SLC_list = [ievt.Ntanks_SLC for ievt in evtZBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtZBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(Ntanks_SLC_list,energy_list,"o",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"N$_{\mathrm{SLC\,hits}}$", fontsize=22)
    self.ax.set_ylabel(r"log$_{10}$(E/eV)", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/nTanksSLCEnergy"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/nTanksSLCEnergy"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')
 
  def plotQtotSLCEnergyZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      Qtot_SLC_list = [np.log10(ievt.Qtot_SLC) for ievt in evtZBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtZBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(Qtot_SLC_list,energy_list,"o",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"log10(Q$_{\mathrm{tot\,SLC}})$", fontsize=22)
    self.ax.set_ylabel(r"log$_{10}$(E/eV)", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/QtotSLCEnergy"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/QtotSLCEnergy"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')



  def plotnTankHLCEnergyZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      Ntanks_HLC_list = [ievt.Ntanks_HLC for ievt in evtZBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtZBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(Ntanks_HLC_list,energy_list,"o",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"N$_{\mathrm{HLC\,hits}}$", fontsize=22)
    self.ax.set_ylabel(r"log$_{10}$(E/eV)", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/nTanksHLCEnergy"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/nTanksHLCEnergy"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')
 
  def plotQtotHLCEnergyZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      Qtot_HLC_list = [np.log10(ievt.Qtot_HLC) for ievt in evtZBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtZBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(Qtot_HLC_list,energy_list,"o",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"log10(Q$_{\mathrm{tot\,HLC}})$", fontsize=22)
    self.ax.set_ylabel(r"log$_{10}$(E/eV)", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/QtotHLCEnergy"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/QtotHLCEnergy"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')




  def plotnTankBothEnergyZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      Ntanks_both_list = [ievt.Ntanks_both for ievt in evtZBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtZBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(Ntanks_both_list,energy_list,"o",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"N$_{\mathrm{both\,hits}}$", fontsize=22)
    self.ax.set_ylabel(r"log$_{10}$(E/eV)", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/nTanksBothEnergy"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/nTanksBothEnergy"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')
 
  def plotQtotBothEnergyZenith(self):
    sin2ZenBins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.822]
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(-1,100,52)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(sin2ZenBins[:-1]):
      ncolor = colorsCustom2[ebin]
      lowEdge_Z = sin2ZenBins[ebin]
      highEdge_Z = sin2ZenBins[ebin+1]
      evtZBin = [ievt for ievt in self.RecoEventList if lowEdge_Z <= np.sin(ievt.zenith_true)**2 < highEdge_Z]
      Qtot_both_list = [np.log10(ievt.Qtot_both) for ievt in evtZBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtZBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtZBin]
      r_list = [ievt.r_diff for ievt in evtZBin]
      chi2_list = [ievt.chi2 for ievt in evtZBin]
      self.ax.plot(Qtot_both_list,energy_list,"o",color=colorsCustom2[ebin],label=r"{0:.1f}$^{{\circ}}$-{1:.1f}$^{{\circ}}$".format(
          np.arcsin(np.sqrt(sin2ZenBins[ebin]))*180.0/np.pi,np.arcsin(np.sqrt(sin2ZenBins[ebin+1]))*180.0/np.pi),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"log10(Q$_{\mathrm{tot\,both}})$", fontsize=22)
    self.ax.set_ylabel(r"log$_{10}$(E/eV)", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/QtotbothEnergy"+plotSuffix+str(self.split)+"Zenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/QtotbothEnergy"+plotSuffix+str(self.split)+"Zenith.pdf",transparent=False,bbox_inches='tight')

  ######################################################################zenith bins######################################
  def plotnTankSLCZenithEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      Ntanks_SLC_list = [ievt.Ntanks_SLC for ievt in evtEBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtEBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(Ntanks_SLC_list,zenithList,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"N$_{\mathrm{SLC\,hits}}$", fontsize=22)
    self.ax.set_ylabel(r"$\theta$$_{\mathrm{true}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/nTanksSLCZenith"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/nTanksSLCZenith"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')
 
  def plotQtotSLCZenithEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      Qtot_SLC_list = [np.log10(ievt.Qtot_SLC) for ievt in evtEBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(Qtot_SLC_list,zenithList,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"log10(Q$_{\mathrm{tot\,SLC}})$", fontsize=22)
    self.ax.set_ylabel(r"$\theta$$_{\mathrm{true}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/QtotSLCZenith"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/QtotSLCZenith"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')



  def plotnTankHLCZenithEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      Ntanks_HLC_list = [ievt.Ntanks_HLC for ievt in evtEBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtEBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(Ntanks_HLC_list,zenithList,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"N$_{\mathrm{HLC\,hits}}$", fontsize=22)
    self.ax.set_ylabel(r"$\theta$$_{\mathrm{true}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/nTanksHLCZenith"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/nTanksHLCZenith"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')
 
  def plotQtotHLCZenithEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      Qtot_HLC_list = [np.log10(ievt.Qtot_HLC) for ievt in evtEBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtEBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(Qtot_HLC_list,zenithList,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"log10(Q$_{\mathrm{tot\,HLC}})$", fontsize=22)
    self.ax.set_ylabel(r"$\theta$$_{\mathrm{true}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/QtotHLCZenith"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/QtotHLCZenith"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')




  def plotnTankBothZenithEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      Ntanks_both_list = [ievt.Ntanks_both for ievt in evtEBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtEBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(Ntanks_both_list,zenithList,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"N$_{\mathrm{both\,hits}}$", fontsize=22)
    self.ax.set_ylabel(r"$\theta$$_{\mathrm{true}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/nTanksBothZenith"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/nTanksBothZenith"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')
 
  def plotQtotBothZenithEnergy(self):
    energyBins = 10**np.linspace(5.0, 8.0, 31)
    energyBinsShort = 10**np.linspace(5, 8.0, 7)
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    # bins = np.linspace(-1,100,102)
    bins = np.linspace(0,700,71)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    for ebin, ebinStart in enumerate(energyBinsShort[:-1]):
      lowEdge_E = energyBinsShort[ebin]
      highEdge_E = energyBinsShort[ebin+1]
      evtEBin = [ievt for ievt in self.RecoEventList if lowEdge_E <= ievt.energy < highEdge_E]
      Qtot_both_list = [np.log10(ievt.Qtot_both) for ievt in evtEBin]
      energy_list = [np.log10(ievt.energy)+9 for ievt in evtEBin]
      openingAngleList = [ievt.opening_angle*180.0/np.pi for ievt in evtEBin]
      zenithList = [ievt.zenith_true*180.0/np.pi for ievt in evtEBin]
      r_list = [ievt.r_diff for ievt in evtEBin]
      chi2_list = [ievt.chi2 for ievt in evtEBin]
      self.ax.plot(Qtot_both_list,zenithList,"o",color=colorsCustom2[ebin],label=r"$10^{{{0:.1f}-{1:.1f}}}$ eV".format(
        np.log10(energyBinsShort[ebin])+9,np.log10(energyBinsShort[ebin+1])+9),lw=2.5,alpha=0.5)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    # self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    # self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"log10(Q$_{\mathrm{tot\,both}})$", fontsize=22)
    self.ax.set_ylabel(r"$\theta$$_{\mathrm{true}}$", fontsize=22)
    # self.ax.set_yscale("log")
    # self.ax.set_xlim(0,10)
    # self.ax.set_ylim(0,100)
    # self.ax.set_aspect('equal')
    self.ax.grid(visible=True,alpha=0.5)
    self.ax.legend(fontsize=12)
    plt.savefig(plotFolder+"/QtotbothZenith"+plotSuffix+str(self.split)+"Energy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/QtotbothZenith"+plotSuffix+str(self.split)+"Energy.pdf",transparent=False,bbox_inches='tight')







  def Finish(self):
    if self.containment:
      self.RecoEventList = containedEvents(self.RecoEventList,410)


    #####################################################
    #####################################################
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    p68 = np.percentile(self.r_diffList,68)
    print(p68)
    bins = np.linspace(-1,2000,2002)
    # self.ax.hist(self.r_diffList,histtype="step",label=r"",lw=2.5)
    self.ax.hist(self.r_diffList,bins=bins,histtype="step",label=r"",lw=2.5)
    self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"abs(r$_{true}$-r$_{reco}$) [m]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_yscale("log")
    self.ax.legend()
    plt.savefig(plotFolder+"/coreDiffHG7Only"+plotSuffix+str(self.split)+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreDiffHG7Only"+plotSuffix+str(self.split)+".pdf",transparent=False,bbox_inches='tight')
    #########################################################
    #############chi2 plot###################################
    # self.plotOpenAngleEnergy()
    # self.plotOpenAngleZenith()
    #########################################
    # self.plotChi2Zenith()
    # self.plotChi2Energy()
    #####################################
    # self.plotCoreRecoX()
    # self.plotCoreRecoY()
    # self.plotCoreRecoXZenith()
    # self.plotCoreRecoYZenith()
    # self.plotCoreRecoRZenith()
    # self.plotCoreRecoREnergy()
    # self.plotPlaneRecoZenZenith()
    # self.plotPlaneRecoAziZenith()
    # self.plotPlaneRecoZenEnergy()
    # self.plotPlaneRecoAziEnergy()
    # self.plotCoreOpenAngleEnergy()
    # self.plotCoreOpenAngleZenith()
    # self.plotCoreChi2Energy()
    # self.plotCoreChi2Zenith()
    # self.plotZenCoreEnergy()
    # self.plotZenOpenAngleEnergy()
    # self.plotZenChi2Energy()
    # self.plotnTankSLCEnergyZenith()
    # self.plotQtotSLCEnergyZenith()
    # self.plotnTankHLCEnergyZenith()
    # self.plotQtotHLCEnergyZenith()
    # self.plotnTankBothEnergyZenith()
    # self.plotQtotBothEnergyZenith()
    self.plotnTankSLCZenithEnergy()
    self.plotQtotSLCZenithEnergy()
    self.plotnTankHLCZenithEnergy()
    self.plotQtotHLCZenithEnergy()
    self.plotnTankBothZenithEnergy()
    self.plotQtotBothZenithEnergy()






tray = I3Tray()
tray.AddModule("I3Reader","reader",
             # filenameList = args.input,
             filenameList = fileList,
             # filename = inFile,
            )

tray.AddModule(test7HG,"7HG",
              streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
              )

# tray.AddModule(excludeITSMT,"notITSMT",
#               streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#               )

tray.AddModule(zenithCheck,"zhits",
              split = "IceTopSplitIncl",
              containment = True,
            # streams = [icetray.I3Frame.DAQ],
            )

# tray.AddModule("I3Writer","i3writer",
#             filename=str(outputDir)+fileName+".i3.gz",
#             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )

tray.Execute()
tray.Finish()