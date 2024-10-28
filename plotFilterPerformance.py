#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio

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
fileDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck_afterFilter/"

# outputDir = "/data/sim/IceTop/2023/generated/untriggered/testFile/splitCheck_afterSplit/"

GCD = "/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoSMTDOMSet.i3.gz"

# plotSuffix = "Only7HG"
plotSuffix = ""
# fileList = sorted(glob.glob(fileDir+"FeDAT000*.i3.*"))[:1]
fileList = sorted(glob.glob(fileDir+"FeDAT000*.i3.*"))
# fileList = ["/home/enpaudel/icecube/triggerStudy/simFiles/dataSetClean1_6/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz"]

def openingAngle(theta1,phi1,theta2,phi2):
  return np.arccos(np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)+np.cos(theta1)*np.cos(theta2))

def test7HG(frame):
  triggerHierarchy = frame["QTriggerHierarchy"]
  icetop7HG = [t for t in triggerHierarchy if (t.key.config_id == 30043 and t.fired)]
  return len(icetop7HG)>0

def excludeITSMT(frame):
  triggerHierarchy = frame["QTriggerHierarchy"]
  itSMT = [t for t in triggerHierarchy if (t.key.config_id == 102 and t.fired)]
  return len(itSMT)<1

class zenithCheck(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)
  def Configure(self):
    self.openingAngleList = []
    self.r_diffList = []
    self.zenithTrueList = []
    self.zenithRecoList = []


  def Physics(self,frame):
    # if int(frame["I3EventHeader"].event_id) == int(8351):
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
    openAngle = openingAngle(zenith_true,azimuth_true,zenith_reco,azimuth_reco)
    # print("zenith True",zenith_reco,zenith_true,np.arcsin(np.sqrt(self.zenithBin[0])),np.arcsin(np.sqrt(self.zenithBin[1])))
    # if np.arcsin(np.sqrt(self.zenithBin[0])) <= zenith_true < np.arcsin(np.sqrt(self.zenithBin[1])) and not np.isnan(openAngle):
    self.openingAngleList.append(openAngle*I3Units.radian/I3Units.degree)
    self.r_diffList.append(r_diff)
    self.zenithTrueList.append(zenith_true*I3Units.radian/I3Units.degree)
    self.zenithRecoList.append(zenith_reco*I3Units.radian/I3Units.degree)

  def Finish(self):
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.openingAngleList = [ielt for ielt in self.openingAngleList if not np.isnan(ielt)]
    bins = np.linspace(-1,100,102)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.openingAngleList,68)
    print(p68)
    self.ax.hist(self.openingAngleList,bins=bins,histtype="step",label=r"",lw=2.5)
    self.ax.axvline(p68,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"p$_{{{:.0f}}}$={:.1f}$^{{\circ}}$".format(68,p68))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\psi$ [$^{\circ}$]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_xlim(0,100)
    self.ax.set_ylim(0.9,5*10**3)
    self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/openAngleAfterFilter"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/openAngleAfterFilter"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')
    plt.close()
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
    plt.savefig(plotFolder+"/coreDiffAfterFilter"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/coreDiffAfterFilter"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')
    plt.close()
    #######################################################
    #######################################################
    self.fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    self.ax = self.fig.add_subplot(gs[0])
    # print(self.zenithDiff,min(self.zenithDiff),max(self.zenithDiff))
    # bins = np.linspace(min(self.zenithDiff),max(self.zenithDiff),80)
    # self.ax.hist(self.zenithDiff,bins=bins,histtype="step")
    # bins = np.linspace(min(self.openingAngleList),max(self.openingAngleList),80)
    self.zenithTrueList = [ielt for ielt in self.zenithTrueList if not np.isnan(ielt)]
    bins = np.linspace(-1,100,102)
    # print(self.openingAngleList,min(self.openingAngleList),max(self.openingAngleList))
    p68 = np.percentile(self.zenithTrueList,68)
    print(p68)
    self.ax.hist(self.zenithTrueList,bins=bins,histtype="step",label=r"",lw=2.5)
    self.ax.axvline(45,ymin=0,ymax=1,color="orange",ls="--",lw=2.5,label=r"$\theta$$_{{{}}}$={:.1f}$^{{\circ}}$".format("th",45))
    self.ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    self.ax.set_xlabel(r"$\theta$ [$^{\circ}$]", fontsize=22)
    self.ax.set_ylabel("count", fontsize=22)
    self.ax.set_xlim(0,100)
    self.ax.set_ylim(0.9,5*10**3)
    self.ax.set_yscale("log")
    self.ax.legend(fontsize=18)
    plt.savefig(plotFolder+"/zenithAfterFilter"+plotSuffix+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/zenithAfterFilter"+plotSuffix+".pdf",transparent=False,bbox_inches='tight')
    plt.close()








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
            # streams = [icetray.I3Frame.DAQ],
            )

# tray.AddModule("I3Writer","i3writer",
#             filename=str(outputDir)+fileName+".i3.gz",
#             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )

tray.Execute()
tray.Finish()