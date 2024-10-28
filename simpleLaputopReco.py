#!/usr/bin/env python3

'''
This example shows how to include IceTop DOMSets to the GCD file in addition to default InIce DOMSets
'''
from icecube.icetray import I3Tray, I3Units
from icecube import dataclasses,icetray,dataio
from icecube.trigger_sim import GetDefaultDOMSets
from icecube.trigger_sim.InjectDefaultDOMSets import InjectDefaultDOMSets
# from icecube.filterscripts.offlineL2 import LaputopStandard
from icecube.toprec import LaputopStandard
from icecube.toprec import LaputopSmallShower
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
parser.add_argument('--output',"-o",type=str,default="../testLapReco.i3.gz",help='output simulation GCD for IceTop')
args = parser.parse_args()

# GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305NoDomSetTankTrig.i3.gz"
GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HGDomsets.i3.gz"

primary = args.input[0].split("/")[-1][0]
# example python reconstructIceTopS125.py -i 
# /home/enpaudel/icecube/triggerStudy/simFiles/dataSetUnique1_6/HeDAT004171GenDetFiltProcUnique.i3.bz2
# -o /home/enpaudel/icecube/triggerStudy/simFiles/dataSetReco/FeDAT000001GenDetFiltProcUniqueCleanVEMEvts.i3.gz

plotFolder = "/home/enpaudel/icecube/filterReco/plots/"

class S125Energy(icetray.I3Module):
  def __init__(self,ctx):
    icetray.I3Module.__init__(self,ctx)

  def Configure(self):
    self.S125 = []
    self.Energy = []

  def DAQ(self,frame):
    """what if there is only one trigger in trigger hierarchy
    """
    pass
  def Physics(self,frame):
    energy = frame["MCPrimary"].energy
    params = I3LaputopParams.from_frame(frame, "LaputopParams")
    lgS125_value = params.value(Par.Log10_S125)
    beta_error = params.error(Par.Beta)
    # print("logS125",lgS125_value)
    # print("logS125",lgS125_value)
    if np.cos(frame["MCPrimary"].dir.zenith) > 0.95 and not (np.isnan(lgS125_value)):
      self.S125.append(lgS125_value)
      # self.Energy.append(np.log10(energy)*I3Units.PeV)
      self.Energy.append(np.log10(energy))

  def Finish(self):
    # print("total number of events",self.totalEvents)
    # print("s125",self.S125)
    # print("Energy",self.Energy)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # h = ax.hist2d(self.S125,self.Energy,bins=1000,density=True, norm=colors.LogNorm())
    H,xedges,yedges = np.histogram2d(self.S125,self.Energy,bins=100)
    H = H.T
    im = ax.imshow(H, interpolation="nearest",origin="lower",aspect="auto",norm=mpl.colors.LogNorm(),extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]])
    cbar = fig.colorbar(im,ax=ax,pad=0.02)
    cbar.set_label(r"count",fontsize=18)
    # cbar.ax.tick_params(labelsize=18)
    # fig.colorbar(h[3], ax=ax)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=18)
    ax.set_xlabel(r"log10(S125)", fontsize=18)
    ax.set_ylabel(r"log10(E0/GeV)", fontsize=18)
    # ax.set_yscale("log")
    ax.set_xlim(0,2.5)
    ax.set_ylim(6,8)
    plt.savefig(plotFolder+"/{0}recoS125.pdf".format(primary),transparent=False,bbox_inches='tight')
    plt.close()



tray = I3Tray()
tray.AddModule("I3Reader","reader",
             filenameList=[GCD]+args.input,
            )

#################################method 1######################
###############################################################
# tray.Add('I3TopRecoCore',
#   DataReadout     = 'OfflineIceTopHLCTankPulses',
#   NTanks          = 7,
#   ShowerCore      = 'ShowerCOG',
#   Verbose         = False,
#   Weighting_Power = 0.5,
#   If              = lambda frame: 'OfflineIceTopHLCTankPulses' in frame
#   )

# tray.Add('I3TopRecoPlane',
#   DataReadout = 'OfflineIceTopHLCTankPulses',
#   ShowerPlane = 'ShowerPlane',
#   Trigger     = 3,
#   Verbose     = False,
#   If          = lambda frame: 'OfflineIceTopHLCTankPulses' in frame
# )


# tray.Add(LaputopStandard,"Laputop",
#   pulses='OfflineIceTopHLCTankPulses'
#          )

# tray.Add(LaputopSmallShower,"LaputopSmall",
#   pulses='OfflineIceTopHLCTankPulses'
#          )
##############################################################
##############################################################


###############################################################
###########################method 1####################################
tray.AddSegment(level3_IceTop,"level3",
  detector="IC86",
  do_select=False,
  isMC=True,
  simulate_background_rate=0,
  add_jitter=False,
  snowLambda=2.6,
  icetopStreamName="IceTopSplit",
  nameChanges=[["QTriggerHierarchy","I3TriggerHierarchy"]],
  ignore_slc_calib=True
  )





# tray.AddModule(S125Energy)
# tray.AddModule("I3Writer","i3writer",
#             filename=args.output,
#             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )
####################################################################
####################################################################

tray.AddModule("I3Writer","i3writer",
            filename=args.output,
            streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
            )

tray.Execute()
tray.Finish()

