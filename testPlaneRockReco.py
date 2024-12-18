#!/usr/bin/env python3

'''
This example shows how to include IceTop DOMSets to the GCD file in addition to default InIce DOMSets
'''
from icecube.icetray import I3Tray, I3Units
from icecube import dataclasses,icetray,dataio


import numpy as np



import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input',"-i",nargs="+",type=str,default="",help="input simulation GCD for IceTop")
args = parser.parse_args()


GCD="/data/user/enpaudel/triggerStudy/simFiles/GeoCalibDetectorStatus_2020.Run135057.Pass2_V0_Snow210305_W7HGDomsets.i3.gz"

def recoCompare(frame):
  if frame.Stop == icetray.I3Frame.Physics:
    rockReco = frame["ITReconstructionStep1a"].dir.zenith
    planeReco = frame["ShowerPlane"].dir.zenith
    lapuReco = frame["Laputop"].dir.zenith
    if not (np.isnan(rockReco) and np.isnan(planeReco) and np.isnan(lapuReco)):
      print(f"zenith reco rock {np.rad2deg(rockReco):.1f} plane {np.rad2deg(planeReco):.1f} lapu {np.rad2deg(lapuReco):.1f}")


tray = I3Tray()
tray.AddModule("I3Reader","reader",
             filenameList=[GCD]+args.input,
            )

tray.AddModule(recoCompare)


# tray.AddModule("I3Writer","i3writer",
#             filename=args.output,
#             streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
#             )

tray.Execute()
tray.Finish()
