#!/usr/bin/env python3

from icecube.icetray import I3Tray, I3Units
from icecube import icetray,dataclasses,dataio

from icecube.phys_services.which_split import which_split

import glob

file_path = "/home/enpaudel/dataExp/run2023/new_filter_test/filtered/"
file_list = sorted(glob.glob(file_path+"PFFilt_PhysicsFiltering_Run00138615_Subrun*.i3.gz"))[:1]
outputDir = "/home/enpaudel/dataExp/run2023/new_filter_test/passedEvents/"
n_incl = 0

def countFilter(frame,filter_name):
    global n_incl
    if frame.Has("OfflineFilterMask"):
      mask = frame["OfflineFilterMask"]
      if mask[filter_name].condition_passed and mask[filter_name].prescale_passed:
        # print("has passed filter",frame["I3EventHeader"].run_id,frame["I3EventHeader"].event_id)
        n_incl += 1


def testFilter(frame,filter_name):
    pass_filt = False
    if frame.Has("OfflineFilterMask"):
      mask = frame["OfflineFilterMask"]
      if mask[filter_name].condition_passed and mask[filter_name].prescale_passed:
        pass_filt = True
    return pass_filt

tray = I3Tray()
tray.AddModule("I3Reader","reader",
              FilenameList=file_list,
              # filenameList=inputList[0],
              # filename=GCD,
              )
InclinedSplit = "IceTopSplitIncl"
tray.AddModule(countFilter,"count",
  filter_name="IceTopInclined_24",
  Streams=[icetray.I3Frame.Physics],
  If=which_split(split_name=InclinedSplit)
  )
tray.AddModule(testFilter,"test",
  filter_name="IceTopInclined_24",
  Streams=[icetray.I3Frame.Physics],
  # If=which_split(split_name=InclinedSplit)
  )

tray.AddModule("I3Writer","i3writer",
          filename=str(outputDir)+"inclFilt138615.i3.gz",
          streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
          DropOrphanStreams=[icetray.I3Frame.DAQ]
          )



tray.Execute()
tray.Finish()

time = 10*60+41.770381 #for run 138615
print(n_incl)
print(n_incl/time,"Hz")