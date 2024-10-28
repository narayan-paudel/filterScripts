#!/usr/bin/env python3



from icecube import icetray

@icetray.traysegment
def apply_cluster_cleaning(tray,name="",
  SubEventStreamName="IceTopSplitIncl",
  InputPulses="IceTopTankPulses",
  OutputPulses="CleanedTankPulses",
  BadTankList="TankPulseMergerExcludedTanksIncl",
  ExcludedTanks="ClusterCleaningExcludedTanksIncl",
  InterStationTimeTolerance=200.0 * I3Units.ns,  # Default
  IntraStationTimeTolerance=200.0 * I3Units.ns,  # Default
  ):
  from icecube.icetray import I3Tray, I3Units
  from icecube import icetray, dataclasses, dataio
  from icecube import topeventcleaning

  import numpy as np
  import glob

  from icecube.phys_services.which_split import which_split

  tray.AddModule("I3TopHLCClusterCleaning", name + "_cluster",
                 SubEventStreamName=SubEventStreamName,
                 InputPulses=InputPulses,
                 OutputPulses=OutputPulses,
                 BadTankList=BadTankList,
                 ExcludedTanks=ExcludedTanks,
                 InterStationTimeTolerance=InterStationTimeTolerance,
                 IntraStationTimeTolerance=IntraStationTimeTolerance,
                 # If=lambda frame: "IceTopTankPulses" in frame and len(frame["IceTopTankPulses"])>0,
                 If=lambda frame: InputPulses in frame,
                 )
  class PhysicsCopyTriggers(icetray.I3ConditionalModule):
      def __init__(self, context):
          icetray.I3ConditionalModule.__init__(self, context)
          self.AddOutBox("OutBox")

      def Configure(self):
          pass

      def Physics(self, frame):
          if frame.Has("QTriggerHierarchy"):
              myth = dataclasses.I3TriggerHierarchy.from_frame(frame,
                                                               "QTriggerHierarchy")
              if frame.Has("TriggerHierarchy"):
                  icetray.logging.log_error("ERROR: PhysicsCopyTriggers: triggers in frame")
              else:
                  frame.Put("TriggerHierarchy", myth)
          else:
              icetray.logging.log_error("Error: PhysicsCopyTriggers: Missing QTriggerHierarchy")
          self.PushFrame(frame)

  # Null split and IceTop split need the triggerHierarchy put back in P frame

  tray.AddModule(PhysicsCopyTriggers, name + "_IceTopTrigCopy",
                 If=which_split(split_name=SubEventStreamName))


