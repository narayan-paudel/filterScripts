#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2024 The IceTray Contributors
#
# SPDX-License-Identifier: BSD-2-Clause

"""
This script will reconstruct the Icetop parameters using the rock_bottom framework
First the shower core is estimated using the center-of-gravity and the direction using a plane
The fit is then performed using the log-log LDF 
It is not supposed to represent an "official" reconstruction, just something
to get you started on your own analysis
"""

## python -------------------------------
import argparse

## IceCube -------------------------------
from icecube.icetray import I3Tray, I3Units
from icecube import icetray, dataclasses, dataio
from icecube import toprec, phys_services
from icecube import rock_bottom
from icecube import gulliver, gulliver_modules, lilliput
from icecube.rock_bottom import RbPair, RbList
from icecube.rock_bottom.modules import WriteSeed

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, nargs="+", help="data files")
parser.add_argument("--gcd", required=True, type=str, help="output file name")
parser.add_argument("--output", type=str, default="Ouput_IceTop_1Step.i3.gz", help="output file")
args = parser.parse_args()

tray = I3Tray()

# ---------------------------------------------------------------------
#                 Reader
# ---------------------------------------------------------------------

tray.Add("I3Reader", FilenameList=[args.gcd] + args.input)

# ---------------------------------------------------------------------
#            IceTop initials for reconstruction
# ---------------------------------------------------------------------

loglog = rock_bottom.ldf.LogLog()
TimeFcn = rock_bottom.showerfront.TimeGauss()
HLCPulses = "OfflineIceTopHLCTankPulses"
SnowService = "SimpleSnow"

tray.AddService(
    "I3SimpleSnowAttenuationServiceFactory",
    SnowService,
    Lambda=2.1,
)

tray.AddModule(
    "I3TopRecoCore",
    "core",
    DataReadout=HLCPulses,
    NTanks=5,
    If=lambda frame: "ShowerCOG" not in frame,
)

tray.AddModule(
    "I3TopRecoPlane",
    "plane",
    DataReadout=HLCPulses,
    If=lambda frame: "ShowerPlane" not in frame,
)

tray.AddModule(
    WriteSeed,
    "WriteSeed",
    ShowerCOG="ShowerCOG",
    ShowerPlane="ShowerPlane",
    OutputName="ShowerSimpleSeedIT",
)

# ---------------------------------------------------------------------
#          Minimizer preparation
# ---------------------------------------------------------------------

tray.AddService(
    "I3GulliverMinuitFactory",
    "Minuit",
    MinuitPrintLevel=-2,
    FlatnessCheck=True,
    Algorithm="SIMPLEX",
    MaxIterations=2500,
    MinuitStrategy=2,
    Tolerance=0.001,
)

# ---------------------------------------------------------------------
#            1 STEP
# ---------------------------------------------------------------------

tray.AddService(
    "LaputopSignalModel",
    "LaputopSignalModel",
    LDF=loglog,
    SnowService=SnowService,
    ParameterNames=RbList([RbPair("slopeLdf", 0)]),
    ParameterValues=[
        2.6,
    ],
    BoundNames=RbList([RbPair("slopeLdf", 0), RbPair("lgSref", 0)]),
    BoundValues=[(2.9, 3.1), (-3.0, 8.0)],
    StepNames=RbList([RbPair("slopeLdf", 0), RbPair("lgSref", 0)]),
    StepValues=[
        0.6,
        1.0,
    ],
)

tray.AddService(
    "I3RbLDFLikelihoodFactory",
    "IceTopLikelihood_LDF",
    DetectorType=rock_bottom.IceTop,
    Model="LaputopSignalModel",
    Pulses1=HLCPulses,
    UseSilent=True,
    MinSignal=0.7,
)

# ---------------------------------------------------------------------
#          Run combined reconstruction
# ---------------------------------------------------------------------

tray.AddService(
    "I3MultiSurfaceSeedServiceFactory",
    "SeedsStep1",
    FirstGuesses=["ShowerSimpleSeedIT"],
    SignalModels=["LaputopSignalModel"],
    ParticleXStep=10.0,
    ParticleYStep=10.0,
    ParticleXRelBounds=[-200.0, 200.0],
    ParticleYRelBounds=[-200.0, 200.0],
)

tray.AddService(
    "I3MultiSurfaceParametrizationFactory",
    "ParamsStep1",
    SeedService="SeedsStep1",
)


tray.AddModule(
    "I3SimpleFitter",
    "CombinedReconstructionStep1",
    SeedService="SeedsStep1",
    Parametrization="ParamsStep1",
    LogLikelihood="IceTopLikelihood_LDF",
    OutputName="ITReconstructionStep1a",
    Minimizer="Minuit",
)

# ---------------------------------------------------------------------
#           Writing
# ---------------------------------------------------------------------

tray.AddModule(
    "I3Writer",
    "i3writer",
    Filename=args.output,
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
)

tray.Execute()
tray.Finish()
