#!/usr/bin/env python3

import os
import glob
import subprocess

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
# print("abs path",ABS_PATH_HERE)

############################################################################
inputPath = "/home/enpaudel/dataExp/run2023/IceTopTrig/"
# inputPath = "/data/ana/CosmicRay/IceTop_level3/sim/IC86.2012/"
# OinputList = sorted(glob.glob(inputPath+"12631/"+"Level3_IC86.2012_12631_*.i3.gz"))[:3000]
# pinputList = sorted(glob.glob(inputPath+"12360/"+"Level3_IC86.2012_12360_*.i3.gz"))[:3000]
# HeinputList = sorted(glob.glob(inputPath+"12630/"+"Level3_IC86.2012_12630_*.i3.gz"))[:3000]
inputList = sorted(glob.glob(inputPath+"PFFilt_PhysicsFiltering_Run00138615_Subrun*.i3.gz"))[:1]
#############################################################################
# print("inputList",inputList)



submitFileName = ABS_PATH_HERE+"tempSubmitInclFilt2023.sub"


def makeSubFile(fileList):
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /data/user/enpaudel/filterReco/scripts/applyInclinedFilterData.sh\n")
  submitFile.write("Log        = /scratch/enpaudel/log/inclFilt$(Process).log\n")
  submitFile.write("Output     = /data/user/enpaudel/filterReco/log/inclFilt$(Process).out\n")
  submitFile.write("Error      = /data/user/enpaudel/filterReco/log/inclFilt$(Process).err\n")
  submitFile.write("request_cpus = 1\n")
  submitFile.write("request_memory = 2GB\n")
  submitFile.write("request_disk = 1GB\n")
  submitFile.write("#request_gpus = 1\n")
  submitFile.write("#should_transfer_files   = IF_NEEDED\n")
  submitFile.write("#when_to_transfer_output = ON_EXIT\n")
  submitFile.write("#notification = Complete\n")
  submitFile.write("#notify_user = <email-address>\n")
  submitFile.write("#priority = <integer>\n")
  submitFile.write("##long job\n")
  priority=9900
  # submitFile.write("+AccountingGroup=\"1_week.$ENV(USER)\" \n\n")
  submitFile.write("priority = {}\n".format(priority))
  submitFile.write("#set arguments to executable\n")
  submitFile.write("arguments = ")
  for ifile in fileList:
    submitFile.write(ifile + " ")
  submitFile.write("\n")
  submitFile.write("queue 1\n")
  submitFile.close()

def submitToCondorFile(fileList):
  makeSubFile(fileList)
  # print(fileList[0])
  runID = fileList[0].split("/")[-1].split("_")[2].split("n00")[1]
  subrunID = fileList[0].split("/")[-1].split("_")[4].split(".")[0]
  # print(runID,subrunID)
  # corsikaID = int(''.join(i for i in corsikaID if i.isdigit()))
  # print("file list",*fileList[:2])
  # print("corsika id",corsikaID,primary)
  subprocess.call(["condor_submit tempSubmitInclFilt2023.sub -batch-name {0}---{1}".format(runID,subrunID)], shell=True)
  subprocess.call(["rm tempSubmitInclFilt2023.sub"], shell=True)

def submitToCondor(fileList,chunk):
  fileChunks = [fileList[i:i + chunk] for i in range(0, len(fileList), chunk)]
  for ifileList in fileChunks:
    print("submitting to condor",ifileList[0])
    submitToCondorFile(ifileList)

submitToCondor(inputList,1)
