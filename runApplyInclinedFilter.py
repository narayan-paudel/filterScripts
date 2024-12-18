#!/usr/bin/env python3

import os
import glob
import subprocess

import re

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))
# ABS_PATH_HERE = "./"
ABS_PATH_HERE += "/"
# print("abs path",ABS_PATH_HERE)

############################################################################
inputPath = "/data/sim/IceTop/2023/generated/untriggered/filterStudy/splitted/"
inputList = sorted(glob.glob(inputPath+"*.i3.gz"))[:]
#############################################################################
# print("inputList",inputList)



submitFileName = ABS_PATH_HERE+"tempSubmitInclFilt.sub"


def makeSubFile(fileList):
  submitFile = open(submitFileName,"w")
  submitFile.write("########################################\n")
  submitFile.write("## submit description file\n")
  submitFile.write("########################################\n\n")
  submitFile.write("Universe   = vanilla\n")
  submitFile.write("Executable = /data/user/enpaudel/filterReco/scripts/applyInclinedFilter.sh\n")
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
  print(fileList[0])
  filename = fileList[0].split("/")[-1]
  primary = filename.split("DAT")[0]
  print(primary)
  runID = re.findall(r'\d+', filename)[0]
  print(runID)

  # subrunID = fileList[0].split("/")[-1].split("_")[4].split(".")[0]
  # print(runID,subrunID)
  # corsikaID = int(''.join(i for i in corsikaID if i.isdigit()))
  # print("file list",*fileList[:2])
  # print("corsika id",corsikaID,primary)
  subprocess.call(["condor_submit tempSubmitInclFilt.sub -batch-name {0}--{1}".format(primary,runID)], shell=True)
  subprocess.call(["rm tempSubmitInclFilt.sub"], shell=True)

def submitToCondor(fileList,chunk):
  fileChunks = [fileList[i:i + chunk] for i in range(0, len(fileList), chunk)]
  for ifileList in fileChunks:
    print("submitting to condor",ifileList[0])
    submitToCondorFile(ifileList)

submitToCondor(inputList,1)