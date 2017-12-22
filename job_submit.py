import sys
import os
import time

theta_range = 90
theta_int = 10
phi_range = 90
phi_int = 10

WCSim_exe_file = "/neut/neut21/tknight/Documents/nuPRISM_Code/WCSim/exe/bin/Linux-g++/WCSim"
job_macro_path = "/neut/neut21/tknight/Documents/nuPRISM_Code/WCSim/job_macros/"

output_path = "/neut/datasrv2a/tknight/WCSim/job_outputs/"
error_path = "/neut/datasrv2a/tknight/WCSim/job_errors/"

WCSim_dir = "/neut/neut21/tknight/Documents/nuPRISM_Code/WCSim/"
source_file = "/neut/neut21/tknight/Documents/nuPRISM_Code/Analysis/Source_At_Start_nuPRISM.sh"

reflector = raw_input("refl or norefl or 3point5inch\n")
while(reflector != "refl" and reflector != "norefl" and reflector != "3point5inch"):
  reflector = raw_input("refl or norefl or 3point5inch\n")

reflector_change = raw_input("If using default reflectors, enter 1. If using different reflectors, enter 2.\n")
while(reflector_change != "1" and reflector_change != "2"):
  reflector_change = raw_input("If using default reflectors, enter 1. If using different reflectors, enter 2.\n")
  
if(reflector_change == "2"):
  reflector_angle = raw_input("Enter the refector angle in deg\n")
  reflector_height = raw_input("Enter the reflector height in mm (note: decimal places are written as a p (ex. 13p5))\n")

count = 0

for theta in range(-theta_range,theta_range+1,theta_int):
#  for phi in range(-90,-80,5):

  print "Job for theta = " + str(theta)
#  scriptname = "/neut/neut21/tknight/Documents/nuPRISM_Code/WCSim/job_macros/x_mPMT_NP1_" + str(theta) + "_" + str(phi) + ".mac"
  if(reflector_change == "2"):
    scriptname = "job_x_run_" + reflector + "_" + reflector_angle + "deg_" + reflector_height + "mm_" + str(theta) + ".pbs"
  else:
    scriptname = "job_x_run_" + reflector + "_" + str(theta) + ".pbs"

  fout = open(scriptname, 'w')
  fout.write("#!/bin/bash\n")

  if(reflector_change == "2"):
    fout.write("#PBS -N tknight_job" + str(count) + "_" + reflector + "_" + reflector_angle + "deg_" + reflector_height + "mm\n")
  else:
    fout.write("#PBS -N tknight_job" + str(count) + "_" + reflector + "\n")
#  fout.write("#PBS -j oe\n")
  if(reflector_change == "2"):
    fout.write("#PBS -o " + output_path + "tknight_job"+ str(count) + "_" + reflector + "_" + reflector_angle + "deg_" + reflector_height + "mm_" + str(theta) )
    fout.write("#PBS -e " + error_path + "tknight_job" + str(count) + "_" + reflector + "_" + reflector_angle + "deg_" + reflector_height + "mm_" + str(theta) )
  else:
    fout.write("#PBS -o " + output_path + "tknight_job"+ str(count) + "_" + reflector + "_" + str(theta))
    fout.write("#PBS -e " + error_path + "tknight_job" + str(count) + "_" + reflector + "_" + str(theta) )
  fout.write("#PBS -l walltime=72:00:00,mem=2000mb,vmem=10000mb\n")
#  fout.write("export DIRNAME='echo $HOST | sed -e 's/neut/data/'`a\n")
  fout.write("cd " + WCSim_dir + "\n")

  fout.write("source " + source_file + "\n")
  
  for phi in range(-phi_range,phi_range+1,phi_int):
    if(reflector_change == "2"):
      fout.write(WCSim_exe_file + " " + job_macro_path + "x_mPMT_NP1_" + reflector + "_" + reflector_angle + "deg_" + reflector_height + "mm_" + str(theta) + "_" + str(phi) + ".mac > /dev/null \n")
    else:
      fout.write(WCSim_exe_file + " " + job_macro_path + "x_mPMT_NP1_" + reflector + "_" + str(theta) + "_" + str(phi) + ".mac > /dev/null \n")
#  fout.write("/neut/neut21/tknight/Documents/nuPRISM_Code/WCSim/exe/bin/Linux-g++/WCSim /neut/neut21/tknight/Documents/nuPRISM_Code/WCSim/job_macros/x_mPMT_NP1_" + reflector + "_" + str(theta) + ".mac\n")
#  fout.write("mv x_mPMT_NP1_refl_*.root /neut/data17/tknight/WCSim/z_NP1_jobs/\n")
  fout.close()
  os.system("qsub " + scriptname + " -q srvq")
  os.remove(scriptname)
  count = count + 1  

