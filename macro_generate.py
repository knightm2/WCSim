import sys
import os
import time
import math

theta_range = 90
theta_int = 10
phi_range = 90
phi_int = 10

s_theta_range = 60
s_theta_int = 4
s_phi_range = 60
s_phi_int = 4

output_filepath = "/neut/datasrv2a/tknight/WCSim/z_NP1_jobs_Dec13_COC/" #directory to which generated root files will be saved

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
  for phi in range(-phi_int,phi_int+1,phi_int):

    print reflector

    if(reflector_change == "2"):
      macroname = "job_macros/x_mPMT_NP1_" + reflector + "_" + reflector_angle + "deg_" + reflector_height + "mm_" + str(theta) + "_" + str(phi) + ".mac"
      print reflector_angle
      print reflector_height
    else:
      macroname = "job_macros/x_mPMT_NP1_" + reflector + "_" + str(theta) + "_" + str(phi) + ".mac"
    fout = open(macroname, 'w')
    
    fout.write("/run/verbose 0\n")
    fout.write("/tracking/verbose 0\n")
    fout.write("/hits/verbose 0\n")

    if(reflector_change == "2"):
      fout.write("/control/execute /neut/neut21/tknight/Documents/nuPRISM_Code/WCSim/mPMT_NP1_" + reflector + "_" + reflector_angle + "deg_" + reflector_height + "mm.mac\n")
    else:
      fout.write("/control/execute /neut/neut21/tknight/Documents/nuPRISM_Code/WCSim/mPMT_NP1_" + reflector + ".mac\n")

    fout.write("/WCSim/PMTQEMethod DoNotApplyQE\n")
    fout.write("/WCSim/PMTCollEff off\n")
    fout.write("/WCSim/SavePi0 false\n")
    fout.write("/DAQ/Digitizer SKI\n")
    fout.write("/DAQ/Trigger NDigits\n")
    fout.write("/control/execute daq.mac\n")
    fout.write("/DarkRate/SetDarkRate 0 kHz\n")
    fout.write("/DarkRate/SetDarkMode 19\n")
    fout.write("/DarkRate/SetDarkHigh 100000\n")
    fout.write("/DarkRate/SetDarkLow 0\n")
    fout.write("/DarkRate/SetDarkWindow 4000\n")
#    fout.write("/mygen/generator laser\n")
    fout.write("/mygen/generator gps\n")
    fout.write("/gps/particle opticalphoton\n")
    fout.write("/gps/ene/mono 3.102 eV\n")
    fout.write("/gps/pos/type Point\n")
#    fout.write("/gps/pos/shape Circle\n")
#    fout.write("/gps/pos/radius 0.001 mm\n")
#    fout.write("/gps/pos/radius 5.0 mm\n") 
    fout.write("/gps/ang/type iso\n")

    print "theta = " + str(theta) + "  |  phi = " + str(phi)

    n_photons = 250    
    count = count + 1
    pi = 4.0*math.atan(1)
    theta_r = theta*(pi/180)
    phi_r = phi*(pi/180)
#    CEN1 = -3.617+math.cos(phi_r)*math.cos(theta_r)
#    CEN1 = -3.7+math.cos(phi_r)*math.cos(theta_r)
    ROC = 1.229 #Radius of curvature of source movement
    CEN1 = -3.7-0.229 + ROC*math.cos(phi_r)*math.cos(theta_r)
    CEN2 = ROC*math.cos(phi_r)*math.sin(theta_r)
    CEN3 = ROC*math.sin(phi_r)
    ROT11 = math.sin(phi_r)*math.cos(theta_r)
    ROT12 = math.sin(phi_r)*math.sin(theta_r)
    ROT13 = -math.cos(phi_r)
    ROT21 = -math.sin(theta_r)/math.cos(theta_r)
    ROT22 = 1
    ROT23 = 0
    DIR1 = -math.cos(phi_r)*math.cos(theta_r)
    DIR2 = -math.cos(phi_r)*math.sin(theta_r)
    DIR3 = -math.sin(phi_r)

    fout.write("/gps/pos/centre " + str(CEN1) + " " + str(CEN2) + " " + str(CEN3) + " m\n")
    fout.write("/gps/pos/rot1 " + str(ROT11) + " " + str(ROT12) + " " + str(ROT13) + "\n")
    fout.write("/gps/pos/rot2 " + str(ROT21) + " " + str(ROT22) + " " + str(ROT23) + "\n")
#    fout.write("/gps/direction " + str(DIR1) + " " + str(DIR2) + " " + str(DIR3) + "\n")
    for s_theta in range(-s_theta_range,s_theta_range+1,s_theta_int):
      for s_phi in range(-s_phi_range,s_phi_range+1,s_phi_int):

        ang_span_theta = float(s_theta_int) #degrees spanned in theta and phi for each cone of light
        ang_span_phi = float(s_phi_int)

        s_theta_r = float(s_theta)*pi/180.0
        s_phi_r = float(s_phi)*pi/180.0

        #scanning across source phi and theta angles 
        s_theta_low = math.acos(-DIR3) + s_theta_r - (ang_span_theta/2.0)*pi/180.0
        s_theta_up = s_theta_low + ang_span_theta*pi/180.0
        s_phi_low = math.atan(DIR2/DIR1) + s_phi_r - (ang_span_phi/2.0)*pi/180.0
        s_phi_up = s_phi_low + ang_span_phi*pi/180.0

        fout.write("/gps/ang/mintheta " + str(s_theta_low) + " rad\n")
        fout.write("/gps/ang/maxtheta " + str(s_theta_up) + " rad\n")
        fout.write("/gps/ang/minphi " + str(s_phi_low) + " rad\n")
        fout.write("/gps/ang/maxphi " + str(s_phi_up) + " rad\n")
#        fout.write("/gps/ang/mintheta -0.573 deg\n")
#        fout.write("/gps/ang/maxtheta 0.573 deg\n")  
#        fout.write("/gps/ang/minphi -0.573 deg\n")
#        fout.write("/gps/ang/maxphi 0.573 deg\n")  
        
        
#        y_distance = y_steps*0.020 #taking each step as 20mm
#        z_distance = z_steps*0.020
#        x_distance = 0
#        l = math.sqrt(y_distance**2 + z_distance**2)
#        if(l<0.254):
#          x_distance = (0.332 - math.sqrt(l**2 + 0.214**2))*math.cos(math.atan(l/0.214))

#        DIR1 = x_distance-math.cos(phi_r)*math.cos(theta_r)
#        DIR1 = -math.cos(phi_r)*math.cos(theta_r)
#        DIR2 = y_distance-math.cos(phi_r)*math.sin(theta_r)
#        DIR3 = z_distance-math.sin(phi_r)
#        LEN = math.sqrt(DIR1**2+DIR2**2+DIR3**2)
#        DIR1 = DIR1/LEN
#        DIR2 = DIR2/LEN
#        DIR3 = DIR3/LEN
    
        ang_cov_factor = (math.cos(s_theta_low)-math.cos(s_theta_up))/(math.cos((90.0-ang_span_theta/2.0)*pi/180.0)-math.cos((90.0+ang_span_theta/2.0)*pi/180.0))   
        n_counts = abs(int(250.0*ang_cov_factor))
        fout.write("\n")
#        fout.write("/gps/direction " + str(DIR1) + " " + str(DIR2) + " " + str(DIR3) + "\n")
        if(reflector_change == "2"):
          fout.write("/WCSimIO/RootFile " + output_filepath + reflector + "_" + str(theta) + "/x_mPMT_NP1_" + reflector + "_" + reflector_angle + "deg_" + reflector_height + "mm_" + str(theta) + "_" + str(phi) + "_" + str(s_theta) + "_" + str(s_phi) + ".root\n")
        else:
          fout.write("/WCSimIO/RootFile " + output_filepath + reflector + "_" + str(theta) + "/x_mPMT_NP1_" + reflector + "_" + str(theta) + "_" + str(phi) + "_" + str(s_theta) + "_" + str(s_phi) + ".root\n")
#        fout.write("/run/beamOn 5000\n")
        fout.write("/run/beamOn " + str(n_counts) + "\n")
        fout.write("\n")

    fout.close()

#    scriptname = "x_mPMT_" + str(theta) + "_" + str(phi) + ".pbs"
#    gout = open(scriptname,'w')
#    gout.write("#!/bin/bash\n")
#    print "Looking at file " + scriptname
#    gout.write("#PBS -N tknight_job" + str(count) + "\n")
#    gout.write("#PBS -j oe\n")
#    gout.write("#PBS -l walltime=48:00:00,mem=2000mb,vmem=2400mb\n")
#    gout.write("export DIRNAME=`echo $HOST | sed -e 's/neut/data/'`a\n")
#    gout.write("source /neut/data21/tknight/Documents/nuPRISM_Code/Analysis/Source_At_Start_nuPRISM.sh")
#    gout.write(






