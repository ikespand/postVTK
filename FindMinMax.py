#!/usr/bin/env python
'''
ABOUT	: A Pyhton 2.7 script to read the data from VTU file from Nektar++ CFD code.
NOTE	: RUN IT WITH PYTHON SHELL IN PARAVIEW OR pvpython/pvbatch from the terminal (only supports serial).
USAGE	: aprun -n 1 pvbatch TimeAvgEuu.py
DATE	: 31.10.2018
AUTHOR	: Sandeep Pandey (sandeep.pandey@ike.uni-stuttgart.de)
'''
import time
import numpy as np
#np.set_printoptions(threshold=np.inf) # To print out whole array, helps in debugging
import math
import glob
#import matplotlib.pyplot as plt
# Import the modules from the paraview
from paraview.simple import *
from vtk.util import numpy_support as npvtk 
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------- MAIN PROGRAM -----------------------------
ExtraInfo = True # !! USER INPUT !! # 
def log(msg,num):
	if ExtraInfo:
		print msg, num
#
# -------------------------START:Reading/reshaping -------------------------
#print'' # Empty line     
print'-----------------------------------------------------------'
print'-------------------START OF THE PROGRAM--------------------'
print'-----------------------------------------------------------'
start_time = time.time()
# List out all files with Channel_al_*.vtu. 
rSTChannelal_vtu = glob.glob('./Channel_al_*.vtu') # INPUT
print'Total number of files to be checked=',len(rSTChannelal_vtu)
wMax=[]
wMin=[]
for file_name in rSTChannelal_vtu:
	print'***********************************************************'
	print'Processing ',file_name # Print out the filename we are currently processing
	VTUFILE = XMLUnstructuredGridReader(FileName=[file_name])
	VTUFILE.PointArrayStatus = ['u', 'v', 'w', 'p']
	SetActiveSource(VTUFILE)
	dataPlane = servermanager.Fetch(VTUFILE)
	w = npvtk.vtk_to_numpy(dataPlane.GetPointData().GetArray('w'))
	log('Maximum streamwise velocity=',np.ndarray.max(w))
	log('Minimum streamwise velocity=',np.ndarray.min(w))
	wMax.append(np.ndarray.max(w))
	wMin.append(np.ndarray.min(w))
	del w
	Delete(VTUFILE)
	del VTUFILE
print'***********************************************************'
np.savetxt('wMax.csv', wMax, delimiter=',')	
np.savetxt('wMin.csv', wMin, delimiter=',')	
elapsed_time = (time.time()-start_time)/60.0
print'Total elapsed time for this code is=', "%.2f" %elapsed_time,'minutes'
print'-----------------------------------------------------------'
print'---------------------END OF THE PROGRAM--------------------'
print'-----------------------------------------------------------'