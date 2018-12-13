#!/usr/bin/env python
'''
ABOUT	: A Pyhton 2.7 script to read the data from 'spectrum' folder generated for different time step. It will do further averaging.
		  Currently, it saves wavenumber and energy in two different csv files, which can be used for contourmap. 
		  Only averaging of Euu is performed, for k and y, it is not needed because they are fix.
USAGE	: aprun -n 1 python2 AveragingSpectrum.py
DATE	: 11.12.2018
AUTHOR	: Sandeep Pandey (sandeep.pandey@ike.uni-stuttgart.de)
'''
from __future__ import print_function
import time
import sys
import numpy as np
#np.set_printoptions(threshold=np.inf) # To print out whole array, helps in debugging
import math
import os
import glob

Dir='spectrum_uu/' # !! USER INPUT !! # 

try:
    # Create target Directory
    os.path.exists(Dir)
    print("Processing files...") 
except OSError:
	print("Directory does not exist")
		
# ----------------------------- MAIN PROGRAM -----------------------------
fileDelta=1		# !! USER INPUT !! # 
startFile=91   	# !! USER INPUT !! # 
# -------------------------START:Reading/reshaping -------------------------
start_time = time.time()
print("-----------------------------------------------------------")
print("-------------------START OF THE PROGRAM--------------------")
print("-----------------------------------------------------------")
numFiles = len(glob.glob(Dir+'*_E.csv'))
print("Number of files for averaging=",numFiles)
tempEuu=np.genfromtxt(Dir+'Channel_al_'+str(startFile)+'_E.csv', delimiter=',') 
Euu=np.zeros(tempEuu.shape)
for i in range(numFiles):
	timeStamp=(fileDelta*i)+startFile						
	filename='Channel_al_%s'% str(timeStamp)
	EuuFile = np.genfromtxt(Dir+filename+'_E.csv', delimiter=',')
	Euu=Euu+EuuFile	
	del EuuFile
	
EuuAvg=Euu/numFiles
currentDir = os.getcwd()
OutFileDir = os.path.join(currentDir,Dir)
OutFile='Avg'
#ksName=os.path.join(OutFileDir, OutFile+'_kxAvg.csv')
EuuName=os.path.join(OutFileDir, OutFile+'_EAvg.csv')
yName=os.path.join(OutFileDir, OutFile+'_yAvg.csv')
#np.savetxt(ksName,kxS, delimiter=',')
np.savetxt(EuuName,EuuAvg, delimiter=',')
#np.savetxt(yName,yS, delimiter=',')
elapsed_time = (time.time()-start_time)/60.0
print("Total elapsed time for this code is=", "%.2f" %elapsed_time,"minutes")
print("-----------------------------------------------------------")
print("---------------------END OF THE PROGRAM--------------------")
print("-----------------------------------------------------------")