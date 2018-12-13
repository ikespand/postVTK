#!/usr/bin/env python
'''
ABOUT	: A Pyhton 2.7 script to read the data from VTU file from Nektar++ CFD code for DNS. 
		  Then it computes the RSS and TKE and finally write down a csv.
NOTE	: RUN IT WITH PYTHON SHELL IN PARAVIEW OR pvpython/pvbatch from the terminal (only supports serial).
USAGE	: aprun -n 1 pvbatch TimeAvgEuu.py
DATE	: 24.10.2018
'''
import time
#np.set_printoptions(threshold=np.inf) # To print out whole array, helps in debugging
import math
import glob
# Import the modules from the paraview
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
# ----------------------------- MAIN PROGRAM -----------------------------                                       
print'------------------------------------------------------------'
print'-------------------START OF THE PROGRAM---------------------'
print'-------------------TURBULENCE STATISTICS--------------------'
print'------------------------------------------------------------'
start_time = time.time()
rSTChannelal_136 = glob.glob('./Channel_al_*.vtu')
numFiles=len(rSTChannelal_136)
print'Total number of files to be processed= ',numFiles,'and now reading them'
rSTChannelal_136 = XMLUnstructuredGridReader(FileName=rSTChannelal_136)
rSTChannelal_136.PointArrayStatus = ['u', 'v', 'w', 'w_y']
print'Reading is completed, now performing post-processing ..'

# create a new 'Calculator'
calculator1 = Calculator(Input=rSTChannelal_136)
calculator1.ResultArrayName = 'uu'
calculator1.Function = 'u*u'

# create a new 'Calculator'
calculator2 = Calculator(Input=calculator1)
calculator2.ResultArrayName = 'vv'
calculator2.Function = 'v*v'

# create a new 'Calculator'
calculator3 = Calculator(Input=calculator2)
calculator3.ResultArrayName = 'ww'
calculator3.Function = 'w*w'

# create a new 'Calculator'
calculator4 = Calculator(Input=calculator3)
calculator4.ResultArrayName = 'vw'
calculator4.Function = 'v*w'

# create a new 'Temporal Statistics'
temporalStatistics1 = TemporalStatistics(Input=calculator4)
# Properties modified on temporalStatistics1
temporalStatistics1.ComputeMinimum = 0
temporalStatistics1.ComputeMaximum = 0
temporalStatistics1.ComputeStandardDeviation = 0

# create a new 'Calculator'
calculator5 = Calculator(Input=temporalStatistics1)
calculator5.ResultArrayName = 'uuprime'
calculator5.Function = 'uu_average-(u_average*u_average)'

# create a new 'Calculator'
calculator6 = Calculator(Input=calculator5)
calculator6.ResultArrayName = 'vvprime'
calculator6.Function = 'vv_average-(v_average*v_average)'

# create a new 'Calculator'
calculator7 = Calculator(Input=calculator6)
calculator7.ResultArrayName = 'wwprime'
calculator7.Function = 'ww_average-(w_average*w_average)'

# create a new 'Calculator'
calculator8 = Calculator(Input=calculator7)
calculator8.ResultArrayName = 'vwprime'
calculator8.Function = 'vw_average-(v_average*w_average)'

# set active source
SetActiveSource(calculator8)
# save data
SaveData('./AvgStats.vtk', proxy=calculator8)
print'Successfully finished!'
elapsed_time = (time.time()-start_time)/60.0
print'Total elapsed time for this code is=', "%.2f" %elapsed_time,'minutes'
print'------------------------------------------------------------'
print'-------------------TURBULENCE STATISTICS--------------------'
print'---------------------END OF THE PROGRAM---------------------'
print'------------------------------------------------------------'