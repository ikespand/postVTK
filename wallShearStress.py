#!/usr/bin/env python
'''
ABOUT	: A Pyhton 2.7 script to calculate the wall shear stress and friction Reynols number from the VTK File. 
		  It needs w_y_average variable which means dw/dy in VTK files.
		  So it is a post-processing tools for VTK which we get from averaging VTUs.
USAGE	: aprun -n 1 pvbatch TimeAvgEuu.py
DATE	: 24.10.2018
'''
import time
import numpy as np
#np.set_printoptions(threshold=np.inf) # To print out whole array, helps in debugging
import math
import glob
# Import the modules from the paraview
from paraview.simple import *
from vtk.util import numpy_support as npvtk 
paraview.simple._DisableFirstRenderCameraReset()

# USER INPUT: Provide the case specific variables
VarName ='w_y_average'
nu = 1.0/182.0
rho = 1.0
H = 1.0
# -----------------------------MY FUNCTION-----------------------------
def MeanValue(Data_y_value):
	Data_y_value=Data_y_value.transpose()
	Data_y_value = Data_y_value[np.argsort(Data_y_value[:, 0])[::1]]
	yy=Data_y_value[:,0]
	values=Data_y_value[:,1]
	sum1=0
	width=0
	j=0
	yUnique=(np.unique(yy))
	yAver=np.zeros((len(yUnique)))
	ValAver=np.zeros((len(yUnique)))
	for i in range(0,len(yy)-1):
		if  yy[i]== yy[i+1]:
			sum1 = sum1+values[i]
			width=width+1
		else:
			sum1 = sum1+values[i]
			aver = sum1/width
			yAver[j]= yy[i]
			ValAver[j]= aver
			j = j+1
			width=1
			sum1=0
	#BELOW-this will take care of last y=1
	yAver[j]= yy[i]
	ValAver[j]= sum1/width		
	#ABOBE-this will take care of last y=1			
	return yAver, ValAver
# ----------------------------- MAIN PROGRAM -----------------------------                                       
print'------------------------------------------------------------'
print'-------------------START OF THE PROGRAM---------------------'
print'------------------------------------------------------------'
start_time = time.time()
avgStatsvtk = LegacyVTKReader(FileNames=['./AvgStats.vtk'])
print'Reading is completed, now performing post-processing ..\n'
dataPlane = servermanager.Fetch(avgStatsvtk)
nbPoints = dataPlane.GetNumberOfPoints()
# Put the points coordinates in x, y and z arrays
x=[]
y=[]
z=[]
for i in range(nbPoints):
	coord = dataPlane.GetPoint(i)
	xx, yy, zz = coord[:3]
	x.append(xx)
	y.append(yy)
	z.append(zz)
	del xx, yy, zz
x=np.array(x)
y=np.array(y)
z=np.array(z)

wAvg = npvtk.vtk_to_numpy(dataPlane.GetPointData().GetArray(VarName))	
wAvg1 = np.zeros((2,nbPoints))
wAvg1=y,wAvg
wAvg1=np.array(wAvg1)
w_average=MeanValue(wAvg1)

dwdy=w_average[1][0] #Taking value at the wall
tau_w=(nu/rho)*dwdy
u_tau=np.sqrt(tau_w/rho)
Re_tau=u_tau*H/nu;
print'Re_tau',Re_tau,'\ntau_w',tau_w,'\nu_tau',u_tau
np.savetxt('w_y_average.csv',w_average, delimiter=',')
print'\nSuccessfully finished!'
elapsed_time = (time.time()-start_time)/60.0
print'Total elapsed time for this code is=', "%.2f" %elapsed_time,'minutes'
print'------------------------------------------------------------'
print'---------------------END OF THE PROGRAM---------------------'
print'------------------------------------------------------------'