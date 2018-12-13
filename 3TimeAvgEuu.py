#!/usr/bin/env python
'''
ABOUT	: A Pyhton 2.7 script to read the data from VTU file from Nektar++ CFD code for channel flow DNS. 
		  After reading, we compute 2-D energy spectra from the fluctuations for a given timestep.
		  Currently, it saves wavenumber and energy in two different csv files, which can be used for contourmap.
NOTE	: RUN IT WITH PYTHON SHELL IN PARAVIEW OR pvpython/pvbatch from the terminal (only supports serial).
USAGE	: aprun -n 1 pvbatch TimeAvgEuu.py
TODO	: 1. A new interface to draw the map from this code itself. 
		  2. parallelize it.
DATE	: 03.08.2018
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
#import matplotlib.pyplot as plt
# Import the modules from the paraview
from paraview.simple import *
from vtk.util import numpy_support as npvtk 
paraview.simple._DisableFirstRenderCameraReset()
newDir='spectrum'
try:
    # Create target Directory
    os.mkdir(newDir)
    print("Directory spectrum  Created ") 
except OSError as e:
    if e.errno != os.errno.EEXIST:
		raise   
		
# ----------------------------- MAIN PROGRAM -----------------------------
HomeMode=384
# NEED EXTRA INFO?
ExtraInfo = False # !! USER INPUT !! # 
def log(msg,num):
	if ExtraInfo:
		print(msg, num)
#
# -------------------------START:Reading/reshaping -------------------------
#print'' # Empty line     
print("-----------------------------------------------------------")
print("-------------------START OF THE PROGRAM--------------------")
print("-----------------------------------------------------------")
start_time = time.time()
# List out all files with RSTChannel-al_*.vtu. 
rSTChannelal_vtu = sys.argv[1:]#glob.glob('./RSTChannel-al_*.vtu') # INPUT
print("Total number of files to be processed=",len(rSTChannelal_vtu))
StartFile = XMLUnstructuredGridReader(FileName=[rSTChannelal_vtu[0]])
SetActiveSource(StartFile)
dataPlane = servermanager.Fetch(StartFile)
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
yUnique=(np.unique(y))
NAX=int(math.floor(HomeMode/2-1))
NAX=NAX-1
Euu = np.zeros(NAX)
Delete(StartFile)
del StartFile
for file_name in rSTChannelal_vtu:
	print("***********************************************************")
	print("***********************************************************")
	print("Processing ",file_name) # Print out the filename we are currently processing
	VTUFILE = XMLUnstructuredGridReader(FileName=[file_name])
	VTUFILE.PointArrayStatus = ['u', 'v', 'w', 'p']
	SetActiveSource(VTUFILE)
	dataPlane = servermanager.Fetch(VTUFILE)
	nbCells = dataPlane.GetNumberOfCells()
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
	x=np.array(x)
	y=np.array(y)
	z=np.array(z)
	Nz=len(np.unique(x))
	Nx=len(np.unique(z))
	Ny=len(np.unique(y))
	log("DomainNx=",(Nx))
	log("DomainNy=",(Ny))
	log("DomainNz=",(Nz))
	kxS = np.zeros(NAX)
	cfmean = np.zeros(NAX)
	kxS = np.zeros(NAX)
	cfmean = np.zeros(NAX)
	print("-----------------------------------------------------------")
	SetActiveSource(VTUFILE)
	# create a new 'Slice'
	slice1 = Slice(SliceType="Plane")
	slice1.SliceType.Origin = [0.0, 0.835164835, 12.5]
	slice1.SliceType.Normal = [0.0, 1.0, 0.0]
	slice1.Triangulatetheslice = 0 #Don't triangulate the slice
	dataPlane = servermanager.Fetch(slice1)
	nbCells = dataPlane.GetNumberOfCells()
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
	# Put the velocity vector U into a numpy array		
	w = npvtk.vtk_to_numpy(dataPlane.GetPointData().GetArray('w'))	
	x=np.array(x)
	y=np.array(y)
	z=np.array(z)
	Nx=len(np.unique(z))
	Nz=len(w)/Nx #len(np.unique(x))
	wMean=np.mean(w)
	wFluc=w-wMean
	# Create a matrix to sort and reshape according to the XZ-plane
	raw=  x,y,z,wFluc # Span, WallNormal, Stream, VelFluc
	raw2=np.array(raw).transpose()
	raw3 = raw2[raw2[:,0].argsort()] 
	raw3 = raw3[raw3[:,2].argsort(kind='mergesort')]
	rawX=np.reshape(raw3[:,0], (Nx, Nz)).transpose()
	rawZ=np.reshape(raw3[:,2], (Nx, Nz)).transpose()
	rawW=np.reshape(raw3[:,3], (Nx, Nz)).transpose()
	# COMPUTE SPECTRA
	xStream=np.mean(rawZ, axis=0)
	n=len(xStream)
	na=n-1
	Nc=int(math.floor(na/2-1))
	per=xStream[n-1]-xStream[0]
	kx=2*np.pi*(np.arange(0,Nc))/per
	kxS[:]=kx[np.arange(1,Nc)]
	cf1 = np.zeros((Nc-1, Nz))
	for ii in range(0,len(rawW[:,1])):
		y=rawW[ii,:]
		cff=np.fft.fft(y,na)/na
		cf=np.multiply(cff,cff.conjugate())
		cf1[:,ii]=cf[np.arange(1,Nc)]
	cfmean[:]=np.pi*np.mean(cf1, axis=1)
	# Delete slice after calculation.
	Delete(slice1)
	del slice1
	Euu=Euu+cfmean
	# Delete the variable, so they don't have any stored values
	del cfmean, raw2, raw3, rawX, rawZ, rawW
	# Delete the file after processing.
	Delete(VTUFILE)
	del VTUFILE
print("***********************************************************")
Euu=Euu/float(len(rSTChannelal_vtu)) # Average of Euu over the time. For kx, not necessary.
OutFile=sys.argv[1].replace('.vtu','')
currentDir = os.getcwd()
OutFileDir = os.path.join(currentDir,newDir)
#kxS=np.array([1, 2, 3])
ksName=os.path.join(OutFileDir, OutFile+'_kx.csv')
EuuName=os.path.join(OutFileDir, OutFile+'_Euu.csv')
np.savetxt(ksName,kxS, delimiter=',')
np.savetxt(EuuName,Euu, delimiter=',')
#np.savetxt(OutFile+'Euu.csv',Euu, delimiter=',')
#np.savetxt(OutFile+'y.csv',yS, delimiter=',')
elapsed_time = (time.time()-start_time)/60.0
print("Total elapsed time for this code is=", "%.2f" %elapsed_time,"minutes")
print("-----------------------------------------------------------")
print("---------------------END OF THE PROGRAM--------------------")
print("-----------------------------------------------------------")