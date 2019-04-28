#!/usr/bin/env python
'''
ABOUT	: A Pyhton 2.7 script to read the data from VTU file from Nektar++ CFD code for porous DNS. 
		  After reading, we compute 2-D energy spectra from the fluctuations for a given timestep.
		  Currently, it saves wavenumber and energy in two different csv files, which can be used for contourmap.
NOTE	: RUN IT WITH PYTHON SHELL IN PARAVIEW OR pvpython/pvbatch from the terminal (only supports serial).
USAGE	: aprun -n 1 pvbatch PoroSpectra.py
DATE	: 20.03.2019
AUTHOR	: Sandeep Pandey (sandeep.pandey@ike.uni-stuttgart.de)
'''
from __future__ import print_function
from __future__ import division
import time
import sys
import numpy as np
import math
import os
import glob
import scipy.interpolate as scipyIn
# Import the modules from the paraview
from paraview.simple import *
from vtk.util import numpy_support as npvtk 
paraview.simple._DisableFirstRenderCameraReset()
	
# -----------------------------START: INPUT PARAMETERS ------------------------------        
# From where y start and where it ends and how many points you want
y_1=10
y_n=20
n_y=200
# Ho many points in the XZ-plane for interpolation
n_x = 600
n_z = 200
# How 'y' should be distributed from the lower wall? 'uniform or 'GP'?
y_type='uniform'
# Name of new directory where all spectrum data will be saved
newDir='spectrum_uu'
# Do you need extra information for e.g. number of cells/points etc?
ExtraInfo = False
# -----------------------------END: INPUT PARAMETERS -------------------------------  
# Create / check new directory
try:
    # Create target Directory
    os.mkdir(newDir)
    print("Directory spectrum  Created ") 
except OSError as e:
    if e.errno != os.errno.EEXIST:
		raise   
# ----------------------------- USER DEFINED FUNCTIONS -----------------------------
def log(msg,num):
	if ExtraInfo:
		print(msg, num)
#
def yArrange(y_1 = 10, y_n=20, n_y=20, y_type='uniform'):
    if y_type == "GP":
        r_y=(y_n/y_1)**(1/(n_y-1))
        yUnique = [y_1 * r_y ** (n - 1) for n in range(1, n_y + 1)]
    else:
        yUnique=np.linspace(y_1, y_n, num=n_y)         
    return yUnique;
#
def generate_grid(x_min = 0, x_max=10, n_x=200,z_min = 0, z_max=25, n_z=200):
    x=np.linspace(x_min, x_max, num=n_x)
    z=np.linspace(z_min, z_max, num=n_z)
    return np.meshgrid(x,z);      

# -------------------------MAIN PROGRAM -------------------------  
print("-----------------------------------------------------------")
print("-------------------START OF THE PROGRAM--------------------")
print("-----------------------------------------------------------")
start_time = time.time()
# List out all files with RSTChannel-al_*.vtu. 
rSTChannelal_vtu =  sys.argv[1:] # glob.glob('./clipped1.vtu') #
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

# Custom functions: yArrange will give us values of y in GP.
yUnique=yArrange(y_1 , y_n, n_y)
x_grid,z_grid=generate_grid(x_min = np.amin(x), x_max = np.amax(x), n_x=n_x, 
							z_min = np.amin(z), z_max = np.amax(z), n_z=n_z )
					
NAX=int(math.floor(len(x_grid[1,:])/2))
NAX=NAX
Eww = np.zeros((NAX, len(yUnique)))
Delete(StartFile)
del StartFile

for file_name in rSTChannelal_vtu:
	print("***********************************************************")
	# Print out the filename we are currently processing
	print("Processing ",file_name) 
	VTUFILE = XMLUnstructuredGridReader(FileName=[file_name])
	VTUFILE.PointArrayStatus = ['u'] #['u', 'v', 'w', 'p']
	SetActiveSource(VTUFILE)
	dataPlane = servermanager.Fetch(VTUFILE)
	nbCells = dataPlane.GetNumberOfCells()
	nbPoints = dataPlane.GetNumberOfPoints()
	log("Total number of cells in the domain", nbCells)
	log("Total number of points in the domain", nbPoints)
	kzS = np.zeros((NAX, len(yUnique)))
	cfmean = np.zeros((NAX, len(yUnique)))
	for jj in range(0,len(yUnique)):
		print("-----------------------------------------------------------")
		print("Calculating for slice=",jj ,"Remaining slice=",len(yUnique)-jj)
		print("-----------------------------------------------------------")
		SetActiveSource(VTUFILE)
		# Create a new 'Slice'
		slice1 = Slice(SliceType="Plane")
		slice1.SliceType.Origin = [0.0, yUnique[jj], 12.5]
		slice1.SliceType.Normal = [0.0, 1.0, 0.0]
		# Don't triangulate the slice
		slice1.Triangulatetheslice = 0 
		dataPlane = servermanager.Fetch(slice1)
		nbCells = dataPlane.GetNumberOfCells()
		nbPoints = dataPlane.GetNumberOfPoints()
		log("nbPoints=",nbPoints)
		# Put the points coordinates in x, y and z arrays
		x_slice=[]
		y_slice=[]
		z_slice=[]
		for i in range(nbPoints):
			coord = dataPlane.GetPoint(i)
			xx_sl, yy_sl, zz_sl = coord[:3]
			x_slice.append(xx_sl)
			y_slice.append(yy_sl)
			z_slice.append(zz_sl)
		x_slice=np.array(x_slice)
		y_slice=np.array(y_slice)
		z_slice=np.array(z_slice)    
		# Put the velocity vector U into a numpy array	
		u = npvtk.vtk_to_numpy(dataPlane.GetPointData().GetArray('u'))
		Velocity = scipyIn.griddata(np.transpose(np.array([x_slice,z_slice])), u, (x_grid, z_grid), method='nearest')
		Nx=len(np.unique(z_grid[:,1]))
		Nz=len(np.unique(x_grid[1,:]))
		uMean=np.mean(Velocity)
		uFluc=Velocity-uMean
		# Compute spectra with FFT
		xStream=np.mean(x_grid, axis=0)
		Lz=xStream[Nx-1]-xStream[0]
		n=len(xStream)
		na=n-1
		Nc=int(math.floor(n/2))
		per=xStream[n-1]-xStream[0]
		kz=2*np.pi*(np.arange(0,Nc))/Lz
		kzS[:,jj]=kz
		cf1 = np.zeros((Nc, Nz))
		for ii in range(0,len(uFluc[:,1])):
			Uzfluc=uFluc[ii,:]
			cff0=(np.fft.fft(Uzfluc,Nx))/Nx
			cff=cff0
			cf=np.multiply(cff,cff.conjugate())
			cf1[:,ii]=cf[np.arange(0,Nc)]
		cfmean[:,jj]=np.pi*np.mean(cf1, axis=1)
		# Delete slice after calculation.
		Delete(slice1)
		del slice1
	Eww=Eww+cfmean
	# Delete the variable, so they don't have any stored values
	del u, uFluc, uMean, Velocity, cfmean
	# Delete the file after processing.
	Delete(VTUFILE)
	del VTUFILE
print("***********************************************************")
# Average of E over the time. For kx, not necessary.
Eww=Eww/float(len(rSTChannelal_vtu)) 
# Few thinsg to save the files in desired location with desired name
OutFile=sys.argv[1].replace('.vtu','')
currentDir = os.getcwd()
OutFileDir = os.path.join(currentDir,newDir)
kName=os.path.join(OutFileDir, OutFile+'_k.csv')
EName=os.path.join(OutFileDir, OutFile+'_E.csv')
yName=os.path.join(OutFileDir, OutFile+'_y.csv')
# Save them as csv
np.savetxt(kName,kzS, delimiter=',')
np.savetxt(EName,Eww, delimiter=',')
np.savetxt(yName,yUnique, delimiter=',')

elapsed_time = (time.time()-start_time)/60.0
print("Total elapsed time for this code is=", "%.2f" %elapsed_time,"minutes")
print("-----------------------------------------------------------")
print("---------------------END OF THE PROGRAM--------------------")
print("-----------------------------------------------------------")