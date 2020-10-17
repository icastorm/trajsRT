###############################################################################
# trajUtils_v.py
# 
# Last update: 6/21/2020 by Isaac Arseneau
#
# Description: Supports the real-time trajDriver
#
# Similar to trajUtils.py but with horizontal vectorization
###############################################################################


### Imports ###
# Imports
from netCDF4 import Dataset
import numpy as np
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim, ll_to_xy, xy_to_ll, xy_to_ll_proj,destagger)
import metpy
from metpy import interpolate as interp
from metpy.units import units
import metpy.calc as mcalc
import math
from datetime import datetime, timedelta
from os import path



# Take a list of x and y values, return lists for 2d analysis
# For 1D inputs of sizes M,N, returns 2 lists of size M*N,M*N
def make_2d(ogX, ogY):
    newX = ogX*len(ogY)
    newY = []
    for i in range(len(ogY)):
        for j in range(len(ogX)):
            newY.append(ogY[i])
    return newX, newY


# Read in the current settings for the domain
# Currently returns latMin,latMax,lonMin,lonMax,roaMin,roaMax
def getSettings(sFile = '/home/iarsenea/trajs/settings.txt'):
	settings = np.loadtxt(sFile,usecols=1,dtype=str)
	return [float(i) for i in settings[0:6]]


# Given the mm dd hh info of the run start time, output the wrf file
# that corresponds to 12 hours later 
def getWrfFile(mm,dd,hh,member,sFile = '/home/iarsenea/trajs/settings.txt'):

	initDate="2020"+mm+dd+hh
	settings = np.loadtxt(sFile,usecols=1,dtype=str)
	basePath = settings[6]
	fPath = basePath+"mem"+str(member)+"/"
	date = datetime.strptime(initDate,"%Y%m%d%H")+timedelta(hours=12.)

	wrfFile = fPath+"wrfout_d01_"+date.strftime("%Y-%m-%d_%H:00:00")

	return wrfFile


# Given the mm dd hh info of the run start time, return the path to
# the wrf_done file, if it were to exist
def getDoneFile(mm,dd,hh,member,sFile = '/home/iarsenea/trajs/settings.txt'):

	initDate="2020"+mm+dd+hh
	settings = np.loadtxt(sFile,usecols=1,dtype=str)
	basePath = settings[6]
	fPath = basePath+"mem"+str(member)+"/"

	doneFile = fPath+"wrf_done"

	return doneFile


# Given the key value of the run, return the available members
def getAvailMems(key):

	key = str(key)

	# First number is id, second is total number of cores being used
	# The first character must always be a number

	# Using 2
	if key == "12":
		return [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,33,34,35,36,37]

	if key == "22":
		return [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,38,39,40,41,42]

	# Some tests
	if key == "1t":
		return [500,1]

	if key == "2t":
		return [5,6,7,8]

	print("No valid key detected.")
	return []

# Read the names of files currently in the linkPath and infer the
# mm dd hh values from those names
def getRecentDate(sFile = '/home/iarsenea/trajs/settings.txt'):

	linkPath = "/home/iarsenea/trajs/trajFiles_links/"

	link = "trajLink_0"

	ogPath = path.realpath(linkPath+link)

	print(ogPath)


# Given the mm dd hh info of the run start time, output the dest traj file
# that corresponds to 12 hours later 
def getTrajFile(mm,dd,hh,member,roa,sFile = '/home/iarsenea/trajs/settings.txt'):

	initDate="2020"+mm+dd+hh
	settings = np.loadtxt(sFile,usecols=1,dtype=str)
	trajPath = settings[7]
	date = datetime.strptime(initDate,"%Y%m%d%H")+timedelta(hours=12.)

	trajFile = trajPath+"traj_"+str(member)+"_"+str(roa)+"_"+date.strftime("%Y-%m-%d_%H:00:00")+"_v"

	return trajFile


# Given the starting value and the number of files to use,
# returns the names of the links to be used for the files
# that will be created on this core
def getLinks(start,nFiles,sFile = '/home/iarsenea/trajs/settings.txt'):

	# Read in the path for the link directory
	settings = np.loadtxt(sFile,usecols=1,dtype=str)
	linkPath = settings[8]

	links = []
	for i in range(nFiles):
		links.append(linkPath+"trajLink_"+str(i+start))

	return links


# Given a wrf file path, open the file, get what is needed
# for the trajs out of it, and then close it.
def getWrfData(wrfPath,dx,dy,dz,levs,debug=True):

	# Open the file
	perDat = Dataset(wrfPath)
	refDat = Dataset("/home/iarsenea/trajs/wrfoutREFd01")

	#print(data.variables)

	# Pull in the values from the base state that we will add to the perturbations
	ref_h = getvar(refDat, "height", units="m", timeidx=0)
	thght = np.asarray(refDat.variables["HGT"])[0] * units.meter
	lats, lons = latlon_coords(ref_h)

	# Pull in the values from the perturbation 
	ph = np.asarray(destagger(perDat.variables["PH"][0], 0)) * units('m^2/s^2')
	phb = np.asarray(destagger(refDat.variables["PHB"][0], 0)) * units('m^2/s^2')
	ua = np.asarray(destagger(perDat.variables["U"][0], 2)) * units('m/s')
	va = np.asarray(destagger(perDat.variables["V"][0], 1)) * units('m/s')

	# Calculate geopotential
	print("Converting from perturbation height to height AGL...")
	geo = ph + phb

	# Convert to height
	hght = mcalc.geopotential_to_height(geo)

	# Convert to height_agl
	h = np.zeros_like(hght)
	for i in range(hght.shape[0]):
	    h[i] = hght[i] - thght
	print("Done.\n")

	# Get the x and y values of the lat and lon coordinates
	x, y = ll_to_xy(refDat, lats, lons)
	x = np.arange(0,np.max(x.data)+1)*dx
	y = np.arange(0,np.max(y.data)+1)*dy

	# Interpolate the winds speeds and temps to the heights specified in levs
	if debug: print("\nInterpolating wind values to every " + str(dz.m) + " meters... \n")
	ua_m = metpy.interpolate.interpolate_1d(levs, h, ua)
	va_m = metpy.interpolate.interpolate_1d(levs, h, va)


	perDat.close()
	refDat.close()


	return h,x,y,ua_m,va_m


# Create the proper file and run the organic trajectory
def calcTrajs_org(startXs,startYs,x,y,ua_m,va_m,levs,dx,dy,dz,roa,destFile):

	# Create an ncfile to write the trajectories to
	outfile = Dataset(destFile, "w", format="NETCDF4")
	dim_level = outfile.createDimension("level", len(levs)-1)
	dim_traj = outfile.createDimension("trajectory", len(startXs))
	var_level = outfile.createVariable("level", "f4", ("trajectory","level"))
	var_x = outfile.createVariable("x", "f4", ("trajectory","level"))
	var_y = outfile.createVariable("y", "f4", ("trajectory","level"))

	# Create a bunch more variables, for later sensitivity data storage
	# Response var 1
	var_uESA_1 = outfile.createVariable("uESA_1", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_1 = outfile.createVariable("vESA_1", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_1 = outfile.createVariable("tESA_1", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_1 = outfile.createVariable("qESA_1", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 2
	var_uESA_2 = outfile.createVariable("uESA_2", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_2 = outfile.createVariable("vESA_2", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_2 = outfile.createVariable("tESA_2", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_2 = outfile.createVariable("qESA_2", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 3
	var_uESA_3 = outfile.createVariable("uESA_3", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_3 = outfile.createVariable("vESA_3", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_3 = outfile.createVariable("tESA_3", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_3 = outfile.createVariable("qESA_3", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 4
	var_uESA_4 = outfile.createVariable("uESA_4", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_4 = outfile.createVariable("vESA_4", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_4 = outfile.createVariable("tESA_4", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_4 = outfile.createVariable("qESA_4", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 5
	var_uESA_5 = outfile.createVariable("uESA_5", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_5 = outfile.createVariable("vESA_5", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_5 = outfile.createVariable("tESA_5", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_5 = outfile.createVariable("qESA_5", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 6
	var_uESA_6 = outfile.createVariable("uESA_6", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_6 = outfile.createVariable("vESA_6", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_6 = outfile.createVariable("tESA_6", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_6 = outfile.createVariable("qESA_6", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 7
	var_uESA_7 = outfile.createVariable("uESA_7", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_7 = outfile.createVariable("vESA_7", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_7 = outfile.createVariable("tESA_7", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_7 = outfile.createVariable("qESA_7", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 8
	var_uESA_8 = outfile.createVariable("uESA_8", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_8 = outfile.createVariable("vESA_8", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_8 = outfile.createVariable("tESA_8", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_8 = outfile.createVariable("qESA_8", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 9
	var_uESA_9 = outfile.createVariable("uESA_9", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_9 = outfile.createVariable("vESA_9", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_9 = outfile.createVariable("tESA_9", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_9 = outfile.createVariable("qESA_9", "f4", ("trajectory","level"), fill_value = -999.0)

	# Response var 10
	var_uESA_10 = outfile.createVariable("uESA_10", "f4", ("trajectory","level"), fill_value = -999.0)
	var_vESA_10 = outfile.createVariable("vESA_10", "f4", ("trajectory","level"), fill_value = -999.0)
	var_tESA_10 = outfile.createVariable("tESA_10", "f4", ("trajectory","level"), fill_value = -999.0)
	var_qESA_10 = outfile.createVariable("qESA_10", "f4", ("trajectory","level"), fill_value = -999.0)

	# Calculate the levels for every single starting point at once
	ogX = startXs * units.meter
	ogY = startYs * units.meter

	# Calculate the trajectory
	xs,ys = trajectory_w_val_v(ogX,ogY,x,y,ua_m,va_m,levs,dx,dy,dz,roa)

	# Write the trajectory to the outfile
	for i in range(len(startXs)):
		var_level[i,:]=levs[:-1]
	var_x[:]=xs
	var_y[:]=ys


	outfile.close()


# Calculate a trajectory given a starting point and the relevant data field
# Returns xs,ys,vals
def trajectory_w_val_v(ogX,ogY,x,y,u,v,levs,dx,dy,dz,ROA,output=False):

	# Create arrays of 0s to store the x and y values in
	xs = np.zeros((len(ogX),len(levs)-1))
	ys = np.zeros_like(xs)

	#print(xs.shape)

	# The first level is just the starting positions
	xs[:,0] = ogX
	ys[:,0] = ogY

	thisX,thisY = ogX,ogY

	#print(np.column_stack((thisX,thisY)).shape)
	#print(u[0].flatten().shape)
	#print(x.shape)
	#print(np.amax(x),np.amax(thisX))

	x,y = np.meshgrid(x,y)

	# Use x and y to generate a coordinate array
	og_coords = np.vstack((x.flatten(),y.flatten())).T
	#print(og_coords.shape)

	print("Calculating balloon trajectories...\n")

	# Iterate over each vertical level, saving the x/y values as you go
	for i in range(len(levs)-2):

		# Step up a level
		thisX,thisY = step(u,v,thisX,thisY,i,og_coords,levs,dx,dy,dz,ROA,output)

		# Save the new values
		xs[:,i+1] = thisX
		ys[:,i+1] = thisY

	print("Done.\n")

	return xs,ys


# Step the trajectory up one level
def step(uVals,vVals,thisX,thisY,iZ,og_coords,levs,dx,dy,dz,roa,output=True):

	# First find the layer average wind by finding the value of the starting point 
	# and the value of the point directly above

	# u/v values of the current level and the level above
	u0 = interp.interpolate_to_points(og_coords/dx,uVals[iZ].flatten(),np.column_stack((thisX,thisY))/dx)
	v0 = interp.interpolate_to_points(og_coords/dy,vVals[iZ].flatten(),np.column_stack((thisX,thisY))/dy)
	u1 = interp.interpolate_to_points(og_coords/dx,uVals[iZ+1].flatten(),np.column_stack((thisX,thisY))/dx)
	v1 = interp.interpolate_to_points(og_coords/dy,vVals[iZ+1].flatten(),np.column_stack((thisX,thisY))/dy)


	u = ((u0+u1)/2)*units.meter/units.second
	v = ((v0+v1)/2)*units.meter/units.second

	# Second, use that value to figure out new x and y locations
	# outX = ((dz/roa)*u)+thisX

	#print(dz)
	#print(roa)
	#print(u[0])
	#print(thisX[0])
	dz_roa = dz/roa

	outX = ((dz_roa)*u)+thisX
	print(iZ)
	#print(outX.shape)
	#print(((dz_roa)*u)[0])
	#print()
	outY = ((dz_roa)*v)+thisY

	if output:
		print("Layer is from " + str(levs[iZ]) + " to " + str(levs[iZ+1]) + ".")
		print("Average u component = " + str(u[0]))
		print("Average v component = " + str(v[0]))
		print("Average wind speed = " + str(((u[0]**2)+(v[0]**2))**.5))
		print("ROA = " + str(roa))
		print("dz = " + str(dz))
		print("Amount of time spent in layer: " + str((dz_roa)))
		print("DX = " + str(((dz_roa)*u)[0]))
		print("Return X is " + str((outX[0])))
		print()

	return outX, outY

# Return the given ROA (x10) scaled from 0 to 1
def scaledROA(roa,minROA=30,maxROA=70):
	return (roa-minROA)/(maxROA-minROA)


# Return the curve used in a numpy traj1d object
def getCurve():
	return np.poly1d([0.66702064,0.79965154,0.83248499,0.29680275])


# Use the curve to estimate the trajectories for the intermediate ROA values.
# The output files will only have levs and x and y values
def calcTrajs_curve(dataS,dataF,newFiles,newROAs):

	# Pull in the information needed for each roa
	newROAs_scaled = [scaledROA(roa) for roa in newROAs]
	#print(newROAs_scaled)
	dataLevs = dataS["level"][:]
	curve = getCurve()

	xS = dataS["x"][:]
	yS = dataS["y"][:]
	xF = dataF["x"][:]
	yF = dataF["y"][:]
	numLevs = xS.shape[1]
	numTrajs = xS.shape[0]

	# Calculate the dx0 and dy0 values
	dx0 = xS - xF
	dy0 = yS - yF

	# Create the ratio array
	ratio = []
	for i in range(len(newROAs)):
		ratio.append(curve(.5-scaledROA(newROAs[i])))
	ratio=np.asarray(ratio)

	# Create each file and perform the calculation on it
	for i in range(len(newFiles)):

		destFile = newFiles[i]
		roa = newROAs[i]

		print("Calculating for " + destFile + "...")

		# Create an ncfile to write the trajectories to
		outfile = Dataset(destFile, "w", format="NETCDF4")
		dim_level = outfile.createDimension("level", numLevs)
		dim_traj = outfile.createDimension("trajectory", numTrajs)
		var_level = outfile.createVariable("level", "f4", ("trajectory","level"))
		var_x = outfile.createVariable("x", "f4", ("trajectory","level"))
		var_y = outfile.createVariable("y", "f4", ("trajectory","level"))

		# Create a bunch more variables, for later sensitivity data storage
		# Response var 1
		var_uESA_1 = outfile.createVariable("uESA_1", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_1 = outfile.createVariable("vESA_1", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_1 = outfile.createVariable("tESA_1", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_1 = outfile.createVariable("qESA_1", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 2
		var_uESA_2 = outfile.createVariable("uESA_2", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_2 = outfile.createVariable("vESA_2", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_2 = outfile.createVariable("tESA_2", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_2 = outfile.createVariable("qESA_2", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 3
		var_uESA_3 = outfile.createVariable("uESA_3", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_3 = outfile.createVariable("vESA_3", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_3 = outfile.createVariable("tESA_3", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_3 = outfile.createVariable("qESA_3", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 4
		var_uESA_4 = outfile.createVariable("uESA_4", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_4 = outfile.createVariable("vESA_4", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_4 = outfile.createVariable("tESA_4", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_4 = outfile.createVariable("qESA_4", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 5
		var_uESA_5 = outfile.createVariable("uESA_5", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_5 = outfile.createVariable("vESA_5", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_5 = outfile.createVariable("tESA_5", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_5 = outfile.createVariable("qESA_5", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 6
		var_uESA_6 = outfile.createVariable("uESA_6", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_6 = outfile.createVariable("vESA_6", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_6 = outfile.createVariable("tESA_6", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_6 = outfile.createVariable("qESA_6", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 7
		var_uESA_7 = outfile.createVariable("uESA_7", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_7 = outfile.createVariable("vESA_7", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_7 = outfile.createVariable("tESA_7", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_7 = outfile.createVariable("qESA_7", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 8
		var_uESA_8 = outfile.createVariable("uESA_8", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_8 = outfile.createVariable("vESA_8", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_8 = outfile.createVariable("tESA_8", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_8 = outfile.createVariable("qESA_8", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 9
		var_uESA_9 = outfile.createVariable("uESA_9", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_9 = outfile.createVariable("vESA_9", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_9 = outfile.createVariable("tESA_9", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_9 = outfile.createVariable("qESA_9", "f4", ("trajectory","level"), fill_value = -999.0)

		# Response var 10
		var_uESA_10 = outfile.createVariable("uESA_10", "f4", ("trajectory","level"), fill_value = -999.0)
		var_vESA_10 = outfile.createVariable("vESA_10", "f4", ("trajectory","level"), fill_value = -999.0)
		var_tESA_10 = outfile.createVariable("tESA_10", "f4", ("trajectory","level"), fill_value = -999.0)
		var_qESA_10 = outfile.createVariable("qESA_10", "f4", ("trajectory","level"), fill_value = -999.0)


		# Calculate the trajectory based on the curve
		xs = np.empty((numTrajs,numLevs))
		ys = np.empty((numTrajs,numLevs))

		# First level is the same
		xs[:,0],ys[:,0] = xS[:,0],yS[:,0]

		# Now do the rest
		xs[:,1:], ys[:,1:] = ratio[i]*dx0[:,1:] + xF[:,1:], ratio[i]*dy0[:,1:] + yF[:,1:]

		# Write the results to file
		var_x[:] = xs
		var_y[:] = ys
		var_level[:] = dataLevs

		outfile.close()



