###############################################################################
# trajDriver_v.py
# 
# Last update: 6/12/2020 by Isaac Arseneau
#
# Description: Activated when the realtime driver detects that the
# ensemble has finished running out to +12 hours (valid at 00z) for
# the 12z daily run. Sets up and manages trajectory runs, both
# organic and curve-generated, while tracking progress and how the
# trajs were initialized in an output file. Additionally, a file
# containing a list of the file paths created is created to be read
# in by the code interpolating target values to the trajs.
#
# Traj files are generated for each member and rate of ascent (ROA),
# and these files contain the x/y (m) and lat/lon locations of each
# trajectory with height in the shape (nTrajs,nLevels).
# Additionally, each file contains 10 empty variables in the same
# shape as above that will store the interpolated target values.
#
# All variables are single-precision floats.
#
# Takes in the multiple index of the symlinks, the mm dd hh of the run start
# time, and the members to use
#
# Similar to trajDriver.py but with horizontal vectorization
###############################################################################


### Imports ###

from netCDF4 import Dataset
import numpy as np
import trajUtils_v as utils
from metpy.units import units
import sys
from os import (symlink, unlink, path)
args = sys.argv

### Setup Parameters ###

# These may change
dz = 50*units.meter # Distance between vertical levels in meters
droa = 10 # m/10s
buff = 10 # Number of wind points to surround the domain with
memGoal = 16 # Total number of members that should be used

# These will likely not change
latMin,latMax,lonMin,lonMax,roaMin,roaMax = utils.getSettings() # Read in the necessary settings
dx,dy = 12000*units.meter,12000*units.meter # Horizontal spacing in (m)
ogX,ogY = list(np.arange(lonMin*dx.m,lonMax*dx.m,dx.m)),list(np.arange(latMin*dy.m,latMax*dy.m,dy.m))
ogX,ogY = utils.make_2d(ogX,ogY)
levs = np.arange(50, 10000+(dz.m*2), dz.m) * units.meter
roas = np.arange(roaMin,roaMax+droa,droa)
organicROAs = (roaMin,roaMax)

print()
print(latMin,latMax,lonMin,lonMax)
print(roas)
print()

# Read in the args
startSym = int(str(args[1])[0]) - 1 # should be 0 (a-1) for the first job, 1 (a-1) for the second, and so on
memKey = str(args[1]) # ab where a is the core number and b is the total number of cores being used
mm = args[2]
dd = args[3]
hh = args[4]

# Get the members to be used for this run
members = utils.getAvailMems(memKey)

# Get the names of the links to be used
symNames = utils.getLinks(startSym*len(roas)*memGoal,len(roas)*memGoal)

# Remove the existing links ahead of time
for link in symNames:
	try:
		unlink(link)
	except:
		print("No link to remove.")

print("Link names:",symNames[0],"...",symNames[-1])


### Run the trajs ###

# For each member, generate the proper files to look for to run the organic trajs
# and do the wrf file operations that need to be done

# Also, track the file # being generated
iLink = 0

# Keep track of the members used
usedMems = [] 

keepGoing = True
while keepGoing:

	# Pop the fist unused/uncrashed member
	mem = members.pop(0)

	# Same wrf file will be used for both organic trajs
	wrfFile = utils.getWrfFile(mm,dd,hh,mem)
	doneFile = utils.getDoneFile(mm,dd,hh,mem)

	# Check to ensure that the wrf file exists.
	# If not, start going through the other members

	# NOTE: If there are many members who crash/are very behind the others
	# and do not put out a wrf_done file, then the script may hang here as 
	# it waits for those files to come in.
	while not path.exists(wrfFile):

		print()
		print("WRF File",wrfFile,"does not exist. Trying extra members...")

		# If the run simply isn't done yet, add this member to the end of the list.

		# If the wrf_done file is there and the output file is not, that means
		# the member has likely crashed and should be skipped.

		# This way, if the wrf member simply isn't to +12 yet and it is needed later,
		# the list will come back around after checking to make sure none of the
		# other available members are already done.

		# Essentially, the fastest "goal#" members are used.

		if not path.exists(doneFile):

			print("This member may not be finished. Adding to end of list...")
			members.append(mem)

		# Pop the next available member
		mem = members.pop(0)

		# Generate the new file locations for the next member to be used
		wrfFile = utils.getWrfFile(mm,dd,hh,mem)
		doneFile = utils.getDoneFile(mm,dd,hh,mem)

	print()
	print("--------------------------------------")
	print("Calculating Trajs using Member: ",mem)
	print("--------------------------------------")
	print()

	# Read in the necessary information from the wrf file
	h,x,y,ua_m,va_m = utils.getWrfData(wrfFile,dx,dy,dz,levs)
	print()
	print(x.shape)

	# Subset the x/y and u/v horizontal locations to speed up interpolation
	x = x[int(lonMin)-buff:int(lonMax)+buff]
	y = y[int(latMin)-buff:int(latMax)+buff]
	ua_m = ua_m[:,int(latMin)-buff:int(latMax)+buff,int(lonMin)-buff:int(lonMax)+buff]
	va_m = va_m[:,int(latMin)-buff:int(latMax)+buff,int(lonMin)-buff:int(lonMax)+buff]
	print(ua_m.shape)
	print(x.shape)

	print()
	print('Calculating organic trajs...')
	print()

	# Run both of the organic trajs needed
	orgFiles = []
	for j in range(len(organicROAs)):

		# Generate the destination file name for the trajs
		trajFile = utils.getTrajFile(mm,dd,hh,mem,organicROAs[j])
		orgFiles.append(trajFile)
		print(wrfFile,trajFile,"\n")

		# Run the trajectory calculation
		utils.calcTrajs_org(ogX,ogY,x,y,ua_m,va_m,levs,dx,dy,dz,(organicROAs[j]/10)*units.meter/units.second,trajFile)

		# Link the new file
		symlink(trajFile, symNames[iLink])
		iLink += 1

	print()
	print('Calculating curve trajs...')
	print()

	# Get the intermediate file names
	interFiles = []
	interROAs = []
	for j in range(len(roas)-2):

		thisROA = roas[j+1]
		interROAs.append(thisROA)

		# Get the proper name of the file to be made
		interFiles.append(utils.getTrajFile(mm,dd,hh,mem,thisROA))

	# Read in the two organic files, slow then fast
	sFile = Dataset(orgFiles[0])
	fFile = Dataset(orgFiles[1])

	# Create and populate each of the new files
	utils.calcTrajs_curve(sFile,fFile,interFiles,interROAs)

	# Link the proper files
	for trajFile in interFiles:

		# Link the new file
		symlink(trajFile, symNames[iLink])
		iLink += 1

	print("Done.")
	print()

	# Close the files used
	sFile.close()
	fFile.close()

	# Record that this traj calc went through
	usedMems.append(mem)

	# Stop if the total number of members used is the proper number
	if len(usedMems) >= memGoal:
		keepGoing = False

print()
print("-----------------------------------------------------------")
print()
print("Done with all members.")
print("Members used:",usedMems)
print()
