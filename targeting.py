### Imports ###

from netCDF4 import Dataset
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
from wrf import xy_to_ll
from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.signal import argrelmax
from functools import partial
import sys

args = sys.argv

date = str(args[1])
print(date)


### Definitions ###

# Returns the list of linked files
def getLinks(nLinks = 160,linkPath = "/home/iarsenea/trajs/trajFiles_links/"):
    links = []
    for i in range(nLinks):
        links.append(linkPath+"trajLink_"+str(i))
    return links


def getLinks_test(nLinks = 160,baseDir = "/home/iarsenea/trajs/trajFiles_links/"):
    links = os.listdir(baseDir)#[::2]
    for i in range(len(links)):
        links[i] = baseDir+links[i]
    return links


def getLocs_test():
    links = getLinks_test()
    data = Dataset(links[0])
    xs = data['x'][:,0]
    ys = data['y'][:,0]
    print(np.amin(np.ma.masked_greater(data['qESA_2'],0)))
    print(np.amax(np.ma.masked_greater(data['qESA_2'],0)))
    
    return xs,ys


def convertLatLon(xs,ys):
    # Convert the launch points to lat/lon
    refDat = Dataset("/home/iarsenea/trajs/wrfoutREFd01")
    lats,lons = xy_to_ll(refDat,np.asarray(xs),np.asarray(ys))
    return(lats,lons)


def centers_to_edges_1d(x):
    dx = np.absolute(x[0]-x[1])/2
    return np.arange(np.min(x)-dx,np.max(x)+dx+(dx/2),2*dx)


# Should only require one read through of the files
# Final output for each response function is four fields: mean maximum value, % maximum values above 90th percentile,
# mean maximum % variance, and % of maximum % variance above 20%
def doCalcs(date,nResps = 10):
    
    # Get the links
    #links = getLinks()
    links = getLinks_test()
    
    # Read in the first dataset to get the proper shape
    data = Dataset(links[0])
    nTrajs = data['x'].shape[0]
    nLevs = data['x'].shape[1]
    
    data.close()
    
    # Create output arrays to store data
    mmVal = np.zeros((nResps,nTrajs)) # mean maximum value
    pmVal = np.zeros((nResps,nTrajs)) # % maximum values above 90th percentile
    mmVar = np.zeros((nResps,nTrajs)) # mean maximum % variance
    pmVar = np.zeros((nResps,nTrajs)) # % of maximum % variance above 20%
    
    print()
    print(mmVal.shape)
    print()
    
    # Get the proper percentile values to compare against
    percVals = getPercs()
    
    # Get the total variance values to compare against
    varVals = getRVar(date)
    
    # Read in every file and do the necessary calculations
    for i in range(len(links)):
        link = links[i]
        data = Dataset(link)
        print(i,link)
        
        # Make an empty array so that the values can be stored
        vals = np.empty((nResps,nTrajs,nLevs),dtype=float)
        print(vals.shape)
        
        # Read the values in and sum over the base variables
        for j in range(nResps):
            
            vals[j] = data["uESA_"+str(j+1)][:]+data["vESA_"+str(j+1)][:]+data["tESA_"+str(j+1)][:]+data["qESA_"+str(j+1)][:]
            
        data.close()
            
        # Take max over the last dimension (levs)
        # Left with (nResps, nTrajs)
        m = np.amax(vals,axis=2)
        
        # Create copies of m to use for calcs
        per90 = np.copy(m)
        per20 = np.copy(m)
        
        mVar = np.copy(m)
        
        print(m.shape)
        
        # Change anything less than the percentile values to 0, more than to 1s
        for j in range(m.shape[0]):
            
            # Percent of max vals > 50th percentile
            per90[j][m[j]<percVals[j+1,0]]=0
            per90[j][m[j]>=percVals[j+1,0]]=1
            
            # Percent of total variance
            mVar[j] = (mVar[j]/(varVals[j]))*100
            print("Response Function:",j+1)
            print("Variance:",varVals[j])
            print("Maximum value:",np.amax(m[j]))
            print("Maximum % variance:",np.amax(mVar[j]))
            print()
            
            # Percent of max vals > 15% of total variance (roughly the 75th percentile)
            per20[j][mVar[j]<15]=0
            per20[j][mVar[j]>=15]=1
            
        
        # Add to the running sums
        mmVal += m
        pmVal += per90
        mmVar += mVar
        pmVar += per20
        
    # Once all files have been read, divide the sums to get the averages
    nLinks = len(links)
    mmVal /= nLinks
    pmVal /= nLinks
    mmVar /= nLinks
    pmVar /= nLinks
        
    return mmVal,pmVal*100,mmVar,pmVar*100


# Given the four main calculate outputs, return as many as 3 locations for each response function to target
# Write the values to file
def chooseTargs(mmVal,pmVal,mmVar,pmVar,date):
    
    lFile = "/home/iarsenea/trajs/targeting/"+date+"/log_"+date+".txt"
    
    clearLog(logFile = lFile)
    log = partial(writeLog,logFile = lFile)
    print(np.amax(pmVar))
    print(np.amax(pmVar[:,478:]))
    
    # Get the location of the launch points
    xs,ys = getLocs_test()
    xs = xs/12000
    ys = ys/12000
    
    lats,lons=convertLatLon(xs,ys)
    
    #lats,lons = convertLatLon(edge_xs,edge_ys)
    
    N = int(np.amax(ys)-np.amin(ys))+1 # yMax-yMin
    M = int(np.amax(xs)-np.amin(xs))+1 # xMax-xMin
    
    lats = np.asarray(lats).reshape(N,M)
    lons = np.asarray(lons).reshape(N,M)
    
    filtmmVal = np.copy(mmVal)
    print(filtmmVal.shape)
    
    targLats = []
    targLons = []
    
    for i in range(10):
        #print()
        #print("Response Function",i+1)
        #print("================================================================================")
        
        log("\nResponse Function "+str(i+1)+"\n")
        log("================================================================================\n")
        
        # Filter points by those w/pmVar >= 75%
        filter1 = np.where(pmVar[i]<75)[:]
        filtmmVal[i][filter1] = 0
    
        # Apply second filter for points by those with pmVal >= 75%
        filter2  = np.where(pmVal[i]<75)[:]
        filtmmVal[i][filter2] = 0
        
        # Convert the values to 2d
        data = filtmmVal[i].reshape(N,M)
        
        # Take the top 3 local mmVal values and return their locations plus some info
        maxima = np.asarray(argrelmax(data,order=20))
        
        # Get the values of the local maxima
        vals = []
        newLats = []
        newLons = []
        for j in range(maxima.shape[1]):
            vals.append(data[maxima[0,j],maxima[1,j]])
            newLats.append(lats[maxima[0,j],maxima[1,j]])
            newLons.append(lons[maxima[0,j],maxima[1,j]])
        
        # Get the top three
        top_3_idx = np.flip(np.argsort(vals)[-15:])
        top_3_values = [vals[i] for i in top_3_idx]
        
        newLats2 = []
        newLons2 = []
        # Get the lats and lons for those values:
        if list(top_3_idx):
            for j in range(len(top_3_idx)):
            #for j in range(len(vals)):
                
                tLat = newLats[top_3_idx[j]]
                tLon = newLons[top_3_idx[j]]

                #tLat = newLats[j]
                #tLon = newLons[j]

                newLats2.append(tLat)
                newLons2.append(tLon)
                
                #print(str(j+1)+":")
                #print()
                #print(top_3_values[j])
                #print(str(tLat)+", "+str(tLon))
                #print("--------------------------------------------------------------------------------")
                
                log(str(j+1)+":\n")
                log(str(top_3_values[j])+"\n")
                #log(str(vals[j])+"\n")
                log(str(tLat)+", "+str(tLon)+"\n")
                log("--------------------------------------------------------------------------------\n")
        
        targLats.append(newLats2)
        targLons.append(newLons2)
        print()
        
    return targLats,targLons


# Given a file to read from returns an array w/shape nResps+1 x nPercs
def getPercs(percFile='/home/iarsenea/trajs/targeting/percs.txt'):
    
    return np.loadtxt(percFile,delimiter=',')


def addLocs():
    reese = (33.593274, -102.029487, "Reese")
    amarillo = (35.192011, -101.837456, "Amarillo")
    clovis = (34.402701, -103.203542, "Clovis")
    wFalls = (33.906197, -98.491838, "Wichita Falls")
    midland = (31.990741, -102.081358, "Midland")
    abilene = (32.442535, -99.735791, "Abilene")
    plainview = (34.192500, -101.708145, "Plainview")
    hobbs = (32.707751, -103.135328, "Hobbs")
    roswell = (33.396840, -104.524914, "Roswell")
    canadian = (35.911075, -100.383835, "Canadian")
    dalhart = (36.052066, -102.517591, "Dalhart")
    floydada = (33.9845, -101.3377, "Floydada")
    
    places = (reese,amarillo,clovis,wFalls,midland,abilene,plainview,
              hobbs,roswell,canadian,dalhart,floydada)
    
    for place in places:
        plt.scatter(place[1],place[0],s=10,color='k')
        plt.text(place[1]+.05,place[0]+.05,place[2],fontsize=9)


# Given a forecast init date string YYYYMMDDHH, get the variance of the response values for each response function.
def getRVar(date):
    rFile = "/home/bancell/enkfALT/src/targout/Rvalues_"+date
    data = Dataset(rFile)
    
    vals = []
    for resp in np.arange(1,11,1):
        #vals.append(float(data["R"+str(resp)][0].data)) ##################
        vals.append(float(data["R"+str(resp)][-1].data))
    
    return vals

        
def writeLog(text,logFile = 'logTest.txt'):
    
    file = open(logFile,'a')
    file.write(text)
    file.close()
    

def clearLog(logFile = 'logTest.txt'):
    
    file = open(logFile,'w')
    file.close()
    
    
    
### Calculations ###

mmVal,pmVal,mmVar,pmVar = doCalcs(date)
targLats,targLons = chooseTargs(mmVal,pmVal,mmVar,pmVar,date)



### Plotting Values ###

# Get the location of the launch points
xs,ys = getLocs_test()
xs = xs/12000
ys = ys/12000

lats,lons=convertLatLon(xs,ys)

N = int(np.amax(ys)-np.amin(ys))+1 # yMax-yMin
M = int(np.amax(xs)-np.amin(xs))+1 # xMax-xMin

lats = np.asarray(lats).reshape(N,M)
lons = np.asarray(lons).reshape(N,M)

# Get data to plot state and province boundaries
states_provinces = cfeature.NaturalEarthFeature( # pulls data from "natural earth" company
        category='cultural',
        name='admin_1_states_provinces_lakes',
        scale='50m', # 1:50 million
        facecolor='none')

proj = ccrs.PlateCarree()
buff = .7

# Save all of the plots to a pdf
pFile = "/home/iarsenea/trajs/targeting/"+date+"/pdf_"+date+".pdf"
with PdfPages(pFile) as pdf:

    ### Plot Response Box ###

    rFile = "/home/bancell/enkfALT/src/targout/subsetTARG_"+date+".txt"

    # Read in the response file
    hour,minLon,maxLon,minLat,maxLat = np.loadtxt(rFile)

    # Plot that sucker
    # Get data to plot state and province boundaries
    states_provinces = cfeature.NaturalEarthFeature( # pulls data from "natural earth" company
            category='cultural',
            name='admin_1_states_provinces_lakes',
            scale='50m', # 1:50 million
            facecolor='none')

    # Set up figure
    fig = plt.figure(figsize=(11,8))

    # Plot the requested background
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    #ax.set_extent((-180, -50, 5, 70))
    ax.set_extent((minLon-5, maxLon+5, minLat-5, maxLat+5))
    ax.coastlines('50m', linewidth=0.55)
    ax.add_feature(states_provinces,edgecolor='black',linewidth=0.45)

    plt.plot((minLon,minLon,maxLon,maxLon,minLon),(minLat,maxLat,maxLat,minLat,minLat),c='r')
    plt.title("Run Date: "+date,loc="Left")
    plt.title(str(hour),loc="Right")
    pdf.savefig()
    plt.close()



    for iRes in range(10):

        #iRes = 0

        # Reshape the z values
        z1 = mmVal[iRes].reshape((N, M))
        z2 = pmVal[iRes].reshape((N, M))
        z3 = mmVar[iRes].reshape((N, M))
        z4 = pmVar[iRes].reshape((N, M))

        
        ### Plot % max above percentile and mean max ###

        fig = plt.figure(figsize=(10,8)) 

        # Plot the requested background
        ax = fig.add_subplot(111,projection=proj)
        ax.set_extent((np.amin(lons)-buff,np.amax(lons)+buff,np.amin(lats)-buff,np.amax(lats)+buff))
        ax.add_feature(states_provinces,edgecolor='k',linewidth=0.45)

        # Set up the gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='gray', linestyle='--')
        gl.xlabel_style = {'size': 7}
        gl.ylabel_style = {'size': 7}
        gl.xlabels_top = False
        gl.ylabels_right = False

        # Plot the values
        levs = np.arange(0,105,10)
        plt.contourf(lons,lats,z2,cmap="Spectral_r",levels=levs,extend="max")
        plt.contour(lons,lats,z1,alpha=.75)
        if np.amax(z1) > 0:
            plt.colorbar()
        
        plt.title("Response Function "+str(iRes+1), loc = "Left")
        plt.title("% Above 50th Percentile Fill, Mean Contour", loc = "Right")
        
        # Plot the three targets
        for j in range(len(targLons[iRes])):
            plt.scatter(targLons[iRes][j],targLats[iRes][j],marker="x",s=50,c='white',alpha=1-((j)/(len(targLons[iRes])*2)))
            plt.text(targLons[iRes][j]+.075,targLats[iRes][j]+.1,str(j+1),c='white',alpha=1-((j)/(len(targLons[iRes])*2)))
        
        addLocs()
        plt.grid(markevery=1)
        pdf.savefig()
        plt.close()
        
        
        ### Plot % % var above 15, and % var ###

        fig = plt.figure(figsize=(10,8)) 

        # Plot the requested background
        ax = fig.add_subplot(111,projection=proj)
        ax.set_extent((np.amin(lons)-buff,np.amax(lons)+buff,np.amin(lats)-buff,np.amax(lats)+buff))
        ax.add_feature(states_provinces,edgecolor='k',linewidth=0.45)

        # Set up the gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='gray', linestyle='--')
        gl.xlabel_style = {'size': 7}
        gl.ylabel_style = {'size': 7}
        gl.xlabels_top = False
        gl.ylabels_right = False

        # Plot the values
        plt.contourf(lons,lats,z4,cmap="Spectral_r",extend='max',levels=levs)
        if np.amax(z4) > 0:
            plt.colorbar()
        plt.contour(lons,lats,z3,cmap="YlGnBu_r",alpha = .8,linewidths=1.75,levels=np.arange(0,120,20))

        plt.title("Response Function "+str(iRes+1),loc="Left")
        plt.title("% > 15% Variance Fill, Mean % Variance Contour", loc="Right")
    
        addLocs()
        plt.grid(markevery=1)
        pdf.savefig()
        plt.close()

