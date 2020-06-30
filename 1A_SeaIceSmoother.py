'''*********************************************
Authors: Alex Crawford
Date Created: 6/9/15
Date Modified: 1/28/19; 8/20/19 edited for Python 3
10/4/19 edited for Arrigo data
Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell of a particular sector each year.

Inputs: 
    concentration threshold (cts) -- value between 0 and 1
    years of interest (ymin, ymax) -- integers
    months for  -- the month on which to start the annual cycle -- March (3) is
        a good idea because it is the closest to the maximum
    
Outputs: A csv file with the concentration threshold and setor noted in the 
    file name. Retreat and advance days are recorded as a "DOY", with 1 being
    the 1st day of January. If retreat below the
    concentration threshold never occurs, the minimum day is recorded instead.
*********************************************'''
# Import clock:
from time import perf_counter as clock
# Start script stopwatch. The clock starts running when time is imported
start = clock()

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import os
import netCDF4 as nc
import numpy as np
import MERRA_Module as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

### Input Variables ###
cts = [0.15] # A number between 0 and 1 for the concentration threshold
n = 5 # Moving Average Size (n = # observation on either side of current day;
        ## so 1 = 3-point, 2 = 5-point, 4 = 9-point)

### Time Variables ###
ymin, ymax = 1998, 2018 # to 1998 to 2002
maxmo = [1,4] # months in which the sea ice maximum may occur
minmo = [8,10] # months in which the sea ice minimum may occur

### Path Variables ###
path = "/Volumes/Prospero/Arrigo"
inpath = path+"/ice"
outpath = path+"/ice_SmoothedMA"+str(n)

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Time Set Up
years = range(ymin,ymax+1)
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13",\
    "14","15","16","17","18","19","20","21","22","23","24","25","26","27",\
    "28","29","30","31"]
    
for y in years:
    Y = str(y)
    
    files = os.listdir(inpath+"/ice_arctic_50N_"+Y)
    files = [f for f in files if (f.endswith('.nc') and (f.startswith('.') == 0))]
    
    # Load example file for dimensions
    ncf = nc.Dataset(inpath+"/ice_arctic_50N_"+Y+"/"+files[0])
    arr = ncf.variables['seaice_conc'][:]
    
    # Generate Mask
    lats = ncf.variables['lat'][:]
    lons = ncf.variables['lon'][:]
    
    if len(lats.shape) == 1:
        lons, lats = np.meshgrid(lons,lats)
    
    validrows, validcols = np.where( (arr <= 1) & (lats >= 45) )
    
    del arr
    ncf.close()

    arr0 = []
    for f in files:
        ncf = nc.Dataset(inpath+"/ice_arctic_50N_"+Y+"/"+f)
        
        arr0.append(ncf.variables['seaice_conc'][:])
        
        ncf.close()
        
    # Smoothing
    arr0 = np.array(arr0)
    arrMA = np.zeros(arr0.shape)*np.nan
    print(" -- Smoothing at " + Y)
    for i in range(len(validrows)):
        arrMA[:,validrows[i],validcols[i]] = md.movingAverage2(arr0[:,validrows[i],validcols[i]],n)
    
    # Write new netcdf file
    fsmooth = f[:29]+"MA"+str(n)+"_"+Y+f[-31:]
    
    ncf1 = nc.Dataset(outpath+"/"+fsmooth, 'w', format='NETCDF4')
    ncf1.createDimension('y', arrMA.shape[1])
    ncf1.createDimension('x', arrMA.shape[2])
    ncf1.createDimension('time', arrMA.shape[0])
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    sicNC = ncf1.createVariable('siconc', np.float32, ('time','y','x',))
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = 'Sea Ice Concentration with 5-Day Moving Average'
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'day of year (Jan 1 = 1, Feb 1 = 32)' #e.g., 'days since 1979-01-01 00:00:00.0'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'      
    sicNC.units = 'Percentage'

    tNC[:] = range(1,arr0.shape[0]+1)
    latNC[:] = lats
    lonNC[:] = lons
    sicNC[:] = arrMA

    ncf1.close()
    
    del arrMA, arr0, validrows, validcols
