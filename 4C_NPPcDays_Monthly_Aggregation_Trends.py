"""
Author: Alex Crawford
Date Created: 30 Dec 2019
Date Modified: 30 Dec 2019
Purpose: Calculate annual number of "continuous NPP days" for Arrigo data (from monthly netcdf files)
 and their trends. Works directly from raw Arrigo data.
"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import os
import netCDF4 as nc
import numpy as np
from scipy import stats
import MERRA_Module as md

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")
# File Variables
nrow = 2325 # Number of rows
ncol = 2014 # Number of columns
nanvalue = -99 # Value for no data
minlat = 66 # minimum latitude of reasonable data
ct = "75"

path = "/Volumes/Prospero/Arrigo"
inpath = path+"/nppcdays/Annual"
outpath = path+"/nppcdays/Monthly"
area1 = "prod_arctic_50N_Reg_v3_A_"
area2 = "prod_arctic_50N_Reg_v3_corr_S_"
ext = ".0_reanal1_d4_c2c90p0_prod.bin"

# Time Variables
mos = list(range(1,12+1))
ymin, ymax = 1998, 2018
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
dpm = [31,28,31,30,31,30,31,31,30,31,30,31]

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

#### Prep with date/location ####
YY = str(ymin)+"_"+str(ymax)
years = list(range(ymin,ymax+1))
lyb = md.leapyearBoolean(years)

lats = np.fromfile(path+"/Projections/lats_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
lons = np.fromfile(path+"/Projections/lons_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))

# Start Monthly Loop
for m in mos:
    
    ncsts, nceds, ncdayss = [], [], []
    
    for y in years:
        # Set up Leap Year
        Y = str(y)
        print(Y+mons[m-1])
 
        # Set up Leap Year       
        dpmly = dpm*1
        if md.leapyearBoolean([y])[0] == 1:
            dpmly[1] = 29
    
        # Set days (accounting for leap years)
        dmax = sum(dpmly[:m])
        dmin = dmax - dpmly[m-1]
        
        # Load files
        ncf = nc.Dataset(inpath+"/nppcdays"+ct+"_arctic_66N_Annual_"+YY+".nc")
        
        # Store 'start' and 'end'
        ast = ncf['ncst'][y-ymin,:,:]
        aed = ncf['nced'][y-ymin,:,:]
        ncf.close()

        # Initiate 
        ncst = np.where((ast <= dmin) & (aed > dmin), dmin, np.where((ast > dmin) & (ast <= dmax), ast, np.nan))
        nced = np.where(np.isnan(ast) == 1, np.nan, np.where(aed < dmax, aed, dmax))
        
        # Calculate the number of continuous NPP days
        ncdays = nced - ncst
        ncdays = np.where( np.isfinite(ncdays) == 1, ncdays, 0 )
        
        ncsts.append(ncst), nceds.append(nced), ncdayss.append(ncdays)
        
    ### Write new netcdf file ###
    fout = "nppcdays"+ct+"_arctic_66N_Month"+mons[m-1]+"_"+YY+".nc"
    
    ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
    ncf1.createDimension('y', lats.shape[0])
    ncf1.createDimension('x', lats.shape[1])
    ncf1.createDimension('time', len(ncsts))
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    ncdaysNC = ncf1.createVariable('nppcdays'+ct, np.float32, ('time','y','x',))
    ncstNC = ncf1.createVariable('ncst', np.float32, ('time','y','x',))
    ncedNC = ncf1.createVariable('nced', np.float32, ('time','y','x',))
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = 'Number of Continuous Days of NPP'
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'years'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'      
    ncdaysNC.units = 'Days (per Month)'
    ncstNC.units = 'NPP Start Day (DOY)'
    ncedNC.units = 'NPP end Day (DOY)'
    
    tNC[:] = years
    latNC[:] = lats
    lonNC[:] = lons
    ncdaysNC[:] = np.array(ncdayss)
    ncstNC[:] = np.array(ncsts)
    ncedNC[:] = np.array(nceds)
    
    ncf1.close()
