"""
Author: Alex Crawford
Date Created: 30 Dec 2019
Date Modified: 30 Dec 2019
Purpose: Calculate estimate of total NPP based on annual NPP Rate and annual 
NPP start and end dates. Requires that NPPcDays and NPPRate be calculated first.
"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import netCDF4 as nc
import numpy as np
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
ct = '75'

path = "/Volumes/Prospero/Arrigo"
var1 = "npprate"
var2 = "nppcdays"
outpath = path+"/npptot/Monthly"+ct

# Time Variables
mos = list(range(5,9+1))
ymin, ymax = 1998, 2018
dmin, dmax = 59, -92 # Limit to March through September
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
dpm = [31,28,31,30,31,30,31,31,30,31,30,31]

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

##################
### TOTAL NPP ####
##################
print(" -- Total Calculation")

# Prep with date/location
YY = str(ymin)+"_"+str(ymax)
years = list(range(ymin,ymax+1))
lyb = md.leapyearBoolean(years)

lats = np.fromfile(path+"/Projections/lats_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
lons = np.fromfile(path+"/Projections/lons_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))

for m in mos:
    M = mons[m-1]
    print(M)
    
    # Load Files
    npprate = nc.Dataset(path+"/"+var1+"/Monthly/"+var1+"_arctic_66N_Month"+M+"_"+YY+".nc")
    nppcdays = nc.Dataset(path+"/"+var2+"/Monthly/"+var2+ct+"_arctic_66N_Month"+M+"_"+YY+".nc")
    
    # Calculate Total NPP by Year
    npptot = np.array(npprate[var1][:]) * np.array(nppcdays[var2+ct][:])
    
    ### Write new netcdf file ###
    fout = "npptot_arctic_66N_Month"+M+"_"+YY+".nc"
    
    ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
    ncf1.createDimension('y', lats.shape[0])
    ncf1.createDimension('x', lats.shape[1])
    ncf1.createDimension('time', len(years))
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    npptotNC = ncf1.createVariable('npptot', np.float32, ('time','y','x',))
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = 'Total Monthly NPP based on the Arrigo Model'
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'years'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'      
    npptotNC.units = 'Total Monthly NPP (mgC/m^2)'
    
    tNC[:] = years
    latNC[:] = lats
    lonNC[:] = lons
    npptotNC[:] = npptot
    
    ncf1.close()
