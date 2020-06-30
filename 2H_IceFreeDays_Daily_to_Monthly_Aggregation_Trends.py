"""
Author: Alex Crawford
Date Created: 4 Oct 2019
Date Modified: 18 Jun 2020
Purpose: Calculate monthly  totals of NPP for Arrigo data and their trends.
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
nmin = 15
ct = '75'

path = "/Volumes/Prospero/Arrigo"
inpath3 = path+"/ice/AdvanceRetreat/siphenology_C"+ct+"_1998-2018.nc"
outpath = path+"/icefreedays/Monthly"+ct

# Time Variables
mos = list(range(5,9+1))
ymin, ymax = 1998, 2018
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
dpm = [31,28,31,30,31,30,31,31,30,31,30,31]

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

## Load Ice days
ice = nc.Dataset(inpath3)
lrd = ice.variables['lrd'][:]
fad = ice.variables['fad'][:]
ice.close()

#### Prep with date/location ####
YY = str(ymin)+"_"+str(ymax)
years = list(range(ymin,ymax+1))

lats = np.fromfile(path+"/Projections/lats_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
lons = np.fromfile(path+"/Projections/lons_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))

b1list, a1list, r1list, p1list, e1list = [], [], [], [], []

#########################
#### MONTHLY VALUES ####
for m in mos[:]:
    M = mons[m-1]
    print(' -'+M)
    
    icefreedays = []
    for y in years:
        # Set up Leap Year       
        dpmly = dpm*1
        if md.leapyearBoolean([y])[0] == 1:
            dpmly[1] = 29
    
        # Set days (accounting for leap years)
        dmax = sum(dpmly[:m])
        dmin = dmax - dpmly[m-1]
        
        # Store 'start' and 'end'
        ast = lrd[y-ymin,:,:]
        aed = fad[y-ymin,:,:]
    
        # Initiate
        ist = np.where((ast <= dmin) & (aed > dmin), dmin, np.where((ast > dmin) & (ast <= dmax), ast, np.nan))
        ied = np.where(np.isnan(ast) == 1, np.nan, np.where(aed < dmax, aed, dmax))
        
        # Calculate the number of continuous ice-free days
        days = ied - ist + 1
        days = np.where( np.isfinite(days) == 1, days, 0 )
        
        icefreedays.append(days)
    
    icefreedays = np.array(icefreedays)

    ### Write new netcdf file ###
    fout = "icefreedays_arctic_66N_Month"+M+"_"+YY+".nc"
    
    ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
    ncf1.createDimension('y', icefreedays.shape[1])
    ncf1.createDimension('x', icefreedays.shape[2])
    ncf1.createDimension('time', icefreedays.shape[0])
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    iceNC = ncf1.createVariable('icefreedays', np.float32, ('time','y','x'))
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = 'Ice-free days on the grid from the Arrigo Model'
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'years'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'      
    iceNC.units = '# of days'

    tNC[:] = years
    latNC[:] = lats
    lonNC[:] = lons
    iceNC[:] = icefreedays

    ncf1.close()
    
    ## Loading necessary if doing trends separately ##
    # ncf1 = nc.Dataset(outpath+"/icefreedays_arctic_66N_Month"+M+"_"+YY+".nc")
    # icefreedays = ncf1.variables['icefreedays'][:]

    #### MONTHLY TRENDS ####
    # print(" -- Trend Calculation")
    # rows, cols = np.where( np.apply_along_axis(np.sum, 0, np.ma.getmask(icefreedays)) > nmin)

    # b1, a1, r1, p1, e1 = np.zeros_like(icefreedays[0,:,:])*np.nan, np.zeros_like(icefreedays[0,:,:])*np.nan, \
    #     np.zeros_like(icefreedays[0,:,:])*np.nan, np.ones_like(icefreedays[0,:,:])*np.nan, np.zeros_like(icefreedays[0,:,:])*np.nan
        
    # for i in range(len(rows)):
    #     ri, ci = rows[i], cols[i]
    #     b1[ri,ci], a1[ri,ci], r1[ri,ci], p1[ri,ci], e1[ri,ci] = stats.linregress(years,icefreedays[:,ri,ci])
    
    # b1list.append(b1), a1list.append(a1), r1list.append(r1), p1list.append(p1), e1list.append(e1)
    

### Write new netcdf file ###
# fout = "icefreedays_arctic_66N_Monthly_Trend_"+YY+".nc"

# ncf3 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
# ncf3.createDimension('y', icefreedays.shape[1])
# ncf3.createDimension('x', icefreedays.shape[2])
# ncf3.createDimension('time', len(mos))

# yNC = ncf3.createVariable('y', np.float32, ('y',))
# xNC = ncf3.createVariable('x', np.float32, ('x',))
# tNC = ncf3.createVariable('time', np.float32, ('time',))

# bNC = ncf3.createVariable('trend', np.float32, ('time','y','x',))
# aNC = ncf3.createVariable('intercept', np.float32, ('time','y','x',))
# rNC = ncf3.createVariable('rsquared', np.float32, ('time','y','x',))
# pNC = ncf3.createVariable('pvalue', np.float32, ('time','y','x',))
# eNC = ncf3.createVariable('stderr', np.float32, ('time','y','x',))
# latNC = ncf3.createVariable('lat', np.float32, ('y','x',))
# lonNC = ncf3.createVariable('lon', np.float32, ('y','x',))

# ncf3.description = 'Trend ('+YY+') in number of open days from the Arrigo Model'
# ncf3.source = 'netCDF4 python module'
# tNC.units = 'months'
# latNC.units = 'degrees north'
# lonNC.units = 'degrees east'      
# bNC.units = 'Trend (per yr) in number of open days'
# aNC.units = 'Intercept of trend in number of open days'
# rNC.units = 'r-squared value for trend'
# pNC.units = 'p-value for trend'
# eNC.units = 'standard error for trend'

# tNC[:] = mos
# latNC[:] = lats
# lonNC[:] = lons
# bNC[:] = np.array(b1list)
# aNC[:] = np.array(a1list)
# rNC[:] = np.array(r1list)**2
# pNC[:] = np.array(p1list)
# eNC[:] = np.array(e1list)

# ncf3.close()