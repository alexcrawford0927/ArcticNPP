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
ct = "75" # concentration threshold for sea ice

path = "/Volumes/Prospero/Arrigo"
inpath = path+"/prod"
icepath = '/Volumes/Prospero/Arrigo/ice/AdvanceRetreat/siphenology_C'+ct+'_1998-2018.nc'
outpath = path+"/nppcdays/Annual"
area1 = "prod_arctic_50N_Reg_v3_A_"
area2 = "prod_arctic_50N_Reg_v3_corr_S_"
ext = ".0_reanal1_d4_c2c90p0_prod.bin"

# Time Variables
mos = list(range(1,12+1))
ymin, ymax = 1998, 2018
dmin, dmax = 59, -92 # Limit to March through September
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

#### Read in Sea Ice Data ####
ice = nc.Dataset(icepath)

sts, eds, ncdayss = [], [], []

for y in years:
    # Set up Files
    Y = str(y)
    print(Y)
    
    if y < 2003:
        inpathY = inpath+"/"+area2+Y
    else:
        inpathY = inpath+"/"+area1+Y
        
    files = os.listdir(inpathY)
    files = [f for f in files if f.startswith(Y)]
    files = files[dmin:dmax] # Limit to months of interest
    
    # Initiate 
    ncst, nced = np.zeros_like(lats)*np.nan, np.zeros_like(lats)*np.nan
    
    # Assign start and end #
    for F in files:
        #print("--"+F[4:7])
        # Load each daily file
        f = int(F[4:7])
        arr = np.fromfile(inpathY+"/"+F,dtype='f').reshape((nrow,ncol))
        
        # Set NaNs
        arr = np.where( (arr == nanvalue) | (lats < minlat), np.nan, f )
        
        # Assign Start (if appicable) and End
        ncst = np.where( (np.isnan(ncst) == 1) & (arr == f), f, ncst )
        nced = np.where( (arr == f) & (np.isfinite(ncst) == 1), f, nced )
        
    # Calculate the number of continuous NPP days (as an intersection with ice sea)
    #### First assumption -- the starts and ends follow NPP only (always open water)
    st = ncst+0
    ed = nced+0
    
    #### Modify that assumption wherever the sea ice retreat occurs after the ncst 
    #######or the sea ice advance occurs before the nced
    lrd = ice['lrd'][y-ymin,:,:]
    fad = ice['fad'][y-ymin,:,:]

    st = np.where(lrd > ncst, lrd, st)
    st = np.where(lrd > nced, np.nan, st)
    ed = np.where(fad < nced, fad, ed)
    ed = np.where(fad < ncst, np.nan, ed)
        
    # The total number of days is inclusive, so add one to the difference
    ncdays = ed - st + 1
    ncdays = np.where( np.isfinite(ncdays) == 1, ncdays, 0 )
    
    sts.append(st), eds.append(ed), ncdayss.append(ncdays)
    
### Write new netcdf file ###
fout = "nppcdays"+ct+"_arctic_66N_Annual_"+YY+".nc"

ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf1.createDimension('y', lats.shape[0])
ncf1.createDimension('x', lats.shape[1])
ncf1.createDimension('time', len(sts))

yNC = ncf1.createVariable('y', np.float32, ('y',))
xNC = ncf1.createVariable('x', np.float32, ('x',))
tNC = ncf1.createVariable('time', np.float32, ('time',))

ncdaysNC = ncf1.createVariable('nppcdays'+ct, np.float32, ('time','y','x',))
ncstNC = ncf1.createVariable('ncst', np.float32, ('time','y','x',))
ncedNC = ncf1.createVariable('nced', np.float32, ('time','y','x',))
latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))

ncf1.description = 'Number of Ice-Free Days between first and last valid NPP'
ncf1.source = 'netCDF4 python module'
tNC.units = 'years'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
ncdaysNC.units = 'Days (per Year)'
ncstNC.units = 'NPP Ice-Free Start Day (DOY)'
ncedNC.units = 'NPP Ice-Free End Day (DOY)'

tNC[:] = years
latNC[:] = lats
lonNC[:] = lons
ncdaysNC[:] = np.array(ncdayss)
ncstNC[:] = np.array(sts)
ncedNC[:] = np.array(eds)

ncf1.close()
