"""
Author: Alex Crawford
Date Created: 6 Jan 2020
Date Modified: 6 Jan 2020
Purpose: Calculates the annual sea surface temperature from monthly values for
open water period.
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
ct = '75'

var = "oisst"
path = "/Volumes/Prospero/Arrigo"
inpath = path+"/"+var+"/Monthly"+ct
wtpath = path+"/nppcdays/Monthly"
outpath = path+"/"+var+"/Annual"+ct

# Time Variables
mos = list(range(5,9+1))
ymin, ymax = 1998, 2018
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
#### Prep with date/location ####
YY = str(ymin)+"_"+str(ymax)
years = list(range(ymin,ymax+1))

lats = np.fromfile(path+"/Projections/lats_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
lons = np.fromfile(path+"/Projections/lons_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))

#### FOR LOOP BY YEAR ####
ylist = []
for y in years:
    Y = str(y)
    print(Y)
    
    # Prep empty lists for each wind variable
    weights, ssts = [], []
    
    for m in mos:
        M = mons[m-1]
        
        # Identify number of ice free days in the month for weights
        wnc = nc.Dataset(wtpath+"/nppcdays"+ct+"_arctic_66N_Month"+M+"_1998_2018.nc")
        weight = wnc['nppcdays'+ct][y-ymin,:,:].data
        weights.append(weight)
        
        # Extract each variable, multiplying by the weights
        mnc = nc.Dataset(inpath+"/"+var+"_arctic_66N_Month"+M+"_1998_2018.nc")
        ssts.append( weight*mnc[var][y-ymin,:,:].data )
    
    # Take the sum for all months
    weightSum = np.apply_along_axis(np.nansum,0,weights)
    sstSum = np.apply_along_axis(np.nansum,0,ssts)
    
    # Divide each wind variable by the sum of the weights
    ylist.append( sstSum/weightSum )

### Write new netcdf file ###
fout = var+"_arctic_66N_Annual_"+YY+".nc"

ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf1.createDimension('y', lats.shape[0])
ncf1.createDimension('x', lats.shape[1])
ncf1.createDimension('time', len(years))

yNC = ncf1.createVariable('y', np.float32, ('y',))
xNC = ncf1.createVariable('x', np.float32, ('x',))
tNC = ncf1.createVariable('time', np.float32, ('time',))

chlNC = ncf1.createVariable(var, np.float32, ('time','y','x',))
latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))

ncf1.description = 'Average sea surface temperature during open water period'
ncf1.source = 'netCDF4 python module'
tNC.units = 'years'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
chlNC.units = 'degrees C'

tNC[:] = years
latNC[:] = lats
lonNC[:] = lons
chlNC[:] = np.array(ylist)

ncf1.close()    
        