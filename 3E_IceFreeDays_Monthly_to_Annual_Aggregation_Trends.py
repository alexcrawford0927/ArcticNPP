"""
Author: Alex Crawford
Date Created: 4 Oct 2019
Date Modified: 18 Jun 2020
Purpose: Calculate annual rate of NPP for Arrigo data (from monthly netcdf files)
 and their trends. Must run the daily to monthly script first. Weighted by the number of 
 nppdays.
"""

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import netCDF4 as nc
import numpy as np
from scipy import stats

'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")
# File Variables
nrow = 2325 # Number of rows
ncol = 2014 # Number of columns
ncvar = 'icefreedays'
ct = '75'

path = "/Volumes/Prospero/Arrigo"
inpath = path+"/ice/AdvanceRetreat/"
outpath = path+"/"+ncvar+"/Annual"+ct

ymin, ymax = 1998, 2018
mos = list(range(5,9+1))
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

### Load sea ice ###
ice = nc.Dataset(inpath+"/siphenology_C"+ct+"_"+str(ymin)+"-"+str(ymax)+".nc")

icefreedays = ice['opc'][:]

#### Write to File ####
fout = ncvar+"_arctic_66N_Annual_"+YY+".nc"

ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf1.createDimension('y', lats.shape[0])
ncf1.createDimension('x', lats.shape[1])
ncf1.createDimension('time', len(years))

yNC = ncf1.createVariable('y', np.float32, ('y',))
xNC = ncf1.createVariable('x', np.float32, ('x',))
tNC = ncf1.createVariable('time', np.float32, ('time',))

arrNC = ncf1.createVariable(ncvar, np.float32, ('time','y','x',))
latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))

ncf1.description = 'Number of Ice Free Days during season of interest'
ncf1.source = 'netCDF4 python module'
tNC.units = 'years'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
arrNC.units = 'days'

tNC[:] = years
latNC[:] = lats
lonNC[:] = lons
arrNC[:] = icefreedays

ncf1.close()

#### Calculate Trend ####
print(" - Calculating Trend")
b, a, r, p, e = np.zeros_like(lats), np.zeros_like(lats), \
    np.zeros_like(lats), np.ones_like(lats), np.zeros_like(lats)

rows, cols = np.where(np.isfinite(icefreedays[0,:,:]))

for i in range(len(rows)):
    ri, ci = rows[i], cols[i]
    b[ri,ci], a[ri,ci], r[ri,ci], p[ri,ci], e[ri,ci] = stats.linregress(years,icefreedays[:,ri,ci])

#### Write to File ####
fout = ncvar+"_arctic_66N_Annual_Trend_"+YY+".nc"

ncf3 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf3.createDimension('y', lats.shape[0])
ncf3.createDimension('x', lats.shape[1])

yNC = ncf3.createVariable('y', np.float32, ('y',))
xNC = ncf3.createVariable('x', np.float32, ('x',))

bNC = ncf3.createVariable('trend', np.float32, ('y','x',))
aNC = ncf3.createVariable('intercept', np.float32, ('y','x',))
rNC = ncf3.createVariable('rsquared', np.float32, ('y','x',))
pNC = ncf3.createVariable('pvalue', np.float32, ('y','x',))
eNC = ncf3.createVariable('stderr', np.float32, ('y','x',))
latNC = ncf3.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf3.createVariable('lon', np.float32, ('y','x',))

ncf3.description = 'Trend ('+YY+') in Annual Number of Ice Free Days in Season of Interest'
ncf3.source = 'netCDF4 python module'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
bNC.units = 'Trend (per yr) in Seasonal Ice Free Days [days] per yr'
aNC.units = 'Intercept of in Seasonal Ice Free Days [days]'
rNC.units = 'r-squared value for trend'
pNC.units = 'p-value for trend'
eNC.units = 'standard error for trend'

latNC[:] = lats
lonNC[:] = lons
bNC[:] = b
aNC[:] = a
rNC[:] = r**2
pNC[:] = p
eNC[:] = e

ncf3.close()
