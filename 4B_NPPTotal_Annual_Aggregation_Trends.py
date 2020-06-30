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
ct="75"

path = "/Volumes/Prospero/Arrigo"
var1 = "npprate"
var2 = "nppcdays"
outpath = path+"/npptot/Annual"+ct

# Time Variables
ymin, ymax = 1998, 2018

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

# Load Files
npprate = nc.Dataset(path+"/"+var1+"/Annual/"+var1+"_arctic_66N_Annual_"+YY+".nc")
nppcdays = nc.Dataset(path+"/"+var2+"/Annual/"+var2+ct+"_arctic_66N_Annual_"+YY+".nc")

# Calculate Total NPP by Year
npptot = np.array(npprate[var1][:]) * np.array(nppcdays[var2+ct][:])

### Write new netcdf file ###
fout = "npptot_arctic_66N_Annual_"+YY+".nc"

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

ncf1.description = 'Total Annual NPP based on the Arrigo Model'
ncf1.source = 'netCDF4 python module'
tNC.units = 'years'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
npptotNC.units = 'Total Annual NPP (mgC/m^2)'

tNC[:] = years
latNC[:] = lats
lonNC[:] = lons
npptotNC[:] = npptot

ncf1.close()
    
########################
### TREND TOTAL NPP ####
########################
#print(" -- Trend Calculation")
#b, a, r, p, e = np.zeros_like(lats), np.zeros_like(lats), \
#    np.zeros_like(lats), np.ones_like(lats), np.zeros_like(lats)
#
#rows, cols = np.where(np.isfinite(npptot[0,:,:]))
#
#for i in range(len(rows)):
#    ri, ci = rows[i], cols[i]
#    
#    xvals = np.array(years)[np.isfinite(npptot[:,ri,ci]) == 1]
#    yvals = npptot[:,ri,ci][np.isfinite(npptot[:,ri,ci]) == 1]
#    
#    b[ri,ci], a[ri,ci], r[ri,ci], p[ri,ci], e[ri,ci] = stats.linregress(xvals,yvals)
#
#### Write new netcdf file ###
#fout = "npptot_arctic_66N_Annual_Trend_"+YY+".nc"
#
#ncf2 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
#ncf2.createDimension('y', lats.shape[0])
#ncf2.createDimension('x', lats.shape[1])
#
#yNC = ncf2.createVariable('y', np.float32, ('y',))
#xNC = ncf2.createVariable('x', np.float32, ('x',))
#
#bNC = ncf2.createVariable('trend', np.float32, ('y','x',))
#aNC = ncf2.createVariable('intercept', np.float32, ('y','x',))
#rNC = ncf2.createVariable('rsquared', np.float32, ('y','x',))
#pNC = ncf2.createVariable('pvalue', np.float32, ('y','x',))
#eNC = ncf2.createVariable('stderr', np.float32, ('y','x',))
#latNC = ncf2.createVariable('lat', np.float32, ('y','x',))
#lonNC = ncf2.createVariable('lon', np.float32, ('y','x',))
#
#ncf2.description = 'Trend ('+YY+') in Total Annual Net Primary Productivity from the Arrigo Model'
#ncf2.source = 'netCDF4 python module'
#latNC.units = 'degrees north'
#lonNC.units = 'degrees east'      
#bNC.units = 'Trend (per yr) in Annual Sum of mg C/m^2'
#aNC.units = 'Intercept of in Annual Sum of mg C/m^2'
#rNC.units = 'r-squared value for trend'
#pNC.units = 'p-value for trend'
#eNC.units = 'standard error for trend'
#
#latNC[:] = lats
#lonNC[:] = lons
#bNC[:] = b
#aNC[:] = a
#rNC[:] = r**2
#pNC[:] = p
#eNC[:] = e
#
#ncf2.close()

##########################
### AVERAGE TOTAL NPP ####
##########################
print(" -- Average Calculation")
### Write new netcdf file ###
fout = "npptot_arctic_66N_Annual_Avg_"+YY+".nc"

ncf3 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf3.createDimension('y', lats.shape[0])
ncf3.createDimension('x', lats.shape[1])

yNC = ncf3.createVariable('y', np.float32, ('y',))
xNC = ncf3.createVariable('x', np.float32, ('x',))

stdNC = ncf3.createVariable('avg', np.float32, ('y','x',))
avgNC = ncf3.createVariable('std', np.float32, ('y','x',))
latNC = ncf3.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf3.createVariable('lon', np.float32, ('y','x',))

ncf3.description = 'Average ('+YY+') in Total Annual Net Primary Productivity from the Arrigo Model'
ncf3.source = 'netCDF4 python module'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
stdNC.units = 'Standard Deviation of Annual Sum of mg C/m^2'
avgNC.units = 'Average ('+YY+') in Annual Sum of mg C/m^2'

latNC[:] = lats
lonNC[:] = lons
stdNC[:] = np.apply_along_axis(np.nanstd,0,npptot)
avgNC[:] = np.apply_along_axis(np.nanmean,0,npptot)

ncf3.close()