"""
Author: Alex Crawford
Date Created: 1 Nov 2019
Date Modified: 19 Dec 2019
Purpose: 
"""

################
# Load Modules & Functions
import numpy as np
import netCDF4 as nc
import os

###############
# Declare Variables

# File Variables
nrow = 2325 # Number of rows
ncol = 2014 # Number of columns
maxval = 10**15
ct='75'

mos = list(range(5,9+1))
ymin, ymax = 1998, 2018
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
dpm = [31,28,31,30,31,30,31,31,30,31,30,31]

path = "/Volumes/Prospero/Arrigo"
inpath = path+"/openwater/Monthly"+ct
outpath = path+"/openwater/Annual"+ct

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

YY = str(ymin)+"_"+str(ymax)
years = range(ymin,ymax+1)

lats = np.fromfile(path+"/Projections/lats_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
lons = np.fromfile(path+"/Projections/lons_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))

siclist = []
for y in years:
    Y = str(y)
    print(Y)
    
    owsum = np.zeros_like(lats)
    for m in mos:
        M = mons[m-1]
        
        ncow = nc.Dataset(inpath+"/openwater_Arctic_66N_Month"+M+"_"+YY+".nc")
        ow = ncow['openwater'][y-ymin,:,:].data
        
        owsum += ow
        
    # Append to list for the given month/year
    siclist.append( owsum )
    
### Write new netcdf file ###
fout = "openwater_arctic_66N_Annual_"+YY+".nc"

ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf1.createDimension('y', lats.shape[0])
ncf1.createDimension('x', lons.shape[1])
ncf1.createDimension('time', ymax-ymin+1)

yNC = ncf1.createVariable('y', np.float32, ('y',))
xNC = ncf1.createVariable('x', np.float32, ('x',))
tNC = ncf1.createVariable('time', np.float32, ('time',))

iceNC = ncf1.createVariable('openwater', np.float32, ('time','y','x'))
latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))

ncf1.description = 'Open water on the grid from the Arrigo Model'
ncf1.source = 'netCDF4 python module'
tNC.units = 'years'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
iceNC.units = 'Sum of Annual %SIC'

tNC[:] = years
latNC[:] = lats
lonNC[:] = lons
iceNC[:] = np.array(siclist)

ncf1.close()
    
    
