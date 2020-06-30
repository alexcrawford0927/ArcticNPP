"""
Author: Alex Crawford
Date Created: 2 Jan 2020
Date Modified: 6 Jan 2020
Purpose: Calculates the high wind fraction, average wind speed and wind 
direction, and average and high-wind direction from component winds for multi-
month period.
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

path = "/Volumes/Prospero/Arrigo"
inpath = path+"/wind4/Monthly"+ct
wtpath = path+"/nppcdays/Monthly"
outpath = path+"/wind4/Annual"+ct

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
hwf, uvavg, uavg, vavg, huavg, hvavg = [], [], [], [], [], []
hWper, hSper, uWper, vSper = [], [], [], []
for y in years:
    Y = str(y)
    print(Y)
    
    # Prep empty lists for each wind variable
    weights, hwfs, uvavgs, uavgs, vavgs, huavgs, hvavgs = [], [], [], [], [], [], []
    hWpers, hSpers, uWpers, vSpers = [], [], [], []
    
    for m in mos:
        M = mons[m-1]
        
        # Identify number of ice free days in the month for weights
        wnc = nc.Dataset(wtpath+"/nppcdays"+ct+"_arctic_66N_Month"+M+"_1998_2018.nc")
        weight = wnc['nppcdays'+ct][y-ymin,:,:].data
        weights.append(weight)
        
        # Extract each wind variable, multiplying by the weights
        mnc = nc.Dataset(inpath+"/wind4_arctic_66N_Month"+M+"_1998_2018.nc")
        hwfs.append( weight*mnc['hwf'][y-ymin,:,:].data )
        uvavgs.append( weight*mnc['uvavg'][y-ymin,:,:].data )
        uavgs.append( weight*mnc['uavg'][y-ymin,:,:].data )
        vavgs.append( weight*mnc['vavg'][y-ymin,:,:].data )
        huavgs.append( weight*mnc['huavg'][y-ymin,:,:].data )
        hvavgs.append( weight*mnc['hvavg'][y-ymin,:,:].data )
        uWpers.append( weight*mnc['uWper'][y-ymin,:,:].data )
        vSpers.append( weight*mnc['vSper'][y-ymin,:,:].data )
        hWpers.append( weight*mnc['hWper'][y-ymin,:,:].data )
        hSpers.append( weight*mnc['hSper'][y-ymin,:,:].data )
        
    # Take the sum for all months
    weightSum = np.apply_along_axis(np.nansum,0,weights)
    hwfSum = np.apply_along_axis(np.nansum,0,hwfs)
    uvavgSum = np.apply_along_axis(np.nansum,0,uvavgs)
    uavgSum = np.apply_along_axis(np.nansum,0,uavgs)
    vavgSum = np.apply_along_axis(np.nansum,0,vavgs)
    huavgSum = np.apply_along_axis(np.nansum,0,huavgs)
    hvavgSum = np.apply_along_axis(np.nansum,0,hvavgs)
    uWperSum = np.apply_along_axis(np.nansum,0,uWpers)
    vSperSum = np.apply_along_axis(np.nansum,0,vSpers)
    hWperSum = np.apply_along_axis(np.nansum,0,hWpers)
    hSperSum = np.apply_along_axis(np.nansum,0,hSpers)
    
    # Divide each wind variable by the sum of the weights
    hwf.append(hwfSum/weightSum)
    uvavg.append(uvavgSum/weightSum)
    uavg.append(uavgSum/weightSum)
    vavg.append(vavgSum/weightSum)
    huavg.append(huavgSum/weightSum)
    hvavg.append(hvavgSum/weightSum)
    uWper.append(uWperSum/weightSum)
    vSper.append(vSperSum/weightSum)
    hWper.append(hWperSum/weightSum)
    hSper.append(hSperSum/weightSum)
    
### Write new netcdf file ###
fout = "wind4_arctic_66N_Annual_"+YY+".nc"

ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf1.createDimension('y', lats.shape[0])
ncf1.createDimension('x', lats.shape[1])
ncf1.createDimension('time', len(years))

yNC = ncf1.createVariable('y', np.float32, ('y',))
xNC = ncf1.createVariable('x', np.float32, ('x',))
tNC = ncf1.createVariable('time', np.float32, ('time',))

hwfNC = ncf1.createVariable('hwf', np.float32, ('time','y','x',))
uvavgNC = ncf1.createVariable('uvavg', np.float32, ('time','y','x'))
uavgNC = ncf1.createVariable('uavg', np.float32, ('time','y','x'))
vavgNC = ncf1.createVariable('vavg', np.float32, ('time','y','x'))
huavgNC = ncf1.createVariable('huavg', np.float32, ('time','y','x'))
hvavgNC = ncf1.createVariable('hvavg', np.float32, ('time','y','x'))
uWperNC = ncf1.createVariable('uWper', np.float32, ('time','y','x'))
vSperNC = ncf1.createVariable('vSper', np.float32, ('time','y','x'))
hWperNC = ncf1.createVariable('hWper', np.float32, ('time','y','x'))
hSperNC = ncf1.createVariable('hSper', np.float32, ('time','y','x'))
latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))

ncf1.description = 'Average wind, wind components, and high wind fraction during open water period'
ncf1.source = 'netCDF4 python module'
tNC.units = 'years'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
hwfNC.units = 'fraction of days'
uvavgNC.units = 'm/s'
uavgNC.units = 'm/s'
vavgNC.units = 'm/s'
huavgNC.units = 'm/s'
hvavgNC.units = 'm/s'
uWperNC.units = 'ratio'
vSperNC.units = 'ratio'
hWperNC.units = 'ratio'
hSperNC.units = 'ratio'

tNC[:] = years
latNC[:] = lats
lonNC[:] = lons
hwfNC[:] = np.array(hwf)
uvavgNC[:] = np.array(uvavg)
uavgNC[:] = np.array(uavg)
vavgNC[:] = np.array(vavg)
huavgNC[:] = np.array(huavg)
hvavgNC[:] = np.array(hvavg)
uWperNC[:] = np.array(uWper)
vSperNC[:] = np.array(vSper)
hWperNC[:] = np.array(hWper)
hSperNC[:] = np.array(hSper)

ncf1.close()    
        