"""
Author: Alex Crawford
Date Created: 14 May 2020
Date Modified: 14 May 2020
               
Purpose: Calculates the percentage of observations of wind coming from 
different directions from component winds by month for only the ice-free period.
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
minuv = 10 # Minimum wind speed in m/s for "high wind events"
minlat = 66
ct = '75'

path = "/Volumes/Prospero/Arrigo"
wpath = "/Volumes/Ferdinand/ERAI/wind/wind_EASE2_4km/10m/Hourly6"
outpath = path+"/wind4/Monthly"+ct
npppath = path+"/nppcdays"

# Time Variables
mos = list(range(9,9+1))
ymin, ymax = 1998, 2018
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
tres = 6 # Temporal resolution of wind data in hours
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

tscale = int(24/tres)

#########################
#### MONTHLY VALUES ####
for m in mos[:]:
    M = mons[m-1]
    print(' -'+M)
    
    # Load Start and End Points for NPP
    ncNPP = nc.Dataset(npppath+"/Monthly/nppcdays"+ct+"_arctic_66N_Month"+M+"_"+YY+".nc")
    ncst = ncNPP['ncst'][:].data
    nced = ncNPP['nced'][:].data
    
    dmin = np.nanmin(ncst)

    uvavgs, uavgs, vavgs, huavgs, hvavgs = [], [], [], [], []
    hwfs, uWpers, vSpers, hWpers, hSpers = [], [], [], [], []
    for y in years[:]:
        Y = str(y)
        yi = y-ymin
        print(' --'+Y)
        
        # Identify valid cells
        st, ed = ncst[yi,:,:].astype(int), nced[yi,:,:].astype(int)
        
        # Load times for Wind Data
        ncf = nc.Dataset(wpath+"/wind10m_EASE4km_"+Y+M+".nc")

        # Prep Arrays
        uvavg, uavg, vavg, huavg, hvavg = np.zeros_like(lats), np.zeros_like(lats), np.zeros_like(lats), np.zeros_like(lats), np.zeros_like(lats)
        hwf, uWper, vSper, hWper, hSper = np.zeros_like(lats), np.zeros_like(lats), np.zeros_like(lats), np.zeros_like(lats), np.zeros_like(lats)
        valid = np.zeros_like(lats)
        
        # Loop through every hour
        print("Loading Hourly Data")
        for h in range(ncf['time'][:].shape[0]):  
            #print(h)
            # Time Adjustment
            d = int(dmin+h/tscale)
            
            # Identify valid cells for this time
            valid = valid + np.where( (np.isfinite(ncst[yi,:,:]) == 0) | (ncst[yi,:,:] > d) | (nced[yi,:,:] <= d) | (lats < minlat), 0, 1)
                        
            # Extract wind data
            v = ncf['v'][h,:,:].data
            u = ncf['u'][h,:,:].data
            
            # Set NaNs for land and sea ice
            v = np.where( (np.isfinite(ncst[yi,:,:]) == 0) | (ncst[yi,:,:] > d) | (nced[yi,:,:] <= d) | (lats < minlat), 0, v ) 
            u = np.where( (np.isfinite(ncst[yi,:,:]) == 0) | (ncst[yi,:,:] > d) | (nced[yi,:,:] <= d) | (lats < minlat), 0, u ) 
            
            # Wind Analysis
            uv = np.sqrt( np.square(u) + np.square(v) )
            hw = uv > minuv
            
            # Aggregate to month by summing
            uvavg = uvavg + uv
            vavg = vavg + v
            uavg = uavg + u
            hwf = hwf + hw
            
            huavg = huavg + np.where(hw, u, 0)
            hvavg = hvavg + np.where(hw, v, 0)
            hWper = hWper + np.where(hw, (u > 0), 0)
            hSper = hSper + np.where(hw, (v > 0), 0)
            uWper = uWper + (u > 0)
            vSper = vSper + (v > 0)
        
        ncf.close()
        
        print("Taking Monthly Averages")
        
        # Wind Analysis
        valid = np.where(valid == 0, np.nan, valid)
        validhwf = np.where(valid == 0, np.nan, hwf)
        
        # Append arrays to lists
        uvavgs.append( uvavg/valid )
        uavgs.append( uavg/valid )
        vavgs.append( vavg/valid ) 
        hwfs.append( hwf/valid)
        huavgs.append( huavg/validhwf )
        hvavgs.append( hvavg/validhwf )
        hWpers.append( hWper/validhwf )
        hSpers.append( hSper/validhwf )
        uWpers.append( uWper/valid)
        vSpers.append( vSper/valid )
        
    ### Write new netcdf file ###
    fout = "wind4_arctic_66N_Month"+M+"_"+YY+".nc"
    
    ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
    ncf1.createDimension('y', lats.shape[0])
    ncf1.createDimension('x', lats.shape[1])
    ncf1.createDimension('time', len(years))
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    uvavgNC = ncf1.createVariable('uvavg', np.float32, ('time','y','x'))
    uavgNC = ncf1.createVariable('uavg', np.float32, ('time','y','x'))
    vavgNC = ncf1.createVariable('vavg', np.float32, ('time','y','x'))
    hwfNC = ncf1.createVariable('hwf', np.float32, ('time','y','x'))
    huavgNC = ncf1.createVariable('huavg', np.float32, ('time','y','x'))
    hvavgNC = ncf1.createVariable('hvavg', np.float32, ('time','y','x'))
    hWperNC = ncf1.createVariable('hWper', np.float32, ('time','y','x'))
    hSperNC = ncf1.createVariable('hSper', np.float32, ('time','y','x'))
    uWperNC = ncf1.createVariable('uWper', np.float32, ('time','y','x'))
    vSperNC = ncf1.createVariable('vSper', np.float32, ('time','y','x'))
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = '''Includes various wind variabiles, including the high-
    wind fraction and the percentage of observations with westerly or northerly 
    wind (including during high-wind events) -- all based on an ice-free period.'''
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'years'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'
    uvavgNC.units = 'm/s'
    uavgNC.units = 'm/s'
    vavgNC.units = 'm/s'
    hwfNC.units = 'ratio'
    huavgNC.units = 'm/s'
    hvavgNC.units = 'm/s'    
    hWperNC.units = 'ratio'
    hSperNC.units = 'ratio'
    uWperNC.units = 'ratio'
    vSperNC.units = 'ratio'
    
    tNC[:] = years
    latNC[:] = lats
    lonNC[:] = lons
    uvavgNC[:] = np.array(uvavgs)
    uavgNC[:] = np.array(uavgs)
    vavgNC[:] = np.array(vavgs)
    hwfNC[:] = np.array(hwfs)
    huavgNC[:] = np.array(huavgs)
    hvavgNC[:] = np.array(hvavgs)
    hWperNC[:] = np.array(hWpers)
    hSperNC[:] = np.array(hSpers)
    uWperNC[:] = np.array(uWpers)
    vSperNC[:] = np.array(vSpers)
    
    ncf1.close()
