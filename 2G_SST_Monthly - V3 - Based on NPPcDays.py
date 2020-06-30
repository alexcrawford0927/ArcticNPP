"""
Author: Alex Crawford
Date Created: 2 Jan 2020
Date Modified: 3 Jan 2020
Purpose: Calculates the average SST by month
for only the ice-free period.
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
var = "oisst"
ext = "_arctic_50N_"+str(ncol)+"x"+str(nrow)+".bin"
nanvalue = -999
minlat = 66
nmin = 15
ct = '75'

path = "/Volumes/Prospero/Arrigo"
npppath = path+"/nppcdays"
inpath = path+"/"+var
outpath = path+"/"+var+"/Monthly"+ct

# Time Variables
mos = list(range(6,9+1))
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

blist, alist, rlist, plist, elist = [], [], [], [], []

#########################
#### MONTHLY VALUES ####
for m in mos[:]:
    M = mons[m-1]
    print(' -'+M)
    
    # Load Start and End Points for NPP
    ncNPP = nc.Dataset(npppath+"/Monthly/nppcdays"+ct+"_arctic_66N_Month"+M+"_"+str(ymin)+"_"+str(ymax)+".nc")
    ncst = ncNPP['ncst'][:].data
    nced = ncNPP['nced'][:].data
    
    nced[np.isfinite(ncst) == 0] = np.nan
    
    mlist = []
    for y in years:
        Y = str(y)
        print (' --'+Y)
        yi = y-ymin
        
        ### IDENTIFY CORRECT FILES FOR MONTH/YEAR COMBO ###
        # Account for leap years
        if m == 1: # No impact on January
            lymin = 0
            lymax = 0
        elif m == 2: # Impact on February only at end
            lymin = 0
            lymax = lyb[y-ymin]
        else: # Impact on beginning and end for other months
            lymin = lyb[y-ymin]
            lymax = lyb[y-ymin]
        
        dmin = sum(dpm[:(m-1)])+lymin+1 # INCLUSIVE
        dmax = sum(dpm[:(m)])+lymax+1 # EXCLUSIVE

        ### START DAILY LOOP ###
        ylist = []
        for d in range(dmin,dmax):
            # Load daily data
            YD = str((y*1000)+d)
            arr = np.fromfile(inpath+"/"+var+"_arctic_50N_"+Y+"/"+YD+"_"+var+ext,dtype='f').reshape((nrow,ncol))
        
            # Set NaNs for land and sea ice
            arr = np.where((arr == nanvalue) | (ncst[yi,:,:] > d) | (nced[yi,:,:] <= d) | (lats < minlat), np.nan, arr) 
            
            # Append to Monthly list
            ylist.append(arr)
        
        # Take average for the month
        mlist.append( np.apply_along_axis(np.nanmean,0,np.array(ylist)) )
    
    marr = np.array(mlist)
        
    ### Write new netcdf file ###
    fout = var+"_arctic_66N_Month"+M+"_"+YY+".nc"
    
    ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
    ncf1.createDimension('y', lats.shape[0])
    ncf1.createDimension('x', lats.shape[1])
    ncf1.createDimension('time', len(years))
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    arrNC = ncf1.createVariable(var, np.float32, ('time','y','x',))
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = 'Average SST during open water period'
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'years'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'      
    arrNC.units = 'Degrees C'
    
    tNC[:] = years
    latNC[:] = lats
    lonNC[:] = lons
    arrNC[:] = marr
    
    ncf1.close()
#    
#    #### MONTHLY TRENDS ####
#    print(" -- Trend Calculation")
#    b, a, r, p, e = np.zeros_like(lats), np.zeros_like(lats), \
#        np.zeros_like(lats), np.ones_like(lats), np.zeros_like(lats)
#    
#    marrfin = np.apply_along_axis(np.sum,0,np.isfinite(marr))
#    
#    rows, cols = np.where(marrfin >= nmin)
#    
#    for i in range(len(rows)):
#        ri, ci = rows[i], cols[i]
#        
#        xvals = np.array(years)[np.isfinite(marr[:,ri,ci]) == 1]
#        yvals = marr[np.isfinite(marr[:,ri,ci]) == 1,ri,ci]
#        
#        b[ri,ci], a[ri,ci], r[ri,ci], p[ri,ci], e[ri,ci] = stats.linregress(xvals,yvals)
#    
#    blist.append(b), alist.append(a), rlist.append(r), plist.append(p), elist.append(e)
#        
#### Write new netcdf file ###
#fout = var+"_arctic_66N_Monthly_Trend_"+YY+".nc"
#
#ncf2 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
#ncf2.createDimension('y', lats.shape[0])
#ncf2.createDimension('x', lats.shape[1])
#ncf2.createDimension('time', len(mos))
#
#yNC = ncf2.createVariable('y', np.float32, ('y',))
#xNC = ncf2.createVariable('x', np.float32, ('x',))
#tNC = ncf2.createVariable('time', np.float32, ('time',))
#
#bNC = ncf2.createVariable('trend', np.float32, ('time','y','x',))
#aNC = ncf2.createVariable('intercept', np.float32, ('time','y','x',))
#rNC = ncf2.createVariable('rsquared', np.float32, ('time','y','x',))
#pNC = ncf2.createVariable('pvalue', np.float32, ('time','y','x',))
#eNC = ncf2.createVariable('stderr', np.float32, ('time','y','x',))
#latNC = ncf2.createVariable('lat', np.float32, ('y','x',))
#lonNC = ncf2.createVariable('lon', np.float32, ('y','x',))
#
#ncf2.description = 'Trend ('+YY+') in Monthly SST'
#ncf2.source = 'netCDF4 python module'
#tNC.units = 'months'
#latNC.units = 'degrees north'
#lonNC.units = 'degrees east'      
#bNC.units = 'Trend (per yr) in Monthly SST (degrees C)'
#aNC.units = 'Intercept of in Monthly SST (degrees C)'
#rNC.units = 'r-squared value for trend'
#pNC.units = 'p-value for trend'
#eNC.units = 'standard error for trend'
#
#tNC[:] = mos
#latNC[:] = lats
#lonNC[:] = lons
#bNC[:] = np.array(blist)
#aNC[:] = np.array(alist)
#rNC[:] = np.array(rlist)**2
#pNC[:] = np.array(plist)
#eNC[:] = np.array(elist)
#
#ncf2.close()