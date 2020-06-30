"""
Author: Alex Crawford
Date Created: 4 Oct 2019
Date Modified: 30 Dec 2019
Purpose: Calculate the number of days with valid satellite observations of NPP
for Arrigo data and their trends.
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

path = "/Volumes/Prospero/Arrigo"
inpath = path+"/prod"
outpath = path+"/nppdays/Monthly"
area1 = "prod_arctic_50N_Reg_v3_A_"
area2 = "prod_arctic_50N_Reg_v3_corr_S_"
ext = ".0_reanal1_d4_c2c90p0_prod.bin"

# Time Variables
mos = list(range(3,10+1))
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
for m in mos[:1]:
    M = mons[m-1]
    print(' -'+M)
        
    sumMlist = []
    
    for y in years:
        Y = str(y)
        
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
        
        dmin = sum(dpm[:(m-1)])+lymin # EXCLUSIVE
        dmax = sum(dpm[:(m)])+lymax # INCLUSIVE
        
        if y < 2003:
            inpathY = inpath+"/"+area2+Y
        else:
            inpathY = inpath+"/"+area1+Y
        
        files = os.listdir(inpathY)
        files = [f for f in files if f.startswith(Y)]
        files = [f for f in files if ( (int(f[:7]) > (y*1000)+dmin) and (int(f[:7]) <= (y*1000)+dmax) ) ]
        
        ### LOAD FILES AND CALCULATE NUMBER OF VALID DAYS PER MONTH ###
        print(" -- Summing "+Y)
        arrMlist = []
        
        for f in files:
            arr = np.fromfile(inpathY+"/"+f,dtype='f').reshape((nrow,ncol))
            arrMlist.append( np.where((arr == nanvalue) | (lats < minlat), np.nan, arr) )
        
        arrM = np.array(arrMlist)
        del arrMlist
        
        sumMlist.append( np.apply_along_axis( np.sum,0,np.isfinite(arrM) ) )
    
    sumM = np.array(sumMlist)
    
    ### Write new netcdf file ###
    fout = "nppdays_arctic_66N_Month"+M+"_"+YY+".nc"
    
    #
    #ncf = nc.Dataset(outpath+"/old"+fout)
    #sumM = ncf.variables['prod'][:]
    #
    
    ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
    ncf1.createDimension('y', sumM.shape[1])
    ncf1.createDimension('x', sumM.shape[2])
    ncf1.createDimension('time', sumM.shape[0])
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    arrNC = ncf1.createVariable('nppdays', np.float32, ('time','y','x',))
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = 'Number of Days with a valid (ice-free, cloud-free) NPP value from the Arrigo Model'
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'years'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'      
    arrNC.units = '# of days'

    tNC[:] = years
    latNC[:] = lats
    lonNC[:] = lons
    arrNC[:] = sumM

    ncf1.close()

    #### MONTHLY TRENDS ####
    print(" -- Trend Calculation")
    b, a, r, p, e = np.zeros_like(sumM[0,:,:]), np.zeros_like(sumM[0,:,:]), \
        np.zeros_like(sumM[0,:,:]), np.ones_like(sumM[0,:,:]), np.zeros_like(sumM[0,:,:])
    
    rows, cols = np.where(np.isfinite(sumM[0,:,:]))
    
    for i in range(len(rows)):
        ri, ci = rows[i], cols[i]
        b[ri,ci], a[ri,ci], r[ri,ci], p[ri,ci], e[ri,ci] = stats.linregress(years,sumM[:,ri,ci])
    
    blist.append(b), alist.append(a), rlist.append(r), plist.append(p), elist.append(e)

### Write new netcdf file ###
fout = "nppdays_arctic_66N_Monthly_Trend_"+YY+".nc"

ncf2 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf2.createDimension('y', sumM.shape[1])
ncf2.createDimension('x', sumM.shape[2])
ncf2.createDimension('time', len(mos))

yNC = ncf2.createVariable('y', np.float32, ('y',))
xNC = ncf2.createVariable('x', np.float32, ('x',))
tNC = ncf2.createVariable('time', np.float32, ('time',))

bNC = ncf2.createVariable('trend', np.float32, ('time','y','x',))
aNC = ncf2.createVariable('intercept', np.float32, ('time','y','x',))
rNC = ncf2.createVariable('rsquared', np.float32, ('time','y','x',))
pNC = ncf2.createVariable('pvalue', np.float32, ('time','y','x',))
eNC = ncf2.createVariable('stderr', np.float32, ('time','y','x',))
latNC = ncf2.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf2.createVariable('lon', np.float32, ('y','x',))

ncf2.description = 'Trend ('+YY+') in Monthly Number of Days with valid NPP from the Arrigo Model'
ncf2.source = 'netCDF4 python module'
tNC.units = 'months'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
bNC.units = 'Trend (per yr) in number of valid NPP days'
aNC.units = 'Intercept of in number of valid NPP days'
rNC.units = 'r-squared value for trend'
pNC.units = 'p-value for trend'
eNC.units = 'standard error for trend'

tNC[:] = mos
latNC[:] = lats
lonNC[:] = lons
bNC[:] = np.array(blist)
aNC[:] = np.array(alist)
rNC[:] = np.array(rlist)**2
pNC[:] = np.array(plist)
eNC[:] = np.array(elist)

ncf2.close()
  