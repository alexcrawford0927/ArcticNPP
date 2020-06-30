"""
Author: Alex Crawford
Date Created: 4 Oct 2019
Date Modified: 4 Oct 2019
Purpose: Calculate monthly  totals of NPP for Arrigo data and their trends.
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
minlat = 66 # minimum latitude of reasonable data
maxNnans = 6 # maximum number of nans allowed in a time series for trend calculations

path = "/Volumes/Prospero/Arrigo"
inpath1 = path+"/nppdays/Monthly"
inpath2 = path+"/prod/Monthly"
outpath = path+"/npprateOnly/Monthly"
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
for m in mos[:]:
    M = mons[m-1]
    print(' -'+M)
    
    # Calculate the average rate of NPP during days with valid data
    print(' -- NPPrate')
    npp = nc.Dataset(inpath2+"/prod_arctic_66N_Month"+M+"_1998_2018.nc")
    nppdays = nc.Dataset(inpath1+"/nppdays_arctic_66N_Month"+M+"_1998_2018.nc")
    npprate = npp.variables['prod'][:] / nppdays.variables['nppdays'][:]
    npp.close(), nppdays.close()
    
    ### Write new netcdf file ###
    fout = "npprate_arctic_66N_Month"+M+"_"+YY+".nc"
    
    ncf1 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
    ncf1.createDimension('y', npprate.shape[1])
    ncf1.createDimension('x', npprate.shape[2])
    ncf1.createDimension('time', npprate.shape[0])
    
    yNC = ncf1.createVariable('y', np.float32, ('y',))
    xNC = ncf1.createVariable('x', np.float32, ('x',))
    tNC = ncf1.createVariable('time', np.float32, ('time',))
    
    arrNC = ncf1.createVariable('npprate', np.float32, ('time','y','x',))
    latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
    lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))
    
    ncf1.description = 'Average Monthly NPP rate (mg C / m^2 / day) during valid days from the Arrigo Model'
    ncf1.source = 'netCDF4 python module'
    tNC.units = 'years'
    latNC.units = 'degrees north'
    lonNC.units = 'degrees east'      
    arrNC.units = 'mg C / m^2 / day'

    tNC[:] = years
    latNC[:] = lats
    lonNC[:] = lons
    arrNC[:] = npprate

    ncf1.close()
    
    ## Loading necessary if doing trends separately ##
    ncf1 = nc.Dataset(outpath+"/npprate_arctic_66N_Month"+M+"_"+YY+".nc")
    npprate = ncf1.variables['npprate'][:]

    #### MONTHLY TRENDS ####
    print(" -- Trend Calculation")
    b, a, r, p, e = np.zeros_like(npprate[0,:,:]), np.zeros_like(npprate[0,:,:]), \
        np.zeros_like(npprate[0,:,:]), np.ones_like(npprate[0,:,:]), np.zeros_like(npprate[0,:,:])
    
    rows, cols = np.where( np.apply_along_axis(np.sum, 0, np.ma.getmask(npprate)) <= maxNnans)
    
    for i in range(len(rows)):
        ri, ci = rows[i], cols[i]
        
        xvals = np.array(years)[np.ma.getmask(npprate[:,ri,ci]) == 0]
        yvals = npprate[:,ri,ci][np.ma.getmask(npprate[:,ri,ci]) == 0]
        
        b[ri,ci], a[ri,ci], r[ri,ci], p[ri,ci], e[ri,ci] = stats.linregress(xvals,yvals)
    
    blist.append(b), alist.append(a), rlist.append(r), plist.append(p), elist.append(e)

### Write new netcdf file ###
fout = "npprate_arctic_66N_Monthly_Trend_"+YY+".nc"

ncf2 = nc.Dataset(outpath+"/"+fout, 'w', format='NETCDF4')
ncf2.createDimension('y', npprate.shape[1])
ncf2.createDimension('x', npprate.shape[2])
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

ncf2.description = 'Trend ('+YY+') in Average Monthly NPP rate for Days with valid NPP from the Arrigo Model'
ncf2.source = 'netCDF4 python module'
tNC.units = 'months'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
bNC.units = 'Trend (per yr) in NPP rate [ (mg C / m^2 / day) ]'
aNC.units = 'Intercept of in NPP rate [ (mg C / m^2 / day) ]'
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

    