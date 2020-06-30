"""
Author: Alex Crawford
Date Created: 21 Jan 2020
Date Modified: 22 Jan 2020
Purpose: Calculates spatio-temporal averages for various regions of the 
"""
'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import netCDF4 as nc
import numpy as np
import pandas as pd

'''#####################
 Declare Variables
#####################'''
# File Variables
nrow = 2325 # Number of rows
ncol = 2014 # Number of columns
nanv = -999 # Value for no data
minlat = 66 # minimum latitude of reasonable data
nanmax = 9999999
ct='75'

rg = '2' # Region Version
op = '6' # Output Version

ymin, ymax = 1998, 2018

reglist = list(range(2,18+1,1)) #[0,100,200,400,800,1600,3200]

path = "/Volumes/Prospero/Arrigo"

'''#####################
Main Analysis
#####################'''
print("Main Analysis")
#### Prep with date/location ####
YY = str(ymin)+"_"+str(ymax)
years = list(range(ymin,ymax+1))

lats = np.fromfile(path+"/Projections/lats_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
lons = np.fromfile(path+"/Projections/lons_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
area = np.fromfile(path+"/Projections/arctic_50N_areas_2014x2325.bin",dtype='f')[1:-1].reshape((nrow,ncol))

# Load mask
regs = pd.read_pickle(path+"/Projections/Regions"+rg+"_4km.pkl")

# Create output data frame
pdf = pd.DataFrame(columns=["Reg","Year","npprate","npptot","chl","oisst",\
                      "icefreedays","hwf","uvavg","uavg","vavg","huavg","hvavg",\
                      "hWper","hSper","uWper","vSper","obsdays","area","opendays"])  
    
# Load NC files
NCnppcdays = nc.Dataset(path+"/nppcdays/Annual/nppcdays"+ct+"_arctic_66N_Annual_"+YY+".nc")
NCicefreedays = nc.Dataset(path+"/icefreedays/Annual"+ct+"/icefreedays_arctic_66N_Annual_"+YY+".nc")
NCnpprate =  nc.Dataset(path+"/npprate/Annual/npprate_arctic_66N_Annual_"+YY+".nc")
NCnpptot = nc.Dataset(path+"/npptot/Annual"+ct+"/npptot_Arctic_66N_Annual_"+YY+".nc")
NCchl = nc.Dataset(path+"/chl/Annual/chl_Arctic_66N_Annual_"+YY+".nc")
NCwind = nc.Dataset(path+"/wind4/Annual"+ct+"/wind4_Arctic_66N_Annual_"+YY+".nc")
NCoisst = nc.Dataset(path+"/oisst/Annual"+ct+"/oisst_Arctic_66N_Annual_"+YY+".nc")
NCow = nc.Dataset(path+"/openwater/Annual"+ct+"/openwater_Arctic_66N_Annual_"+YY+".nc")

for y in years:
    print("  " + str(y))
    # Extract nupy arrays
    nppcdays = NCnppcdays['nppcdays'+ct][y-ymin,:,:].data
    icefreedays = NCicefreedays['icefreedays'][y-ymin,:,:].data
    opendays = NCow['openwater'][y-ymin,:,:].data
    npprate = NCnpprate['npprate'][y-ymin,:,:].data
    npprate = np.where(npprate > nanmax, np.nan, npprate) # set NaNs for npprate
    npptot = NCnpptot['npptot'][y-ymin,:,:].data
    npptot = np.where(npptot > nanmax, np.nan, npptot) # set NaNs for npprate
    chl = NCchl['chl'][y-ymin,:,:].data
    oisst = NCoisst['oisst'][y-ymin,:,:].data
    hwf = NCwind['hwf'][y-ymin,:,:].data
    uvavg = NCwind['uvavg'][y-ymin,:,:].data
    uavg = NCwind['uavg'][y-ymin,:,:].data
    vavg = NCwind['vavg'][y-ymin,:,:].data
    huavg = NCwind['huavg'][y-ymin,:,:].data
    hvavg = NCwind['hvavg'][y-ymin,:,:].data
    hWper = NCwind['hWper'][y-ymin,:,:].data
    hSper = NCwind['hSper'][y-ymin,:,:].data
    uWper = NCwind['uWper'][y-ymin,:,:].data
    vSper = NCwind['vSper'][y-ymin,:,:].data
    
    for reg in reglist: # For each region
        # Identify subset of gridcells that are within the region
        regi = np.where(regs == reg) 
        regarea = np.sum( area[regi] )
        
        # Take sum for # of days
        avgdays = np.nansum( icefreedays[regi]*area[regi] ) / regarea
        obsdays = np.sum( nppcdays[regi]*area[regi] ) / regarea
        owdays = np.sum( opendays[regi]*area[regi] ) / regarea
        
        # Take an average that is weighted by both the number of opendays and the area of the grid cell
        avgnpp = np.nansum( npprate[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(npprate[regi])*opendays[regi]*area[regi] )
        sumnpp = np.nansum(npptot[regi]*area[regi])
        avgchl = np.nansum( chl[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(chl[regi])*opendays[regi]*area[regi] )
        avgoisst = np.nansum( oisst[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(oisst[regi])*opendays[regi]*area[regi] )
        avghwf = np.nansum( hwf[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(hwf[regi])*opendays[regi]*area[regi] )
        avguvavg = np.nansum( uvavg[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(uvavg[regi])*opendays[regi]*area[regi] )
        avguavg = np.nansum( uavg[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(uavg[regi])*opendays[regi]*area[regi] )
        avgvavg = np.nansum( vavg[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(vavg[regi])*opendays[regi]*area[regi] )
        avghuavg = np.nansum( huavg[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(huavg[regi])*opendays[regi]*area[regi] )
        avghvavg = np.nansum( hvavg[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(hvavg[regi])*opendays[regi]*area[regi] )
        avghWper = np.nansum( hWper[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(hWper[regi])*opendays[regi]*area[regi] )
        avghSper = np.nansum( hSper[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(hSper[regi])*opendays[regi]*area[regi] )
        avguWper = np.nansum( uWper[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(uWper[regi])*opendays[regi]*area[regi] )
        avgvSper = np.nansum( vSper[regi]*opendays[regi]*area[regi] ) / np.sum( np.isfinite(vSper[regi])*opendays[regi]*area[regi] )

        # Append to dataframe
        row = pd.DataFrame([{"Reg":reg,"Year":y,"npprate":avgnpp,"npptot":sumnpp,"icefreedays":avgdays,\
                       "chl":avgchl,"oisst":avgoisst,"hwf":avghwf,"uvavg":avguvavg,\
                       "uavg":avguavg,"vavg":avgvavg,"huavg":avghuavg,"hvavg":avghvavg,\
                   "hWper":avghWper,"hSper":avghSper,"uWper":avguWper,"vSper":avgvSper,\
                   "obsdays":obsdays,"area":regarea,"opendays":owdays},])
        pdf = pdf.append(row,ignore_index=True)

# Close netcdf files
NCnppcdays.close(), NCicefreedays.close(), NCnpprate.close(), NCnpptot.close(), NCchl.close(), NCoisst.close(), NCwind.close()

# Write to file
pdf.to_csv(path+"/Regional"+ct+"/ArrigoData_Regional_Annual_v"+op+".csv",index=False)
