'''
Author: Alex Crawford
Date Created: 15 Aug 2018
Date Modified: 27 Nov 2018

Aggregates post-MIZ bloom npp, chl, sst, par for each grid cell for each 8-day period,
stores results as large pickled numpy arrays.
'''
'''############
Set up Modules
############'''
from osgeo import gdal, gdalnumeric
import os
import numpy as  np
from netCDF4 import Dataset
import MERRA_Module as md
import pandas as pd

'''############
Declare Variables
############'''
path = "/Volumes/Prospero/MODIS"
inpath0 = path+"/VGPM/VGPM_PostBloom"
inpath1 = path+"/VGPM/VGPM8_EASE2_25km"
inpath2 = path+"/par/par8_EASE2_25km"
inpath3 = path+"/chl/chl8_EASE2_25km"
inpath4 = path+"/sst/sst8_EASE2_25km"
inpath5 = "/Volumes/Ferdinand/ERAI/wind/uvAb_EASE2_25km/10m/Value/Hourly6"
inpath6 = "/Volumes/Ferdinand/ERAI/wind/u_EASE2_25km/10m/Value/Hourly6"
inpath7 = "/Volumes/Ferdinand/ERAI/wind/v_EASE2_25km/10m/Value/Hourly6"

path2 = "/Volumes/Prospero/SeaIce/AdvanceRetreat_EASE2/C10"
outpath = path+"/VGPM/VGPM_PostBloom8/loccsv_sw4"

nan2 = -99
dtype = gdal.GDT_Float32
latmin = 60
LatN = "/Volumes/Prospero/Projections/EASE2_N60_25km_Lats.tif"
LonN = "/Volumes/Prospero/Projections/EASE2_N60_25km_Lons.tif"

dpf = 8. # The number of days per file
bloomdays = 20 # Minimum number of days for the initial MIZ bloom
windmin = 10 # Minimum wind speed needed to be classified as a "wind event"
ymin, ymax = 2011, 2014
startdoy = 65 # Should match the sea ice data (e.g., if sea ice starts March 1, 
                # then start as close to March 1st as possible)
startmonth = 3 # Should match the sea ice data
refdate = [1900,1,1,0,0,0]
fdays = np.array([65,73,81,89,97,105,113,121,129,137,145,153,161,169,177,185,193,\
                  201,209,217,225,233,241,249,257,265,273,281,289,297,305,313,321,\
                  329,337,345,353,361,366,374,382,390,398,406,414,422])

mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
days = [ "0"+str(i) for i in range(1,10)] + [str(i) for i in range(10,32) ]
hours = [ "0"+str(i)+"00" for i in range(10)] + [str(i)+"00" for i in range(10,24) ]

from mpl_toolkits.basemap import interp

# Load New Projection Info
refout = gdal.Open("/Volumes/Ferdinand/Projections/EASE2_N60_25km_Lats.tif")
gt =  refout.GetGeoTransform()
xout = np.arange(gt[0],refout.RasterXSize*gt[1]+gt[0],gt[1])
yout = np.arange(gt[3],refout.RasterYSize*gt[5]+gt[3],gt[5])

# Create 2-D versions of x and y grids
x2out = np.transpose( np.repeat(xout,refout.RasterYSize).reshape(refout.RasterYSize,refout.RasterXSize) )
y2out = np.repeat(yout,refout.RasterXSize).reshape(refout.RasterYSize,refout.RasterXSize)
   
# Load Current Projection Info
refin = gdal.Open("/Volumes/Ferdinand/Projections/EASE2_N0_100km_Lats.tif")
gt =  refin.GetGeoTransform()
xin = np.arange(gt[0],refin.RasterXSize*gt[1]+gt[0],gt[1])
yin = np.arange(gt[3],refin.RasterYSize*gt[5]+gt[3],gt[5])
    
'''################
Main Analysis
################'''
ref = gdal.Open(LatN)
lats = gdalnumeric.LoadFile(LatN)
lons = gdalnumeric.LoadFile(LonN)

#######################################
### Identify Open Post-Bloom Period ###
#######################################
# Load LRD to find valid ocean gridcells
lrds = []
for y in range(ymin,ymax+1):
    lrd = gdalnumeric.LoadFile(path2+"/LRD_EASE2/LastRetreatDay_"+str(y)+".tif")[49:-50,49:-50]
    lrds.append( np.where( (lrd <= nan2) | (lats < latmin), np.nan, lrd) )
    
lrdfin = np.isfinite(lrds)
lrdfinsum = np.apply_along_axis(np.sum,0,lrdfin)
rows_, cols_ = np.where(lrdfinsum > 0)

# Loop through each location, making a pandas dataframe for which rows are times
for i in range(11492,len(rows_)):#len(rows_)): 0, 4000, 8000, ... 11935
    ri = rows_[i]
    ci = cols_[i]
    
    pdf = pd.DataFrame(columns=['lat','lon','year','doy','opendays','npp','chl','par','sst','hwe','u','v'])
    
    print ri, "  ", ci
    for y in range(ymin,ymax+1):
        Y = str(y)
        dpy = 365 + md.leapyearBoolean([y])[0]
        
        # Load open water period
        lrd = lrds[y-ymin][ri,ci]

        edoy = gdalnumeric.LoadFile(inpath0+"/enddoy/enddoy_"+Y+".tif")[ri,ci]
        sdoy = gdalnumeric.LoadFile(inpath0+"/startdoy/startdoy_"+Y+".tif")[ri,ci]
        
        edate = md.timeAdd([y,1,0,0,0,0],[0,0,edoy,0,0,0])
        sdate = md.timeAdd([y,1,0,0,0,0],[0,0,sdoy,0,0,0])
        
        if (sdoy < 0) | (edoy < 0) | (sdoy >= edoy):
            continue
        
        else:
            darr = fdays[np.where((fdays >= sdoy) & (fdays < edoy))]
            
            # Prep Lists
            npp, chl, sst, par, hwe, u, v, hu, hv, uv = [], [], [], [], [], [], [], [], [], []
            
            for d in range(len(darr)):
                if darr[d] < 10:
                    D = "00"+str(darr[d])
                elif darr[d] < 100:
                    D = "0"+str(darr[d])
                else:
                    D = str(darr[d])
                
                # Load NPP, Chl, SST, PAR files
                try:
                    npp.append( gdalnumeric.LoadFile(inpath1+"/"+Y+"/VGPM."+Y+D+"_EASE2_25km_sw_to_md.tif")[ri,ci] )
                    par.append( gdalnumeric.LoadFile(inpath2+"/"+Y+"/par."+Y+D+"_EASE2_25km_sw_to_md.tif")[ri,ci] )
                    chl.append( gdalnumeric.LoadFile(inpath3+"/"+Y+"/chl."+Y+D+"_EASE2_25km_sw_to_md.tif")[ri,ci] )
                    sst.append( gdalnumeric.LoadFile(inpath4+"/"+Y+"/sst."+Y+D+"_EASE2_25km_sw_to_md.tif")[ri,ci] )
                except:
                    if (darr[d]-365) < 10:
                        D2 = "00"+str(darr[d]-365)
                    elif (darr[d]-365) < 100:
                        D2 = "0"+str(darr[d]-365)
                    else:
                        D2 = str(darr[d]-365)
                    
                    if y == 2002:
                        npp.append( gdalnumeric.LoadFile(inpath1+"/"+str(y+1)+"/VGPM_"+str(y+1)+D2+"_EASE2_25km.tif")[ri,ci] )
                        par.append( gdalnumeric.LoadFile(inpath2+"/"+str(y+1)+"/par."+str(y+1)+D2+"_EASE2_25km.tif")[ri,ci] )
                        chl.append( gdalnumeric.LoadFile(inpath3+"/"+str(y+1)+"/chl."+str(y+1)+D2+"_EASE2_25km.tif")[ri,ci] )
                        sst.append( gdalnumeric.LoadFile(inpath4+"/"+str(y+1)+"/sst."+str(y+1)+D2+"_EASE2_25km.tif")[ri,ci] )
                    else:
                        npp.append( gdalnumeric.LoadFile(inpath1+"/"+str(y+1)+"/VGPM."+str(y+1)+D2+"_EASE2_25km_sw_to_md.tif")[ri,ci] )
                        par.append( gdalnumeric.LoadFile(inpath2+"/"+str(y+1)+"/par."+str(y+1)+D2+"_EASE2_25km_sw_to_md.tif")[ri,ci] )
                        chl.append( gdalnumeric.LoadFile(inpath3+"/"+str(y+1)+"/chl."+str(y+1)+D2+"_EASE2_25km_sw_to_md.tif")[ri,ci] )
                        sst.append( gdalnumeric.LoadFile(inpath4+"/"+str(y+1)+"/sst."+str(y+1)+D2+"_EASE2_25km_sw_to_md.tif")[ri,ci] )
                
            # Load Wind Data
            alist, ulist, vlist, tlist = [], [], [], []
            if edoy <= dpy:
                for m in range(sdate[1],edate[1]+1):
                    nc = Dataset(inpath5+"/ubAb10m_EASE25km_"+Y+mons[m-1]+".nc")
                    ncu = Dataset(inpath6+"/u10m_EASE25km_"+Y+mons[m-1]+".nc")
                    ncv = Dataset(inpath7+"/v10m_EASE25km_"+Y+mons[m-1]+".nc")
                    
                    alist.append(nc.variables['uvab'][:,ri,ci])
                    ulist.append(ncu.variables['u'][:,ri,ci])
                    vlist.append(ncv.variables['v'][:,ri,ci])
                    tlist = tlist + list(nc.variables['time'][:])
                    nc.close()
            else:
                for m in range(sdate[1],12+1):
                    nc = Dataset(inpath5+"/ubAb10m_EASE25km_"+Y+mons[m-1]+".nc")
                    ncu = Dataset(inpath6+"/u10m_EASE25km_"+Y+mons[m-1]+".nc")
                    ncv = Dataset(inpath7+"/v10m_EASE25km_"+Y+mons[m-1]+".nc")
                    
                    alist.append(nc.variables['uvab'][:,ri,ci])
                    ulist.append(ncu.variables['u'][:,ri,ci])
                    vlist.append(ncv.variables['v'][:,ri,ci])
                    tlist = tlist + list(nc.variables['time'][:])
                    nc.close()
            
                for m in range(1,edate[1]+2):
                    nc = Dataset(inpath5+"/ubAb10m_EASE25km_"+str(y+1)+mons[m-1]+".nc")
                    ncu = Dataset(inpath6+"/u10m_EASE25km_"+str(y+1)+mons[m-1]+".nc")
                    ncv = Dataset(inpath7+"/v10m_EASE25km_"+str(y+1)+mons[m-1]+".nc")
                    
                    alist.append(nc.variables['uvab'][:,ri,ci])
                    ulist.append(ncu.variables['u'][:,ri,ci])
                    vlist.append(ncv.variables['v'][:,ri,ci])
                    tlist = tlist + list(nc.variables['time'][:])
                    nc.close()
            
            arr = np.concatenate(alist,axis=0)
            uarr = np.concatenate(ulist,axis=0)
            varr = np.concatenate(vlist,axis=0)
            tarr = np.array(tlist)
            
            if y == 2003:
                tarr[np.where(tarr == 0)] = np.arange(909024,909054,6)
            elif y == 2007:
                tarr[np.where(tarr == 0)] = np.arange(940440,940464,6)
            
            del alist, m, tlist, ulist, vlist
            
            # Identify high wind events
            arr2 = np.where(arr > windmin, 1, np.nan)
            
            # Append Wind Data to Lists
            for d in range(len(darr)):
                # Create indices for determining start and end times for wind
                sd = md.timeAdd([y,1,0,0,0,0],[0,0,darr[d],0,0,0])
                ed = md.timeAdd([y,1,0,0,0,0],[0,0,darr[d]+8,0,0,0])
                sti = np.where(tarr == md.daysBetweenDates(refdate,sd)*24)[0][0]
                eti = np.where(tarr == md.daysBetweenDates(refdate,ed)*24)[0][0]
                    
                # Sum wind event occurrences, converting to units of "event days"
                hwe.append( np.nansum(arr2[sti:eti])/4. )
                hu.append( np.nanmean(uarr[sti:eti]*arr2[sti:eti]) )
                hv.append( np.nanmean(varr[sti:eti]*arr2[sti:eti]) )
                
                u.append( np.mean(uarr[sti:eti]) )
                v.append( np.mean(varr[sti:eti]) )
                uv.append( np.mean(arr[sti:eti]) )
                
            # Add a new row!
            for d in range(len(darr)):
                newrow = pd.DataFrame([{'lat':lats[ri,ci],'lon':lons[ri,ci],\
                                    'year':y,'doy':darr[d],'opendays':darr[d]-lrd,\
                                    'npp':npp[d],'chl':chl[d],'sst':sst[d],'par':par[d],\
                                    'hwe':hwe[d],'hu':hu[d], 'hv':hv[d], 'uv':uv[d], \
                                    'u':u[d], 'v':v[d]},])
                pdf = pdf.append(newrow,ignore_index=True)
        
    pdf.to_csv(outpath+"/PostBloom_8Day_Values_Row"+str(ri)+"_Col"+str(ci)+".csv", index=False)