'''
Author: Alex Crawford
Date Created: 13 Aug 2018
Date Modified: 13 Sep 2018

Aggregates post-MIZ bloom npp, chl, sst, par for each grid cell for each year,
stores results as geotiffs
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
inpath1 = path+"/VGPM/vgpm8_EASE2_25km"
inpath2 = path+"/par/par8_EASE2_25km"
inpath3 = path+"/chl/chl8_EASE2_25km"
inpath4 = path+"/sst/sst8_EASE2_25km"
inpath5 = "/Volumes/Ferdinand/ERAI/wind/uvAb_EASE2_25km/10m/Value/Hourly6"
inpath6 = "/Volumes/Ferdinand/ArcticCyclone/detection10_7E/CycloneFields" # path+"/VGPM/countA8" # 
inpath7 ="/Volumes/Miranda/OceanCurrents/BeringStrait/BeringStrait_Monthly_HeatTempVol_2017_Z.csv"

path2 = "/Volumes/Miranda/SeaIce/AdvanceRetreat_EASE2/C10"
outpath = "/Volumes/Prospero/MODIS/VGPM/VGPM_PostBloom"
outpath6 = path+"/VGPM/countA8"

nan2 = -99
nanlt = 0
dtype = gdal.GDT_Float32
latmin = 60
LatN = "/Volumes/Prospero/Projections/EASE2_N60_25km_Lats.tif"

dpf = 8. # The number of days per file
bloomdays = 20 # Minimum number of days for the initial MIZ bloom
windmin = 10 # Minimum wind speed needed to be classified as a "wind event"
ymin, ymax = 1998, 2003
startdoy = 65 # Should match the sea ice data (e.g., if sea ice starts March 1, 
                # then start as close to March 1st as possible)
startmonth = 3 # Should match the sea ice data
refdate = [1900,1,1,0,0,0]

mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

'''################
Main Analysis
################'''
ref = gdal.Open(LatN)
lats = gdalnumeric.LoadFile(LatN)

for y in range(ymin,ymax+1):
    Y = str(y)
    dpy = 365 + md.leapyearBoolean([y])[0]
    
    #######################################
    ### Identify Open Post-Bloom Period ###
    #######################################
    # Load LRD and FAD to obtain open water period
    lrd = gdalnumeric.LoadFile(path2+"/LRD_EASE2/LastRetreatDay_"+Y+".tif")[49:-50,49:-50]
    fad = gdalnumeric.LoadFile(path2+"/FAD_EASE2/FirstAdvanceDay_"+Y+".tif")[49:-50,49:-50]
    lrd = np.where( (lrd <= nan2) | (lats < latmin), np.nan, lrd)
    fad = np.where( (fad <= nan2) | (lats < latmin), np.nan, fad)
    
    # Identify valid grid cells for year
    start = np.ceil((lrd + bloomdays - 1)/dpf)-np.ceil((startdoy-1)/dpf)
    # Subtract 1 to convert to an index; also subtract to account for non-calendar year start date
    
    # If there's a wrap-around to the next calendar year, one of the files may have fewer days
    end = np.where( fad <= dpy, np.round(fad/dpf), 1 + np.round(fad/dpf - (dpy/dpf % 1)) ) - np.ceil((startdoy-1)/dpf)
    # Do not subtract 1 because Python is an exclusive subsetter; subtract to account for non-calendar year start date
    
    rows, cols = np.where(start < end)
    
    ########################
    ### Calculate OP NPP ###
    ########################
    print Y + ": NPP"
    # Load NPP file lists
    flist0 = [f for f in os.listdir(inpath1+"/"+Y) if ( (int(f[9:12]) >= startdoy) and f.startswith('.') == 0 and f.endswith(".tif") )]
    flist1 = [f for f in os.listdir(inpath1+"/"+str(y+1)) if ( (int(f[9:12]) < startdoy) and f.startswith('.') == 0 and f.endswith(".tif") )]
    flist0.sort()
    flist1.sort()
    
    fdays = [int(f[9:12]) for f in flist0] + [int(f[9:12])+dpy for f in flist1]

    # Load NPP files
    alist = []
    for f in flist0[:-1]:
        # Multiply by 8 to get 8-day total instead of a rate
        alist.append( dpf * gdalnumeric.LoadFile(inpath1+"/"+Y+"/"+f) )
    # Append final file, which might only has 5 days (or 6 in leap years...)
    lday = int(f[8:11])
    if (dpy - lday) < dpf:
        alist.append( (dpy - lday) * gdalnumeric.LoadFile(inpath1+"/"+Y+"/"+flist0[-1]) )
    else:
        alist.append( dpf * gdalnumeric.LoadFile(inpath1+"/"+Y+"/"+flist0[-1]) )
    for f in flist1:
        # Multiply by 8 to get 8-day total instead of a rate
        alist.append( dpf * gdalnumeric.LoadFile(inpath1+"/"+str(y+1)+"/"+f) )
    
    arr = np.array(alist)
    arr = np.where(arr < nanlt, np.nan, arr)
    del alist, flist0, flist1, lday
    
    # Sum NPP
    npp = np.zeros_like(arr[0])*np.nan
    sdoys = np.zeros_like(arr[0])*np.nan
    edoys = np.zeros_like(arr[0])*np.nan
    endnpp = np.zeros_like(arr[0])*np.nan
    for i in range(len(rows)): # units mg C / m^2 / "Period"
        npp[rows[i],cols[i]] = np.nansum(arr[int(start[rows[i],cols[i]]):int(end[rows[i],cols[i]]),rows[i],cols[i]])
        if np.sum(np.isfinite(arr[int(start[rows[i],cols[i]]):int(end[rows[i],cols[i]]),rows[i],cols[i]])) > 0:
            sdoys[rows[i],cols[i]] = fdays[int(start[rows[i],cols[i]])]
            edoys[rows[i],cols[i]] = fdays[np.where(np.isfinite(arr[:int(end[rows[i],cols[i]]),rows[i],cols[i]]) == 1)[0][-1]+1]
            endnpp[rows[i],cols[i]] = np.where(np.isfinite(arr[:int(end[rows[i],cols[i]]),rows[i],cols[i]]) == 1)[0][-1]
    
    md.writeNumpy_gdalObj(sdoys,outpath+"/startdoy/startdoy_"+Y+".tif",ref,dtype)
    md.writeNumpy_gdalObj(edoys,outpath+"/enddoy/enddoy_"+Y+".tif",ref,dtype)
    md.writeNumpy_gdalObj(edoys-sdoys,outpath+"/pbap/pbap_"+Y+".tif",ref,dtype)
    md.writeNumpy_gdalObj(fad-sdoys,outpath+"/pbop/pbop_"+Y+".tif",ref,dtype)
    
    md.writeNumpy_gdalObj(npp,outpath+"/VGPM/VGPM_"+Y+".tif",ref,dtype)
    md.writeNumpy_gdalObj(npp/(edoys-sdoys),outpath+"/VGPMrate/VGPMrate_"+Y+".tif",ref,dtype)
    del arr, npp, i
    
    # All following variables are summed/averaged only for the period that is both post-bloom and during non-nan npp
    rows, cols = np.where(start < endnpp)
    
    ########################
    ### Calculate OP PAR ###
    ########################
    print Y + ": PAR"
    # Load PAR file lists
    flist0 = [f for f in os.listdir(inpath2+"/"+Y) if ( (int(f[8:11]) >= startdoy) and f.startswith('.') == 0 and f.endswith(".tif") )]
    flist1 = [f for f in os.listdir(inpath2+"/"+str(y+1)) if ( (int(f[8:11]) < startdoy) and f.startswith('.') == 0 and f.endswith(".tif") )]
    flist0.sort()
    flist1.sort()
    
    # Load PAR files
    alist = []
    for f in flist0[:-1]:
        # Multiply by 8 to get 8-day total instead of a rate
        alist.append( dpf * gdalnumeric.LoadFile(inpath2+"/"+Y+"/"+f) )
    # Append final file, which might only has 5 days (or 6 in leap years...)
    lday = int(f[8:11])
    if (dpy - lday) < dpf:
        alist.append( (dpy - lday) * gdalnumeric.LoadFile(inpath2+"/"+Y+"/"+flist0[-1]) )
    else:
        alist.append( dpf * gdalnumeric.LoadFile(inpath2+"/"+Y+"/"+flist0[-1]) )
    for f in flist1:
        # Multiply by 8 to get 8-day total instead of a rate
        alist.append( dpf * gdalnumeric.LoadFile(inpath2+"/"+str(y+1)+"/"+f) )
    
    arr = np.array(alist)
    arr = np.where(arr < nanlt, np.nan, arr)
    del alist, flist0, flist1, lday

    # Sum PAR
    par = np.zeros_like(arr[0])*np.nan
    for i in range(len(rows)): # units mol photons / m^2 / "Period"
        par[rows[i],cols[i]] = np.nansum(arr[int(start[rows[i],cols[i]]):int(endnpp[rows[i],cols[i]])+1,rows[i],cols[i]])
    
    md.writeNumpy_gdalObj(par,outpath+"/par/par_"+Y+".tif",ref,dtype)
    md.writeNumpy_gdalObj(par/(edoys-sdoys),outpath+"/parrate/parrate_"+Y+".tif",ref,dtype)
    
    del par, arr, i
    
    ########################
    ### Calculate OP Chl ###
    ########################
    print Y + ": Chl"
    # Load Chl file lists
    flist0 = [f for f in os.listdir(inpath3+"/"+Y) if ( (int(f[8:11]) >= startdoy) and f.startswith('.') == 0 and f.endswith(".tif") )]
    flist1 = [f for f in os.listdir(inpath3+"/"+str(y+1)) if ( (int(f[8:11]) < startdoy) and f.startswith('.') == 0 and f.endswith(".tif") )]
    flist0.sort()
    flist1.sort()
    
    # Load Chl files
    alist = []
    for f in flist0: # These are concentrations, not rates, so raw values are fine
        alist.append( gdalnumeric.LoadFile(inpath3+"/"+Y+"/"+f) )
    for f in flist1:
        alist.append( gdalnumeric.LoadFile(inpath3+"/"+str(y+1)+"/"+f) )
    
    arr = np.array(alist)
    arr = np.where(arr < nanlt, np.nan, arr)
    del alist, flist0, flist1

    # Average Chl
    chl = np.zeros_like(arr[0])*np.nan
    for i in range(len(rows)): # units mg / m^3
        chl[rows[i],cols[i]] = np.nanmean(arr[int(start[rows[i],cols[i]]):int(endnpp[rows[i],cols[i]])+1,rows[i],cols[i]])
    
    md.writeNumpy_gdalObj(chl,outpath+"/chl/chl_"+Y+".tif",ref,dtype)
    
    del chl, arr, i
    
    ########################
    ### Calculate OP SST ###
    ########################
    print Y + ": SST"
    # Load SST file lists
    flist0 = [f for f in os.listdir(inpath4+"/"+Y) if ( (int(f[8:11]) >= startdoy) and f.startswith('.') == 0 and f.endswith(".tif") )]
    flist1 = [f for f in os.listdir(inpath4+"/"+str(y+1)) if ( (int(f[8:11]) < startdoy) and f.startswith('.') == 0 and f.endswith(".tif") )]
    flist0.sort()
    flist1.sort()
    
    # Load SST files
    alist = []
    for f in flist0: # These are concentrations, not rates, so raw values are fine
        alist.append( gdalnumeric.LoadFile(inpath4+"/"+Y+"/"+f) )
    for f in flist1:
        alist.append( gdalnumeric.LoadFile(inpath4+"/"+str(y+1)+"/"+f) )
    
    arr = np.array(alist)
    arr = np.where(arr == nan2, np.nan, arr)
    del alist, flist0, flist1

    # Average SST
    sst = np.zeros_like(arr[0])*np.nan
    for i in range(len(rows)): # units K
        sst[rows[i],cols[i]] = np.nanmean(arr[int(start[rows[i],cols[i]]):int(endnpp[rows[i],cols[i]])+1,rows[i],cols[i]])
    
    md.writeNumpy_gdalObj(sst,outpath+"/sst/sst_"+Y+".tif",ref,dtype)
    
    del sst, arr, i
    
    ###############################
    ### Calculate OP Windy Days ###
    ###############################
    print Y + ": Wind"
    # Load Wind files
    alist = []
    tlist = []
    for m in range(startmonth,12+1) + range(1,startmonth):
        if m < 10:
            M = "0"+str(m)
        else:
            M = str(m)
        if m >= startmonth:
            nc = Dataset(inpath5+"/ubAb10m_EASE25km_"+Y+M+".nc")
        else:
            nc = Dataset(inpath5+"/ubAb10m_EASE25km_"+str(y+1)+M+".nc") 
        alist.append(nc.variables['uvab'][:])
        tlist = tlist + list(nc.variables['time'][:])
        
        nc.close()
    
    arr = np.concatenate(alist,axis=0)
    tarr = np.array(tlist)
    if y == 2003:
        tarr[np.where(tarr == 0)] = np.arange(909024,909054,6)
    elif y == 2007:
        tarr[np.where(tarr == 0)] = np.arange(940440,940464,6)
    
    del alist, M, m, tlist
    
    # Identify high wind events
    arr2 = np.where(arr > windmin, 1, 0)
    
    # Count wind events during period of interest
    hwe = np.zeros_like(arr2[0,:,:])*np.nan
    uvAb = np.zeros_like(arr2[0,:,:])*np.nan
    for i in range(len(rows)):
        # Create indices for determining start and end times for wind
        sdate = md.timeAdd([y,1,0,0,0,0],[0,0,sdoys[rows[i],cols[i]],0,0,0])
        edate = md.timeAdd([y,1,0,0,0,0],[0,0,edoys[rows[i],cols[i]],0,0,0])
        sti = np.where(tarr == md.daysBetweenDates(refdate,sdate)*24)[0][0]
        eti = np.where(tarr == md.daysBetweenDates(refdate,edate)*24)[0][0]
        
        # Sum wind event occurrences, converting to units of "event days"
        hwe[rows[i],cols[i]] = np.sum(arr2[sti:eti,rows[i],cols[i]])/4.
        uvAb[rows[i],cols[i]] = np.mean(arr[sti:eti,rows[i],cols[i]])
  
    md.writeNumpy_gdalObj(hwe,outpath+"/hwe"+str(windmin)+"/hwe"+str(windmin)+"_"+Y+".tif",ref,dtype)
    md.writeNumpy_gdalObj(hwe/(edoys-sdoys),outpath+"/hwf"+str(windmin)+"/hwf"+str(windmin)+"_"+Y+".tif",ref,dtype)
    md.writeNumpy_gdalObj(uvAb,outpath+"/uvAb/uvAb_"+Y+".tif",ref,dtype)
    
    del arr2, hwe, i, arr
    
    ##################################
    ### Calculate OP Cyclone Count ###
    ##################################
    print Y + ": Cyclone Count"
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

    # Load Cyclone Area Fields
    alist = []
    tlist = []
    for m in range(startmonth,12+1) + range(1,startmonth):
        if m >= startmonth:
            flist = os.listdir(inpath6+"/"+Y+"/"+months[m-1])
            flist = [f for f in flist if (f.endswith(".pkl") and f.startswith("CF"))]
            for f in flist:
                try:
                    alist.append( md.unpickle(inpath6+"/"+Y+"/"+f) )
                except:
                    cf = md.unpickle(inpath6+"/"+Y+"/"+months[m-1]+"/"+f)
                    alist.append( interp(np.flipud(cf.fieldAreas),xin,np.flipud(yin),x2out,y2out,order=0) )
                    md.pickle(alist[-1],outpath6+"/"+Y+"/"+f)
                tlist.append([int(f[2:6]),int(f[6:8]),int(f[8:10]),int(f[11:13]),0,0])
        else:
            flist = os.listdir(inpath6+"/"+str(y+1)+"/"+months[m-1])
            flist = [f for f in flist if (f.endswith(".pkl") and f.startswith("CF"))]
            for f in flist:
                try:
                    alist.append( md.unpickle(inpath6+"/"+str(y+1)+"/"+f) )
                except:
                    cf = md.unpickle(inpath6+"/"+str(y+1)+"/"+months[m-1]+"/"+f)
                    alist.append( interp(np.flipud(cf.fieldAreas),xin,np.flipud(yin),x2out,y2out,order=0) )
                    md.pickle(alist[-1],outpath6+"/"+Y+"/"+f)
                tlist.append([int(f[2:6]),int(f[6:8]),int(f[8:10]),int(f[11:13]),0,0]) 
    
    arr = np.array(alist)
    
    del alist, m, flist, cf
    
    # Sum the number of cyclone events during period of interest
    countA = np.zeros_like(arr[0,:,:])*np.nan
    for i in range(len(rows)):
        # Create indices for start and end time
        sdate = md.timeAdd([y,1,0,0,0,0],[0,0,sdoys[rows[i],cols[i]],0,0,0])
        edate = md.timeAdd([y,1,0,0,0,0],[0,0,edoys[rows[i],cols[i]],0,0,0])
        sti = [t for t in range(len(tlist)) if tlist[t] == sdate][0]
        eti = [t for t in range(len(tlist)) if tlist[t] == edate][0]
        
        # Sum cyclone event occurrences, converting to units of "event days"
        countA[rows[i],cols[i]] = np.sum(arr[sti:eti,rows[i],cols[i]])/4.
    
    md.writeNumpy_gdalObj(countA,outpath+"/countA/countA_"+Y+".tif",ref,dtype)
    md.writeNumpy_gdalObj(countA/(edoys-sdoys),outpath+"/caf/caf_"+Y+".tif",ref,dtype)
    
    #########################################################
    ### Calculate Bering Heat Flux (temporally-speaking) ###
    #########################################################
    print Y + ": BHF"
    pdf = pd.read_csv(inpath7)
    
    bvf = np.zeros_like(lats)*np.nan
    for i in range(len(rows)):
        # Create indices for start and end time
        sdate = md.timeAdd([y,1,0,0,0,0],[0,0,sdoys[rows[i],cols[i]],0,0,0])
        edate = md.timeAdd([y,1,0,0,0,0],[0,0,edoys[rows[i],cols[i]],0,0,0])
        
        if sdate[0] == edate[0]:
            bvf[rows[i],cols[i]] = np.mean( pdf['Volume'].loc[(pdf['Year'] == sdate[0]) & (pdf['Month'] >= sdate[1]) & (pdf['Month'] <= edate[1])] )
        else:
            bvf[rows[i],cols[i]] = np.mean( pdf['Volume'].loc[( (pdf['Year'] == sdate[0]) & (pdf['Month'] >= sdate[1])) | ((pdf['Year'] == edate[0]) & (pdf['Month'] <= edate[1]))] )
    
    md.writeNumpy_gdalObj(bvf,outpath+"/bvf/bvf_"+Y+".tif",ref,dtype)
