'''
Author: Alex Crawford
Date Created: 14 Sep 2018
Date Modified: 28 Nov 2018

Reads csv files designed by gridcell and aggregates them to regions.
Only aggregates for 8-day periods that have valid data for npp.
'''
'''############
Set up Modules
############'''
import pandas as pd
from osgeo import gdal, gdalnumeric
import numpy as np
import scipy.stats as stats
import MERRA_Module as md

'''###########
Declare Variables
###########'''
r = "3" # "2" or "3"
v = "4" # "", "2", "3", or "4"
nv = -99

path = "/Volumes/Prospero"
inpath = path+"/MODIS/VGPM/VGPM_PostBloom8/loccsv"
outpath = path+"/MODIS/VGPM/VGPM_PostBloom8/regcsv"+r+"_sw"

maskN = path+"/Boundaries/EASE2_SSMI_Regions"+r+"_N60.tif"

regs = [3,8,9,10,11,12,13,16,18,19] #[11,12] #  [3,8,9,10,11,12,13,16] #

if r == "2":
    REGS = ["","","","Bering Sea","","","","","Barents Sea","Kara Sea","Laptev Sea",\
    "East Siberian Sea","Chukchi Sea","Beaufort Sea","","","Baffin Bay"]
elif r == "3":
    REGS = ["","","","Bering Sea","","","","","Barents Sea","Kara Sea","Laptev Sea",\
    "E East Siberian Sea","N Chukchi Sea","Beaufort Sea","","","Baffin Bay",\
    "", "W East Siberian Sea","S Chukchi Sea"]

'''################
Main Analysis
################'''

mask = gdalnumeric.LoadFile(maskN)

for reg in regs:
    print REGS[reg] + ": Loading"

    rows, cols = np.where(mask == reg)
    
    ### LOAD FILES IN THIS REGION ###
    pdf = pd.DataFrame(columns=['lat','lon','year','doy','opendays',\
                    'npp','chl','par','sst','hwe','hu','hv','uv','u','v'])
    rdf = pd.DataFrame(columns=['year','doy','opendays','area',\
                    'npp','chl','par','sst','hwe','hu','hv','uv','u','v'])

    for i in range(len(rows)):
        try:
            gdf = pd.read_csv(inpath+"_sw"+v+"/PostBloom_8Day_Values_Row"+str(rows[i])+"_Col"+str(cols[i])+".csv")
        except:
            gdf = pd.DataFrame(columns=['lat','lon','year','doy','opendays',\
                    'npp','chl','par','sst','hwe','hu','hv','uv','u','v'])
        try:
            hdf = pd.read_csv(inpath+v+"/PostBloom_8Day_Values_Row"+str(rows[i])+"_Col"+str(cols[i])+".csv")
        except:
            hdf = pd.DataFrame(columns=['lat','lon','year','doy','opendays',\
                    'npp','chl','par','sst','hwe','hu','hv','uv','u','v'])
        
        gdf = gdf.append(hdf,ignore_index=True)
        
        if len(gdf) > 0:
            pdf = pdf.append(gdf,ignore_index=True)
    
    ### AGGREGATE GRIDCELLS FOR EACH 8-DAY PERIOD ###
    print REGS[reg] + ": Aggregating"
    for co in pdf.columns: # Declare NaNs
        pdf.loc[pdf[co] == nv,co] = np.nan
    
    pdf = pdf.loc[np.isfinite(pdf['npp']) == 1] # Limit to valid npp weeks
    
    for y in np.unique(pdf['year']):
        ydf = pdf.loc[pdf['year'] == y]
        
        for w in np.unique(ydf['doy']):
            wdf = ydf.loc[ydf['doy'] == w]
            newrow = pd.DataFrame([{'year':y,'doy':w,'opendays':np.nanmean(list(wdf['opendays'])),\
                            'area':len(wdf),'npp':np.sum(list(wdf['npp'])),'chl':np.nanmean(list(wdf['chl'])),\
                            'par':np.nanmean(list(wdf['par'])),'sst':np.nanmean(list(wdf['sst'])),\
                            'hwe':np.mean(list(wdf['hwe'])),'uv':np.mean(list(wdf['uv'])),\
                            'u':np.mean(list(wdf['u'])),'v':np.mean(list(wdf['v'])),\
                            'hu':np.nanmean(list(wdf['hu'])),'hv':np.nanmean(list(wdf['hv'])),\
                            },])
            rdf = rdf.append(newrow)
            del newrow, wdf
        
        del ydf
    
    # Write to File
#    rdf.to_csv(outpath+"/PostBloom_8Day_"+REGS[reg]+".csv",index=False)
