'''*********************************************
Authors: Alex Crawford
Date Created: 9/6/19
Date Modified: 
Purpose: To calculate the day of the year on which sea ice retreats and
advances in each grid cell of a particular sector each year.

Inputs: 
    concentration threshold (cts) -- value between 0 and 1
    years of interest (ymin, ymax) -- integers
    months for  -- the month on which to start the annual cycle -- March (3) is
        a good idea because it is the closest to the maximum
    
Outputs: A csv file with the concentration threshold and setor noted in the 
    file name. Retreat and advance days are recorded as a "DOY", with 1 being
    the 1st day of January. If retreat below the
    concentration threshold never occurs, the minimum day is recorded instead.
*********************************************'''
# Import clock:
from time import perf_counter as clock
# Start script stopwatch. The clock starts running when time is imported
start = clock()

'''*******************************************
Set up Modules
*******************************************'''
print("Importing Modules")

import os
import netCDF4 as nc
import numpy as np
def leapyearBoolean(years):
    '''
    Given a list of years, this function will identify which years are leap 
    years and which years are not. Returns a list of 0s (not a leap year) and 
    1s (leap year). Each member of the year list must be an integer or float.
    
    Requires numpy.
    '''
    ly = [] # Create empty list
    for y in years: # For each year...
        if (y%4 == 0) and (y%100 != 0): # If divisible by 4 but not 100...
            ly.append(1) # ...it's a leap year
        elif y%400 == 0: # If divisible by 400...
            ly.append(1) # ...it's a leap year
        else: # Otherwise...
            ly.append(0) # ...it's NOT a leap year
    
    return ly

def daysBetweenDates(date1,date2,lys=1):
    '''
    Calculates the number of days between date1 (inclusive) and date2 (exclusive)
    when given dates in list format [year,month,day,hour,minute,second] or 
    [year,month,day]. Works even if one year is BC (BCE) and the other is AD (CE). 
    If hours are used, they must be 0 to 24. Requires numpy.
    
    date1 = the start date (earlier in time; entire day included in the count if time of day not specified)\n
    date2 = the end date (later in time; none of day included in count unless time of day is specified)
    '''
    db4 = [0,31,59,90,120,151,181,212,243,273,304,334] # Number of days in the prior months
    
    if date1[0] == date2[0]: # If the years are the same...
        # 1) No intervening years, so ignore the year value:        
        daysY = 0
    
    else: # B) If the years are different...
        
        # 1) Calculate the total number of days based on the years given:
        years = range(date1[0],date2[0]) # make a list of all years to count
        years = [yr for yr in years if yr != 0]
        if lys==1:
            lyb = leapyearBoolean(years) # annual boolean for leap year or not leap year
        else:
            lyb = [0]
        
        daysY = 365*len(years)+np.sum(lyb) # calculate number of days
    
    if lys == 1:
        ly1 = leapyearBoolean([date1[0]])[0]
        ly2 = leapyearBoolean([date2[0]])[0]
    else:
        ly1, ly2 = 0, 0
    
    # 2) Calcuate the total number of days to subtract from start year
    days1 = db4[date1[1]-1] + date1[2] -1 # days in prior months + prior days in current month - the day you're starting on
    # Add leap day if appropriate:    
    if date1[1] > 2:
        days1 = days1 + ly1
    
    # 3) Calculate the total number of days to add from end year
    days2 = db4[date2[1]-1] + date2[2] - 1 # days in prior months + prior days in current month - the day you're ending on
    # Add leap day if appropriate:    
    if date2[1] > 2:
        days2 = days2 + ly2
        
    # 4) Calculate fractional days (hours, minutes, seconds)
    day1frac, day2frac = 0, 0
    
    if len(date1) == 6:
        day1frac = (date1[5] + date1[4]*60 + date1[3]*3600)/86400.
    elif len(date1) != 3:
        raise Exception("date1 does not have the correct number of values.")
    
    if len(date2) == 6:
        day2frac = (date2[5] + date2[4]*60 + date2[3]*3600)/86400.
    elif len(date2) != 3:
        raise Exception("date2 does not have the correct number of values.")
    
    # 5) Final calculation
    days = daysY - days1 + days2 - day1frac + day2frac
    
    return days
    
'''*******************************************
Declare Variables
*******************************************'''
print("Declaring Variables")

### Input Variables ###
nrow = 2325 # Number of rows
ncol = 2014 # Number of columns
ct = 0.75 # A number between 0 and 1 (or 0 and 100) for the concentration threshold
n = 5 # Moving Average Size (n = # observation on either side of current day;
        ## so 1 = 3-point, 2 = 5-point, 4 = 9-point)

### Time Variables ###
maxmo = [1,4] # months in which the sea ice maximum may occur
minmo = [8,10] # months in which the sea ice minimum may occur

### Path Variables ###
path = "/Volumes/Prospero/Arrigo/"
inpath = path+"/ice/SmoothedMA"+str(n)
outpath = path+"/ice/AdvanceRetreat"

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")
# Time Set Up
mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
days = ["01","02","03","04","05","06","07","08","09","10","11","12","13",\
    "14","15","16","17","18","19","20","21","22","23","24","25","26","27",\
    "28","29","30","31"]
    
files = os.listdir(inpath)
files = [f for f in files if f.startswith('.') == 0]
files.sort()

CT = str(int(ct*100))
lats = np.fromfile(path+"/Projections/lats_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
lons = np.fromfile(path+"/Projections/lons_arctic_50N_2014x2325_4f.flat",dtype='f').reshape((nrow,ncol))
    
lrdL, frdL, ladL, fadL, minL, mndL, opcL, opL = [], [], [], [], [], [], [], []

for fi in range(len(files)): # Loop through each year
    print("Starting " + files[fi])
        
    # Load indices for dates
    y = int(files[fi][33:37]) # The Current Year
    
    maxi0 = daysBetweenDates([y,1,1],[y,maxmo[0],1], 1)
    maxi1 = daysBetweenDates([y,1,1],[y,maxmo[1]+1,1], 1)
    
    mini0 = daysBetweenDates([y,1,1],[y,minmo[0],1], 1)
    mini1 = daysBetweenDates([y,1,1],[y,minmo[1]+1,1], 1)
    
    maxi2 = daysBetweenDates([y,1,1],[y+1,maxmo[0],1], 1)
    maxi3 = daysBetweenDates([y,1,1],[y+1,maxmo[1]+1,1], 1)
    
    # Current Year
    ncf = nc.Dataset(inpath+"/"+files[fi])
    arr = ncf.variables['siconc'][:]
    
    # Next Year
    try:
        ncf2 = nc.Dataset(inpath+"/"+files[fi+1])
        arr2 = ncf2.variables['siconc'][:((maxi3-arr.shape[0])+1),:,:]
        arr = np.concatenate((arr, arr2), 0)
        del arr2
    except:
        arr = np.concatenate((arr, arr[:((maxi3-arr.shape[0])+1),:,:]), 0)
        print(" -- Note, this is the last year, so any result greater than 365 is invalid")
    
    # Calculate minimum & maximum value for year
    Maxes1 = np.amax(arr[maxi0:maxi1],0)
    Mins = np.amin(arr[mini0:mini1],0)
    Maxes2 = np.amax(arr[maxi2:maxi3],0)
    
    ### Prep Outputs ###
    mndArr, mxdArr = np.zeros_like(Mins)*np.nan, np.zeros_like(Mins)*np.nan
    frdArr, lrdArr = np.zeros_like(Mins)*np.nan, np.zeros_like(Mins)*np.nan
    fadArr, ladArr = np.zeros_like(Mins)*np.nan, np.zeros_like(Mins)*np.nan
    opArr, opcArr = np.zeros_like(Mins)*np.nan, np.zeros_like(Mins)*np.nan
    
    ### Calculate Retreat and Advance Events ###
    validcells = np.where( (np.isnan(Mins) == 0) & (np.isnan(Maxes1) == 0) & (np.isnan(Maxes2) == 0) )
    for i in range(len(validcells[0])):
        r,c = validcells[0][i], validcells[1][i] # Assign row and column
        
        # Calculate index for minimum & maximum
        MaxesI1 = int(np.median(np.where(arr[maxi0:maxi1,r,c] == Maxes1[r,c]))) # Gives first occurrence of maximum if multiples present
        MinsI = mini0 + int(np.median(np.where(arr[mini0:mini1,r,c] == Mins[r,c]))) # Gives first occurrence of minimum if multiples present
        MaxesI2 = maxi2 + int(np.median(np.where(arr[maxi2:maxi3,r,c] == Maxes2[r,c]))) # Gives first occurrence of maximum if multiples present
        
        # Store Minimum Day
        mxdArr[r,c] = MaxesI1 + 1
        mndArr[r,c] = MinsI + 1
        mxd2rr = MaxesI2 + 1
        
        # If it's always above the concentration threshold...
        if Mins[r,c] >= ct: 
            opArr[r,c], opcArr[r,c] = 0, 0
            
        # If it's never above the concentration threshold... 
        elif (Maxes1[r,c] < ct) & (Maxes2[r,c] < ct):
            opArr[r,c], opcArr[r,c] = 365, 365
            
        # Otherwise...
        else:
            above = np.where(arr[:,r,c] >= ct)[0] # Indices above concentration
            below = np.where(arr[:,r,c] < ct)[0] # Indices below concentration
            
            # First Retreat Day
            # First index after Maxes1 and before/on Mins for which concentration is below threshold
            try:
                frdArr[r,c] = below[np.where((below <= MinsI) & (below > MaxesI1))][0] + 1
            except:
                frdArr[r,c] = np.nan
            
            # Last Retreat Day
            # Last index after Maxes1 and before/on Mins for which concentration is below threshold
            try:
                lrdArr[r,c] = above[np.where((above < MinsI) & (above >= MaxesI1))][-1] + 1
            except:
                lrdArr[r,c] = np.nan

            # First Advance Day
            # First index after Mins and before/on Maxes2 for which concentration is above threshold
            try: 
                fadArr[r,c] = above[np.where((above > MinsI) & (above <= MaxesI2))][0] + 1
            except:
                fadArr[r,c] = np.nan
            
            # Last Advance Day
            # Last index after Mins anbd before/on Maxes2 for which concentration is below threshold
            try:
                ladArr[r,c] = below[np.where((below >= MinsI) & (below < MaxesI2))][-1] + 1
            except:
                ladArr[r,c] = np.nan
            
            # Open Water Periods
            if (Maxes1[r,c] < ct)  & (Maxes2[r,c] >= ct): # When it starts below threshold but ends above
                opArr[r,c] =  ladArr[r,c] - mxdArr[r,c] #np.min([365, ladArr[r,c] - MaxesI1])
                opcArr[r,c] = fadArr[r,c] - mxdArr[r,c] #np.min([365, fadArr[r,c] - MaxesI1])
                
                lrdArr[r,c] = np.nan
                frdArr[r,c] = np.nan
                
            elif (Maxes1[r,c] >= ct)  & (Maxes2[r,c] < ct): # When it starts above threshold but ends below
                opArr[r,c] =  mxd2rr - frdArr[r,c] #np.min([365, MaxesI2 - frdArr[r,c]])
                opcArr[r,c] = mxd2rr - lrdArr[r,c] #np.min([365, MaxesI2 - lrdArr[r,c]])
                
                fadArr[r,c] = np.nan
                ladArr[r,c] = np.nan
                
            else: # Simple Case
                opArr[r,c] =  ladArr[r,c] - frdArr[r,c] #np.min([365, ladArr[r,c] - frdArr[r,c]])
                opcArr[r,c] =  fadArr[r,c] - lrdArr[r,c] #np.min([365, fadArr[r,c] - lrdArr[r,c]])
    
    # Append to lists
    minL.append(Mins)
    mndL.append(mndArr)
    lrdL.append(lrdArr)
    ladL.append(ladArr)
    frdL.append(frdArr)
    fadL.append(fadArr)
    opcL.append(opcArr)
    opL.append(opArr)
    
    # Remove objects from memory
    del frdArr, lrdArr, fadArr, ladArr, opArr, opcArr
    del MinsI, validcells, above, below, MaxesI1, MaxesI2, mxdArr, mxd2rr, mndArr

### Write Outputs ###
outName = "siphenology_C"+CT+"_"+str(int(files[0][33:37]))+"-"+str(int(files[-1][33:37]))+".nc"

ncf1 = nc.Dataset(outpath+"/"+outName, 'w', format='NETCDF4')
ncf1.createDimension('y', Mins.shape[0])
ncf1.createDimension('x', Mins.shape[1])
ncf1.createDimension('time', int(files[-1][33:37])-int(files[0][33:37])+1 )

yNC = ncf1.createVariable('y', np.float32, ('y',))
xNC = ncf1.createVariable('x', np.float32, ('x',))
tNC = ncf1.createVariable('time', np.float32, ('time',))

minNC = ncf1.createVariable('min', np.float32, ('time','y','x',))
mndNC = ncf1.createVariable('mnd', np.float32, ('time','y','x',))
lrdNC = ncf1.createVariable('lrd', np.float32, ('time','y','x',))
ladNC = ncf1.createVariable('lad', np.float32, ('time','y','x',))
frdNC = ncf1.createVariable('frd', np.float32, ('time','y','x',))
fadNC = ncf1.createVariable('fad', np.float32, ('time','y','x',))
opcNC = ncf1.createVariable('opc', np.float32, ('time','y','x',))
opNC = ncf1.createVariable('op', np.float32, ('time','y','x',))

latNC = ncf1.createVariable('lat', np.float32, ('y','x',))
lonNC = ncf1.createVariable('lon', np.float32, ('y','x',))

ncf1.description = 'Sea Ice Concentration with 5-Day Moving Average'
ncf1.source = 'netCDF4 python module'
tNC.units = 'years' #e.g., 'days since 1979-01-01 00:00:00.0'
latNC.units = 'degrees north'
lonNC.units = 'degrees east'      
minNC.units = 'percentage'
mndNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
lrdNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
frdNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
ladNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
fadNC.units = 'day of year (1st of Jan = 1, 1st of Feb = 32)'
opcNC.units = 'days'
opNC.units = 'days'
   

tNC[:] = np.arange(int(files[0][33:37]),int(files[-1][33:37])+1,1)
latNC[:] = ncf.variables['lat'][:]
lonNC[:] = ncf.variables['lon'][:]
minNC[:] = np.array(minL)
mndNC[:] = np.array(mndL)
lrdNC[:] = np.array(lrdL)
ladNC[:] = np.array(ladL)
frdNC[:] = np.array(frdL)
fadNC[:] = np.array(fadL)
opcNC[:] = np.array(opcL)
opNC[:] = np.array(opL)

ncf1.close()
ncf.close()

now = clock()
print(" -- Completed, Elapsed Time: " + str(now-start))
    
print("Complete.")
