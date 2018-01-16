"""
Scripts reads in sea ice thickness data from CS-2 and PIOMAS for 
exploratory data analysis
 
Notes
-----
    Author : Zachary Labe
    Date   : 15 January 2018
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import datetime
from scipy.interpolate import griddata as g
import pandas as pd

### Define directories
directorydata1 = '/surtsey/zlabe/seaice_obs/PIOMAS/'  
directorydata2 = '/surtsey/zlabe/seaice_obs/cs2/' 
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Projects/CAAthickness/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Interpolate PIOMAS/CS-2 data sets - %s----' % titletime) 

### Alott time series
yearmin = 2011
yearmax = 2017
years = np.arange(yearmin,yearmax+1,1)
yearsp = np.arange(1979,2017+1,1)
months = [r'April']

### Read in data
lonvals = []
latvals = []
sitc = []
situ = []
for i in range(years.shape[0]):
    filename = 'corrected_FYI_SIT.xlsx'
    dataq = pd.read_excel(directorydata2 + filename,sheetname='%s' % years[i],
                         header=0,parse_cols=range(0,4))
    data = dataq.as_matrix()
    latq = data[:,0] # latitudes
    lonq = data[:,1] # longitudes
    situq = data[:,2] # uncorrected CS-2 FYI thickness
    sitcq = data[:,3] # corrected CS-2 FYI thickness
    
    latvals.append(latq)
    ### Change CS2 longitude coordinates for consistency with PIOMAS
    lon36 = np.where(lonq <= 0.)[0]
    lonq[lon36] = lonq[lon36] + 360.
    lonvals.append(lonq)
    
    situ.append(situq)
    sitc.append(sitcq)

    print('Completed: Read in %s data!' % years[i])  
    
### Appends lat/lon
lonvalsn = np.append(lonvals,lonvals)
latvalsn = np.append(latvals,latvals)
              
def piomasReader(directory,types,years,actual):
    """
    Function reads PIOMAS binary and converts to standard numpy array.

    Parameters
    ----------
    directory : string
        working directory for stored PIOMAS files
    types : string
        thickness or snow
    years : integers
        years for data files
    actual : True or False
        actual thickness or heff

    Returns
    -------
    lats : 2d array
        latitudes
    lons : 2d array
        longitudes
    var : 4d array [year,month,lat,lon]
        sea ice thickness (m) or snow cover (m)

    Usage
    -----
    lats,lons,var = piomasReader(directory,types,years)
    """    
    print('\n>>> Using piomasReader function!')
    
    if types == 'thickness':
        filename = 'Thickness/piomas_regrid1x1_sit_LENS_19792017.nc'
        data = Dataset(directory + filename)
        lat = data.variables['lat'][:,:]
        lon = data.variables['lon'][:,:]
        var = data.variables['sit'][:,:,:,:]
        data.close()
        
        if actual == True:
            filename = 'SeaIceConcentration/piomas_regrid1x1_sic_LENS_19792017.nc'
            data = Dataset(directory + filename)
            lat = data.variables['lat'][:,:]
            lon = data.variables['lon'][:,:]
            sic = data.variables['sit'][:,:,:,:]
            data.close()
            
            var = var/sic
            print('Completed: SIT/SIC!')
        
    elif types == 'snow':
        filename = 'SnowCover/piomas_regrid1x1_snc_LENS_19792017.nc'
        data = Dataset(directory + filename)
        lat = data.variables['lat'][:,:]
        lon = data.variables['lon'][:,:]
        var = data.variables['sit'][:,:,:,:]
        data.close()
            
    var[np.where(var == 0.)] = np.nan
    
    print('Completed: PIOMAS data read!')
    return lat,lon,var

### Call functions to read in PIOMAS 1x1 degree data    
lat1,lon1,sit1 = piomasReader(directorydata1,'thickness',years,True)
lat1,lon1,snc1 = piomasReader(directorydata1,'snow',years,False)

### Select month for data analysis (April)
sitpa = sit1[:,3,:,:]
sncpa = snc1[:,3,:,:]

### Select years for data analysis 
yearq = np.where(yearsp >= 2011)[0]
sitpa = sitpa[yearq,:,:]
sncpa = sncpa[yearq,:,:]

#### Find PIOMAS values of CS2 lat/lons    
interpsit = []
interpsnow = []
for i in range(years.shape[0]):    
    psit = g((lon1.flatten(),lat1.flatten()),
             sitpa[i,:,:].flatten(),(lonvals[i],latvals[i]))  
    psnc = g((lon1.flatten(),lat1.flatten()),
         sncpa[i,:,:].flatten(),(lonvals[i],latvals[i]))  
             
    interpsit.append(psit)
    interpsnow.append(psnc)
    
    print('\nCompleted: Interpolated year %s!' % years[i])

###########################################################################
###########################################################################
###########################################################################
### Save interpolated PIOMAS data
directorytext = '/home/zlabe/Documents/Projects/CAAthickness/Data/'

for i in range(years.shape[0]):
    np.savetxt(directorytext + 'PIOMAS_CS2UC_sit_04_%s.txt' % years[i],
               np.c_[latvals[i],lonvals[i],interpsit[i],interpsnow[i],
                situ[i],sitc[i]],fmt='%1.3f',
                header='### Lats,Lons,PIOMAS SIT,PIOMAS Snow,CS2U SIT,CS2C SIT',
                comments='################## FYI, April %s values\n' % years[i])
    print('Completed: Created text file of interpolation year %s!' % years[i])

print('Completed: Script done!')