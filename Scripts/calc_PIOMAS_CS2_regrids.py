"""
Scripts reads in sea ice thickness data from CS-2 and PIOMAS for 
exploratory data analysis
 
Notes
-----
    Author : Zachary Labe
    Date   : 8 August 2017
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import datetime
from scipy.interpolate import griddata as g

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
print '\n' '----Interpolate PIOMAS/CS-2 data sets - %s----' % titletime 

### Alott time series
yearmin = 2011
yearmax = 2017
years = np.arange(yearmin,yearmax+1,1)
yearsp = np.arange(1979,2017+1,1)
months = [r'April']

### Read in CryoSat-2 (CS2) data
lonvals = []
latvals = []
snowvals = []
sitvals = []
for i in xrange(years.shape[0]):
    ### Call data files
    datacs = np.genfromtxt(directorydata2 + 'cs2_04_%s.txt' % years[i],
                           unpack=False,usecols=[0,1,2,3,4],skip_header=1)

    ### Read CS2 data in separate arrays
    latcs = datacs[:,1] # latitudes
    latvals.append(latcs)
    
    loncs = datacs[:,2] # longitudes
    ### Change CS2 longitude coordinates for consistency with PIOMAS
    lon36 = np.where(loncs <= 0.)[0]
    loncs[lon36] = loncs[lon36] + 360.
    lonvals.append(loncs)
    
    snccs = datacs[:,3] # Warren snow depth
    snowvals.append(snccs)
    
    sitcs = datacs[:,4] # CS2 sea ice thickness  
    sitvals.append(sitcs)
    
    print '\nCompleted: Read in CS2 data for April %s!' % years[i]
              
def piomasReader(directory,types,years):
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
    print '\n>>> Using piomasReader function!'
    
    if types == 'thickness':
        filename = 'Thickness/piomas_regrid1x1_sit_LENS_19792017.nc'
        data = Dataset(directory + filename)
        lat = data.variables['lat'][:,:]
        lon = data.variables['lon'][:,:]
        var = data.variables['sit'][:,:,:,:]
        data.close()
        
    elif types == 'snow':
        filename = 'SnowCover/piomas_regrid1x1_snc_LENS_19792017.nc'
        data = Dataset(directory + filename)
        lat = data.variables['lat'][:,:]
        lon = data.variables['lon'][:,:]
        var = data.variables['sit'][:,:,:,:]
        data.close()
            
    var[np.where(var == 0.)] = np.nan
    
    print 'Completed: PIOMAS data read!'
    return lat,lon,var

### Call functions to read in PIOMAS 1x1 degree data    
lat1,lon1,sit1 = piomasReader(directorydata1,'thickness',years)
lat1,lon1,snc1 = piomasReader(directorydata1,'snow',years)

### Select month for data analysis (April)
sitpa = sit1[:,3,:,:]
sncpa = snc1[:,3,:,:]

### Select years for data analysis 
yearq = np.where(yearsp >= 2011)[0]
sitpa = sitpa[yearq,:,:]
sncpa = sncpa[yearq,:,:]

#### Find PIOMAS values of CS2 lat/lons    
interpsit = []
interpsnc = []
for i in xrange(years.shape[0]):    
    psit = g((lon1.flatten(),lat1.flatten()),
             sitpa[i,:,:].flatten(),(lonvals[i],latvals[i]))  
    psnc = g((lon1.flatten(),lat1.flatten()),
             sncpa[i,:,:].flatten(),(lonvals[i],latvals[i]))
             
    interpsit.append(psit)
    interpsnc.append(psnc)
    
    print '\nCompleted: Interpolated year %s!' % years[i]

###########################################################################
###########################################################################
###########################################################################
### Save interpolated PIOMAS data
directorytext = '/home/zlabe/Documents/Projects/CAAthickness/Data/'

for i in xrange(years.shape[0]):
    np.savetxt(directorytext + 'PIOMAS_sit_04_%s.txt' % years[i],
               np.c_[latvals[i],lonvals[i],interpsit[i],interpsnc[i],
                sitvals[i],snowvals[i]],fmt='%1.3f',
                header='### Lats,Lons,PIOMAS SIT,PIOMAS Snow Depth,' \
                'CS2 SIT,Warren Snow Depth',
                comments='################## FYI, April %s values\n' % years[i])
    print 'Completed: Created text file of interpolation year %s!' % years[i]

print 'Completed: Script done!'