"""
Scripts reads in sea ice thickness data from CS-2 and PIOMAS for 
exploratory data analysis
 
Notes
-----
    Author : Zachary Labe
    Date   : 22 July 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_SeaIceThick_PIOMAS as CT
import read_PiomasArea as CA
import statsmodels.api as sm
from mpl_toolkits.basemap import Basemap
from openpyxl import load_workbook

### Define directories
directorydata1 = '/surtsey/zlabe/seaice_obs/PIOMAS/'  
directorydata2 = '/surtsey/zlabe/seaice_obs/cs2/' 
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceVariability/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calc sea ice thickness data sets - %s----' % titletime 

### Alott time series
yearmin = 2010
yearmax = 2017
years = np.arange(yearmin,yearmax+1,1)
months = [r'April']

### Call functions
lats,lons,sitp = CT.readPiomas(directorydata1,years,0.15)
area = CA.readPiomasArea(directorydata1)

### Read in other data - column 2 is warren snow
wb = load_workbook(directorydata2 + 'cryosat_FYI.xlsx',read_only=True)
ws = wb.get_sheet_by_name('cryosat')

# latitudes
lats = np.array([r[0].value for r in ws.iter_rows()])
latsq = np.asarray(map(lambda x: float(x),lats[1:]))

# longitudes
lons = np.array([r[1].value for r in ws.iter_rows()])
lonsq = np.asarray(map(lambda x: float(x),lons[1:]))

# sea ice thickness from cryosat 
sitc = np.array([r[3].value for r in ws.iter_rows()])
sitcq = np.asarray(map(lambda x: float(x),sitc[1:]))

### Retrieve April SIT
sitpa = sitp[:,4,:,:]

###########################################################################
###########################################################################
###########################################################################
### Calculating temporal sit
def weightThick(var,area):
    """
    Area weights sit array 4d [year,month,lat,lon] into [year,month]
    """
    sityr = np.empty((var.shape[0],var.shape[1]))
    for i in xrange(var.shape[0]):
        for j in xrange(var.shape[1]):
            varq = var[i,j,:,:]
            mask = np.isfinite(varq) & np.isfinite(area)
            varmask = varq[mask]
            areamask = area[mask]
            sityr[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted SIT average!' 
    return sityr
     
#sitave = weightThick(sit,area)
#sitave = weightThick(sit,area)

    
