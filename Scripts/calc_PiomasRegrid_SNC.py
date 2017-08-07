"""
Function regrids Piomas SIT onto selected grid. Also, a test plotter
is provided to check the regridding
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 21 October 2016
    
Usage
-----
    sitnew,lats,lons = regridPiomas(directory,sit,lat1,lon1,lat2,lon2)
"""

print '\n>>> Using regridPiomas function!'

import numpy as np
import read_SnowCover_PIOMAS as CS
from netCDF4 import Dataset
from scipy.interpolate import griddata as g
import matplotlib.colors as c
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import nclcmaps as ncm

### Define directories
directorydata1 = '/surtsey/zlabe/seaice_obs/PIOMAS/' 
directorydata2 = '/surtsey/ypeings/CESM_large_ensemble/SIT/interp_1deg/'

### Alott time series
yearmin = 1979
yearmax = 2017
years = np.arange(yearmin,yearmax+1,1)

### Call functions
lat1,lon1,snc = CS.readPiomas(directorydata1,years,0.0)

### Read new lat/lon grid
files = 'CESM_large_ensemble/SIT/interp_1deg/'
filename = directorydata2 +'b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc'
data = Dataset(filename)
lat2 = data.variables['lat'][:]
lon2 = data.variables['lon'][:]
data.close()

if lat1.ndim == 1:
    if lon1.ndim == 1:
        lon1,lat1 = np.meshgrid(lon1,lat1)
        print 'Made meshgrid of new lats/lons!'

if lat2.ndim == 1:
    if lon2.ndim == 1:
        lon2,lat2 = np.meshgrid(lon2,lat2)
        print 'Made meshgrid of new lats/lons!'

def regrid(lat1,lon1,lat2,lon2,var,years):
    """
    Interpolated on selected grid. Reads PIOMAS in as 4d with 
    [year,month,lat,lon]
    """
    
    varn_re = np.reshape(var,(var.shape[0],var.shape[1],(120*360)))   
    
    varn = np.empty((var.shape[0],var.shape[1],lat2.shape[0],lon2.shape[1]))
    
    print 'Completed: Start regridding process:'
    
    for i in xrange(varn.shape[0]):
        for j in xrange(varn.shape[1]):
            z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,j,:],(lat2,lon2),method='linear')
            varn[i,j,:,:] = z
        print 'Completed: Year %s Regridding---' % (years[i])
    return varn

def netcdfPiomas(lats,lons,var,directory):
    print '\n>>> Using netcdfPIOMAS function!'
    
    name = 'SnowCover/piomas_regrid1x1_snc_LENS_19792017.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'PIOMAS SNC from 1979-2017 ' \
                        'interpolated on selected grid from LENS (1x1deg)'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('months',var.shape[1])
    ncfile.createDimension('lat',var.shape[2])
    ncfile.createDimension('lon',var.shape[3])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    months = ncfile.createVariable('months','f4',('months'))
    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
    varns = ncfile.createVariable('sit','f4',('years','months','lat','lon'))
    
    ### Units
    varns.units = 'meters'
    ncfile.title = 'PIOMAS SNC on LENS Grid 1x1 degree'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'PIOMAS, University of Washington'
    ncfile.references = 'Zhang and Rothrock [2003]'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    months[:] = list(xrange(var.shape[1]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print '*Completed: Created netCDF4 File!'
    
sncn = regrid(lat1,lon1,lat2,lon2,snc,years)
netcdfPiomas(lat2,lon2,sncn,directorydata1)
       
### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)

cmap = ncm.cmap('MPL_PuRd')

sncn[np.where(sncn==0.0)] = np.nan

var = sncn[-2,3,:,:]*100.

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(0,40,10)
values = np.arange(0,30.1,1)

cs = m.contourf(lon2,lat2,var[:,:],
                values,extend='max',latlon=True)
cs1 = m.contour(lon2,lat2,var[:,:],
                barlim,linewidths=0.2,colors='k',
                linestyles='-',latlon=True)                       
        
cs.set_cmap(cmap)

cbar = m.colorbar(cs,location='right',pad='10%',drawedges=False)
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))  
cbar.ax.tick_params(axis='x', size=.1)
cbar.set_label(r'\textbf{snow depth (cm)}')

fig.suptitle(r'\textbf{Testing PIOMAS Regrid}')

plt.savefig('/home/zlabe/Desktop/' + 'testing_regrid_PIOMAS_snownew.png',dpi=300)
'Completed: Script done!'