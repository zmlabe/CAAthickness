"""
Scripts reads in sea ice thickness data from CS-2 corrected/uncorrected
and plots a comparison over 2011-2017
 
Notes
-----
    Author : Zachary Labe
    Date   : 15 January 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as c
import datetime
import statsmodels.api as sm
from mpl_toolkits.basemap import Basemap
import statsmodels.api as sm
import scipy.stats as sts
import nclcmaps as ncm
import pandas as pd
import cmocean

### Define directories
directorydata = '/surtsey/zlabe/seaice_obs/cs2/'
directoryfigure = '/home/zlabe/Desktop/CS2PIOMAS/Corrected/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calc sea ice thickness data sets - %s----\n' % titletime) 

### Alott time series
yearmin = 2011
yearmax = 2017
years = np.arange(yearmin,yearmax+1,1)
months = [r'April']

### Read in data
lonvals = []
latvals = []
sitc = []
situ = []
for i in range(years.shape[0]):
    filename = 'corrected_FYI_SIT.xlsx'
    dataq = pd.read_excel(directorydata + filename,sheetname='%s' % years[i],
                         header=0,parse_cols=range(0,4))
    data = dataq.as_matrix()
    latq = data[:,0] # latitudes
    lonq = data[:,1] # longitudes
    situq = data[:,2] # uncorrected CS-2 FYI thickness
    sitcq = data[:,3] # corrected CS-2 FYI thickness
    
    latvals.append(latq)
    lonvals.append(lonq)
    situ.append(situq)
    sitc.append(sitcq)

    print('Completed: Read in %s data!' % years[i])  
    
### Appends lat/lon
lonvalsn = np.append(lonvals,lonvals)
latvalsn = np.append(latvals,latvals)

############################################################################
############################################################################
############################################################################
### Plot comparison between corrected and uncorrected SIT in CS2
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']})     
    
fig = plt.figure(figsize=(10,4))
for i in range(14):
    ax = plt.subplot(2,7,i+1)
    
    x1 = lonvalsn[i]
    y1 = latvalsn[i]
    
    if i < 7:
        comp = situ[i]
    else:
        comp = sitc[i-7]
    
    style = 'polar'
    if style == 'ortho':
        m = Basemap(projection='ortho',lon_0=-90,
                    lat_0=70,resolution='l',round=True)
    elif style == 'polar':
        m = Basemap(projection='npstere',boundinglat=56,lon_0=270,
                    resolution='l',round =True)
        
    m.drawmapboundary(fill_color='darkgrey',linewidth=0.7)
    m.drawcoastlines(color='darkgrey',linewidth=0.2)
    m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')
    
    #### Plot filled contours 
    x,y = m(x1,y1)
    comp[np.isnan(comp)]=0.0
    cs = m.hexbin(x,y,C=comp,vmin = 0,vmax = 3) 
    
    barlim = np.arange(0,4,3) 
    cs.set_cmap(cmocean.cm.thermal)
    m.fillcontinents(color='dimgrey')
    
cbar_ax = fig.add_axes([0.313,0.12,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)                      
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim)))    
cbar.set_label(r'\textbf{Thickness (m)}',color='k',labelpad=0.1)
cbar.ax.tick_params(axis='x', size=.001)
cbar.outline.set_edgecolor('dimgrey')

ax.annotate(r'\textbf{CS2$_u$}',
        xy=(0, 0),xytext=(0.09,0.7),xycoords='figure fraction',
        fontsize=25,color='k',rotation=90)

ax.annotate(r'\textbf{CS2$_c$}',
    xy=(0, 0),xytext=(0.09,0.363),xycoords='figure fraction',
    fontsize=25,color='k',rotation=90)
    
fig.subplots_adjust(hspace=-0.2)
fig.subplots_adjust(wspace=0)    

plt.savefig(directoryfigure + 'CS2_correctUncorrect_compare.png',dpi=300)
