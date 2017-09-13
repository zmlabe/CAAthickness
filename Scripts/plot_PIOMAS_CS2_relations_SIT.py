"""
Scripts reads in sea ice thickness data from CS-2 and interpolated PIOMAS 
for exploratory data analysis
 
Notes
-----
    Author : Zachary Labe
    Date   : 11 August 2017
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

### Define directories
directorydata = '/home/zlabe/Documents/Projects/CAAthickness/Data/'
directoryfigure = '/home/zlabe/Desktop/CS2PIOMAS/Thickness/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calc sea ice thickness data sets - %s----\n' % titletime 

### Alott time series
yearmin = 2011
yearmax = 2017
years = np.arange(yearmin,yearmax+1,1)
months = [r'April']

### Read in data
lonvals = []
latvals = []
sitp= []
snowp = []
sitc = []
snowc = []
for i in xrange(years.shape[0]):
    filename = 'PIOMAS_sit_04_%s.txt' % years[i]
    data = np.genfromtxt(directorydata + filename,unpack=False,
                         usecols=[0,1,2,3,4,5],skip_header=2)
    latq = data[:,0]
    lonq = data[:,1]
    sitpq = data[:,2]
    snowpq = data[:,3]
    sitcq = data[:,4]
    snowcq = data[:,5] 
    
    lonvals.append(lonq)
    latvals.append(latq)
    sitp.append(sitpq)
    snowp.append(snowpq)
    sitc.append(sitcq)
    snowc.append(snowcq)

    print 'Completed: Read in %s data!' % years[i]   

### Calculate linear relationship
def calcStats(varx,vary,years):
    timex = np.arange(varx.shape[0])
    timey = np.arange(vary.shape[0])
    
    mask = np.isfinite(varx) & np.isfinite(vary)
    slope, intercept, r_value, p_value, std_err = \
                                    sts.linregress(varx[mask],vary[mask])
    line = slope*timex + intercept

    print 'Completed: Calculate statistics between data sets!'
    return timex,timey,line,slope,r_value

###########################################################################
###########################################################################
###########################################################################
### Assess statistics for sea ice thickness 
timexi = []
timeyi = []
linei = []
slopei = []
r_valuei = []
for i in xrange(years.shape[0]):    
    timexiq,timeyiq,lineiq,slopeiq,r_valueiq = calcStats(sitp[i],
                                                         sitc[i],years)    
    timexi.append(timexiq)
    timeyi.append(timeyiq)
    linei.append(lineiq)
    slopei.append(slopeiq)
    r_valuei.append(r_valueiq)
    
### Plot figures
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 0))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 
        
###########################################################################
### Plot absolute differences
limit = np.arange(-1,1.1,0.1)
barlim = np.arange(-1,2,1)        
        
fig = plt.figure()
for i in xrange(years.shape[0]):
    diffsit = sitc[i] - sitp[i]
    diffsit[np.isnan(diffsit)]=0.0

    ax = plt.subplot(2,4,i+1)    
    m = Basemap(projection='npstere',boundinglat=59,lon_0=270,
                resolution='l',round =True,area_thresh=1000.)              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='k',linewidth=0.2)
    m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')
    
    x,y = m(lonvals[i],latvals[i])
    cs = m.hexbin(x,y,C=diffsit,vmin = -1,vmax = 1)
    
    cmap = ncm.cmap('NCV_blu_red')            
    cs.set_cmap(cmap)  
    
    ax.annotate(r'\textbf{%s}' % years[i], xy=(0, 0), 
                 xytext=(0.7,0.97),xycoords='axes fraction',
                fontsize=15,color='dimgrey',rotation=0)
    
cbar_ax = fig.add_axes([0.313,0.15,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=False)                      
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))    
cbar.set_label(r'\textbf{Difference (m)}',color='dimgrey')
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(hspace=-0.3)
plt.subplots_adjust(wspace=0.01)
    
plt.savefig(directoryfigure + 'diff_SIT.png',dpi=300)

###########################################################################
### Plot CryoSat-2 magnitude
limit = np.arange(0,6,1)
barlim = np.arange(0,6,1)        
        
fig = plt.figure()
for i in xrange(years.shape[0]):
    varsit = sitc[i]
    varsit[np.isnan(varsit)]=0.0

    ax = plt.subplot(2,4,i+1)    
    m = Basemap(projection='npstere',boundinglat=59,lon_0=270,
                resolution='l',round =True,area_thresh=1000.)              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='k',linewidth=0.2)
    m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')
    
    x,y = m(lonvals[i],latvals[i])
    cs = m.hexbin(x,y,C=varsit,vmin = 0,vmax = 5)
               
    cs.set_cmap('cubehelix')  
    
    ax.annotate(r'\textbf{%s}' % years[i], xy=(0, 0), 
                 xytext=(0.7,0.97),xycoords='axes fraction',
                fontsize=15,color='dimgrey',rotation=0)
    
cbar_ax = fig.add_axes([0.313,0.15,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)                      
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))    
cbar.set_label(r'\textbf{Thickness (m)}',color='dimgrey')
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(hspace=-0.3)
plt.subplots_adjust(wspace=0.01)
    
plt.savefig(directoryfigure + 'cs2_SIT.png',dpi=300)

###########################################################################
### Plot PIOMAS magnitude
limit = np.arange(0,6,1)
barlim = np.arange(0,6,1)        
        
fig = plt.figure()
for i in xrange(years.shape[0]):
    varsit = sitp[i]
    varsit[np.isnan(varsit)]=0.0

    ax = plt.subplot(2,4,i+1)    
    m = Basemap(projection='npstere',boundinglat=59,lon_0=270,
                resolution='l',round =True,area_thresh=1000.)              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='k',linewidth=0.2)
    m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')
    
    x,y = m(lonvals[i],latvals[i])
    cs = m.hexbin(x,y,C=varsit,vmin = 0,vmax = 5)
               
    cs.set_cmap('cubehelix')  
    
    ax.annotate(r'\textbf{%s}' % years[i], xy=(0, 0), 
                 xytext=(0.7,0.97),xycoords='axes fraction',
                fontsize=15,color='dimgrey',rotation=0)
    
cbar_ax = fig.add_axes([0.313,0.15,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)                      
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))    
cbar.set_label(r'\textbf{Thickness (m)}',color='dimgrey')
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(hspace=-0.3)
plt.subplots_adjust(wspace=0.01)
    
plt.savefig(directoryfigure + 'piomas_SIT.png',dpi=300)

###########################################################################
### Plot absolute differences
limit = np.arange(-1,1.1,0.1)
barlim = np.arange(-1,2,1)        
        
fig = plt.figure()
for i in xrange(years.shape[0]):    
    ax = plt.subplot(2,4,i+1,aspect='equal')               
    adjust_spines(ax, ['left', 'bottom'])            
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_color('dimgrey')
    ax.spines['bottom'].set_color('dimgrey')
    ax.tick_params(axis='both',direction='out',length=4,width=2,
                   which='major',pad=4,color='dimgray')
    
    varx = sitc[i]
    vary = sitp[i]

    plt.plot(timexi[i],linei[i],linewidth=1.5,linestyle='-',
             color='m',zorder=3)
    plt.plot(timexi[i],timeyi[i],color='k',linestyle='-',linewidth=2,
             zorder=2)
    plt.scatter(varx,vary,s=9,color='dodgerblue',edgecolor='darkblue',
                linewidth=0.2,alpha=0.7,zorder=1)
    plt.xlim([0,4])
    plt.ylim([0,4])
    plt.xticks(np.arange(0,5,1),map(str,np.arange(0,5,1)),fontsize=7)
    plt.yticks(np.arange(0,5,1),map(str,np.arange(0,5,1)),fontsize=7)
    
    ax.annotate(r'\textbf{%s, R$^{2}$=%s}' \
                % (years[i],round(r_valuei[i]**2,2)),
                xy=(0, 0),xytext=(0.05,1.02),xycoords='axes fraction',
                fontsize=9,color='dimgrey',rotation=0)
                
    ax.annotate(r'\textbf{CryoSat-2 [m]}',
            xy=(0, 0),xytext=(0.395,0.1),xycoords='figure fraction',
            fontsize=14,color='k',rotation=0)
    ax.annotate(r'\textbf{PIOMAS [m]}',
        xy=(0, 0),xytext=(0.05,0.64),xycoords='figure fraction',
        fontsize=14,color='k',rotation=90)
    
    plt.subplots_adjust(hspace=0.)
    plt.subplots_adjust(wspace=0.3)
    
plt.savefig(directoryfigure + 'scatter_SIT.png',dpi=300)