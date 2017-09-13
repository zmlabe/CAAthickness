"""
Scripts reads in snow depth data from CS-2 and interpolated PIOMAS 
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
directoryfigure = '/home/zlabe/Desktop/CS2PIOMAS/Snow/'

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
    snowpq = data[:,3] * 100.
    sitcq = data[:,4]
    snowcq = data[:,5] * 100.
    
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
timexs = []
timeys = []
lines = []
slopes = []
r_values = []
for i in xrange(years.shape[0]):    
    timexiq,timeyiq,lineiq,slopeiq,r_valueiq = calcStats(snowc[i],
                                                         snowp[i],years)    
    timexs.append(timexiq)
    timeys.append(timeyiq)
    lines.append(lineiq)
    slopes.append(slopeiq)
    r_values.append(r_valueiq)
    
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
barlim = np.arange(-40,41,20)        
        
fig = plt.figure()
for i in xrange(years.shape[0]):
    diffsnc = snowc[i] - snowp[i]
    diffsnc[np.isnan(diffsnc)]=0.0

    ax = plt.subplot(2,4,i+1)    
    m = Basemap(projection='npstere',boundinglat=59,lon_0=270,
                resolution='l',round =True,area_thresh=1000.)              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='k',linewidth=0.2)
    m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')
    
    x,y = m(lonvals[i],latvals[i])
    cs = m.hexbin(x,y,C=diffsnc,vmin = -40,vmax = 40)
    
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
cbar.set_label(r'\textbf{Difference (cm)}',color='dimgrey')
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(hspace=-0.3)
plt.subplots_adjust(wspace=0.01)
    
plt.savefig(directoryfigure + 'diff_snc.png',dpi=300)

###########################################################################
### Plot CryoSat-2 magnitude
limit = np.arange(0,41,20)
barlim = np.arange(0,41,20)        
        
fig = plt.figure()
for i in xrange(years.shape[0]):
    varsnc = snowc[i]
    varsnc[np.isnan(varsnc)]=0.0

    ax = plt.subplot(2,4,i+1)    
    m = Basemap(projection='npstere',boundinglat=59,lon_0=270,
                resolution='l',round =True,area_thresh=1000.)              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='k',linewidth=0.2)
    m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')
    
    x,y = m(lonvals[i],latvals[i])
    cs = m.hexbin(x,y,C=varsnc,vmin = 0,vmax = 40)
               
    cs.set_cmap('cool')  
    
    ax.annotate(r'\textbf{%s}' % years[i], xy=(0, 0), 
                 xytext=(0.7,0.97),xycoords='axes fraction',
                fontsize=15,color='dimgrey',rotation=0)
    
cbar_ax = fig.add_axes([0.313,0.15,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)                      
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))    
cbar.set_label(r'\textbf{Snow Depth (cm)}',color='dimgrey')
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(hspace=-0.3)
plt.subplots_adjust(wspace=0.01)
    
plt.savefig(directoryfigure + 'cs2_snc.png',dpi=300)

###########################################################################
### Plot PIOMAS magnitude
limit = np.arange(0,41,20)
barlim = np.arange(0,41,20)       
        
fig = plt.figure()
for i in xrange(years.shape[0]):
    varsnc = snowp[i]
    varsnc[np.isnan(varsnc)]=0.0

    ax = plt.subplot(2,4,i+1)    
    m = Basemap(projection='npstere',boundinglat=59,lon_0=270,
                resolution='l',round =True,area_thresh=1000.)              
    m.drawmapboundary(fill_color='white',color='dimgray',linewidth=0.7)
    m.drawcoastlines(color='k',linewidth=0.2)
    m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')
    
    x,y = m(lonvals[i],latvals[i])
    cs = m.hexbin(x,y,C=varsnc,vmin = 0,vmax = 40)
               
    cs.set_cmap('cool')  
    
    ax.annotate(r'\textbf{%s}' % years[i], xy=(0, 0), 
                 xytext=(0.7,0.97),xycoords='axes fraction',
                fontsize=15,color='dimgrey',rotation=0)
    
cbar_ax = fig.add_axes([0.313,0.15,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)                      
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))    
cbar.set_label(r'\textbf{Snow Depth (cm)}',color='dimgrey')
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(hspace=-0.3)
plt.subplots_adjust(wspace=0.01)
    
plt.savefig(directoryfigure + 'piomas_snc.png',dpi=300)

###########################################################################
### Plot scatter
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
    
    varx = snowc[i]
    vary = snowp[i]

    plt.plot(timexs[i],lines[i],linewidth=1.5,linestyle='-',
             color='firebrick',zorder=3)
    plt.plot(timexs[i],timeys[i],color='k',linestyle='-',linewidth=2,
             zorder=2)
    plt.scatter(varx,vary,s=9,color='mediumslateblue',edgecolor='indigo',
                linewidth=0.2,alpha=0.7,zorder=1)
    plt.xlim([0,40])
    plt.ylim([0,40])
    plt.xticks(np.arange(0,41,10),map(str,np.arange(0,41,10)),fontsize=7)
    plt.yticks(np.arange(0,41,10),map(str,np.arange(0,41,10)),fontsize=7)
    
    ax.annotate(r'\textbf{%s, R$^{2}$=%s}' \
                % (years[i],round(r_values[i]**2,2)),
                xy=(0, 0),xytext=(0.05,1.02),xycoords='axes fraction',
                fontsize=9,color='dimgrey',rotation=0)
                
    ax.annotate(r'\textbf{Warren Snow [cm]}',
            xy=(0, 0),xytext=(0.3933,0.1),xycoords='figure fraction',
            fontsize=14,color='k',rotation=0)
    ax.annotate(r'\textbf{PIOMAS [cm]}',
        xy=(0, 0),xytext=(0.05,0.64),xycoords='figure fraction',
        fontsize=14,color='k',rotation=90)
    
    plt.subplots_adjust(hspace=0.)
    plt.subplots_adjust(wspace=0.3)
    
plt.savefig(directoryfigure + 'scatter_snc.png',dpi=300)