"""
Scripts reads in sea ice thickness data from CS-2 corrected/uncorrected
and plots a trend analysis over the 2011-2017 period (April)
 
Notes
-----
    Author : Zachary Labe
    Date   : 16 January 2018
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as c
import datetime
from mpl_toolkits.basemap import Basemap
import nclcmaps as ncm
import cmocean
import scipy.stats as sts

### Define directories
directorydata = '/home/zlabe/Documents/Projects/CAAthickness/Data/'
directoryfigure = '/home/zlabe/Desktop/CS2PIOMAS/Corrected/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Calc sea ice thickness trends - %s----\n' % titletime) 

### Alott time series
yearmin = 2011
yearmax = 2017
years = np.arange(yearmin,yearmax+1,1)
yearsq = np.append(years,years)
months = [r'April']

### Read in data
lonvals = []
latvals = []
sitp= []
situ = []
sitc = []
for i in range(years.shape[0]):
    filename = 'PIOMAS_CS2UC_sit_04_%s.txt' % years[i]
    data = np.genfromtxt(directorydata + filename,unpack=False,
                         usecols=[0,1,2,3,4],skip_header=2)
    latq = data[:,0]
    lonq = data[:,1]
    sitpq = data[:,2]
    situq = data[:,3]
    sitcq = data[:,4]
    
    lonvals.append(lonq)
    latvals.append(latq)
    sitp.append(sitpq)
    situ.append(situq)
    sitc.append(sitcq)

    print('Completed: Read in %s data!' % years[i])   
    
### Appends lat/lon
lonvalsn = np.append(lonvals,lonvals)
latvalsn = np.append(latvals,latvals)

### Conduct a weighted average
def calc_weightedAve(var,lats):
    """
    Area weights sit array 5d [ens,year,month,lat,lon] into [ens,year,month]
    
    Parameters
    ----------
    var : 5d,4d,3d,1d array of a gridded variable
    lats : 2d array of latitudes
    
    Returns
    -------
    meanvar : weighted average for 3d,2d,1d array

    Usage
    -----
    meanvar = calc_weightedAve(var,lats)
    """
    print('\n>>> Using calc_weightedAve function!')
    
    ### Import modules
    import numpy as np
    
    ### Calculate weighted average for various dimensional arrays
    if var.ndim == 5:
        meanvar = np.empty((var.shape[0],var.shape[1],var.shape[2]))
        for ens in range(var.shape[0]):
            for i in range(var.shape[1]):
                for j in range(var.shape[2]):
                    varq = var[ens,i,j,:,:]
                    mask = np.isfinite(varq) & np.isfinite(lats)
                    varmask = varq[mask]
                    areamask = np.cos(np.deg2rad(lats[mask]))
                    meanvar[ens,i,j] = np.nansum(varmask*areamask) \
                                        /np.sum(areamask)  
    elif var.ndim == 4:
        meanvar = np.empty((var.shape[0],var.shape[1]))
        for i in range(var.shape[0]):
            for j in range(var.shape[1]):
                varq = var[i,j,:,:]
                mask = np.isfinite(varq) & np.isfinite(lats)
                varmask = varq[mask]
                areamask = np.cos(np.deg2rad(lats[mask]))
                meanvar[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
    elif var.ndim == 3:
        meanvar = np.empty((var.shape[0],var.shape[1]))
        for i in range(var.shape[0]):
            varq = var[i,:,:]
            mask = np.isfinite(varq) & np.isfinite(lats)
            varmask = varq[mask]
            areamask = np.cos(np.deg2rad(lats[mask]))
            meanvar[i] = np.nansum(varmask*areamask)/np.sum(areamask)
    elif var.ndim == 1:
        meanvar = np.empty((var.shape[0]))
        varq = var[:]
        mask = np.isfinite(varq) & np.isfinite(lats)
        varmask = varq[mask]
        areamask = np.cos(np.deg2rad(lats[mask]))
        meanvar = np.nansum(varmask*areamask)/np.sum(areamask)
    else:
        ValueError('Variable has the wrong dimensions!')
     
    print('Completed: Weighted variable average!')
    
    print('*Completed: Finished calc_weightedAve function!')
    return meanvar

### Calculate average over the entire Arctic for each data set
meanp = np.empty((len(sitc)))
meanu = np.empty((len(sitc)))
meanc = np.empty((len(sitc)))
for i in range(meanp.shape[0]):
    meanp[i] = calc_weightedAve(sitp[i],latvals[i])
    meanu[i] = calc_weightedAve(situ[i],latvals[i])
    meanc[i] = calc_weightedAve(sitc[i],latvals[i])
meansit = [meanp,meanu,meanc]

### Calculate linear trends
def calc_trendLine(varx): 
    """
    Calculates linear regression
    
    Parameters
    ----------
    var : 1d array
    
    Returns
    -------
    line : values for plotting linear line
    slope : slope of linear regression

    Usage
    -----
    line,slope = calc_weightedAve(var,lats)
    """
    print('\n>>> Using calc_trendLine function!')
    
    ### Create array of values for regression
    timey = np.arange(len(varx))

    ### Scipy.stats linear regression
    slope,intercept,r_value,p_value,std_error = sts.linregress(timey,varx)
    line = slope*timey + intercept
    
    print('Completed: Linear regression!')
    
    print('*Completed: Finished calc_trendLine function!')
    return line,slope

### Calculate trend line for each data set
linesitp,slopesitp = calc_trendLine(meanp)
linesitu,slopesitu = calc_trendLine(meanu)
linesitc,slopesitc = calc_trendLine(meanc)
linesit = [linesitp,linesitu,linesitc]
    
###############################################################################
###############################################################################
###############################################################################
### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 5))
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
        
fig = plt.figure()
ax = plt.subplot(111)

varnames = [r'PIOMAS',r'CS2$_u$',r'CS2$_c$']

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='darkgrey')

color=iter(cmocean.cm.thermal(np.linspace(0.1,0.8,len(meansit))))
for i in range(len(meansit)):
    c=next(color)
    plt.plot(meansit[i],linewidth=2.5,color=c,alpha=1,linestyle='-',
             marker='o',markersize=7,label=r'\textbf{%s}' % varnames[i])
    plt.plot(linesit[i],linewidth=0.9,color=c,alpha=1,linestyle='--')

plt.ylabel(r'\textbf{Sea Ice Thickness (m)}',color='k',fontsize=13)
plt.yticks(np.arange(0,2.1,0.5),list(map(str,np.arange(0,2.1,0.5))))
plt.ylim([0,2])

ax.yaxis.grid(zorder=1,color='dimgrey',alpha=0.6)
    
xlabels = list(map(str,years))
plt.xticks(np.arange(0,7,1),xlabels)
plt.xlim([0,6])

plt.legend(shadow=False,fontsize=13,loc='lower center',
           fancybox=True,frameon=False,ncol=4)
    
plt.savefig(directoryfigure + 'meanSITtrend_04_1117.png',dpi=300)