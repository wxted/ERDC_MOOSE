from MOOSEplot import MeteoPlot as mp
from matplotlib import pyplot as plt
import numpy as np

## Define list of 8 stations to use in histogram analysis
stns=['BXK','PHX','CGZ','GEU','FFZ','GXF','IGM','HII']
n_bins=25 ## Number of bins to use in histogram

## set the time-range (start, end)
time_range=['2014-07-03_14:00:00','2014-07-04_05:00:00']


## Define figure and subplot
fig=plt.figure(figsize=(12,10))
ax=plt.subplot(1,1,1)

## loop to load station data for EACH station.
for sdx, s in enumerate(stns):

    ## Define a NEW object for each station with the time range specified above
    mplot=mp.MeteoPlot(stn=s,input_path='./METARS/',timerange=time_range)

    ## access and set wrf_attrs dictionary to load wrf data
    mplot.wrf_attrs['path']='./WRF_POST/'
    mplot.wrf_attrs['pfx']='wrfoutput_post_chemPOST_d03'

    ### load WRF data into the "wrfdf" dataframe
    wrfdf=mplot.load_wrf_data(variable='TMP')
    mdata=mplot.data ## access the suface observations DATA frame directly

    ## THIS Second of code is used to linearly interpolate the station data to WRF time coordinate
    ## Helpful, since stations don't ALWAYS report on the hour, and some surface data reports more than once per hour
    x = np.asarray(wrfdf['TIME'],dtype=np.float64) ## Set WRF time coordinate to a numpy float
    xp = np.asarray(mdata['TIME'],dtype=np.float64) ## Set station time coordinate to numpy float
    fp = np.asarray(mdata['TMP']) ## Grab the station temperature data
    ts = np.interp(x,xp,fp) ## Numpy interp function to interpolate the station temperature data to the WRF data

    ## start building a list of "difference" values that make up the histogram
    if sdx == 0: ## if this is the first station, define the list
        diff=list(wrfdf['TMP'].values-ts) ## model - obs

    else: ## otherwise, add the NEW difference values to the list to expand the list
        diff=diff+list(wrfdf['TMP'].values-ts)

## plot the histogram of differences.
ax.hist(diff, bins=n_bins,width=0.6,edgecolor='k',density=True,color='#FF8A8A')

## ALL CODE BELOW THIS LINE USES MATPLOTLIB TO PLOT THE FIGURE ANNOTATIONS AND
## THE PERCENTILE LINES --> Not required but helps add additional information to the figure
## SEE: https://matplotlib.org/3.1.1/tutorials/text/annotations.html

## 25th percentile
ax.axvline(np.percentile(diff,25),ls='--',color='b',lw=0.9)
ax.annotate("25$^{th}$: %.1f"%np.percentile(diff,25),
            xy=(np.percentile(diff,25), 0.15), xycoords='data',
            xytext=(np.percentile(diff,25)-2.5, 0.12), textcoords='data',
            size=13, va="center", ha="center",
            arrowprops=dict(arrowstyle="simple",
                            connectionstyle="arc3,rad=-0.2"),
            )


## 75th percentile
ax.axvline(np.percentile(diff,75),ls='--',color='b',lw=0.9)
ax.annotate("75$^{th}$: %.1f"%np.percentile(diff,75),
            xy=(np.percentile(diff,75), 0.15), xycoords='data',
            xytext=(np.percentile(diff,75)+2.5, 0.16), textcoords='data',
            size=13, va="center", ha="center",
            arrowprops=dict(arrowstyle="simple",
                            connectionstyle="arc3,rad=-0.2"),
            )

## Median
ax.axvline(np.median(diff),ls='--',color='k',lw=1.6)
ax.annotate("median: %.1f"%np.median(diff),
            xy=(np.median(diff), 0.19), xycoords='data',
            xytext=(np.median(diff)-1.25, 0.16), textcoords='data',
            size=13, va="center", ha="center",
            arrowprops=dict(arrowstyle="simple",
                            connectionstyle="arc3,rad=-0.2"),
            )

### finally, plot the title, labels, and add a grid to the plot.
ax.set_title("Histogram of Model - Observed Temperature \nbetween %s - %s"%(time_range[0],time_range[1]),loc='left')
ax.set_ylabel("Fraction of data")
ax.set_xlabel("Temperature [C]")
ax.grid()
plt.show()
