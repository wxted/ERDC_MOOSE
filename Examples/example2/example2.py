from MOOSEplot import MeteoPlot as mp
from matplotlib import pyplot as plt

## Initialize the meteoplot object with the PHX identifer for "Phoenix Airport"
mplot=mp.MeteoPlot(stn='PHX',input_path='./',timerange=['2014-07-02_12:00:00','2014-07-04_12:00:00'])

## access the weather_keys list and update it to included Dust (DU) and Dust Storm (DS)
mplot.weather_keys=['RA','TS','DU','DS']
mplot.print_wx_attrs() ## Print out attributes of the weather keys (e.g., colors)

## access the wind_attrs dictionary and set specific keywords controlling the number of barbs to plot, and their length
mplot.wind_attrs['skip']=1
mplot.wind_attrs['barb_kwargs']['length']=4.3

## Define the variable dictionary with the Pandas header name from the .csv file, and the ylabel to use on the figure.

## initialize figure and subplot
fig=plt.figure(figsize=(8,12))
ax=plt.subplot(3,1,1) ##(num rows, num columns, current figure)

## add the temperature time-series
mplot.add_subplot(ax,'TMP',ylabel='Temperature',color='#d73027',ls='--',lw=2)

## Access the wrf_attrs dictionary and define the path, and file prefix
mplot.wrf_attrs['path']='/home/Desktop/WRF_POST/'
mplot.wrf_attrs['pfx']='wrfoutput_post_chemPOST_d03'

## Load the WRF data from the "var_dict" defined above into the local "wrfdf" dataframe
wrfdf=mplot.load_wrf_data(variable='TMP')

## plot the WRF data on the subplot.
mplot.plot_wrf_data(ax,wrfdf,variable='TMP',marker='o',ls='none',color='k')
## Set a few plot attributes
ax.set_xticklabels([])
ax.grid()
ax.set_title("Temperature [C]",loc='left')

## define second subplot
ax=plt.subplot(3,1,2)
## Add the observed weather from the weather keys
mplot.add_subplot(ax,'weather')
ax.grid()
ax.set_xticklabels([])
ax.set_title("Observed Weather",loc='left')

## define 3rd and final subplot
ax=plt.subplot(3,1,3)
## Set dataname to "wind" and plot the observed windspeed
mplot.add_subplot(ax,'wind',ylabel='Wind Speed')

## Load WRF wind speed into "wrfdf" dataframe
wrfdf=mplot.load_wrf_data(variable='wind')
## Plot WRF windspeed, NOTE does not plot barbs.
mplot.plot_wrf_data(ax,wrfdf,variable='wind',marker='o',ls='none',color='k')

ax.grid()
ax.set_title("Windspeed and observed Wind Barbs",loc='left')
plt.suptitle(mplot.meta_name) ## Plot overarching title with full station name from the .csv file.
plt.show()
