from MOOSEplot import geoplot as gp
from MOOSEplot import plotCommon as pc
from datetime import datetime
from matplotlib import pyplot as plt

## First define a StandardPlot object from GeoPlot.
Gplot=gp.StandardPlot()

## Load and Parse the namelist.
## Note, you will likely need to change the wrf_path in the namelist
## Such that it works with this namelist.

Gplot.parse_namelist(namelist_file='namelist_ex1.txt')

## Optional, uncomment the below line to print out all of the namelist values.
#Gplot.print_namelist()

## Define the map projection based on current namelist values.
Gplot.define_projection()

#Define Figure and subplot
figure=plt.figure(figsize=(12,9))
#Note that the Cartopy framework is already imported as part of GeoPlot, and the projection is already defined as a cartopy projection
ax=figure.add_subplot(1,1,1,projection=Gplot.projection)

# Plot the model data on the defined subplot.
Gplot.makePlot(ax)

## Add the metars to the data
## Note that numerous arguments are required to match the plotting conventions defined in the StandardPlot
## METAR .csv files are located in the local SW_METARS directory.  Observed Temperature (in C) is defined under the "TMP" title
## Min,Max, and color map are manually set to match the values in the namelist.
## Time is set to the namelist time-step and format, and the surface station lat/lon is saved under the latitude, and longitude headers in the csv file.
## Vectors is set equal to false, so wind barbs are NOT plotted.
pc.add_metars(ax,path='SW_METARS/',vari='TMP',vmin=10.,vmax=45.,cmap='jet',edgecolor='k',marker='o',time=Gplot.namelist_dictionary['time_step'],
              tfmt=Gplot.namelist_dictionary['time_format'],lon='longitude',lat='latitude',zorder=8,vectors=False)


## Finally, define a time-string using the datetime function and print out a title.
timestr=datetime.strptime(Gplot.namelist_dictionary['time_step'],Gplot.namelist_dictionary['time_format']).strftime('%Y%m%d %H UTC')
plt.title(timestr,loc='right')


Gplot.show() ## Show the figure onto the screen
