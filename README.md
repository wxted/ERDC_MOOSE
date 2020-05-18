# ERDC-MOOSE

## The ERDC Model Output and Observations Scripting Engine (ERDC-MOOSE) 
 
ERDC-MOOSE is a set of python scripts designed specifically to more easily handle WRF model output and ERA-INTERIM reanalysis
gridded data sets to generate both quick look plots and publication quality figures.  Key capabilities listed below:

See Ted Letchers ERDC Technical Report: 'xxxxxx'

-   WrfPost: A post processing package that subsets raw wrfout files into a single netCDF file.
-   EraPost: A post processing package that subsets raw ERA-Interim datafiles in grib2 format from NCARs RDA.
-   GeoPlot: A geographic plotting package that allows users to set standard plotting options, that generates a geographic figure
-   MeteoPlot: A function to allow for loading and timeseries plotting of surface weather observations and model data


## Installation:

In the command line type:

    `python setup install`