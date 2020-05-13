# coding=utf-8
"""This file contains the object class and function definition for MeteoPlot.
   files into a single netcdf file with user specified variables.

   Note: This script is very specifically designed to handle surface data related to dust variables.
   It is designed to take data and output it into a format common to the WRF-Chem output such that it can be
   quickly and easily compared using the plotting functions within this larger package.

   Defining class requires a station input that has a .csv station file associated with it created by
   the GetStn class within the this post-process package.  Also takes optional time boundaries.

   Once class is set, a variety of functions allows for fast plotting of meteorlogical time-series variables, and
   addition of WRF, or other gridded met-data.  Works best if the gridded model data has been post processed through
   the post_processing classes in this package.

   Example:

       mplot = MeteoPlot(stn='LEB',input_path='path_to/station_data/',timerange='Auto')
       mplot.add_subplot(

       Inputs:
        - stn = Prefix for met station (without preceeding "K" if applicable)
        - figsize = figure size --> tuple (length , height)
        - input_path = path where .csv METAR files are stored locally.
        - subplots = intial number of subplots in the figure, assumes they are vertically stacked.
        - timerange (default='Auto') --> optional list string times with a start (min) time and end (max) time.
        - tfmt = time format associated with string times in 'timerange'
        - pvar = determine which pressure variable to use in the METAR data (ALTM, or PMSL) ALTM usually works best.


    Coding format follows PEP8 (ish) rules, but not perfectly.
"""
import numpy as np
from matplotlib import pyplot as plt

from datetime import datetime, timedelta
from MOOSEplot import plotCommon as pc

import pandas as pd
import glob as glob
import sys

class MeteoPlot:
    def __init__(self, stn='LEB', figsize=(7, 11), input_path='./', subplots=4, timerange='Auto', pvar='ALTM',
                 tfmt='%Y-%m-%d_%H:%M:%S'):
        """ Initialize Meteogram plot. -- Very specific file, created from meteogram post processing script."""

        try:
            self.file = glob.glob(input_path + '/*%s*' % stn)[0]
        except:
            print("There are no files for station %s in folder %s" % (stn, input_path))
            print("Please check your paths and rerun the script.")
            sys.exit()

        self.TMPunitsAllowed = ('c', 'f', 'k', 'celsius', 'kelvin', 'fahrenheit')
        self.PRESunitsAllowed = ('inhg', 'mb', 'hpa', 'pa', 'pascals')
        self.VISunitsAllowed = ('feet', 'miles', 'mi', 'm', 'meters', 'kilometers', 'km', 'sm', 'ft')
        self.SPDunitsAllowed = ('m/s', 'mph', 'knots', 'kt', 'ms-1', 'mps')
        self.VOLunitsAllowed = ('kg m^-3','kg/m3','kg/m^3','g m^-3','g/m3','g/m^3','g cm^-3','g/cm3','g/cm^3',
                                'ug m^-3','ug/m3','ug/m^3','ug cm^-3','ug/cm3','g/cm^3','micrograms','ppmv','ppm')

        self.pvar = pvar
        self.stnID=stn
        self.weather_keys = weather_keys = ['RA', 'TS', 'SN', 'UP']
        self.weather_col = {'RA': ['#66c2a4', '#2ca25f', '#006d2c'],
                            'SN': ['#ece7f2', '#a6bddb', '#2b8cbe'],
                            'FG': 2 * ['#bdbdbd'],
                            'HZ': 2 * ['#fdae6b'], 'BR': 2 * ['#bdbdbd'],
                            'BD': 2 * ['#d8b365'], 'DU': 2 * ['#d8b365'], 'DS': 2 * ['#8c510a'],
                            'BLSN': ['#fde0dd', '#fa9fb5', '#c51b8a'],
                            'UP': ['#ffeda0', '#feb24c', '#f03b20'], 'TS': 2 * ['#2166ac']}

        # GET THE META-DATA AT THE TOP OF THE FILE
        with open(self.file, 'r') as f:
            self.meta_name = f.readline()
            UnitName = f.readline()

        print("Reading Data For: %s" % self.meta_name)

        self.__UnitDict = {}  # Create unit dictionary
        for i in UnitName.split('|'):
            self.__UnitDict[i.split(':')[0]] = i.split(':')[1].strip()

        df = pd.read_csv(self.file, header=2, parse_dates=['TIME'])
        df = df[np.isnan(df[pvar]) == False].reset_index()  # remove any values with bad altimiter

        self.UnitDict = UD(self.__UnitDict)
        ## quick sanity check to see if timerange is a list or not.
        if not isinstance(timerange, list):
            timerange = 'auto'
        else:
            # ASSUME IT'S A LIST OF TIME VALUES!
            # Reset the dateframe to focus solely on the time range specifeid in the function.
            tmin = pd.to_datetime(datetime.strptime(timerange[0], tfmt))
            tmax = pd.to_datetime(datetime.strptime(timerange[1], tfmt))
            mask = (df['TIME'] > tmin) & (df['TIME'] <= tmax)
            df = df.loc[mask]

        self.data = df.reset_index()

        self.latname='latitude'
        self.lonname='longitude'

        self.wrflatname='Latitude'
        self.wrflonname='Longitude'

        self.wrf_attrs={'path':'./','pfx':'','post':True,'tmp_name':'2 metre temperature',
                        'u_name':'10 metre U wind component',
                        'v_name':'10 metre V wind component',
                        'metVis_name':'Visibility (Surface)',
                        'extcof55_name':'EXTCOF55 (Surface)',
                        'aod_name':'Aerosol Optical Depth',
                        'pm10_name':'PM10 (Surface)'}

        self.wind_attrs={'barbs':True,'skip':3,'spd_name':'SPD','dir_name':'DIR','barb_kwargs':{},'spd_kwargs':{}}

        self.lat = self.data[self.latname][0]
        self.lon = self.data[self.lonname][0]


    def __repr__(self):
        return (f'{self.__class__.__name__}('
                f'Station Id {self.stnID!r}, Description={self.meta_name!r})')

    def __len__(self):
        return len(self.data[self.pvar])

    def print_wind_attrs(self):
        """Help on built-in function print_wx_attrs in module print_wind_attrs:

                print_wind_attrs(...)
                    print_wind_attrs()

            User-callable function that prints out a list of wind attributes that are read through the add_subplot function when
            dataname = wind

            attributes:
                - barbs (bool) = Plot wind barbs on wind speed time-series
                - skip (integer) = controls how many barbs are plotted -> higher number = fewer barbs
                - spd_name (string) = name of windspeed variable in observational dataframe (defaults to SPD)
                - dir_name (string) = name of meteorlogical wind direction in observational data frame (defaults to DIR)
                - barb_kwargs (dictionary) = dictionary with keyword arguments specific to Matplotlibs "barb" function
                    - see: https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.barbs.html
                - spd_kwargs (dictionary) = dictionary for keyword arguments for the "wind speed" time-series

            Returns
            -------
            nothing.
        """
        for key in self.wind_attrs:
            print('%s --> %s'%(key,self.wind_attrs[key]))

    def print_wx_attrs(self,addons=['-','','+']):
        """Help on built-in function print_wx_attrs in module MeteoPlot:

                 print_wx_attrs(...)
                     print_wx_attrs()

            User-callable function that prints out a list of weather report
            attributes

            Returns
            -------
            nothing.
        """
        for key in self.weather_keys:
            print("Weather :%s"%key)
            if len(set(self.weather_col[key])) > 1:
                print("  -Multiple levels-")
                for idx,i in enumerate(set(self.weather_col[key])):
                    print('-- %s -> color = %s'%(addons[idx],i))
            else:
                '  - color = %s'%self.weather_col[key][0]


    def QuickPlot(self):
        """Help on built-in function QuickPlot in module MeteoPlot:

                 QuickPlot(...)
                     QuickPlot()

            User-callable function that generates a quick look figure with time-series of
                - Temperature
                - Weather Reports
                - Wind speed and Wind Vectors
                - Visibility

            Returns
            -------
            matplotlib figure
        """
        figure=plt.figure(figsize=(10,7))
        ax=plt.subplot(4,1,1)
        self.add_subplot(ax, marker='o', ls='none', color='k')
        plt.grid()
        ax.set_xticklabels([])
        plt.title("Meteogram for: %s" % self.meta_name)

        ax = plt.subplot(4, 1, 2)
        self.add_subplot(ax,'weather')
        plt.grid()
        ax.set_xticklabels([])

        ax = plt.subplot(4, 1, 3)
        self.add_subplot(ax,'wind', ylabel ='Wind Speed')
        plt.grid()
        ax.set_xticklabels([])

        ax = plt.subplot(4, 1, 4)
        self.add_subplot(ax,'VIS', ylabel= 'Visibility')
        plt.grid()

        return figure


    def add_subplot(self,ax,variable='TMP',**kwargs):
        """Help on built-in function add_subplot in module MeteoPlot:

                 add_subplot(...)
                     add_subplot(ax[,variable],**kwargs)

            User-callable function to plot a time-series of meteorological data onto a subplot with
            optional keyword arguments to determine plot styling.

            Additional special keyword arguments, specific to "add_subplot":
                - kind: (string: line,bar,scatter), controls type of plot
                - ylabel: (string), adds ylabel to subplot

            Parameters
            -------
            ax: matplotlib subplot object, subplot to plot time-series to
            variable (optional): string, name of variable header in surface data pandas dataframe
            ** kwargs (optional): Any keyword argument allowed by matplotlib subplot + special keywords listed above

            Returns
            -------
            nothing.
        """
        ## SET SPECIAL KWARGS OPTIONS HERE ##
        kind = kwargs['kind'] if 'kind' in kwargs else 'line' ## get "kind" keyword argument
        ylabel = kwargs['ylabel'] if 'ylabel' in kwargs else '' ## get "ylabel" keyword argument
        special_keywords=('kind','ylabel')

        for i in special_keywords:
            if i in kwargs: del kwargs[i]

        if kind.lower() not in ('line','bar','scatter'):
            print("!!WARNING!! Plot type is not allowed: allowable types: line,bar,scatter")
            kind='line'

        if variable.lower() in ['weather','wx']:
            self.plot_weather_types(ax, self.data, self.weather_keys, self.weather_col)
        elif variable.lower() in ['wind']:
            self.plot_wind(ax,self.data,self.UnitDict,ylabel=ylabel,kind=kind)
        else:
            self.plot_regular_variable(ax, self.data, self.UnitDict, vari=variable,
                                       ylabel=ylabel,kind=kind, **kwargs)


    def plot_wrf_data(self,ax,wrfdf,variable='TMP',**kwargs):
        """Help on built-in function plot_wrf_data in module MeteoPlot:

                 plot_wrf_data(...)
                     add_subplot(ax,wrfdf,[,variable],**kwargs)

            User-callable function to plot a time-series of WRF-simulated data onto a subplot with
            optional keyword arguments to determine plot styling.  WRF-Simulated data must be stored within a
            pandas dataframe created using the MeteoPlot function "load_wrf_data"

            Additional special keyword arguments, specific to "add_subplot":
                - kind: (string: line,bar,scatter), controls type of plot
                - ylabel: (string), adds ylabel to subplot

            Parameters
            -------
            ax: matplotlib subplot object, subplot to plot time-series to
            wrfdf: pandas dataframe, dataframe storing WRF data, created using "load_wrf_data" function
            variable (optional): string, name of variable header in wrfdf pandas dataframe
            ** kwargs (optional): Any keyword argument allowed by matplotlib subplot + special keywords listed above

            Returns
            -------
            nothing.
        """
        if variable.lower() in ('wind','wnd','spd','spd10m','speed'):
            variable='SPD'
        elif variable.lower() in ('tmp','temp','temperature','t2m','tmp2m'):
            variable ='TMP'
        elif variable.lower() in ('vis','visibilty'):
            variable='VIS'
        elif variable.lower() in ('dust','dustvis'):
            variable='DUSTVIS'
        elif variable.lower() in ('aod','optical depth'):
            variable='AOD'

        ## SET SPECIAL KWARGS OPTIONS HERE ##
        kind = kwargs['kind'] if 'kind' in kwargs else 'line' ## get "kind" keyword argument
        special_keywords=('kind','ylabel')

        for i in special_keywords:
            if i in kwargs: del kwargs[i]

        if kind.lower() not in ('line','bar','scatter'):
            print("!!WARNING!! Plot type is not allowed: allowable types: line,bar,scatter")
            kind='line'

        if kind == 'line':
            ax.plot(wrfdf['TIME'], wrfdf[variable], **kwargs)
        elif kind == 'bar':
            ax.bar(wrfdf['TIME'], wrfdf[variable], **kwargs)
        elif kind == 'scatter':
            ax.scatter(wrfdf['TIME'], wrfdf[variable], **kwargs)

        return

    def load_wrf_data(self,variable='TMP'):
        """Help on built-in function plot_wrf_data in module MeteoPlot:

                 load_wrf_data(...)
                     load_wrf_data([variable])

            User-callable function to load data from the WRF netcdf file into a pandas data frame
            that can be either passed to the "plot_wrf_data" function, or accessed directly by user outside of
            the MeteoPlot Framework.

            Parameters
            -------

            variable (optional): string, name of variable header in wrfdf pandas dataframe:
                -Note that the variable header name must match one of the following allowed headers:

                - 'TMP': 'tmp','temp','temperature','t2m','tmp2m'
                - 'VIS': 'vis','visibilty'
                - 'DUSTVIS': 'dustvis','dust'
                - 'SPD': 'wind','wnd','spd','spd10m','speed'
                - 'AOD': 'aod', 'aerosol optical depth'
                - 'PM10': 'pm10'

            Returns
            -------
            pandas dataframe with the WRF data header, and an additional 'TIME' variable with the wrf time stored as
            a pandas datetime array.

        """
        from netCDF4 import Dataset
        wrffiles=sorted(glob.glob(self.wrf_attrs['path']+'*%s*'%self.wrf_attrs['pfx']))

        if len(wrffiles) == 0:
            # Check and see if there are any files, if not kick out of this function.
            print("There are no wrf files with the prefix %s in folder %s"%(self.wrf_attrs['path'],
                                                                            self.wrf_attrs['pfx']))
            return
        init=True
        for wdx,w in enumerate(wrffiles):
            wrfdata=Dataset(w,'r')

            if self.wrf_attrs['post'] == True: # This is a post processed file, assume time read in as such:
                times = [pd.to_datetime(datetime.strptime(wrfdata.variables['Time'].units,
                                                          wrfdata.variables['Time'].units) +
                                        timedelta(hours=float(ss))) for ss in wrfdata.variables['Time'][:]]
            else: # Assume it's a normal WRF file, and times gets read in differently.
                times = [pd.to_datetime(datetime.strptime(b"".join(tm).decode("utf-8")), '%Y-%m-%d_%H:%M:%S') for tm in
                            wrfdata.variables['Times']]

            if np.min(times) > np.max(self.data['TIME']):
                # Special case, min wrf time, is greater than max dataframe time (i.e., outside of bounds)
                # Close file and move along.
                wrfdata.close()
                continue

            if np.max(times) < np.min(self.data['TIME']):
                # Special case, max wrf time, is less than min dataframe time (i.e., outside of bounds)
                # Close file and move along.
                wrfdata.close()
                continue

            ## otherwise, just load the whole thing and sort it out later....?

            if variable.lower() in ('tmp','temp','temperature','t2m','tmp2m'):
                variable='TMP'
                if self.wrf_attrs['post'] == False:
                    if self.wrf_attrs['tmp_name'] !='T2':
                        ## Quick warning for user.
                        print("Warning, WRF tmp_name variable is not T2, but 'post' is False"
                              "Assuming this is an unprocessed WRF file and getting T2 from variables.")

                    vari=wrfdata.variables['T2'][:]

                else:
                    try:
                        vari=wrfdata.variables[self.wrf_attrs['tmp_name']][:]
                    except:
                        print("%s does not appear to be in this %s, skipping!"%(self.wrf_attrs['tmp_name'],w))
                        return

                ## Convert WRF Temperature Units (Kelvin) to Units stored in UnitDict
                vari=self.__ConvertTmpUnits(vari,'K',self.UnitDict[variable])

            if variable.lower() in ('vis','visibilty'):
                variable='VIS'
                if self.wrf_attrs['post'] == False:
                    print("Visibility is not an output associated with unprocessed WRF files"
                          "and as of right now, I am not able to compute visibility from unprocessed WRF files"
                          "in this function.  Please post-process the WRF files using wrfpost and try again."
                          "Skipping Visibility....")
                    continue
                else:
                    try:
                        vari=wrfdata.variables[self.wrf_attrs['metVis_name']][:]
                        variunits=wrfdata.variables[self.wrf_attrs['metVis_name']].units
                    except:
                        print("%s does not appear to be in this %s, skipping!"%(self.wrf_attrs['metVis_name'],w))
                        return
                    ## Convert WRF Visibility Units to Units stored in UnitDict
                    ## WRF visibility units come from the netCDF units attribute.
                    vari = self.__ConvertLenUnits(vari, variunits, self.UnitDict[variable])

            if variable.lower() in ['wx','weather']:
                print("This is not performed in this function, please call the function wrf_wx "
                      "(which doesn't exist yet in MOOSE!")

            if variable.lower() in ('dustvis','dust'):
                variable='DUSTVIS'
                if self.wrf_attrs['post'] == False:
                    if self.wrf_attrs['extcof55_name'] != 'EXTCOF55':
                        print("Warning, WRF 'extcof55_name' variable is not EXTCOF55, but 'post' is False"
                              "Assuming this is an unprocessed WRF file and getting EXTCOF55 from variables.")
                        SFCEXTCOF=wrfdata.variables['EXTCOF55'][:,0,:] #surface level only.
                else:
                    try:
                        SFCEXTCOF=wrfdata.variables[self.wrf_attrs['extcof55_name']][:]
                    except:
                        print("%s does not appear to be in this %s, skipping!"%(self.wrf_attrs['extcof55_name'],w))
                        return

                vari = 3.912 / SFCEXTCOF  ## in kilometers.
                ## Convert WRF dust visibility units to units stored in UnitDict
                vari = self.__ConvertLenUnits(vari, 'km', self.UnitDict['VIS'])


            if variable.lower() in ('aod','aerosol optical depth'):
                variable='AOD'
                if self.wrf_attrs['post'] == False:
                    print("AOD is not an output associated with unprocessed WRF files"
                          "and as of right now, I am not able to compute AOD from unprocessed WRF files"
                          "in this function.  Please post-process the WRF files using wrfpost and try again."
                          "Skipping AOD....")
                    continue
                else:
                    try:
                        vari=wrfdata.variables[self.wrf_attrs['aod_name']][:]
                    except:
                        print("%s does not appear to be in this %s, skipping!"%(self.wrf_attrs['aod_name'],w))
                        return

            if variable.lower() in ('pm10'):
                variable='PM10'
                if self.wrf_attrs['post'] == False:
                    try:
                        vari = wrfdata.variables['PM10'][:,0,:]
                    except:
                        print("PM10 is not in this %s file, are you sure it's a WRF-Chem file?"
                              "Skipping PM10 ...")
                else:
                    try:
                        vari=wrfdata.variables[self.wrf_attrs['pm10_name']][:]
                    except:
                        print("%s does not appear to be in this %s, skipping!"%(self.wrf_attrs['aod_name'],w))
                        return

                vari = self.__ConvertVolConcUnits(vari, 'ug/m^3', self.UnitDict['PM10'])

            ## otherwise, just load the whole thing and sort it out later....?
            if variable.lower() in ('wind','wnd','spd','spd10m','speed'):
                variable='SPD'
                if self.wrf_attrs['post'] == False:
                    if self.wrf_attrs['u_name'] !='U10':
                        ## Quick warning for user.
                        print("Warning, WRF u_name variable is not U10, but 'post' is False"
                              "Assuming this is an unprocessed WRF file and getting U10/V10 from variables.")
                    u=wrfdata.variables['U10'][:]
                    v=wrfdata.variables['V10'][:]
                else:
                    try:
                        u=wrfdata.variables[self.wrf_attrs['u_name']][:]
                        v=wrfdata.variables[self.wrf_attrs['v_name']][:]
                    except:
                        print("%s does not appear to be in this %s, skipping!"%(self.wrf_attrs['u_name'],w))
                        return

                vari=np.sqrt(u**2.+v**2.)
                ## Convert data from wrf-speed (Assuming m/s) to Units stored in Variable Dictionary.
                vari = self.__ConvertSpdUnits(vari, 'mps', self.UnitDict[variable])


            if init == True:
                if self.wrf_attrs['post'] == False:
                    print("Defaulting to XLONG,XLAT as lat/lon names for unprocessed file...")
                    wrflon=wrfdata.variables['XLONG'][0,:].squeeze()
                    wrflat=wrfdata.variables['XLAT'][0,:].squeeze()
                else:
                    wrflon = wrfdata.variables[self.wrflonname][:].squeeze()
                    wrflat = wrfdata.variables[self.wrflatname][:].squeeze()

                x,y=pc.find_nearest_latlon(wrflon,wrflat,self.lon,self.lat)

                all_vari=list(vari[:,x,y].squeeze())
                all_times=times

                init=False
            else:
                all_vari=all_vari+list(vari[:,x,y].squeeze())
                all_times=all_times+times

            wrfdata.close()

        ## Make dictionary, then make dataframe! ##

        try:
            wrf_data_dict={'TIME':all_times,variable:all_vari}
            wrfdf=pd.DataFrame.from_dict(wrf_data_dict).sort_values(by='TIME')

            mask = (wrfdf['TIME'] > np.min(self.data['TIME'])) & (wrfdf['TIME'] <= np.max(self.data['TIME']))
            wrfdf = wrfdf.loc[mask]

        except:
            print("Something went wrong, possibly with the time matching ..."
                      "in anycase, I can't do the WRF plotting, so I'm not going to.")
            return
        return wrfdf

    def __ConvertLenUnits(self, data, inunit, newunit):
        """Help on built-in function plot_wrf_data in module _ConvertLenUnits:

                 _ConvertLenUnits(...)
                     _ConvertLenUnits(data, inunit, newunit)

            "private" function to convert units of length for a dataframe.  Called from the user-callable
            Convert Units function

            Parameters
            -------

            data: numpy array, array containing the data to be converted:
            inunit: string, unit corresponding to the units in the input data
            newunit: string, unit corresponding to the units the data will be converted to


            Returns
            -------
            numpy data array converted from inunit to outunit.
        """

        if inunit == 'm':  ## CONVERT meters to KM
            Ubase = data
        elif inunit == 'km':
            Ubase = data * 1000.
        elif inunit == 'sm':
            Ubase = data * 1609.34
        else:
            Ubase = data * 0.3048

        if newunit == 'm':  ## CONVERT meters to meters
            data = Ubase
        elif newunit == 'km':
            data = Ubase / 1000.
        elif newunit == 'sm':
            data = Ubase / 1609.34
        else:
            data = Ubase / 0.3048

        Ubase = None
        return data

    def __ConvertTmpUnits(self, data, inunit, newunit):
        """Help on built-in function _ConvertTmpUnits in module MeteoPlot:

                 _ConvertTmpUnits(...)
                     _ConvertTmpUnits(data, inunit, newunit)

            "private" function to convert units of temperature for a dataframe.  Called from the user-callable
            Convert Units function

            Parameters
            -------

            data: numpy array, array containing the data to be converted:
            inunit: string, unit corresponding to the units in the input data
            newunit: string, unit corresponding to the units the data will be converted to


            Returns
            -------
            numpy data array converted from inunit to outunit.
        """

        if inunit == 'F':  ## CONVERT F to K
            Ubase = (data - 32.) * 5. / 9. + 273.15
        elif inunit == 'C':
            Ubase = (data) + 273.15
        else:
            Ubase = data

        ## NOW COVERT KELVIN TO OUTUNIT!
        if newunit == 'F':  ## CONVERT F to K
            data = (Ubase - 273.15) * 9. / 5. + 32.
        elif newunit == 'C':
            data = Ubase - 273.15
        else:
            data = Ubase

        Ubase = None
        return data


    def __ConvertVolConcUnits(self, data, inunit, newunit):
        """Help on built-in function __ConvertVolConcUnits in module MeteoPlot:

                         _ConvertVolConcUnits(...)
                             _ConvertVolConcUnits(data, inunit, newunit)

                    "private" function to convert units of volume concentration for a dataframe.  Called from the user-callable
                    Convert Units function

                    Parameters
                    -------

                    data: numpy array, array containing the data to be converted:
                    inunit: string, unit corresponding to the units in the input data
                    newunit: string, unit corresponding to the units the data will be converted to


                    Returns
                    -------
                    numpy data array converted from inunit to outunit.
                """
        ## CONVERT EVERYTHING TO kg m^-3
        if inunit == 'kg/m^3':
            Ubase = data
        elif inunit == 'ug/m^3':
            Ubase = data / 1E9
        elif inunit == 'g/m^3':
            Ubase = data / 1000.
        elif inunit == 'g/cm^3':
            Ubase = data * 1000.
        elif inunit == 'ug/cm^3':
            Ubase = data * 0.001
        elif inunit == 'mg/cm^3':
            Ubase = data
        elif inunit == 'mg/m^3':
            Ubase = data * 1E-6
        else: ## ppmv assumed.
            Ubase = data * 0.001

        ## NOW COVERT VOLUME IN KG/m^3 TO OUTUNIT!
        if newunit == 'kg/m^3':
            data = Ubase
        elif newunit == 'ug/m^3':
            data = Ubase * 1.E9
        elif newunit == 'g/m^3':
            data = Ubase * 1000.
        elif newunit == 'g/cm^3':
            data = Ubase / 1000.
        elif newunit == 'ug/cm^3':
            data = Ubase / 0.001
        elif newunit == 'mg/cm^3':
             data = Ubase
        elif newunit == 'mg/m^3':
            date = Ubase / 1E-6
        else: ## ppmv assumed.
            data = Ubase / 0.001

        Ubase = None
        return data


    def __ConvertSpdUnits(self, data, inunit, newunit):
        """Help on built-in function _ConvertSpdUnits in module MeteoPlot:

                 _ConvertSpdUnits(...)
                     __ConvertSpdUnits(data, inunit, newunit)

            "private" function to convert units of speed for a dataframe.  Called from the user-callable
            Convert Units function

            Parameters
            -------

            data: numpy array, array containing the data to be converted:
            inunit: string, unit corresponding to the units in the input data
            newunit: string, unit corresponding to the units the data will be converted to


            Returns
            -------
            numpy data array converted from inunit to outunit.
        """

        if inunit == 'm/s':
            Ubase = data
        elif inunit == 'kts':
            Ubase = data / 1.94384
        else:
            Ubase = data / 2.23693

        ## NOW COVERT KELVIN TO OUTUNIT!
        if newunit == 'm/s':  ## CONVERT F to K
            data = Ubase
        elif newunit == 'kts':
            data = Ubase * 1.94384
        else:
            data = Ubase * 2.23693

        Ubase = None
        return data

    def __ConvertPrsUnits(self, data, inunit, newunit):
        """Help on built-in function _ConvertPrsUnits in module MeteoPlot:

                 _ConvertPrsUnits(...)
                     _ConvertPrsUnits(data, inunit, newunit)

            "private" function to convert units of pressure for a dataframe.  Called from the user-callable
            Convert Units function

            Parameters
            -------

            data: numpy array, array containing the data to be converted:
            inunit: string, unit corresponding to the units in the input data
            newunit: string, unit corresponding to the units the data will be converted to


            Returns
            -------
            numpy data array converted from inunit to outunit.
        """
        if inunit == 'inHg':  ## CONVERT F to K
            Ubase = data * 3386.39
        elif inunit == 'hPa':
            Ubase = data * 100.
        else:
            Ubase = data

        ## NOW COVERT KELVIN TO OUTUNIT!
        if newunit == 'inHg':  ## CONVERT F to K
            data = Ubase / 3386.39
        elif newunit == 'hPa':
            data = Ubase / 100.
        else:
            data = Ubase

        Ubase = None
        return data

    def ConvertUnits(self, value, newunit):
        """Help on built-in function ConvertUnits in module MeteoPlot:

                 ConvertUnits(...)
                     ConvertUnits(value, newunit)

            User-callable function to convert units of a variable with the MeteoPlot data dataframe and __UnitDict.


            Parameters
            -------
            value: string, corresponding to the header name for dataframe that will be converted.
                Note value must correspond to one of the following:
                - 'VIS': Visibility (allowed units for VIS are in the VISunitsAllowed dictionary)
                - 'TMP': Temperature (allowed units for TMP are in the TMPunitsAllowed dictionary)
                - 'DEW': Dew Point Temperature (allowed units for DEW are in the TMPunitsAllowed dictionary)
                - 'ALTM' / 'PMSL': Mean Sea level pressure (allowed units for PMSL are in the PRESunitsAllowed dictionary)
                - 'SPD': Wind Speed (allowed units for vis are in the SPDunitsAllowed dictionary)
            newunit: string, unit corresponding to the units the data will be converted to


            Returns
            -------
            nothing, but converts units in the attribute "data" and changes the values of the private "__UnitDict"
            dictionary.

            Example Usage
            -------------
            MeteoPlot.ConvertUnits('TMP','F') # Converts units to degrees Fahrenheit.
        """
        ## VISIBILITY
        if value in ('VIS'):
            if newunit.lower() in self.VISunitsAllowed:
                if newunit.lower() in ('miles', 'mi', 'sm'):
                    newunit = "sm"
                elif newunit.lower() in ('m', 'meters'):
                    newunit = 'm'
                elif newunit.lower() in ('kilometers', 'km'):
                    newunit = 'km'
                else:
                    newunit = 'ft'
                inunit = self.__UnitDict[value]
                if newunit == inunit:
                    return
                else:
                    self.data[value] = self.__ConvertLenUnits(self.data[value].values, inunit, newunit)
            else:
                print("Unit %s not in allowed visibility units: %s" % (newunit, self.VISunitsAllowed))
                return

        ## TEMPERATURE
        if value in ('TMP', 'DEW'):
            if newunit.lower() in self.TMPunitsAllowed:
                if newunit.lower() in ('f', 'fahrenheit'):
                    newunit = "F"
                elif newunit.lower() in ('c', 'celsius'):
                    newunit = 'C'
                else:
                    newunit = 'K'
                inunit = self.__UnitDict[value]
                if newunit == inunit:
                    return
                else:
                    self.data[value] = self.__ConvertTmpUnits(self.data[value].values, inunit, newunit)
            else:
                print("Unit %s not in allowed temperature units: %s" % (newunit, self.TMPunitsAllowed))
                return

        if value in ('ALTM', 'PMSL'):
            if newunit.lower() in self.PRESunitsAllowed:
                if newunit.lower() in ('hpa', 'mb'):
                    newunit = "hPa"
                elif newunit.lower() in ('inhg'):
                    newunit = 'inHg'
                else:
                    newunit = 'Pa'
                inunit = self.__UnitDict[value]
                if newunit == inunit:
                    return
                else:
                    self.data[value] = self.__ConvertPrsUnits(self.data[value].values, inunit, newunit)
            else:
                print("Unit %s not in allowed pressure units: %s" % (newunit, self.PRESunitsAllowed))
                return

        if value in ('SPD'):
            if newunit.lower() in self.SPDunitsAllowed:
                if newunit.lower() in ('ms-1', 'm/s', 'mps'):
                    newunit = "m/s"
                elif newunit.lower() in ('kt', 'kts'):
                    newunit = 'kts'
                else:
                    newunit = 'mph'
                inunit = self.__UnitDict[value]
                if newunit == inunit:
                    return
                else:
                    self.data[value] = self.__ConvertSpdUnits(self.data[value].values, inunit, newunit)
            else:
                print("Unit %s not in allowed wind speed units: %s" % (newunit, self.SPDunitsAllowed))
                return

        if value in ('PM10'):
            if newunit.lower() in self.VOLunitsAllowed:
                if newunit.lower() in ('kg m^-3', 'kg/m3', 'kg/m^3'):
                    newunit = 'kg/m^3'
                elif newunit.lower() in ('g m^-3', 'g/m3', 'g/m^3'):
                    newunit = 'g/m^3'
                elif newunit.lower() in ('g cm^-3', 'g/cm3', 'g/cm^3'):
                    newunit = 'g/cm^3'
                elif newunit.lower() in ('ug m^-3', 'ug/m3', 'ug/m^3'):
                    newunit = 'ug/m^3'
                elif newunit.lower() in ('ug cm^-3', 'ug/cm3', 'ug/cm^3'):
                    newunit = 'ug/cm^3'
                elif newunit.lower() in ('mg cm^-3', 'mg/cm3', 'mg/cm^3'):
                    newunit = 'mg/cm^3'
                elif newunit.lower() in ('mg m^-3', 'mg/m3', 'mg/m^3'):
                    newunit = 'mg/m^3'
                elif newunit.lower() in ('ppmv', 'ppm'):
                    newunit = 'ppm'

                inunit = self.__UnitDict[value]
                if newunit == inunit:
                    return
                else:
                    self.data[value] = self.__ConvertVolConcUnits(self.data[value].values, inunit, newunit)
            else:
                print("Unit %s not in allowed wind speed units: %s" % (newunit, self.VOLunitsAllowed))
                return



        self.__UnitDict[value] = newunit
        self.UnitDict = UD(self.__UnitDict)

    def plot_regular_variable(self,ax,df,unit_dict,vari='TMP',kind='line',ylabel='',**kwargs):
        """Help on built-in function plot_regular_variable in module MeteoPlot:

                 plot_regular_variable(...)
                     plot_regular_variable(ax,df,unit_dict[,vary,kind,ylabel])

            Internal helper function to simplfy plotting regular variables during the add_subplot function.
            Variables here are mostly passed from the add_subplot function to this function which actually performs
            the plotting.

            Parameters
            ----------
            ax: matplotib subplot object, subplot
            df: pandas dataframe, dataframe containing meteorlogical information
            unit_dict: dictionary, unit_dictionary (__UnitDict).

            Returns
            -------
            nothing.
            """

        if kind == 'line':
            ax.plot(df['TIME'],df[vari],**kwargs)
        elif kind == 'bar':
            ax.bar(df['TIME'], df[vari], **kwargs)
        elif kind == 'scatter':
            ax.scatter(df['TIME'], df[vari], **kwargs)
        ax.set_ylabel("%s (%s)"%(ylabel,unit_dict[vari]),fontsize=12)
        ax.set_xlim(np.min(df['TIME']),np.max(df['TIME']))
        return

    def plot_wind(self,ax,df,unit_dict,ylabel='',kind='line'):
        """Help on built-in function plot_wind in module MeteoPlot:

             plot_wind(...)
                 plot_wind(ax,df,unit_dict,[,kind,ylabel])

            Internal helper function to simplfy plotting wind speed and wind vectors on a subplot if
            the "vari" supplied to add_subplot is equal to 'wind'

            Uses the MeteoPlot Attribute: wind_attrs to set wind keywords and modify vectors.

            Parameters
            ----------
            ax: matplotib subplot object, subplot
            df: pandas dataframe, dataframe containing meteorlogical information
            unit_dict: dictionary, unit_dictionary (__UnitDict).

            Returns
            -------
            nothing.
        """
        if kind == 'line':
            ax.plot(df['TIME'],df[self.wind_attrs['spd_name']],**self.wind_attrs['spd_kwargs'])
        elif kind == 'bar':
            ax.bar(df['TIME'],df[self.wind_attrs['spd_name']],**self.wind_attrs['spd_kwargs'])
        elif kind == 'scatter':
            ax.scatter(df['TIME'], df[self.wind_attrs['spd_name']], **self.wind_attrs['spd_kwargs'])

        ax.set_ylabel("%s (%s)"%(ylabel,unit_dict[self.wind_attrs['spd_name']]),fontsize=12)
        ax.set_xlim(np.min(df['TIME']),np.max(df['TIME']))

        if self.wind_attrs['barbs'] == True:
            u,v=pc.convert_met_2_uv(df[self.wind_attrs['dir_name']],df[self.wind_attrs['spd_name']]) # Convert speed and direction to u and v
            u,v=np.ma.masked_invalid(u),np.ma.masked_invalid(v)

            ax.barbs(df['TIME'][::self.wind_attrs['skip']].values,
                     np.nanmax(df[self.wind_attrs['spd_name']].values)*np.ones_like(df[self.wind_attrs['spd_name']][::self.wind_attrs['skip']])
             ,u[::self.wind_attrs['skip']]*2.28,v[::self.wind_attrs['skip']]*2.28,**self.wind_attrs['barb_kwargs'])

    def plot_weather_types(self,ax, df, weather_keys, color_dict):
        """Help on built-in function plot_weather_types in module MeteoPlot:

             plot_weather_types(...)
                 plot_weather_types(ax,df,weather_keys,color_dict)

            Internal helper function to simplfy plotting METAR weather reports on a subplot if
            the "vari" supplied to add_subplot is equal to 'weather' or 'wx'

            Parameters
            ----------
            ax: matplotib subplot object, subplot
            df: pandas dataframe, dataframe containing meteorlogical information
            weather_keys: list, list of METAR weather codes (e.g., ['BR','RA'])
            color_dict: dictionary, dictionary with METAR keys and string color values (e.g., {'BR':'blue','RA':'green'})

            Returns
            -------
            nothing.
        """
        weather_key_dict = {}
        c = 2
        for i in weather_keys:
            weather_key_dict[i] = c
            c += 1

        weather_data = df['Weather']
        weather_data = np.ma.masked_where(weather_data == 'NaN', weather_data)
        weather_values = []
        for i in weather_data:
            try:  ##sort of a roundabout bassackwards way, but it works I suppose?
                check = np.isnan(i)
                weather_values.append(-1)
            except:
                weather_values.append(i)

        for idx, i in enumerate(weather_values[:-1]):
            dt = df['TIME'][idx + 1] - df['TIME'][idx]
            if i == -1:
                continue  ## don't do anything -> go to next cycle.
            wxdata = i.split()  ## incase there are more than 1 weather values.
            for z in wxdata:  ## loop through all weather reports here.
                mark_idx = 1
                if '-' in z:
                    mark_idx = 0
                elif '+' in z:
                    mark_idx = 2
                for k in weather_keys:
                    if k in z:
                        if 'BLSN' in z and k == 'SN':
                            continue
                        ax.fill_between([df['TIME'][idx] - dt / 2., df['TIME'][idx + 1] - dt / 2.],
                                        weather_key_dict[k] - 0.35, weather_key_dict[k] + 0.35
                                        , color=color_dict[k][mark_idx])
        ax.set_ylim(1, len(weather_key_dict) + 2)
        ax.set_yticks([weather_key_dict[i] for i in weather_key_dict])
        ax.set_yticklabels([i for i in weather_key_dict])
        ax.set_xlim(np.min(df['TIME']), np.max(df['TIME']))
        ax.set_ylabel("Weather", fontsize=12)
        return

class UD(dict):
    """This is a constructor to help manage the UnitDict attributed of the MeteoPlot function.
        It's key purpose is to set the functions that access the hidden __UnitDict dictionary and preclude users from
        editing the Unit Dictionary manually, such that the units can only be edited in the 'ConvertUnits' function along
        with the actual data such that the units listed in the unit dictionary are always correct."""
    def __getitem__(self,key):
        return dict.__getitem__(self,key)
    def __setitem__(self, key, value):
        print("You can't set these keys manually, please use the ConvertUnits Function")
        return
    def __delitem__(self, key):
        dict.__delitem__(self,key)
    def __iter__(self):
        return dict.__iter__(self)
    def __len__(self):
        return dict.__len__(self)
    def __contains__(self, x):
        return dict.__contains__(self,x)