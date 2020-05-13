# coding=utf-8
"""This script calls a "namelist" and post processes UK Met Office (UKMO) netcdf
   files into a single netcdf file with user specified variables.

   Note: This script is very specifically designed to handle surface data related to dust variables.
   It is designed to take data and output it into a format common to the WRF-Chem output such that it can be
   quickly and easily compared using the plotting functions within this larger package.

   Coding format follows PEP8 (ish) rules, but not perfectly.
"""
import numpy as np
from matplotlib import pyplot as plt

from cartopy import crs as ccrs
from cartopy import feature as cfeature
from netCDF4 import Dataset

from datetime import datetime, timedelta
from ERDCMOOSE.plotting import plotCommon as pc
import glob as glob
import sys


class StandardPlot:
    """ Class Constructor for the ERDC-MOOSE standard plotting object contained within GeoPlot."""
    def __init__(self):
        """ This function defines a standard plotting function to plot data from a netCDF file."""
        try:
            ## TRY TO IMPORT LOCAL STAGE 4 DATA SET!
            import pygrib as pg
            self.stg4_special=['stg4_pcp']
        except:
            print("Either I Cannot find stage 4 input support, OR, I cannot find pygrib")
            print("This means, no stage 4 precipitation plotting is possible")
            print("I will let you run everything else though...")
            self.stg4_special=[]

        self.namelist_dictionary={'geo_box':[-125,-65,20,60],'vars':['10m Wind speed','Mean sea level pressure'],
                                  'levels':[500,1000],'var_range':[[25,80],[20,80]],'over_land':1,
                                  'time_step':'2019-02-23 18:00:00','time_format':'%Y-%m-%d %H:%M:%S',
                                  'projection':'LambertConformal','proj_kwargs':{},
                                  'datapath':'/p/work/tletcher/WRF_POST/','smooth_field':[-1,2],
                                  'colorbar':True,'colormap':'jet','plot_winds':[True,False],
                                  'wind_kwargs':{},'barb_skip':[5],'vector_type':'barb',
                                  'contour_kwargs1':{},'contour_kwargs2':{},'contour_kwargs3':{},
                                  'print_vars':True,'pfx':'MARCH','figsize':'(10,6)','mslp_hilo':True,
                                  'narr_folder':'./','shade_type':'pocolormesh','linspace':[True,False],
                                  'extend':'neither','show_title':True,'valid_title':False,'convert_temp':['None'],
                                  'shade_first':True,'stg4_folder':'/p/work/tletcher/STAGE4_DATA/20190224/',
                                  'stg4_inc':'06h','pcp_accum_times':'2019-02-23 18:00:00','clabels':[False,True],
                                  'clabfmt':['%i','%i'],'clabbox':[False,False],
                                  'features':['states','lakes','borders'],'map_resolution':'50m'}


        ## Define allowable keys for certain namelist options!
        self.allowable_temp_conversions = ['None', 'C', 'F', 'K']
        self.allowable_extends = ["neither", "both", "min", "max"]
        self.allowable_projections=['LambertConformal','PlateCarree']
        self.allowable_shade_types=['contourf','pcolormesh']
        self.allowable_map_resolutions = ['110m','50m','10m']

        # Here are a few initial fields that are quasi-hard coded here.
        # These correspond to variable names for lat,lon,u,v, etc, for geographic coordinates, and
        # for special variables.  Note that if you use the post processing parts of this package
        # to post process the ERA or WRF data, they will be saved to netCDF files with the following names.

        self.lonName='Longitude'
        self.latName='Latitude'
        self.projection=ccrs.PlateCarree()

        self.U10='10 metre U wind component'
        self.V10='10 metre V wind component'
        self.U='U velocity'
        self.V='V velocity'


    def KtoC(self,tmp):
        """Help on built-in function KtoC in module GeoPlot:

        KtoC(...)
            KtoC(TK)

            Returns temperature in degrees Celsius given an input temperature Kelvin.

            Assumes freezing point is 273.15 Kelvin.

            Parameters
            ----------
            TK : number, Temperature in Kelvin

            Returns
            -------
            Temperature in Celsius
        """
        return tmp - 273.16

    def KtoF(self,tmp):
        """Help on built-in function KtoF in module GeoPlot:

        KtoF(...)
            KtoF(TK)

            Returns temperature in degrees Fahrenheit given an input temperature Kelvin.

            Assumes freezing point is 273.15 Kelvin.

            Parameters
            ----------
            TK : number, Temperature in Kelvin

            Returns
            -------
            Temperature in Fahrenheit
        """
        tmp = tmp - 273.16
        return 9. / 5. * tmp + 32.

    def print_namelist(self):
        """Help on built-in function print_namelist in module GeoPlot:

        print_namelist(...)
            print_namelist()

            Prints out the keys currently listed in the namelist_dictionary with currently set values.


            Returns
            -------
            nothing.
        """

        print("--------------------------------------------")
        print("|         Printing Namelist Values         |")
        print("--------------------------------------------")
        for n in self.namelist_dictionary:
            if type(self.namelist_dictionary[n]) is list:
                print("- %s" % n)
                for i in self.namelist_dictionary[n]:
                    print("  - %s" % i)
            else:
                print("- %s : %s" % (n, self.namelist_dictionary[n]))

    def parse_namelist(self,namelist_file='Reanalysis_namelist.txt'):
        """Help on built-in function parse_namelist in module GeoPlot:

        parse_namelist(...)
            parse_namelist(namelist_file)

            Resets the GeoPlot object's namelist_dictionary from values parsed out of a user input text file

            In general, the namelist is organized using comma and colon separation to guide parsing function
            commas typically separate different variable options, and colons are used to separate keyword arguments and
            list object values.  See Letcher et al., 2020 for a detailed description of how GeoPlot works with namelists.

            Parameters
            ----------
            namelist_file : filepath, Path to a text file with namelist information

            Returns
            -------
            Nothing, but resets an intrinsic GeoPlot dictionary.
        """


        import os
        cwd = os.getcwd()

        lookforkeys = ['geo_box', 'vars', 'levels', 'var_range', 'over_land', 'time_step',
                       'data_path', 'projection', 'proj_kwargs', 'time_format', 'smooth_field', 'colorbar',
                       'colormap', 'plot_winds', 'wind_kwargs', 'barb_skip', 'vector_type',
                       'contour_kwargs1','contour_kwargs2','contour_kwargs3',
                       'print_vars', 'pfx', 'figsize', 'print_opts',
                       'pcp_accum_times', 'stg4_folder', 'stg4_inc', 'narr_folder', 'shade_type',
                       'track_low', 'track_int', 'show_title', 'valid_title', 'linspace',
                       'convert_temp', 'extend', 'clabels', 'clabfmt', 'clabbox', 'shade_first','features','mslp_hilo']

        try:
            namelistfile = glob.glob(cwd + '/' + namelist_file)[0]
        except:
            self.print_fatal("NAMELIST FILE: %s does not exist, check your paths are try again!")


        list_keys = ['vars', 'smooth_field', 'plot_winds', 'levels', 'geo_box', 'var_range',
                     'pcp_accum_times', 'pfx', 'linspace', 'convert_temp',
                     'clabels', 'clabfmt', 'clabbox','features']

        bool_keys=['colorbar','plot_winds','print_vars','mslp_hilo','linspace','shade_first','clabels','clabbox']
        ## OKAY, DEFAULTS SET! ##
        with open(namelistfile) as fp:
            for line in fp:
                var_key = line.split('=')[0].strip(' \t\n\r')
                if var_key in lookforkeys and '#' not in var_key:  ## if not a comment!
                    values = line.split('=')[1].split('#')[0].strip(' \t\n\r').split(',')
                    if 'kwargs' in var_key:  ## Special case, key word argument dictionary
                        contour_plots = line.split('=')[1].split('#')[0].strip(' \t\n\r').split('|')
                        for cdx, cc in enumerate(contour_plots):
                            self.namelist_dictionary[var_key][cdx] = {}
                            values = cc.split(",")
                            arg_keys = [k.split(':')[0] for k in values]
                            arg_vals = [k.split(':')[1] for k in values]
                            for adx, a in enumerate(arg_keys):
                                try:
                                    self.namelist_dictionary[var_key][cdx][a] = float(arg_vals[adx])
                                except:
                                    self.namelist_dictionary[var_key][cdx][a] = arg_vals[adx]
                    elif var_key == 'var_range':
                        ## Special case 2: variable range nested list!
                        values = [[float(y) for y in x.split(':')] for x in values if x != '']
                        for vdx, v in enumerate(values):
                            if len(v) == 1:
                                self.print_fatal([
                                                "YOU HAVE ONLY CHOSEN A SINGLE VALUE FOR YOUR VARIABLE RANGE FOR ONE OF YOUR VARIABLES",
                                                "I REQUIRE AT LEAST 2 VALUES, A 3rd VALUE MAY BE SUPPLIED TO INDICATE THE INTERVAL!",
                                                "PLEASE CHECK YOUR NAMELIST AND TRY AGAIN!"])
                            if len(v) == 2:
                                values[vdx] = v + [41]  ## add "41" to values to set contour interval.

                        self.namelist_dictionary[var_key] = values
                    else:
                        values = [x for x in values if x != '']
                        if var_key == 'figsize':
                            ## ANOTHER SPECIAL CASE: FIGSIZE GETS COVERTED TO A TUPLE!
                            self.namelist_dictionary[var_key] = (float(values[0]), float(values[1]))
                            continue
                        if var_key == 'features':
                            ## ANOTHER SPECIAL CASE: CONVERT ALL FEATURES TO LOWECASE!
                            self.namelist_dictionary[var_key]=[i.lower() for i in values]
                            continue
                        if len(values) > 1:
                            if var_key in list_keys:
                                try:
                                    self.namelist_dictionary[var_key] = [float(i) for i in values]
                                except:
                                    self.namelist_dictionary[var_key] = values

                                if var_key in bool_keys:
                                    self.namelist_dictionary[var_key]=\
                                        [bool(int(i)) for i in self.namelist_dictionary[var_key]]
                            else:
                                self.print_warning(["THIS IS NOT A LIST VARIABLE",
                                               "SETTING KEY VALUE TO FIRST OPTION"])
                                try:
                                    self.namelist_dictionary[var_key] = float(values[0])
                                except:
                                    self.namelist_dictionary[var_key] = values[0]

                                if var_key in bool_keys:
                                    self.namelist_dictionary[var_key]=bool(int(self.namelist_dictionary[var_key]))
                        else:
                            if var_key in list_keys:
                                try:
                                    self.namelist_dictionary[var_key] = [float(i) for i in values]
                                except:
                                    self.namelist_dictionary[var_key] = values

                                if var_key in bool_keys:
                                    self.namelist_dictionary[var_key] = \
                                        [bool(int(i)) for i in self.namelist_dictionary[var_key]]
                            else:
                                try:
                                    self.namelist_dictionary[var_key] = float(values[0])
                                except:
                                    self.namelist_dictionary[var_key] = values[0]

                                if var_key in bool_keys:
                                    self.namelist_dictionary[var_key] = bool(int(self.namelist_dictionary[var_key]))

    def print_warning(self,message):

        """Help on built-in function print_warning in module GeoPlot:

        print_warning(...)
            print_warning()

            Simple helper function to decorate text surrounding a warning message.

            Returns
            -------
            nothing.
        """
        print("-----------------")
        print("|    WARNING    |")
        if type(message) is not list:
            print(message)
        else:
            for m in message:
                print(m)
        print("-----------------")

    def print_fatal(self,message):

        """Help on built-in function print_fatal in module GeoPlot:

        print_fatal(...)
            print_fatal()

            Simple helper function to decorate text surrounding a fatal error message.

            Returns
            -------
            nothing.
        """
        print("-----------------")
        print("|  FATAL ERROR  |")
        if type(message) is not list:
            print(message)
        else:
            for m in message:
                print(m)
        print("-----------------")
        sys.exit(1)

    def define_projection(self):
        """Help on built-in function define_projection in module GeoPlot:

        define_projection(...)
            define_projection()

            A function to read in Cartopy projection related arguments attached to the GeoPlot namelist_dictionary
            and uses them to build a Cartopy projection object.  ERDC-MOOSE V 1.0 only allows
            Lambert Conformal and Cylindrical map projections.

            By Default a cylindrical projection is used, a Lambert projection can be defined if
            namelist_dictionary['projection'] is equal to one of the following (case-insensitive):
                - lcc
                - lambert
                - lambertconformal

            Returns
            -------
            nothing, but modifies the existing GeoPlot "projection" attribute

            Example Usage
            -------------
            Gplot=GeoPlot()
            ....
            Gplot=define_projection()
            ax=plt.subplot(111,projection=Gplot.projection)
            ....
        """
        try:
            proj_kwargs=self.namelist_dictionary['proj_kwargs'][0]
        except:
         proj_kwargs=self.namelist_dictionary['proj_kwargs']
        print(proj_kwargs)

        if self.namelist_dictionary['projection'].lower() in ['lcc','lambert','lambertconformal']:
            self.projection=ccrs.LambertConformal(**proj_kwargs)
        else:
            self.projection=ccrs.PlateCarree(**proj_kwargs)


    def CheckNamelist(self):
        """Help on built-in function CheckNamelist in module GeoPlot:

        CheckNamelist(...)
            CheckNamelist()

            Simple helper function to make sure namelist values contained within the namelist_dictionary are valid

            Returns
            -------
            nothing.
        """

        ## CHECK TEMPERATURE CONVERSIONS! ##
        for idx, i in enumerate(self.namelist_dictionary['convert_temp']):
            if 'temperature' in self.namelist_dictionary['vars'] and i not in self.allowable_temp_conversions:
                self.print_warning(["TEMPERATURE CONVERSION CHOICE %s NOT ALLOWED" % i,
                               'ALLOWED CONVERSIONS:'] + ["   --%s" % j for j in self.allowable_temp_conversions] +
                              ['DEFAULTING TO "None"'])
                self.namelist_dictionary['convert_temp'][idx] = 'None'

        ## A FEW ADDITIIONAL CHECKS! ##
        if len(self.namelist_dictionary['linspace']) < len(self.namelist_dictionary['vars']):
            self.print_fatal(["LINESPACE LIST IS NOT THE SAME SIZE AS THE NUMBER OF VARIABLES",
                         "PLEASE CHECK YOUR NAMELIST AND TRY AGAIN!"])

        if len(self.namelist_dictionary['var_range']) < len(self.namelist_dictionary['var_range']):
            self.print_fatal(["VAR_RANGE LIST IS NOT THE SAME SIZE AS THE NUMBER OF VARIABLES",
                         "PLEASE CHECK YOUR NAMELIST AND TRY AGAIN!"])

        if self.namelist_dictionary['shade_type'] not in self.allowable_shade_types:
            self.print_warning(["%s is not an allowable shade type" % self.namelist_dictionary['shade_type'],
                           "Defaulting to 'pcolormesh'"])
            self.namelist_dictionary['shade_type'] = 'pcolormesh'

        ## CHECK CONTOURF EXTEND INPUT! ##
        if self.namelist_dictionary['shade_type'] == 'contourf' and self.namelist_dictionary['extend'] not in self.allowable_extends:
            self.print_warning(["CONTOURF EXTENSION CHOICE %s NOT ALLOWED" % i,
                           'ALLOWED EXTENSION OPTIONS:'] + ["   --%s" % j for j in self.allowable_extends] +
                          ['DEFAULTING TO "neither"'])
            self.namelist_dictionary['extend'] = 'neither'


        if self.namelist_dictionary['map_resolution'] not in self.allowable_map_resolutions:
            self.print_warning(['MAP RESOLUTION CHOICE %s NOT ALLOWED!'%i,
                                'ALLOWED MAP RESOLUTIONS:'] + ["   --%s" % j for j in self.allowable_map_resolutions] +
                          ['DEFAULTING TO "50m"'])
            self.namelist_dictionary['map_resolution'] = '50m'

    def makePlot(self,ax,fix_lon_pos=True):

        """Help on built-in function makePlot in module GeoPlot:

        makePlot(...)
            makePlot(ax,[, fix_lon_pos,])

            This function handles the majority of the file IO, and plotting including:
                a) Reading in the NetCDF file, choosing the approrite vertical level and time step as defined in
                   GeoPlot namelist_dictionary.
                b) Making any required unit conversions for temperature, and merging U/V wind components into speed
                c) Contour filling and plotting and geographic overlay
                d) Any automatic titling or colorbars

            Parameters
            ----------
            ax : matplotlib subplot, a defined matplotlib subplot object with a cartopy projection.
            fix_lon_pos : Boolean (optional), A flag that will adjust longitude to run from 0-360 instead of -180 to 180
                useful for domains that straddle the date line.

            Notes
            ----------
            Assumes data is in NetCDF format

            Assumes that the latitude / longitude variables are in decimal degrees and
            that the NetCDF variable names are set to the GeoPlot attributes: latName and lonName.

            Assumes that the time coordinate information is stored in the "Time" variable within the NetCDF file.
            If the "Time" variable has the attribute "calendar" the time coordinate information will be determined using
            num2date, otherwise it will use datetime.  All time information is converted to datetime.

            Can convert the units of temperature only, if the namelist variable "convert_temp" is set to "C" or "F"

            Accumulated precpitation variables are accumuated over the time-period specified by the
            namelist_dictionary['pcp_accum_times'] variable.

            Mean Sea Level Pressure is special and largly ignores namelist contour attributes.

            Returns
            -------
            nothing, but it adds data to a pre-existing subplot, which can be displayed or saved to an image file

            Example Usage
            -------------
            Gplot=GeoPlot()
            ....
            Gplot=define_projection()
            ax=plt.subplot(111,projection=Gplot.projection)
            Gplot.makePlot(ax)
            plt.show()
            ....
        """

        from netCDF4 import Dataset,num2date
        from datetime import datetime, timedelta

        self.CheckNamelist()
        title=''
        show_hls=True
        var_range=self.namelist_dictionary['var_range']
        if type(self.namelist_dictionary['pfx']) == str:
            ##critical check, if the namelist prefix is a string, then make it a list
            ##so it behaves correctly in the below loop.
            self.namelist_dictionary['pfx']=[self.namelist_dictionary['pfx']]
        for pdx, p in enumerate(self.namelist_dictionary['pfx']):
            c_count = 1
            w_count = 0

            wskp = int(self.namelist_dictionary['barb_skip'])
            barb_kwargs=self.namelist_dictionary['wind_kwargs'][0]
            if pdx == 1:
                print("Resetting a couple of things to avoid overlap...")
                self.namelist_dictionary['colorbar'] = False
                show_hls = False

            filepath = '%s/*%s*.nc' % (self.namelist_dictionary['data_path'], p)

            ## Check if file exists.
            try:
                datafile = glob.glob(filepath)[0]
                print("Opening %s" % datafile)
            except:
                self.print_fatal("NO Files with path %s found.... exiting." % filepath)

            ## Open NETCDF file and get times, lat/lon
            ncdata = Dataset(datafile, 'r')  ## open for appending,
            ## DO NOT NEED TO DO THE DESCRIPTION STUFF! or dimenion setting! ##
            ## CHECK IF WRF FILE! ##
            wrf_special = []
            if 'TITLE' in ncdata.ncattrs() and 'WRF' in ncdata.TITLE:
                pcp_varnames = ['Gridscale Precipitation rate', 'Gridscale Snowfall rate']
                if all(item in ncdata.variables for item in pcp_varnames) == True:
                    wrf_special = wrf_special + ['WRF Accumulated Precipitation', 'WRF Accumulated Snow']
                if 'Convective Precipitation rate' in ncdata.variables:
                    ## Set flag to include convective precip in wrf accumulated precip variable
                    ## Doesn't matter if cumulus parameterization is off in model.
                    wrf_convective_pcp = True
                else:
                    wrf_convective_pcp = False

            try:
                times = num2date(ncdata.variables['Time'][:],units=ncdata.variables['Time'].units,
                                                             calendar=ncdata.variables['Time'].calendar)
                times = [datetime.strptime(datetime.strftime(i, '%Y-%m-%d %H:%M:%S'), '%Y-%m-%d %H:%M:%S')
                             for i in times]
            except:
                times = [datetime.strptime(ncdata.variables['Time'].units, ncdata.variables['Time'].units) +
                         timedelta(hours=float(ss)) for ss in ncdata.variables['Time'][:]]

            time_choice_strp = datetime.strptime(self.namelist_dictionary['time_step'],
                                                 self.namelist_dictionary['time_format'])
            time_indx = np.argmin(np.abs(np.array(times) - time_choice_strp))

            try:
                levels = ncdata.variables['Levels'][:].squeeze()
            except:
                print('WARNING, I CANNOT FIND THE VARIABLE "LEVELS" in the netcdf FILE')
                print("Assuming this file doesn't have any levels, assuming this is a surface file.")
                print("Setting dummy array to levels=[1000]")
                levels=np.array([1000])
            lons, lats = ncdata.variables[self.lonName][:].squeeze(), ncdata.variables[self.latName][:].squeeze()
            if fix_lon_pos == True:
                lons = np.ma.masked_less(lons, 0.0).filled(360. + lons)

            lev_idx = []
            for vv in self.namelist_dictionary['levels']:
                if vv not in levels:
                    self.print_warning(["Variable levels is not in dataset",
                                   "Defaulting to nearest vertical level."])
                lev_idx.append(np.argmin(np.abs(levels - vv)))

            ##Build list of avaiable variables from netCDF file.
            skp_list = ['Time', self.lonName, self.latName, 'Levels']
            special = ['Wind speed', '10m Wind speed', 'NARR_pcp'] + self.stg4_special + wrf_special
            var_list = []
            for vdx, v in enumerate(ncdata.variables):
                if vdx == 0 and self.namelist_dictionary['print_vars'] == True:
                    print("--------------------------")
                    print("Printing avaialble variables")
                    print("--------------------------")
                if v not in skp_list:
                    if self.namelist_dictionary['print_vars'] == True:

                    var_list.append(v)
            if self.namelist_dictionary['print_vars'] == True:
                print("--------------------------")

            if self.namelist_dictionary['print_vars'] == True:
                print("--------------------------")
                print("Printing avaiable special variables")
                print("--------------------------")
                for v in special:
                    print(" -- %s " % (v))

            pc.AddMap(ax, extent=self.namelist_dictionary['geo_box'],
                       features=self.namelist_dictionary['features'],
                       reso=self.namelist_dictionary['map_resolution'])
            ## initialize Z order
            if self.namelist_dictionary['over_land'] == True:
                zorder = 3
            else:
                zorder = 2


            for vdx,vari in enumerate(self.namelist_dictionary['vars']):
                print("---------------------------")
                print("Working on Variable %s"%vari)
                print("---------------------------")
                if self.namelist_dictionary['plot_winds'][vdx] == True:
                    if vari == '10m Wind speed':
                        UU=ncdata.variables['10 metre U wind component'][time_indx,:]
                        VV=ncdata.variables['10 metre V wind component'][time_indx,:]
                    else:
                        print("Collecting Wind data for level %i hPa"%self.namelist_dictionary['levels'][vdx])
                        if self.namelist_dictionary['levels'][vdx] == 1000:
                            print("Assuming winds at 1000 hPa are 10m winds....")
                            UU=ncdata.variables[self.U10][time_indx,:]
                            VV=ncdata.variables[self.V10][time_indx,:]
                        else:
                            UU=ncdata.variables[self.U][time_indx,lev_idx[vdx],:]
                            VV=ncdata.variables[self.V][time_indx,lev_idx[vdx],:]

                if vari not in var_list and vari not in special:
                    print("Your Variable is not in the list! Defaulting to 2M Temperature")
                    vari='2 metre temperature'


                # THIS NEXT SOMEWHAT UNWIELDY LOOKING BLOCK OF CODE ACTS TO PLOT "Special Variables"
                # that can't just be pulled directly from the netCDF file.
                # this includes catches for both observational precipitation from (stage 4 data set)
                # and NARR precpitation.
                # Special variables:
                #   - Wind Speed
                #   - 10m Wind Speed
                #   - stg4_pcp
                #   - NARR_pcp
                #   - WRF_pcp
                if vari == 'Wind speed':
                    # Special case here, if we want to plot windspeed, but we have already loaded
                    # u and v wind components at this level from first step, don't take the time and effort
                    # to reload them just for the sake of slightly cleaner code.
                    if self.namelist_dictionary['plot_winds'] == True:
                        ##Don't reload the data if we don't have to...
                        var_data=np.sqrt(UU**2.+ VV**2.)
                    else:
                        ## otherwise, you gotta load the data.
                        var_data=np.sqrt(ncdata.variables[self.U][time_indx,lev_idx[vdx],:]**2.+ \
                                 ncdata.variables[self.V][time_indx,lev_idx[vdx],:]**2.)
                    units='m s^-1'
                    if vdx == 0:
                        title='%i hPa Wind Speed'%self.namelist_dictionary['levels'][vdx]
                    else:
                        title+='\n %i hPa Wind Speed'%self.namelist_dictionary['levels'][vdx]
                elif vari == '10m Wind speed':
                    var_data=np.sqrt(ncdata.variables[self.U10][time_indx,:]**2.+ \
                             ncdata.variables[self.V10][time_indx,:]**2.)
                    units='m s^-1'
                    if vdx == 0:
                        title='10m Wind Speed'
                    else:
                        title+='\n 10m Wind Speed'

                ## SPECIAL BLOCK OF SPECIAL VALUES -- > PRECIPITATION

                elif vari == 'stg4_pcp': ## "Observed" Stage 4 precip
                    lons,lats,var_data,title,units=\
                        pc.acc_precip(self.namelist_dictionary,vdx=vdx,src_name='Stage 4',title=title)

                elif vari == 'NARR_pcp': # "Observed" NARR Precip
                    lons,lats,var_data,title,units=\
                        pc.acc_precip(self.namelist_dictionary,vdx=0,src_name='NARR',title=title)
                elif vari == 'WRF Accumulated Precipitation': ## WRF total Precipitation
                    var_data, title, units = pc.wrf_acc_pcp(self.namelist_dictionary,times,time_indx,ncdata,
                                                            vdx=vdx,title=title)

                elif vari == 'WRF Accumulated Snow':

                    var_data, title, units = pc.wrf_acc_snow(self.namelist_dictionary, times, time_indx, ncdata,
                                                            vdx=vdx, title=title)
                else:
                    units=ncdata.variables[vari].units
                    ## NOT USING A SPECIAL VARIABLE! ###
                    lons,lats=ncdata.variables['Longitude'][:].squeeze(),ncdata.variables['Latitude'][:].squeeze()
                    ## IF NOT A SPECIAL VARIABLE !!#
                    if 'z' not in ncdata.variables[vari].dimensions:
                        #This is a single level (surface) variable, don't use L_index!
                        if np.ndim(ncdata.variables[vari][:]) < 3: ## THEN THIS IS A Time-Invarient Variable!
                            var_data = ncdata.variables[vari][:]
                        else:
                            var_data=ncdata.variables[vari][time_indx,:]
                    else: #otherwise, use the level index.
                        if np.ndim(ncdata.variables[vari][:]) < 4:  ## THEN THIS IS A Time-Invarient Variable!
                            var_data = ncdata.variables[vari][lev_idx[vdx], :]
                        else:
                            var_data=ncdata.variables[vari][time_indx,lev_idx[vdx],:]
                        #print(np.shape(var_data),np.nanmax(var_data),np.nanmin(var_data))

                    ## SPECIAL CASE --> CONVERT TEMPERATURE!
                    if 'temperature' in vari.lower():
                        if self.namelist_dictionary['convert_temp'][vdx] == 'F':
                            print("Converting Temperature to F")
                            var_data=pc.KtoF(var_data)
                            units='F'
                        elif self.namelist_dictionary['convert_temp'][vdx] == 'C':
                            print("Converting Temperature to C")
                            var_data=pc.KtoC(var_data)
                            units='C'
                    if vari != 'Mean sea level pressure':
                        if vdx == 0:
                            title=ncdata.variables[vari].longname+' [%s]'%units
                        else:
                            title='\n %s'%ncdata.variables[vari].longname+' [%s]'%units

                    units=ncdata.variables[vari].units
                if self.namelist_dictionary['smooth_field'][vdx] > 0:
                    print("Smoothing data field with sigma of %i"%self.namelist_dictionary['smooth_field'][vdx])
                    var_data=pc.smooth_field(var_data,self.namelist_dictionary['smooth_field'][vdx])

                if self.namelist_dictionary['mslp_hilo'] == True and vari=='Mean sea level pressure' and pdx==0:
                    print("Special Plot: MSLP with hi's and lo's")
                    mslp_title=pc.plot_mslp(ax,lons,lats,var_data/100.,times[time_indx],hsize=50*wskp/1.8,lsize=25*wskp/1.8,zorder=6,hl=show_hls)
                    units='hPa'
                    if vdx == 0:
                        title=ncdata.variables[vari].longname+' [%s]'%units
                    else:
                        title+='\n %s'%ncdata.variables[vari].longname+' [%s]'%units
                else:
                    if vdx == 0 and self.namelist_dictionary['shade_first'] == True:
                        ## SHADE FIRST VARIABLE! ##
                        title+=' (shaded)'
                        if self.namelist_dictionary['linspace'][vdx] == True:
                            clevels=np.linspace(var_range[vdx][0],var_range[vdx][1],var_range[vdx][2])
                        else:
                            clevels=np.arange(var_range[vdx][0],var_range[vdx][1],var_range[vdx][2])
                        if self.namelist_dictionary['shade_type'] == 'contourf':
                            self.im=plt.contourf(lons,lats,var_data,cmap=self.namelist_dictionary['colormap'],zorder=zorder,levels=clevels,
                                transform=ccrs.PlateCarree(),extend=self.namelist_dictionary['extend'])
                        else:
                            self.im=plt.pcolormesh(lons,lats,np.ma.masked_less(var_data,var_range[vdx][0]),cmap=self.namelist_dictionary['colormap'],
                                zorder=zorder,vmin=var_range[vdx][0],vmax=var_range[vdx][1],transform=ccrs.PlateCarree())
                        if self.namelist_dictionary['colorbar'] == True:
                            cbar=plt.colorbar(self.im,fraction=0.046, pad=0.04)
                            cbar.set_label('[%s]'%units,fontsize=12.)
                    else:
                        if pdx == 0:
                            if self.namelist_dictionary['linspace'][vdx] == True:
                                clevels=np.linspace(var_range[vdx][0],var_range[vdx][1],var_range[vdx][2])
                            else:
                                clevels=np.arange(var_range[vdx][0],var_range[vdx][1],var_range[vdx][2])

                            contour_kwargs=self.namelist_dictionary['contour_kwargs%i'%c_count][0]
                            print(contour_kwargs)
                            C=plt.contour(lons,lats,var_data,zorder=zorder,levels=clevels,transform=ccrs.PlateCarree(),
                                          **contour_kwargs)
                            c_count=c_count+1 ## Move on to next variable --> ONLY IF ENOUGH ARGUMENTS! ##
                            if self.namelist_dictionary['clabels'][vdx] == True:
                                clabels=plt.clabel(C,fmt=self.namelist_dictionary['clabfmt'][vdx])
                                if self.namelist_dictionary['clabbox'][vdx] == True:
                                    for t in clabels:
                                        t.set_bbox({'facecolor': 'white', 'pad': 4})
                                        t.set_fontweight('heavy')
                        title+=' (contoured)'

                    if self.namelist_dictionary['plot_winds'][vdx] == True and pdx == 0:
                        lons,lats=ncdata.variables[self.lonName][:].squeeze(),ncdata.variables[self.latName][:].squeeze()
                        if self.namelist_dictionary['vector_type']=='barb':
                            plt.barbs(lons[::wskp,::wskp],lats[::wskp,::wskp],UU[::wskp,::wskp],VV[::wskp,::wskp],
                                     transform=ccrs.PlateCarree(),zorder=zorder+1,**barb_kwargs)
                        elif self.namelist_dictionary['vector_type']=='quiver':
                            plt.quiver(lons[::wskp,::wskp],lats[::wskp,::wskp],UU[::wskp,::wskp],VV[::wskp,::wskp],
                                       transform=ccrs.PlateCarree(),zorder=zorder+1,**barb_kwargs)
                        else:
                            plt.streamplot(lons, lats, UU, VV, transform=ccrs.PlateCarree(),zorder=zorder+1,**barb_kwargs)

                zorder=zorder+1

        self.current_title=title
        if self.namelist_dictionary['show_title'] == True:
            plt.title(title, loc='left',fontweight='bold')
        if self.namelist_dictionary['valid_title'] == True:
            plt.title('VALID: %s'%times[time_indx].strftime('%c'), loc='right')
