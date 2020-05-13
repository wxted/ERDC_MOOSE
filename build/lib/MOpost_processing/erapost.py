# coding=utf-8
"""This script calls a "namelist" and post processes the grib2 ERA data downloaded from
   the NCAR RDA. These files come in a somewhat unwieldy format with multiple files containing different variables
   and model levels.  Class is similar to WrfPost and allows for the creation and storage of these data into a
   simpler to read netCDF file.
   files into a single netcdf file with user specified variables.

   I added "appending" options make it easy to add a variable to an existing
   file so that you don't need to recompute the full 3D atmosphere everytime you
   want to change the file.  Outputs variable names such that they are equal to the
   WMO grib-2 naming standard (or at least tries to ...)
   which makes it easy to plot reanalysis grib2 data
   and the WRF output with the same plotting script.
   REQUIRES that supplmentary scripts: wrf_synoptic.py and thermodynamics.py are located
   in the same folder as this script to operate.

   !! WARNING !!
   Compared to other modules within this package, this is VERY Unforgiving.  It does not pre-checks for correct times
   and levels prior to running through the files.  This is mostly because of the way the ERA-Interim data is stored in
   the grib2 format.  Specifically, certain variables are stored with specific file prefixes, and each file only contains
   a single time step.  Because the data is spread out over numerous files doing a full check to ensure all time steps and
   levels match each other is burdensome and I don't want to do it.  The best way to avoid these issues is simply to
   download all the ERA-data you need at one time, for all variables and time steps, and place them in the same directory.
   Avoid downloading some data using one subset, and other data using a different subset. This will lead to problems.
   !!!!!!!!!!!!

   Coding format follows PEP8 (ish) rules, but not perfectly.

   Classes: EraPost -> Defines object with functions for post processing grib2 ERA data.  Requires "pygrib" library.
"""

import numpy as np
from ERDCMOOSE.functions import common as helpers
from netCDF4 import Dataset
from datetime import datetime, timedelta
import glob as glob
import pygrib as pg


class EraPost:
    def __init__(self):
        """ERA Post is a class object that defines the namelist and class objects"""

        self.FilePfxDict = {'uv': ['U velocity', 'V velocity'],
                            'sc': ['Geopotential', 'Temperature', 'Relative humidity',
                                   '10 metre U wind component','10 metre V wind component',
                                   '2 metre temperature','2 metre dewpoint temperature','Mean sea level pressure',
                                   'Total column water vapour','Total cloud cover','Total column water'],
                            'sc.10u': ['10 metre U wind component'],
                            'sc.10v': ['10 metre V wind component'],
                            'sc.2t': ['2 metre temperature'],
                            'sc.2d': ['2 metre dewpoint temperature'],
                            'sc.msl': ['Mean sea level pressure'],
                            'sc.tcwv': ['Total column water vapour'],
                            'sc.tcc': ['Total cloud cover'],
                            'sc.tcw_': ['Total column water']}

        self.atmos_description_dict = {"U velocity": 'Wind Component in the Zonal (x) Direction',
                                       "V velocity": 'Wind Component in the Meridonal (y) Direction',
                                       'Geopotential': 'Geopotential Height',
                                       'Temperature': 'Temperature',
                                       'Relative humidity': 'Relative Humidity'}

        self.sfc_description_dict={'10 metre U wind component':'10 meter zonal (x) wind component',
                      '10 metre V wind component':'10 meter meridonal (y) wind component',
                      '2 metre temperature':'2 AGL Temperature',
                      '2 metre dewpoint temperature':'2 AGL Dewpoint Temperature',
                      'Mean sea level pressure':'Mean sea level pressure',
                      'Total column water vapour':'Total column water vapour (precipitable water)',
                      'Total column water':'Total column water (Vertically Integrated Liquid)',
                      'Total cloud cover':'Total Cloud Cover'}

        self.namelist_dictionary = {'levels': [1000, 850, 700, 750, 600, 500, 400, 300, 200, 100, 50],
                                    'atmos_vars': ['U velocity', 'V velocity', 'Geopotential', 'Temperature',
                                                   'Relative humidity'],
                                    'sfc_vars': ['10 metre U wind component', '10 metre V wind component',
                                                 '2 metre temperature', '2 metre dewpoint temperature',
                                                 'Mean sea level pressure', 'Total column water vapour'],
                                    'append': False, 'era_directory': './', 'output_path': './',
                                    'out_pfx': './'}

        self.namelist_dictionary['time_format'] = '%Y-%m-%d_%H:%M:%S'
        self.namelist_dictionary['time_min'] = '2019-02-24_12:00:00'
        self.namelist_dictionary['time_max'] = '2019-02-25_00:00:00'

        self.namelist_dictionary['output_path'] = self.namelist_dictionary['out_pfx'] + 'era_post.nc'

    def print_namelist(self):
        """A simple function to print out all of the namelist values!"""

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

    def AtmosPost(self):
        from datetime import datetime, timedelta
        """This function does all the heavy lifting for the atmospheric variable reanalysis."""
        data_dict = {}

        create_data_dict = True
        levels = self.namelist_dictionary['levels']

        all_var_list = []
        for i in self.namelist_dictionary['atmos_vars']:
            for key in self.FilePfxDict:
                if i in self.FilePfxDict[key]:
                    all_var_list.append(i)
                    break

        time_start = datetime.strptime('00:00:00 01-01-1900', '%H:%M:%S %d-%m-%Y')
        outputfile = self.namelist_dictionary['output_path']

        if self.namelist_dictionary['append'] == False:
            outputdata = Dataset(outputfile, 'w')  ## open up netCDF file for write.
            outputdata.Description = """NETCDF output from ERA-Interim Reanalysis grib2
                 atmospheric data interpolated to pressure levels.  Generated from the
                 European Centre for Medium-Range Weather Forecasts. Data available from NCAR at:
                 https://rda.ucar.edu/#!lfd?nb=y&b=proj&v=ECMWF%20Interim%20Reanalysis"""
            outputdata.History = "File Created %s" % datetime.now().strftime('%c')

            ## some projection information ##
            ## Note that this is FIXED to the ERA-Interim grid.
            outputdata.projection = 'Gaussian Latitude/Longitude Grid'
            outputdata.latitudeOfFirstGridPointInDegrees = 89.463
            outputdata.latitudeOfLastGridPointInDegrees = -89.463
            outputdata.longitudeOfFirstGridPointInDegrees = 0.0
            outputdata.longitudeOfLastGridPointInDegrees = 359.297

            ## Create Dimensions ##
            outputdata.createDimension('z', len(levels))
            outputdata.createDimension('time', None)

            levs = outputdata.createVariable('Levels', 'f4', 'z')
            time = outputdata.createVariable('Time', 'f4', 'time')

            time.units = "Hours since %s" % time_start.strftime('%H:%M:%S %d-%m-%Y')
            time.description = "Hours since %s" % time_start.strftime('%H:%M:%S %d-%m-%Y')

            levs[:] = np.array(levels)
            levs.units = 'isobaricInhPa'
            levs.description = "Vertical Pressure Coodinate"

        else:
            outputdata = Dataset(outputfile, 'r+')  ## open up netCDF file for write.

        for vdx, v in enumerate(self.namelist_dictionary['atmos_vars']):
            for key in self.FilePfxDict:
                if v in self.FilePfxDict[key]:
                    print("WORKING ON %s" % v)
                    # Variable is in this dictionary key --> do your magic!
                    gribfiles = sorted(glob.glob(self.namelist_dictionary['era_directory'] + '*ei*.pl.*%s*' % key))
                    # Now loop through grib files ...
                    tints = []
                    times = []

                    save_indices = []

                    tmin = datetime.strptime(self.namelist_dictionary['time_min'],
                                             self.namelist_dictionary['time_format'])
                    tmax = datetime.strptime(self.namelist_dictionary['time_max'],
                                             self.namelist_dictionary['time_format'])
                    for gbx, gb in enumerate(gribfiles):
                        grb = pg.open(gb)
                        save_time = True
                        # NOW LOOP THROUGH LEVELS! #
                        for ldx, lev in enumerate(levels):
                            try:
                                cdata = grb.select(parameterName=v, level=lev)[0]
                            except:
                                self.print_warning("%s is not a valid level, skipping" % lev)
                                continue

                            if create_data_dict == True:
                                lats, lons = cdata.latlons()  ## assume lat lon constant.=
                                outputdata.createDimension('east_west', lons.shape[1])
                                if self.namelist_dictionary['append'] == False:
                                    outputdata.createDimension('north_south', lons.shape[0])
                                    longitude = outputdata.createVariable('Longitude', 'f4', ('north_south', 'east_west'))
                                    latitude = outputdata.createVariable('Latitude', 'f4', ('north_south', 'east_west'))
                                    longitude[:] = lons
                                    latitude[:] = lats

                                    latitude.units = "Degrees North"
                                    longitude.units = "Degrees East"

                                for vari in all_var_list:
                                    data_dict[vari] = [np.zeros([len(gribfiles), len(levels)] + list(np.shape(lats)))]

                                create_data_dict = False  ## turn off create!

                            if save_time == True:
                                tnow = datetime.strptime(str(cdata['dataDate']) + "{:04d}".format(cdata['dataTime']),
                                                         '%Y%m%d%H%M')

                                ## Check here to see if this time is within the time_range ##
                                ## should be a simple, if between conditional. ##
                                if tmin <= tnow <= tmax:
                                    tints.append((tnow - time_start).total_seconds() / 3600.)
                                    times.append(
                                        datetime.strptime(str(cdata['dataDate']) + "{:04d}".format(cdata['dataTime']),
                                                          '%Y%m%d%H%M'))
                                    save_time = False  # Turn off Save time.
                                    save_indices.append(gbx)
                                else:
                                    break

                            data_dict[v][0][gbx, ldx, :] = cdata.values
                            if len(data_dict[v]) == 1:
                                data_dict[v].append(cdata['parameterUnits'])
                                data_dict[v].append(self.atmos_description_dict[i])

                        grb.close()
                    data_dict[v][0] = data_dict[v][0][tuple(save_indices), :]
        for i in data_dict:
            ncoutvar = outputdata.createVariable(i, 'f4', ('time', 'z', 'north_south', 'east_west'))
            ncoutvar.units = data_dict[i][1]
            ncoutvar.longname = data_dict[i][2]
            if i == 'Geopotential':
                ncoutvar[:] = data_dict[i][0][:] / 9.81  ## divide by gravity to get meters!
            else:
                ncoutvar[:] = data_dict[i][0][:]  ## Otherwise keep data as is.
            time[:] = tints

        outputdata.close()

    def SfcPost(self):
        from datetime import datetime, timedelta
        """This function does all the heavy lifting for the Surface variable reanalysis."""
        data_dict = {}

        create_data_dict = True

        levels = self.namelist_dictionary['levels']

        times = []
        all_var_list = []
        for i in self.namelist_dictionary['sfc_vars']:
            for key in self.FilePfxDict:
                if i in self.FilePfxDict[key]:
                    all_var_list.append(i)
                    break

        time_start = datetime.strptime('00:00:00 01-01-1900', '%H:%M:%S %d-%m-%Y')
        outputfile = self.namelist_dictionary['output_path']

        if self.namelist_dictionary['append'] == False:
            outputdata = Dataset(outputfile, 'w')  ## open up netCDF file for write.
            outputdata.Description = """NETCDF output from ERA-Interim Reanalysis grib2
                 atmospheric data interpolated to pressure levels.  Generated from the
                 European Centre for Medium-Range Weather Forecasts. Data available from NCAR at:
                 https://rda.ucar.edu/#!lfd?nb=y&b=proj&v=ECMWF%20Interim%20Reanalysis"""
            outputdata.History = "File Created %s" % datetime.now().strftime('%c')

            ## some projection information ##
            ## Note that this is FIXED to the ERA-Interim grid.
            outputdata.projection = 'Gaussian Latitude/Longitude Grid'
            outputdata.latitudeOfFirstGridPointInDegrees = 89.463
            outputdata.latitudeOfLastGridPointInDegrees = -89.463
            outputdata.longitudeOfFirstGridPointInDegrees = 0.0
            outputdata.longitudeOfLastGridPointInDegrees = 359.297

            ## Create Dimensions ##
            outputdata.createDimension('z', len(levels))
            outputdata.createDimension('time', None)

            levs = outputdata.createVariable('Levels', 'f4', 'z')
            time = outputdata.createVariable('Time', 'f4', 'time')

            time.units = "Hours since %s" % time_start.strftime('%H:%M:%S %d-%m-%Y')
            time.description = "Hours since %s" % time_start.strftime('%H:%M:%S %d-%m-%Y')

            levs[:] = np.array(levels)
            levs.units = 'isobaricInhPa'
            levs.description = "Vertical Pressure Coordinate"

        else:
            outputdata = Dataset(outputfile, 'r+')  ## open up netCDF file for write.

        save_time = True
        save_indices=[]
        times=[]
        tints=[]
        for vdx, v in enumerate(self.namelist_dictionary['sfc_vars']):
            for key in self.FilePfxDict:
                if v in self.FilePfxDict[key]:
                    # Variable is in this dictionary key --> do your magic!
                    gribfiles = glob.glob(self.namelist_dictionary['era_directory'] + '*ei*.sfc.*%s*' % key)
                    # Now loop through grib files ...

                    tmin = datetime.strptime(self.namelist_dictionary['time_min'],
                                             self.namelist_dictionary['time_format'])
                    tmax = datetime.strptime(self.namelist_dictionary['time_max'],
                                             self.namelist_dictionary['time_format'])
                    for gbx, gb in enumerate(gribfiles):
                        grb = pg.open(gb)
                        cdata = grb.select(parameterName=v)
                        ctimes = []
                        if self.namelist_dictionary['append'] == True: ## appending to existing file!
                            available_times = [datetime.strptime(outputdata.variables['Time'].units, \
                                                "Hours since %s"%time_start.strftime('%H:%M:%S %d-%m-%Y'))+ \
                                                timedelta(hours=float(ss)) for ss in outputdata.variables['Time'][:]]

                            ndata = cdata[0]
                            ctime = datetime.strptime(
                                str(ndata['dataDate']) + "{:04d}".format(ndata['dataTime']), '%Y%m%d%H%M')
                            if ctime in available_times:
                                tdx=available_times.index(ctime)
                                ## now loop through matching time index list.
                                if create_data_dict == True:
                                    lats, lons = ndata.latlons()  ## assume lat lon constant.
                                    for i in all_var_list:
                                        data_dict[i] = [np.zeros([len(available_times)] + list(np.shape(lats)))]
                                    create_data_dict = False

                                data_dict[v][0][tdx, :] = ndata.values
                                if len(data_dict[v]) == 1:
                                    data_dict[v].append(ndata['parameterUnits'])
                                    data_dict[v].append(self.sfc_description_dict[v])
                        else:
                            if create_data_dict == True:
                                lats, lons = cdata.latlons()  ## assume lat lon constant.=
                                outputdata.createDimension('west_east', lons.shape[1])
                                outputdata.createDimension('south_north', lons.shape[0])
                                longitude = outputdata.createVariable('Longitude', 'f4', ('south_north', 'west_east'))
                                latitude = outputdata.createVariable('Latitude', 'f4', ('south_north', 'west_east'))
                                longitude[:] = lons
                                latitude[:] = lats

                                latitude.units = "Degrees North"
                                longitude.units = "Degrees East"

                                for vari in all_var_list:
                                    data_dict[vari] = [np.zeros([len(gribfiles)] + list(np.shape(lats)))]

                                create_data_dict = False  ## turn off create!

                            for tdx in range(len(cdata)):
                                if save_time == True:
                                    tnow = datetime.strptime(str(cdata['dataDate']) + "{:04d}".format(cdata['dataTime']),
                                                             '%Y%m%d%H%M')

                                    ## Check here to see if this time is within the time_range ##
                                    ## should be a simple, if between conditional. ##
                                    if tmin <= tnow <= tmax:
                                        tints.append((tnow - time_start).total_seconds() / 3600.)
                                        times.append(
                                            datetime.strptime(str(cdata['dataDate']) + "{:04d}".format(cdata['dataTime']),
                                                              '%Y%m%d%H%M'))
                                        save_time = False  # Turn off Save time.
                                        save_indices.append(gbx)
                                    else:
                                        continue

                            save_time = False
                            data_dict[v][0][gbx, :] = cdata.values
                            if len(data_dict[v]) == 1:
                                data_dict[v].append(cdata['parameterUnits'])
                                data_dict[v].append(self.atmos_description_dict[i])

                            grb.close()
                            data_dict[v][0] = data_dict[v][0][tuple(save_indices), :]

        for i in data_dict:
            if i in outputdata.variables:
                # VARIABLE ALREADY EXISTS IN NETCDF FILE --> Overwrite it!###
                print("Variable already exists!  Overwriting ...")
                ncoutvar = outputdata.variables[i]
            else:
                ncoutvar = outputdata.createVariable(i, 'f4', ('time', 'north_south', 'east_west'))
            ncoutvar.units = data_dict[i][1]
            ncoutvar.longname = data_dict[i][2]
            ncoutvar[:] = data_dict[i][0][:]
        if self.namelist_dictionary['append']==False:
            time[:]=tints

        outputdata.close()

    def print_warning(self, message):
        print("-----------------")
        print("|    WARNING    |")
        if type(message) is not list:
            print(message)
        else:
            for m in message:
                print(m)
        print("-----------------")

    def ParseNamelist(self, namelist_file='era_post_namelist.txt'):
        """This is the function that parses the namelist and puts it into a neat dictionary
            that is used in the script.  Some pretty stupid rules are:
            NO QUOTES in the namelist!  No "" or '', I know, it seems counter intuitive to
            send strings with spaces without quotation marks, but if you include quotes anywhere
            in an uncommented section of the namelist, you're going to have a bad time.
            Speaking of comments: us a # to indicate a comment.  Function will ignore
            anything written on a line following a #.  Comma separate things, you don't need to include a
            comma at the end of a line, but it also doesn't hurt.  Can write anything into the namelist you want,
            but if it's not in the list: "lookforkeys" it won't be read into the dictionary.
            the 'list_keys' list indicates that the function is looking for multiple values,
            if it doesn't find multiple values, it will just return a list with length 1."""

        lookforkeys = ['levels', 'atmos_vars', 'sfc_vars', 'time_max', 'append',
                       'era_directory', 'out_pfx', 'time_format', 'output_path', 'time_min']

        list_keys = ['levels', 'atmos_vars', 'sfc_vars']

        try:
            namelistfile = glob.glob(namelist_file)[0]
        except RuntimeError:
            print("NAMELIST FILE: %s does not exist, check your paths are try again!" % namelist_file)

        # OKAY, DEFAULTS SET! ##
        with open(namelistfile) as fp:
            for line in fp:
                var_key = line.split('=')[0].strip(' \t\n\r')
                if var_key in lookforkeys and '#' not in var_key:  # if not a comment!
                    values = line.split('=')[1].split('#')[0].strip(' \t\n\r').split(',')
                    values = [x for x in values if x != '']
                    if len(values) > 1:
                        if var_key in list_keys:
                            self.namelist_dictionary[var_key] = values
                        else:
                            print("WARNING, THIS IS NOT A LIST VARIABLE")
                            print("SETTING KEY VALUE TO FIRST OPTION")
                            self.namelist_dictionary[var_key] = values[0]
                    else:
                        if var_key in list_keys:
                            self.namelist_dictionary[var_key] = values
                        else:
                            self.namelist_dictionary[var_key] = values[0]

                    print(var_key, self.namelist_dictionary[var_key])

        return self.namelist_dictionary
