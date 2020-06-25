# coding=utf-8
"""This script calls a "namelist" and post processes the raw wrfout NETCDF
   files into a single netcdf file with user specified variables.
   Note: This is not a replacement for the much faster and more efficient
   nco software for basic time and variable slicing.
   However, this script is useful for interpolating to common pressure levels,
   computing and saving derived variables (e.g., mslp or precipitable water),
   and it also does time and variable slicing.  It can, however, take a long time
   to run depending on how many grid points you're processing (e.g., I've seen it take
   upwards of 7 minutes PER vertically varying variable).  To mitigate this
   I added "appending" options make it easy to add a variable to an existing
   file so that you don't need to recompute the full 3D atmosphere everytime you
   want to change the file.  Outputs variable names such that they are equal to the
   WMO grib-2 naming standard (or at least tries to ...)
   which makes it easy to plot reanalysis grib2 data
   and the WRF output with the same plotting script.
   REQUIRES that supplmentary scripts: wrf_synoptic.py and thermodynamics.py are located
   in the same folder as this script to operate.
   Coding format follows PEP8 (ish) rules, but not perfectly.
"""

import numpy as np
from MOOSEfunctions import common as helpers
from netCDF4 import Dataset
from datetime import datetime, timedelta
import glob as glob

class WrfPost:
    """The MOOSE V 1.0 WrfPost class constructor

       WRF Post performs post processing of raw WRF output from wrfout_d0X_YYYY-mm-dd_HH:MM:SS files and
       outputs post processed and time-subset data to a NetCDF file with variable names matching WMO variable
       naming convention.

       Can process a subset of WRF-Chem air-quality variables.

       Usage:
            - wp = WrfPost()
            - wp.parse_namelist()
            - wp.AtmosPost()

        To print available output variables you can use:

            - wp.PrintPostVars()
        """
    def __init__(self):
        # First Define the "namelist_dictionary" with default values.

        self.namelist_dictionary = {'levels': [1000, 850, 700, 750, 600, 500, 400, 300, 200, 100, 50],
                                    'atmos_vars': ['U velocity', 'V velocity', 'Geopotential', 'Temperature',
                                                   'Relative humidity', 'Visibility'],
                                    'sfc_vars': ['10 metre U wind component', '10 metre V wind component',
                                                 '2 metre temperature', '2 metre dewpoint temperature',
                                                 'Mean sea level pressure', 'Total column water vapour'],
                                    'static_vars': ['EROD', 'Terrain', 'Land Cover'],
                                    'chem_vars': ['EXTCOF55', 'TOT_DUST'], 'append': 0, 'wrf_directory': './',
                                    'wrf_pfx': 'test'}

        self.namelist_dictionary['output_path'] = './'
        self.namelist_dictionary['time_format'] = '%Y-%m-%d_%H:%M:%S'
        self.namelist_dictionary['time_min'] = '2019-02-24_12:00:00'
        self.namelist_dictionary['time_max'] = '2019-02-25_00:00:00'
        self.namelist_dictionary['scale_uv'] = 0
        self.namelist_dictionary['nest'] = 'd01'
        self.namelist_dictionary['dust_scheme'] = 1

        self.namelist_dictionary['scale_uv'] = bool(self.namelist_dictionary['scale_uv'])
        self.namelist_dictionary['append'] = bool(self.namelist_dictionary['append'])

        # Initialize WRFTimes and TimeIndexDict variables.
        self.WRFTimes = None  # Will be set to list type with time integers
        self.TimeIndexDict = None  # Will be dictionary type with filename keys, and index values
        self.TimeStart = datetime.strptime('00:00:00 01-01-1900', '%H:%M:%S %d-%m-%Y')


    def __repr__(self):
        outfile = self.namelist_dictionary['output_path'] + 'wrfoutput_post_%s.nc' \
                % (self.namelist_dictionary['wrf_pfx'])
        indir=self.namelist_dictionary['wrf_directory']
        return (f'{self.__class__.__name__}('
                f'WRF data directory {indir!r}, output file={outfile!r})')


    def print_namelist(self):
        """Help on built-in function print_namelist in module WrfPost:

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

    def standard_sfc_var(self, vari, data, df):
        """Help on built-in function print_namelist in module WrfPost:

        standard_sfc_var(...)
            standard_sfc_var(variable,data,df)

            Simple helper function to cut down on repetitive code when loading a standard WRF variable.
            Function manages the possibility that wrfout files were saved with either a single, or multiple time steps.

            Parameters
            ----------
            vair : string, name of variable in netCDF file
            data : NetCDF object created using the netCDF4 python package, input wrfout data
            df   : string, path for netCDF file, used to identify the time-index to build the time-slice


            Returns
            -------
            nothing.
        """
        if np.shape(data.variables[vari][:])[0] > 1:
            out = data.variables[vari][slice(*self.TimeIndexDict[df]), :]
        else:
            out = data.variables[vari][:]

        return out


    def PrintPostVars(self):
        """Help on built-in function PrintPostVars in module WrfPost:

        PrintPostVars(...)
            PrintPostVars()

            User-callable function that will print out the full list of variables and descriptions that can be computed or saved
            during each of the WrfPost functions:
                - StaticPost
                - AtmosPost
                - ChemPostAtmos
                - SfcPost
                - ChemPostSfc

            Returns
            -------
            nothing.
        """
        static_description_dict = {'3D Erodibility': "Surface Soil Erodibility Scaling Factor (sand/clay/silt)",
                                   '2D Erodibility': "Surface Soil Erodibility Scaling Factor Specific to Dust Scheme",
                                   'Analog Dust': "Analog Dust Potential from ERDC-Geo",
                                   'Terrain': 'Surface Elevation',
                                   'Dust Potential': 'Dust Lofting Potential',
                                   'Land Cover': 'Surface Land Category',
                                   'Sand Fraction': 'Fraction of surface soil that is sand',
                                   'Clay Fraction': 'Fraction of surface soil that is clay',
                                   'Soil Category': 'Dominant top-layer soil category',
                                   'Vegetation Fraction': 'Percent of grid cell that is considered vegetated'}

        atmos_description_dict = {"U velocity": 'Wind Component in the Zonal (x) Direction',
                                  "V velocity": 'Wind Component in the Meridonal (y) Direction',
                                  "W velocity": 'Wind Component in the Vertical (z) Direction',
                                  'Geopotential': 'Geopotential Height',
                                  'Temperature': 'Temperature',
                                  'Potential Temperature': 'Potential Temperature',
                                  'Eq Potential Temperature' : 'Equivalent Potential Temperature (Stull, 1988 Approximation)',
                                  'Relative humidity': 'Relative Humidity',
                                  'Cloud water mixing ratio': 'Cloud water mass mixing ratio',
                                  'Rain mixing ratio': 'Rain water mass mixing ratio',
                                  'Cloud Ice mixing ratio': 'Cloud Ice mass mixing ratio',
                                  'Snow mixing ratio': 'Snow mass mixing ratio',
                                  'Cloud Fraction': 'Sub grid cloud fraction',
                                  'Visibility': "Estimated visibility",
                                  'Radar Reflectivity': "Radar Reflectivity"
                                  }

        sfc_description_dict = {'10 metre U wind component': '10 meter zonal (x) wind component',
                                '10 metre V wind component': '10 meter meridonal (y) wind component',
                                '2 metre temperature': '2 AGL Temperature',
                                '2 metre dewpoint temperature': '2 AGL Dewpoint Temperature',
                                'Mean sea level pressure': 'Mean sea level pressure',
                                'Surface pressure' : 'Surface Pressure',
                                'Total column water vapour': 'Total column water vapour (precipitable water)',
                                'Gridscale Precipitation rate': 'Grid-scale precipitation rate',
                                'Convective Precipitation rate': 'Precipitation rate from the convective parameterization',
                                'Planetary boundary layer height': "Height of the boundary layer top",
                                'Visibility (Surface)': "Horizontal Visibility",
                                'Snow cover fraction': 'Fractional Snow Cover within Gridcell',
                                'Gridscale Snowfall rate': 'Grid-scale Snowfall Rate',
                                'Friction Velocity': 'Surface Friction Velocity',
                                'z0': 'Surface Roughness Length',
                                'Outgoing Longwave Radiation': 'Upwelling Longwave Radiation At the Top of the Atmosphere',
                                '1km AGL Reflectivity': "Simulated Radar Reflectivity at 1km Above Ground Level",
                                'Soil Moisture': 'Soil Moisture'}

        chem_description_dict = {'Black Carbon': "Black Carbon Concentration 1",
                                 'SO2': "Sulfur Dioxide Concentration",
                                 'Organic Carbon': "Organic Carbon Concentration 1",
                                 'Total Dust': "Total Dust Concentration",
                                 'EXTCOF55': "Extinction Coefficient at 550 nm",
                                 'PM 2.5': "Aerosol Particulate Matter ( d <2.5 microns) Dry-mass Concentration",
                                 'PM 10': "Aerosol Particulate Matter ( d <10 microns) Dry-mass Concentration",
                                 'DMS': 'Dimethyl sulfide ([CH_3]_2S) Concentration',
                                 'Sea Salt': 'Total Sea Salt Concentration'}

        chemsfc_description_dict = {'Black Carbon (Surface)': "Black Carbon Concentration 1",
                                 'SO2 (Surface)': "Sulfur Dioxide Concentration",
                                 'Organic Carbon (Surface)': "Organic Carbon Concentration 1",
                                 'Total Dust (Surface)': "Total Dust Concentration at the surface",
                                 'EXTCOF55 (Surface)': "Extinction Coefficient at 550 nm",
                                 'Aerosol Optical Depth': 'Aerosol Optical Depth at 550nm',
                                 'Ideal Dust': 'Idealized Dust Flux',
                                 'Dust Scaling Factor': 'Dust Scaling Factor Computed from Static '
                                                        'Analog Dust Map and Idealize Dust Flux',
                                 'Total Dust Emission': "Total Dust Emission",
                                 'Aerosol Visibility (Surface)':
                                        "Surface Visibility (550 nm) Considering Aerosol Particles from WRF-Chem Only",
                                 'PM 2.5 (Surface)': "Aerosol Particulate Matter ( d <2.5 microns) Dry-mass Concentration",
                                 'PM 10 (Surface)': "Aerosol Particulate Matter ( d <10 microns) Dry-mass Concentration",
                                 'DMS (Surface)': 'Dimethyl sulfide ([CH_3]_2S) Concentration',
                                 'Sea Salt (Surface)': 'Total Sea Salt Concentration'}

        print("---------------------------")
        print("Printing Out Allowable Variables")
        print("")
        print("StaticPost")
        print("-----------")
        for s in static_description_dict:
            print('  > %s - %s'%(s, static_description_dict[s]))
        print("")
        print("AtmosPost")
        print("-----------")
        for s in atmos_description_dict:
            print('  > %s - %s'%(s, atmos_description_dict[s]))

        print("SfcPost")
        print("-----------")
        for s in sfc_description_dict:
            print('  > %s - %s'%(s, sfc_description_dict[s]))

        print("ChemPostAtmos")
        print("-----------")
        for s in chem_description_dict:
            print('  > %s - %s'%(s, chem_description_dict[s]))

        print("ChemPostSfc")
        print("-----------")
        for s in chemsfc_description_dict:
            print('  > %s - %s'%(s, chemsfc_description_dict[s]))

    def SubsetTime(self):

        """Help on built-in function print_namelist in module SubsetTime:

                SubsetTime(...)
                    SubsetTime()

                    This function is designed to take the user input time boundaries defined in the
                    namelist dictionary: time-max,time_min. It replaces the None Type WRFTimes variable with a list of datetime
                    objects and indices from the wrf-files time variable.  Called automatically as part of all post processing functions

                    Returns
                    -------
                    nothing, but sets the WrfPost attributes "WRFTimes", "TimeIndexDict"
                """

        # Initialize start time
        # Then get all WRF datafiles.
        datafiles = sorted(glob.glob(self.namelist_dictionary['wrf_directory']
                                     + 'wrfout_%s*' % self.namelist_dictionary['nest']))

        time_bounds = self.namelist_dictionary['time_min'], self.namelist_dictionary['time_max']
        time_bounds_strp = sorted([datetime.strptime(i, '%Y-%m-%d_%H:%M:%S') for i in time_bounds])

        # Because, there may be several files that fall outside of the max/min times, need to define a new list:
        #   "final_files" that reflects only the files that contain data within the time boundaries.
        final_files = []
        # Now, we need to open each file individually to get the time information from them.
        print("Getting Time Information from %i files")
        time_inds = {}  # Define dictionary to put indices in -> dictionary keys are file names, and values are
        # start/ending indices
        wrftimeints = []  # create an empty list in case there are no files!
        for fdx, f in enumerate(datafiles):
            cfile = Dataset(f, 'r')

            # Get times from within WRF file.
            wrftimes = [datetime.strptime(b"".join(tm).decode("utf-8"), '%Y-%m-%d_%H:%M:%S') for tm in
                        cfile.variables['Times']]

            # If no times in the file are useful.
            if time_bounds_strp[1] < np.min(wrftimes) or time_bounds_strp[0] > np.max(wrftimes):
                print("Entire file outside of time bounds, skipping file.")
                continue

            if time_bounds_strp[0] < np.min(wrftimes) and time_bounds_strp[1] > np.max(wrftimes):
                # Case:, end time bound is past the end time index, start time is less than start index
                # KEEP ENTIRE ARRAY!
                time_inds[f] = [0, len(wrftimes)]
            elif time_bounds_strp[0] > np.min(wrftimes) and time_bounds_strp[1] > np.max(wrftimes):
                # Case: start bounds is outside of dataset end bound is not, and end time is outside
                time_inds[f] = [np.argmin(np.abs(np.array(wrftimes) - time_bounds_strp[0])), len(wrftimes)]
            elif time_bounds_strp[0] < np.min(wrftimes) and time_bounds_strp[1] < np.max(wrftimes):
                # Case: start bounds is outside of dataset end bound is within it.
                time_inds[f] = [0, np.argmin(np.abs(np.array(wrftimes) - time_bounds_strp[1]))]
            else:
                # ONLY remaining option, both are within! ##
                time_inds[f] = [np.argmin(np.abs(np.array(wrftimes) - time_bounds_strp[0])),
                                np.argmin(np.abs(np.array(wrftimes) - time_bounds_strp[1]))]

            final_files.append(f)
            # Save WRF-Times as an integer index rather than a datetime object, this variable will be saved as
            # a list of integers
            wrftimes = [(tt - self.TimeStart).total_seconds() / 3600. for tt in wrftimes[slice(*time_inds[f])]]

            # simple try statement that decides whether to intialize the array
            if fdx == 0:
                print("Initializing WRF time list...")
                wrftimeints = wrftimes
            else:
                wrftimeints = wrftimeints + wrftimes

            cfile.close()

        self.WRFTimes = wrftimeints
        self.TimeIndexDict = time_inds

    def StaticPost(self, wrfinput=True):
        """Help on built-in function print_namelist in module StaticPost:

        StaticPost(...)
            StaticPost([wrfinput])

            This function creates the static (2D) data arrays that get output to the post processed output file.
            Function parses the namelist_dictionary attribute, reads in the static arrays from the WRF data
            and saves the data to the post processed netCDF file. If append = 0, will created a new file, otherwise will
            save to an existing post processed file.  Note that if append = 1, all data will be overwritten.

            Parameters
            ----------
            wrfinput : Boolean (optional), Optional flag, if True, all static data will come for a "wrfinput" instead
                of a "wrfout" file.  This is important if the user is trying to build static data from variables
                that are only saved to the wrfinput intitialization file during real.exe (e.g., EROD).  Note, that this
                option has a minor impact on the metadata saved to the post processed file, as there are minor differences
                in metadata between "wrfinput" and "wrfout" files.

            Returns
            ----------
            Nothing.

        """

        static_description_dict = {'3D Erodibility': "Surface Soil Erodibility Scaling Factor (sand/clay/silt)",
                                   '2D Erodibility': "Surface Soil Erodibility Scaling Factor Specific to ",
                                   'Analog Dust': "Analog Dust Potential from ERDC-Geo",
                                   'Terrain': 'Surface Elevation',
                                   'Dust Potential' : 'Dust Lofting Potential',
                                   'Land Cover': 'Surface Land Category',
                                   'Sand Fraction': 'Fraction of surface soil that is sand',
                                   'Clay Fraction': 'Fraction of surface soil that is clay',
                                   'Soil Category': 'Dominant top-layer soil category',
                                   'Vegetation Fraction': 'Percent of grid cell that is considered vegetated'}

        static_unit_dict = {'3D Erodibility': '-','2D Erodibility': '-', 'Terrain': 'm', 'Land Cover': 'integer category',
                            'Sand Fraction': 'fractional', 'Clay Fraction': 'fractional',
                            'Dust Potential':'mg m^-2 s^-1','Vegetation Fraction': '%',
                            'Soil Category': 'integer category','Analog Dust': 'mg m^-2 s^-1'}

        if wrfinput:
            fpfx_str = 'wrfinput'
        else:
            fpfx_str = 'wrfout'
        # Check if file exists.
        try:
            datafile = sorted(glob.glob(self.namelist_dictionary['wrf_directory']
                                        + '%s_%s*' % (fpfx_str, self.namelist_dictionary['nest'])))[0]
            print("Opening %s to gather metadata," % datafile)
        except RuntimeError:
            print("NO %s nest wrfout files in folder %s found.... exiting." %
                  (self.namelist_dictionary['nest'], self.namelist_dictionary['wrf_directory']))

        metadata = Dataset(datafile, 'r')
        meta_data_dict = GetMetaData(metadata,wrfinput=wrfinput)
        # GET LAT / LON AND ADD TO FILE! ##
        lats = metadata.variables['XLAT'][0, :].squeeze()
        lons = metadata.variables['XLONG'][0, :].squeeze()
        metadata.close()

        # define output file from the namelist dictionary values.
        # noinspection PyPep8
        outputfile = self.namelist_dictionary['output_path'] + 'wrfoutput_post_%s.nc' \
                     % (self.namelist_dictionary['wrf_pfx'])

        levels = [float(l) for l in self.namelist_dictionary['levels']]
        if not self.namelist_dictionary['append']:  # If not appending.

            # Check if times already exists, if not, get (with) the times!
            if self.WRFTimes is None:
                print("WRF Time values are not defined... Defining them.")
                self.SubsetTime()

            outputdata = Dataset(outputfile, 'w')  # open up netCDF file for write.

            outputdata.Description = """NetCDF output from wrfoutput
                        atmospheric data interpolated to pressure levels."""

            outputdata.History = "File Created %s" % datetime.now().strftime('%c')

            outputdata.setncatts(meta_data_dict)

            # Create Dimensions #
            outputdata.createDimension('time', None)
            outputdata.createDimension('z', len(levels))
            time = outputdata.createVariable('Time', 'f4', 'time')

            time.units = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')
            time.description = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')

            time[:] = self.WRFTimes

            levs = outputdata.createVariable('Levels', 'f4', 'z')
            levs.units='hPa'
            levs.description='Model Pressure Levels'
            levs[:] = levels

            outputdata.createDimension('east_west', lons.shape[1])
            outputdata.createDimension('south_north', lons.shape[0])
            longitude = outputdata.createVariable('Longitude', 'f4', ('south_north', 'east_west'))
            latitude = outputdata.createVariable('Latitude', 'f4', ('south_north', 'east_west'))
            longitude[:] = lons
            latitude[:] = lats

            latitude.units = "Degrees North"
            longitude.units = "Degrees East"
        else:
            print("Appending to %s" % outputfile)
            outputdata = Dataset(outputfile, 'r+')  # open up netCDF file for appending

        print("---------------------------")
        print("netCDF output file %s Opened" % outputfile)
        print("Doing Static Data...")
        print(" -- Doing variables: ")
        for vv in self.namelist_dictionary['static_vars']:
            print("    -- %s" % vv)

        # SIMPLE ONLY OPEN FIRST FILE AND GET THAT DATA ... ###

        ncdata = Dataset(datafile, 'r')  # open for reading

        for vv in self.namelist_dictionary['static_vars']:
            if vv == '3D Erodibility':
                if 'EROD' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable EROD is not in this wrfout file, skipping!")

                elif '3D Erodibility' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting values ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['EROD'][:]
                else:
                    # START FROM SCRATCH
                    # ULTRA SPECIAL CASE, NEED TO MAKE NEW DIMENSION! ##
                    outputdata.createDimension('dust_erosion_category', 3)
                    ncoutvar = outputdata.createVariable(vv, 'f4', ('dust_erosion_category',
                                                                    'south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['EROD'][:]


            if vv == '2D Erodibility':
                if 'EROD' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable EROD is not in this wrfout file, skipping!")
                else:
                    ### NOTE THIS NEEDS A SPECIAL FUNCTION HERE TO ACCOUNT FOR THE FACT THAT ###
                    ### EDUST HAS DIFFERENT UNITS FOR ALL DUST SCHEMES.
                    ### THIS IS CONTROLLED BY THE USER IN THE NAMELIST WITH THE FLAG: dust_scheme
                    ### the dust_scheme flag corresponds to the wrf options
                    ### -> 1 = GOCART
                    ### -> 3 = AFWA
                    ### -> 4 = UoC
                    ### If TOT_EDUST is in wrfout file, then it MUST BE AFWA, will adjust accordingly
                    ### but return a warning.
                    ### Otherwise, it will use value in the namelist.

                    ### FIRST CHECK IF DUST SCHEME IS EVEN AN ALLOWABLE NUMBER, if it is not, throw warning and set to 1
                    if int(self.namelist_dictionary['dust_scheme']) not in (1, 3, 4):
                        print("!!!WARNING!!!")
                        print("dust_scheme %i is not allowed, the following schemes are allowed:" %
                              int(self.namelist_dictionary['dust_scheme']))
                        print(" dust_scheme = 1 -> GOCART (First Erodibility Category)")
                        print(" dust_scheme = 3 -> AFWA-GOCART (Sum of Erodibility Categories)")
                        print(" dust_scheme = 4 -> UoC (1 Where Sum of Erodibility Categories > 0.0001)")
                        print("Defaulting to GOCART (dust_scheme = 1")
                        self.namelist_dictionary['dust_scheme'] = 1
                        print("")
                        estring='GOCART'

                    if int(self.namelist_dictionary['dust_scheme']) == 1:
                        print(" dust_scheme = 1 -> Mapping 2D Erodibility assuming GOCART Dust Scheme")
                        estring='GOCART'
                        EROD=ncdata.variables['EROD'][:].squeeze()[0,:]
                    elif int(self.namelist_dictionary['dust_scheme']) == 3:
                        print(" dust_scheme = 3 -> Mapping 2D Erodibility assuming AFWA-GOCART Dust Scheme")
                        estring = 'AFWA-GOCART'
                        EROD=np.sum(ncdata['EROD'][:].squeeze(),axis=0)
                    else:
                        print(" dust_scheme = 4 -> Mapping 2D Erodibility assuming UoC Dust Scheme")
                        estring = 'UoC'
                        mask=np.ma.masked_less(np.sum(ncdata['EROD'][:].squeeze(),axis=0),0.0001).mask
                        EROD=np.ma.masked_array(np.ones_like(ncdata['EROD'][:].squeeze()[0,:]),mask=mask).filled(0.0)

                    if '2D Erodibility' in outputdata.variables:
                        # HEY, it's already here!
                        print("Variable already exists!  Overwriting values ...")
                        ncoutvar.units = static_unit_dict[vv]
                        ncoutvar.longname = static_description_dict[vv]
                        ncoutvar[:] = EROD
                    else:
                        # START FROM SCRATCH
                        ncoutvar = outputdata.createVariable(vv, 'f4', ('south_north', 'east_west'))
                        ncoutvar.units = static_unit_dict[vv]
                        ncoutvar.longname = static_description_dict[vv]+'%s Dust Scheme'%estring
                        ncoutvar[:] = EROD

            if vv == 'Soil Category':
                if 'SOILCTOP' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable SOILCTOP is not in this wrfinput or wrfouput file, skipping!")

                elif 'Soil Category' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting values ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    soil_cat=ncdata.variables['SOILCTOP'][:].squeeze()
                    ncoutarray=np.zeros((soil_cat.shape[1],soil_cat.shape[2]))

                    ## NEED to Get index of max value to get "dominant category"
                    for idx in range(soil_cat.shape[1]):
                        for jdx in range(soil_cat.shape[2]):
                            ncoutarray[idx, jdx] = np.argmax(soil_cat[:, idx, jdx]) + 1.  ## add one to index.
                    ncoutvar[:]=ncoutarray
                else:
                    # START FROM SCRATCH
                    ncoutvar = outputdata.createVariable(vv, 'f4', ('south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    soil_cat = ncdata.variables['SOILCTOP'][:].squeeze()
                    ncoutarray = np.zeros((soil_cat.shape[1], soil_cat.shape[2]))
                    ## NEED to Get index of max value to get "dominant category"
                    for idx in range(soil_cat.shape[1]):
                        for jdx in range(soil_cat.shape[2]):
                            ncoutarray[idx, jdx] = np.argmax(soil_cat[:, idx, jdx]) + 1.  ## add one to index.
                    ncoutvar[:] = ncoutarray

            if vv == 'Sand Fraction':
                if 'SANDFRAC' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable SANDFRAC is not in this wrfinput or wrfouput file, skipping!")

                elif 'Sand Fraction' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting values ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['SANDFRAC'][:]
                else:
                    # START FROM SCRATCH
                    ncoutvar = outputdata.createVariable(vv, 'f4', ('south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['SANDFRAC'][:]


            if vv == 'Clay Fraction':
                if 'CLAYFRAC' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable CLAYFRAC is not in this wrfinput or wrfouput file, skipping!")

                elif 'Clay Fraction' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting values ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['CLAYFRAC'][:]
                else:
                    # START FROM SCRATCH
                    ncoutvar = outputdata.createVariable(vv, 'f4', ('south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['CLAYFRAC'][:]

            if vv == 'Analog Dust':
                if 'DUST_POTENT' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable DUST_POTENT is not in this wrfinput or wrfouput file, skipping!")

                elif 'Analog' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting values ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['DUST_POTENT'][:]
                else:
                    # START FROM SCRATCH
                    ncoutvar = outputdata.createVariable(vv, 'f4', ('south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['DUST_POTENT'][:]

            if vv == 'Terrain':
                if 'HGT' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable HGT is not in this wrfout file, skipping!")
                elif 'Terrain' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['HGT'][0, :]
                else:
                    # START FROM SCRATCH
                    ncoutvar = outputdata.createVariable(vv, 'f4', ('south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['HGT'][0, :]

            if vv == 'Vegetation Fraction':
                print("Please Note: Vegetation Fraction is time-variable."
                      "However, it varies slowly and for most simulations, can be considered static.")
                if 'VEGFRA' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable VEGFRA is not in this wrfout file, skipping!")
                elif 'Vegetation Fraction' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['VEGFRA'][0, :]
                else:
                    # START FROM SCRATCH
                    ncoutvar = outputdata.createVariable(vv, 'f4', ('south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['VEGFRA'][0, :]

            if vv == 'Dust Potential':
                if 'DUST_POTENT' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable DUST_POTENT is not in this wrfout file, skipping!")
                elif 'Dust Potential' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['DUST_POTENT'][0, :]
                else:
                    # START FROM SCRATCH
                    ncoutvar = outputdata.createVariable(vv, 'f4', ('south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['DUST_POTENT'][0, :]

            if vv == 'Land Cover':
                if 'LU_INDEX' not in ncdata.variables:
                    print("!!!WARNING!!!")
                    print("Variable LU_INDEX is not in this wrfout file, skipping!")
                elif 'Land Cover' in outputdata.variables:
                    # HEY, it's already here!
                    print("Variable already exists!  Overwriting ...")
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['LU_INDEX'][0, :]
                else:
                    # START FROM SCRATCH
                    ncoutvar = outputdata.createVariable(vv, 'i2', ('south_north', 'east_west'))
                    ncoutvar.units = static_unit_dict[vv]
                    ncoutvar.longname = static_description_dict[vv]
                    ncoutvar[:] = ncdata.variables['LU_INDEX'][0, :]


        ncdata.close()
        outputdata.close()
        print("- Finished with static data ... ")

    def ChemPostAtmos(self):
        """Help on built-in function print_namelist in module ChemPostAtmos:

        ChemPostAtmos(...)
            ChemPostAtmos()

            This function creates the static (4D) data arrays that get output to the post processed output file.
            Function parses the namelist_dictionary attribute, reads in chemistry variables from the WRF-Chem data
            determined by the namelist_dictionary key: "chem_vars"
            and saves the data to the post processed netCDF file. If append = 0, will created a new file, otherwise will
            save to an existing post processed file.  Note that if append = 1, all data will be overwritten.

            Returns
            ----------
            Nothing.

        """

        chem_description_dict = {'Black Carbon': "Black Carbon Concentration 1",
                                 'SO2': "Sulfur Dioxide Concentration",
                                 'Organic Carbon': "Organic Carbon Concentration 1",
                                 'Total Dust': "Total Dust Concentration",
                                 'EXTCOF55': "Extinction Coefficient at 550 nm",
                                 'PM 2.5': "Aerosol Particulate Matter ( d <2.5 microns) Dry-mass Concentration",
                                 'PM 10': "Aerosol Particulate Matter ( d <10 microns) Dry-mass Concentration",
                                 'DMS': 'Dimethyl sulfide ([CH_3]_2S) Concentration',
                                 'Sea Salt': 'Total Sea Salt Concentration'}

        chem_unit_dict = {'Black Carbon': 'ug kg^-1 dry air', 'EXTCOF55': 'km^-1',
                          'Orgain Carbon': 'ug kg^-1 dry air', 'Total Dust': 'ug m^-3',
                          'SO2': 'ppmv','PM 10':'ug m^-3','PM 2.5':'ug m^-3','DMS':'ppmv','Sea Salt':'ug kg^-1 dry air'}

        data_dict = {}  # Create one large data dictionary that is used to store data temporarily before outputting
        # it to the post processed NetCDF file.
        # Check if file exists.
        try:
            datafile = sorted(glob.glob(self.namelist_dictionary['wrf_directory']
                                        + 'wrfout_%s*' % self.namelist_dictionary['nest']))[0]
            print("Opening %s to gather metadata," % datafile)
        except RuntimeError:
            print("NO %s nest wrfout files in folder %s found.... exiting." %
                  (self.namelist_dictionary['nest'], self.namelist_dictionary['wrf_directory']))

        metadata = Dataset(datafile, 'r')
        meta_data_dict = GetMetaData(metadata)
        # GET LAT / LON AND ADD TO FILE! ##
        lats = metadata.variables['XLAT'][0, :].squeeze()
        lons = metadata.variables['XLONG'][0, :].squeeze()
        metadata.close()

        levels = [float(l) for l in self.namelist_dictionary['levels']]
        # Check if times already exists, if not, get (with) the times!
        if self.WRFTimes is None:
            print("WRF Time values are not defined... Defining them.")
            self.SubsetTime()

        # define output file from the namelist dictionary values.
        outputfile = self.namelist_dictionary['output_path'] + 'wrfoutput_post_%s.nc' \
                     % (self.namelist_dictionary['wrf_pfx'])

        if not self.namelist_dictionary['append']:  # If not appending.
            outputdata = Dataset(outputfile, 'w')  # open up netCDF file for write.

            outputdata.Description = """NetCDF output from wrfoutput
                        atmospheric data interpolated to pressure levels."""

            outputdata.History = "File Created %s" % datetime.now().strftime('%c')

            outputdata.setncatts(meta_data_dict)

            # Create Dimensions #
            outputdata.createDimension('time', None)
            outputdata.createDimension('z', len(levels))
            time = outputdata.createVariable('Time', 'f4', 'time')

            time.units = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')
            time.description = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')

            time[:] = self.WRFTimes

            levs = outputdata.createVariable('Levels', 'f4', 'z')
            levs.units = 'hPa'
            levs.description = 'Model Pressure Levels'
            levs[:] = levels

            outputdata.createDimension('east_west', lons.shape[1])
            outputdata.createDimension('south_north', lons.shape[0])
            longitude = outputdata.createVariable('Longitude', 'f4', ('south_north', 'east_west'))
            latitude = outputdata.createVariable('Latitude', 'f4', ('south_north', 'east_west'))
            longitude[:] = lons
            latitude[:] = lats

            latitude.units = "Degrees North"
            longitude.units = "Degrees East"
        else:
            print("Appending to %s" % outputfile)
            outputdata = Dataset(outputfile, 'r+')  # open up netCDF file for appending

        print("---------------------------")
        print("netCDF output file %s/wrfoutput_post_%s.nc Opened for writing/appending"
              % (self.namelist_dictionary['wrf_directory'], self.namelist_dictionary['wrf_pfx']))
        print("Total # of WRF Files included in post process %i" % len(self.TimeIndexDict))
        print("Total # of WRF time steps %i" % len(self.WRFTimes))
        print("Doing 4D Surface data ")
        print(" -- Doing variables: ")
        for vv in self.namelist_dictionary['chem_vars']:
            print("    -- %s" % vv)
        print("---------------------------")

        # Okay, I think we're finally ready to do this! ##
        tind = 0  # set inital index position for FULL data array.
        # Define a local list of levels to help readability

        for dfx, df in enumerate(self.TimeIndexDict):
            print("Grabbing Data from: %s" % df)
            ncdata = Dataset(df, 'r')  # open for reading
            # Define pressure, since it's used a lot. --> Note that it is in pascals.
            if np.shape(ncdata.variables['P'][:])[0] > 1:
                wrf_pressure = (
                        ncdata.variables['P'][slice(*self.TimeIndexDict[df]), :]
                        + ncdata.variables['PB'][slice(*self.TimeIndexDict[df]), :]).squeeze()
            else:
                wrf_pressure = (
                        ncdata.variables['P'][:]
                        + ncdata.variables['PB'][:])

            for v in self.namelist_dictionary['chem_vars']:
                if v == 'EXTCOF55':  # Do Yellow-band Extinction Coefficient
                    if not all(item in ncdata.variables for item in ['EXTCOF55']):
                        print("!!!WARNING!!!")
                        print("EXTCOF55 are not in thie datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Yellow-band Extinction Coefficient ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    EXT = ncdata.variables['EXTCOF55']
                    try:
                        'bottom_top' in EXT.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()
                    zind = EXT.dimensions.index('bottom_top')

                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(EXT[:])[0] > 1:
                        EXT = EXT[slice(*self.TimeIndexDict[df]), :]
                    else:
                        EXT = EXT[:]

                    EXTInterp = helpers.p_interp(EXT, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:EXTInterp.shape[0] + tind, :] = EXTInterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Yellow-band Extinction: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Black Carbon':  # Do Black Carbon Concentration

                    if not all(item in ncdata.variables for item in ['BC1']):
                        print("!!!WARNING!!!")
                        print("BC1 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Black Carbon Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    BC = ncdata.variables['BC1']

                    try:
                        'bottom_top' in BC.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = BC.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(BC[:])[0] > 1:
                        BC = BC[slice(*self.TimeIndexDict[df]), :]
                    else:
                        BC = BC[:]

                    BCInterp = helpers.p_interp(BC, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:BCInterp.shape[0] + tind, :] = BCInterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Black Carbon Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Organic Carbon':  # Do Organic Carbon Concentration

                    if not all(item in ncdata.variables for item in ['OC1']):
                        print("!!!WARNING!!!")
                        print("OC1 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Organic Carbon Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    OC = ncdata.variables['OC1']

                    try:
                        'bottom_top' in OC.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = OC.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(OC[:])[0] > 1:
                        OC = OC[slice(*self.TimeIndexDict[df]), :]
                    else:
                        OC = OC[:]

                    OCInterp = helpers.p_interp(OC, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:OCInterp.shape[0] + tind, :] = OCInterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Organic Carbon Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Total Dust':  # Do Total Dust Concentration
                    add_dust_bins = False ## initialize a "false" value for adding dust bins.
                    if not all(item in ncdata.variables for item in ['TOT_DUST']):
                        print("!!!WARNING!!!")
                        print("TOT_DUST is not in this datafile,"
                              "Checking if Dust bins 1-5 exist in file...")
                        if not all(item in ncdata.variables for item in ['DUST_1','DUST_2','DUST_3','DUST_4','DUST_5']):
                            print("No dust bins 1-5 either, so I definitely cannot do %s" % v)
                            continue
                        else:
                            print("Found dust bins 1-5, I'll add these to get TOT_DUST!")
                            add_dust_bins=True

                    prg_time = datetime.now()
                    print ("Interpolating Total Dust Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]

                    if add_dust_bins == False:
                        TDUST = ncdata.variables['TOT_DUST']
                    else:
                        TDUST = ncdata.variables['DUST_1'] ## Dummy variable to get the dimensions

                    try:
                        'bottom_top' in TDUST.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TDUST.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if add_dust_bins == False:
                        if np.shape(TDUST[:])[0] > 1:
                            TDUST = TDUST[slice(*self.TimeIndexDict[df]), :]
                        else:
                            TDUST = TDUST[:]
                    else:
                        ## IF ADDING DUST BINS! ##
                        if np.shape(TDUST[:])[0] > 1:
                            TDUST = ncdata.variables['DUST_1'][slice(*self.TimeIndexDict[df]), :] + \
                              ncdata.variables['DUST_2'][slice(*self.TimeIndexDict[df]), :] + \
                              ncdata.variables['DUST_3'][slice(*self.TimeIndexDict[df]), :] + \
                              ncdata.variables['DUST_4'][slice(*self.TimeIndexDict[df]), :] + \
                              ncdata.variables['DUST_5'][slice(*self.TimeIndexDict[df]), :]
                        else:
                            TDUST = ncdata.variables['DUST_1'][:] + \
                              ncdata.variables['DUST_2'][:] + \
                              ncdata.variables['DUST_3'][:] + \
                              ncdata.variables['DUST_4'][:] + \
                              ncdata.variables['DUST_5'][:]

                    TDUSTInterp = helpers.p_interp(TDUST, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:TDUSTInterp.shape[0] + tind, :] = TDUSTInterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Total Dust Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Sea Salt':  # Do Sea-Salt Concentration
                    if not all(item in ncdata.variables for item in ['SEAS_1','SEAS_2','SEAS_3','SEAS_4']):
                        print("!!!WARNING!!!")
                        print("SEAS_1-4 are not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Total Sea Salt Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]

                    TSALT = ncdata.variables['SEAS_1'] ## Dummy variable to get the dimensions

                    try:
                        'bottom_top' in TSALT.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TSALT.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                        ## IF ADDING DUST BINS! ##
                    if np.shape(TSALT[:])[0] > 1:
                        TSALT = ncdata.variables['SEAS_1'][slice(*self.TimeIndexDict[df]), :] + \
                          ncdata.variables['SEAS_2'][slice(*self.TimeIndexDict[df]), :] + \
                          ncdata.variables['SEAS_3'][slice(*self.TimeIndexDict[df]), :] + \
                          ncdata.variables['SEAS_4'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        TSALT = ncdata.variables['SEAS_1'][:] + \
                          ncdata.variables['SEAS'][:] + \
                          ncdata.variables['SEAS'][:] + \
                          ncdata.variables['SEAS'][:]

                    TSALTInterp = helpers.p_interp(TSALT, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:TSALTInterp.shape[0] + tind, :] = TSALTInterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Total Sea Salt Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'SO2':  # Do Sulfur Dioxide Concentration
                    if not all(item in ncdata.variables for item in ['so2']):
                        print("!!!WARNING!!!")
                        print("so2 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Sulfur Dioxide Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    so2 = ncdata.variables['so2']

                    try:
                        'bottom_top' in so2.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = so2.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(so2[:])[0] > 1:
                        so2 = so2[slice(*self.TimeIndexDict[df]), :]
                    else:
                        so2 = so2[:]

                    so2Interp = helpers.p_interp(so2, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:so2Interp.shape[0] + tind, :] = so2Interp  # FINISHED WITH Chemistry VARIABLES!

                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()

                    print("Time to Interpolate Sulfur Dioxide Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'DMS':  # Do Dimethyl Sulfide Concentration

                    if not all(item in ncdata.variables for item in ['dms']):
                        print("!!!WARNING!!!")
                        print("dms is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Dimethyl Sulfide Concentration....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    dms = ncdata.variables['dms']

                    try:
                        'bottom_top' in dms.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = dms.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(dms[:])[0] > 1:
                        dms = dms[slice(*self.TimeIndexDict[df]), :]
                    else:
                        dms = dms[:]

                    dmsInterp = helpers.p_interp(dms, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:dmsInterp.shape[0] + tind, :] = dmsInterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate DMS Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'PM 10':  # Do PM 10 Concentration

                    if not all(item in ncdata.variables for item in ['PM10']):
                        print("!!!WARNING!!!")
                        print("PM10 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating PM 10 Aerosol Particulate Matter Concentration....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    PM10 = ncdata.variables['PM10']

                    try:
                        'bottom_top' in PM10.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = PM10.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(PM10[:])[0] > 1:
                        PM10 = PM10[slice(*self.TimeIndexDict[df]), :]
                    else:
                        PM10 = PM10[:]

                    PM10Interp = helpers.p_interp(PM10, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:PM10Interp.shape[0] + tind, :] = PM10Interp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate PM10 Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'PM 2.5':  # Do PM 2.5 Concentration

                    if not all(item in ncdata.variables for item in ['PM2_5_DRY']):
                        print("!!!WARNING!!!")
                        print("PM2_5_DRY is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating PM 2.5 Dry Aerosol Particulate Matter Concentration....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    PM25 = ncdata.variables['PM2_5_DRY']

                    try:
                        'bottom_top' in PM25.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = PM25.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(PM25[:])[0] > 1:
                        PM25 = PM25[slice(*self.TimeIndexDict[df]), :]
                    else:
                        PM25 = PM25[:]

                    PM25Interp = helpers.p_interp(PM25, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:PM25Interp.shape[0] + tind, :] = PM25Interp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate PM 2.5 Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

            tind += ncdata.variables['Times'][slice(*self.TimeIndexDict[df]), :].shape[
                0]  # add time length to time index!
            ncdata.close()

        for i in data_dict:
            if i in outputdata.variables:
                # VARIABLE ALREADY EXISTS IN NETCDF ###
                print("Variable already exists!  Overwriting ...")
                ncoutvar = outputdata.variables[i]
            else:
                # otherwise, create the variable.
                ncoutvar = outputdata.createVariable(i, 'f4', ('time', 'z', 'south_north', 'east_west'))

            ncoutvar.units = data_dict[i][1]
            ncoutvar.longname = data_dict[i][2]
            ncoutvar[:] = data_dict[i][0][:]
        outputdata.close()
        print("finished output atmospheric chemistry data!")

    def ChemPostSfc(self):
        """Help on built-in function print_namelist in module ChemPostSfc:

        ChemPostSfc(...)
            ChemPostSfc()

            This function creates the (3D) data arrays that get output to the post processed output file.
            Function parses the namelist_dictionary attribute, reads in chemistry variables from the WRF-Chem data
            determined by the namelist_dictionary key: "chem_vars"
            and saves the data to the post processed netCDF file. If append = 0, will created a new file, otherwise will
            save to an existing post processed file.  Note that if append = 1, all data will be overwritten.

            Note that for instances where the variable is 4D atmospheric variable and not a specialized surface variable,
            the surface variable is simply the lowest-level atmospheric variable.

            Returns
            ----------
            Nothing.

        """

        chem_description_dict = {'Black Carbon (Surface)': "Black Carbon Concentration 1",
                                 'SO2 (Surface)': "Sulfur Dioxide Concentration",
                                 'Organic Carbon (Surface)': "Organic Carbon Concentration 1",
                                 'Total Dust (Surface)': "Total Dust Concentration at the surface",
                                 'EXTCOF55 (Surface)': "Extinction Coefficient at 550 nm",
                                 'Aerosol Visibility (Surface)':
                                     "Surface Visibility (550 nm) Considering Aerosol Particles from WRF-Chem Only",
                                 'Aerosol Optical Depth': 'Aerosol Optical Depth at 550nm',
                                 'Ideal Dust' : 'Idealized Dust Flux',
                                 'Dust Scaling Factor': 'Dust Scaling Factor Computed from Static '
                                                        'Analog Dust Map and Idealize Dust Flux',
                                 'Total Dust Emission': "Total Dust Emission",
                                 'PM 2.5 (Surface)': "Aerosol Particulate Matter ( d <2.5 microns) Dry-mass Concentration",
                                 'PM 10 (Surface)': "Aerosol Particulate Matter ( d <10 microns) Dry-mass Concentration",
                                 'DMS (Surface)': 'Dimethyl sulfide ([CH_3]_2S) Concentration',
                                 'Sea Salt (Surface)': 'Total Sea Salt Concentration'}

        chem_unit_dict = {'Black Carbon (Surface)': 'ug kg^-1 dry air', 'EXTCOF55 (Surface)': 'km^-1',
                          'Orgain Carbon (Surface)': 'ug kg^-1 dry air', 'Total Dust (Surface)': 'ug m^-3',
                          'SO2 (Surface)': 'ppmv', 'Aerosol Optical Depth': '-','Ideal Dust':'g m^-2 s^-1',
                          'Total Dust Emission':'g m^-2 s^-1','Aerosol Visibility (Surface)':'km',
                          'PM 10 (Surface)': 'ug m^-3','Dust Scaling Factor':'-',
                          'PM 2.5 (Surface)':'ug m^-3','DMS (Surface)':'ppmv',
                          'Sea Salt (Surface)':'ug kg^-1 dry air'}

        data_dict = {}  # Create one large data dictionary that is used to store data temporarily before outputting
        # it to the post processed NetCDF file.
        # Check if file exists.
        try:
            datafile = sorted(glob.glob(self.namelist_dictionary['wrf_directory']
                                        + 'wrfout_%s*' % self.namelist_dictionary['nest']))[0]
            print("Opening %s to gather metadata," % datafile)
        except RuntimeError:
            print("NO %s nest wrfout files in folder %s found.... exiting." %
                  (self.namelist_dictionary['nest'], self.namelist_dictionary['wrf_directory']))

        metadata = Dataset(datafile, 'r')
        meta_data_dict = GetMetaData(metadata)
        # GET LAT / LON AND ADD TO FILE! ##
        lats = metadata.variables['XLAT'][0, :].squeeze()
        lons = metadata.variables['XLONG'][0, :].squeeze()
        metadata.close()

        levels = [float(l) for l in self.namelist_dictionary['levels']]
        # Check if times already exists, if not, get (with) the times!
        if self.WRFTimes is None:
            print("WRF Time values are not defined... Defining them.")
            self.SubsetTime()

        # define output file from the namelist dictionary values.
        outputfile = self.namelist_dictionary['output_path'] + 'wrfoutput_post_%s.nc' \
                     % (self.namelist_dictionary['wrf_pfx'])

        if not self.namelist_dictionary['append']:  # If not appending.
            outputdata = Dataset(outputfile, 'w')  # open up netCDF file for write.

            outputdata.Description = """NetCDF output from wrfoutput
                        atmospheric data interpolated to pressure levels."""

            outputdata.History = "File Created %s" % datetime.now().strftime('%c')

            outputdata.setncatts(meta_data_dict)

            # Create Dimensions #
            outputdata.createDimension('time', None)
            outputdata.createDimension('z', len(levels))
            time = outputdata.createVariable('Time', 'f4', 'time')

            time.units = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')
            time.description = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')

            time[:] = self.WRFTimes

            levs = outputdata.createVariable('Levels', 'f4', 'z')
            levs.units = 'hPa'
            levs.description = 'Model Pressure Levels'
            levs[:] = levels

            outputdata.createDimension('east_west', lons.shape[1])
            outputdata.createDimension('south_north', lons.shape[0])
            longitude = outputdata.createVariable('Longitude', 'f4', ('south_north', 'east_west'))
            latitude = outputdata.createVariable('Latitude', 'f4', ('south_north', 'east_west'))
            longitude[:] = lons
            latitude[:] = lats

            latitude.units = "Degrees North"
            longitude.units = "Degrees East"
        else:
            print("Appending to %s" % outputfile)
            outputdata = Dataset(outputfile, 'r+')  # open up netCDF file for appending

        print("---------------------------")
        print("netCDF output file %s/wrfoutput_post_%s.nc Opened for writing/appending"
              % (self.namelist_dictionary['wrf_directory'], self.namelist_dictionary['wrf_pfx']))
        print("Total # of WRF Files included in post process %i" % len(self.TimeIndexDict))
        print("Total # of WRF time steps %i" % len(self.WRFTimes))
        print("Doing 3D Surface data ")
        print(" -- Doing variables: ")
        for vv in self.namelist_dictionary['chem_vars']:
            print("    -- %s" % vv)
        print("---------------------------")

        # Okay, I think we're finally ready to do this! ##
        tind = 0  # set initial index position for FULL data array.
        # Define a local list of levels to help readability

        for dfx, df in enumerate(self.TimeIndexDict):
            print("Grabbing Data from: %s" % df)
            ncdata = Dataset(df, 'r')  # open for reading

            for v in self.namelist_dictionary['chem_vars']:
                if v == 'EXTCOF55 (Surface)':  # Do Yellow-band Extinction Coefficient
                    if not all(item in ncdata.variables for item in ['EXTCOF55']):
                        print("!!!WARNING!!!")
                        print("EXTCOF55 is not in thie datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    EXT = ncdata.variables['EXTCOF55']
                    try:
                        'bottom_top' in EXT.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()
                    zind = EXT.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(EXT[:])[0] > 1:
                        EXT = EXT[slice(*self.TimeIndexDict[df]), :]
                    else:
                        EXT = EXT[:]

                    data_dict[v][0][tind:EXT.shape[0] + tind, :] = EXT[:, 0, :]  # Only need first vertical level
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Yellow-band Extinction: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))


                if v == 'Aerosol Visibility (Surface)':  # Do Yellow-band Extinction Coefficient
                    if not all(item in ncdata.variables for item in ['EXTCOF55']):
                        print("!!!WARNING!!!")
                        print("EXTCOF55 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    EXT = ncdata.variables['EXTCOF55']
                    try:
                        'bottom_top' in EXT.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()
                    zind = EXT.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(EXT[:])[0] > 1:
                        EXT = EXT[slice(*self.TimeIndexDict[df]), :]
                    else:
                        EXT = EXT[:]

                    VIS=3.921/EXT
                    data_dict[v][0][tind:VIS.shape[0] + tind, :] = VIS[:, 0, :]  # Only need first vertical level
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Surface Aerosol Visibility: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Black Carbon (Surface)':  # Do Black Carbon Concentration

                    if not all(item in ncdata.variables for item in ['BC1']):
                        print("!!!WARNING!!!")
                        print("BC1 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    BC = ncdata.variables['BC1']

                    try:
                        'bottom_top' in BC.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = BC.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(BC[:])[0] > 1:
                        BC = BC[slice(*self.TimeIndexDict[df]), :]
                    else:
                        BC = BC[:]

                    data_dict[v][0][tind:BC.shape[0] + tind, :] = BC[:, 0, :]  # Only need first vertical level
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Black Carbon Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Organic Carbon (Surface)':  # Do Organic Carbon Concentration

                    if not all(item in ncdata.variables for item in ['OC1']):
                        print("!!!WARNING!!!")
                        print("OC1 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Organic Carbon Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    OC = ncdata.variables['OC1']

                    try:
                        'bottom_top' in OC.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = OC.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(OC[:])[0] > 1:
                        OC = OC[slice(*self.TimeIndexDict[df]), :]
                    else:
                        OC = OC[:]

                    data_dict[v][0][tind:OC.shape[0] + tind, :] = OC[:, 0, :]  # Only need first vertical level
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Organic Carbon Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Total Dust (Surface)':  # Do Total Dust Concentration
                    add_dust_bins = False  ## initialize a "false" value for adding dust bins.
                    if not all(item in ncdata.variables for item in ['TOT_DUST']):
                        print("!!!WARNING!!!")
                        print("TOT_DUST is not in this datafile,"
                              "Checking if Dust bins 1-5 exist in file...")
                        if not all(item in ncdata.variables for item in
                                   ['DUST_1', 'DUST_2', 'DUST_3', 'DUST_4', 'DUST_5']):
                            print("No dust bins 1-5 either, so I definitely cannot do %s" % v)
                            continue
                        else:
                            print("Found dust bins 1-5, I'll add these to get TOT_DUST!")
                            add_dust_bins = True

                    prg_time = datetime.now()
                    print ("Interpolating Total Dust Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    if add_dust_bins == False:
                        TDUST = ncdata.variables['TOT_DUST']
                    else:
                        TDUST = ncdata.variables['DUST_1']  ## Dummy variable to get the dimensions

                    try:
                        'bottom_top' in TDUST.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TDUST.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if add_dust_bins == False:
                        if np.shape(TDUST[:])[0] > 1:
                            TDUST = TDUST[slice(*self.TimeIndexDict[df]), :]
                        else:
                            TDUST = TDUST[:]
                    else:
                        ## IF ADDING DUST BINS! ##
                        if np.shape(TDUST[:])[0] > 1:
                            TDUST = ncdata.variables['DUST_1'][slice(*self.TimeIndexDict[df]), :] + \
                                    ncdata.variables['DUST_2'][slice(*self.TimeIndexDict[df]), :] + \
                                    ncdata.variables['DUST_3'][slice(*self.TimeIndexDict[df]), :] + \
                                    ncdata.variables['DUST_4'][slice(*self.TimeIndexDict[df]), :] + \
                                    ncdata.variables['DUST_5'][slice(*self.TimeIndexDict[df]), :]
                        else:
                            TDUST = ncdata.variables['DUST_1'][:] + \
                                    ncdata.variables['DUST_2'][:] + \
                                    ncdata.variables['DUST_3'][:] + \
                                    ncdata.variables['DUST_4'][:] + \
                                    ncdata.variables['DUST_5'][:]

                    data_dict[v][0][tind:TDUST.shape[0] + tind, :] = TDUST[:, 0, :]  # Only need first vertical level
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Total Dust Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Sea Salt (Surface)':  # Do Total Sea Salt Concentration
                    if not all(item in ncdata.variables for item in ['SEAS_1', 'SEAS_2', 'SEAS_3', 'SEAS_4']):
                        print("!!!WARNING!!!")
                        print("SEAS_1-4 are not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Total Sea Salt Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    TSALT = ncdata.variables['SEAS_1']  ## Dummy variable to get the dimensions

                    try:
                        'bottom_top' in TSALT.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TSALT.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    ## IF ADDING DUST BINS! ##
                    if np.shape(TSALT[:])[0] > 1:
                        TSALT = ncdata.variables['SEAS_1'][slice(*self.TimeIndexDict[df]), :] + \
                                ncdata.variables['SEAS_2'][slice(*self.TimeIndexDict[df]), :] + \
                                ncdata.variables['SEAS_3'][slice(*self.TimeIndexDict[df]), :] + \
                                ncdata.variables['SEAS_4'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        TSALT = ncdata.variables['SEAS_1'][:] + \
                                ncdata.variables['SEAS_2'][:] + \
                                ncdata.variables['SEAS_3'][:] + \
                                ncdata.variables['SEAS_4'][:]

                    data_dict[v][0][tind:TSALT.shape[0] + tind, :] = TSALT[:, 0, :]  # Only need first vertical level
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Total Sea Salt Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))


                if v == 'Total Dust Emission':  # Do Total Dust Emission
                    ### NOTE THIS NEEDS A SPECIAL FUNCTION HERE TO ACCOUNT FOR THE FACT THAT ###
                    ### EDUST HAS DIFFERENT UNITS FOR ALL DUST SCHEMES.
                    ### THIS IS CONTROLLED BY THE USER IN THE NAMELIST WITH THE FLAG: dust_scheme
                    ### the dust_scheme flag corresponds to the wrf options
                        ### -> 1 = GOCART
                        ### -> 3 = AFWA
                        ### -> 4 = UoC
                    ### If TOT_EDUST is in wrfout file, then it MUST BE AFWA, will adjust accordingly
                    ### but return a warning.
                    ### Otherwise, it will use value in the namelist.

                    ### FIRST CHECK IF DUST SCHEME IS EVEN AN ALLOWABLE NUMBER, if it is not, throw warning and set to 1
                    if int(self.namelist_dictionary['dust_scheme']) not in (1,3,4):
                        print("!!!WARNING!!!")
                        print("dust_scheme %i is not allowed, the following schemes are allowed:"%
                              int(self.namelist_dictionary['dust_scheme']))
                        print(" dust_scheme = 1 -> GOCART (converts kg to g m^-2 s^-1")
                        print(" dust_scheme = 3 -> AFWA-GOCART (no unit conversion")
                        print(" dust_scheme = 4 -> UoC (converts ug m^-2 s^-1 to g m^-2 s^-1)")
                        print("Defaulting to GOCART (dust_scheme = 1")
                        self.namelist_dictionary['dust_scheme']=1
                        print("")

                    add_dust_bins = False  ## initialize a "false" value for adding dust bins.
                    if not all(item in ncdata.variables for item in ['TOT_EDUST']):
                        print("!!!WARNING!!!")
                        print("TOT_EDUST is not in this datafile,"
                              "Checking if Dust bins 1-5 exist in file...")
                        if not all(item in ncdata.variables for item in
                                   ['EDUST1', 'EDUST2', 'EDUST3', 'EDUST4', 'EDUST5']):
                            print("No dust bins 1-5 either, so I definitely cannot do %s" % v)
                            continue
                        else:
                            print("Found dust bins 1-5, I'll add these to get ETOT_DUST!")
                            add_dust_bins = True

                    prg_time = datetime.now()
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    if add_dust_bins == False:
                        print("Since TOT_EDUST exists in this file, I assume it's the AFWA dust emission scheme")
                        print("Note, as of WRF-Chem V 4.1, TOT_EDUST was only available in the AFWA scheme.")
                        if self.namelist_dictionary['dust_scheme'] != 3:
                            print("WARNING, dust_scheme in the namelist is NOT equal to 3 (AFWA)")
                            print("I'm overriding this and setting it to 3, but I wanted to make you aware...")
                            self.namelist_dictionary['dust_scheme']=3
                        TDUST = ncdata.variables['TOT_EDUST']
                    else:
                        TDUST = ncdata.variables['EDUST1']  ## Dummy variable to get the dimensions

                    if add_dust_bins == False:
                        if np.shape(TDUST[:])[0] > 1:
                            TDUST = TDUST[slice(*self.TimeIndexDict[df]), :]
                        else:
                            TDUST = TDUST[:]
                    else:
                        ## IF ADDING DUST BINS! ##
                        if np.shape(TDUST[:])[0] > 1:
                            TDUST = ncdata.variables['EDUST1'][slice(*self.TimeIndexDict[df]), :].squeeze() + \
                                    ncdata.variables['EDUST2'][slice(*self.TimeIndexDict[df]), :].squeeze() + \
                                    ncdata.variables['EDUST3'][slice(*self.TimeIndexDict[df]), :].squeeze() + \
                                    ncdata.variables['EDUST4'][slice(*self.TimeIndexDict[df]), :].squeeze() + \
                                    ncdata.variables['EDUST5'][slice(*self.TimeIndexDict[df]), :].squeeze()
                        else:
                            TDUST = ncdata.variables['EDUST1'][:].squeeze() + \
                                    ncdata.variables['EDUST2'][:].squeeze() + \
                                    ncdata.variables['EDUST3'][:].squeeze() + \
                                    ncdata.variables['EDUST4'][:].squeeze() + \
                                    ncdata.variables['EDUST5'][:].squeeze()

                    ## DO UNIT CONVERSION! ##

                    if int(self.namelist_dictionary['dust_scheme']) ==1:
                        print("dust_scheme = 1 (GOCART)")
                        print("converting GOCART EDUST from kg to g m^-2 s^-1")
                        DX=float(ncdata.DX)
                        DY=float(ncdata.DY)
                        dxdy=DX*DY
                        DT=float(ncdata.DT)
                        TDUST = TDUST*1000./(DT*dxdy)

                    elif int(self.namelist_dictionary['dust_scheme']) ==4:
                        print('dust_scheme = 4 (UoC)')
                        print("converting UoC EDUST from ug m^-2 s^-1 to g m^-2 s^-1")
                        TDUST=TDUST*1.E-6 ## micro grams to grams

                    elif int(self.namelist_dictionary['dust_scheme']) == 3:
                        print("dust_scheme = 3 (AFWA)")
                        print("NO Conversion needed.")

                    data_dict[v][0][tind:TDUST.shape[0] + tind, :] = TDUST[:]  # Only need first vertical level
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Total Dust Emission: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'SO2 (Surface)':  # Do Sulfur Dioxide Concentration
                    if not all(item in ncdata.variables for item in ['so2']):
                        print("!!!WARNING!!!")
                        print("so2 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Sulfur Dioxide Concentration ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    so2 = ncdata.variables['so2']

                    try:
                        'bottom_top' in so2.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = so2.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(so2[:])[0] > 1:
                        so2 = so2[slice(*self.TimeIndexDict[df]), :]
                    else:
                        so2 = so2[:]

                    data_dict[v][0][tind:so2.shape[0] + tind, :] = so2[:, 0, :]  # Only need first vertical level

                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Sulfur Dioxide Concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'PM 10 (Surface)':  # Do Surface PM10
                    if not all(item in ncdata.variables for item in ['PM10']):
                        print("!!!WARNING!!!")
                        print("PM10 is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating PM10 values ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    pm10 = ncdata.variables['PM10']

                    try:
                        'bottom_top' in pm10.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = pm10.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(pm10[:])[0] > 1:
                        pm10 = pm10[slice(*self.TimeIndexDict[df]), :]
                    else:
                        pm10 = pm10[:]

                    data_dict[v][0][tind:pm10.shape[0] + tind, :] = pm10[:, 0, :]  # Only need first vertical level

                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get PM 2.5 concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'PM 2.5 (Surface)':  # Do Surface PM2.5
                    if not all(item in ncdata.variables for item in ['PM2_5_DRY']):
                        print("!!!WARNING!!!")
                        print("PM2_5_DRY is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating PM 2.5 values ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    pm25 = ncdata.variables['PM2_5_DRY']

                    try:
                        'bottom_top' in pm25.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = pm25.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(pm25[:])[0] > 1:
                        pm25 = pm25[slice(*self.TimeIndexDict[df]), :]
                    else:
                        pm25 = pm25[:]

                    data_dict[v][0][tind:pm25.shape[0] + tind, :] = pm25[:, 0, :]  # Only need first vertical level

                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get PM 2.5 concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'DMS (Surface)':  # Do Surface Dimethyl Sulfide
                    if not all(item in ncdata.variables for item in ['dms']):
                        print("!!!WARNING!!!")
                        print("dms is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating dms values ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    dms = ncdata.variables['dms']

                    try:
                        'bottom_top' in dms.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = dms.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(dms[:])[0] > 1:
                        dms = dms[slice(*self.TimeIndexDict[df]), :]
                    else:
                        dms = dms[:]

                    data_dict[v][0][tind:dms.shape[0] + tind, :] = dms[:, 0, :]  # Only need first vertical level

                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get DMS concentration: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))


                if v == 'Ideal Dust':  # Do Ideal Dust Concentration
                    if not all(item in ncdata.variables for item in ['DUST_IDEAL']):
                        print("!!!WARNING!!!")
                        print("DUST_IDEAL is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    dust_id = ncdata.variables['DUST_IDEAL']

                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(dust_id[:])[0] > 1:
                        dust_id = dust_id[slice(*self.TimeIndexDict[df]), :]
                    else:
                        dust_id = dust_id[:]

                    data_dict[v][0][tind:dust_id.shape[0] + tind, :] = dust_id[:]  # Only need first vertical level

                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Idealized Dust: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))


                if v == 'Dust Scaling Factor':  # Do Dust Scaling Factor
                    ## Note that this function requires dust potential / Analog Dust has already been added to the output
                    ## file.  Super special variable, that will only ever be called by a select few people.
                    ## Think of it as something off of the "secret" menu at McDonalds: The assumption here is that if you
                    ## are calling this variable, you know how it's computed and how to use the post processor.
                    ## Requires DUST_IDEAL is in wrfout data, and Dust Analog has already been saved to the post processed file
                    if not all(item in ncdata.variables for item in ['DUST_IDEAL']):
                        print("!!!WARNING!!!")
                        print("DUST_IDEAL is not in this datafile, cannot do %s" % v)
                        continue

                    if not all(item in outputdata.variables for item in ['Analog Dust', '2D Erodibility']):
                        print("!!!WARNING!!!")
                        print("Analog Dust or 2D Erodibility have not been saved to %s, cannot do %s" % (outputfile,v))
                        continue

                    print("Special Variable From ERDC-GEO: "
                          "                       Analog Dust\n"
                          "Dust Scaling Factor =  -----------\n"
                          "                    Idealized Dust Flux")
                    print("NOTE:  Filling values outside of ERDC-GEO with AFWA Erodibility Scaling Factor")
                    prg_time = datetime.now()
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    dust_id = ncdata.variables['DUST_IDEAL']

                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(dust_id[:])[0] > 1:
                        dust_id = dust_id[slice(*self.TimeIndexDict[df]), :]
                    else:
                        dust_id = dust_id[:]

                    dust_potent2D=outputdata.variables['Analog Dust'][:]
                    erod2D=outputdata.variables['2D Erodibility'][:]

                    mask=np.ma.masked_less_equal(dust_id,1E-6).mask ## need to save mask for use in a minute
                    ## take care of pesky divide by zero warnings and convert to milligrams.
                    dust_id=np.ma.masked_less_equal(dust_id,1E-6).filled(1E-6)*1000.

                    ##Need to match static array to time-variable array sizes:
                    ##Note: Not the most efficient method, revisit in future versions?
                    dust_potent=np.zeros_like(dust_id)
                    erod3D=np.zeros_like(dust_id)
                    for i in range(dust_potent.shape[0]):
                        dust_potent[i,:]=dust_potent2D[:]
                        erod3D[i,:]=erod2D[:]

                    output=np.ma.masked_where(dust_potent[:]>9.9,dust_potent[:]/dust_id[:]).filled(erod3D)
                    ##Now reapply the masked less than 1E-6 ideal dust mask defined above

                    output=np.ma.masked_array(output,mask=mask).filled(erod3D[:])


                    data_dict[v][0][tind:dust_id.shape[0] + tind, :] = output[:]

                    output=None
                    dust_potent=None
                    dust_id=None
                    dust_potent2D=None
                    erod2D=None
                    erod3D=None

                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Get Dust Scaling Factor: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Aerosol Optical Depth':  # Do Total AOD

                    if not all(item in ncdata.variables for item in ['EXTCOF55', 'PH', 'PHB']):
                        print("!!!WARNING!!!")
                        print("EXTCOF55,PH,PHB are not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Doing Aerosol Optical Depth... ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]
                    EXT = ncdata.variables['EXTCOF55']

                    try:
                        'bottom_top' in EXT.dimensions
                    except KeyError:
                        print("Utoh, %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = EXT.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(EXT[:])[0] > 1:
                        HGT = ncdata.variables['PH'][slice(*self.TimeIndexDict[df]), :] + \
                              ncdata.variables['PHB'][slice(*self.TimeIndexDict[df]), :]
                        EXT = EXT[slice(*self.TimeIndexDict[df]), :]
                    else:
                        HGT = ncdata.variables['PH'][:] + ncdata.variables['PHB'][:]
                        EXT = EXT[:]
                    # ---> GEO POTENTIAL (m2/s2) --> height (m) --> hgt (km)
                    HGT = (HGT / helpers.constants['g']) / 1000.

                    # AOD = ∫ ext*dz ---> Extinction Coeff (km-1) * dz (km) -- > sum = column integrated
                    WRFAOD = np.sum(EXT * np.diff(HGT, axis=zind), axis=zind)
                    data_dict[v][0][tind:WRFAOD.shape[0] + tind, :] = WRFAOD  # FINISHED WITH ATMOS VARIABLES!
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(chem_unit_dict[v])
                        data_dict[v].append(chem_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate AOD: %.2f seconds" % ((end_prg_time - prg_time).total_seconds()))

            tind += ncdata.variables['Times'][slice(*self.TimeIndexDict[df]), :].shape[
                0]  # add time length to time index!
            ncdata.close()

        for i in data_dict:
            if i in outputdata.variables:
                # VARIABLE ALREADY EXISTS IN NETCDF FILE --> Overwrite it!###
                print("Variable already exists!  Overwriting ...")
                ncoutvar = outputdata.variables[i]
            else:
                # otherwise, create the variable.
                ncoutvar = outputdata.createVariable(i, 'f4', ('time', 'south_north', 'east_west'))

            ncoutvar.units = data_dict[i][1]
            ncoutvar.longname = data_dict[i][2]
            ncoutvar[:] = data_dict[i][0][:]
        outputdata.close()
        print("finished output atmospheric chemistry data!")

    def AtmosPost(self):
        """Help on built-in function print_namelist in module AtmosPost:

        AtmosPost(...)
            AtmosPost()

            This function creates the (4D) data arrays that get output to the post processed output file.
            Function parses the namelist_dictionary attribute, reads in chemistry variables from the WRF data
            determined by the namelist_dictionary key: "atmos_vars"
            and saves the data to the post processed netCDF file. If append = 0, will created a new file, otherwise will
            save to an existing post processed file.  Note that if append = 1, all data will be overwritten.

            Note that for instances where the variable is 4D atmospheric variable and not a specialized surface variable,
            the surface variable is simply the lowest-level atmospheric variable.

            Returns
            ----------
            Nothing.

        """

        # Define specific atmospheric variables ##
        atmos_description_dict = {"U velocity": 'Wind Component in the Zonal (x) Direction',
                                  "V velocity": 'Wind Component in the Meridonal (y) Direction',
                                  "W velocity": 'Wind Component in the Vertical (z) Direction',
                                  'Geopotential': 'Geopotential Height',
                                  'Temperature': 'Temperature',
                                  'Potential Temperature': 'Potential Temperature',
                                  'Eq Potential Temperature': 'Equivalent Potential Temperature (Stull, 1988 Approximation)',
                                  'Relative humidity': 'Relative Humidity',
                                  'Cloud water mixing ratio': 'Cloud water mass mixing ratio',
                                  'Rain mixing ratio': 'Rain water mass mixing ratio',
                                  'Cloud Ice mixing ratio': 'Cloud Ice mass mixing ratio',
                                  'Snow mixing ratio': 'Snow mass mixing ratio',
                                  'Cloud Fraction': 'Sub grid cloud fraction',
                                  'Visibility': "Estimated visibility",
                                  'Radar Reflectivity': "Radar Reflectivity"
                                  }

        atmos_unit_dict = {"U velocity": 'm s^-1', "V velocity": 'm s^-1', "W velocity": 'm s^-1',
                           'Geopotential': 'm^2 s^-2',
                           'Temperature': 'K','Potential Temperature': 'K',
                           'Eq Potential Temperature': 'K','Relative humidity': '%',
                           'Cloud water mixing ratio': 'kg kg^-1',
                           'Rain mixing ratio': 'kg kg^-1',
                           'Cloud Ice mixing ratio': 'kg kg^-1',
                           'Snow mixing ratio': 'kg kg^-1',
                           'Cloud Fraction': '-', 'Visibility': 'km',
                           'Radar Reflectivity': "dBZ"}

        data_dict = {}  # Create one large data dictionary that is used to store data temporarily before outputting
        # it to the post processed NetCDF file.
        # Check if file exists.
        try:
            datafile = sorted(glob.glob(self.namelist_dictionary['wrf_directory']
                                        + 'wrfout_%s*' % self.namelist_dictionary['nest']))[0]
            print("Opening %s to gather metadata," % datafile)
        except RuntimeError:
            print("NO %s nest wrfout files in folder %s found.... exiting." %
                  (self.namelist_dictionary['nest'], self.namelist_dictionary['wrf_directory']))

        metadata = Dataset(datafile, 'r')
        meta_data_dict = GetMetaData(metadata)
        # GET LAT / LON AND ADD TO FILE! ##
        lats = metadata.variables['XLAT'][0, :].squeeze()
        lons = metadata.variables['XLONG'][0, :].squeeze()
        metadata.close()

        # Check if times already exists, if not, get (with) the times!
        if self.WRFTimes is None:
            print("WRF Time values are not defined... Defining them.")
            self.SubsetTime()

        levels = [float(l) for l in self.namelist_dictionary['levels']]
        # define output file from the namelist dictionary values.
        outputfile = self.namelist_dictionary['output_path'] + 'wrfoutput_post_%s.nc' \
                     % (self.namelist_dictionary['wrf_pfx'])

        if not self.namelist_dictionary['append']:  # If not appending.
            outputdata = Dataset(outputfile, 'w')  # open up netCDF file for write.

            outputdata.Description = """NetCDF output from wrfoutput
                 atmospheric data interpolated to pressure levels."""

            outputdata.History = "File Created %s" % datetime.now().strftime('%c')

            outputdata.setncatts(meta_data_dict)

            # Create Dimensions #
            outputdata.createDimension('time', None)
            outputdata.createDimension('z', len(levels))
            time = outputdata.createVariable('Time', 'f4', 'time')

            time.units = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')
            time.description = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')

            time[:] = self.WRFTimes

            levs = outputdata.createVariable('Levels', 'f4', 'z')
            levs.units = 'hPa'
            levs.description = 'Model Pressure Levels'
            levs[:] = levels

            outputdata.createDimension('east_west', lons.shape[1])
            outputdata.createDimension('south_north', lons.shape[0])
            longitude = outputdata.createVariable('Longitude', 'f4', ('south_north', 'east_west'))
            latitude = outputdata.createVariable('Latitude', 'f4', ('south_north', 'east_west'))
            longitude[:] = lons
            latitude[:] = lats

            latitude.units = "Degrees North"
            longitude.units = "Degrees East"
        else:
            print("Appending to %s" % outputfile)
            outputdata = Dataset(outputfile, 'r+')  # open up netCDF file for appending

        print("---------------------------")
        print("netCDF output file %s/wrfoutput_post_%s.nc Opened for writing/appending"
              % (self.namelist_dictionary['wrf_directory'], self.namelist_dictionary['wrf_pfx']))
        print("Total # of WRF Files included in post process %i" % len(self.TimeIndexDict))
        print("Total # of WRF time steps %i" % len(self.WRFTimes))
        print("Doing 4D Surface data ")
        print(" -- Doing variables: ")
        for vv in self.namelist_dictionary['atmos_vars']:
            print("    -- %s" % vv)
        print("---------------------------")
        # Okay, I think we're finally ready to do this! ##
        tind = 0  # set inital index position for FULL data array.
        # Define a local list of levels to help readability
        for dfx, df in enumerate(self.TimeIndexDict):
            print("Grabbing Data from: %s" % df)
            ncdata = Dataset(df, 'r')  # open for reading
            # Define pressure, since it's used a lot. --> Note that it is in pascals.
            # NOTE --> Need conditional statement to account for the possibility that the time dimension is
            # a single dimension and, the slicing method will not load the data.
            if np.shape(ncdata.variables['P'][:])[0] > 1:
                wrf_pressure = (
                        ncdata.variables['P'][slice(*self.TimeIndexDict[df]), :]
                        + ncdata.variables['PB'][slice(*self.TimeIndexDict[df]), :]).squeeze()
            else:
                wrf_pressure = (
                        ncdata.variables['P'][:]
                        + ncdata.variables['PB'][:])

            for v in self.namelist_dictionary['atmos_vars']:
                if self.namelist_dictionary['append'] and v in outputdata.variables:
                    # ONLY DO Something if you are not appending and the variable doesn't already exist!
                    print("This variable already in file, skipping...")
                    continue

                if v == 'U velocity':  # DO U WIND!!
                    if 'U' not in ncdata.variables:
                        print("!!!WARNING!!!")
                        print("U Is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()

                    print ("Interpolating U Velocity ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    UU = ncdata.variables['U']
                    try:
                        'bottom_top' in UU.dimensions
                    except KeyError:
                        print("Utoh ... %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = UU.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(UU[:])[0] > 1:
                        UU = (UU[slice(*self.TimeIndexDict[df]), :, :, :-1]
                              + UU[slice(*self.TimeIndexDict[df]), :, :, 1:]) / 2.  # UNSTAGGER
                        if self.namelist_dictionary['scale_uv']:
                            VV = ncdata.variables['V'][slice(*self.TimeIndexDict[df]), :]
                            VV = (VV[:, :, :-1, :] + VV[:, :, 1:, :]) / 2.  # UNSTAGGER

                            for j in range(UU.shape[1]):
                                UU[:, j, :] = (
                                        UU[:, j, :] * ncdata.variables['COSALPHA'][slice(*self.TimeIndexDict[df]), :]
                                        - VV[:, j, :] * ncdata.variables['SINALPHA'][slice(*self.TimeIndexDict[df]),
                                                        :])
                    else:
                        UU = (UU[:, :, :, :-1] + UU[:, :, :, 1:]) / 2.  # UNSTAGGER
                        if self.namelist_dictionary['scale_uv']:
                            VV = ncdata.variables['V'][:]
                            VV = (VV[:, :, :-1, :] + VV[:, :, 1:, :]) / 2.  # UNSTAGGER

                            for j in range(UU.shape[1]):
                                UU[:, j, :] = (UU[:, j, :] * ncdata.variables['COSALPHA'][:]
                                               - VV[:, j, :] * ncdata.variables['SINALPHA'][:])

                    UIterp = helpers.p_interp(UU, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:UIterp.shape[0] + tind, :] = UIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])

                    end_prg_time = datetime.now()
                    print("Time to Interpolate U velocity: %.2f seconds" % ((end_prg_time - prg_time).total_seconds()))

                if v == 'V velocity':  # DO V WIND!!

                    if 'V' not in ncdata.variables:
                        print("!!!WARNING!!!")
                        print("V Is not in this datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()

                    print ("Interpolating V Velocity ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    VV = ncdata.variables['V']

                    try:
                        'bottom_top' in VV.dimensions
                    except KeyError:
                        print("Utoh %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = VV.dimensions.index('bottom_top')
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(VV[:])[0] > 1:
                        VV = (VV[slice(*self.TimeIndexDict[df]), :, :-1, :]
                              + VV[slice(*self.TimeIndexDict[df]), :, 1:, :]) / 2.  # UNSTAGGER

                        if self.namelist_dictionary['scale_uv']:
                            UU = ncdata.variables['U'][slice(*self.TimeIndexDict[df]), :]
                            UU = (UU[:, :, :, :-1] + UU[:, :, :, 1:]) / 2.  # UNSTAGGER x dimension.
                            for j in range(VV.shape[1]):
                                VV[:, j, :] = (
                                        VV[:, j, :] * ncdata.variables['COSALPHA'][slice(*self.TimeIndexDict[df]),
                                                      :]
                                        + UU[:, j, :] * ncdata.variables['SINALPHA'][slice(*self.TimeIndexDict[df]),
                                                        :]
                                )

                    else:
                        VV = (VV[:, :, :-1, :]
                              + VV[:, :, 1:, :]) / 2.  # UNSTAGGER

                        if self.namelist_dictionary['scale_uv']:
                            UU = ncdata.variables['U'][:]
                            UU = (UU[:, :, :, :-1] + UU[:, :, :, 1:]) / 2.  # UNSTAGGER x dimension.
                            for j in range(VV.shape[1]):
                                VV[:, j, :] = (VV[:, j, :] * ncdata.variables['COSALPHA'][:]
                                               + UU[:, j, :] * ncdata.variables['SINALPHA'][:])

                    VIterp = helpers.p_interp(VV, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:VIterp.shape[0] + tind, :] = VIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])

                    end_prg_time = datetime.now()
                    print("Time to Interpolate V velocity: %.2f seconds" % ((end_prg_time - prg_time).total_seconds()))

                if v == 'Temperature':  # DO Temperature!!

                    if 'T' not in ncdata.variables:
                        print("!!!WARNING!!!")
                        print("T Is not in this datafile, cannot do %s" % v)
                        continue
                    prg_time = datetime.now()
                    print ("Interpolating Temperature ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    TT = ncdata.variables['T']

                    try:
                        'bottom_top' in TT.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TT.dimensions.index('bottom_top')
                    # Convert to temperature

                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(TT[:])[0] > 1:
                        TT = helpers.Theta2TmpK(TT[slice(*self.TimeIndexDict[df]), :] + 300., wrf_pressure / 100.)
                    else:
                        TT = helpers.Theta2TmpK(TT[:] + 300., wrf_pressure / 100.)
                    TIterp = helpers.p_interp(TT, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:TIterp.shape[0] + tind, :] = TIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Temperature: %.2f seconds" % ((end_prg_time - prg_time).total_seconds()))

                if v == 'Potential Temperature':  # DO Potential Temperature!!

                    if 'T' not in ncdata.variables:
                        print("!!!WARNING!!!")
                        print("T Is not in this datafile, cannot do %s" % v)
                        continue
                    prg_time = datetime.now()
                    print ("Interpolating Potential Temperature ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    TT = ncdata.variables['T']

                    try:
                        'bottom_top' in TT.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TT.dimensions.index('bottom_top')
                    # Convert to temperature

                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(TT[:])[0] > 1:
                        TT = TT[slice(*self.TimeIndexDict[df]), :] + 300.
                    else:
                        TT = TT[:] + 300.
                    TIterp = helpers.p_interp(TT, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:TIterp.shape[0] + tind, :] = TIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Potential Temperature: %.2f seconds" % ((end_prg_time - prg_time).total_seconds()))

                if v == 'Eq Potential Temperature':  # DO Potential Temperature!!

                    if not all(item in ncdata.variables for item in ['T', 'QVAPOR']):
                        print("!!!WARNING!!!")
                        print("T and QVAPOR are not in thie datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Equivalent Potential Temperature ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    TT = ncdata.variables['T']

                    try:
                        'bottom_top' in TT.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TT.dimensions.index('bottom_top')
                    # Convert to temperature

                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.
                    if np.shape(TT[:])[0] > 1:
                        TT = helpers.Theta2TmpK(TT[slice(*self.TimeIndexDict[df]), :] + 300., wrf_pressure / 100.)
                        QV = ncdata.variables['QVAPOR'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        TT = helpers.Theta2TmpK(TT[:] + 300., wrf_pressure / 100.)
                        QV = ncdata.variables['QVAPOR'][:]

                    THETAE=(TT+helpers.constants['Lv']/helpers.constants['Cp']*QV)*\
                       (helpers.constants['pref']/(wrf_pressure / 100.))**(helpers.constants['Rd']
                                                                           /helpers.constants['Cp'])

                    TIterp = helpers.p_interp(THETAE, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:TIterp.shape[0] + tind, :] = TIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Equivalent Potential Temperature: %.2f seconds" % ((end_prg_time - prg_time).total_seconds()))


                if v == 'Reflectivity':  # DO Reflectivity!!

                    if 'REFL_10CM' not in ncdata.variables:
                        print("!!!WARNING!!!")
                        print("REFL_10CM Is not in this datafile, cannot do %s" % v)
                        continue
                    prg_time = datetime.now()
                    print ("Interpolating Reflectivity ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    REF = ncdata.variables['REFL_10CM']

                    try:
                        'bottom_top' in REF.dimensions
                    except KeyError:
                        print("Utoh. %s does not have a valid vertical coordinate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = REF.dimensions.index('bottom_top')
                    if np.shape(REF[:])[0] > 1:
                        REF = REF[slice(*self.TimeIndexDict[df]), :]
                    else:
                        REF = REF[:]
                    # NOTE --> Need conditional statement to account for the possibility that the time dimension is
                    # a single dimension and, the slicing method will not load the data.

                    # Convert to temperature
                    REFinterp = helpers.p_interp(REF, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:REFinterp.shape[0] + tind, :] = REFinterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Reflectivity: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Relative humidity':  # DO RELATIVE Humidity!!

                    if not all(item in ncdata.variables for item in ['T', 'QVAPOR']):
                        print("!!!WARNING!!!")
                        print("T and QVAPOR are not in thie datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Relative Humidity ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    TT = ncdata.variables['T']
                    QV = ncdata.variables['QVAPOR']

                    try:
                        'bottom_top' in TT.dimensions
                        'bottom_top' in QV.dimensions
                    except KeyError:
                        print("Utoh %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TT.dimensions.index('bottom_top')

                    if np.shape(TT[:])[0] > 1:
                        TT = helpers.Theta2TmpK(TT[slice(*self.TimeIndexDict[df]), :] + 300.,
                                                wrf_pressure / 100.)  # Convert to temperature
                        ES = helpers.VaporPressure(TT, phase="liquid")  # GET Saturation Vapor Pressure
                        E = helpers.MixR2VaporPress(QV[slice(*self.TimeIndexDict[df]), :],
                                                    wrf_pressure/100.)  # GET Vapor Pressure

                    else:
                        TT = helpers.Theta2TmpK(TT[:] + 300., wrf_pressure / 100.)  # Convert to temperature
                        ES = helpers.VaporPressure(TT, phase="liquid")  # GET Saturation Vapor Pressure
                        E = helpers.MixR2VaporPress(QV[:], wrf_pressure/100.)  # GET Vapor Pressure

                    RH = (E / ES) * 100.  # Get relative humidity in %

                    RHIterp = helpers.p_interp(RH, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:RHIterp.shape[0] + tind, :] = RHIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Relative Humidity: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Visibility':  # DO VISIBILITY

                    if not all(
                            item in ncdata.variables for item in ['T', 'QVAPOR', 'QCLOUD', 'QICE', 'QRAIN', 'QSNOW']):
                        print("!!!WARNING!!!")
                        print("T, QVAPOR,QCLOUD,QICE,QRAIN,QSNOW are not in thie datafile, cannot do %s" % v)
                        continue

                    prg_time = datetime.now()
                    print ("Interpolating Visibility ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    TT = ncdata.variables['T']
                    QV = ncdata.variables['QVAPOR']
                    QC = ncdata.variables['QCLOUD']
                    QI = ncdata.variables['QICE']
                    QR = ncdata.variables['QRAIN']
                    QS = ncdata.variables['QSNOW']

                    try:
                        'bottom_top' in TT.dimensions
                        'bottom_top' in QV.dimensions
                        'bottom_top' in QC.dimensions
                        'bottom_top' in QI.dimensions
                        'bottom_top' in QR.dimensions
                        'bottom_top' in QS.dimensions
                    except KeyError:
                        print("Utoh %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = TT.dimensions.index('bottom_top')

                    if np.shape(TT[:])[0] > 1:
                        TT = helpers.Theta2TmpK(TT[slice(*self.TimeIndexDict[df]), :] + 300.,
                                                wrf_pressure / 100.)  # Convert to temperature
                        QV = QV[slice(*self.TimeIndexDict[df]), :]
                        QC = QC[slice(*self.TimeIndexDict[df]), :]
                        QI = QI[slice(*self.TimeIndexDict[df]), :]
                        QR = QR[slice(*self.TimeIndexDict[df]), :]
                        QS = QS[slice(*self.TimeIndexDict[df]), :]

                    else:
                        TT = helpers.Theta2TmpK(TT[:] + 300., wrf_pressure / 100.)  # Convert to temperature

                        QV = QV[:]
                        QC = QC[:]
                        QI = QI[:]
                        QR = QR[:]
                        QS = QS[:]

                    rho_a = helpers.AirDensity(TT, wrf_pressure / 100., QV, moist=True)
                    VIS = helpers.wrf_visibility(rho_a, QC, QV, QI, QR, QS)

                    VISIterp = helpers.p_interp(VIS, wrf_pressure, lev=levels, zind=zind)
                    data_dict[v][0][tind:VISIterp.shape[0] + tind, :] = VISIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])
                    end_prg_time = datetime.now()

                    # Clear memory by setting (probably) large arrays = None.
                    TT = None
                    QC = None
                    QV = None
                    QI = None
                    QR = None
                    QS = None

                    print("Time to Interpolate Visibility: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'Geopotential':  # DO Geopotential Height!!
                    if not all(item in ncdata.variables for item in ['PH', 'PHB']):
                        print("!!!WARNING!!!")
                        print("PH and PHB are not in thie datafile, cannot do %s" % v)
                        continue
                    prg_time = datetime.now()
                    print ("Interpolating Geopotential height ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    HH = ncdata.variables['PH']

                    try:
                        'bottom_top_stag' in HH.dimensions  # Remember that geopotential is staggered in z.
                    except KeyError:
                        print("WARNING %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = HH.dimensions.index('bottom_top_stag')
                    if np.shape(HH[:])[0] > 1:
                        HH = (HH[slice(*self.TimeIndexDict[df]), :]
                              + ncdata.variables['PHB'][slice(*self.TimeIndexDict[df]), :]).squeeze()
                        HH = (HH[:, :-1, :] + HH[:, 1:, :]) / 2.
                    else:
                        HH = (HH[:] + ncdata.variables['PHB'][:])
                        HH = (HH[:, :-1, :] + HH[:, 1:, :]) / 2.

                    HIterp = helpers.p_interp(HH, wrf_pressure, lev=levels, zind=zind)
                    HIterp = HIterp / helpers.constants['g']

                    data_dict[v][0][tind:HIterp.shape[0] + tind, :] = HIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate Geopotential Height: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))

                if v == 'W velocity':  # DO Geopotential Height!!
                    if 'W' not in ncdata.variables:
                        print("!!!WARNING!!!")
                        print("W Is not in thie datafile, cannot do %s" % v)
                        continue
                    prg_time = datetime.now()
                    print ("Interpolating Vertical (W) Velocity ....")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes), len(levels)] + list(np.shape(lons)))]
                    WW = ncdata.variables['W']

                    try:
                        'bottom_top_stag' in WW.dimensions
                    except KeyError:
                        print("Utoh %s does not have a valid vertical coordiate!" % v)
                        print("You must choose a variable with the vertical coordinate: bottom_top")
                        print("Something appears wrong with wrfout file ... Exiting script")
                        ncdata.close()

                    zind = WW.dimensions.index('bottom_top_stag')

                    if np.shape(WW[:])[0] > 1:
                        WW = (WW[slice(*self.TimeIndexDict[df]), :])
                    else:
                        WW = WW[:]
                    WW = (WW[:, :-1, :] + WW[:, 1:, :]) / 2.

                    WIterp = helpers.p_interp(WW, wrf_pressure, lev=levels, zind=zind)

                    data_dict[v][0][tind:WIterp.shape[0] + tind, :] = WIterp
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(atmos_unit_dict[v])
                        data_dict[v].append(atmos_description_dict[v])
                    end_prg_time = datetime.now()
                    print("Time to Interpolate W Velocity: %.2f seconds" % ((end_prg_time - prg_time).total_seconds()))
            # add time length to time index!
            tind += ncdata.variables['Times'][slice(*self.TimeIndexDict[df]), :].shape[0]
            ncdata.close()

        for i in data_dict:
            if i in outputdata.variables:
                # VARIABLE ALREADY EXISTS IN NETCDF ###
                print("Variable already exists!  Overwriting ...")
                ncoutvar = outputdata.variables[i]
            else:
                ncoutvar = outputdata.createVariable(i, 'f4', ('time', 'z', 'south_north', 'east_west'))
            ncoutvar.units = data_dict[i][1]
            ncoutvar.longname = data_dict[i][2]
            ncoutvar[:] = data_dict[i][0][:]

        outputdata.close()
        print("finished output atmospheric data!")

    def SfcPost(self):
        """Help on built-in function print_namelist in module SfcPost:

                SfcPost(...)
                    SfcPost()

                    This function creates the (3D) data arrays that get output to the post processed output file.
                    Function parses the namelist_dictionary attribute, reads in surface variables from the WRF data
                    determined by the namelist_dictionary key: "sfc_vars"
                    and saves the data to the post processed netCDF file. If append = 0, will created a new file, otherwise will
                    save to an existing post processed file.  Note that if append = 1, all data will be overwritten.

                    Note that for instances where the variable is 4D atmospheric variable and not a specialized surface variable,
                    the surface variable is simply the lowest-level atmospheric variable.

                    Note that in the namelist dictionary a special value of sfc_vars = 'all' will save ALL available
                    surface variables to the post processed file.

                    Returns
                    ----------
                    Nothing.

                """

        sfc_description_dict = {'10 metre U wind component': '10 meter zonal (x) wind component',
                                '10 metre V wind component': '10 meter meridonal (y) wind component',
                                '2 metre temperature': '2 AGL Temperature',
                                '2 metre dewpoint temperature': '2 AGL Dewpoint Temperature',
                                'Mean sea level pressure': 'Mean sea level pressure',
                                'Surface pressure': 'Surface Pressure',
                                'Total column water vapour': 'Total column water vapour (precipitable water)',
                                'Gridscale Precipitation rate': 'Grid-scale precipitation rate',
                                'Convective Precipitation rate': 'Precipitation rate from the convective parameterization',
                                'Planetary boundary layer height': "Height of the boundary layer top",
                                'Visibility (Surface)': "Horizontal Visibility",
                                'Snow cover fraction': 'Fractional Snow Cover within Gridcell',
                                'Gridscale Snowfall rate': 'Grid-scale Snowfall Rate',
                                'Friction Velocity': 'Surface Friction Velocity',
                                'z0': 'Surface Roughness Length',
                                'Outgoing Longwave Radiation': 'Upwelling Longwave Radiation At the Top of the Atmosphere',
                                '1km AGL Reflectivity' : "Simulated Radar Reflectivity at 1km Above Ground Level",
                                'Normal Adjusted Friction Velocity':
                                    'Normalized Surface Friction Velocity Scaled by Windspeed',
                                'Soil Moisture': 'Volumetric Soil Moisture'}

        sfc_unit_dict = {'10 metre U wind component': 'm s ^-1', '10 metre V wind component': 'm s ^-1',
                         '2 metre temperature': 'K', '2 metre dewpoint temperature': 'K',
                         'Mean sea level pressure': 'Pa',
                         'Total column water vapour': 'kg m^-2', 'Gridscale Precipitation rate': 'kg m^-2 s^-1',
                         'Convective Precipitation rate': 'kg m^-2 s^-1', 'Planetary boundary layer height': 'm',
                         'Visibility (Surface)': 'km', 'Snow cover fraction': '-',
                         'Gridscale Snowfall rate': 'kg m^-2 s^-1',
                         'Friction Velocity': 'm s^-1', 'z0': 'm', 'Outgoing Longwave Radiation': 'W m^-2',
                         '1km AGL Reflectivity':'dBZ','Normal Adjusted Friction Velocity':'m/s','Surface pressure':'Pa',
                         'Soil Moisture':'m^3 m^-3'}

        data_dict = {}
        create_data_dict = True

        if 'all' in self.namelist_dictionary['sfc_vars']:
            self.namelist_dictionary['sfc_vars']=[v for v in sfc_description_dict.keys()]

        ## Check if file exists.
        try:
            datafile = sorted(glob.glob(self.namelist_dictionary['wrf_directory']
                                        + 'wrfout_%s*' % self.namelist_dictionary['nest']))[0]
            print("Opening %s to gather metadata," % datafile)
        except RuntimeError:
            print("NO %s nest wrfout files in folder %s found.... exiting." %
                  (self.namelist_dictionary['nest'], self.namelist_dictionary['wrf_directory']))

        metadata = Dataset(datafile, 'r')
        meta_data_dict = GetMetaData(metadata)
        # GET LAT / LON AND ADD TO FILE! ##
        lats = metadata.variables['XLAT'][0, :].squeeze()
        lons = metadata.variables['XLONG'][0, :].squeeze()
        metadata.close()

        # Check if times already exists, if not, get (with) the times!
        if self.WRFTimes is None:
            print("WRF Time values are not defined... Defining them.")
            self.SubsetTime()

        levels = [float(l) for l in self.namelist_dictionary['levels']]
        # define output file from the namelist dictionary values.
        outputfile = self.namelist_dictionary['output_path'] + 'wrfoutput_post_%s.nc' \
                     % (self.namelist_dictionary['wrf_pfx'])

        if not self.namelist_dictionary['append']:  # If not appending.
            outputdata = Dataset(outputfile, 'w')  # open up netCDF file for write.

            outputdata.Description = """NetCDF output from wrfoutput
                         atmospheric data interpolated to pressure levels."""

            outputdata.History = "File Created %s" % datetime.now().strftime('%c')

            outputdata.setncatts(meta_data_dict)

            # Create Dimensions #
            outputdata.createDimension('time', None)
            outputdata.createDimension('z', len(levels))
            time = outputdata.createVariable('Time', 'f4', 'time')

            time.units = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')
            time.description = "Hours since %s" % self.TimeStart.strftime('%H:%M:%S %d-%m-%Y')

            time[:] = self.WRFTimes

            levs = outputdata.createVariable('Levels', 'f4', 'z')
            levs.units = 'hPa'
            levs.description = 'Model Pressure Levels'
            levs[:] = levels

            outputdata.createDimension('east_west', lons.shape[1])
            outputdata.createDimension('south_north', lons.shape[0])
            longitude = outputdata.createVariable('Longitude', 'f4', ('south_north', 'east_west'))
            latitude = outputdata.createVariable('Latitude', 'f4', ('south_north', 'east_west'))
            longitude[:] = lons
            latitude[:] = lats

            latitude.units = "Degrees North"
            longitude.units = "Degrees East"
        else:
            print("Appending to %s" % outputfile)
            outputdata = Dataset(outputfile, 'r+')  # open up netCDF file for appending

        print("---------------------------")
        print("netCDF output file %s/wrfoutput_post_%s.nc Opened for writing/appending"
              % (self.namelist_dictionary['wrf_directory'], self.namelist_dictionary['wrf_pfx']))
        print("Total # of WRF Files included in post process %i" % len(self.TimeIndexDict))
        print("Total # of WRF time steps %i" % len(self.WRFTimes))
        print("Doing 3D Surface data ")
        print(" -- Doing variables: ")
        for vv in self.namelist_dictionary['sfc_vars']:
            print("    -- %s" % vv)
        print("---------------------------")

        # Okay, I think we're finally ready to do this! ##
        tind = 0  # set inital index position for FULL data array.
        # Define a local list of levels to help readability
        for dfx, df in enumerate(self.TimeIndexDict):
            print("Grabbing Data from: %s" % df)
            ncdata = Dataset(df, 'r')  # open for reading

            if np.shape(ncdata.variables['P'][:])[0] > 1:
                wrf_pressure = (
                        ncdata.variables['P'][slice(*self.TimeIndexDict[df]), 0, :]
                        + ncdata.variables['PB'][slice(*self.TimeIndexDict[df]), 0, :]).squeeze()
            else:
                wrf_pressure = (
                        ncdata.variables['P'][:, 0, :]
                        + ncdata.variables['PB'][:, 0, :]).squeeze()

            for v in self.namelist_dictionary['sfc_vars']:

                if v == 'Mean sea level pressure':
                    if not all(item in ncdata.variables for item in ['PH', 'PHB', 'T', 'QVAPOR']):
                        print("!!!WARNING!!!")
                        print("PH,PHB,T, and QVAPOR are not in this datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Sea Level Pressure")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    if np.shape(ncdata.variables['PH'][:])[0] > 1:
                        HGT = (ncdata.variables['PH'][slice(*self.TimeIndexDict[df]), 0, :]
                               + ncdata.variables['PHB'][slice(*self.TimeIndexDict[df]), 0, :]
                               + ncdata.variables['PH'][slice(*self.TimeIndexDict[df]), 1, :]
                               + ncdata.variables['PHB'][slice(*self.TimeIndexDict[df]), 1, :]) / 2.

                        T = ncdata.variables['T'][slice(*self.TimeIndexDict[df]), 0, :] + 300.
                        Q = ncdata.variables['QVAPOR'][slice(*self.TimeIndexDict[df]), 0, :]
                    else:
                        HGT = (ncdata.variables['PH'][:, 0, :] + ncdata.variables['PHB'][:, 0, :] +
                               ncdata.variables['PH'][:, 1, :] + ncdata.variables['PHB'][:, 1, :]) / 2.
                        T = ncdata.variables['T'][slice(*self.TimeIndexDict[df]), 0, :] + 300.
                        Q = ncdata.variables['QVAPOR'][slice(*self.TimeIndexDict[df]), 0, :]

                    T = helpers.Theta2TmpK(T, wrf_pressure.squeeze() / 100.)
                    MSLP = helpers.mslp(wrf_pressure.squeeze()/ 100.,T.squeeze(), HGT.squeeze()/9.81, Q.squeeze())

                    data_dict[v][0][tind:MSLP.shape[0] + tind, :] = MSLP
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == '10 metre U wind component':
                    if not all(item in ncdata.variables for item in ['U10', 'V10', 'COSALPHA', 'SINALPHA']):
                        print("!!!WARNING!!!")
                        print("U10,V10, COSALPHA,and SINALPHA, are not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing 10 meter U wind")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    if np.shape(ncdata.variables['U10'][:])[0] > 1:
                        UU = ncdata.variables['U10'][slice(*self.TimeIndexDict[df]), :]
                        if self.namelist_dictionary['scale_uv']:  ## IF correcting map factor.
                            VV = ncdata.variables['V10'][slice(*self.TimeIndexDict[df]), :]
                            UU = UU * ncdata.variables['COSALPHA'][slice(*self.TimeIndexDict[df]), :] \
                                 - VV * ncdata.variables['SINALPHA'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        UU = ncdata.variables['U10'][:]
                        if self.namelist_dictionary['scale_uv']:  ## IF correcting map factor.
                            VV = ncdata.variables['V10'][:]
                            UU = UU * ncdata.variables['COSALPHA'][:] - VV * ncdata.variables['SINALPHA'][:]

                    data_dict[v][0][tind:UU.shape[0] + tind, :] = UU
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == '10 metre V wind component':
                    if not all(item in ncdata.variables for item in ['U10', 'V10', 'COSALPHA', 'SINALPHA']):
                        print("!!!WARNING!!!")
                        print("U10,V10, COSALPHA,and SINALPHA, are not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing 10 meter V wind")

                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    if np.shape(ncdata.variables['V10'][:])[0] > 1:
                        VV = ncdata.variables['V10'][slice(*self.TimeIndexDict[df]), :]
                        if self.namelist_dictionary['scale_uv']:  ## IF correcting map factor.
                            UU = ncdata.variables['U10'][slice(*self.TimeIndexDict[df]), :]
                            VV = VV * ncdata.variables['COSALPHA'][slice(*self.TimeIndexDict[df]), :] \
                                 + UU * ncdata.variables['SINALPHA'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        VV = ncdata.variables['V10'][:]
                        if self.namelist_dictionary['scale_uv']:  ## IF correcting map factor.
                            UU = ncdata.variables['U10'][:]
                            VV = VV * ncdata.variables['COSALPHA'][:] + UU * ncdata.variables['SINALPHA'][:]

                    data_dict[v][0][tind:VV.shape[0] + tind, :] = VV
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == '2 metre temperature':
                    if not all(item in ncdata.variables for item in ['T2']):
                        print("!!!WARNING!!!")
                        print("T2 is not in this datafile, cannot do %s" % v)
                        continue
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    print("Grabbing 2 meter temperature")
                    TT = self.standard_sfc_var('T2', ncdata, df)

                    data_dict[v][0][tind:TT.shape[0] + tind, :] = TT
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Soil Moisture':
                    if not all(item in ncdata.variables for item in ['SMOIS', 'ZS']):
                        print("!!!WARNING!!!")
                        print("Variable SMOIS and ZS are not in this datafile, cannot do %s" % v)
                        continue

                    ## SOIL MOISTURE IS SPECIAL, WILL LOOP THROUGH ALL HEIGHTS AND MAKE SOIL
                    if np.ndim(ncdata.variables['ZS'][:]) > 1:
                        ZS = ncdata['ZS'][0, :].squeeze() * 100.  ## Time / Soil Height in cm
                    else:
                        ZS = ncdata['ZS'][:].squeeze() * 100.
                    if dfx == 0:
                        ## NOW NEED TO MAKE NEW VARIABLES FOR EACH SOIL LEVEL.
                        for z in ZS:
                            vsm = '%s at %i cm below surface'%(v,z)
                            data_dict[vsm] = [np.zeros([len(self.WRFTimes)] +list(np.shape(lons)))]


                    print("Grabbing Soil Moisture")
                    SM = self.standard_sfc_var('SMOIS', ncdata, df)

                    for zdx, z in enumerate(ZS):
                        vsm = '%s at %i cm below surface' % (v, z)
                        data_dict[vsm][0][tind:SM.shape[0] + tind, :] = SM[:,zdx,:]
                        if len(data_dict[vsm]) == 1:
                            data_dict[vsm].append(sfc_unit_dict[v])
                            data_dict[vsm].append(sfc_description_dict[v]+ '%i cm below the surface' %z)

                if v == '2 metre dewpoint temperature':
                    if not all(item in ncdata.variables for item in ['Q2', 'PSFC']):
                        print("!!!WARNING!!!")
                        print("Q2 and PSFC are not in this datafile, cannot do %s" % v)
                        continue

                    print("Grabbing 2 meter dew point temperature")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    if np.shape(ncdata.variables['Q2'][:])[0] > 1:
                        Q2 = ncdata.variables['Q2'][slice(*self.TimeIndexDict[df]), :]
                        press = ncdata.variables['PSFC'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        Q2 = ncdata.variables['Q2'][:]
                        press = ncdata.variables['PSFC'][:]

                    E2 = helpers.MixR2VaporPress(Q2, press/100.)
                    TD = helpers.DewPoint(E2)

                    data_dict[v][0][tind:TT.shape[0] + tind, :] = TD
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Total column water vapour':
                    if not all(item in ncdata.variables for item in ['P', 'PB', 'PSFC', 'QVAPOR']):
                        print("!!!WARNING!!!")
                        print("P and PSFC, PB, and QVAPOR are not in this datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Total column water vapor (a.k.a., precipitable water)")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    if np.shape(ncdata.variables['P'][:])[0] > 1:
                        PRES1 = ncdata.variables['P'][slice(*self.TimeIndexDict[df]), :] \
                                + ncdata.variables['PB'][slice(*self.TimeIndexDict[df]), :]  # Add base-state pressure
                        PSFC = ncdata.variables['PSFC'][slice(*self.TimeIndexDict[df]), :]
                        PTOTAL = np.concatenate((PSFC[:, None, :, :], PRES1),
                                                axis=1)  # ADD SURFACE PRESSURE TO Pressure
                        DP = np.diff(PTOTAL, axis=1)  # GET DP, used in all calculations!
                        QV = ncdata.variables['QVAPOR'][slice(*self.TimeIndexDict[df]), :]  # WATER VAPOR MIXING RATIO!
                    else:
                        PRES1 = ncdata.variables['P'][:] + ncdata.variables['PB'][:]  # Add base-state pressure
                        PSFC = ncdata.variables['PSFC'][:]
                        PTOTAL = np.concatenate((PSFC[:, None, :, :], PRES1),
                                                axis=1)  # ADD SURFACE PRESSURE TO Pressure
                        DP = np.diff(PTOTAL, axis=1)  # GET DP, used in all calculations!
                        QV = ncdata.variables['QVAPOR'][:]  # WATER VAPOR MIXING RATIO!

                    pwater = -1. / (9.81 * 1000.) * np.sum(QV * DP, axis=1) * 1000.  # To get mm
                    data_dict[v][0][tind:pwater.shape[0] + tind, :] = pwater
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])



                if v == '1km AGL Reflectivity':

                    prg_time=datetime.now()
                    ## 1km above ground reflectivity.
                    if not all(item in ncdata.variables for item in
                               ['PHB', 'PH', 'REFL_10CM']):
                        print("!!!WARNING!!!")
                        print("PHB', 'PH', 'REFL_10CM are not in this datafile"
                              ", cannot do %s" % v)
                        continue
                    print("Computing Reflectivity at 1km AGL")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    REF=ncdata.variables['REFL_10CM']
                    zind = REF.dimensions.index('bottom_top')

                    if np.shape(ncdata.variables['PHB'][:])[0] > 1:
                        REF = ncdata.variables['REFL_10CM'][slice(*self.TimeIndexDict[df]), :, :]
                        HGT = ncdata.variables['PH'][slice(*self.TimeIndexDict[df]), :] + \
                              ncdata.variables['PHB'][slice(*self.TimeIndexDict[df]), :]

                    else:
                        REF = ncdata.variables['REFL_10CM'][:, :, :]
                        HGT = ncdata.variables['PH'][:, :] + \
                              ncdata.variables['PHB'][:, :]

                    HGT = (HGT[:, :-1, :] + HGT[:, 1:, :]) / 2.  ## Unstagger the Height variable
                    HGT=HGT / helpers.constants['g'] ## Convert to geometric height.

                    REF1km = helpers.z_interp(REF, HGT, lev=1000, zind=zind)
                    data_dict[v][0][tind:REF1km.shape[0] + tind, :] = REF1km.squeeze()
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                    end_prg_time = datetime.now()
                    print("Time to Interpolate 1km Reflectivity: %.2f seconds" % (
                        (end_prg_time - prg_time).total_seconds()))
                if v == 'Visibility (Surface)':

                    if not all(item in ncdata.variables for item in
                               ['P', 'PB', 'T', 'QVAPOR', 'QCLOUD', 'QICE', 'QRAIN', 'QSNOW']):
                        print("!!!WARNING!!!")
                        print("P, PB, T, QCLOUD, QRAIN, QSNOW, QICE, and and QVAPOR are not in this datafile"
                              ", cannot do %s" % v)
                        continue
                    print("Computing Visibility at first model level")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    if np.shape(ncdata.variables['P'][:])[0] > 1:
                        PRES = ncdata.variables['P'][slice(*self.TimeIndexDict[df]), 0, :] + ncdata.variables['PB'][
                                                                                    slice(*self.TimeIndexDict[df]), 0,
                                                                                    :]  # Add base-state pressure
                        TT = ncdata.variables['T'][slice(*self.TimeIndexDict[df]), 0, :] + 300.
                        QV = ncdata.variables['QVAPOR'][slice(*self.TimeIndexDict[df]), 0, :]
                        QC = ncdata.variables['QCLOUD'][slice(*self.TimeIndexDict[df]), 0, :]
                        QI = ncdata.variables['QICE'][slice(*self.TimeIndexDict[df]), 0, :]
                        QS = ncdata.variables['QSNOW'][slice(*self.TimeIndexDict[df]), 0, :]
                        QR = ncdata.variables['QRAIN'][slice(*self.TimeIndexDict[df]), 0, :]

                    else:
                        PRES = ncdata.variables['P'][:, 0, :] + ncdata.variables['PB'][:, 0,
                                                                :]  # Add base-state pressure
                        TT = ncdata.variables['T'][:, 0, :] + 300.
                        QV = ncdata.variables['QVAPOR'][:, 0, :]
                        QC = ncdata.variables['QCLOUD'][:, 0, :]
                        QI = ncdata.variables['QICE'][:, 0, :]
                        QS = ncdata.variables['QSNOW'][:, 0, :]
                        QR = ncdata.variables['QRAIN'][:, 0, :]

                    TMP = helpers.Theta2TmpK(TT, PRES, pref=100000.)  # Convert Theta to Tempearture
                    Tvirtual = helpers.VirtualTemp(TMP, QV)  # GETS VIRTUAL TEMPERATURE
                    rho_a = helpers.AirDensity(TMP, PRES, QV)
                    VIS = helpers.wrf_visibility(rho_a, QC, QV, QI, QR, QS)
                    data_dict[v][0][tind:VIS.shape[0] + tind, :] = VIS
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Gridscale Precipitation rate':
                    if not all(item in ncdata.variables for item in ['RAINNC']):
                        print("!!!WARNING!!!")
                        print("RAINNC is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Gridscale Precipitation rate --> Warning assumes output frequency is consant!")

                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    SF = ncdata.variables['RAINNC']
                    time_locator = SF.dimensions.index('Time')

                    if np.shape(SF[:])[0] > 0:
                        PCP = ncdata.variables['RAINNC'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        PCP = ncdata.variables['RAINNC'][:]

                    data_dict[v][0][tind:PCP.shape[0] + tind, :] = PCP
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Gridscale Snowfall rate':
                    if not all(item in ncdata.variables for item in ['SNOWNC']):
                        print("!!!WARNING!!!")
                        print("SNOWNC is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Gridscale Snowfall rate --> Warning assumes output frequency is constant!")

                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    SF = ncdata.variables['SNOWNC']
                    time_locator = SF.dimensions.index('Time')

                    if np.shape(SF[:])[0] > 1:
                        PCP = ncdata.variables['SNOWNC'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        PCP = ncdata.variables['SNOWNC'][:]

                    data_dict[v][0][tind:PCP.shape[0] + tind, :] = PCP
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Convective Precipitation rate':
                    if not all(item in ncdata.variables for item in ['RAINC']):
                        print("!!!WARNING!!!")
                        print("RAINC is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Convective Precipitation rate --> Warning assumes output frequency is consant!")

                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    RAIN = ncdata.variables['RAINC']
                    time_locator = RAIN.dimensions.index('Time')
                    if np.shape(RAIN[:])[0] > 0:
                        PCP = ncdata.variables['RAINC'][slice(*self.TimeIndexDict[df]), :]
                    else:
                        PCP = ncdata.variables['RAINC'][:]

                    data_dict[v][0][tind:PCP.shape[0] + tind, :] = PCP
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Planetary boundary layer height':
                    if not all(item in ncdata.variables for item in ['PBLH']):
                        print("!!!WARNING!!!")
                        print("PBLH is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Planetary boundary layer height")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    PBL = self.standard_sfc_var('PBLH', ncdata, df)
                    data_dict[v][0][tind:PBL.shape[0] + tind, :] = PBL
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])


                if v == 'Surface pressure':
                    if not all(item in ncdata.variables for item in ['PSFC']):
                        print("!!!WARNING!!!")
                        print("PSFC is not in this datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Surface Pressure")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    PSFC = self.standard_sfc_var('PSFC', ncdata, df)
                    data_dict[v][0][tind:PSFC.shape[0] + tind, :] = PSFC
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Snow cover fraction':
                    if not all(item in ncdata.variables for item in ['SNOWC']):
                        print("!!!WARNING!!!")
                        print("SNOWC is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Snow cover fraction")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    SNOWC = self.standard_sfc_var('SNOWC', ncdata, df)

                    data_dict[v][0][tind:SNOWC.shape[0] + tind, :] = SNOWC
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Friction Velocity':
                    if not all(item in ncdata.variables for item in ['UST']):
                        print("!!!WARNING!!!")
                        print("UST is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Friction Velocity")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    UST = self.standard_sfc_var('UST', ncdata, df)

                    data_dict[v][0][tind:UST.shape[0] + tind, :] = UST
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'z0':
                    if not all(item in ncdata.variables for item in ['ZNT']):
                        print("!!!WARNING!!!")
                        print("ZNT is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Surface Roughness")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    ZNT = self.standard_sfc_var('ZNT', ncdata, df)

                    data_dict[v][0][tind:ZNT.shape[0] + tind, :] = ZNT
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Normal Adjusted Friction Velocity':
                    if not all(item in ncdata.variables for item in ['TEDS_UST']):
                        print("!!!WARNING!!!")
                        print("TEDS_UST is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Normal Adjusted Friction Velocity")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    UST = self.standard_sfc_var('TEDS_UST', ncdata, df)

                    data_dict[v][0][tind:UST.shape[0] + tind, :] = UST
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'z0':
                    if not all(item in ncdata.variables for item in ['ZNT']):
                        print("!!!WARNING!!!")
                        print("ZNT is not in thie datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Surface Roughness")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    ZNT = self.standard_sfc_var('ZNT', ncdata, df)

                    data_dict[v][0][tind:ZNT.shape[0] + tind, :] = ZNT
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

                if v == 'Outgoing Longwave Radiation':
                    if not all(item in ncdata.variables for item in ['OLR']):
                        print("!!!WARNING!!!")
                        print("OLR is not in this datafile, cannot do %s" % v)
                        continue
                    print("Grabbing Outgoing Longwave Radiation")
                    if dfx == 0:
                        data_dict[v] = [np.zeros([len(self.WRFTimes)] + list(np.shape(lons)))]

                    OLR = self.standard_sfc_var('OLR', ncdata, df)

                    data_dict[v][0][tind:OLR.shape[0] + tind, :] = OLR
                    if len(data_dict[v]) == 1:
                        data_dict[v].append(sfc_unit_dict[v])
                        data_dict[v].append(sfc_description_dict[v])

            tind += ncdata.variables['Times'][slice(*self.TimeIndexDict[df]), :].shape[
                0]  ## add time lenght to time index!
            ncdata.close()

        for i in data_dict:

            if i in ['Gridscale Snowfall rate','Convective Precipitation rate','Gridscale Precipitation rate']:
                if len(self.WRFTimes) > 1:
                    denom=(self.WRFTimes[1]-self.WRFTimes[0])*3600.
                    data_dict[i][0][:]=np.gradient(data_dict[i][0][:])[0]/denom
                else:
                    print("!!! WARNING !!!")
                    print("You are only outputting 1 time, therefore I have nothing to compare the wrf accumulated"
                          "%s variable against to get a rate, therefore %s will be an accumulated precip from "
                          "the beginning of the simulation, and not a rate."%(i,i))

            if i in outputdata.variables:
                ### VARIABLE ALREADY EXISTS IN NETCDF ###
                print("Variable already exists!  Overwriting ...")
                ncoutvar = outputdata.variables[i]
            else:
                if i == 'Soil Moisture':
                    ## special case, need soil height dimension added!
                    ncoutvar = outputdata.createVariable(i, 'f4', ('time', 'soil_layer','south_north', 'east_west'))
                else:
                    ncoutvar = outputdata.createVariable(i, 'f4', ('time', 'south_north', 'east_west'))
            ncoutvar.units = data_dict[i][1]
            ncoutvar.longname = data_dict[i][2]
            ncoutvar[:] = data_dict[i][0][:]

        outputdata.close()
        append = 1
        return append
        print("finished output Surface data!")

    def ParseNamelist(self, namelist_file='wrf_post_namelist.txt'):

        """Help on built-in function print_namelist in module ParseNamelist:

                ParseNamelist(...)
                    ParseNamelist([namelist_file])

                    This is a user-callable function that parses the namelist and puts it into the "namelist_dictionary"
                    attribute

                    Notes/Rules for namelist formatting:
                    NO QUOTES in the namelist!  No "" or '', I know, it seems counter intuitive to
                    send strings with spaces without quotation marks, but if you include quotes anywhere
                    in an uncommented section of the namelist, you're going to have a bad time.

                    Namelist comments are designiated by "#" This function will ignore
                    anything written on a line following a #.

                    All list values (e.g., variables) should be comma separated. separate things.

                    You can write anything into the namelist you want,

                    Only values that are containted within the fixed "lookforkeys" will be read into the dictionary.
                    the 'list_keys' list indicates that the function is looking for multiple values,
                    if it doesn't find multiple values for list keys, it will just return a list with length 1.

                    lookforkeys:
                        - levels (list): pressure levels (in hPa) to interpolate atmospheric data to
                        - atmos_vars (list): standard "vanilla" WRF output variable names for 4D atmospheric data
                        - sfc_vars (list): standard "vanilla" WRF output variable names for 3D surface data
                        - chem_vars (list): WRF-Chem output variable names for both 4D atmospheric and 3D surface chem data
                        - static_vars (list): static variable names

                        - wrf_directory: Path to folder where raw wrfout files are stored
                        - wrf_pfx: User defined prefix to the output post processed NetCDF file
                        - time_format: time format (e.g., %Y-%m-%d_%H:%M:%S) in python
                            strftime (https://strftime.org) format corresponding to "time_min" and "time_max" namelist values
                        - output_path: Path to folder where the wrf post processed file will be saved
                        - time_min: lower time boundary for time subsetting
                        - time_max: higher time boundary for time subsetting
                        - scale_uv: Boolean Flag (0=False,1=True) to convert the grid-relative u/v wind components to Earth relative u/v
                            Important for high-latitude
                        - nest: nest identifier to separate the raw wrfinput and wrfout files
                        - append: Boolean Flag (0=False,1=True) to determine if the data should be added to a preexisting
                            post process file, or if a new post processed file should be created
                        - dust_scheme: Special flag only relevant when outputing "Total Dust Emission" chem var.
                            Matches wrf-chem namelist option.  Required to ensure that the Total Dust Emission units are
                            correct.

                    Parameters
                    ----------
                    namelist_dictionary: string (optional), path to namelist.txt file with the namelist information
                        default value: wrf_post_namelist.txt

                    Returns
                    ----------
                    namelist_dictionary

                """

        lookforkeys = ['levels', 'atmos_vars', 'sfc_vars', 'chem_vars',
                       'wrf_directory', 'wrf_pfx', 'time_format', 'output_path', 'time_min',
                       'time_max', 'scale_uv', 'nest', 'append', 'static_vars','dust_scheme']

        list_keys = ['levels', 'atmos_vars', 'sfc_vars', 'chem_vars', 'static_vars']

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



def GetMetaData(metadata,wrfinput=False):
    """Help on built-in function print_namelist in module GetMetaData:

        GetMetaData(...)
            GetMetaData([wrfinput])

            This is a helper function to populate the metadata dictionary for the wrf post processed file.
            Automatically called when creating a post processed file.

            Parameters
            -----------
            wrfinput: Boolean (optional): Flag indicating the input data is from a wrfinput file

            Returns
            ----------
            dictionary contating metadata information that will be saved to the post processed netcdf file.
    """
    meta_data_dict = {}

    meta_data_dict['DX'] = metadata.DX
    meta_data_dict['DY'] = metadata.DY
    meta_data_dict['TITLE'] = metadata.TITLE
    meta_data_dict['MP_PHYSICS'] = metadata.MP_PHYSICS

    meta_data_dict["RA_LW_PHYSICS"] = metadata.RA_LW_PHYSICS
    meta_data_dict["RA_SW_PHYSICS"] = metadata.RA_SW_PHYSICS
    meta_data_dict["SF_SFCLAY_PHYSICS"] = metadata.SF_SFCLAY_PHYSICS
    meta_data_dict["SF_SURFACE_PHYSICS"] = metadata.SF_SURFACE_PHYSICS
    meta_data_dict["BL_PBL_PHYSICS"] = metadata.BL_PBL_PHYSICS
    meta_data_dict["CU_PHYSICS"] = metadata.CU_PHYSICS
    # meta_data_dict["NX"]=metadata.WEST-EAST_GRID_DIMENSION
    # meta_data_dict["NY"]=metadata.SOUTH-NORTH_GRID_DIMENSION
    # meta_data_dict["NZ"]=metadata.BOTTOM-TOP_GRID_DIMENSION
    meta_data_dict['MAP_PROJ_CHAR'] = metadata.MAP_PROJ_CHAR
    meta_data_dict['CEN_LAT'] = metadata.CEN_LAT
    meta_data_dict["CEN_LON"] = metadata.CEN_LON
    meta_data_dict["STAND_LON"] = metadata.STAND_LON
    meta_data_dict['DT'] = metadata.DT
    meta_data_dict['LAND_USE_CLASS'] = metadata.MMINLU

    if wrfinput == False:
        # note these meta-data files are not in wrfinput files, so they are only called with wrfinput == false.
        meta_data_dict['RADT'] = metadata.RADT
        meta_data_dict['RAIN_BUCKET'] = metadata.BUCKET_MM
        meta_data_dict['ENERGY_BUCKET'] = metadata.BUCKET_J


    return meta_data_dict
