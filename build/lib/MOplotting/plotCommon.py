# coding=utf-8
"""A file with functions that are common to most, if not all plotting classes within geoplot.
   Admittedly, this whole class of functions could use a lot more documentation."""

from cartopy import crs as ccrs
from cartopy import feature as cfeature
import numpy as np
from matplotlib import pyplot as plt

import glob as glob
import sys
import os

def KtoC(tmp):
    """Standard Kelvin to Celsius Temperature conversion"""
    return tmp - 273.16

def KtoF(tmp):
    """Standard Kelvin to Fahrenheit Temperature conversion"""
    tmp = tmp - 273.16
    return 9. / 5. * tmp + 32.

def AddMetars(ax,path='./',stns='all',time='2019-02-24 18:00',tfmt='%Y-%m-%d %H:%M',pvar='ALTM',
               vari='Temperature',vectors=True,spd='SPD',drct='DIR',lon='LON',lat='LAT',**kwargs):
    """This function adds meters data to an already specific map.  Designed to be used as a stand alone function.

       looks for station data that has been formatted and produced using the get metar part of this script
       in a specified folder.  If stns == 'all', it will grab all stations within that folder, otherwise it takes a
       "list" of stations to use.

       Assumes longitude and latitude are in decimal degrees.
       This function performs NO unit conversions in MOOSE V 1.0, future version may include unit support.

    Parameters
        ----------
        ax : matplotlib subplot object with a cartopy projection attachment, subplot to plot METARS on.
        path (default: './'): string, path to METAR DATA
        stns (default: 'all'): string ('all') or list of stations to plot from within the 'path' folder
        time (default: '2019-02-24 18:00'): string, plot METAR data nearest to this time
        tfmt (default: '%Y-%m-%d %H:%M'): string, format for time variable following strftime convention
        pvar (default: 'ALTM'): string,  header for pressure variable used to identify missing data in METAR .csv file
        vari (default: 'Temperature'): string, header for variable to be plotted
        vectors (default: True): Boolean flag, if True, plots wind vectors in addition to variable
        spd (default: 'SPD'): string, header for windspeed variable in .csv file (only used if vectors = True)
        drct (default: 'DIR'): string, header for wind direction in .csv file (only used if vectors = True)
        lon (default: 'LON'): string, header for longitude data in .csv file
        lat (default: 'LAT'): string, header for latitude data in .csv file

        **kwargs : standard matplotlib keywork arguments accepted by "scatter"

    Example usage:
        ----------
        ax=plt.subplot(111,projection=ccrs.PlateCarree())
        AddMetars(ax,path='METARS/',stns='all',vari='TMP',vmin=-10, vmax=30,cmap='jet',vectors =False)

        Will scatter temperature observations on a map scaled from -10 to 30 (presumably in degrees Celsius)
        using the jet colormap.  Will not plot wind vectors.

    """
    import pandas as pd
    from datetime import datetime


    print("LOADING METAR DATA!")
    pdtime=pd.to_datetime(datetime.strptime(time,tfmt))
    stnfiles = glob.glob(path + '*.csv')
    if stns.lower() != 'all': #IF YOU'RE NOT LOOKING FOR ALL STATIONS
        if type(stns) is not list: #if stations is not a list (e.g., a single station)
            stns=[stns] #make it a list?

        ## Basically, keep file "i" if any part of the file "i" is in the stns list.
        stnfiles=[i for i in stnfiles if any(i.split('/')[-1] in j for j in stns)]

    lons,lats=[],[]
    uu,vv=[],[]
    outvar=[]
    for fds, fx in enumerate(stnfiles):
    # GET THE META-DATA AT THE TOP OF THE FILE
        with open(fx, 'r') as f:
            meta_name = f.readline()
            unit_name = f.readline()

        print("Reading Data For: %s" % meta_name)

        df = pd.read_csv(fx, header=2, parse_dates=['TIME'])
        df = df[np.isnan(df[pvar]) == False].reset_index()  # remove any values with bad altimiter

        df.set_index('TIME',inplace=True)


        idx=np.argmin(np.abs(df.index-pdtime))
        data=df.loc[df.index == df.index[idx]]

        if vectors == True:
            spd_val=data[spd]
            drct_val=data[drct]

            u,v=convert_met_2_uv(drct_val,spd_val)
            uu.append(u)
            vv.append(v)

        lons.append(data[lon][0])
        lats.append(data[lat][0])
        outvar.append(data[vari][0])

        df = None  ## release memory, just in case.

    lons,lats=np.array(lons).squeeze(),np.array(lats).squeeze()
    uu,vv=np.array(uu).squeeze(),np.array(vv).squeeze()
    outvar=np.array(outvar).squeeze()

    if vectors == True:
        ax.barbs(lons,lats,uu,vv,length=6.5,color='r',transform=ccrs.PlateCarree(),zorder=10)

    ax.scatter(lons,lats, c=outvar, transform=ccrs.PlateCarree(), **kwargs)



    return

        ## first set input time to pandas datetime

def AddMap(ax, extent=[0, 360, -90, 90], zorder=1, features=[],reso='50m'):
    """Function to add a map with geographic features to a figure."""
    ax.add_feature(cfeature.LAND, zorder=zorder + 1)
    ax.add_feature(cfeature.OCEAN, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale(reso), zorder=8)
    if 'borders' in features:
        ax.add_feature(cfeature.BORDERS.with_scale(reso), linestyle=':', zorder=8)
    if 'states' in features:
        ax.add_feature(cfeature.STATES.with_scale(reso), linestyle=':', zorder=8)
    if 'lakes' in features:
        ax.add_feature(cfeature.LAKES.with_scale(reso), alpha=1.0, zorder=zorder + 1)
    ax.set_extent(extent, ccrs.PlateCarree())


def plot_maxmin_points(ax, lon, lat, data, extrema, nsize, symbol, color='k',
                       plotValue=True, transform=ccrs.PlateCarree(), zorder=8):
    """
    This function will find and plot relative maximum and minimum for a 2D grid.
    This function is a modified reproduction of an example from
    Unidata's MetPy: https://unidata.github.io/MetPy/latest/index.html:
    Specifically this one:
    https://unidata.github.io/python-gallery/examples/HILO_Symbol_Plot.html#sphx-glr-examples-hilo-symbol-plot-py
    """
    from scipy.ndimage.filters import maximum_filter, minimum_filter

    if (extrema == 'max'):
        data_ext = maximum_filter(data, nsize, mode='nearest')
    elif (extrema == 'min'):
        data_ext = minimum_filter(data, nsize, mode='nearest')
    else:
        raise ValueError('Value for hilo must be either max or min')

    mxy, mxx = np.where(data_ext == data)

    for i in range(len(mxy)):
        ax.text(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]], symbol, color=color, size=24,
                clip_on=True, horizontalalignment='center', verticalalignment='center',
                transform=transform, zorder=5)
        ax.text(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]],
                '\n' + str(np.int(data[mxy[i], mxx[i]])),
                color=color, size=12, clip_on=True, fontweight='bold',
                horizontalalignment='center', verticalalignment='top', transform=transform, zorder=zorder)


def plot_mslp(ax, lons, lats, mslp, time, hsize=50, lsize=25, zorder=5, hl=True):
    """Plots Mean Sea Level Pressure and Finds max/min points."""
    # Plot MSLP
    kw_clabels = {'fontsize': 11, 'inline': True, 'inline_spacing': 5, 'fmt': '%i',
                  'rightside_up': True, 'use_clabeltext': True}
    clevmslp = np.arange(800., 1120., 4)
    cs2 = ax.contour(lons, lats, mslp, clevmslp, colors='k', linewidths=1.25,
                     linestyles='solid', transform=ccrs.PlateCarree(), zorder=5)
    plt.clabel(cs2, **kw_clabels)

    # Use definition to plot H/L symbols
    if hl == True:
        plot_maxmin_points(ax, lons, lats, mslp, 'max', hsize, symbol='H', color='b', transform=ccrs.PlateCarree(),
                           zorder=zorder)
        plot_maxmin_points(ax, lons, lats, mslp, 'min', lsize, symbol='L', color='r', transform=ccrs.PlateCarree(),
                           zorder=zorder)

    # Put on some titles
    title = 'MSLP (hPa) (contour)'
    return title


def smooth_field(data, sigma):
    """Super simple function to smooth a 2D field."""
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(data, sigma=sigma)


def acc_precip(namelist_dictionary, vdx=0, src_name='Stage 4', title='', time_indx=1):
    """Helper function to read in precipitation data from NARR, WRF or STAGE 4!
       This function handles all the common steps, such as defining the time boundaries, and handling the titles
    """
    from datetime import datetime, timedelta

    if type(namelist_dictionary['pcp_accum_times']) is not list or \
            len(namelist_dictionary['pcp_accum_times']) == 1:
        print("""Assuming %s precip accumulation is over single time
            step ending on:%s""" % (src_name, namelist_dictionary['pcp_accum_times']))

        if 'stage 4' in src_name.lower():
            intime = datetime.strptime(namelist_dictionary['pcp_accum_times'],
                                       namelist_dictionary['time_format'])
        else:  ## assume it's NARR
            intime = [datetime.strptime(namelist_dictionary['pcp_accum_times'],
                                        namelist_dictionary['time_format']) - timedelta(hours=3),
                      datetime.strptime(namelist_dictionary['pcp_accum_times'],
                                        namelist_dictionary['time_format'])]

        if 'stage 4' in src_name.lower():
            lons, lats, var_data, actual_st, endtime, units, step = \
                pull_stg4_data_single_time(
                    namelist_dictionary['stg4_folder'], intime,
                    namelist_dictionary['stg4_inc'])
        if vdx == 0:
            title = """%s Accumulated Precipitation %s \n starting at %s and ending %s""" % \
                    (src_name, units,actual_st.strftime(namelist_dictionary['time_format']),
                     endtime.strftime(namelist_dictionary['time_format']))
        else:
            title += """\n %s Accumulated Precipitation %s \n starting at %s and ending %s""" % \
                     (src_name, units,actual_st.strftime(namelist_dictionary['time_format']),
                      endtime.strftime(namelist_dictionary['time_format']))
    else:
        print("""Assuming stage 4 precip accumulation is accumulated between %s and %s """ % (
            namelist_dictionary['pcp_accum_times'][0],
            namelist_dictionary['pcp_accum_times'][1]))

        intime = [datetime.strptime(namelist_dictionary['pcp_accum_times'][0],
                                    namelist_dictionary['time_format']),
                  datetime.strptime(namelist_dictionary['pcp_accum_times'][1],
                                    namelist_dictionary['time_format'])]

        if 'stage 4' in src_name.lower():
            lons, lats, var_data, actual_st, endtime, units, step = \
                pull_stg4_data_accum_times(
                    namelist_dictionary['stg4_folder'], intime,
                    namelist_dictionary['stg4_inc'])

        else:  ## ASSUME IT's NARR
            lons, lats, var_data, actual_st, endtime, units = \
                read_NARR_pcp(namelist_dictionary['narr_folder'], intime)

        if vdx == 0:
            title = """%s Accumulated Precipitation %s \n starting at %s and ending %s""" % \
                    (src_name, units, actual_st.strftime(namelist_dictionary['time_format']),
                     endtime.strftime(namelist_dictionary['time_format']))
        else:
            title += """\n %s Accumulated Precipitation %s \n starting at %s and ending %s""" % \
                     (src_name, units, actual_st.strftime(namelist_dictionary['time_format']),
                      endtime.strftime(namelist_dictionary['time_format']))

        return lons, lats, var_data, title, units


def pull_stg4_data_single_time(stg4path, time, increment='06h'):
    """Function to load stage4 data for a single time into an array
    and returns for easy plotting"""
    import pygrib as pg
    from datetime import datetime, timedelta

    stg4name = 'Total Precipitation'
    files = sorted(glob.glob(stg4path + '*%s*' % increment))
    ## SORT FILES! ##
    if len(files) == 0:
        ## UTOH! NO FILES!
        print("NO FILES IN THIS DIRECTORY:")
        print("%s" % stg4path)
        print("CHECK YOUR DIRECTORY AND TRY AGAIN PLEASE!")

    for fdx, f in enumerate(files):
        ## LOOK FOR CORRECT TIME FIRST! ##
        pdata = pg.open(f)

        vari = pdata.select(name=stg4name)[0]
        starttime = vari.validDate
        endtime = starttime + timedelta(hours=int(vari.endStep))
        ## CHECK TIMES ##
        if time == endtime:
            # GREAT, YOU DID IT!
            print("FOUND TIME in STG4 DATA!")
            step = int(vari.endStep)
            units = vari.units
            lats, lons = vari.latlons()
            pcp = vari.values
            pdata.close()
            return lons, lats, pcp, starttime, endtime, units, step
        elif starttime < time < endtime:
            ## TIME IS WITHIN THE 2 TIME INCREMENTS!##
            print("""Hmm, your chosen time %s is between the starttime %s
            and endtime %s, I'll return the data, but the accumulated
            precpitation ending on %s may not be exactly what
            you're looking for""" % (time.strftime('%Y-%m-%d %H:%M:%S'),
                                     starttime.strftime('%Y-%m-%d %H:%M:%S'),
                                     endtime.strftime('%Y-%m-%d %H:%M:%S'),
                                     endtime.strftime('%Y-%m-%d %H:%M:%S')))

            step = int(vari.endStep)
            units = vari.units
            lats, lons = vari.latlons()
            pcp = vari.values
            pdata.close()
            return lons, lats, pcp, starttime, endtime, units, step

        else:
            pdata.close()

    print("""I COULD NOT FOR THE LIFE OF ME FIND A STAGE 4
        PRECIPITATION FILE WITH THE TIME YOU ARE LOOKING FOR!
        PLEASE CHECK YOUR TIMES AND YOUR FILES AND TRY AGAIN!""")

    return


def pull_stg4_data_accum_times(stg4path, times, increment='06h', return_all=False):
    """Function to load stage4 data for a single time into an array
    and returns for easy plotting"""
    import pygrib as pg
    from datetime import datetime, timedelta

    stg4name = 'Total Precipitation'
    files = sorted(glob.glob(stg4path + '*%s*' % increment))
    ## SORT FILES! ##
    if len(files) == 0:
        ## UTOH! NO FILES!
        print("NO FILES IN THIS DIRECTORY:")
        print("%s" % stg4path)
        print("CHECK YOUR DIRECTORY AND TRY AGAIN PLEASE!")

    if type(times) is not list:
        print("YOUR TIME INPUT IS NOT A LIST...")
        print("ARE YOU SURE YOU didn't want to call: 'pull_stg4_data_single_time'?")
        print("EXITING...")
        sys.exit(1)
    if len(times) == 1:
        print("YOUR TIME LIST ONLY HAS ONE VALUE!")
        print("ARE YOU SURE YOU didn't want to call: 'pull_stg4_data_single_time'?")
        print("EXITING...")
        sys.exit(1)

    ## UP START TIME BY INCREMENT TO ACCOUNT FOR FACT THAT FILE IS END-TIME ACC
    input_st = times[0] + timedelta(hours=int(increment[:2]))
    input_et = times[-1]

    start = False
    for fdx, f in enumerate(files):
        ## LOOK FOR CORRECT TIME FIRST! ##
        pdata = pg.open(f)

        vari = pdata.select(name=stg4name)[0]
        #starttime = vari.validDate
        starttime=datetime.strptime(f.split('.')[-2],'%Y%m%d%H')
        endtime = starttime
        ## CHECK TIMES FIRST LOOK FOR START TIME! ##
        if start == False:
            if input_st == starttime:
                ## PERFECT MATCH, GREAT! ##
                print("FOUND TIME in STG4 DATA! %s" % starttime.strftime('%Y-%m-%d-%H:%M:%S'))
                start = True  ## Turn on start flag
                step = int(vari.endStep)
                units = vari.units
                lats, lons = vari.latlons()
                pcp = vari.values
                if return_all == True:
                    ## Set up dummy axis for concatenation
                    pcp = pcp[None, :, :]
                pdata.close()
                actual_st = starttime - timedelta(hours=int(increment[:2]))
            elif starttime < input_st < endtime:
                ## NOT A PERFECT MATCH, but OKAY!
                print("""Hmm, your chosen starttime %s is between the starttime %s
                and endtime %s, we'll go with it, but may not be exactly what
                you're looking for""" % (input_st.strftime('%Y-%m-%d %H:%M:%S'),
                                         starttime.strftime('%Y-%m-%d %H:%M:%S'),
                                         endtime.strftime('%Y-%m-%d %H:%M:%S')))
                start = True
                step = int(vari.endStep)
                units = vari.units
                lats, lons = vari.latlons()
                pcp = vari.values
                if return_all == True:
                    ## Set up dummy axis for concatenation
                    pcp = pcp[None, :, :]
                pdata.close()
                actual_st = starttime - timedelta(hours=int(increment[:2]))
            else:
                ## KEEP LOOKING! ##
                pdata.close()
        else:  ## START IS TRUE, LOAD ALL THE DATA UNTIL YOU HIT THE END!
            if input_et == endtime:
                # GREAT, YOU DID IT!
                print("FOUND TIME in STG4 DATA! %s" % starttime.strftime('%Y-%m-%d-%H:%M:%S'))
                step = int(vari.endStep)
                units = vari.units
                lats, lons = vari.latlons()
                if return_all == True:
                    ## CONCATENATE TO ARRAY INSTEAD OF ADD!
                    pcp = np.concatenate((pcp, vari.values), axis=0)
                else:
                    # otherwise, just add it up.
                    pcp = pcp + vari.values
                pdata.close()
                return lons, lats, pcp, actual_st, endtime, units, step
            elif starttime < input_et < endtime:
                ## TIME IS WITHIN THE 2 TIME INCREMENTS!##
                print("""Hmm, your chosen time %s is between the starttime %s
                and endtime %s, I'll return the data, but the accumulated
                precpitation ending on %s may not be exactly what
                you're looking for""" % (input_et.strftime('%Y-%m-%d %H:%M:%S'),
                                         starttime.strftime('%Y-%m-%d %H:%M:%S'),
                                         endtime.strftime('%Y-%m-%d %H:%M:%S'),
                                         endtime.strftime('%Y-%m-%d %H:%M:%S')))

                step = int(vari.endStep)
                units = vari.units
                lats, lons = vari.latlons()
                if return_all == True:
                    ## CONCATENATE TO ARRAY INSTEAD OF ADD!
                    pcp = np.concatenate((pcp, vari.values), axis=0)
                else:
                    # otherwise, just add it up.
                    pcp = pcp + vari.values
                pdata.close()
                print(actual_st, endtime)
                return lons, lats, pcp, actual_st, endtime, units, step

            else:
                ## IF NOT AT END TIME, just add to it.
                if return_all == True:
                    ## CONCATENATE TO ARRAY INSTEAD OF ADD!
                    pcp = np.concatenate((pcp, vari.values), axis=0)
                else:
                    # otherwise, just add it up.
                    pcp = pcp + vari.values
                pdata.close()

    print("""I COULD NOT FOR THE LIFE OF ME FIND A STAGE 4
        PRECIPITATION FILE WITH THE TIME YOU ARE LOOKING FOR!
        PLEASE CHECK YOUR TIMES AND YOUR FILES AND TRY AGAIN!""")


def read_NARR_pcp(narrpath, intimes):
    """helper function to read in netcdf NARR precpitation data from:
       https://www.esrl.noaa.gov/psd/data/gridded/data.narr.monolevel.html.
       I cannot guarantee that NARR data from another source will work here...
    """
    from datetime import datetime, timedelta
    from netCDF4 import Dataset

    tfmt = '%Y-%m-%d %H:%M:%S'
    narr_start = datetime.strptime('1800-01-01 00:00:00', '%Y-%m-%d %H:%M:%S')

    intimes[0] = intimes[0] + timedelta(hours=3)  ## NEED TO ADD 3 HOURS  to account for time-ending adjustment.

    narr_files = sorted(glob.glob(narrpath + 'apcp*.nc'))
    print("Found %i NARR Precipitation Files... checking times for matches in dataset." % len(narr_files))
    init_data_array = True  ## flag indicated we are intializating the data array
    for ndx, n in enumerate(narr_files):
        cfile = Dataset(n, 'r')
        times = cfile.variables['time']
        tinit, tfin = narr_start + timedelta(hours=times[0]), narr_start + timedelta(hours=times[-1])

        if intimes[-1] < tinit:
            ## Final time choice is before FIRST time in NARR file, continue!
            print("File outside of time boundaries, skipping....")
            cfile.close()
            continue

        if intimes[0] > tfin:
            ## initial time choice is after LAST time in NARR file, continue!
            print("File outside of time boundaries, skipping....")
            cfile.close()
            continue

        ## okay, made it this far ... handle special cases.
        timelist = [narr_start + timedelta(hours=i) for i in times]
        if intimes[-1] > tfin:
            ## Special case one, final time is outside of time bounds.
            t1idx = np.argmin(np.abs(np.array(timelist) - intimes[0]))
            t2idx = -1
        elif intimes[0] < tinit:
            ## Special case two, intial time is outside of timebounds:
            t1idx = 0
            t2idx = np.argmin(np.abs(np.array(timelist) - intimes[-1]))
        else:
            ##both times should be in same file...
            t1idx = np.argmin(np.abs(np.array(timelist) - intimes[0]))
            t2idx = np.argmin(np.abs(np.array(timelist) - intimes[-1]))

        if init_data_array == True:  # intialize array!
            outlon = cfile.variables['lon'][:]
            outlon = np.ma.masked_less_equal(outlon, 0).filled(outlon + 360.)
            outlat = cfile.variables['lat'][:]
            units = cfile.variables['apcp'].units
            outpcp = np.nansum(cfile.variables['apcp'][t1idx:t2idx + 1, :], axis=0)
            outst = timelist[t1idx] - timedelta(
                hours=3)  ## Adjust back to title reflects actual accumulation increment.
            outet = timelist[t2idx]
            init_data_array = False
        else:
            outpcp = outpcp + np.nansum(cfile.variables['apcp'][t1idx:t2idx + 1, :], axis=0)
            outet = timelist[t2idx]
        cfile.close()

    print(outst, outet, intimes)
    return outlon, outlat, outpcp, outst, outet, units

def wrf_acc_pcp(namelist_dictionary, times, time_indx, ncdata,vdx=0,title=''):
    from datetime import datetime, timedelta
    ## GETTING WRF ACCUMULATED PRECIPTIATION FROM PRECIP RATE VAR!
    dt = (np.array(times)[1] - np.array(times)[0]).total_seconds()
    if type(namelist_dictionary['pcp_accum_times']) is not list or \
            len(namelist_dictionary['pcp_accum_times']) == 1:
        print("""Assuming stage 4 precip accumulation is over single time
                            WRF output step ending on:%s""" % namelist_dictionary['pcp_accum_times'])

        if time_indx == 0:
            print("Come on, you know you can't have any accumulated precipitation at the zero time, in fact,"
                  "Just to be safe, I'm going to exit out of the script here so you can take a good look at the namelist"
                  "and make sure you have it configured correctly.")
            sys.exit()
        ti1 = time_indx - 1
        ti2 = time_indx

    else:
        wrf1 = datetime.strptime(namelist_dictionary['pcp_accum_times'][0], namelist_dictionary['time_format'])
        wrf2 = datetime.strptime(namelist_dictionary['pcp_accum_times'][1], namelist_dictionary['time_format'])
        ti1 = np.argmin(np.abs(np.array(times) - wrf1))
        ti2 = np.argmin(np.abs(np.array(times) - wrf2))

    ## INCLUDE BOTH SNOW AND RAIN!
    wrf_total_pcp = np.nansum(ncdata.variables['Gridscale Precipitation rate'][ti1:ti2 + 1, :] * dt, axis=0)
    if 'Convective Precipitation rate' in ncdata.variables:
        ## ADD CONVECTIVE PRECIP
        wrf_total_pcp = wrf_total_pcp + \
                        np.nansum(ncdata.variables['Convective Precipitation rate'][ti1:ti2 + 1, :] * dt, axis=0)

    ## NOW CONVERT RATE TO ACCUMULATED...
    units = 'kg m^-2'

    if vdx == 0:
        title = """WRF Accumulated Precipitation [%s] \n starting at %s and ending %s""" % \
                (units, times[ti1].strftime(namelist_dictionary['time_format']),
                 times[ti2].strftime(namelist_dictionary['time_format']))
    else:
        title += """\n WRF Accumulated Precipitation [%s] \n starting at %s and ending %s""" % \
                (units, times[ti1].strftime(namelist_dictionary['time_format']),
                 times[ti2].strftime(namelist_dictionary['time_format']))


    return wrf_total_pcp,title,units



def wrf_acc_snow(namelist_dictionary, times, time_indx, ncdata,vdx=0,title=''):
    """Helper function to accumulate WRF snowfall between two times given a precipitation rate."""
    from datetime import datetime, timedelta
    ## GETTING WRF ACCUMULATED PRECIPTIATION FROM PRECIP RATE VAR!
    dt = (np.array(times)[1] - np.array(times)[0]).total_seconds()
    if type(namelist_dictionary['pcp_accum_times']) is not list or \
            len(namelist_dictionary['pcp_accum_times']) == 1:
        print("""Assuming stage 4 precip accumulation is over single time
                            WRF output step ending on:%s""" % namelist_dictionary['pcp_accum_times'])

        if time_indx == 0:
            print("Come on, you know you can't have any accumulated precipitation at the zero time, in fact,"
                  "Just to be safe, I'm going to exit out of the script here so you can take a good look at the namelist"
                  "and make sure you have it configured correctly.")
            sys.exit()
        ti1 = time_indx - 1
        ti2 = time_indx

    else:
        wrf1 = datetime.strptime(namelist_dictionary['pcp_accum_times'][0], namelist_dictionary['time_format'])
        wrf2 = datetime.strptime(namelist_dictionary['pcp_accum_times'][1], namelist_dictionary['time_format'])
        ti1 = np.argmin(np.abs(np.array(times) - wrf1))
        ti2 = np.argmin(np.abs(np.array(times) - wrf2))

    ## INCLUDE BOTH SNOW AND RAIN!
    wrf_total_pcp = np.nansum(ncdata.variables['Gridscale Snowfall rate'][ti1:ti2 + 1, :] * dt, axis=0)

    ## NOW CONVERT RATE TO ACCUMULATED...
    units = 'kg m^-2'

    if vdx == 0:
        title = """WRF Accumulated Snow [%s] \n starting at %s and ending %s""" % \
                (units, times[ti1].strftime(namelist_dictionary['time_format']),
                 times[ti2].strftime(namelist_dictionary['time_format']))
    else:
        title += """\n WRF Accumulated Snow [%s] \n starting at %s and ending %s""" % \
                (units, times[ti1].strftime(namelist_dictionary['time_format']),
                 times[ti2].strftime(namelist_dictionary['time_format']))


    return wrf_total_pcp,title,units



def convert_met_2_uv(drct,spd):
    """Converts a meteorological windspeed and direction into u/v vector components.
       Units are equal to whatever the speed units read in are.
    """
    drct = drct.values*np.pi/180.
    u = -spd.values * np.sin(drct)
    v = -spd.values * np.cos(drct)
    return u, v

def find_nearest_latlon(lon_array,lat_array,lon_point,lat_point):
    """Finds nearest grid-cell (indices) from WRF lat/lon outputs to an input lat/lon coordinate pair."""
    latmin=np.abs(lat_array-lat_point)
    lonmin=np.abs(lon_array-lon_point)
    c=np.maximum(lonmin,latmin)
    x, y = np.where(c == np.min(c))
    return x,y