"""
ERDC NWP COMMON FUNCTIONS

This file contains a number of different functions that are useful both for the
main post processing and visualization classes contained within this package, and
as stand alone functions for the end user.
"""
import numpy as np

# Define All Global Constants as a dictionary with values, and a
# second dictionary with descriptions and units of each constant.
constants = {
    'Rd': 287.05, 'Rv': 461.51,
    'Cp': 1004.6, 'Cw': 4218.0,
    'epsilon': 0.622,
    'rho_w': 1000., 'rho_i': 917.0,
    'g': 9.81,
    'Lv': 2.5E6, 'Lf': 3.34E5,
    'Sb': 5.67E-8,
    'Pref': 1000.
}

constant_description = {
    'Rd': ' Dry Gas Constant, J kg^{-1} K^{-1}', 'Rv': ' Moist Gas Constant, J kg^{-1} K^{-1}',
    'Cp': 'Specific Heat of Dry Air J kg^{-1} K^{-1}', 'Cw': 'Specific Heat of Water J kg^{-1} K^{-1}',
    'epsilon': 'Rd/Rv --> Ratio between the dry air and moist air gas constants',
    'rho_w': 'Density of Water kg m^{-3}', 'rho_i': 'Density of Ice kg m^{-3}',
    'g': 'Gravitational constant m s{^-2}',
    'Lv': 'Latent Heat of Vaporization J kg^{-1}', 'Lf': 'Latent Heat of Fusion J kg^{-1}',
    'Sb': 'Stefan-Boltzmann constant W m^{-2} K^{-4}',
    'Pref': 'Reference Pressure for Potential Temperature (hPa)'
}

#   Common Meteorological Functions
# -----------------------------------
#    - Potential Temperature
#    - Air Density
#    - Virtual Temperature
#    - Vapor Pressure Approximation
#    - Dew point temperature (given vapor pressure)
#    - Wet Bulb temperature approximation
#    - Dew point temperature (given mixing ratio)
# -----------------------------------


def print_constants():
    """This is a simple function to print out available constants along with a brief description and units

        To use the constants dictionary: simply load the constants: from functions.common import constants

        Rd= constants['Rd'] #Dry gas constant for example.
    """
    print("-------------------------")
    print("Printing Constants")
    print("-------------------------")
    for c in constants:
        print("%s = %.2f | %s"%(c,constants[c],constant_description[c]))
    print("-------------------------")

def Theta2TmpK(theta, P, pref=constants['Pref']):
    """Converts potential temperature (theta) at a given pressure to temperature

    Inputs:
        theta = potential temperature (K)
        P = pressure (hPa)
        pref (optional) = reference level pressure (hPa)

    Returns:
        temperature (K)
    """

    return theta * (P / pref) ** (constants['Rd'] / constants['Cp'])


def TmpK2Theta(T, P, pref=constants['Pref']):
    """Converts temperature at a given pressure to potential temperature (theta)

    Inputs:
        T = temperature  (K)
        P = pressure (hPa)
        Pref (optional) = reference level pressure (hPa)

    Returns:
         Theta (K)
    """

    return T * (pref / P) ** (constants['Rd'] / constants['Cp'])


def VirtualTemp(T, Q):
    """Computes Virtual Temperature (Tv)

    Inputs:
        T: temperature (K)
        Q: water vapor mass mixing Ratio (kg/kg)

    Returns:
        Tv: Virtual temperature (K)
    """
    return T * (1.0 + 0.6 * Q)


def AirDensity(T, P, Q, moist=True):
    """Computes the density of moist air

    Inputs:
        T = temperature (K)
        P = pressure (hPa)
        Q = water vapor mass mixing ratio (kg/kg)
        moist (optional) = flag to compute moist air density (default = True)

    Returns:
        rho_air (kg/m^3)
    """
    P = P * 100.  # Convert Pressure in hPa to pressure in Pa
    # Boolean flag for adjusting the air temperature for moisture content (recommended)
    if moist:
        Tv = VirtualTemp(T, Q)
    else:
        Tv = T
    return P / (constants['Rd'] * Tv)


def VaporPressure(T, phase="liquid"):
    """Compute Water vapor pressure over liquid water or ice.
       Note: This function is ported over and modified from Thomas Chubb's "Skew-T.py package.

    Inputs:
        T = dew point temperature (K) or (for saturation vapor pressure) temperature (K)
        phase (optional): ['liquid'],'ice'. If 'liquid', do simple dew point. If 'ice',
        return saturation vapour pressure as follows:

            Tc>=0: es = es_liquid
            Tc <0: es = es_ice

    Returns:
         vapor pressure (VP) (Pa)

    SOURCE: http://cires.colorado.edu/~voemel/vp.html (#2:
    CIMO guide (WMO 2008), modified to return values in Pa)

    This formulation is chosen because of its appealing simplicity,
    but it performs very well with respect to the reference forms
    at temperatures above -40 C.
    """

    Tc = T - 273.15  # convert temperature in K to temperature in C

    over_liquid = 6.112 * np.exp(17.67 * Tc / (Tc + 243.12)) * 100.
    over_ice = 6.112 * np.exp(22.46 * Tc / (Tc + 272.62)) * 100.

    if phase.lower() not in ['liquid', 'ice']:
        print("I don't know what you're trying to do with the phase")
        print("I'm setting the default phase to liquid")
        phase = 'liquid'

    # if/then conditional function to return vapor pressure with respect to ice
    # where the temperature array is < 0 C.
    if phase == "liquid":
        return over_liquid
    elif phase == "ice":
        return np.ma.where(Tc < 0, over_ice, over_liquid)


def DewPoint(VP):
    """ Computes (approximates) dew point temperature
        using Bolton's (1980, MWR, p1047) formula.

    Inputs:
        VP = vapor pressure (Pa)

    Returns:
        dew point temperature (Td) (K)
    """

    ln_ratio = np.log(VP / 611.2)
    Td = ((17.67 - ln_ratio) * 273.15 + 243.5 * ln_ratio) / (17.67 - ln_ratio)
    return Td


def WetBulb(T, RH):
    """Computes (approximates) Wet-Bulb Temperature from Stull (2011).
       Note: This function is ported over and modified from Thomas Chubb's "Skew-T.py package.

    Inputs:
        T = temperature (K)
        RH = Relative Humidity (%)
    Returns:
        Twb (K)
    """
    Tc = T - 273.  # Convert Temperature from Kelvin to degrees Celcius

    Tw = Tc * np.arctan(0.151977 * (RH + 8.313659) ** 0.5) + \
         np.arctan(Tc + RH) - np.arctan(RH - 1.676331) + \
         0.00391838 * RH ** 1.5 * np.arctan(0.023101 * RH) - \
         4.686035

    Tw = Tw + 273.15  # Convert Temperature from degrees Celcius to Kelvin

    return Tw

def MixR2VaporPress(Q, P):
    """Computes Vapor Pressure From water vapor mixing ratio and atmospheric pressure.
       Note: This function is ported over and modified from Thomas Chubb's "Skew-T.py package.

    Inputs:
        Q = water vapor mass mixing ratio (kg/kg)
        P = atmospheric pressure (hPa)

    returns
        vapor pressure (VP) (Pa)
    """
    P = P * 100.  # Convert pressure from hPa to Pa

    return Q * P / (constants['epsilon'] + Q)


#  Supplementary functions used in WRF post processing
# -----------------------------------
#    - Mean Sea Level Pressure approximation
#    - pressure level interpolation
#    - Visibility from hydrometeor mass mixing ratios.
# -----------------------------------

def mslp(Ps, T, HGT, Q):
    """
        This function approximates the mean sea level pressure reduction using the simple method
        that assumes an isothermal temperature profile "through the ground"
        Note, this probably has pretty significant inaccuracies for elevations above 2000 meters.

        mslp = Ps * e^(g*HGT/(Rd*T))

        Inputs:
            T = temperature (K)
            Ps = surface pressure (hPa)
            HGT = surface elevation (terrain height) (m)
            Q = water vapor mass mixing ratio (kg/kg)

        Outputs:
            MSLP = Mean Sea Level Pressure (Pa)
    """

    TV = VirtualTemp(T, Q)
    Ps = Ps * 100.  # Convert pressure in hPa to pressure in Pa

    MSLP = Ps * np.exp((constants['g'] * HGT) / (287.05 * TV))

    return MSLP

def p_interp(vari, P, zind=1, lev=850, masked=True):
    """This Function interpolates a variable to a constant pressure surface.
       The Masked option indicates to mask out where the pressure surface is below the ground,
       i.e. level > p_sfc.

       Note: This function uses the Numpy "Interp" function which only works on 1D arrays.
       Therefore, it typically runs SLOWLY over large data-sets, since it has to manually loop through all
       x,y,t grid points.

       Future versions will include a "multiprocessing" option.

       Inputs:
            vari = variable to interpolate (e.g., temperature) [Must vary in z]
            P = pressure (Pa) [Must vary in z]
            lev (default = 850) = pressure level (or list of pressure levels to interpolate to) (hPa)
            zind (default = 1) = array index of height coordinate (remember Python indexing starts at zero)
            masked (default = True) = flag to mask the variable where the variable is masked at
                locations where the pressure level is greater then the maximum pressure
                (i.e., it's below the ground)

        Returns:
    """

    P = P / 100.  # Convert pressure to hPa to match the pressure levels

    if type(lev) is not list:
        # If level is not a list (i.e, it's a simple number), convert it to a list.
        lev = [lev]

    # Sort list incase some idiot fed in a list that was all out of order ... Ted.
    lev = sorted(lev)[::-1]

    # This is an important step where the variable is "rolled" such that the z-index is in the first column.
    vari = np.rollaxis(vari, zind)
    P = np.rollaxis(P, zind)

    # Set up the end results to be shape (len(lev), other dimensions)
    shape = list(np.shape(vari))
    shape[0] = len(lev)
    result = np.ma.empty(shape)

    # Enumerate over ALL dimensions that are not the 1st (vertical dimension)
    for idx, x in np.ndenumerate(vari[0, :]):
        index = tuple([slice(None, None, None)] + list(idx))
        i0 = tuple([0] + list(idx))
        # Get pressure and Variable at current index.
        Pnow = np.log10(P[index])
        Vnow = vari[index]
        # Organize such that pressure increases (requirement for numpy)
        # re-align the variable to the standard decreasing pressure levels [::-1]
        r = np.interp(np.log10(lev[::-1]), Pnow[::-1], Vnow[::-1])[::-1]
        if masked:
            # if masked, mask out anywhere the pressure level is greater than the max pressure.
            r = np.ma.masked_where(np.array(lev) > P[i0], r)
        result[index] = r

    # return the rolled variable back to to it's original shape by moving the zero axis to the input z-index.
    return np.rollaxis(result, 0, zind + 1)


def z_interp(vari, Z, zind=1, lev=1000):
    """This Function interpolates a variable to a constant height surface.

       Note: This function uses the Numpy "Interp" function which only works on 1D arrays.
       Therefore, it typically runs SLOWLY over large data-sets, since it has to manually loop through all
       x,y,t grid points.

       Future versions will include a "multiprocessing" option.

       Inputs:
            vari = variable to interpolate (e.g., temperature) [Must vary in z]
            Z = Height (meters) [Must vary in z]
            lev (default = 1000) = height level (or list of height levels to interpolate to) (m)
            zind (default = 1) = array index of height coordinate (remember Python indexing starts at zero"

        Returns: height interpolated variable

    """

    if type(lev) is not list:
        # If level is not a list (i.e, it's a simple number), convert it to a list.
        lev = [lev]

    # Sort list incase some idiot fed in a list that was all out of order ... Ted.
    lev = sorted(lev)

    # This is an important step where the variable is "rolled" such that the z-index is in the first column.
    vari = np.rollaxis(vari, zind)
    Z = np.rollaxis(Z, zind)

    # Set up the end results to be shape (len(lev), other dimensions)
    shape = list(np.shape(vari))
    shape[0] = len(lev)
    result = np.ma.empty(shape)

    # Enumerate over ALL dimensions that are not the 1st (vertical dimension)
    for idx, x in np.ndenumerate(vari[0, :]):
        index = tuple([slice(None, None, None)] + list(idx))
        i0 = tuple([0] + list(idx))
        # Get pressure and Variable at current index.
        Znow = Z[index]
        Vnow = vari[index]
        r = np.interp(lev, Znow, Vnow)

        result[index] = r

    # return the rolled variable back to to it's original shape by moving the zero axis to the input z-index.
    return np.rollaxis(result, 0, zind + 1)

def get_volume(Qcld, Qvapor, Qice, Qsnow, Qrain, rho_a):
    """Compute volume of the combined air / hydrometeor medium"""
    vol = (1. + Qvapor) / rho_a + (Qcld + Qrain) / constants['rho_w'] + (Qice + Qsnow) / constants['rho_i']
    return vol


def get_conc(hydrometeor, volume):
    """Converts hydrometeor mass mixing ratio (kg kg^1) to mass concentration (g m^-3) from volume."""
    conc = hydrometeor / volume
    return np.ma.masked_less(conc * 1000., 0).filled(0.0)


def get_vis(data_dict, volume, bsnow=False):
    """This function
       estimates visibility from a known extinction formula for different hydrometeor species
       and concentrations (g m^-3)  Returns a visibility based on the 2% visual contrast threshold
       (Koschmeider Forumla) in meters.
       This function is ported over from the Unified Post Processor Code:
            https://dtcenter.org/community-code/unified-post-processor-upp

       Based on the formulas detailed in Warner and Stoelinga 1999:
            https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281999%29038%3C0385%3ANMSMSO%3E2.0.CO%3B2

       Inputs:
            data_dict = dictionary of gridded data with hydrometeror mass mixing ratios of: Rain, Snow,
                Cloud Water, and Cloud Ice (kg/kg).

            volume = combined air/hydrometeor volume (m^3)
            bsnow (default = False) = Either false, or array of blowing snow mixing ratio (kg/kg)
                Added for a specific project, and most likely not applicable to you.

       Returns:
            visibility (km) --> automatically masked where greater than 24km.
       """

    # This block defines all functions required to get visibility from WRF data following the UPP method.

    Extinction_dictionary = {'QCLOUD': [144.7, 0.8800],
                             'QICE': [327.8, 1.0],
                             'QRAIN': [2.24, 0.75],
                             'QSNOW': [10.36, 0.7776]}  # First number = coefficient, second = exponent

    ext_coeff = 1.0E-5
    for i in Extinction_dictionary:
        conc = get_conc(data_dict[i], volume)
        ext_coeff += (Extinction_dictionary[i][0] * conc ** Extinction_dictionary[i][1])

    if bsnow:
        ext_coeff += bsnow / 1000.
    vis = np.ma.masked_greater(3.912 / ext_coeff, 24).filled(24.0)
    return vis


def wrf_visibility(rho_a, QCLOUD, QVAPOR, QICE, QRAIN, QSNOW):
    """This is the parent visibility function wrapper that defines the data-dictionary and combined air/hydrometeor
           volume which is passed to the "get_vis" function.  See "get_vis", "get_volume" and "get_conc" for more
           specific information on the computation of visibility.

        This function is ported over from the Unified Post Processor Code:
            https://dtcenter.org/community-code/unified-post-processor-upp

         Based on the formulas detailed in Warner and Stoelinga 1999:
            https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281999%29038%3C0385%3ANMSMSO%3E2.0.CO%3B2

        Inputs:
            rho_a =  air density (kg m^-3)
            QCLOUD = cloud water mass mixing ratio (kg/kg)
            QICE   = cloud ice mass mixing ratio (kg/kg)
            QSNOW  = snow water mass mixing ratio (kg/kg)
            QRAIN  = rain water mass mixing ratio (kg/kg)
            QVAPOR = water vapor mass mixing ratio (kg/kg)

            Returns:
            Visibility (km)
        """

    # Define dictionary to hold hydrometeor arrays
    Data_Dictionary = {'QCLOUD': np.squeeze(QCLOUD), 'QVAPOR': np.squeeze(QVAPOR), 'QICE': np.squeeze(QICE),
                       'QRAIN': np.squeeze(QRAIN), 'QSNOW': np.squeeze(QSNOW)}

    # get the combine air / hydrometor volume
    volume = get_volume(Data_Dictionary['QCLOUD'], Data_Dictionary['QVAPOR'],
                        Data_Dictionary['QICE'], Data_Dictionary['QSNOW'], Data_Dictionary['QRAIN'],
                        rho_a)

    # Compute visibility
    vis = get_vis(Data_Dictionary, volume)
    return vis
