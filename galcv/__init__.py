# This script contains all of the code necessary for galcv.

# We want a function that returns the cosmic variance values for input apparent magnitudes

import numpy as np
import pandas as pd
import os
from scipy import integrate
from scipy import interpolate

# A dictionary of the parameters for which I have the exact interp files
fitParams = dict(mag=np.array([22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34]), z=np.array([5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15.]), zW=np.array([0.1, 0.15, 0.25, 0.5, 0.75, 1., 1.5, 2.]))

# This is a tuple of int and float types to check inputs against
intOrFloat = (int, np.int8, np.int16, np.int32, np.int64, float, np.float16, np.float32, np.float64)

# How many decimal places to round the outputs to
roundTo = 4


def getcv(mag, area, z, zW=1., appOrAbs='apparent', CMF_method='nu-scaling', interpWarning=1, goFast='No'):
    '''
    This function returns relative cosmic variance results. This function is a wrapper function for formatting. The actual calculation happens in singlecv()

    Parameters
    -------------------------
    mag : int, float, list, or numpy.ndarray
        The magnitude(s) to consider. This must be in rest-UV (1500 - 2800 Angstroms) AB magnitude (default apparent magnitude; see appOrAbs param)
    area : int or float
        The area of a survey in arcmin^2 (square survey pattern only)
    z : int or float
        The central redshift of the survey
    zW : int or float
        The width of the redshift bin the survey is considering. Default is 1.
    appOrAbs: 'apparent' or 'absolute'
        Whether the mag input(s) are in apparent magnitudes or in absolute magnitudes
    CMF_method: 'nu-scaling' or 'PS-scaling'
        The method used for generating the conditional mass function. See Trapp & Furlanetto (2020) for details.
    interpWarning: int or float
        Flag for displaying interpolation warning message. 0 for no message, 1 for short message (Default), 2 for long message

    Returns
    -------------------------
    A Python list of cosmic variance values of the same length as the mag input
    '''

    # My code uses apparent magnitudes, so it they are given in absolute magnitudes, convert to apparent first
    if appOrAbs == 'absolute':
        mag = absToApp(Mabs=mag, z=z)

    # Check to make sure the keywords have the correct formats
    if goFast == 'No':
        checkVars(mag=mag, area=area, z=z, zW=zW, appOrAbs=appOrAbs, CMF_method=CMF_method, interpWarning=interpWarning)

    # Now, if mag is just an int or float, return an int or float
    if isinstance(mag, intOrFloat):
        return singlecv(mag=mag, area=area, z=z, zW=zW, CMF_method=CMF_method, interpWarning=interpWarning, goFast=goFast)

    else:  # else, return a list of the cv values
        answer = [singlecv(mag=a_mag, area=area, z=z, zW=zW, CMF_method=CMF_method, interpWarning=interpWarning, goFast=goFast) for a_mag in mag]
        if any([a_answer == np.nan for a_answer in answer]):
            if goFast == 'No':
                print('\nSome mag values are too bright to estimate cosmic variance. Those values are returned as np.nan objects.')
        return answer

    # If we want other variables able to be passed as an array, we need a nest of
    # if statements for the ways to properly call singlecv()


def singlecv(mag, area, z, zW, CMF_method, interpWarning, goFast):
    '''
    This function returns relative cosmic variance results by reading in interp files of cosmic variance for parameters near those given, and interpolating between them.

    Parameters
    -------------------------
    mag : int or float
        The magnitudes to consider. The type of magnitude is determined by the magType argument
    area : int or float
        The area of a survey in arcmin^2
    z : int or float
        The central redshift of the survey
    zW : int or float
        The width of the redshift bin the survey is considering. Default is 1
    interpWarning: bolean
        Flag for displaying interpolation warning message
    '''

    kwargs_fitMatches = dict(mag=mag, z=z, zW=zW)

    # For each parameter, check to see if the provided one matches one of the fit files exactly
    isExact, interpBetween = fitMatches(kwargs_fitMatches, fitParams)

    # isExact is a bolean that tells us if the parameters ALL have exact fit files for them
    # if so, this function becomes very simple.
    if isExact:
        # Read in the fit files and output the answer
        thecv = readcv(mag=mag, area=area, z=z, zW=zW, CMF_method=CMF_method)
        return round(thecv, roundTo)
        ##############
    else:  # Gotta do some interpolation
        if interpWarning == 1:
            print('Parameter combination requires interpolation.')
        elif interpWarning == 2:
            print('This code does not provide results for the specific combination of \'z\' and \'zW\' parameters. In this case, we interpolate between nearest valid parameters, which are:')
            print('mag: ', fitParams['mag'])
            print('z: ', fitParams['z'])
            print('zW: ', fitParams['zW'])
            print('We interpolate with log10(cosmic variance) between parameters, as cosmic variance is close to a powerlaw with these parameters.')
            print('To shorten this message, set interpWarning = 1. To stop this message, set interpWarning = 0')

        # interpBetween has all of the combination of parameters we need.
        # Create empty arrays to hold the
        mags = np.zeros(2)
        zs = np.zeros((2, 2))
        zWs = np.zeros((2, 2, 2))
        zW_epcvs = np.zeros((2, 2, 2))

        # Read in all of the fit files adjacent to the desired parameters.
        # If one or more of the parameters are an exact match for a fit file, read it in twice
        for i in range(2):
            i_mag = interpBetween['mag'][i]
            mags[i] = i_mag
            for j in range(2):
                j_z = interpBetween['z'][j]
                zs[i][j] = j_z
                for k in range(2):
                    k_zW = interpBetween['zW'][k]

                    zWs[i][j][k] = k_zW
                    zW_epcvs[i][j][k] = readcv(mag=i_mag, area=area, z=j_z, zW=k_zW, CMF_method=CMF_method, verbose=False)
        if np.any(np.isnan(zW_epcvs)):  # If any of the epcv values required are np.nan, can't do the interpolation
            if goFast == 'No':
                print('Apparent magnitude {:.2f} is too bright for a cosmic variance interpolation estimate at this area, z, and zW'.format(mag))
            return np.nan

        # Now we interpolate between constant mag and z, but differing zW epcvs
        z_epcvs = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                j_epcvs = zW_epcvs[i][j]
                j_zWs = zWs[i][j]
                z_epcvs[i][j] = interpcv(zW, ixs=j_zWs, iys=j_epcvs)

        # now interpolate between constant mag but differing values of z
        mag_epcvs = np.zeros(2)
        for i in range(2):
            i_epcvs = z_epcvs[i]
            i_zs = zs[i]
            mag_epcvs[i] = interpcv(z, ixs=i_zs, iys=i_epcvs)

        # Finally, interpolate between magnitudes
        final_epcvs = interpcv(mag, ixs=mags, iys=mag_epcvs)
        return round(final_epcvs, roundTo)


def fitMatches(kwargs, fitParams):
    '''
    This function takes two dictionaries, one is the parameters that are being called, the other is the parameters for which I have exact fit files. If I have exact fit files for all parameters, then isExact is returned True, and interpBetween is returned None. Otherwise, isExact is returned False and interpBetween is a dict of the parameters of fit files that I should use to interp between

    Parameters
    -------------------------
    kwargs : dict
        A dictionary of all the parameters passed to singlecv
    fitParams : dict
        A dictionary of all the parameters for which I have exact fit files

    '''
    mag = kwargs['mag']
    mags = fitParams['mag']
    z = kwargs['z']
    zs = fitParams['z']
    zW = kwargs['zW']
    zWs = fitParams['zW']

    # If all of the parameters have exact fit files, return True and None
    if (mag in mags) and (z in zs) and (zW in zWs):
        return (True, None)
    else:
        # Check to see if mag is one of the ones I have an fit file for.
        # If so, return itself twice for the singlecv function to use
        if mag in mags:
            interpBetween = dict(mag=[mag, mag])
        else:  # Otherwise, append the closest values on either side
            maglo = mags[mags < mag].max()
            maghi = mags[mags > mag].min()
            interpBetween = dict(mag=[maglo, maghi])
        # Now the same for the z parameters
        if z in zs:
            interpBetween['z'] = [z, z]
        else:
            zlo = zs[zs < z].max()
            zhi = zs[zs > z].min()
            interpBetween['z'] = [zlo, zhi]
        # Finally the zW param
        if zW in zWs:
            interpBetween['zW'] = [zW, zW]
        else:
            zWlo = zWs[zWs < zW].max()
            zWhi = zWs[zWs > zW].min()
            interpBetween['zW'] = [zWlo, zWhi]

        return (False, interpBetween)


def readcv(mag, area, z, zW, CMF_method, verbose=True):
    '''
    This function reads in the cv fits and returns the cv value, or np.nan if there is no estimate

    Parameters
    -------------------------
    mag : int or float
        The apparent magnitude to consider
    area : int or float
        The survey area in arcmin^2
    z : int or float
        The redshift
    zW : int or float
        The redshift bin width
    '''

    if mag not in fitParams['mag']:
        raise Exception('Error, readcv should only be recieving magnitudes I know are in the fit files.')

    # Which CMF method file to read in
    if CMF_method == 'nu-scaling':
        CMF_string = 'stdScale'
    elif CMF_method == 'PS-scaling':
        CMF_string = 'PSfrac'

    # Read in the file of fits
    THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
    my_file = os.path.join(THIS_FOLDER, 'varepsilonCV_' + CMF_string + '_z{:.0f}_zW{:.2f}.pkl'.format(z, zW))
    dffit = pd.read_pickle(my_file)

    mapps = dffit['mapp'].values
    Psis = dffit['Psi'].values
    gammas = dffit['gamma'].values
    bs = dffit['b'].values
    minAs = dffit['minA'].values
    maxAs = dffit['maxA'].values

    if mag not in mapps:
        raise Exception('Error, the mag supplied to readcv is not in the mapps array from the fit file, but it should be.')

    # What index is the mag in the mapps array?
    whichIndex = np.where(mapps == mag)[0][0]
    # Check to see if there is an estimate at this survey area,
    # and also make sure there is an estimate at this magnitude at all
    if (area >= minAs[whichIndex]) and ~np.isnan(minAs[whichIndex]):
        return 10**log10Eps(area, Psis[whichIndex], gammas[whichIndex], bs[whichIndex])
    else:
        if verbose:
            print('apparent magnitude {:.2f} is too bright to calculate cosmic variance.'.format(mag))
            i = 1
            while (mag + i) < max(mapps):  # Finding the brightest mag that has an estimate at this area
                checkIdx = np.where(mapps == (mag + i))[0][0]
                if (area >= minAs[checkIdx]) and minAs[checkIdx] != -9999:
                    print('minimum magnitude for this area, z, and zW: mag = ', (mag + i))
                    break
                i += 1
            print('Returning np.nan object.')
        return np.nan


def interpcv(x, ixs, iys):
    '''
    This function takes 2 x values (ixs) and 2 y values (iys), interpolates a line between them (in log10(y)), and then evaluates x on that line

    Parameters
    -------------------------
    x : int or float
        The value to evaluate at
    ixs : int or float
        The x values for the interpolation
    iys : int or float
        The epcv values
    '''

    if ixs[0] > ixs[1]:
        raise Exception('The two ixs put into interpcv should be sorted low to high')
    if x < ixs[0] or x > ixs[1]:
        raise Exception('The x value given to interpcv should be between the two ixs values')

    iys = np.log10(iys)
    rise = iys[1] - iys[0]
    run = ixs[1] - ixs[0]
    # If the same x and y values are given, just return one of the y values.
    if (rise == 0) or (run == 0):
        return 10**iys[1]
    else:
        slope = rise / run
        ans = iys[0] + slope * (x - ixs[0])
        return 10**ans


def inRange(value, theRange):
    '''
    This function checks if the value is within the range theRange

    Parameters
    -------------------------
    value : int or float
        The numerical value
    theRange: list
        The lower and upper bounds of the range
    '''
    if not isinstance(value, intOrFloat):
        raise Exception('The argument of inRange() must be a single int or float')

    if type(theRange) != list:
        raise Exception('theRange must be a list with two values')
    elif len(theRange) != 2:
        raise Exception('theRange must be a list with two values')

    if (value >= theRange[0]) and (value <= theRange[1]):
        return True
    else:
        return False


def checkVars(mag, area, z, zW, appOrAbs, CMF_method, interpWarning):
    # Check if the variables have the correct dtype and are in the correct ranges. If not, raise an exception.
    magRange = [min(fitParams['mag']), max(fitParams['mag'])]
    areaRange = [1, 3.16e4]
    zRange = [min(fitParams['z']), max(fitParams['z'])]
    zWRange = [min(fitParams['zW']), max(fitParams['zW'])]

    # Check the mag variable
    if type(mag) in (list, np.ndarray):
        if any([not isinstance(a_mag, intOrFloat) for a_mag in mag]):
            raise Exception('All values in the list/array \'mag\' must be a int or float')
        elif any([not inRange(a_mag, magRange) for a_mag in mag]):
            raise Exception('at least one mag value is outside of apparent mag range: {}'.format(magRange))

    elif isinstance(mag, intOrFloat):
        if not inRange(mag, magRange):
            raise Exception('mag value is outside of apparent mag range: {}'.format(magRange))
    else:
        raise Exception('\'mag\' must be a float, int, list, or numpy array')

    # Now the area variable
    if isinstance(area, intOrFloat):
        if not inRange(area, areaRange):
            raise Exception('area value outside of area range: {}'.format(areaRange))
    else:
        raise Exception('area must be a float or an int')

    # Now the z variable
    if isinstance(z, intOrFloat):
        if not inRange(z, zRange):
            raise Exception('z value outside of z range: {}'.format(zRange))
    else:
        raise Exception('z must be a float or an int')

    # Now the zW variable
    if isinstance(zW, intOrFloat):
        if not inRange(zW, zWRange):
            raise Exception('zW value outside of zW range: {}'.format(zWRange))
    else:
        raise Exception('zW must be a float or an int')

    # Now the appOrAbs variable
    if type(appOrAbs) != str:
        raise Exception('appOrAbs must be \'apparent\' or \'absolute\'')

    # Now the CMF_method variable
    if type(CMF_method) != str:
        raise Exception('CMF_method must be \'nu-scaling\' or \'PS-scaling\'')

    # Now for interpWarning
    if isinstance(interpWarning, intOrFloat):
        if interpWarning not in [0, 1, 2]:
            raise Exception('interpWarning be equal to 0, 1, or 2')
    else:
        raise Exception('interpWarning must be int or float')


def log10Eps(area, Psi, gamma, b):
    '''
    This function returns log10(epsilon_cv) given the fit parameters Psi, gamma, and b

    Parameters
    -------------------------
    area : int or float
        The area of the survey
    Psi, gamma, b: int or float
        The parameters of the varepsiloncv fit
    '''
    return Psi * area**gamma + b


def lincv(mass, area, z, zW=1, message='yes'):
    '''
    Warning! Use with caution and only if outside the bounds of 'galcv.getcv()'. This function is designed to be used at larger areas and larger masses (brighter galaxies) than galcv.getcv(). For additional questions see the README on GitHub. This function returns the 1-sigma linear approximation of cosmic variance for haloes of the chosen mass (in solar masses) in the chosen volume at the chosen redshift. Also, this method assumes the survey volume is a sphere. If your survey volume is actually very elongated in some direction, this method will overestimate cosmic variance.

    Parameters
    -------------------------
    mass : int or float or array-like
        Mass of a halo (in units of solar mass)
    area : int or float
        Survey area in square arcminutes
    z : int or float
        Central redshift (range: 4 - 15)
    zW : int or float
        Redshift width (default = 1; range: 0.1 - 2)
    message : 'yes' or 'no'
        Whether or not to print the warning message

    Returns
    -------------------------
    A NumPy list of cosmic variance values of the same length as the mass input
    '''
    if message == 'yes':
        print('Warning! Use with caution and only if outside the bounds of \'galcv.getcv()\'. This function is designed to be used at larger areas and larger masses (brighter galaxies) than galcv.getcv(). For additional questions see the README on GitHub. This function returns the 1-sigma linear approximation of cosmic variance for haloes of the chosen mass (in solar masses) in the chosen volume at the chosen redshift. Also, this method assumes the survey volume is a sphere. If your survey volume is actually very elongated in some direction, this method will overestimate cosmic variance.')

    tckS = np.load(os.path.dirname(os.path.abspath(__file__)) + '/sigma0fM_interp.npy', allow_pickle=True)
    # sigma associated with the halo mass
    sigM = 10**interpolate.splev(np.log10(mass), tckS) * growthFac(z=z)

    alpha = -1.36
    beta = -1.14

    # Calculating the Trac linear bias
    bTR = 1
    bTR += (alpha / delCrit0(z=z)) * (sigM / 2.54)**alpha / (1 + (sigM / 2.54)**alpha)
    bTR -= (2 * beta) / (sigM**2 * delCrit0(z=z))

    CRITDENMSOLMPC = 2.7755e11
    # Calculating the sigma associated with the survey area
    s_sr = arcmin2ToSr(arcmin2=area)
    s_Vol = volCM(z1=z - 0.5 * zW, z2=z + 0.5 * zW, Omega=s_sr)
    s_RegMass = s_Vol * CRITDENMSOLMPC * om0hh
    sigM_Reg = 10**interpolate.splev(np.log10(s_RegMass), tckS) * growthFac(z=z)

    # The cosmic variance is the survey-wide 1-sigma fluctuation times the bias factor
    if z < 4 or z > 15:
        print('\n\n\nAnother Warning! You have chosen a redshift outside the recommended range of 4 - 15. I\'ll still output an answer but be wary.')
    if area < 3e4:
        print('\n\n\nAnother Warning! This area is covered in getcv(); you may be able to use that instead.')
    if zW/z > 0.25:
        print('\n\n\nAnother Warning! The redshift bin width is getting large compared to the redshift.')
    if zW > 2:
        print('\n\n\nAnother Warning! The redshift bin width is getting large. Recommend to keep it below 2')
    if np.any(mass > 1e14):
        print('\n\n\nAnother Warning! There is a very massive halo here! Are you sure that is what you wanted?')
    if np.any(mass < 4e8):
        print('\n\n\nAnother Warning! There is a very small-mass halo here... Are you sure that is what you wanted?')
    return bTR * sigM_Reg

###Now to define various cosmology equations for use in the main methods

def absToApp(Mabs='ER', z='ER'):
    '''
    This function converts absolute magnitude into apparent

    Parameters
    -------------------------
    Mabs : int or float or array-like
        absolute magnitudes
    z : int or float
        redshift
    '''

    DL_parsec = lum_dist(z=z)
    # Convert to parsecs
    DL_parsec = DL_parsec * 1.e6

    mapp = Mabs
    mapp = mapp + 5.0 * np.log10(DL_parsec / 10)
    mapp = mapp - 2.5 * np.log10(1.0 + z)

    return mapp


def lum_dist(z='ER', **kwargs):
    '''
    This function returns the luminosity distance to some redshift

    Parameters
    -------------------------
    z : int or float
        redshift
    '''
    ans = comv_dist(z=z) * (1.0 + z)

    return ans

HPARAM = 0.678  # hubble constant today/100
OMEGA0 = 0.308  # matter fraction today
OMEGANU = 0.0  # radiation fraction today
LAMBDA0 = 0.692  # dark energy fraction today
SPEEDOFLIGHTMPCPERSEC = 9.71561e-15
TCMB = 2.726  # CMB temp today
thetaCMB = TCMB / 2.7
om0hh = OMEGA0 * HPARAM * HPARAM
zEquality = 2.50e4 * om0hh * pow(thetaCMB, -4.0) - 1.0


def comv_dist(z='ER', **kwargs):
    '''
    This function returns the comoving distance to some redshift

    Parameters
    -------------------------
    z : int or float
        redshift
    '''

    def wrapper_comv_dist(x):
        # What is the hubble constant today in 1/s
        Hnot = 100.0 * HPARAM * 3.24078e-20
        Hofz = Hnot * np.sqrt(OMEGA0 * (1. + x)**3. + OMEGANU * (1. + x)**4. + LAMBDA0)
        ansWR = SPEEDOFLIGHTMPCPERSEC / Hofz
        return ansWR

    ans, abserr = integrate.quad(wrapper_comv_dist, 0.0, z)

    if abs(abserr / ans) > 1e-4:
        print('Warning! Comoving distance calculation err is high')
        print('err/value = ' + str(abs(abserr / ans)))

    return ans


def growthFac(z):

    omZ = omegaZ(z=z)
    lamZ = lambdaZ(z=z)

    D = ((1.0 + zEquality) / (1.0 + z) * 5.0 * omZ / 2.0 * pow(pow(omZ, 4.0 / 7.0) - lamZ + (1.0 + omZ / 2.0) * (1.0 + lamZ / 70.0), -1.0))
    D = D / ((1.0 + zEquality) * 5.0 / 2.0 * OMEGA0 * pow(pow(OMEGA0, 4.0 / 7.0) - LAMBDA0 + (1.0 + OMEGA0 / 2.0) * (1.0 + LAMBDA0 / 70.0), -1.0))
    return D


def omegaZ(z):

    z = z + 1.
    temp = OMEGA0 * z**3.0
    temp = temp / (temp + LAMBDA0 + (1. - OMEGA0 - LAMBDA0) * z**2)

    return temp


def lambdaZ(z):
    z = z + 1.
    temp = LAMBDA0
    temp = temp / (temp + OMEGA0 * z**3 + ((1. - OMEGA0 - LAMBDA0) * z**2))
    return temp


def delCrit0(z):

    t1 = 0.15 * (12 * np.pi)**(2. / 3.)
    omZ = 0 + omegaZ(z=z)

    if abs(omZ - 1) < 1.0e-5:
        return t1
    elif abs(LAMBDA0) < 1.0e-5:
        t1 = t1 * omZ**0.0185
        return t1
    else:
        t1 = t1 * omZ**0.0055
        return t1


def arcmin2ToSr(arcmin2):

    steradians = arcmin2 * (1. / (60. * 180. / np.pi))**2.

    return steradians


def volCM(z1='ER', z2='ER', Omega='ER'):

    if z1 == 0:
        print('Warning 8623478: You can\'t have z1=0. Pick some small value instead if you must.')

    # First, interpolate the comoving radius between the redshift bounds
    the_dz = 0.1
    the_zs = np.arange(z1, z2 + 2 * the_dz, the_dz)
    the_rcoms = np.zeros(len(the_zs))
    for i in range(len(the_zs)):
        the_rcoms[i] = comv_dist(the_zs[i])
    tckRCOM = interpolate.splrep(the_zs, the_rcoms, k=1, s=0)

    def wrapper_volCM(x):
        Hnot = 100.0 * HPARAM * 3.24078e-20
        Hofz = Hnot * np.sqrt(OMEGA0 * (1. + x)**3. + OMEGANU * (1. + x)**4. + LAMBDA0)
        c_over_Hofz = SPEEDOFLIGHTMPCPERSEC / Hofz

        ansWR = c_over_Hofz
        ansWR = ansWR * interpolate.splev(x, tckRCOM)**2

        return ansWR

    ans, abserr = integrate.quad(wrapper_volCM, z1, z2, limit=100)

    if abs(abserr / ans) > 1e-4:
        print('Warning! Volum (comoving) calculation err is high')
        print('err/value = ' + str(abs(abserr / ans)))

    return ans * Omega
