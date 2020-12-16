# galcv

This package provides predictions of cosmic variance for the high-z UV luminosity function (UVLF) of galaxies. The methods for this code are described in Trapp & Furlanetto (2020).

This package returns a linear approximation of the relative cosmic variance of the UVLF (1 sigma) for the following parameter ranges:

#### Apparent rest-UV AB magnitude: 22 -> 34

#### Redshift: 5 -> 15

#### Survey Area \(sqr arcmin\): 1 -> 31640

**Note:** Cosmic variance is not available for *all* combinations of these parameters, even within these ranges. This occurs most often at the lowest survey areas and brightest apparent magnitudes. The code will print a warning if this is the case.

---
**ANNOUNCEMENT:** The new method `lincv()` is similar to getcv but can be used at wider parameter ranges: Redshift 4 -> 15; Survey Area \(sqr arcmin\) 1e4 -> 1e8, and halo masses higher than 1e13 solar masses (abs mag < -24). However, it is only designed for limited cases and is a function of halo mass, not galaxy luminosity. See below for more details or email me at atrapp@astro.ucla.edu
---

---
# Installation and Use

The *simplest* way to install and use `galcv` is through 'pip' in a python environment:
```python
> pip install galcv
```

The package can then be imported in any python environment or in a script using:
```python
> import galcv
```

There are currently two user-facing functions: `getcv()` is the main function; `lincv()` is a new addition that is more approximate but can be used with larger volumes, more massive galaxies, and goes down to z = 4. The rest of the functions are intended for internal use. Example use:
```python
> galcv.getcv(mag=[30,29,28], area=100, z=9)
> [0.178, 0.208, 0.245]

> galcv.lincv(mass=[1e11,1e12,1e13], area=1e7, z=4)
> array([0.0004263 , 0.00072753, 0.00142989])
```

`getcv()` takes three required parameters (mag, area, z), and has four default parameters (zW, appOrAbs, CMF_method, interpWarning). The following is the docstring for `getcv()` that explains the inputs and output:

```
This function returns relative cosmic variance results. This function is a wrapper function for formatting. The actual calculation happens in singlecv()

Parameters
-------------------------
mag : int, float, list, or numpy.ndarray
    The magnitude(s) to consider. This must be in APPARENT rest-UV (1500 - 2800 Angstroms) AB magnitude
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
```

`lincv()` takes three required parameters (mass, area, z), and has two default parameters (zW, message). `lincv()` simply outputs the linear halo bias function (derived from Trac et al. 2015 halo mass function) multiplied by the 1-sigma rms fluctuation of the dark matter density field on the scale of the survey. **If you want to use `lincv()`: be careful! It is much more approximate and is only appropriate in limited circumstances. Email atrapp@astro.ucla.edu or comment on GitHub to see if this is actually useful for your purposes.** The following is the docstring for `lincv()` that explains the inputs and output:

```
Warning! Use with caution and only if outside the bounds of 'galcv.getcv()'. This function is designed to be used at larger areas and larger masses (brighter galaxies) than galcv.getcv(). In these regions, Poisson noise SHOULD be dominating anyway. For additional questions please comment in the GitHub. This function returns the 1-sigma linear approximation of cosmic variance for haloes of the chosen mass (in solar masses) in the chosen volume at the chosen redshift. Note: you must use your own halo-mass to luminosity relation if you want to connect to the UV luminosity function. Also, this method assumes the survey volume is a sphere. If your survey volume is actually very elongated in some direction, this method will overestimate cosmic variance.

Parameters
-------------------------
mass : int or float or array-like of ints and floats
    Mass of a halo (in units of solar mass)
area : int or float
    Survey area in square arcminutes
z : int or float
    Central redshift
zW : int or float
    Redshift width (default = 1)
message : 'yes' or 'no'
    Whether or not to print the warning message

Returns
-------------------------
A NumPy list of cosmic variance values of the same length as the mass input
```
    
---
# Alternate Installation and Use Methods

If 'pip' is not working, or you would prefer to run the code yourself, you may clone the github repo and run the \_\_init\_\_.py script (in the /galcv folder) in a python environment. You will then have access to the `getcv()` function.

In fact, all the code needs to run is the \_\_init\_\_.py script along with all of the .pkl files that are in the /galcv folder.

If you would like to use the `getcv()` function in your script (without installing it with pip and importing it), you may do the following:
- Copy the \_\_init\_\_.py file into the same directory as your script
- Also copy *all* of the .pkl files from the /galcv folder to that same directory
- At the beginning of your script, include the line:
```python
from __init__ import *
```
- You should then be able to use `getcv()` or `lincv()` in that script.

---
# Links

GitHub: https://github.com/adamtrapp/galcv

PyPi: https://pypi.org/project/galcv/
