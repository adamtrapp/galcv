# galcv

### WARNING: This package is still in its alpha stage

This package provides predictions of cosmic variance for the high-z UV luminosity function (UVLF) of galaxies. The methods for this code are described in Trapp & Furlanetto (2020, in prep.).

This package provides the relative cosmic variance of the UVLF for the following parameter ranges:

#### Apparent rest-UV AB magnitude: 22 -> 34

#### Redshift: 5 -> 15

#### Survey Area \(sqr arcmin\): 1 -> 31640

**Note:** Cosmic variance is not available for *all* combinations of these parameters, even within these ranges. This occurs most often at the lowest survey areas and brightest apparent magnitudes. The code will print a warning if this is the case.

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

There is currently one user-facing function: `getcv()`. The rest of the functions are intended for internal use. Example use:
```python
> galcv.getcv(mag=[30,29,28], area=100, z=9)
> [0.178, 0.208, 0.245]
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
- You should then be able to use `getcv()` in that script.

---
# Links

GitHub: https://github.com/adamtrapp/galcv

PyPi: https://pypi.org/project/galcv/
