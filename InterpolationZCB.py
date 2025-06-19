#Naive way to create a ZCB that Match Mkt spine points

import numpy as np
from scipy.interpolate import splrep, splev, 

def P0TModel(t,ti,ri,method):
    rInterp = method(ti,ri)
    r = rInterp(t)
    return np.exp(-r*t)

def linear_interpolation(ti,ri):
    interpolator = lambda t: np.interp(t, ti, ri)
    return interpolator

ri   = np.array(...)/100
mat = np.array([1.0,2.0,3.0,5.0,7.0,10.0,20.0,30.0])

method = linear_interpolation

P0T = lambda t: P0TModel(t,mat,ri,method)
