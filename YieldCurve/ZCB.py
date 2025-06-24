import numpy as np
from scipy.interpolate import splrep, splev, interp1d

def P0TModel(t,ti,ri,method):
    rInterp = method(ti,ri)
    r = rInterp(t)
    return np.exp(-r*t)
    
def linear_interpolation(ti,ri):
    interpolator = lambda t: np.interp(t, ti, ri)
    return interpolator
    
def scipy_1d_interpolate(ti, ri):
    interpolator = lambda t: interp1d(ti, ri, kind='cubic',fill_value="extrapolate")(t)
    return interpolator    
        
def spline_interpolate(ti,ri):
    interpolator = splrep(ti, ri, s=0.01)
    interp = lambda t: splev(t,interpolator)
    return interp

########## excel data
data = xl("RawData!V4:AI110")
data = data.to_numpy()

ri   = np.array(data[0,5:14])/100
mat = np.array([1.0,2.0,3.0,5.0,7.0,10.0,20.0,30.0,100.0])

method = scipy_1d_interpolate

P0T = lambda t: P0TModel(t,mat,ri,method)
