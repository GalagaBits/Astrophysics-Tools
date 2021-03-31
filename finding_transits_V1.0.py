#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 20:31:07 2021

@author: dealderod
"""

import numpy as np
import pandas as pd
from astropy import units as u

data = pd.read_csv('/Users/dealderod/Documents/Number10_1.csv', header=None, delimiter=',')

xval = np.array(data.iloc[:,0])    

yval = np.array(data.iloc[:,1])

#function from unutbu on Stack
def nearest_index(array, value):
    array = np.asarray(array)
    idx = np.abs(array - value).argmin()
    return idx


def transit_std(flux, t1, t2):
    t_start, t_end = (nearest_index(xval, t1),
                      nearest_index(xval, t2))
    transit = flux[t_start:t_end]
    return np.std(transit)

def gaussian(xval, a, sigma, time):
    return (a*(np.exp(-np.log(2)*(((xval-time)**2)/((0.5*(sigma*np.sqrt(8*np.log(2))))**2)))))


#Interactive Gaussian Fit of Transits
'''
import pyspeckit
import pylab as pl

sp = pyspeckit.Spectrum(data=yval, xarr=xval)

fig = pl.figure(figsize=(8,6))
sp.plotter(color='black',xmin=0.0e7, xmax=0.010e7 ,ymin=0.99900, ymax=1.00075,figure=fig)

pl.xlabel('Time (s)', fontname="Times New Roman")
pl.ylabel('Normalized Flux',fontname="Times New Roman")

# sp.baseline(interactive=True, subtract=False)

# sp.specfit(interactive=True)
'''

#Overplot Transits with Phase shifts

import pylab as pl
from astropy.timeseries import TimeSeries


orbital_period = 27.7435 #days


ts = TimeSeries(time_start='2016-03-22T12:30:31',
                time_delta=7.1889609587E+01 * u.s,
                data={'normalized flux': yval})

ts_folded = ts.fold(period=orbital_period*u.s) 

a_1,a_2,a_3,a_4,a_5      = -0.00061,-0.00057,-0.00062,-0.00052,-0.00056
w_1,w_2,w_3,w_4,w_5      =3.78e3,3e3,3.3e3,3.1e3,3.1e3
x = np.linspace(-100000*7.1889609587E+01,100000*7.1889609587E+01,10000000)
pl.plot(x, gaussian(x,a_1,w_1,0), color="blue")
pl.plot(x, gaussian(x,a_2,w_2,0), color="red")
pl.plot(x, gaussian(x,a_3,w_3,0), color="yellow")
pl.plot(x, gaussian(x,a_4,w_4,0), color="orange")
pl.plot(x, gaussian(x,a_5,w_5,0), color="green")

pl.xlabel('Time (s)', fontname="Times New Roman")
pl.ylabel('Normalized Flux',fontname="Times New Roman")

pl.show()
