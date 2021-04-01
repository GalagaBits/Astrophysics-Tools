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

a_1,a_2,a_3,a_4,a_5      = -0.00061,-0.00057,-0.00062,-0.00052,-0.00056
w_1,w_2,w_3,w_4,w_5      =3.78e3,3e3,3.3e3,3.1e3,3.1e3
x = np.linspace(-100000*7.1889609587E+01,100000*7.1889609587E+01,10000000)


pl.plot(x, gaussian(x,a_1,w_1,0), color="blue")
pl.plot(x, gaussian(x,a_2,w_2,0), color="red")
pl.plot(x, gaussian(x,a_3,w_3,0), color="yellow")
pl.plot(x, gaussian(x,a_4,w_4,0), color="orange")
pl.plot(x, gaussian(x,a_5,w_5,0), color="green")

pl.title('Period Folding of Gaussian Fit Transits', fontname="Times New Roman")
pl.xlabel('Time (s)', fontname="Times New Roman")
pl.ylabel('Normalized Flux',fontname="Times New Roman")
pl.show()

#Overplot Observations
#'''
import pylab as pl
from astropy.timeseries import TimeSeries


transit1 = yval[nearest_index(xval, -15300):nearest_index(xval, 24700)]
transit2 = yval[nearest_index(xval, 2380000):nearest_index(xval, 2420000)]
transit3 = yval[nearest_index(xval, 4780000):nearest_index(xval, 4820000)]
transit4 = yval[nearest_index(xval, 7170000):nearest_index(xval, 7210000)]
transit5 = yval[nearest_index(xval, 9570000):nearest_index(xval, 9610000)]

orbital_period = 27.7435 #days


ts1 = TimeSeries(time_start='2000-01-01T00:00:00.000',
                time_delta=7.1889609587E+01 * u.s,
                data={'normalized flux': transit1})


ts2 = TimeSeries(time_start='1999-12-31T19:45:00.000',
                 time_delta=7.1889609587E+01 * u.s,
                 data={'normalized flux': transit2})

ts3 = TimeSeries(time_start='2000-01-28T14:51:38.000',
                time_delta=7.1889609587E+01 * u.s,
                data={'normalized flux': transit3})

ts4 = TimeSeries(time_start='2000-02-25T07:03:16.000',
                time_delta=7.1889609587E+01 * u.s,
                data={'normalized flux': transit4})

ts5 = TimeSeries(time_start='2000-03-24T20:32:54.000',
                time_delta=7.1889609587E+01 * u.s,
                data={'normalized flux': transit5})



ts1_folded = ts1.fold(period=2397038*u.s, epoch_time='2000-01-01T01:16:00')
ts2_folded = ts2.fold(period=2397038*u.s, epoch_time='2000-01-28T19:06:38')
ts3_folded = ts3.fold(period=2397038*u.s, epoch_time='2000-02-25T12:57:16')
ts4_folded = ts4.fold(period=2397038*u.s, epoch_time='2000-03-24T06:47:54')
ts5_folded = ts5.fold(period=2397038*u.s, epoch_time='2000-04-21T00:38:32')


pl.plot(ts1_folded.time.jd, ts1_folded['normalized flux'], 'rx', markersize=2, drawstyle='steps-pre',label="Transit 1")
pl.plot(ts2_folded.time.jd, ts2_folded['normalized flux'], 'gx', markersize=2, drawstyle='steps-pre',label="Transit 2")
pl.plot(ts3_folded.time.jd, ts3_folded['normalized flux'], 'yx', markersize=2, drawstyle='steps-pre',label="Transit 3")
pl.plot(ts4_folded.time.jd, ts4_folded['normalized flux'], 'bx', markersize=2, drawstyle='steps-pre',label="Transit 4")
pl.plot(ts4_folded.time.jd, ts4_folded['normalized flux'], 'mx', markersize=2, drawstyle='steps-pre',label="Transit 5")

#pl.title('Period Folding of Gaussian Fit Transits', fontname="Times New Roman")
pl.xlabel('Phase (s)', fontname="Times New Roman")
pl.ylabel('Normalized Flux',fontname="Times New Roman")
pl.legend()
pl.savefig("/Users/dealderod/Documents/GitHub/Astrophysics-Tools/Plots/Transits_Folding")

pl.show() 
#'''



















