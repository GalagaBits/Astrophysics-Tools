#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 20:31:07 2021

@author: dealderod
"""

import numpy as np
import pyspeckit
import pandas as pd
from astropy import units as u

data = pd.read_csv('/Users/dealderod/Documents/Number10_1.csv', header=None, delimiter=',')

xval = np.array(data.iloc[:,0])    

yval = np.array(data.iloc[:,1])

print(xval)
sp = pyspeckit.Spectrum(data=yval, xarr=xval)

#
flux=yval
t1 = 0
t2 = 9.12E3


#function pulled from the user "unutbu" on overflow on Stack
def nearest_index(array, value):
    array = np.asarray(array)
    idx = np.abs(array - value).argmin()
    return idx


def transit_std(flux, t1, t2):
    t_start, t_end = (nearest_index(xval, t1),
                      nearest_index(xval, t2))
    transit = flux[t_start:t_end]
    return np.std(transit)


#Plot the First Transit

import pylab as pl

# fig = pl.figure(figsize=(8,6))
# sp.plotter(color='black',xmin=0.0e7, xmax=0.010e7 ,ymin=0.99900, ymax=1.00075,figure=fig)

# pl.xlabel('Time (s)', fontname="Times New Roman")
# pl.ylabel('Normalized Flux',fontname="Times New Roman")
# pl.show()
# # sp.baseline(interactive=True, subtract=False)

# sp.specfit(interactive=True)


#Overplot Transits with Phase shifts

'''  
Helpful Resources

# Use this link from astropy: https://docs.astropy.org/en/stable/timeseries/

https://doi.org/10.1051/0004-6361/201834672  
https://www.qmul.ac.uk/spa/media/school-of-physics/outreach/research-in-schools/PlanetHuntingWithPythonManual-StudentVersion-UpdatedDataDownloadInstructions.pdf

To find the Phase shift, find the orbital period, divde it my 2, and substract 
it by all of the values in the set.

Notes to self
 - Splice and make sets for each transit
 - convert the seconds to days using astropy
 - Plot the points of observation and gaussians for each plot 

'''
import pylab as pl
from astropy.timeseries import TimeSeries

'''
From astropy...

Parameters

    datanumpy ndarray, dict, list, Table, or table-like object, optional

        Data to initialize time series. This does not need to contain the times,
        which can be provided separately, but if it does contain the times they 
        should be in a column called 'time' to be automatically recognized.
        
    fold(period=None, epoch_time=None, epoch_phase=0, wrap_phase=None,
         normalize_phase=False)[source]
'''

transit1 = yval[nearest_index(xval, -15300):nearest_index(xval, 24700)]
transit2 = yval[nearest_index(xval, 2380000):nearest_index(xval, 2420000)]
# transit3 = yval[nearest_index(xval, 4780000):nearest_index(xval, 4820000)]
# transit4 = yval[nearest_index(xval, 7170000):nearest_index(xval, 7210000)]
# transit5 = yval[nearest_index(xval, 9570000):nearest_index(xval, 9610000)]

orbital_period = 27.7435 #days


ts1 = TimeSeries(time_start='2000-01-01T00:00:00',
                time_delta=7.1889609587E+01 * u.s,
                data={'normalized flux': transit1})


ts2 = TimeSeries(time_start='2000-02-25T11:41:16.000',
                 time_delta=7.1889609587E+01 * u.s,
                 data={'normalized flux': transit2})

# ts3 = TimeSeries(time_start='2000-02-25T11:41:16.000',
#                 time_delta=7.1889609587E+01 * u.s,
#                 data={'normalized flux': transit3})

# ts4 = TimeSeries(time_start='2000-02-25T11:41:16.000',
#                 time_delta=7.1889609587E+01 * u.s,
#                 data={'normalized flux': transit4})

# ts5 = TimeSeries(time_start='2000-02-25T11:41:16.000',
#                 time_delta=7.1889609587E+01 * u.s,
#                 data={'normalized flux': transit5})



ts1_folded = ts1.fold(period=2397038*u.s, epoch_time='2000-01-28T17:50:38.000')
ts2_folded = ts2.fold(period=2397038*u.s, epoch_time='2000-02-25T11:41:16.000')
# ts3_folded = ts3.fold(period=2397038*u.s, epoch_time='2000-03-24T05:31:54.000')
# ts4_folded = ts4.fold(period=2397038*u.s, epoch_time='2000-04-20T23:22:32.000')
# ts5_folded = ts5.fold(period=2397038*u.s, epoch_time='2000-05-18T17:13:10.000')

def gaussian(xval, a, sigma, time):
    return (a*(np.exp(-np.log(2)*(((xval-time)**2)/((0.5*(sigma*np.sqrt(8*np.log(2))))**2)))))

pl.plot(ts1_folded.time.jd, gaussian(ts1_folded.time.jd,-0.00057,3e3,0), color="blue")

# # pl.plot(ts1_folded.time.jd, ts1_folded['normalized flux'], 'k.', markersize=1, drawstyle='steps-pre')
# pl.plot(ts2_folded.time.jd, ts2_folded['normalized flux'], 'b-', markersize=1, drawstyle='steps-pre')
# # pl.plot(ts3_folded.time.jd, ts3_folded['normalized flux'], 'g-', markersize=1, drawstyle='steps-pre')
# # pl.plot(ts4_folded.time.jd, ts4_folded['normalized flux'], 'r-', markersize=1, drawstyle='steps-pre')
# # pl.plot(ts4_folded.time.jd, ts4_folded['normalized flux'], 'y-', markersize=1, drawstyle='steps-pre')

# pl.show() 
