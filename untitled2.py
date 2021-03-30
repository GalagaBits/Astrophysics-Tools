#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 19:15:21 2021

@author: dealderod
"""

import pylab as pl
import pandas as pd
import numpy as np
from astropy.timeseries import TimeSeries
from astropy import units as u

from astropy.timeseries import BoxLeastSquares

data = pd.read_csv('/Users/dealderod/Documents/Number10_1.csv', header=None, delimiter=',')

xval = np.array(data.iloc[:,0])    
yval = np.array(data.iloc[:,1])

#function pulled from the user "unutbu" on overflow on Stack
def nearest_index(array, value):
    array = np.asarray(array)
    idx = np.abs(array - value).argmin()
    return idx

def gaussian(xval, a, sigma, time):
    return (a*(np.exp(-np.log(2)*(((xval-time)**2)/((0.5*(sigma*np.sqrt(8*np.log(2))))**2)))))


transit2 = yval[nearest_index(xval, 2380000):nearest_index(xval, 2420000)]

ts = TimeSeries(time_start='2000-02-25T11:41:16',
                 time_delta=7.1889609587E+01 * u.s,
                 data={'normalized flux': transit2})

periodogram = BoxLeastSquares.from_timeseries(ts)  
results = periodogram.autopower(0.2 * u.day)  
best = np.argmax(results.power)  
period = results.period[best]  

transit_time = results.transit_time[best]  




ts2_folded = ts2.fold(period=2397038*u.s, epoch_time='2000-01-01T00:00:00')

pl.plot(ts2_folded.time.jd, ts2_folded['normalized flux'], 'k.', markersize=1, drawstyle='steps-pre')

