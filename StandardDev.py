#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 20:31:07 2021

@author: dealderod
"""

import numpy as np
import pyspeckit
import pandas as pd

data = pd.read_csv('/Users/dealderod/Documents/Number10_1.csv', header=None, delimiter=',')

xval = np.array(data.iloc[:,0])    

yval = np.array(data.iloc[:,1])

sp = pyspeckit.Spectrum(data=yval, xarr=xval)

#
flux=yval
t1 = 0
t2 = 9.12E3


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


import pylab as pl

fig = pl.figure(figsize=(8,6))
sp.plotter(color='black',xmin=0.0e7, xmax=0.010e7 ,ymin=0.99900, ymax=1.00075,figure=fig)

pl.xlabel('Time (s)', fontname="Times New Roman")
pl.ylabel('Normalized Flux',fontname="Times New Roman")

sp.baseline(interactive=True, subtract=False)

sp.specfit(interactive=True)


# def transit_std(flux, t1, t2):
#     t_start, t_end = (find_nearest(flux, t1),
#                       find_nearest(flux, t2))
#     transit = data[t_start:t_end]
#     np.std(transit)
