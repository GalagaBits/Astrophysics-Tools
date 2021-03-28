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

# sp = pyspeckit.Spectrum(data=yval, xarr=xval)

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


#Plot the First Transit
import pylab as pl

fig = pl.figure(figsize=(8,6))
sp.plotter(color='black',xmin=0.0e7, xmax=0.010e7 ,ymin=0.99900, ymax=1.00075,figure=fig)

pl.xlabel('Time (s)', fontname="Times New Roman")
pl.ylabel('Normalized Flux',fontname="Times New Roman")

sp.baseline(interactive=True, subtract=False)

sp.specfit(interactive=True)


#Overplot Transits with Phase shifts

'''  
Helpful Resources

https://doi.org/10.1051/0004-6361/201834672  
https://www.qmul.ac.uk/spa/media/school-of-physics/outreach/research-in-schools/PlanetHuntingWithPythonManual-StudentVersion-UpdatedDataDownloadInstructions.pdf

To find the Phase shift, find the orbital period, divde it my 2, and substract 
it by all of the values in the set.

Notes to self
 - Splice and make sets for each transit
 - convert the seconds to days using astropy
 - Plot the points of observation and gaussians for each plot 

'''

orbital_period = 597.8 #days

def gaussian(xval, a, FWHM, rest_f):
    return (a*(np.exp(-np.log(2)*(((xval-rest_f)**2)/((0.5*FWHM)**2)))))

a_1,a_2,a_3,a_4,a_5      = -0.00061,-0.00057,-0.00062,-0.00052,-0.00056
w_1,w_2,w_3,w_4,w_5      =3.78e3,3e3,3.3e3,3.1e3,3.1e3
x = np.linspace(-5,5,10)

pl.plot(x, gaussian(x,a_1,w_2,0), color="blue")
# pl.plot(x, gaussian(x,a_2,w_2,0), color="yellow")
# pl.plot(x, gaussian(x,a_3,w_3,0), color="orange")
# pl.plot(x, gaussian(x,a_4,w_4,0), color="red")
# pl.plot(x, gaussian(x,a_5,w_5,0), color="green")
pl.show()
# def transit_std(flux, t1, t2):
#     t_start, t_end = (find_nearest(flux, t1),
#                       find_nearest(flux, t2))
#     transit = data[t_start:t_end]
#     np.std(transit)
