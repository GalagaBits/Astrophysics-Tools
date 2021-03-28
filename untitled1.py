#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 17:02:15 2021

@author: dealderod
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:35:58 2021

@author: dealderod
"""
import numpy as np
import pyspeckit
import pandas as pd

data = pd.read_csv('/Users/dealderod/Documents/Number10_1.csv', header=None, delimiter=',')

print(data.tail(6))

xval = np.array(data.iloc[:,0])
    

yval = np.array(data.iloc[:,1])

print(xval,yval)

sp = pyspeckit.Spectrum(data=yval, xarr=xval)


import pylab as pl

fig = pl.figure(figsize=(8,6))
sp.plotter(color='black',xmin=0.0e7, xmax=0.010e7 ,ymin=0.99900, ymax=1.00075,figure=fig)

pl.xlabel('Time (s)', fontname="Times New Roman")
pl.ylabel('Normalized Flux',fontname="Times New Roman")

sp.baseline(interactive=True, subtract=False)

sp.specfit(interactive=True)
