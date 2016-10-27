# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:33:17 2016

@author: Eugene
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as ml
import numpy as np

n = 1e5
x = y = np.linspace(-5, 5, 100)
X, Y  = np.meshgrid(x, y)
Z1 = ml.bivariate_normal(X, Y, 2, 2, 0, 0)
Z2 = ml.bivariate_normal(X, Y, 4, 1, 1, 1)
ZD = Z2 - Z1
x = X.ravel()
y = Y.ravel()
z = ZD.ravel()
gridsize=30
plt.subplot(111)

# if 'bins=None', then color of each hexagon corresponds directly to its count
# 'C' is optional--it maps values to x-y coordinates; if 'C' is None (default) then 
# the result is a pure 2D histogram 

plt.hexbin(x, y, C=z, gridsize=gridsize, cmap=cm.jet, bins=None)
plt.axis([x.min(), x.max(), y.min(), y.max()])

cb = plt.colorbar()
cb.set_label('mean value')
plt.show()   
