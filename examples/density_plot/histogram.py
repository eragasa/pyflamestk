# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 13:40:47 2016

@author: Eugene
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

xedges = [0, 1, 1.5, 3, 5]
yedges = [0, 2, 3, 4, 6]

x = np.random.normal(3, 1, 10000)
y = np.random.normal(1, 1, 10000)
H, xedges, yedges = np.histogram2d(y, x, bins=(xedges, yedges))

fig = plt.figure(figsize=(7, 3))
ax = fig.add_subplot(131)
ax.set_title('imshow: equidistant')
im = plt.imshow(H, interpolation='nearest', origin='low',
                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
                
                
ax = fig.add_subplot(132)
ax.set_title('pcolormesh: exact bin edges')
X, Y = np.meshgrid(xedges, yedges)
ax.pcolormesh(X, Y, H)
ax.set_aspect('equal')        
                
ax = fig.add_subplot(133)
ax.set_title('NonUniformImage: interpolated')
im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
im.set_data(xcenters, ycenters, H)
ax.images.append(im)
ax.set_xlim(xedges[0], xedges[-1])
ax.set_ylim(yedges[0], yedges[-1])
ax.set_aspect('equal')
plt.show()