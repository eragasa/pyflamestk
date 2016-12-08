import numpy as np
from scipy import stats

def measure(n):
    "Measurement model, return two coupled measurements."
    m1 = np.random.normal(size=n)
    m2 = np.random.normal(scale=0.5, size=n)
    return m1+m2, m1-m2


# generate some random 2d data    
m1, m2 = measure(2000)
xmin = m1.min()
xmax = m1.max()
ymin = m2.min()
ymax = m2.max()

# perform kernel density estimate on data
X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([m1, m2])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)


# resample
n1, n2 = kernel.resample(size=2000)
nvalues = np.vstack([n1,n2])
nkernel = stats.gaussian_kde(values)
nZ = np.reshape(kernel(positions).T, X.shape)

import matplotlib.pyplot as plt

plt.subplot(1,2,1)
plt.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
           extent=[xmin, xmax, ymin, ymax])
plt.plot(m1, m2, 'k.', markersize=2)
plt.axis([xmin,xmax,ymin,ymax])
#plt.set_xlim([xmin, xmax])
#plt.set_ylim([ymin, ymax])

plt.subplot(1,2,2)
plt.imshow(np.rot90(nZ), cmap=plt.cm.gist_earth_r,
           extent=[xmin, xmax, ymin, ymax])
plt.plot(n1, n2, 'k.', markersize=2)
plt.axis([xmin,xmax,ymin,ymax])
#plt.set_xlim([xmin, xmax])
#plt.set_ylim([ymin, ymax])

plt.show()

