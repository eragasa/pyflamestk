def f1(x):
    x1 = x[0]
    x2 = x[1]
    return x1**4+x2**4-x1**2+x2**2 - 10*x1*x2 + 0.25*x1 +20
    
def f2(x):
    x1 = x[0]
    x2 = x[1]
    return (x1-1)**2+x2**2

def C(x):
    f_1 = f1(x)
    f_2 = f2(x)
    C = w1*f_1 + w2*f_2
    
    return C
    
w1 = 1
w2 = 0

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

x0 = np.array([1,1])
res1 = optimize.fmin_cg(C, x0, full_output=True, retall=True)
    
results = []
weights = np.arange(1,10,0.1)
weights = np.hstack((weights,1/weights))
x01s = np.arange(-2,2,.1)
x02s = np.arange(-2,2,.1)


for x01 in x01s:
    for x02 in x02s:
        x0 = np.array([x01,x02])
        for i, w2 in enumerate(weights):
            result = optimize.fmin_cg(C, x0, full_output=True, retall=True)
            x1 = result[0][0]
            x2 = result[0][1]
            new_res = [w1,w2,x1,x2,f1([x1,x2]),f2([x1,x2])]
            if i == 0:
                res = new_res
            else:
                res = np.vstack((res,new_res))
            results.append(result)
        plt.scatter(res[:,4],res[:,5],s=1)
plt.xlabel('f1')
plt.ylabel('f2')
plt.show()

#%%
n_iterations = 0        
for res in results:
   n_iterations += res[3]
print("n_iterations = {}".format(n_iterations))
