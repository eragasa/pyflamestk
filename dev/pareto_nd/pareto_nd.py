import numpy as np

def pareto_frontier_multi(myArray):
    # Sort on first dimension
    myArray = myArray[myArray[:,0].argsort()]
    # Add first row to pareto_frontier
    pareto_frontier = myArray[0:1,:]
    # Test next row against the last row in pareto_frontier
    for row in myArray[1:,:]:
        if sum([row[x] >= pareto_frontier[-1][x]
                for x in range(len(row))]) == len(row):
            # If it is better on all features add the row to pareto_frontier
            pareto_frontier = np.concatenate((pareto_frontier, [row]))
    return pareto_frontier

def create_pareto_set(dataset):
    pareto_set = []
    for point in dataset:
        is_pareto_point = False
        if len(pareto_set)==0:
            pareto_set.append(point)

        for pareto_point in pareto_set:
            # check if new point is a pareto point
            for i,val in enumerate(pareto_point):
                point[i] > val
                is_pareto_point = True
            # check if pareto point is dominated
        if is_pareto_point:
           pareto_set.append(point)

        # check if any existing pareto point is dominated
                        
    return pareto_set
    
pareto_set = []
def update(p):
  
  if any(q > p for q in pareto_set):
    return
  for q in [q for q in pareto_set if p > q]:
    pareto_set.remove(q)
  pareto_set.append(p)
  print("p:",pareto_set)
  
my_array_1 = np.random.normal(0,1,size=[1,3])
n_simulations = 10000
for i in range(n_simulations):
  my_array_1 = np.concatenate((my_array_1,
                               np.random.normal(0,1,size=[1,3])))
pareto_set = create_pareto_set(dataset)
#pareto_frontier = pareto_frontier_multi(my_array_1)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(my_array_1[:,0], my_array_1[:,1], my_array_1[:,2])
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(pareto_set[:,0],
            pareto_set[:,1],
            pareto_set[:,2])
plt.show()