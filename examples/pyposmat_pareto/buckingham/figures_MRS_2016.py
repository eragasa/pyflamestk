# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 12:44:54 2016

@author: Eugene
"""
import pyflamestk.pyposmat
import pyflamestk.pareto
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#------------------------------------------------------------------------------
# figure 1

x_label = 'MgO_NaCl_latt_abserr'
y_label = 'MgO_NaCl_c11_abserr'

axes    = [0., 0.04, 0, 30]

x_idx   = all_names.index(x_label)
y_idx   = all_names.index(y_label)
x_ps_idx = qois.index(x_label)
y_ps_idx = qois.index(y_label)
pareto_front = pareto_frontier_2d(np_all_sims[:,x_idx],
                                  np_all_sims[:,y_idx],
                                  maxX = False, maxY= False)
plt.figure()
plt.scatter(np_all_sims[:,x_idx],
            np_all_sims[:,y_idx])
plt.scatter(pareto_set[:,x_ps_idx],
            pareto_set[:,y_ps_idx],color='y')
plt.plot(pareto_front[0],pareto_front[1],)
plt.axis(axes)
#plt.axis([min(pareto_front[0]),
#          max(pareto_front[0]),
#          min(pareto_front[1]),
#          max(pareto_front[1])])
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.show()