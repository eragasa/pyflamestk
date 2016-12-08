import pyflamestk.eam
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

eam_filename = 'Ni99.eam.alloy'

def func_pairr_effective_charge(r, a1, a2, rc):
    # Ref: Foiles PhysRevB 32 3409 (1985)
    # "Application of the Embedded atom method to liquid transition method
    print(a1,a2,rc)
    myr = r.tolist()
    phi = []
    for r in myr:
        if r < rc:
            z_sq = a1 * (rc-r)**3.0 + a2 * (rc- r)**4
            phi.append(z_sq)
        else:
            phi.append(0)
    return np.array(phi)
#    return np.array(A * q**2./r)

def func_embedding_universal(rho, F0, p, q, F1):
  # Reference:  Johnson and Oh.  J. Mat. Resrch. 1989
  # Analytic embedded method model for bcc metals
  #print('F0',F0,'p',p,'q',q,'F1',F1)
  #print("{:18.8f} {:18.8f} {:18.8f} {:18.8f} {:18.8f}".format(rho,F0,p,q,F1))
  print(F0, p, q, F1)
  return F

def func_rho_3s_metal(r, a0, a1, a2, c):
  # Reference: P. Mitev et al, Sim in Mat Sci and Engineering, 2006
  # Embedded atom method potentials employing a faithful density representation
  # Note: this functional is based on the 3S-electronic wave function
  print(a0,a1,a2,c)
  rho = (a0 + a1 * r + a2 * r**2)**2 * np.exp(-r/c)
  return rho
  
  
eam_file = pyflamestk.eam.EamFileReader(eam_filename)
eam_file.read()

plt.figure(1)
plt.plot(eam_file.fembed['Ni'][:,0],eam_file.fembed['Ni'][:,1])
plt.show()

plt.figure(2)
plt.plot(eam_file.fdens['Ni'][:,0],eam_file.fdens['Ni'][:,1])
plt.show()

plt.figure(3)
plt.plot(eam_file.effpot['NiNi'][:,0],eam_file.effpot['NiNi'][:,1])
plt.show()

curve_fitter = pyflamestk.eam.EamCurveFitter(element_names = ['Ni'])





effpot_r   = eam_file.effpot['NiNi'][:,0]
effpot     = eam_file.effpot['NiNi'][:,1]
p0 = [0.070937, 0.146031, 3.0045]
opt_param_pair,  opt_cov_pair  = scipy.optimize.curve_fit(func_pair_effective_charge,
                                                          effpot_r.tolist(), 
                                                          effpot.tolist(),
                                                          p0 = p0,
                                                          maxfev=200000)
fembed_rho = eam_file.fembed['Ni'][:,0]
fembed     = eam_file.fembed['Ni'][:,1]
p0         = [1,2,1,1]
opt_param_fembed,  opt_cov_fembed  = scipy.optimize.curve_fit(func_embedding_universal,
                                                          fembed_rho.tolist(), 
                                                          fembed.tolist(),
                                                          p0 = p0,
                                                          maxfev=20000)
fdens_r = eam_file.fdens['Ni'][:,0]
fdens   = eam_file.fdens['Ni'][:,1]
opt_param_fdens,  opt_cov_fdens  = scipy.optimize.curve_fit(func_rho_3s_metal,
                                                          fdens_r.tolist(), 
                                                          fdens.tolist(),
                                                          maxfev=20000)


phi_fit = func_pair_effective_charge(eam_file.effpot['NiNi'][:,0],*opt_param_pair)
fembed_fit = func_embedding_universal(eam_file.fembed['Ni'][:,0],*opt_param_fembed)
fdens_fit = func_rho_3s_metal(eam_file.fdens['Ni'][:,0],*opt_param_fdens)
plt.figure(1)
plt.plot(eam_file.fembed['Ni'][:,0],eam_file.fembed['Ni'][:,1])
plt.plot(eam_file.fembed['Ni'][:,0],fembed_fit)
plt.show()

plt.figure(2)
plt.plot(eam_file.fdens['Ni'][:,0],eam_file.fdens['Ni'][:,1])
plt.plot(eam_file.fdens['Ni'][:,0],fdens_fit)
plt.show()

plt.figure(3)
plt.plot(eam_file.effpot['NiNi'][:,0],eam_file.effpot['NiNi'][:,1])
plt.plot(eam_file.effpot['NiNi'][:,0],phi_fit)

plt.show()

print("optimal_pair:",opt_param_pair)
print("optimal_fembed:",opt_param_fembed)
print("opt_param_fdens:",opt_param_fdens)
