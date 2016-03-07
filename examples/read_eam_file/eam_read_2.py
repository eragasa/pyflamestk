import pyflamestk.eam
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------
# pairr functions
# pairr are pair potentials multiplied by r
# pairr(r) = pair(r) * r
#------------------------------------------------------------------------------

def func_pairr_effective_charge(r, a1, a2, rc):
    # Ref: Foiles PhysRevB 32 3409 (1985)
    # "Application of the Embedded atom method to liquid transition method
    myr = r.tolist()
    phi = []
    for r in myr:
        if r < rc:
            z_sq = a1 * (rc-r)**3.0 + a2 * (rc- r)**4
            phi.append(z_sq)
        else:
            phi.append(0)
    return np.array(phi)

#------------------------------------------------------------------------------
# embedding functionls
#------------------------------------------------------------------------------

def func_embedding_universal(rho, F0, p, q, F1):
  # Reference:  Johnson and Oh.  J. Mat. Resrch. 1989
  # Analytic embedded method model for bcc metals
  #print('F0',F0,'p',p:w,'q',q,'F1',F1)
  #print("{:18.8f} {:18.8f} {:18.8f} {:18.8f} {:18.8f}".format(rho,F0,p,q,F1))
  vals = F0 *(q/(q-p) * (rho+1.e-30)**p - p/(q-p) * (rho+1.e-30) ** q) + F1 * rho
  return vals
  
def func_embedding_bjs(rho,F0,gamma,F1):
    ''' implements BJS functional form of embedding function '''
    F0,gamma,F1 = params
    embedvals = F0 * ( 1. - gamma*log(r+1.e-30) ) * r**gamma + F1 * r
    return embedvals

def func_embedding_foiles(rho, Ecoh,am,r0,rho0,lambda0,lambdarose,De,rp,rd):
    ''' implements Foiles-style embeding function
        (i.e. density and pair potential forced to match Rose EOS)
        parameter list:
        p[0]   p[1]     p[2] p[3] p[4]    p[5]       p[6]   p[7] p[8]
        E_coh  a(morse) r0  rho0 lambda0 lambdarose De      rp   rd  '''
    Ecoh,am,r0,rho0,lambda0,lambdarose,De,rp,rd = params
    embedvals = empty_like(rho)
    k=0
    for rhostar in rho:
        #solve density for lattice constant (a) where density (rhostar) is found
        rhop = (rho0,r0,lambda0,rd,rhostar)
        a = brentq(rhofxn,0.,10000.,rhop,xtol=1.0e-8) #lattice constant where rhostar is found
        #find E_Rose for lattice constant
        astar = (a-r0*sqrt(2.)) / (lambdarose*sqrt(2.)*r0)
        Erose = Ecoh*(1+astar)*exp(-astar)
        #find pair potential for lattice constant
        pp = (De,am,r0,rp)
        Epot = 12.*fpair(a/sqrt(2.),pp)+6.*fpair(a,pp)+24.*fpair(sqrt(1.5)*a,pp)+12.*fpair(sqrt(2.)*a,pp)+24.*fpair(sqrt(2.5)*a,pp)+8.*fpair(sqrt(3.)*a,pp)
        #calculate value of embedding fxn
        embedvals[k] = Erose - 0.5*Epot
        k += 1

#------------------------------------------------------------------------------
# election density functions
#------------------------------------------------------------------------------
def func_edensity_3s_metal(r, a0, a1, a2, c):
  # Reference: P. Mitev et al, Sim in Mat Sci and Engineering, 2006
  # Embedded atom method potentials employing a faithful density representation
  # Note: this functional is based on the 3S-electronic wave function
  rho = (a0 + a1 * r + a2 * r**2)**2 * np.exp(-r/c)
  return rho

def func_edensity_exponential(r, lambda0, r0, rd):
    ''' implements Exponent ial e- density functional form '''
    rho0,lambda0,r0,rd = params
    densvals = rho0 * exp( -1.*( r - r0 )/lambda0 ) * psi( (r-rd)/(globalcutoff-rd))
    return densvals

def func_edensity_exponential2(a,rho0,r0,lambda0,rd,rhostar):
    ''' calculates ideal e- density based on exponential functional form
        data input format:  rho0  r0    lambda0  rd   rhostar '''
    return rho0*(12.*exp(-(a/sqrt(2.)-r0)/lambda0) * psi( (a/sqrt(2.)-rd) / (globalcutoff-rd) )
            + 6.*exp(-(a-r0)/lambda0) * psi( (a-rd) / (globalcutoff-rd) )
            + 24.*exp(-(a*sqrt(1.5)-r0)/lambda0) * psi( (a*sqrt(1.5)-rd) / (globalcutoff-rd) )
            + 12.*exp(-(a*sqrt(2.)-r0)/lambda0) * psi( (a*sqrt(2.)-rd) / (globalcutoff-rd) )
            + 24.*exp(-(a*sqrt(2.5)-r0)/lambda0) * psi( (a*sqrt(2.5)-rd) / (globalcutoff-rd) )
            + 8.*exp(-(a*sqrt(3.)-r0)/lambda0) * psi( (a*sqrt(3.)-rd) / (globalcutoff-rd) )
            ) - rhostar

#------------------------------------------------------------------------------
# cutoff functions
#------------------------------------------------------------------------------
  
eam_filename1 = 'Ni99.eam.alloy'
eam_filename2 = 'Ragasa.eam.alloy'
  
eam_file1 = pyflamestk.eam.EamFileReader(eam_filename1)
eam_file1.read()

curve_fitter = pyflamestk.eam.EamCurveFitter(element_names = ['Ni'])

Nrho = eam_file1.Nrho
Nr   = eam_file1.Nr
drho = eam_file1.drho
dr   = eam_file1.dr

effpot_r   = eam_file1.effpot['NiNi'][:,0]
effpot     = eam_file1.effpot['NiNi'][:,1]
p0 = [0.070937, 0.146031, 3.0045]

opt_param_pair,  opt_cov_pair  = scipy.optimize.curve_fit(func_pairr_effective_charge,
                                                          effpot_r.tolist(), 
                                                          effpot.tolist(),
                                                          p0 = p0,
                                                          maxfev=200000)
fembed_rho = eam_file1.fembed['Ni'][:,0]
fembed     = eam_file1.fembed['Ni'][:,1]
p0         = [1,2,1,1]
opt_param_fembed,  opt_cov_fembed  = scipy.optimize.curve_fit(func_embedding_universal,
                                                          fembed_rho.tolist(), 
                                                          fembed.tolist(),
                                                          p0 = p0,
                                                          maxfev=20000)
fdens_r = eam_file1.fdens['Ni'][:,0]
fdens   = eam_file1.fdens['Ni'][:,1]
opt_param_fdens,  opt_cov_fdens  = scipy.optimize.curve_fit(func_edensity_3s_metal,
                                                          fdens_r.tolist(), 
                                                          fdens.tolist(),
                                                          maxfev=20000)

fpair_fit  = func_pairr_effective_charge(eam_file1.effpot['NiNi'][:,0],*opt_param_pair)
fembed_fit = func_embedding_universal(eam_file1.fembed['Ni'][:,0],*opt_param_fembed)
fdens_fit  = func_edensity_3s_metal(eam_file1.fdens['Ni'][:,0],*opt_param_fdens)

eam_file2  = pyflamestk.eam.EamFileWriter(title    = "effective_charge/5smetal/universal",
                                          elements = ['Ni'],
                                          filename = eam_filename2,
                                          Nrho   = Nrho,
                                          drho   = drho,
                                          Nr     = Nr,
                                          dr     = dr)

eam_file2.compute_embedding(func_embedding_universal,
                            opt_param_fembed)
eam_file2.compute_density(func_edensity_3s_metal,
                          opt_param_fdens)
eam_file2.compute_pairr(func_pairr_effective_charge,
                       opt_param_pair)                                          
eam_file2.write()

eam_file3 = pyflamestk.eam.EamFileReader(eam_filename2)                                          
eam_file3.read()

plt.figure(1)
plt.plot(eam_file1.fembed['Ni'][:,0],  eam_file1.fembed['Ni'][:,1])
plt.plot(eam_file2.rho ,eam_file2.elements.embed)
plt.plot(eam_file3.fembed['Ni'][:,0],  eam_file3.fembed['Ni'][:,1])
#plt.plot(eam_file1.fembed['Ni'][:,0],  fembed_fit)
plt.show()

plt.figure(2)
plt.plot(eam_file1.fdens['Ni'][:,0], eam_file1.fdens['Ni'][:,1])
plt.plot(eam_file3.fdens['Ni'][:,0], eam_file3.fdens['Ni'][:,1])
#plt.plot(eam_file1.fdens['Ni'][:,0],  fdens_fit)
plt.show()

plt.figure(3)
plt.plot(eam_file1.effpot['NiNi'][:,0],eam_file1.effpot['NiNi'][:,1])
plt.plot(eam_file3.effpot['NiNi'][:,0],eam_file3.effpot['NiNi'][:,1])
#plt.plot(eam_file1.effpot['NiNi'][:,0],fpair_fit)

plt.show()

print("optimal_pair:",opt_param_pair)
print("optimal_fembed:",opt_param_fembed)
print("opt_param_fdens:",opt_param_fdens)
