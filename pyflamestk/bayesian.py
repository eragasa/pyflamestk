
class BayesianEstimator:
    
    def __init__(self):
        self._theta
        self._q
        self._rho_t = None          # prior distribution, \rho(t)
        self._rho_q                 # marginal likelyhood
        self._rho_q_t               # likelihood
        self._rho_t_q               # posterior distribution
        self._marginal_likelihood
        
    @property
    def posterior_distribution(self): return self._rho_t_q
    
    @posterior_distribution.setter
    def posterior_distribution(self, rho_t_q):
        self._rho_t_q = rho_t_q
    
    @property
    def prior_distribution(self): return self._rho_t
    
    @prior_distribution.setter
    def prior_distribution(self, rho_t):
        self._rho_t = rho_t 
        
    @property
    def marginal_likelihood(self): return self._rho_q
        
    @marginal_likelihood.setter
    def marginal_likelihood(self, rho_q):
        self._rho_q = rho_q
        
    @property
    def likelihood_function(self): return self._rho_q_t
        
    @likelihood_function.setter
    def likelihood_function(self, rho_q_t):
        self._rho_q_t = rho_q_t
        
    
    
