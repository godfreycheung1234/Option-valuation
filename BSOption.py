import numpy as np
import scipy.stats as st

class BSOption:
    def __init__(self, option_type, S, K, T, r, sigma):

        # Convert option_type to lowercase to make it case insensitive
        self.option_type = option_type.lower()
        self.S = S
        self.K = np.array(K).reshape([len(K),1])
        self.T = T
        self.r = r
        self.sigma = sigma
        
        if self.option_type not in ['call', 'put']:
            raise ValueError("Option type must be 'call' or 'put'.")
            
    def d1(self):
        return (np.log(self.S / self.K) + (self.r + 0.5 * np.power(self.sigma,2.0)) * (self.T)) / (self.sigma * np.sqrt(self.T))
    
    def d2(self):
        return self.d1() - self.sigma * np.sqrt(T)
        
    def price(self):

        if self.option_type == 'call':
            return (self.S * st.norm.cdf(self.d1()) - 
                    self.K * np.exp(-self.r * self.T) * st.norm.cdf(self.d2()))

        elif self.option_type == 'put':
            return (self.K * np.exp(-self.r * self.T) * st.norm.cdf(-self.d2()) - 
                    self.S * st.norm.cdf(-self.d1()))

    def delta(self):
        if self.option_type == 'call':
            return st.norm.cdf(self.d1())
            
        elif self.option_type == 'put':
            return (st.norm.cdf(self.d1())-1.0)
            
    def gamma(self):
        return st.norm.pdf(self.d1()) / (self.S * self.sigma * np.sqrt(self.T))

    def vega(self):
        return self.S * st.norm.pdf(self.d1()) * np.sqrt(self.T)

    def theta(self):

        if self.option_type == 'call':
            return (-self.S * st.norm.pdf(self.d1()) * self.sigma / (2 * np.sqrt(self.T)) -
                    self.r * self.K * np.exp(-self.r * self.T) * st.norm.cdf(self.d2()))

        elif self.option_type == 'put':
            return (-self.S * st.norm.pdf(self.d1()) * self.sigma / (2 * np.sqrt(self.T)) +
                    self.r * self.K * np.exp(-self.r * self.T) * st.norm.cdf(-self.d2()))

    def rho(self):

        if self.option_type == 'call':
            return self.K * self.T * np.exp(-self.r * self.T) * st.norm.cdf(self.d2())

        elif self.option_type == 'put':
            return -self.K * self.T * np.exp(-self.r * self.T) * st.norm.cdf(-self.d2())
