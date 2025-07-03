import numpy as np
import enum 
import scipy.stats as st

def Annuity(Tm,Tn,f,P0T):
    n = (Tn - Tm) * f + 1
    ti_grid = np.linspace(Tm,Tn,int(n))
    tau = ti_grid[1]- ti_grid[0]
    temp =0
    for (idx,ti) in enumerate(ti_grid):
        if ti>Tm:
            temp = temp + tau * P0T(ti)
    
    return temp

class OptionType(enum.Enum):
    RECEIVER = 1.0
    PAYER = -1.0

def BachelierSwaption(CP, S,K,sigma,Tm,Tn,f,Annuity,P0T):
    A = Annuity(Tm,Tn,f,P0T)
    d = (S-K) / (sigma*np.sqrt(Tm))

    value = 0.0

    if CP==OptionType.PAYER:
        value = A * ((S-K) * st.norm.cdf(d) + sigma * np.sqrt(Tm) * st.norm.pdf(d) )
    elif CP==OptionType.RECEIVER:
        value = A * ((K-S) * st.norm.cdf(-d) + sigma * np.sqrt(Tm) * st.norm.pdf(-d) )

    return value

