"""

Original from "Financial Engineering" course by L.A. Grzelak
Modified the code by using without integration approach, as mentioned in the course

"""
## HW Start
def f0T(t,P0T):
    dt = 0.0001    
    expr = - (np.log(P0T(t+dt))-np.log(P0T(t-dt)))/(2*dt)
    return expr

def HW_theta(lambd,eta,P0T):
    dt = 0.001    
    theta = lambda t: 1.0/lambd * (f0T(t+dt,P0T)-f0T(t-dt,P0T))/(2.0*dt) + f0T(t,P0T) + eta*eta/(2.0*lambd*lambd)*(1.0-np.exp(-2.0*lambd*t))
    return theta

def HW_A_Integrate(lambd,eta,P0T,T1,T2): ### Original Code
    tau = T2-T1
    zGrid = np.linspace(0.0,tau,250)
    B_r = lambda tau: 1.0/lambd * (np.exp(-lambd *tau)-1.0)
    theta = HW_theta(lambd,eta,P0T)    
    temp1 = lambd * integrate.trapz(theta(T2-zGrid)*B_r(zGrid),zGrid)  
    temp2 = eta*eta/(4.0*np.power(lambd,3.0)) * (np.exp(-2.0*lambd*tau)*(4*np.exp(lambd*tau)-1.0) -3.0) + eta*eta*tau/(2.0*lambd*lambd)
    return temp1 + temp2

def HW_A_Closed(lambd,eta,P0T,T1,T2): ### Derived and consistent with wiki Closed form, Checked same
    return -f0T(T1,P0T)* HW_B(lambd,eta,T1,T2) +np.log(P0T(T2))-np.log(P0T(T1)) - eta*eta/(4.0*lambd) *HW_B(lambd,eta,T1,T2) *HW_B(lambd,eta,T1,T2) *(1- np.exp(-2*lambd * T1))

def HW_A(lambd,eta,P0T,T1,T2):
    return -f0T(T1,P0T)* HW_B(lambd,eta,T1,T2) +np.log(P0T(T2))-np.log(P0T(T1)) - eta*eta/(4.0*lambd) *HW_B(lambd,eta,T1,T2) *HW_B(lambd,eta,T1,T2) *(1- np.exp(-2*lambd * T1))

def HW_B(lambd,eta,T1,T2):
    return 1.0/lambd *(np.exp(-lambd*(T2-T1))-1.0)

def HW_ZCB(lambd,eta,P0T,T1,T2,rT1):
    B_r = HW_B(lambd,eta,T1,T2)
    A_r = HW_A(lambd,eta,P0T,T1,T2)
    return np.exp(A_r + B_r *rT1)

##HW End
