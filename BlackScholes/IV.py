import scipy.optimize as optimize

def IV_BlackSwpation(CP,marketPrice,S,K,Tm,Tn,f,Annuity,P0T):
    sigmaGrid = np.linspace(0.0,0.5,1000)
    optPriceGrid = BlackSwaption(CP,S,K,sigmaGrid,Tm,Tn,f,Annuity,P0T)
    sigmaInitial = np.interp(marketPrice,optPriceGrid,sigmaGrid)

    func = lambda sigma: np.power(BlackSwaption(CP,S,K,sigma,Tm,Tn,f,Annuity,P0T) - marketPrice, 1.0)
    impliedVol = optimize.newton(func, sigmaInitial, tol=1e-15)

    return impliedVol
