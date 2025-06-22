import scipy.optimize as optimize

class SwaptionType(enum.Enum):
    RECIEVER = 1.0
    PAYER = -1.0

def PaymentDates(Tm,Tn,f):
    tau = 1/f
    n = 0
    n = (Tn - Tm)/tau + 1
    ti_grid = np.linspace(Tm,Tn,int(n))
    return ti_grid

def SpotSwap(lambd,eta,P0T,Tm,Tn,K,f=4, rStar = 0.03):
    #rStar is the one for Jamishidian’s trick
    #rStar = 0.03 is dummy value

    ti_grid=[]
    ti_grid = PaymentDates(Tm,Tn,f)
    temp=np.zeros([len(ti_grid), 1])

    for (idx,ti) in enumerate(ti_grid):        
        if ((ti>Tm) and (ti<Tn)):
            temp[idx] = K * (1/f) * HW_ZCB(lambd,eta,P0T,Tm,ti,rStar)
        elif (ti == Tn):
            temp[idx] = (1+ K/f) * HW_ZCB(lambd,eta,P0T,Tm,ti,rStar)
    
    swap = 1 - np.sum(temp)
    return swap

def findrStar(lambd,eta,P0T,Tm,Tn,K,f=4):
    A = lambda x: SpotSwap(lambd,eta,P0T,Tm,Tn,K,f,x)
    result = optimize.newton(A,0.01)
    return result

def SwaptionHW1F_Pricing(SwaptionType,lambd,eta,P0T,Tm,Tn,K, f=4):
    ti_grid=[]
    ti_grid = PaymentDates(Tm,Tn,f)
    temp=np.zeros([len(ti_grid), 1])

    rStar = findrStar(lambd,eta,P0T,Tm,Tn,K)
    annuity = 0

    for (idx,ti) in enumerate(ti_grid):     

        K_i =  HW_ZCB(lambd,eta,P0T,Tm,ti,rStar)
        # Future Annuity, Price at 0, Start at Tm+1, Last payment at Tn
        annuity = annuity + (1/f) * P0T(ti) 

        if ((ti>Tm) and (ti<Tn)):
            temp[idx] = K * (1/f) * HW_ZCB_CallPutPrice(OptionType.PUT,K_i,lambd,eta,P0T,Tm,ti)

        elif (ti == Tn):
            temp[idx] = (1+ K/f) * HW_ZCB_CallPutPrice(OptionType.PUT,K_i,lambd,eta,P0T,Tm,ti)

        value = np.sum(temp)


    if SwaptionType == SwaptionType.PAYER:
        return value

    elif SwaptionType == SwaptionType.RECIEVER:
        return value + K * annuity + P0T(Tn) - P0T(Tm)
