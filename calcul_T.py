import scipy.integrate as sp
import math


def calcul_T(T0,T1,p,Cp):
    """calcul de la fonction primitive de Cp(T), discrétisée avec un pas de 1/p
    sur l'intervalle [T0,T1], entiers, p entier"""

    R1=[T0]
    R2=[sp.quad(Cp,T0,T0+(1/p))[0]]#liste vide contenat les futurs couple (T,F(T))
    for i in range(1,(T1-T0)*p):
        R1.append(T0+i*(1/p))#incrémente T
        R2.append(R2[-1]+sp.quad(Cp,T0+i*(1/p),T0+(i+1)*(1/p))[0])#camcul intég
    #R1.pop(0)#on enlève les premiers termes nuls
    #R2.pop(0)
    return [R1,R2]#R1 température, R2 valeur de l'intégrale 
    
    
def f(T):# à definir plus précisément, fonction Cp(T)
    a = lambda T: exp(3090/T)    
    
    return 3.5-0.000028*T+0.0000000224*(T**2)+(3000*3000*a(T))/(T*(a(T)-1)*(a(T)-1))

import matplotlib.pyplot as pypl
def trace(T0,T1,p,Cp):
    """trace l'intégrale des CP(T)dT"""
    R1=calcul_T(T0,T1,p,Cp)[0]
    R2=calcul_T(T0,T1,p,Cp)[1]
    pypl.plot(R1,R2)
    pypl.show()
    return R2[-1],R2[0]
    
