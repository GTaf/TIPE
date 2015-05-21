from calcul_T import *
import scipy
from math import *
from scipy import interpolate
import numpy as np


def gamma_temperature(fluide,alpha):
    """ Prend en entrée la nature du fluide : 0=air, 1=kerosene, 2= melange air-kerosene,
    la richesse alpha (valant 0 pour l'air et 1 pour le kerosene pur),
    la vitesse du son a.
    Renvoi gamma(T)"""
    a = lambda T: exp(3090/T)
    if fluide==0:
        Cpr=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a(T))/(T*T*(a(T)-1)*(a(T)-1))
        gamma=lambda T: Cpr(T)/(Cpr(T)-1)
    if fluide==1:
        Cpr=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
        gamma=lambda T: Cpr(T)/(Cpr(T)-1)
    if fluide==2:
        Cpra=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a(T))/(T*T*(a(T)-1)*(a(T)-1))
        Cprk=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
        Cpr=lambda T: (Cpra(T)+alpha*Cprk(T))/(alpha+1)
        gamma=lambda T: Cpr(T)/(Cpr(T)-1)
    return gamma
    
def Cp_temperature(fluide,alpha):
    """ Prend en entrée la nature du fluide : 0=air, 1=kerosene, 2= melange air-kerosene,
    la richesse alpha (valant 0 pour l'air et 1 pour le kerosene pur),
    la vitesse du son a. RIEN A VOIR
    Renvoi Cp(T)"""
    a = lambda T: exp(3090/T)
    if fluide==0:
        Cpr=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a(T))/(T*T*(a(T)-1)*(a(T)-1))
    if fluide==1:
        Cpr=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
    if fluide==2:
        Cpra=lambda T: 3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a(T))/(T*T*(a(T)-1)*(a(T)-1))
        Cprk=lambda T: -0.0000018373*T*T+0.00801994*T+4.47659
        Cpr=lambda T: (Cpra(T)+alpha*Cprk(T))/(alpha+1)
    return Cpr
#calcul des fonctions utilisés dans les blocs du turboréacteur



a = lambda T: exp(3090/T)

"""integrale de gamma sur gamma-1, pour l'air pur"""
fp1=lambda T: gamma_temperature(0,0)(T)/(T*(gamma_temperature(0,0)(T)-1)) #fonction gamma/T*(gamma-1) de l'air
F1=calcul_T(20,10000,1,fp1) #calcul de la primitive de gamma/T*(gamma-1) pour l'air
    
    
f1 = interpolate.interp1d(F1[0],F1[1])#fonction primitive
finv1=interpolate.interp1d(F1[1],F1[0]) #reciproque

""" integrale des cp"""
Cpra=lambda T: (3.5-0.000028*T+0.0000000224*T*T+(3090*3090*a(T))/(T*T*(a(T)-1)*(a(T)-1)))*298 #air
Cprk=lambda T: (-0.0000018373*T*T+0.00801994*T+4.47659)*298 #kerosene
F2a = calcul_T(10,10000,1,Cpra)#passage long
F2k = calcul_T(10,10000,1,Cprk)

f2a = interpolate.interp1d(F2a[0],F2a[1])#fonction primitive
finv2a=interpolate.interp1d(F2a[1],F2a[0]) #fonction l'inverse

f2k = interpolate.interp1d(F2k[0],F2k[1])#fonction primitive
finv2k=interpolate.interp1d(F2k[1],F2k[0]) #fonction l'inverse


    
def turboreacteur(T1,P1,ts,tcbp,tchp,tt,rs,rcbp,rchp,rtbp,rthp,alpha,lamb,WA,WF,VA):
    """Modélisation d'un turboréacteur
    Hypothèses : transformations isentropiques dans com  """          
    
    f2 = lambda T : (f2a(T)+alpha*f2k(T))/(1+alpha)
    F2 = (np.array(F2a)+alpha*np.array(F2k))/(1+alpha)
    finv2=interpolate.interp1d(F2[1],F2[0])
            
    print("T1 = ",T1)
    #Soufflante, obtenu par integration(Laplace adapté)
    T2=finv1(f1(T1)+log(ts)) #log correspond au logarithme néperien
    print("T2 = ",T2)
    #Compresseur BP, obtenu par integration(Laplace adapté)
    T3=finv1(f1(T2)+log(tcbp))#log correspond au logarithme néperien
    print("T3 = ",T3)
        
    #Compresseur HP, obtenu par integration(Laplace adapté)
    T4=finv1((f1(T3)+log(tchp)))#log correspond au logarithme néperien
        
    #Chambre de combustion, 1er principe thermochimie
    DfCO2=394000#DfCO2
    DfH2O=280000#DfH20
    DfCH4=88000
    avanct=WF/(0.012+4*0.001)
    Df=avanct*(DfCO2+2*DfH2O-DfCH4)
    
    print("T4 = ",T4)
    
    T5=finv2(f2(T4)+(Df/(WA+WF)))
    
    print("Température chambre = ",T5)
    
    #Turbine HP, obtenu par equilibre HP
    T6=finv2(f2(T5)-(f2(T4)-f2(T3))*WA/(WA+WF)*rthp*rchp) #prise en compte rendement turbine+compresseur
    print("T6 = ",T6)
    
    #Turbine BP, obtenu par equilibre BP
    T7=finv2(f2(T6)-(f2(T3)-f2(T1))*WA/(WA+WF)*rtbp*rcbp-(f2(T2)-f2(T1))*lamb*WA/(WA+WF)*rtbp*rs) #prise en compte turbine+compresseur+souflante
    print("T7 = ",T7)
    
    #Mélangeur, application 1er principe
    T8=finv2(((WA+WF)*f2(T7)+lamb*WA*f2(T2))/(WF+WA+WA*lamb))
    print("T8 = ",T8)    
    
    #Tuyère, application 1er principe
    T9=finv1(f1(T8)+log(tt))
    C9=sqrt(2*(-f2(T9)+f2(T8)))
    
    #Rendement
    Pcin=(((1+lamb)*WA+WF)*C9**2-(1+lamb)*WA*VA**2)/2
    Pth=(WA+WF)*(f2(T5)-f2(T4))
    Rendement=Pcin/Pth
    
    return Rendement,T5
    
    
