from Modèle import *
from random import *

#fonction donnant la température et la pression en fonction de l'altitude    
ISA_temp = interpolate.interp1d([0,11000,20000,32000],[288,216.5,216.5,228.5])
ISA_P = interpolate.interp1d([0,11000,20000,32000],[101325,22632,5474.9,868.02])


def monte_carlo(N):
    """Prend en entree le nombre d'essais N
    Retourne une liste contenant N listes de paramètres d'entrees
    Chaque liste correspond a une experience, les parametres etant choisis au 
    hasard"""
    
    resultat = open("resulat.txt","a")#cré fichier resultat
    for i in range (N):
        experience=[]
        
        #Altitude
        z=10000#randint(0,32000)
        
        TC = ISA_temp(z)
        PC= ISA_P(z)
        
        #Nombre de Mach
        M=0.9#choice([0,0.1,0.3,0.5,0.7,0.8,0.85,0.9,0.95,0.99])#choix au hasard
        #Temperature totale
        TT=TC*(1+0.2*M**2)
        #Pression totale
        PT=PC*(1+0.2*M**2)**3.5
        #Taux de compression
        ts=1.6#randint(1,3)#soufflante
        tcbp=3#randint(1,3)#compresseur bp
        tchp=3#randint(1,3)#compresseur hp
        #Rendement
        rs=0.7#uniform(0.5,1)
        rcbp=0.7#uniform(0.5,1)
        rchp=0.7#uniform(0.5,1)
        rtbp=0.7#uniform(0.5,1)
        rthp=0.7#uniform(0.5,1)
        #Alpha
        alpha=uniform(0.06,0.06)
        tt=1.5#turbine
        #Coefficient de partage du flux
        lamb=randint(1,10)
        #Flux
        WA=1
        WF=alpha*WA
        #Vitesse
        VA=M*1224#(1.4*237*TC)**0.5      
        
        experience=[TT,PT,ts,tcbp,tchp,tt,rs,rcbp,rchp,rtbp,rthp,alpha,lamb,WA,WF,VA,turboreacteur(TT,PT,ts,tcbp,tchp,tt,rs,rcbp,rchp,rtbp,rthp,alpha,lamb,WA,WF,VA)]
        resultat.write(str(experience)+'\n')#ajoute au fichier
    resultat.close()
