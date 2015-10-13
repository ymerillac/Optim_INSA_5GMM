"""NACA contains a function to calculate a NACA profile"""

import matplotlib.pyplot as plt
import sys
from numpy import *

def create(m,p,t,C=1,n=150):
    """ NACA calculate the NACA profile of an airfoil where the inputs are:
            - m is the maximum camber as percentage of the chord
            - p is the location of maximum camber from the airfoil leading edge in tens of percents of the chord
            - t describing maximum thickness of the airfoil as percent of the chord 
            - C is the chord length
            - n is the number of points"""
    try:
        if not ( 0<n and m<1 and p<1 and t<1 and n<=173): #and 0<=m and 0<=p and 0<=t
            raise ValueError
    except ValueError:
        sys.exit("In NACA.create, all inputs must be positive and (m,p,t) must be less than 1. Number of point (n)  must be less than 173")
        

    p=abs(p)
    t=abs(t)
    
    x=linspace(0,C,n)

    Yc=zeros(n)
    Teta=zeros(n)

    Yt=(t*C/0.2)*(0.2969*( (1./C)*x )**(0.5)-0.1260*( (1./C)*x ) -0.3516*( (1./C)*x )**2 +0.2843*( (1./C)*x )**3 -0.1015*( (1./C)*x )**4)

    if p>0:
        k=0
        while x[k]<p*C:
            k=k+1
            
        Yc[0:k]=m*x[0:k]/(p**2)*(2*p - x[0:k]/C)
        Yc[k:n]=m*(C-x[k:n])/((1-p)**2)*(1-2*p+x[k:n]/C)

        Teta[0:k]=arctan( 2*m/(p**2)*(p -x[0:k]/C) )
        Teta[k:n]=arctan( 2*m/((1-p)**2)*(p -x[k:n]/C) )
    

    Yu=Yc+Yt*cos(Teta)
    Yl=Yc-Yt*cos(Teta)

    Xu=x-Yt*sin(Teta)
    Xl=x-Yt*sin(Teta)


    return Xu,Xl,Yu,Yl


# tests:
if __name__=="__main__":
    Xu,Xl,Yu,Yl=create(m=0.08,p=0.4,t=0.15,C=1,n=150)
    plt.plot(Xu,Yu,"b",Xl,Yl,"b")
    plt.ylim(-0.5,0.5)
    plt.xlim(-0.,1.)
    plt.show()
