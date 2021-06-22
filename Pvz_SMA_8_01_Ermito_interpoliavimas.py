
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
import tkinter as tk
from tkinter import * 
from PyFunkcijos import *
import math
import time
import sympy


# ---------- Ermito daugianariai -----------
def Hermite(X,j,x):
    N=size(X);
    L=Lagrange(X,j,x); 
    DL=D_Lagrange(X,j,X[j])
    U=(1-2*DL*(x-X[j]))*np.square(L)
    V=(x-X[j])*np.square(L);
    return U,V

# --------- Lagranzo daugianario isvestine pagal x --------------
def D_Lagrange(X,j,x):
    N=size(X)
    DL=np.zeros(x.shape, dtype=np.double); # DL israiskos skaitiklis
    for i in range(0,N): # ciklas per atmetamus narius
        if (i == j) : continue
        Lds=np.ones(x.shape,  dtype=np.double)
        for k in range(0,N):
            if ((k != j) and (k != i)): Lds=Lds*(x-X[k]);
        DL=DL+Lds
    
    Ldv=np.ones(x.shape,  dtype=np.double);   # DL israiskos vardiklis 
    for k in range(0,N):
         if (k != j): Ldv=Ldv*(X[j]-X[k]) 

    DL=DL/Ldv
    return DL

#-------- Lagranzo daugianaris -------------
def Lagrange(X,j,x):
    N=size(X);
    L=np.ones(x.shape,  dtype=np.double);
    for k in range(0,N) :
        if (j != k):  L=L*((x - X[k]) / (X[j] - X[k]))
    return L 
#----------------------------------

#
#*******************  Programa ************************************ 
#

T=ScrollTextBox(100,20) # sukurti teksto isvedimo langa
T.insert(tk.END,"Ermito interpoliavimas\n");
x,f,df=sympy.symbols(('x','f','df'))
f=1/(1+5*x**2)
df=f.diff(x)
xrange=np.array([-np.pi,np.pi])

nP=5      # interpoliavimo mazgu skaicius
SpausdintiMatrica(nP,T," interpoliavimo mazgu skaicius nP=")
nnn=1600  # vaizdavimo tasku skaicius
SpausdintiMatrica(nnn,T," vaizdavimo tasku skaicius  nnn=")

X=np.linspace(xrange[0],xrange[1],nP)
xxx=np.linspace(xrange[0],xrange[1],nnn)

SpausdintiMatrica(X,T,"X=")
F=np.zeros(X.shape,dtype=np.double)
DF=np.zeros(X.shape,dtype=np.double)
ffg=np.zeros(xxx.shape,dtype=np.double)
for i in range (0,nP): F[i]=f.subs(x,X[i]);DF[i]=df.subs(x,X[i])
for i in range (0,nnn): ffg[i]=f.subs(x,xxx[i]);
SpausdintiMatrica(F,T,"F=")

fff=np.zeros(xxx.shape,dtype=np.double)

# braizomos priklausomybes X(t),  Y(t):
fig1=plt.figure(1,figsize=plt.figaspect(2.5));plt.title("Ermito interpoliavimas")
ax1=fig1.add_subplot(3, 1, 1);ax1.set_xlabel('x',fontsize=14);ax1.set_ylabel('Ui',fontsize=14);ax1.grid(color='k', linestyle='-', linewidth=0.5)
ax2=fig1.add_subplot(3, 1, 2);ax2.set_xlabel('x',fontsize=14);ax2.set_ylabel('Vi',fontsize=14);ax2.grid(color='k', linestyle='-', linewidth=0.5)
ax3=fig1.add_subplot(3, 1, 3);ax3.set_xlabel('x',fontsize=14);ax3.set_ylabel('F',fontsize=14);ax3.grid(color='k', linestyle='-', linewidth=0.5)


for j in range(0,nP):
    uuu,vvv=Hermite(X,j,xxx)
    rgb=jetColormapValueToRGB(double(j/(nP-1)))
    Ut,=ax1.plot(xxx,uuu,color=rgb,label=('U'+"%2d"%(j+1)))
    Vt,=ax2.plot(xxx,vvv,color=rgb,label=('V'+"%2d"%(j+1)))
    plt.draw(); plt.pause(1);
    fff=fff+F[j]*uuu+DF[j]*vvv

ax1.legend(fontsize=14);ax2.legend(fontsize=14);
ax1.plot(X,np.zeros(X.shape),'ko'); ax2.plot(X,np.zeros(X.shape),'ko');plt.draw()

FtP,=ax3.plot(X,F,'bo') 
Ft=ax3.plot(xxx,fff,'b',label='interpoliuota funkcija');
Fgt=ax3.plot(xxx,ffg,'r',linestyle='--',label='duotoji funkcija');
ax3.legend(fontsize=14);plt.draw()
plt.show();

