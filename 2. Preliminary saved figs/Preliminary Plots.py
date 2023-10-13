# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 01:00:31 2021

@author: imman
"""

#%%

import matplotlib.pyplot as plt
import scipy as sp2
import numpy as sp
import scipy.integrate as spi
import scipy.special as sps
import random
from scipy.stats import invgauss
import numpy as np

#%%
data=sp.loadtxt("Experiment_Data copy.csv",  skiprows=7,delimiter=',', unpack=True)

#%%

#amplitude
A=data[0]

#no source callibration
cno=data[1]

#Na-22 callibration
cNa=data[2]

#Mn-54 callibration
cMn=data[3]

#Cs-137 callibration
cCs=data[4]

#Am callibration
cAm=data[5]

#Cs-137 scattering 1
sCs1=data[6]

#Cs-137 background 1
bCs1=data[7]

#Cs-137 scattering 2
sCs2=data[8]

#Cs-137 background 2
bCs2=data[9]

#Cs-137 scattering 3
sCs3=data[10]

#Cs-137 background 3
bCs3=data[11]

#%%

from scipy.optimize import curve_fit
import scipy.optimize as op
import scipy.integrate as spi
plt.style.use('seaborn-deep')
#['seaborn-bright', 'seaborn-deep', 'seaborn-paper', 'seaborn-poster',
#'Solarize_Light2', 'bmh', 'seaborn-white', 'seaborn-darkgrid', 'seaborn-dark', 'classic', 
#'dark_background', 'fivethirtyeight', 'fast', 'seaborn-dark-palette', 'seaborn-notebook', 
#'seaborn-pastel', 'seaborn-ticks', 'seaborn-talk', 'seaborn-colorblind', 'seaborn-muted', 
#'seaborn', '_classic_test', 'ggplot', 'seaborn-whitegrid', 'grayscale']
params = {   
    'axes.labelsize': 20,
   'font.size': 20,
   'legend.fontsize': 20,
   'xtick.labelsize': 20,
   'ytick.labelsize': 20,
   'figure.figsize': [10,8]
   } 
plt.rcParams.update(params)



#%% define the functions for fitting the peaks
#define voigt profile
def voi(x,s,g,a,u,m,b):
    z=((x-u)+1j*g)/(s*np.sqrt(2))
    w=np.exp(-z**2)*sps.erfc(-1j*z)             # -?erfc = 1 - erf = complementary error function
    return a*np.real(w)/(s*np.sqrt(2*np.pi))+m+b*x
x=np.linspace(-10,10,1000)
#plt.plot(x,voi(x,1,1,1,1,1,0.001))   # ?How did you get voigt width params?


#s,g width, s normal distribution width, g width of the Lorentzian
#a height multiply
#u peak frequency
#m height shifting
#b gradient

#define normal distribution
def nor(x,u,s,a,m,b):
    return a*np.exp(-(x-u)**2/s**2)+m+b*x #normal distribution plus a linear function

def nnor(x,u,s,a,m):
    return a*np.exp(-(x-u)**2/s**2)+m




#%%
    
# fit for Na peak
cNa0=cNa

'''
plt.plot(A,cNa0,'o')
plt.errorbar(A,cNa0,fmt='o',yerr=(sp.array(cNa0)**0.5),capsize=3)
plt.legend(['cNa0 freq Amplitude plot, x untrimmed'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')    
plt.show()
'''


Ap=[]
cNap=[]
for i in range(0,len(A)):
    if 120<A[i]<180:
        Ap.append(A[i])
        cNap.append(cNa0[i])
        
   
        
fitcNaV,covcNa = curve_fit(voi, Ap, cNap,[5,5,800,150,25,-0.5],sigma=1/(sp.array(cNap)**0.5))
print(fitcNaV)

for i in range(0,6):
    element = 'Na'
    covval = covcNa[i,i]
    cov = str(covval)
    err = str(covval**0.5)
    if i == 0:
        print('cov ' + element + ' Gauss-width is ' + cov)
        print('err ' + element + ' Gauss-width is ' + err)
    elif i == 1:
        print('cov ' + element + ' Lorentz-width is ' + cov)
        print('err ' + element + ' Lorentz-width is ' + err)
    elif i == 2:
        print('cov ' + element + ' Peak-Height is ' + cov)
        print('err ' + element + ' Peak-Height is ' + err)
    elif i == 3:
        print('cov ' + element + ' Peak-freq or Gauss-mean is ' + cov)
        print('err ' + element + ' Peak-freq or Gauss-mean is ' + err)
    elif i == 4:
        print('cov ' + element + ' Height-Shifting is ' + cov)
        print('err ' + element + ' Height-Shifting is ' + err)
    elif i == 5:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    elif i == 6:
        print('cov ' + element + ' Unknown Param is ' + cov)
        print('err ' + element + ' Unknown Param ' + err) #covariance, error is estimated by standard deviation
    

Ax=sp.linspace(0,500,1000)

#plt.plot(Ax,voi(Ax,*fitcNav))

fitcNa,covcNa = curve_fit(nor, Ap, cNap,[150,5,500,0,0])
print(fitcNa)

for i in range(0,5):
    element = 'Na'
    covval = covcNa[i,i]
    cov = str(covval)
    err = str(covval**0.5)
    if i == 0:
        print('cov ' + element + ' Gauss-mean is ' + cov)
        print('err ' + element + ' Gauss-mean is ' + err)
    elif i == 1:
        print('cov ' + element + ' Sigma :- Gauss-width is ' + cov)
        print('err ' + element + ' Sigma :- Gauss-width is ' + err)
    elif i == 2:
        print('cov ' + element + ' Peak-Height is ' + cov)
        print('err ' + element + ' Peak-Height is ' + err)
    elif i == 3:
        print('cov ' + element + ' Gauss shift Constant is ' + cov)
        print('err ' + element + ' Gauss shift Constant is ' + err)
    elif i == 4:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)#covariance, error is estimated by standard deviation
    
Ax=sp.linspace(0,500,1000)



name = 'Callibration Na0 Amplitude v freq plot'
plt.plot(Ap,cNap,'o') 
plt.plot(Ax,voi(Ax,*fitcNaV))
plt.plot(Ax,nor(Ax,*fitcNa))
plt.errorbar(Ap,cNap,fmt='o',yerr=(sp.array(cNap)**0.5),capsize=3)
plt.legend(['Callibration data Na-22 Amplitude v freq plot','Na0 Voigt fit plot'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.title('Na peak at 149.87±0.13')
#plt.savefig('2 '+name)
plt.show()


#plt.grid()
#plt.xlabel('amplitude',size=8)
#plt.ylabel('frequency',size=8)
#plt.title('Na peak',size=10)
#plt.show()
# peak position 149.87±0.13

#%%

# fit for Mn peak
cMn0=cMn-cno
#plt.errorbar(A,cMn0,fmt='o',yerr=np.abs((cMn0))**0.5,capsize=3)

Ap=[]
cMnp=[]
for i in range(0,len(A)):
    if 185<A[i]<300:
        Ap.append(A[i])
        cMnp.append(cMn0[i])
        
#Check for errorbars below


'''
errcC=[]
for i in range(len(cMnp)):
    if cMnp[i]>=20:
        errcC.append(cMnp[i]**0.5)
    else:
        errcC.append(0)
print(errcC)
'''
        
        
fitcMnV,covcMn = curve_fit(voi, Ap, cMnp,[8.49590479e+00,-2.36642502e+00,4.47826625e+02,2.38801836e+02,1.32477875e+01,-4.66356655e-02],sigma=1/(np.abs(sp.array(cMnp))**0.5))
print(fitcMnV)

for i in range(0,6):
    element = 'Mn'
    covval = covcMn[i,i]
    cov = str(covval)
    err = str(covval**0.5)
    if i == 0:
        print('cov ' + element + ' Gauss-width is ' + cov)
        print('err ' + element + ' Gauss-width is ' + err)
    elif i == 1:
        print('cov ' + element + ' Lorentz-width is ' + cov)
        print('err ' + element + ' Lorentz-width is ' + err)
    elif i == 2:
        print('cov ' + element + ' Peak-Height is ' + cov)
        print('err ' + element + ' Peak-Height is ' + err)
    elif i == 3:
        print('cov ' + element + ' Peak-freq or Gauss-mean is ' + cov)
        print('err ' + element + ' Peak-freq or Gauss-mean is ' + err)
    elif i == 4:
        print('cov ' + element + ' Height-Shifting is ' + cov)
        print('err ' + element + ' Height-Shifting is ' + err)
    elif i == 5:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    elif i == 6:
        print('cov ' + element + ' Unknown Param is ' + cov)
        print('err ' + element + ' Unknown Param ' + err)

    

Ax=sp.linspace(0,500,1000)



fitcMn,covcMn = curve_fit(nor, Ap, cMnp,[240,5,50,0,0],sigma=1/(np.abs(sp.array(cMnp))**0.5))
print(fitcMn)

for i in range(0,5):
    element = 'Mn'
    covval = covcMn[i,i]
    cov = str(covval)
    err = str(covval**0.5)
    if i == 0:
        print('cov ' + element + ' Gauss-mean is ' + cov)
        print('err ' + element + ' Gauss-mean is ' + err)
    elif i == 1:
        print('cov ' + element + ' Sigma :- Gauss-width is ' + cov)
        print('err ' + element + ' Sigma :- Gauss-width is ' + err)
    elif i == 2:
        print('cov ' + element + ' Peak-Height is ' + cov)
        print('err ' + element + ' Peak-Height is ' + err)
    elif i == 3:
        print('cov ' + element + ' Gauss shift Constant is ' + cov)
        print('err ' + element + ' Gauss shift Constant is ' + err)
    elif i == 4:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    
Ax=sp.linspace(0,500,1000)


# peak for Mn at 238.80±0.28


name = 'Callibration Mn0 Amplitude v freq plot'  
plt.errorbar(Ap,cMnp,fmt='o',yerr=np.abs(sp.array(cMnp))**0.5,capsize=3) 
plt.plot(Ap,cMnp,'o')
plt.plot(Ax,voi(Ax,*fitcMnV))
plt.plot(Ax,nor(Ax,*fitcMn))
plt.legend(['Callibration Mn0 Amplitude v freq plot','Mn0 Voigt fit Amplitude vs. freq(A) plot', 'Mn0 Normal-fit Amplitude vs. freq(A) plot'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.title('Mn peak at 238.80±0.28',size=10)
#plt.savefig('3 '+name)
plt.grid()
plt.show()



#%%
# fit for Cs peak
cCs0=cCs-cno
#plt.errorbar(A,cCs0,fmt='o',yerr=cCs0**0.5,capsize=3)

Ap=[]
cCsp=[]
for i in range(0,len(A)):
    if 150<A[i]<220:
        Ap.append(A[i])
        cCsp.append(cCs0[i])

#Check for errorbars below


errcC=[]
for i in range(len(cCsp)):
    if cCsp[i]>=1000:
        errcC.append(cCsp[i]**0.5)
    else:
        errcC.append(0)
print(errcC)
        
        
        
fitcCsV,covcCs = curve_fit(voi, Ap, cCsp,[ 5.70486142e+00 , 2.24738185e+00 , 2.53071838e+04,  1.92045165e+02,
  5.26171410e+02, -2.64437092e+00],sigma=1/((sp.array(errcC)**0.5+1)))
print(fitcCsV)


for i in range(0,6):
    element = 'Cs'
    covval = covcCs[i,i]
    cov = str(covval)
    err = str(covval**0.5)
    if i == 0:
        print('cov ' + element + ' Gauss-width is ' + cov)
        print('err ' + element + ' Gauss-width is ' + err)
    elif i == 1:
        print('cov ' + element + ' Lorentz-width is ' + cov)
        print('err ' + element + ' Lorentz-width is ' + err)
    elif i == 2:
        print('cov ' + element + ' Peak-Height is ' + cov)
        print('err ' + element + ' Peak-Height is ' + err)
    elif i == 3:
        print('cov ' + element + ' Peak-freq or Gauss-mean is ' + cov)
        print('err ' + element + ' Peak-freq or Gauss-mean is ' + err)
    elif i == 4:
        print('cov ' + element + ' Height-Shifting is ' + cov)
        print('err ' + element + ' Height-Shifting is ' + err)
    elif i == 5:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    elif i == 6:
        print('cov ' + element + ' Unknown Param is ' + cov)
        print('err ' + element + ' Unknown Param ' + err)
    
Ax=sp.linspace(0,500,1000)


fitcCs,covcCs = curve_fit(nor, Ap, cCsp,[192,5,1300,0,0],sigma=1/((sp.array(errcC)**0.5+1)))
print(fitcCs)

for i in range(0,5):
    element = 'Cs'
    covval = covcCs[i,i]
    cov = str(covval)
    err = str(covval**0.5)
    if i == 0:
        print('cov ' + element + ' Gauss-mean is ' + cov)
        print('err ' + element + ' Gauss-mean is ' + err)
    elif i == 1:
        print('cov ' + element + ' Sigma :- Gauss-width is ' + cov)
        print('err ' + element + ' Sigma :- Gauss-width is ' + err)
    elif i == 2:
        print('cov ' + element + ' Peak-Height is ' + cov)
        print('err ' + element + ' Peak-Height is ' + err)
    elif i == 3:
        print('cov ' + element + ' Gauss shift Constant is ' + cov)
        print('err ' + element + ' Gauss shift Constant is ' + err)
    elif i == 4:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    
Ax=sp.linspace(0,500,1000)


# peak for Cs at 192.00±0.05



name = 'callibration Cs0 Amplitude v freq(A) plot'
plt.errorbar(Ap,cCsp,fmt='o',yerr=sp.array(cCsp)**0.5,capsize=3)
plt.plot(Ap,cCsp,'o')
plt.plot(Ax,voi(Ax,*fitcCsV))
plt.plot(Ax,nor(Ax,*fitcCs))
plt.legend(['callibration Cs0 Amplitude v freq(A) plot', 'Cs0 Voigt fit', 'Cs0 Normal-fit'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.title('Cs peak at 192.00±0.05',size=10)
#plt.savefig('4 '+name)
plt.show()

#%%

# fit for Am peak
cAm0=cAm-cno
#plt.errorbar(A,cAm0,fmt='o',yerr=cAm0**0.5,capsize=3)

Ap=[]
cAmp=[]
for i in range(0,len(A)):
    if 13<A[i]<50:
        Ap.append(A[i])
        cAmp.append(cAm0[i])
        
#plot errorbars later 

errcC=[]
for i in range(len(cAmp)):
    if cAmp[i]>=1000:
        errcC.append(cAmp[i]**0.5)
    else:
        errcC.append(0)
print(errcC)
        
        
fitcAmV,covcAm = curve_fit(voi, Ap, cAmp,[ 1.50901365e+00,  8.91068122e-02,  3.51782797e+04,  2.17232445e+01,
  2.45872987e+01, -6.67251889e-01],sigma=1/((sp.array(errcC)**0.5+1)))
print(fitcAmV)


for i in range(0,6):
    element = 'Am'
    covval = covcAm[i,i]
    cov = str(covval)
    err = str(covval**0.5)
    if i == 0:
        print('cov ' + element + ' Gauss-width is ' + cov)
        print('err ' + element + ' Gauss-width is ' + err)
    elif i == 1:
        print('cov ' + element + ' Lorentz-width is ' + cov)
        print('err ' + element + ' Lorentz-width is ' + err)
    elif i == 2:
        print('cov ' + element + ' Peak-Height is ' + cov)
        print('err ' + element + ' Peak-Height is ' + err)
    elif i == 3:
        print('cov ' + element + ' Peak-freq or Gauss-mean is ' + cov)
        print('err ' + element + ' Peak-freq or Gauss-mean is ' + err)
    elif i == 4:
        print('cov ' + element + ' Height-Shifting is ' + cov)
        print('err ' + element + ' Height-Shifting is ' + err)
    elif i == 5:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    elif i == 6:
        print('cov ' + element + ' Unknown Param is ' + cov)
        print('err ' + element + ' Unknown Param ' + err)
    
Ax=sp.linspace(0,50,1000)

fitcAm,covcAm = curve_fit(nor, Ap, cAmp,[22,5,9000,0,0],sigma=1/((sp.array(errcC)**0.5+1)))
print(fitcAm)

for i in range(0,5):
    element = 'Am'
    covval = covcAm[i,i]
    cov = str(covval)
    err = str(covval**0.5)
    if i == 0:
        print('cov ' + element + ' Gauss-mean is ' + cov)
        print('err ' + element + ' Gauss-mean is ' + err)
    elif i == 1:
        print('cov ' + element + ' Sigma :- Gauss-width is ' + cov)
        print('err ' + element + ' Sigma :- Gauss-width is ' + err)
    elif i == 2:
        print('cov ' + element + ' Peak-Height is ' + cov)
        print('err ' + element + ' Peak-Height is ' + err)
    elif i == 3:
        print('cov ' + element + ' Gauss shift Constant is ' + cov)
        print('err ' + element + ' Gauss shift Constant is ' + err)
    elif i == 4:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    
Ax=sp.linspace(0,50,1000)



name = 'Americium-241 Amplitude v freq. plot'
plt.plot(Ap,cAmp,'bo',ms=13)   
plt.errorbar(Ap,cAmp,fmt='bo',yerr=sp.array(cAmp)**0.5,capsize=13)
plt.plot(Ax, voi(Ax,*fitcAmV), 'm-', ms=13)
plt.plot(Ax,nor(Ax,*fitcAm), 'g-', ms=13)
plt.grid()
plt.legend(['Callibration Data','Am-241 Voigt fit plot', 'Am-241 Gauss-fit'])
plt.xlabel('Amplitude', size = 21)
plt.ylabel('Frequency', size = 21)
plt.xticks(size=18)
plt.yticks(size=18)
plt.xlim(12,40)
plt.title('Am-241 peak at 21.720±0.005',size=21)
#plt.savefig('Americium')
plt.show()


#%%
#fitting for the conversion

#amplitude scale:
Amplitude=np.array([149.94,238.8,192.00,21.72])  #Array of Peak Amplitude Freq.
erramplitude=np.array([0.11,0.28,0.05,0.01])

#the emission peaks
energy=np.array([511,821,662,59.6])


#linear fit for converting amplitude to energy
fit = np.polyfit(Amplitude,energy,1,w=1/erramplitude)
pc=np.poly1d(fit)
T=np.linspace(0,250,1000)

name =  'Peak Amplitude Frequencies Vs Known Peak Energies of Elements'
plt.errorbar(Amplitude, energy, xerr=erramplitude,fmt='o')
plt.plot(T,pc(T))
plt.legend(['Peak Amplitude Frequencies Vs Known Peak Energies of Elements','Best-fit Amplitude vs Energy Correlation'])
plt.xlabel('Amplitude')
plt.ylabel('Energy')
plt.title('Energy vs Amplitude (Conversion & Proportionality)')
plt.savefig('6 '+name)
plt.grid()
plt.show()

print('Amplitude Energy Gradient is ', fit)

expenergy=pc(Amplitude)


#average perventage error
delta=expenergy-energy
print('Delta Expected Energy - Measured is ', delta)

a=0
for i in range(len(delta)):
    a+=(delta[i]/energy[i])**2
    
#print(a)
toerror=a**0.5
print('total percentage error from fitting the points',toerror)





