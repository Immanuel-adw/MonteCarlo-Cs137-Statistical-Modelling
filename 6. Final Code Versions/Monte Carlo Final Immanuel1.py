#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 11:50:33 2021

@author: Ginny
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
    'axes.labelsize': 10,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'figure.figsize': [7.9,8]
   } 
plt.rcParams.update(params)

#%%
'''
#overall plot for peaks


name = 'Background Radiation'
plt.plot(A,cno,'o')
plt.legend(['Background Radiation'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('A) '+name)
plt.show()


name = 'Na radiation'
plt.plot(A,cNa,'o')
plt.legend(['Na radiation - Amplitude v FreqCount'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('B) '+name)
plt.show()

name = 'Mn radiation'
plt.plot(A,cMn,'o')
plt.legend(['Mn'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('C) '+name)
plt.show()

name = 'Cs radiation'
plt.plot(A,cCs,'o')
plt.legend(['Cs radiation - Amplitude v FreqCount'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('D) '+name)
plt.show()

namee = 'Am radiation'
plt.plot(A,cAm,'o')
plt.legend(['Am radiation - Amplitude v FreqCount'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('E) '+name)
plt.show()

#plt.legend(['no','Na','Mn','Cs','Am'])
'''

#%% define the functions for fitting the peaks
#define voigt profile
def voi(x,s,g,a,u,m,b):
    z=((x-u)+1j*g)/(s*np.sqrt(2))
    w=np.exp(-z**2)*sps.erfc(-1j*z)             # -?erfc = 1 - erf = complementary error function
    return a*np.real(w)/(s*np.sqrt(2*np.pi))+m+b*x
x=np.linspace(-10,10,1000)
plt.plot(x,voi(x,1,1,1,1,1,0.001))   # ?How did you get voigt width params?


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
    
#if the peak cannot be fitted with voigt profile.

#%%
    
# fit for Na peak
cNa0=cNa-cno
plt.plot(A,cNa0,'o')
plt.errorbar(A,cNa0,fmt='o',yerr=(sp.array(cNa0)**0.5),capsize=3)
plt.legend(['cNa0 freq Amplitude plot, x untrimmed'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')    
plt.show()



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
plt.legend(['Callibration Na0 Amplitude v freq plot','Na0 Voigt fit plot'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.title('Na peak at 149.87±0.13')
plt.savefig('2 '+name)
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
plt.savefig('3 '+name)
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
plt.savefig('4 '+name)
plt.show()

#plt.grid()


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
    
Ax=sp.linspace(0,500,1000)

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
    
Ax=sp.linspace(0,250,1000)



name = 'Callibration Am Amplitude v freq(A) plot' 
plt.errorbar(Ap,cAmp,fmt='o',yerr=sp.array(cAmp)**0.5,capsize=3)
plt.plot(Ap,cAmp,'bo')  
plt.plot(Ax, voi(Ax,*fitcAmV))
plt.plot(Ax,nor(Ax,*fitcAm))
plt.grid()
plt.legend(['Callibration Am0 Amplitude v freq(A) plot','Am0 Voigt fit plot', 'Am0 Gauss-fit'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.title('Am peak at 21.72±0.01',size=10)
plt.grid()
plt.savefig('Americium')
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



#%%>>>>>>>energy and peak width at 0 degrees

# peak width is proportional to the square root of energy

E=sp.array([149.81,238.80,192,21.72])
errE=sp.array([0.13,0.28,0.05,0.01])
W=sp.array([7.93,10.49,9.58,2.2]) 
errW=sp.array([0.17,0.43,0.18,0.01]) # How was this calculated ?

errEsqrt=0.5*E**(-0.5)*errE



fit = sp.polyfit(W,E**0.5,1,w=1/errW)
p=sp.poly1d(fit)
T=sp.linspace(2,12,100)


name = 'Peak-Width vs Root-Energy(amplitude)  at Zero degrees'
plt.plot(T,p(T))
plt.errorbar(W,E**0.5,fmt='o',capsize=3,yerr=errEsqrt,xerr=errW)
plt.xlabel('Width W',size=10)
plt.ylabel('square root of energy $E^{0.5}$',size=10)
plt.title('Peak-Width vs $E^{1/2}$')
plt.savefig(name)
plt.grid()
plt.show()


#%% >>>>>fitting for the position of back scattering peaks

# fit for Cs peak
cCs0=cCs-cno
#plt.errorbar(A,cCs0,fmt='o',yerr=cCs0**0.5,capsize=3)

Ap=[]
cCsp=[]
for i in range(0,len(A)):
    if 50<A[i]<100:
        Ap.append(A[i])
        cCsp.append(cCs0[i])
        
#Check for Errorbar below       

fitcCs,covcCs = curve_fit(nor, Ap, cCsp,[55,5,300,0,0])
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

plt.plot(Ax,nor(Ax,*fitcCs))

# peak for Cs at 192.00±0.05

name = 'Cs 0 degrees fitting for back-scattering Peaks within the Detector'
plt.errorbar(Ap,cCsp,fmt='o',yerr=sp.array(cCsp)**0.5,capsize=3)
plt.plot(Ax,nor(Ax,*fitcCs))
plt.xlabel('amplitude',size=8)
plt.ylabel('frequency',size=8)
plt.title('Cs - Detector back-scattering Peaks',size=10)
plt.grid()
plt.show()

#peak at 192, main peak 0 degree
#peak at 122.89, back-scattering give max energy electron, which de-excite
#peak at 59.46 back-scattered photon
#peak at 27.01 characteristic x-ray of iodine in the scintillation material
#peak at 11.47 characteristic x-ray of lead on the equipment

#back scattering is more dominant than other angles---Klein-Nishina scattering cross-section



#%%#energy peak-width relationship for Cs-Cesium at 0,20,30 & 45 deg
E=sp.array([192,59.46,27.01,11.47])
errE=sp.array([0.05,0.64,0.17,0.10])
W=sp.array([9.45,6.66,3.49,2.14])
errW=sp.array([0.09,1.06,0.29,0.16])

errEsqrt=0.5*E**(-0.5)*errE


uncert = errEsqrt * errW

fit = sp.polyfit(W,E**0.5,1,w=1/uncert)
p=sp.poly1d(fit)
T=sp.linspace(2,12,100)


name = 'Root-Energy v. Peak-Width for Cs-Cesium at 0,20,30 & 45 deg.'
plt.plot(T,p(T))
plt.errorbar(W,E**0.5,fmt='o',capsize=3,yerr=errEsqrt,xerr=errW)
plt.xlabel('Peak-Width W',size=10)
plt.ylabel('Square root of Energy $E^{0.5}$',size=10)
plt.title(name)
plt.savefig(name)
plt.grid()
plt.show()

#peak at 122.89 is excluded since it has different broadening mechanism Klein-Nishina scattering. The photons emitted by electron de-excitation can be in any direction, introducing another peak broadeing step.


#%%
# COMPTON TESTING

AA=pc(A)

cCs0=cCs-cno  #
#plt.errorbar(A,cCs0,fmt='o',yerr=cCs0**0.5,capsize=3)


Apdud = []
sCspdud = []
Ap=[]
sCsp=[]
for i in range(0,len(A)):
    if 500<AA[i]<750:
        Ap.append(AA[i])
        sCsp.append(cCs0[i])
    else:
        Apdud.append(0)
        sCspdud.append(0)

#Check for errorbars below


errcC=[]
for i in range(len(sCsp)):
        errcC.append(sCsp[i]**0.5)
    
print(errcC)
        
plt.plot(Ap,sCsp,'o')
plt.show()        
        
fitcCsV,covcCs = curve_fit(voi, Ap, sCsp,[ 82, 81,  1400,  680,
  0, 0], sigma = 1/np.array(errcC)+1)
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


fitcCs,covcCs = curve_fit(nor, Ap, sCsp,[682,80,1400,0,0],sigma=1/((sp.array(errcC)+1)))
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
    
Ax=sp.linspace(500,900,1000)


# peak for Cs at 192.00±0.05



name = 'callibration Cs0 Amplitude v freq(A) plot'
plt.errorbar(Ap,sCsp,fmt='bo',yerr=sp.array(sCsp)**0.5,capsize=3)
#plt.plot(Ap,cCsp,'o')
plt.plot(Ax,voi(Ax,*fitcCsV))
plt.plot(Ax,nor(Ax,*fitcCs))
plt.legend(['Compton Test at 0 Degrees - Cs0 Amplitude v freq(A) plot', 'Cs0 Voigt fit', 'Cs0 Normal-fit'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.title('Cs peak at 0 Degrees',size=10)
plt.savefig('4 '+name)
plt.show()

#plt.grid()




#%%
#plt.plot(A,cCs/10)


#plt.errorbar(AA,sCs1-bCs1,fmt='o',yerr=np.abs(sCs1-bCs1)**0.5,capsize=3)

A1=[]
sCsp1=[]
for i in range(0,len(A)):
    if 400<AA[i]<900:
        A1.append(AA[i])
        sCsp1.append(sCs1[i]-bCs1[i])
print(sCsp1)       
# See errorbars below

'''
#group the data
gA1=[]
gsCsp1=[]
for i in range(0,len(A1)-8):
    gA1.append(sp.mean([A1[i],A1[i+1],A1[i+2],A1[i+3],A1[i+4],A1[i+5],A1[i+6],A1[i+7]]))
    gsCsp1.append(sp.mean([sCsp1[i],sCsp1[i+1],sCsp1[i+2],sCsp1[i+3],sCsp1[i+4],sCsp1[i+5],sCsp1[i+6],sCsp1[i+7]]))
plt.errorbar(gA1,gsCsp1,fmt='o',yerr=sp.array(gsCsp1)**0.5)


    '''    

    
Ax=sp.linspace(0,1750,1000)


fitsCsV,covsCs = curve_fit(voi, A1, sCsp1,[ 4.19738061e+01, 4.14325173e+01,  3.31997291e+02,  5.94834892e+02,
  8.53485283e-01, -5.99601714e-04], sigma = 1/(np.abs(sp.array(sCsp1))**0.5 +0.4))
print(fitsCsV)

for i in range(0,6):
    element = 'Cs'
    covval = covsCs[i,i]
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
    
    
Ax=sp.linspace(0,1750,1000)

fitsCs,covsCs = curve_fit(nor, A1, sCsp1,[ 5.95159322e+02,  3.73930801e+01,  9.32464930e+00, -3.72857332e-01,
  5.19971127e-04], sigma = 1/((np.abs(sp.array(sCsp1))**0.5)+0.004))
print(fitsCs)

for i in range(0,4):
    element = 'Cs'
    covval = covsCs[i,i]
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
        print('cov ' + element + ' Peak-freq or Gauss-mean is ' + cov)
        print('err ' + element + ' Peak-freq or Gauss-mean is ' + err)
    elif i == 4:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)



name = 'Cs-Experiment peak for 20 degrees at 610_01±6_80 MeV'
name2 = 'Cs-Experiment peak for 20 degrees at 610.01±6.80 MeV'
plt.errorbar(A1,sCsp1,fmt='gx',yerr=1/np.abs(sp.array(sCsp1))**0.5,capsize=3,ms=1)
#plt.plot(A1,sCsp1,'gx', ms=1.5)
plt.plot(Ax,voi(Ax,*fitsCsV), 'm-')
plt.plot(Ax,nor(Ax,*fitsCs), 'b-')
plt.legend(['Energy vs FreqCount 20deg Cs plot', 'Voigt-fit', 'Gauss-Fit'])
plt.xlabel('Energy MeV')
plt.ylabel('Frequency (Counts)')
plt.title(name2)
plt.savefig(name)
plt.show()


#plt.plot(A,sCs2,'o')
#plt.plot(A,sCs3,'o')
#plt.legend(['20','30','45'])


print(0.49471538412305943**0.5+6.10)
print('the peak for 20 degree is at 610.01±6.80')





#%%

#plt.errorbar(AA,sCs2-bCs2,fmt='o',yerr=(np.abs(sCs2-bCs2)**0.5),capsize=3)

A2=[]
sCsp2=[]
for i in range(0,len(AA)):
    if 400<AA[i]<750:
        A2.append(AA[i])
        sCsp2.append(sCs2[i]-bCs2[i])

# See Errorbars below


'''
gA2=[]
gsCsp2=[]
for i in range(0,len(A2)-8):
    gA2.append(sp.mean([A2[i],A2[i+1],A2[i+2],A2[i+3],A2[i+4],A2[i+5],A2[i+6],A2[i+7]]))
    gsCsp2.append(sp.mean([sCsp2[i],sCsp2[i+1],sCsp2[i+2],sCsp2[i+3],sCsp2[i+4],sCsp2[i+5],sCsp2[i+6],sCsp2[i+7]]))
plt.plot(gA2,gsCsp2,'o')
'''


fitsCs2V,covsCs2 = curve_fit(voi, A2, sCsp2,[10,10,20,550,0,0], sigma = 1/(sp.sqrt(sp.array(np.abs(sCsp2)))+0.03))
print(fitsCs2V)

for i in range(0,6):
    element = 'Cs'
    covval = covsCs2[i,i]
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
    
Ax=sp.linspace(0,1750,1000)

        
fitsCs2,covsCs2 = curve_fit(nor, A2, sCsp2,[570,5,20,0,-0.05], sigma = 1/sp.sqrt(sp.array(np.abs(sCsp2))))
print(fitsCs2)

for i in range(0,4):
    element = 'Cs'
    covval = covsCs2[i,i]
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
        print('cov ' + element + ' Peak-freq or Gauss-mean is ' + cov)
        print('err ' + element + ' Peak-freq or Gauss-mean is ' + err)
    elif i == 4:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    


name = 'Cs-Experiment peak 30 degrees at 549_13±6_25 MeV' 
name2 = 'Cs-Experiment peak 30 degrees at 549.13±6.25 MeV' 
plt.errorbar(A2,sCsp2,fmt='yx',yerr=sp.sqrt(sp.array(np.abs(sCsp2))),capsize=3, ms=1.8)    
plt.plot(A2,sCsp2,'yx', ms=1)
plt.plot(Ax,voi(Ax,*fitsCs2V), 'm-')
plt.plot(Ax,nor(Ax,*fitsCs2), 'b-')
plt.legend(['Energy vs FreqCount 30deg Cs plot', 'Voigt-fit', 'Gauss-Fit'])
plt.xlabel('Energy MeV')
plt.ylabel('Frequency (Counts)')
plt.title(name2)
plt.savefig(name)
plt.show()
#plt.grid()


'''
plt.plot(sp.linspace(0,995,200),ampd)
plt.plot(gA2,gsCsp2,'o')

'''
        
print(0.5842594969439249**0.5+5.49)
print('the peak for 30 degree is at 549.13±6.25')


#%%
#plt.errorbar(AA,sCs3-bCs3,fmt='o',yerr=np.abs(sCs3-bCs3)**0.5,capsize=3)

A3=[]
sCsp3=[]
for i in range(0,len(AA)):
    if 400<AA[i]<750:
        A3.append(AA[i])
        sCsp3.append(sCs3[i]-bCs3[i])
        
# See Errorbars later 

        
gA3=[]
gsCsp3=[]
for i in range(0,len(A3)-8):
    gA3.append(sp.mean([A3[i],A3[i+1],A3[i+2],A3[i+3],A3[i+4],A3[i+5],A3[i+6],A3[i+7]]))
    gsCsp3.append(sp.mean([sCsp3[i],sCsp3[i+1],sCsp3[i+2],sCsp3[i+3],sCsp3[i+4],sCsp3[i+5],sCsp3[i+6],sCsp3[i+7]]))
#plt.plot(gA3,gsCsp3,'o')
        


fitsCs3V,covsCs3 = curve_fit(voi, A3, sCsp3,[50,50,30,470,0,0], sigma=1/(np.abs(sCsp3)**0.5))
print(fitsCs2)

for i in range(0,6):
    element = 'Cs'
    covval = covsCs3[i,i]
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
     
        

fitsCs3,covsCs3 = curve_fit(nor, A3, sCsp3,[500,3,30,0,-0.05], sigma=1/((np.abs(sCsp3)**0.5)))
print(fitsCs3)

for i in range(0,4):
    element = 'Cs'
    covval = covsCs3[i,i]
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
        print('cov ' + element + ' Peak-freq or Gauss-mean is ' + cov)
        print('err ' + element + ' Peak-freq or Gauss-mean is ' + err)
    elif i == 4:
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
 
    

Ax=sp.linspace(400,800,1000)
        
    
name = 'Cs-Experiment peak 45 degrees at 470_21±5_41 MeV'
name2 = 'Cs-Experiment peak 45 degrees at 470.21±5.41 MeV'
plt.errorbar(A3,sCsp3,yerr=np.abs(sCsp3)**0.5, fmt = 'gx',capsize=3, markersize = 1)    
plt.plot(A3,sCsp3,'gx', ms=2)  
#plt.errorbar(gA3,gsCsp3,yerr=np.abs(gsCsp3)**0.5,fmt = '', capsize=3)    
#plt.plot(gA3,gsCsp3,'bx', markersize = 2)
plt.plot(Ax,voi(Ax,*fitsCs3V), 'm-')
plt.plot(Ax,nor(Ax,*fitsCs3), 'b-')

plt.legend(['Energy vs FreqCount 45deg Cs plot', 'Voigt-fit', 'Gauss-Fit'])

#plt.legend(['E vs Freq Smoothed - bins-of-8', 'Voigt-fit', 'Gauss-Fit'])

plt.xlabel('Energy MeV')
plt.ylabel('Frequency (Counts)')
plt.title(name2)
plt.savefig(name)
plt.show()

    
print(0.5035778448149948**0.5+4.70)
print('the peak for 45 degree is at 470.21±5.41')



#%%
'''
Eth=np.array([662,609.41,549.43,470.41])*1.6e-16
errEth=np.array([0.1,2.09,3.18,1.98])*1.6e-16
ang=np.array([0,20,30,45])
'''


Eth=np.array([662,598.55,549.86,475.25])*1.6e-16
errEth=np.array([0.15,3.54, 8.73,4.12])*1.6e-16
ang=np.array([0,20,30,45])


fit = np.polyfit(np.cos(ang*np.pi/180),1/Eth,1,w=1/(errEth/Eth**2))
#fit, cov = np.polyfit(np.cos(ang*np.pi/180),1/Eth,1,w=1/(errEth))
p=np.poly1d(fit)
#[y, delta] = np.polyval(p,np.cos(ang*np.pi/180))
T=np.linspace(0.67,1.05,1000)
print('fit 1/E vs Costheta is ', fit)
#print('delta S is ', S)
#print('cov 1/E vs Costheta is ', cov1)

name = 'Best-fit line for Inverse of Peak-Energy vs Cos-theta - Cs plot'
plt.errorbar(np.cos(ang*np.pi/180),1/Eth,yerr=errEth/Eth**2,fmt='o')
plt.plot(T,p(T))
plt.legend(['Cos-theta vs 1/Peak-Energy Cs plot', 'Best-fit line for Cos-theta vs 1/Peak-Energy Cs plot'])
plt.grid()
plt.xlabel('Cos-theta')
plt.ylabel('1/Peak-Energy')
plt.savefig('Energy Amplitude Conversion')
plt.show()

#plt.plot(np.cos(ang*np.pi/180),p(np.cos(ang*np.pi/180)))

dif=(1/Eth)-p(np.cos(ang*np.pi/180))
#print(dif)

var=0
for i in range(len(dif)):
    var+=(dif[i]/p(np.cos(ang*np.pi/180))[i])**2
    
print('variance',var*fit,'sd',(var*abs(fit))**0.5)

#%%

#Monte Carlo Simulation 0

#source 
E0=662 # keV#input photon energy

#target
erm=511 #keV #electron rest mass
sa=0 #degree #scattering angle

#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

#ideal signal

def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))

Ei=E(E0) #ideal energy
Ai=Ei*gain # ideal analog signal
print('ideal energy',Ei)


#first set


n=9000 #the number of points

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))
    

# noise
s=res*Ai/2
print('s',s)
noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
#noise=invgauss.ppf(sp.array(randseed),s,loc=0,scale=1)


# outcome bin
Ob=(Ai+noise)/maxsig*2**(bitdepth)
Obf=sp.floor(Ob)
#print(Obf)


#outcome set

#count the amplitude distribution
num=105 #number of intervals
ampd=sp.array([0]*num)

for i in range(0,n):
    for j in range(0,num):
        if 5*j<=Obf[i]<5*j+5:
            ampd[j]+=1
        
        
#print(ampd)
#print( sp.linspace(0,1000,201))   

plt.errorbar(sp.linspace(0,520,105)/0.512,ampd,fmt='o',yerr=sp.sqrt(ampd),capsize=3)
plt.errorbar(A*3.53-17,cCs-cno,fmt='o',yerr=sp.sqrt(np.abs(cCs-cno)))
#for a<Obf<=b:


#fitting the curve
fitcCs,covcCs = curve_fit(nnor, (A*3.53-17), (cCs-cno),[600,5,1000,100])
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,1750,1000)
plt.plot(Ax,nnor(Ax,*fitcCs))


fitcCs1,covcCs1 = curve_fit(nnor, sp.linspace(0,520,105)/0.512,ampd,[660,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)

for i in range(0,4):
    print(covcCs1[i,i])
    
plt.plot(Ax,nnor(Ax,*fitcCs1))

plt.title('Cs Monte-Carlo Frequency-energy graph at 0 degree')
plt.xlabel('Energy, keV')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('ZERO deg')
plt.show()

#%%

#Monte Carlo Simulation 20

#source 
E0=662 # keV#input photon energy

#target
erm=511 #keV #electron rest mass
sa=20 #degree #scattering angle

#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

#ideal signal

def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))

Ei=E(E0) #ideal energy
Ai=Ei*gain # ideal analog signal
print('ideal energy',Ei)


#first set

n=45 #the number of points

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))
    
# noise
s=res*Ai/2
noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
#noise=invgauss.ppf(sp.array(randseed),s,loc=0,scale=1)

# outcome bin
Ob=(Ai+noise)/maxsig*2**(bitdepth)
Obf=sp.floor(Ob)
#print(Obf)

#outcome set

#count the amplitude distribution
num=105 #number of intervals
ampd=sp.array([0]*num)

for i in range(0,n):
    for j in range(0,num):
        if 5*j<=Obf[i]<5*j+5:
            ampd[j]+=1
              
#print(ampd)
#print( sp.linspace(0,1000,201))   
#for a<Obf<=b:


#Fitting the Normal Curve for the Experiment
fitcCs,covcCs = curve_fit(nnor, sp.array(A1),sCsp1,[615,10,50,10])
print('experiment, U,S,A,M',fitcCs)
# A*3.53 - 17.15 is amplitude to energy conversion, A is the background radiation
for i in range(0,4):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,1750,1000)

#Fitting the Voigt Curve for the Experiment
xx1=sp.linspace(0,520,105)/0.512   
fitcCsV,covcCsV = curve_fit(voi, sp.array(A1),sCsp1,[10,10,30,600,0,0])
print('experiment, U,S,A,M',fitcCs)
for i in range(0,4):
    print(covcCs[i,i])  
Ax=sp.linspace(250,780,1000)
#plt.plot(Ax,voi(Ax,*fitcCsV))


#Fitting the Normal Curve for the Simulation 
xx1=sp.linspace(0,520,105)/0.512
fitcCs1,covcCs1 = curve_fit(nnor, xx1,ampd,[595,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)
for i in range(0,4):
    print(covcCs1[i,i])
   
plt.plot(sp.array(A1),sCsp1,'bo') #Experiment data 20 degrees Cs 
#plt.plot([],[],'bo') #Experiment data 20 degrees Cs
plt.plot(sp.linspace(0,520,105)/0.512,ampd,'mo') #simulation data 20 degrees Cs 
plt.errorbar(sp.array(A1),sCsp1,fmt='bo',yerr=sp.sqrt(np.abs(sCsp1)), capsize=3)
plt.errorbar(sp.linspace(0,520,105)/0.512, ampd,fmt='mo',yerr=sp.sqrt(ampd),capsize=3)
plt.legend(['Experiment-data', 'Simulation-data'])
plt.title('Cs Monte-Carlo Frequency-energy graph at 20 degree')
plt.grid()
plt.show()


plt.plot(Ax,nnor(Ax,*fitcCs), 'b-')  # Experiment-fit at 20 degrees
plt.plot(Ax,nnor(Ax,*fitcCs1), 'm-') # Simulation-fit at 20 degrees
plt.legend(['Experiment-fit', 'Simulation-fit'])
plt.title('Cs Monte-Carlo Frequency-energy graph at 20 degree')
plt.xlabel('Energy, keV')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('TWENTY deg')
plt.show()



#%%

#Monte Carlo Simulation 30

#source 
E0=662 # keV#input photon energy

#target
erm=511 #keV #electron rest mass
sa=30 #degree #scattering angle

#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

#ideal signal

def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))

Ei=E(E0) #ideal energy
Ai=Ei*gain # ideal analog signal
print('ideal energy',Ei)


#first set


n=35 #the number of points

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))
    

# noise
s=res*Ai/2
noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
#noise=invgauss.ppf(sp.array(randseed),s,loc=0,scale=1)


# outcome bin
Ob=(Ai+noise)/maxsig*2**(bitdepth)
Obf=sp.floor(Ob)
#print(Obf)


#outcome set

#count the amplitude distribution
num=105 #number of intervals
ampd=sp.array([0]*num)

for i in range(0,n):
    for j in range(0,num):
        if 5*j<=Obf[i]<5*j+5:
            ampd[j]+=1
        
        
#print(ampd)
#print( sp.linspace(0,1000,201))   
#for a<Obf<=b:


fitcCs,covcCs = curve_fit(nnor, sp.array(A2),sCsp2,[570,5,20,10]) # I added arg 2
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])


fitcCs1,covcCs1 = curve_fit(nnor, sp.linspace(0,520,105)/0.512,ampd,[595,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)

for i in range(0,4):
    print(covcCs1[i,i])

Ax=sp.linspace(0,1750,1000)


'''
plt.plot(sp.array(A2),sCsp2,'bo') #Experiment data 30 degrees Cs 
plt.plot(sp.linspace(0,520,105)/0.512,ampd,'mo') #simulation data 30 degrees Cs 
plt.plot(Ax,nnor(Ax,*fitcCs), 'b-')  # Experiment-fit at 30 degrees
plt.plot(Ax,nnor(Ax,*fitcCs1), 'm-') # Simulation-fit at 30 degrees
plt.errorbar(sp.array(A2),sCsp2,fmt='bo',yerr=sp.sqrt(np.abs(sCsp2)))
plt.errorbar(sp.linspace(0,520,105)/0.512, ampd,fmt='mo',yerr=sp.sqrt(ampd),capsize=3)
plt.legend(['Experimental-data', 'Simulation-data','Experimental-fit', 'Simulation-fit'])
plt.title('Cs Monte-Carlo Frequency-energy graph at 30 degree')
plt.grid()
plt.show()
'''


plt.plot(sp.array(A2),sCsp2,'bo') #Experiment data 30 degrees Cs 
#plt.plot([],[],'bo') #Experiment data 30 degrees Cs
plt.plot(sp.linspace(0,520,105)/0.512,ampd,'mo') #simulation data 30 degrees Cs 
plt.errorbar(sp.array(A2),sCsp2,fmt='bo',yerr=sp.sqrt(np.abs(sCsp2)), capsize=3)
plt.errorbar(sp.linspace(0,520,105)/0.512,ampd,fmt='mo',yerr=sp.sqrt(ampd),capsize=3)
plt.legend(['Experiment-data', 'Simulation-data'])
plt.title('Cs Monte-Carlo Frequency-energy graph at 30 degree')
plt.grid()
plt.show()

plt.plot(Ax,nnor(Ax,*fitcCs), 'b-')  # Experiment-fit at 30 degrees
plt.plot(Ax,nnor(Ax,*fitcCs1), 'm-') # Simulation-fit at 30 degrees
plt.legend(['Experiment-fit', 'Simulation-fit'])
plt.title('Cs Monte-Carlo Frequency-energy graph at 30 degree')
plt.xlabel('Energy, keV')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('THIRTY deg')
plt.show()


#%%

#Monte Carlo Simulation 45

#source 
E0=662 # keV#input photon energy

#target
erm=511 #keV #electron rest mass
sa=45 #degree #scattering angle

#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

#ideal signal

def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))

Ei=E(E0) #ideal energy
Ai=Ei*gain # ideal analog signal
print('ideal energy',Ei)


#first set


n=90 #the number of points

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))
    
# noise
s=res*Ai/2
noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
#noise=invgauss.ppf(sp.array(randseed),s,loc=0,scale=1)


# outcome bin
Ob=(Ai+noise)/maxsig*2**(bitdepth)
Obf=sp.floor(Ob)
#print(Obf)


#outcome set

#count the amplitude distribution
num=105 #number of intervals
ampd=sp.array([0]*num)

for i in range(0,n):
    for j in range(0,num):
        if 5*j<=Obf[i]<5*j+5:
            ampd[j]+=1
               
#print(ampd)
#print( sp.linspace(0,1000,201))   

#for a<Obf<=b:

fitcCs,covcCs = curve_fit(nnor, sp.array(A3),sCsp3,[500,25,25,0])
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])

fitcCs1,covcCs1 = curve_fit(nnor, sp.linspace(0,520,105)/0.512,ampd,[500,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)

for i in range(0,4):
    print(covcCs1[i,i])
    
Ax=sp.linspace(0,1750,1000)

plt.plot(sp.array(A3),sCsp3,'bo') #Experiment data 45 degrees Cs 
#plt.plot([],[],'bo') #Experiment data 45 degrees Cs
plt.plot(sp.linspace(0,520,105)/0.512,ampd,'mo') #simulation data 45 degrees Cs 
plt.errorbar(sp.array(A1),sCsp1,fmt='bo',yerr=sp.sqrt(np.abs(sCsp1)), capsize=3)
plt.errorbar(sp.linspace(0,520,105)/0.512, ampd,fmt='mo',yerr=sp.sqrt(ampd),capsize=3)
plt.legend(['Experiment-data', 'Simulation-data'])
plt.title('Cs Monte-Carlo Frequency-energy graph at 45 degree')
plt.grid()
plt.show()

plt.plot(Ax,nnor(Ax,*fitcCs), 'b-')  # Experiment-fit at 45 degrees
plt.plot(Ax,nnor(Ax,*fitcCs1), 'm-') # Simulation-fit at 45 degrees
plt.legend(['Experiment-fit', 'Simulation-fit'])
plt.title('Cs Monte-Carlo Frequency-energy graph at 45 degree')
plt.xlabel('Energy, keV')
plt.ylabel('Frequency')
plt.grid()
plt.savefig('FOURTY FIVE deg')
plt.show()


#%%
#probability distribution over the angle
def sec(th,Eo): # define Klein-Nishina cross section
    a=1/(1+Eo/511*(1-th))
    return 0.5*(2.8179e-15)**2*(a)**2*((a+1/a)-2*(1-th**2))/1.3350666520239514e-29 #normalised

def E(t,EE):
    return EE/(1+EE/511*(1-t))        

  
norm=spi.quad(sec,0,sp.pi,args=662)
x=sp.arange(-1,1,0.01) # the angle


print(norm)
y=sec(x,662)

plt.plot(x,y)   
plt.xlabel('Scattering angle')
plt.ylabel('Cross-section')
plt.legend(['Cos(theta)-Scattering angle v. Cross-Section'])
plt.title('Klein-Nishina - Cos(theta) v. cross section ')
plt.grid()
plt.show()


plt.plot((180/np.pi)*np.arccos(x),E(x,662))   
plt.xlabel('Scattering angle')
plt.ylabel('Scattered Photon energy | keV')
plt.legend(['Scattering angle v. Photon Energy'])
plt.title('Scattering angle v. Photon Energy')
plt.grid()
plt.show()


plt.plot(E(x,662),y) 
plt.xlabel('Scattered photon energy')
plt.ylabel('Cross-section')
plt.legend(['Scattered photon Energy vs Cross-Section'])
plt.title('Klein-Nishina Scattered photon Energy v. Cross section')
plt.grid()
plt.show()

plt.plot(612-E(x,662),y)
plt.xlabel('Scattered electron energy')
plt.ylabel('Cross-section')
plt.legend(['Scattered electron Energy vs Cross-Section'])
plt.title('Klein-Nishina Scattered electron Energy v. Cross section')
plt.grid()
plt.show()

#%%
# th is costheta
def sec(th,Eo):      # this formula to go in lab report 
    a=1/(1+Eo/511*(1-th))
    return 0.5*(2.8179e-15)**2*(a)**2*((a+1/a)-(1-th**2))/1.3350666520239514e-29 #normalised

def E(t,EE):
    return EE/(1+EE/511*(1-t))        
  
norm=spi.quad(sec,-1,1,args=662)
x=sp.arange(-1,1,0.01)

En=sp.arange(min(E(x,662)),max(E(x,662)),0.1)
theta=1-(662/En-1)*511/662   # Costheta - Comes from Eq a defined above in sec but...
                              #May be necessary to draw Geometry in lab report 
y=sec(theta,662)

#plt.plot(En,y)

ener=662
#create a distribution from 0 to 1000 for...
# Scattered Photon Energies..
dispa=sp.arange(0,1000,0.1)     # 0.1 Scaling - Dist. p for a (incoming Energies) 
newdist=[0]*len(dispa)       # newdist - Dist of Photon Cross-section after Scattering
for i in range(int(min(E(x,662))/0.1),int(max(E(x,662))/0.1)): # 0.1 Scaling
    #newdist[i]+=y[i-int(min(E(x,662))/0.1)] 
    if min(E(x,ener))<=dispa[i]<max(E(x,ener)):
        newdist[i]+=y[i-int((min(E(x,ener)))/0.1)] 



                           
'''      
dispa=sp.arange(-300,1000,1)
newdist=[0]*len(dispa)
for i in range(len(dispa)):
    if min(E(x,ener))<=dispa[i]<max(E(x,ener)):
        newdist[i]+=y[i-int((min(E(x,ener))+300))-1] 
'''   
     
#plt.plot(dispa,newdist) #scattered photon Energies vs Cross section

                    #plt.plot(662-dispa,newdist)

#scattered electron distribution

disnew=[]   # Rev Counting from 1000 to 0 APPEND Inverse Energies - List
            # For Electron Energy Distributions (e inverse of photons Energies)
            # As sumtotal Energy equals 662

#Flipping the Cross section list 
for i in range(len(newdist)):
    disnew.append(newdist[len(newdist)-i-1])
    
print ('len dispa vs len disnew ', len(dispa), len(disnew) )
#plt.plot(dispa,disnew)


# Counting downwards from..
# Points where Cross-section is Non-zero
diste0=[]    # Include all numbers but 0 
for i in range(len(disnew)):
    if disnew[i]!=0:
        diste0.append(disnew[i])
for i in range(len(dispa)-len(diste0)):
    diste0.append(0)

diste=[0]*300
for i in range(len(diste0)-300):
    diste.append(diste0[i])
print('len of dist. of e- Energies and dispa resp. are ', len(diste),len(dispa))





plt.plot(dispa,newdist) # Scattered Photon Cross-section
plt.plot(dispa,sp.array(diste))  # Scattered Electron Cross-section
plt.xlabel('Incident photon Energy, MeV')
plt.ylabel('Klein-Nishina Cross section')
plt.legend(['Incident photon Energy vs Cross-Section'])
plt.title('Klein-Nishina Incident Photon E v. Cross-section')
plt.savefig('Compton SHELF')
plt.grid()
plt.show()


plt.plot(dispa,nor(dispa,11.5*3.53-17,(11.5*3.53-17)**1.5/343,1,0,0))  #Characteristic iodine Peak
plt.xlabel('Incident photon Energies, MeV')
plt.ylabel('Frequency')
plt.legend(['Incident photon Energies - Iodine / (MeV)'])
plt.title('Characteristic Iodine Peak - Normal Dist.')
plt.grid()
plt.show()


plt.plot(dispa,nor(dispa,28.5*3.53-17,(28.5*3.53-17)**1.5/343,1,0,0))   #Unknown Peak
plt.xlabel('Incident photon Energy, MeV')
plt.ylabel('Frequency')
plt.legend(['Incident photon Energy - Unknown'])
plt.title('Unknown Peak - Normal Dist.')
plt.grid()
plt.show()


plt.plot(dispa,nor(dispa,662,662*0.075,1.5,0,0)) #Original Monte-Carlo Normal Dist. simulation
#plt.errorbar(A*3.53-17,cCs-cno,fmt='yo',yerr=sp.sqrt(np.abs(cCs-cno)))
plt.xlabel('Incident photon Energy, MeV')
plt.ylabel('Frequency')
plt.legend(['Incident photon Energy - Cs 137 Peak'])
plt.title('Predominant Cs-137 Peak - Normal Dist.')
plt.grid()
plt.show()


# Find out why it is multiplies by exp and Reason for...
#Values within Exponent - is it for smearing?
# Check again reason for raised power by 1.5/343 - try sqroot(mean)
backscadel=sp.array(newdist)+sp.array(diste)*sp.exp(0.001*dispa-1)+nor(dispa,662,662**1.5/343,0.0,0,0)


plt.plot(dispa,backscadel*2000)
plt.xlabel('Energy, MeV')
plt.ylabel('Frequency')
plt.legend(['Cs 137 - Compton Analysis Compton Scatter Model']) # at 0 degree measurement
#plt.title('Model Backscattering Cs 137 - Compton Analysis Backscatter Model.')
plt.title('Backscadel Model')
plt.grid()
plt.show()

#normpeak1=nor(dispa,662,res*Ai/2/0.512,2,0,0)


# Convolution Not used ! ! !
#cov=np.convolve(nor(dispa,662,662*0.075,1,0,0),sp.array(newdist)+sp.array(diste)*sp.exp(0.001*dispa-1)+nor(dispa,662,662*0.075,1,0,0))

#plt.plot(sp.linspace(0,1000,len(cov)),cov)

#cov1=np.convolve(nor(dispa,11.5,11.5**1.5/343,2,0,0),cov)
#plt.plot(sp.linspace(0,1000,len(cov1)),cov1)
#plt.plot(dispa,sp.array(newdist)*sp.array(newdiste)*10)

'''
plt.plot(x,y)   


plt.xlabel('Scattering angle')
plt.ylabel('Cross section')
plt.grid()

plt.show()

plt.plot(E(x,612),y) 
plt.xlabel('Scattered photon energy')
plt.ylabel('Cross section')
plt.grid()

plt.show()

plt.plot(612-E(x,612),y)
plt.xlabel('Scattered electron energy')
plt.ylabel('Cross section')
plt.grid()

'''

#%%


for i in range(len(diste[3000:])):
    if diste[3000:][i]==0:
        cut=i+3000-1
        print('diste Zero point is ', i+3000)
        break
      
#cutting=dispa[777] 
cutting=dispa[cut]
print('cuttoff photon energy-value dispa[777] point is, ', cutting)

newdistex=[]

for i in range(len(newdist)):
    #if dispa[i]<=477:
    if dispa[i]<=cutting:
        newdistex.append(newdist[i])
    else:
        newdistex.append(0)
      
#plt.plot(dispa,newdistex)

count=0

k = 3000
for i in range(len(newdistex[:k])):
    if newdistex[:k][i]==0:
        count+=1

print('newdistex is extension of Cross-section Newdist-Photons')        
print('count of Zeros in newdistex is ',count)
print('dispa count is index (from newdistex Zeros count is)', dispa[count]  )      

disteex=[]
for i in range(len(diste)):
    if dispa[i]<=dispa[count]:  # if Photon Energy is less than intro cutoff point..# append Zero
        disteex.append(0)
    else:
        disteex.append(diste[i]) # Else append the original Diste Cross-section for electron
 

# Heaviside Step-Function
# What are the Reasons for the Values 0.5, 0.25 and 0 - could merely state as an Observation
shift=[]
n = 0.000
t = 500
m = 602
for i in range(len(dispa)):
    
    if t<dispa[i]<=m:
        shift.append(0.055-n*((dispa[i]-t))**2)
    else:
        shift.append(0)
 


Photopeak1 = (nor(dispa,662,662*0.075,1.4,0,0)).tolist()
Photopeak = []
for i in range(len(dispa)):
    if m<dispa[i]:
        Photopeak.append(Photopeak1[i])
    else:
        Photopeak.append(0)
        

shift2=[]
n = 0.0000
o = 102
p = 170
for i in range(len(dispa)): 
    if o<dispa[i]<=p:
        shift2.append(0.1-(n*(dispa[i]-p))**2)
    else:
        shift2.append(0)             
shift2 = np.array(shift2)



Pb = nor(dispa,77.66,3*77.67**0.5,0.435,0,0)
Pb0 = []
for i in range(len(dispa)):
    if dispa[i]<o:
        Pb0.append(Pb[i])
    else:
        Pb0.append(0)
Pb0=np.array(Pb0)

print('len photopeak', len(Photopeak1))
print('len shift', len(shift))
Photopeak = np.array(Photopeak)


Pb1 = nor(dispa,75,3*75**0.5,0.435,0,0)
Pb2 = nor(dispa,73,75**0.5,0.435,0,0)
Pb3 = nor(dispa,85,75**0.5,0.435,0,0)

# Why is ar multiplicative-height 0.5 ?
# What do newdist and diste represent
# Why are there 2 Steps in the Step-Function?
# Why do we cuttoff Newdist and diste at the points where we cut them off?
# Is Backscadel ever used again?
#sim=sp.array(newdistex)+sp.array(disteex)+nor(dispa,11.5*3.53-17,(11.5*5)**1.5/343,1,0,0)+nor(dispa,27.01*3.53-17,(27.01*3.53-17)**1.5/343,0.5,0,0)+sp.array(shift)
#sim3 = (sp.array(newdistex))+(sp.array(disteex))+nor(dispa,11.5*3.53-17,(11.5*5)**1.5/343,0.68,0,0) + nor(dispa,662,662*0.075,1.4,0,0)
sim2_unscaled = shift2 + sp.array(shift) + (sp.array(newdistex))+(sp.array(disteex))+nor(dispa,11.5*3.53-17,(11.5*3.53-17)**0.5,0.68,0,0) + Photopeak + Pb0 +sp.array(shift)
sim2 = 1000*sim2_unscaled

plt.plot(dispa,1000*np.array(newdistex))
plt.plot(dispa,1000*np.array(disteex))
plt.plot(dispa,sim2)
plt.xlabel('Energy, KeV')
plt.ylabel('Frequency')
plt.legend(['Photon E-Cross-Section','Electron E-Cross-Section','Cs 137 - Compton Analysis Compton Scatter Model']) # at 0 degree measurement
plt.title('Model Backscattering Cs 137 - Compton Analysis Backscatter Model.')
plt.grid()
plt.show()

nK = 380
nq = 12000
plt.errorbar(A[0:nK]*3.53-17,(cCs-cno)[0:nK],fmt='yo',yerr=sp.sqrt(np.abs((cCs-cno)[0:nK])))
plt.plot(A[0:nK]*3.53-17,(cCs-cno)[0:nK],'yo')
plt.plot(dispa[0:nq],(sim2)[0:nq])
plt.xlabel('Energy, KeV')
plt.ylabel('Frequency')
plt.legend(['Measured Data','Simulated Model Cs 137']) # at 0 degree measurement
plt.title('Data vs Monte-Carlo Generated Model for Cs 137 - Compton Analysis')
plt.savefig('COMPTON MODEL - Montecarlo Theory Plot')
plt.grid()
plt.show()




#%% 

#Monte Carlo Numerical - Immanuel

# Compton Shelf

E0=662 # keV#input photon energy
#target
erm=511 #keV #electron rest mass
sa=0 #degree #scattering angle
#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution 

#ideal signal

def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))


Ei=sp.array(dispa)#ideal energy
Ai=Ei*gain # ideal480 analog signals
num=105
ampd1=sp.array([0]*num)

r = 0.1 # 90/14 #21/90 # 14/90 * 1.5

#for l in range(300,len(dispa)):
print('len dispa ', len(dispa))
#the number of points
for l in range(98,610):     # Check here for 10 factor
    l = int(l*10)
    #print('running ', l)
    n=int(sim2[l]*r) #
    #print('sim l is ', sim2[l])
    #print('n is ', n)
    # random seeds
    randseed=[]
    for i in range(0,n):
        randseed.append(random.uniform(0, 1))
    # noise
    #s=res*Ai[l]/2
    #s=(1.26683723*Ei[l]**0.5+1.92435425)*gain  # What are these?
    s = (Ei[l])**0.5*gain
    #s = 2*gain
    #s=  Ei[l]**5
    noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))

    # outcome bin
    Ob=(Ai[l]+noise)/maxsig*2**(bitdepth)
    Obf=sp.floor(Ob)
    #print(Obf)

    #outcome set
    #count the amplitude distribution
    num=105 #number of intervals
    ampdin1 = sp.array([0]*num)
    ampdin1 = ampdin1.tolist()
    index_Obf1 = ampdin1.index(max(ampdin1))
    for i in range(0,n):
        for j in range(0,num):
            if 5*j<=Obf[i]<5*j+5:
                ampd1[j]+=1
    #for j in range(0,num):
    #    if 5*j<=index_Obf1<5*j+5:
    #        ampd1[j]+=1
    #l = int(l/10)-200
    #ampd1[l]+=A1




erm=511 #keV #electron rest mass
sa=0 #degree #scattering angle
#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

               
#Peaks 
def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))

Ep1 = sp.array([int((11.5*3.53)-17), 74, 85,  662])     
#Ep  = 10*E(Ep1*0.1)      
Ai=Ep1*gain
           
ampd=sp.array([0]*num)
#the number of points
    #Cgood = (cCs-cno)
for l in range(4):
    if l == 0:
        r = 2 #np.sqrt((sim2[int(Ep1[l])]/1400)*90/14)
        n= int(sim2[int(Ep1[l]*10)]*r) #
        #print('sim Ep1 l is ', sim2[int(Ep1[l]*9)])
        #print('n is ', n)
        # random seeds
        randseed=[]
        for i in range(0,n):
            randseed.append(random.uniform(0, 1))
            # noise
            #s=res*Ai[l]/2
            #s=(1.26683723*(Ep1[l])**0.5+1.92435425)*gain  # What are these?
            s = (Ep1[l])**0.5*gain
            noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
            # outcome bin
            Ob=(Ai[l]+noise)/maxsig*2**(bitdepth)
            Obf=sp.floor(Ob)
            #print(Obf)
            #outcome set
            #count the amplitude distribution
            num=105 #number of intervals
            for j in range(0,num):
                if 5*j<=Obf[i]<5*j+5:
                    ampd[j]+=1
    if l == 1:
        r = 1.4 #np.sqrt((sim2[int(Ep1[l])]/1400)*90/14)
        n=int(sim2[int(Ep1[l]*10)]*20) #
        #print('sim Ep1 l is ', sim2[int(Ep1[l]*9)])
        #print('n is ', n)
        # random seeds
        randseed=[]
        for i in range(0,n):
            randseed.append(random.uniform(0, 1))
            # noise
            #s=res*Ai[l]/2
            #s=(1.26683723*(Ep1[l])**0.5+1.92435425)*gain  # What are these?
            s = (Ep1[l])**0.5*gain
            noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
            # outcome bin
            Ob=(Ai[l]+noise)/maxsig*2**(bitdepth)
            Obf=sp.floor(Ob)
            #print(Obf)
            #outcome set
            #count the amplitude distribution
            num=105 #number of intervals
            for j in range(0,num):
                if 5*j<=Obf[i]<5*j+5:
                    ampd[j]+=1
    if l == 2:
        r = 1 #np.sqrt((sim2[int(Ep1[l])]/1400)*90/14)
        n=int(sim2[int(Ep1[l]*10)]*20) #
        #print('sim Ep1 l is ', sim2[int(Ep1[l]*9)])
        #print('n is ', n)
        # random seeds
        randseed=[]
        for i in range(0,n):
            randseed.append(random.uniform(0, 1))
            # noise
            #s=res*Ai[l]/2
            #s=(1.26683723*(Ep1[l])**0.5+1.92435425)*gain  # What are these?
            s = (Ep1[l])**0.5*gain
            noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
            # outcome bin
            Ob=(Ai[l]+noise)/maxsig*2**(bitdepth)
            Obf=sp.floor(Ob)
            #print(Obf)
            #outcome set
            #count the amplitude distribution
            num=105 #number of intervals
            for j in range(0,num):
                if 5*j<=Obf[i]<5*j+5:
                    ampd[j]+=1
    if l == 3:
        r = np.sqrt((sim2[int(Ep1[l])]/1400)*90/14)
        n=int(sim2[int(Ep1[l]*10)]) #
        #print('sim Ep1 l is ', sim2[int(Ep1[l]*9)])
        #print('n is ', n)
        # random seeds
        randseed=[]
        for i in range(0,n):
            randseed.append(random.uniform(0, 1))
            # noise
            #s=res*Ai[l]/2
            #s=(1.26683723*(Ep1[l])**0.5+1.92435425)*gain  # What are these?
            s = (Ep1[l])**0.5*gain
            noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
            # outcome bin
            Ob=(Ai[l]+noise)/maxsig*2**(bitdepth)
            Obf=sp.floor(Ob)
            #print(Obf)
            #outcome set
            #count the amplitude distribution
            num=105 #number of intervals
            for j in range(0,num):
                if 5*j<=Obf[i]<5*j+5:
                    ampd[j]+=1
        
        
        
ampd = 0.1 * ampd

#%%

#first fit
# keV#input photon energy

E0 = 662
#target
erm=511 #keV #electron rest mass
sa=0 #degree #scattering angle
#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

#ideal signal
def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))
Ei=E(E0) #ideal energy
Ai=Ei*gain # ideal analog signal
print('ideal energy',Ei)

 
#first set
# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))    
#print('randseed Final Monte Carlo is ', randseed)
 
# noise         # study inverse error function
s=res*Ai/2    # s may represent Triangle area under the signal 
print('s',s)
noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
#noise=invgauss.ppf(sp.array(randseed),s,loc=0,scale=1)

# outcome bin
Ob=(Ai+noise)/maxsig*2**(bitdepth)
Obf=sp.floor(Ob)
#print(Obf)


#outcome set


#count the amplitude distribution
num=105 #number of intervals
ampd2=sp.array([0]*num)

for i in range(0,n):
    for j in range(0,num):
        if 5*j<=Obf[i]<5*j+5:
            ampd2[j]+=1


plt.plot(sp.linspace(0,520,105)/0.512,((ampd)+ampd1),'o')
nK = 380
plt.errorbar(A[0:nK]*3.53-17,(cCs-cno)[0:nK],fmt='yo',yerr=sp.sqrt(np.abs((cCs-cno)[0:nK])))
plt.xlabel('Energy, KeV')
plt.ylabel('Frequency')
plt.legend(['Measured Data','Simulated Model Cs 137']) # at 0 degree measurement
plt.title('Immanuels Numerical Monte-Carlo Test')
plt.savefig('Immanuels Numerical Monte-Carlo Test')
plt.show()






#%%
#E0=ener=511 #main peak energy
E0 = 662
n=8000 #the number of points

 


def sec(th,Eo):
    a=1/(1+Eo/511*(1-th))
    return 0.5*(2.8179e-15)**2*(a)**2*((a+1/a)-(1-th**2))/1.3350666520239514e-29 #normalised

 

def E(t,EE):
    return EE/(1+EE/511*(1-t))        

 

ener=662

 

print('norm',spi.quad(sec,-1,1,args=ener))
x=sp.arange(-1,1.01,0.01)

 

En=sp.arange(min(E(x,ener)),max(E(x,ener)),1)
theta=1-(ener/En-1)*511/ener
y=sec(theta,ener)

 
print(max(E(x,ener)))


#create a distribution from 0 to 1000

dispa=sp.arange(-300,1750,1)
newdist=[0]*len(dispa)
for i in range(len(dispa)):
    if min(E(x,ener))<=dispa[i]<max(E(x,ener)):
        newdist[i]+=y[i-int((min(E(x,ener))+300))-1]    


#scattered electron distribution


#plt.plot(662-dispa,newdist)

 
disnew=[]
for i in range(len(newdist)):
    disnew.append(newdist[len(newdist)-i-1])
#plt.plot(dispa,disnew)
 

diste0=[]
for i in range(len(disnew)):
    if disnew[i]!=0:
        diste0.append(disnew[i])
for i in range(len(dispa)-len(diste0)):
    diste0.append(0)

    
diste=[0]*300
for i in range(len(diste0)-300):
    diste.append(diste0[i])

 
print(len(diste),len(dispa))



for i in range(len(diste[300:])):
    if diste[300:][i]==0:
        cut=i+300-1
        break
        
cutting=dispa[cut] 
print('cut point number dispa for newdist-photon cross-sections list ', cutting)

 

newdistex=[]

 

for i in range(len(newdist)):
    if dispa[i]<=cutting:
        newdistex.append(newdist[i])
    else:
        newdistex.append(0)
        
plt.plot(dispa,newdistex)

 

count=0   #Zeros count before photon cross-sections starts 

for i in range(len(newdistex[:500])):
    if newdistex[:500][i]==0:
        count+=1
        
print('count',count)
print(dispa[count]  )      


disteex=[]
for i in range(len(diste)):
    if dispa[i]<=dispa[count]:
        disteex.append(0)
    else:
        disteex.append(diste[i])
        
shift=[]
for i in range(len(dispa)):
    if dispa[i]<=cutting:
        shift.append(0.5)
    elif cutting<dispa[i]<=ener:
        shift.append(0.25)
    else:
        shift.append(0)
        
'''
plt.plot(En,y)*2000
plt.plot(dispa,sp.array(newdist))
plt.plot(dispa,sp.array(diste))
'''

#plt.plot(dispa,disteex)
sim=sp.array(newdistex)+sp.array(disteex)+nor(dispa,27.6,(27.6)**0.5,0.5,0,0)+nor(dispa,75,(75)**0.5,0.5,0,0)+sp.array(shift)
plt.plot(dispa,sim)
nK = 380
#plt.errorbar(A[0:nK]*3.53-17,(cCs-cno)[0:nK],fmt='yo',yerr=sp.sqrt(np.abs((cCs-cno)[0:nK])))
plt.xlabel('Energy, KeV')
plt.ylabel('Frequency')
plt.legend(['Measured Data','Simulated Model Cs 137']) # at 0 degree measurement
plt.title('Wangs Data vs Monte-Carlo Generated Model for Cs 137 - Compton Analysis')
plt.savefig('Wangs Montecarlo Model Theory Plot')
plt.show()


#first fit
# keV#input photon energy


#target
erm=511 #keV #electron rest mass
sa=0 #degree #scattering angle

 
#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution


#ideal signal

 
def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))


Ei=E(E0) #ideal energy
Ai=Ei*gain # ideal analog signal
print('ideal energy',Ei)

 


#first set

 

 

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))
    
#print('randseed Final Monte Carlo is ', randseed)
 
# noise         # study inverse error function
s=res*Ai/2    # s may represent Triangle area under the signal 
print('s',s)
noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
#noise=invgauss.ppf(sp.array(randseed),s,loc=0,scale=1)


# outcome bin
Ob=(Ai+noise)/maxsig*2**(bitdepth)
Obf=sp.floor(Ob)
#print(Obf)


#outcome set


#count the amplitude distribution
num=105 #number of intervals
ampd1=sp.array([0]*num)

 

for i in range(0,n):
    for j in range(0,num):
        if 5*j<=Obf[i]<5*j+5:
            ampd1[j]+=1
        
        
#print(ampd)
#print( sp.linspace(0,1000,201))   

 

#plt.errorbar(sp.linspace(0,520,105)/0.512,ampd1,fmt='o',yerr=sp.sqrt(ampd),capsize=3)
#plt.errorbar(A*3.53-17,cCs-cno,fmt='o',yerr=sp.sqrt(np.abs(cCs-cno)))
#for a<Obf<=b:

 
 

 

#

 

E0=662 # keV#input photon energy

 

#target
erm=511 #keV #electron rest mass
sa=45 #degree #scattering angle

 

#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

 

#ideal signal

 

def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))

 


Ei=sp.array(dispa)#ideal energy
Ai=Ei*gain # ideal analog signals
#print('ideal energy',Ei)

 


#first set

 
#Monte Carlo Numerical
ampd=sp.array([0]*num)
#the number of points
for l in range(300,len(dispa)):
    n=int(sim[l]*50) #
    # random seeds
    randseed=[]
    for i in range(0,n):
        randseed.append(random.uniform(0, 1))

    # noise
    #s=res*Ai[l]/2
    s=(1.26683723*Ei[l]**0.5+1.92435425)*gain  # What are these?
    noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))

    # outcome bin
    Ob=(Ai[l]+noise)/maxsig*2**(bitdepth)
    Obf=sp.floor(Ob)
    #print(Obf)

    #outcome set
    #count the amplitude distribution
    num=105 #number of intervals
    for i in range(0,n):
        for j in range(0,num):
            if 5*j<=Obf[i]<5*j+5:
                ampd[j]+=1




#plt.errorbar((A[:300]*3.53-17),(cCs[:300]-cno[:300]),fmt='o',yerr=sp.sqrt(np.abs(cCs[:300]-cno[:300]))) 
#plt.errorbar((A[:300]*3.53-17),(cNa[:300]-cno[:300]),fmt='o',yerr=sp.sqrt(np.abs(cNa[:300]-cno[:300]))) 

plt.plot(sp.linspace(0,520,105)/0.512,(ampd+ampd1)/6,'o')
nK = 380
plt.errorbar(A[0:nK]*3.53-17,(cCs-cno)[0:nK],fmt='yo',yerr=sp.sqrt(np.abs((cCs-cno)[0:nK])))
plt.xlabel('Energy, KeV')
plt.ylabel('Frequency')
plt.legend(['Measured Data','Simulated Model Cs 137']) # at 0 degree measurement
plt.title('Wangs Numerical Monte-Carlo Test')
plt.savefig('Wangs Numerical Monte-Carlo Test')
plt.show()

# Generated Model for Cs 137 - Compton Analysis

#plt.errorbar((A[:300]*3.53-17),(cMn[:300]-cno[:300]),fmt='o',yerr=sp.sqrt(np.abs(cMn[:300]-cno[:300]))) 
#plt.errorbar((A[:300]*3.53-17),(cAm[:300]-cno[:300]),fmt='o',yerr=sp.sqrt(np.abs(cAm[:300]-cno[:300]))) 

 

 

 


#print(ampd)
#plt.show()
#print( sp.linspace(0,1000,201))   

#plt.plot(sp.linspace(0,520,105)/0.512,ampd,'o')
#plt.plot(sp.array(gA3),gsCsp3,'o')
#for a<Obf<=b:


#%%

'''
#Monte Carlo Simulation 0

#source 
E0=662 # keV#input photon energy

#target
erm=511 #keV #electron rest mass
sa=0 #degree #scattering angle

#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

#ideal signal

def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))

Ei=E(E0) #ideal energy
Ai=Ei*gain # ideal analog signal
print('ideal energy',Ei)


#first set


n=9000 #the number of points

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))
    

# noise
s=res*Ai/2
print('s',s)
noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
#noise=invgauss.ppf(sp.array(randseed),s,loc=0,scale=1)


# outcome bin
Ob=(Ai+noise)/maxsig*2**(bitdepth)
Obf=sp.floor(Ob)
#print(Obf)


#outcome set

#count the amplitude distribution
num=105 #number of intervals
ampd1=sp.array([0]*num)

for i in range(0,n):
    for j in range(0,num):
        if 5*j<=Obf[i]<5*j+5:
            ampd1[j]+=1
        
        
#print(ampd)
#print( sp.linspace(0,1000,201))   

plt.errorbar(sp.linspace(0,520,105)/0.512,ampd1,fmt='o',yerr=sp.sqrt(ampd1),capsize=3)
plt.errorbar(A*3.53-17,cCs-cno,fmt='o',yerr=sp.sqrt(np.abs(cCs-cno)))
#for a<Obf<=b:


#fitting the curve
fitcCs,covcCs = curve_fit(nnor, A*3.53-17, cCs-cno,[600,5,1000,100])
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,1750,1000)
plt.plot(Ax,nnor(Ax,*fitcCs))


fitcCs1,covcCs1 = curve_fit(nnor, sp.linspace(0,520,105)/0.512,ampd,[660,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)

for i in range(0,4):
    print(covcCs1[i,i])
    
plt.plot(Ax,nnor(Ax,*fitcCs1))

plt.title('Frequency-energy graph at 0 degree')
plt.show()


#%%
#Monte Carlo Simulation 45
plt.plot(dispa,backscadel*200)
#source 
E0=662 # keV#input photon energy

#target
erm=511 #keV #electron rest mass
sa=45 #degree #scattering angle

#detector
gain=5e-3
bitdepth=9 #bits
maxsig=5
res=0.075 #resolution

#ideal signal

def E(EE):
    return EE/(1+EE/erm*(1-sp.cos(sa/180*sp.pi)))


Ei=sp.array(dispa)#ideal energy
Ai=Ei*gain # ideal analog signal
print('ideal energy',Ei)


#first set

ampd=sp.array([0]*num)
#the number of points
for l in range(300,len(dispa)):
    n=int(sim[l]*50) #
    # random seeds
    randseed=[]
    for i in range(0,n):
        randseed.append(random.uniform(0, 1))


    # noise
    s=res*Ai[l]/2
    noise=-2**0.5*s*sps.erfinv(1-2*sp.array(randseed))
    #noise=invgauss.ppf(sp.array(randseed),s,loc=0,scale=1)


    # outcome bin
    Ob=(Ai[l]+noise)/maxsig*2**(bitdepth)
    Obf=sp.floor(Ob)
    #print(Obf)


    #outcome set

    #count the amplitude distribution
    num=105 #number of intervals
    

    for i in range(0,n):
        for j in range(0,num):
            if 5*j<=Obf[i]<5*j+5:
                ampd[j]+=1

        
#print(ampd)
#print( sp.linspace(0,1000,201))   
plt.plot(sp.linspace(0,520,105)/0.512,ampd+ampd1,'o')
plt.show()
#plt.plot(sp.linspace(0,520,105)/0.512,ampd,'o')
#plt.plot(sp.array(gA3),gsCsp3,'o')
#for a<Obf<=b:

'''


"""
fitcCs,covcCs = curve_fit(nnor, sp.array(gA3),gsCsp3,[500,10,50,10])
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,1750,1000)
plt.plot(Ax,nnor(Ax,*fitcCs))


fitcCs1,covcCs1 = curve_fit(nnor, sp.linspace(0,520,105)/0.512,ampd,[500,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)

for i in range(0,4):
    print(covcCs1[i,i])
    
plt.plot(Ax,nnor(Ax,*fitcCs1))
plt.title('Frequency-energy graph at 45 degree')

"""



'''
plt.errorbar(A*3.53-17,cCs-cno,fmt='o',yerr=sp.sqrt(np.abs(cCs-cno)))
#plt.errorbar(sp.linspace(0,520,105)/0.512,ampd1,fmt='o',yerr=sp.sqrt(ampd),capsize=3)
'''
