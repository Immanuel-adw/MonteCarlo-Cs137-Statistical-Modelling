#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 13:45:19 2021

@author: Ginny
"""

#%%
import matplotlib.pyplot as plt
import scipy as sp
import scipy.integrate as spi
import scipy.special as sps
from scipy.optimize import curve_fit
import scipy.optimize as op
import numpy as np

#%% import data
data=sp.loadtxt("Experiment_Data copy.csv",  skiprows=7,delimiter=',', unpack=True)

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

#%%  formatting

plt.style.use('seaborn-deep')
params = {   
    'axes.labelsize': 10,
   'font.size': 8,
   'legend.fontsize': 10,
   'xtick.labelsize': 9,
   'ytick.labelsize': 9,
   'figure.figsize': [8.9,7]
   } 
plt.rcParams.update(params)

#%%
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
sp    
#if the peak cannot be fitted with voigt profile.
    
#%% fit for Na peak
cNa0=cNa-cno
plt.plot(A,cNa0,'o')
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
     
     

guess_params_Na = [5,5,800,150,25,-0.5]
g = guess_params_Na
#s,g width, s normal distribution width, g width of the Lorentzian
#a height multiply
#u peak frequency
#m height shifting
#b gradient
       
fitcNa,covcNa = curve_fit(voi, Ap, cNap, g)
print(fitcNa)

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
    
Ax=np.linspace(0,500,1000)

name = 'Callibration Na0 Amplitude v freq plot'
plt.plot(Ap,cNap,'o') 
plt.plot(Ax,voi(Ax,*fitcNa))
plt.legend(['Callibration Na0 Amplitude v freq plot','Na0 Voigt fit plot'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('2 '+name)
plt.show()

# peak position 149.94±0.11
#%% fit for Mn peak
cMn0=cMn-cno
#plt.plot(A,cMn0,'o')

Ap=[]
cMnp=[]
for i in range(0,len(A)):
    if 190<A[i]<300:
        Ap.append(A[i])
        cMnp.append(cMn0[i])


guess_params_Mn = [5,5,50,250,5,-0.5]
g = guess_params_Mn
#s,g width, s normal distribution width, g width of the Lorentzian
#a height multiply
#u peak frequency
#m height shifting
#b gradient
      
fitcMn,covcMn = curve_fit(voi, Ap, cMnp, g)
print(fitcMn)

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
    
Ax=np.linspace(0,500,1000)

name = 'Callibration Mn0 Amplitude v freq plot'   
plt.plot(Ap,cMnp,'o')
plt.plot(Ax,voi(Ax,*fitcMn))
plt.legend(['Callibration Mn0 Amplitude v freq plot','Mn0 Voigt fit Amplitude vs. freq(A) plot'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('3 '+name)
plt.show()
# peak for Mn at 238.80±0.08

#%% fit for Cs peak
cCs0=cCs-cno
plt.plot(A,cCs0,'o')
plt.show()

Ap=[]
cCsp=[]
for i in range(0,len(A)):
    if 150<A[i]<220:
        Ap.append(A[i])
        cCsp.append(cCs0[i])

   

guess_params_Mn = [190,5,1000,100,-0.5]
g = guess_params_Mn
#u peak frequency
#s width, s normal distribution width,
#a height multiply
#m height shifting
#b gradient

   
fitcCs,covcCs = curve_fit(nor, Ap, cCsp,[190,5,1000,100,-0.5])
print(fitcCs)

for i in range(0,4):
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
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    
Ax=np.linspace(0,500,1000)

name = 'callibration Cs0 Amplitude v freq(A) plot'
plt.plot(Ap,cCsp,'o')
plt.plot(Ax,nor(Ax,*fitcCs))
plt.legend(['callibration Cs0 Amplitude v freq(A) plot', 'Cs0 Voigt fit'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('4 '+name)
plt.show()
# peak for Cs at 192.00±0.02

#%% fit for Am peak
cAm0=cAm-cno
plt.plot(A,cAm0,'o')
plt.show()

Ap=[]
cAmp=[]
for i in range(0,len(A)):
    if 13<A[i]<50:
        Ap.append(A[i])
        cAmp.append(cAm0[i])
    
        
fitcAm,covcAm = curve_fit(nor, Ap, cAmp,[25,5,8000,10,0.05])
print(fitcAm)

#u peak frequency
#s width, s normal distribution width,
#a height multiply
#m height shifting
#b gradient
for i in range(0,4):
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
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
        
        
Ax=np.linspace(0,500,1000)

name = 'Callibration Am0 Amplitude v freq(A) plot'   
plt.plot(Ap,cAmp,'o')
plt.plot(Ax,nor(Ax,*fitcAm))
plt.legend(['Callibration Am0 Amplitude v freq(A) plot','Am0 Voigt fit plot'])
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('5 '+name)
plt.show()

# peak for Mn at 21.72±0.02

#%%
#fitting for the conversion

#amplitude scale:
Amplitude=np.array([149.94,238.8,192.00,21.72])  #Array of Peak Amplitude Freq.
erramplitude=np.array([0.11,0.28,0.05,0.01])

#the emission peaks
energy=np.array([511,821,662,59.6])



#linear fit for converting amplitude to energy
fit = np.polyfit(Amplitude,energy,1,w=1/erramplitude)
p=np.poly1d(fit)
T=np.linspace(0,250,1000)

name =  'Peak Amplitude Frequencies Vs Known Peak Energies of Elements'
plt.errorbar(Amplitude, energy, xerr=erramplitude,fmt='o')
plt.plot(T,p(T))
plt.legend(['Peak Amplitude Frequencies Vs Known Peak Energies of Elements','Best-fit Amplitude vs Energy Correlation'])
plt.xlabel('Amplitude')
plt.ylabel('Energy')
plt.savefig('6 '+name)
plt.show()

print('Amplitude Energy Gradient is ', fit)

expenergy=p(Amplitude)


#average perventage error
delta=expenergy-energy
print('Delta Expected Energy - Measured is ', delta)

a=0
for i in range(len(delta)):
    a+=(delta[i]/energy[i])**2
    
#print(a)
toerror=a**0.5
print('total percentage error from fitting the points',toerror)

plt.grid()
plt.show()

#%%
#plt.plot(A,cCs/10)

AA=p(A)     #  Scintillator Values Calibrated to represent Energies in MeV
plt.plot(AA,sCs1,'o')



A1=[]
sCsp1=[]
for i in range(0,len(A)):
    if 450<AA[i]<750:
        A1.append(AA[i])
        sCsp1.append(sCs1[i])
        
name = 'Energy vs FreqCount 20deg Cs plot'
plt.plot(A1,sCsp1,'o')
plt.legend(['Energy vs FreqCount 20deg Cs plot'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.savefig('7 '+name)
plt.show()

#group the data - Creating Bins of 8 
gA1=[]
gsCsp1=[]
for i in range(0,len(A1)-8):
    gA1.append(np.mean([A1[i],A1[i+1],A1[i+2],A1[i+3],A1[i+4],A1[i+5],A1[i+6],A1[i+7]]))
    gsCsp1.append(np.mean([sCsp1[i],sCsp1[i+1],sCsp1[i+2],sCsp1[i+3],sCsp1[i+4],sCsp1[i+5],sCsp1[i+6],sCsp1[i+7]]))



fitsCs,covsCs = curve_fit(nor, gA1, gsCsp1,[600,10,50,10,0.0])
print(fitsCs)

for i in range(0,4):
    element = 'Cs 20deg'
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
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    
Ax=np.linspace(0,1750,1000)

 
name = 'Energy vs FreqCount 20deg Cs plot - bins of 8'
plt.plot(gA1,gsCsp1,'o')       
plt.plot(Ax,nor(Ax,*fitsCs))
plt.legend(['Energy vs FreqCount 20deg Cs plot - bins of 8', 'Gauss-Fit - Energy vs FreqCount 20deg Cs plot'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.savefig('8 '+name)
plt.show()

#plt.plot(A,sCs2,'o')
#plt.plot(A,sCs3,'o')
#plt.legend(['20','30','45'])

print(0.49471538412305943**0.5+6.10)
print('the peak for 20 degree is at 610.01±6.80')

#%%
plt.plot(AA,sCs2,'o')
plt.legend(['Energy vs FreqCount 30deg Cs plot - untrimmed x-values'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.show()

A2=[]
sCsp2=[]
for i in range(0,len(AA)):
    if 460<AA[i]<1750:
        A2.append(AA[i])
        sCsp2.append(sCs2[i])

plt.plot(A2,sCsp2,'o')
plt.legend(['Energy vs FreqCount 30deg Cs plot'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.show()

gA2=[]
gsCsp2=[]
for i in range(0,len(A2)-8):
    gA2.append(np.mean([A2[i],A2[i+1],A2[i+2],A2[i+3],A2[i+4],A2[i+5],A2[i+6],A2[i+7]]))
    gsCsp2.append(np.mean([sCsp2[i],sCsp2[i+1],sCsp2[i+2],sCsp2[i+3],sCsp2[i+4],sCsp2[i+5],sCsp2[i+6],sCsp2[i+7]]))

        
fitsCs2,covsCs2 = curve_fit(nor, gA2, gsCsp2,[570,5,20,0,-0.05])
print(fitsCs2)

for i in range(0,4):
    element = 'Cs 30deg'
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
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    
Ax=np.linspace(0,1750,1000)

name =  'Energy vs FreqCount 30deg Cs plot - bins of 8'
plt.plot(gA2,gsCsp2,'o')
plt.plot(Ax,nor(Ax,*fitsCs2))  
plt.plot(A2,sCsp2,'o')
plt.legend(['Energy vs FreqCount 30deg Cs plot - bins of 8', 'Gaussian Dist. Fit'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.savefig('9 '+name)
plt.show()

#plt.plot(A,sCs2,'o')
#plt.plot(A,sCs3,'o')
#plt.legend(['20','30','45'])

print(0.5842594969439249**0.5+5.49)
print('the peak for 30 degree is at 549.13±6.25')
#%%
plt.plot(AA,sCs3,'o')


A3=[]
sCsp3=[]
for i in range(0,len(AA)):
    if 390<AA[i]<1710:
        A3.append(AA[i])
        sCsp3.append(sCs3[i])
plt.plot(A3,sCsp3,'o')
plt.legend(['Energy vs FreqCount 45deg Cs plot - untrimmed x-values'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.show()
        

plt.plot(A3,sCsp3,'o')
plt.legend(['Energy vs FreqCount 45deg Cs plot'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.show()


gA3=[]
gsCsp3=[]
for i in range(0,len(A3)-8):
    gA3.append(np.mean([A3[i],A3[i+1],A3[i+2],A3[i+3],A3[i+4],A3[i+5],A3[i+6],A3[i+7]]))
    gsCsp3.append(np.mean([sCsp3[i],sCsp3[i+1],sCsp3[i+2],sCsp3[i+3],sCsp3[i+4],sCsp3[i+5],sCsp3[i+6],sCsp3[i+7]]))
 
        

     
fitsCs3,covsCs3 = curve_fit(nor, gA3, gsCsp3,[500,3,20,0,-0.05])
print(fitsCs3)

for i in range(0,4):
    element = 'Cs 45deg'
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
        print('cov ' + element + ' Gradient is ' + cov)
        print('err ' + element + ' Gradient is ' + err)
    
Ax=np.linspace(0,1750,1000)


name =  'Energy vs FreqCount 45deg Cs plot - bins of 8'
plt.plot(A3,sCsp3,'o')
plt.plot(gA3,gsCsp3,'o')
plt.plot(Ax,nor(Ax,*fitsCs3))  
plt.legend(['Energy vs FreqCount 45deg Cs plot','Energy vs FreqCount 45deg Cs plot - bins of 8', 'Gaussian Dist. Fit'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.savefig('10 '+name)
plt.show()

#plt.plot(A,sCs2,'o')
#plt.plot(A,sCs3,'o')


print(0.5035778448149948**0.5+4.70)
print('the peak for 45 degree is at 470.21±5.41')
#%%
Eth=np.array([662,609.41,549.43,470.41])*1.6e-16
errEth=np.array([0.1,2.09,3.18,1.98])*1.6e-16
ang=np.array([0,20,30,45])



fit= np.polyfit(np.cos(ang*np.pi/180),1/Eth,1,w=1/(errEth/Eth**2))
p=np.poly1d(fit)
T=np.linspace(0.67,1.05,1000)

name = 'Best-fit line for Inverse of Peak-Energy vs Cos-theta - Cs plot'
plt.errorbar(np.cos(ang*np.pi/180),1/Eth,yerr=errEth/Eth**2,fmt='o')
plt.plot(T,p(T))
print(fit)
plt.legend(['Cos-theta vs 1/Peak-Energy Cs plot', 'Best-fit line for Cos-theta vs 1/Peak-Energy Cs plot'])
plt.xlabel('Cos-theta')
plt.ylabel('1/Peak-Energy')
plt.savefig('11 '+name)
plt.show()

#plt.plot(np.cos(ang*np.pi/180),p(np.cos(ang*np.pi/180)))

dif=1/Eth-p(np.cos(ang*np.pi/180))
print(dif)

var=0
for i in range(len(dif)):
    var+=dif[i]**2
    
print('variance',var,'sd',var**0.5)
#%%
