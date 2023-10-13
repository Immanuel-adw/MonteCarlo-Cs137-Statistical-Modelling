# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 08:46:10 2021

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
#Amplitude to Energy conversion

#amplitude scale:
Amplitude=np.array([149.94,238.8,192.00,21.72])  #Array of Peak Amplitude Freq.
erramplitude=np.array([0.11,0.28,0.05,0.01])

#the emission peaks
energy=np.array([511,821,662,59.6])


#linear fit for converting amplitude to energy
fit = np.polyfit(Amplitude,energy,1,w=1/erramplitude)
pc=np.poly1d(fit)
T=np.linspace(0,250,1000)


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
#plt.savefig('4 '+name)
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


#group the data
gA1=[]
gsCsp1=[]
for i in range(0,len(A1)-8):
    gA1.append(sp.mean([A1[i],A1[i+1],A1[i+2],A1[i+3],A1[i+4],A1[i+5],A1[i+6],A1[i+7]]))
    gsCsp1.append(sp.mean([sCsp1[i],sCsp1[i+1],sCsp1[i+2],sCsp1[i+3],sCsp1[i+4],sCsp1[i+5],sCsp1[i+6],sCsp1[i+7]]))
plt.errorbar(gA1,gsCsp1,fmt='o',yerr=sp.array(gsCsp1)**0.5)

   

    
Ax=sp.linspace(0,1750,1000)


fitsCsV,covsCs = curve_fit(voi, A1, sCsp1,[ 4.19738061e+01, 4.14325173e+01,  3.31997291e+02,  5.94834892e+02, 8.53485283e-01, -5.99601714e-04], sigma = 1/(np.abs(sp.array(sCsp1))**0.5 +0.4))
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
#plt.savefig(name)
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
#plt.savefig(name)
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

plt.xlabel('Energy keV')
plt.ylabel('Frequency (Counts)')
plt.title(name2)
#plt.savefig(name)
plt.show()

    
print(0.5035778448149948**0.5+4.70)
print('the peak for 45 degree is at 470.21±5.41')



#%%

Eth=np.array([662,609.41,549.43,470.41])*1.6e-16
errEth=np.array([0.1,2.09,3.18,1.98])*1.6e-16
ang=np.array([0,20,30,45])


'''
Eth=np.array([662,596.04,549.865,475.25])*1.6e-16
errEth=np.array([0.15,0.838, 8.728,4.12])*1.6e-16
ang=np.array([0,20,30,45])
'''

fit = np.polyfit(np.cos(ang*np.pi/180),1/Eth,1,w=1/(errEth/Eth**2))
#fit, cov = np.polyfit(np.cos(ang*np.pi/180),1/Eth,1,w=1/(errEth))
p=np.poly1d(fit)
#[y, delta] = np.polyval(p,np.cos(ang*np.pi/180))
T=np.linspace(0.67,1.05,1000)
print('fit 1/E vs Costheta is ', fit)
#print('delta S is ', S)
#print('cov 1/E vs Costheta is ', cov1)

name = 'Best-fit line for Inverse of Peak-Energy vs Cos-theta - Cs plot'
plt.errorbar(np.cos(ang*np.pi/180),1/Eth,yerr=errEth/Eth**2,fmt='bo', ms=3, capsize = 5)
plt.plot(T,p(T), 'g-')
plt.legend(['Cos-theta vs 1/Peak-Energy', 'Best-fit line'])
plt.grid()
plt.xlabel('Cos-theta')
plt.ylabel('1/Peak-Energy | $keV^{-1}$')
plt.title('Verifying Compton-Scattering', size=10  )
plt.savefig('Verifying Compton-Scattering')
plt.show()

#plt.plot(np.cos(ang*np.pi/180),p(np.cos(ang*np.pi/180)))

dif=(1/Eth)-p(np.cos(ang*np.pi/180))
#print(dif)

var=0
for i in range(len(dif)):
    var+= fit *(dif[i]/p(np.cos(ang*np.pi/180))[i])**2
var = sum((errEth/Eth**2))   
print('percentage variance',var,'percentage sd',(abs(var)**0.5))


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
s= (Ei**0.5)*gain
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

#plt.errorbar(sp.linspace(0,520,105)/0.512,ampd,fmt='o',yerr=sp.sqrt(ampd),capsize=3)
#plt.errorbar(A*3.53-17,cCs-cno,fmt='o',yerr=sp.sqrt(np.abs(cCs-cno)))
#for a<Obf<=b:

#fitting the curve
fitcCs,covcCs = curve_fit(nnor, (A*3.53-17), (cCs-cno),[600,5,1000,100])
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,1750,1000)

fitcCs1,covcCs1 = curve_fit(nnor, sp.linspace(0,520,105)/0.512,ampd,[660,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)

for i in range(0,4):
    print(covcCs1[i,i])


plt.plot(Ax[0:600],nnor(Ax,*fitcCs)[0:600])    
plt.plot(Ax[0:600],nnor(Ax,*fitcCs1)[0:600])
plt.title('Cs Monte-Carlo Frequency-energy graph at 0 degree')
plt.legend(['Experiment-fit', 'Simulation-fit'])
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
r = 9000/1380
n= 47  #int(9*r)    #45 #the number of points

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))

    
 
# noise
s=res*Ai/2
s= (Ei**0.5)*gain
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


n= 40 #int(6*900/1400) #35 #the number of points

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))
 
      
# noise
s=res*Ai/2
s= (Ei**0.5)*gain
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


n= 90 #int(16.25*9000/1400) #90 #the number of points

# random seeds
randseed=[]
for i in range(0,n):
    randseed.append(random.uniform(0, 1))

  
# noise
#s=res*Ai/2
s= (Ei**0.5)*gain
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
