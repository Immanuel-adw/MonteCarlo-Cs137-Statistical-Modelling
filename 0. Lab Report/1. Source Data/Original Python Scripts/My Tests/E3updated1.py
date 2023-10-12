#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 11:50:33 2021

@author: Ginny
"""
#%%

import matplotlib.pyplot as plt
import scipy as sp
import scipy.integrate as spi
import scipy.special as sps
import random
from scipy.stats import invgauss
import numpy as np

#%%
data=sp.loadtxt("/Users/apple/Desktop/Lab E3/callibration py/Experiment_Data copy.csv",  skiprows=7,delimiter=',', unpack=True)

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
    'axes.labelsize': 40,
   'font.size': 20,
   'legend.fontsize': 24,
   'xtick.labelsize': 24,
   'ytick.labelsize': 24,
   'figure.figsize': [20,12]
   } 
plt.rcParams.update(params)

#%%

#overall plot for peaks
#plt.plot(A,cno,'o')
plt.plot(A,cNa-cno,'o')
plt.plot(A,cMn-cno,'o')
#plt.plot(A,cCs,'o')
#plt.plot(A,cAm,'o')
plt.legend(['no','Na','Mn','Cs','Am'])

#%%

#define voigt profile
def voi(x,s,g,a,u,m,b):
    z=((x-u)+1j*g)/(s*sp.sqrt(2))
    w=sp.exp(-z**2)*sps.erfc(-1j*z)
    return a*sp.real(w)/(s*sp.sqrt(2*sp.pi))+m+b*x
x=sp.linspace(-10,10,1000)
plt.plot(x,voi(x,1,1,1,1,1,0.001))

#u peak
#s,g width, s normal distribution width, g width of the Lorentzian
#a height multiply
#m height shifting

#define normal distribution
def nor(x,u,s,a,m,b):
    return a*sp.exp(-(x-u)**2/s**2)+m+b*x

def nnor(x,u,s,a,m):
    return a*sp.exp(-(x-u)**2/s**2)+m
    

#%%
    
# fit for Na peak
cNa0=cNa-cno
plt.errorbar(A,cNa0,fmt='o',yerr=(sp.array(cNa0)**0.5),capsize=3)
    
plt.plot(A,cNa0,'o')

Ap=[]
cNap=[]
for i in range(0,len(A)):
    if 120<A[i]<180:
        Ap.append(A[i])
        cNap.append(cNa0[i])
        
plt.errorbar(Ap,cNap,fmt='o',yerr=(sp.array(cNap)**0.5),capsize=3)
    
        
fitcNa,covcNa = curve_fit(voi, Ap, cNap,[5,5,800,150,25,-0.5],sigma=1/(sp.array(cNap)**0.5))
print(fitcNa)

for i in range(0,6):
    print(covcNa[i,i]) #covariance, error is estimated by standard deviation
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,voi(Ax,*fitcNa))

fitcNa,covcNa = curve_fit(nor, Ap, cNap,[150,5,500,0,0])
print(fitcNa)

for i in range(0,5):
    print(covcNa[i,i]) #covariance, error is estimated by standard deviation
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,nor(Ax,*fitcNa))



plt.grid()
plt.xlabel('amplitude',size=30)
plt.ylabel('frequency',size=30)
plt.title('Na peak',size=40)

# peak position 149.87±0.13

#%%

# fit for Mn peak
cMn0=cMn-cno
plt.errorbar(A,cMn0,fmt='o',yerr=np.abs((cMn0))**0.5,capsize=3)

Ap=[]
cMnp=[]
for i in range(0,len(A)):
    if 185<A[i]<300:
        Ap.append(A[i])
        cMnp.append(cMn0[i])
        
plt.errorbar(Ap,cMnp,fmt='o',yerr=np.abs(sp.array(cMnp))**0.5,capsize=3)

errcC=[]
for i in range(len(cMnp)):
    if cMnp[i]>=20:
        errcC.append(cMnp[i]**0.5)
    else:
        errcC.append(0)
print(errcC)
        
        
fitcMn,covcMn = curve_fit(voi, Ap, cMnp,[8.49590479e+00,-2.36642502e+00,4.47826625e+02,2.38801836e+02,1.32477875e+01,-4.66356655e-02],sigma=1/(np.abs(sp.array(cMnp))**0.5))
print(fitcMn)

for i in range(0,6):
    print(covcMn[i,i])
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,voi(Ax,*fitcMn))


fitcMn,covcMn = curve_fit(nor, Ap, cMnp,[240,5,50,0,0],sigma=1/(np.abs(sp.array(cMnp))**0.5))
print(fitcMn)

for i in range(0,5):
    print(covcMn[i,i])
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,nor(Ax,*fitcMn))

# peak for Mn at 238.80±0.28


plt.xlabel('amplitude',size=30)
plt.ylabel('frequency',size=30)
plt.title('Mn peak',size=40)
plt.grid()

#%%
# fit for Cs peak
cCs0=cCs-cno
plt.errorbar(A,cCs0,fmt='o',yerr=cCs0**0.5,capsize=3)

Ap=[]
cCsp=[]
for i in range(0,len(A)):
    if 150<A[i]<220:
        Ap.append(A[i])
        cCsp.append(cCs0[i])
        
plt.errorbar(Ap,cCsp,fmt='o',yerr=sp.array(cCsp)**0.5,capsize=3)

errcC=[]
for i in range(len(cCsp)):
    if cCsp[i]>=1000:
        errcC.append(cCsp[i]**0.5)
    else:
        errcC.append(0)
print(errcC)
        
        
        
fitcCs,covcCs = curve_fit(voi, Ap, cCsp,[ 5.70486142e+00 , 2.24738185e+00 , 2.53071838e+04,  1.92045165e+02,
  5.26171410e+02, -2.64437092e+00],sigma=1/((sp.array(errcC)**0.5+1)))
print(fitcCs)


for i in range(0,6):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,voi(Ax,*fitcCs))

fitcCs,covcCs = curve_fit(nor, Ap, cCsp,[192,5,1300,0,0],sigma=1/((sp.array(errcC)**0.5+1)))
print(fitcCs)

for i in range(0,5):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,nor(Ax,*fitcCs))

# peak for Cs at 192.00±0.05

plt.xlabel('amplitude',size=30)
plt.ylabel('frequency',size=30)
plt.title('Cs peak',size=40)
plt.grid()

#%%

# fit for Am peak
cAm0=cAm-cno
plt.errorbar(A,cAm0,fmt='o',yerr=cAm0**0.5,capsize=3)

Ap=[]
cAmp=[]
for i in range(0,len(A)):
    if 13<A[i]<50:
        Ap.append(A[i])
        cAmp.append(cAm0[i])
        
plt.errorbar(Ap,cAmp,fmt='o',yerr=sp.array(cAmp)**0.5,capsize=3)

errcC=[]
for i in range(len(cAmp)):
    if cAmp[i]>=1000:
        errcC.append(cAmp[i]**0.5)
    else:
        errcC.append(0)
print(errcC)
        
        
fitcAm,covcAm = curve_fit(voi, Ap, cAmp,[ 1.50901365e+00,  8.91068122e-02,  3.51782797e+04,  2.17232445e+01,
  2.45872987e+01, -6.67251889e-01],sigma=1/((sp.array(errcC)**0.5+1)))
print(fitcAm)

for i in range(0,4):
    print(covcAm[i,i])
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,voi(Ax,*fitcAm))


fitcAm,covcAm = curve_fit(nor, Ap, cAmp,[22,5,9000,0,0],sigma=1/((sp.array(errcC)**0.5+1)))
print(fitcAm)

for i in range(0,5):
    print(covcAm[i,i])
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,nor(Ax,*fitcAm))



# peak for Mn at 21.72±0.01
plt.xlabel('amplitude',size=30)
plt.ylabel('frequency',size=30)
plt.title('Am peak',size=40)
plt.grid()
#%%

#fitting for the conversion

#amplitude scale:
Amplitude=sp.array([149.94,238.8,192.00,21.72])
erramplitude=sp.array([0.13,0.31,0.06,0.02])

#the emission peaks
energy=sp.array([511,821,662,59.6])

plt.errorbar(Amplitude, energy, xerr=erramplitude,fmt='o',capsize=3)

#linear fit for converting amplitude to energy
fit = sp.polyfit(Amplitude,energy,1,w=1/erramplitude)
p=sp.poly1d(fit)
T=sp.linspace(0,250,1000)
plt.plot(T,p(T))

print(fit)

expenergy=p(Amplitude)


#average perventage error
delta=expenergy-energy
print(delta)

a=0
for i in range(len(delta)):
    a+=(delta[i]/energy[i])**2
    
print(a)
totalerror=a**0.5
print('total percentage error',totalerror)

plt.grid()

#%%>>>>>>>energy and peak width

'''
E=sp.array([149.81,238.80,192,21.72])
errE=sp.array([0.13,0.28,0.05,0.01])
W=sp.array([5.38,8.49,5.7,1.51])
errW=sp.array([0.42,1.09,0.23,0.03])

errEsqrt=0.5*E**(-0.5)*errE

plt.errorbar(W,E,fmt='o',capsize=3,yerr=errE,xerr=errW)

fit = sp.polyfit(W,E,1,w=errW)
p=sp.poly1d(fit)
T=sp.linspace(0,10,100)
#plt.plot(T,p(T))

plt.grid()

'''
# peak width is proportional to the square root of energy

E=sp.array([149.81,238.80,192,21.72])
errE=sp.array([0.13,0.28,0.05,0.01])
W=sp.array([7.93,10.49,9.58,2.2])
errW=sp.array([0.17,0.43,0.18,0.01])

errEsqrt=0.5*E**(-0.5)*errE

plt.errorbar(W,E**0.5,fmt='o',capsize=3,yerr=errEsqrt,xerr=errW)

fit = sp.polyfit(W,E**0.5,1,w=errW)
p=sp.poly1d(fit)
T=sp.linspace(2,12,100)
plt.plot(T,p(T))
plt.xlabel('Width W',size=30)
plt.ylabel('square root of energy $E^{0.5}$',size=30)

plt.grid()

#%% >>>>>fitting for the position of back scattering peaks

# fit for Cs peak
cCs0=cCs-cno
plt.errorbar(A,cCs0,fmt='o',yerr=cCs0**0.5,capsize=3)

Ap=[]
cCsp=[]
for i in range(0,len(A)):
    if 50<A[i]<100:
        Ap.append(A[i])
        cCsp.append(cCs0[i])
        
plt.errorbar(Ap,cCsp,fmt='o',yerr=sp.array(cCsp)**0.5,capsize=3)
        
        

fitcCs,covcCs = curve_fit(nor, Ap, cCsp,[55,5,300,0,0])
print(fitcCs)

for i in range(0,5):
    print(sp.sqrt(covcCs[i,i]))
    
Ax=sp.linspace(0,500,1000)

plt.plot(Ax,nor(Ax,*fitcCs))

# peak for Cs at 192.00±0.05

plt.xlabel('amplitude',size=30)
plt.ylabel('frequency',size=30)
plt.title('Cs peak',size=40)
plt.grid()

#peak at 192, main peak 0 degree
#peak at 122.89, back-scattering give max energy electron, which de-excite
#peak at 59.46 back-scattered photon
#peak at 27.01 characteristic x-ray of iodine in the scintillation material
#peak at 11.47 characteristic x-ray of lead on the equipment

#back scattering is more dominant than other angles---Klein-Nishina scattering cross-section



#%%#energy peak relationship 
E=sp.array([192,59.46,27.01,11.47])
errE=sp.array([0.05,0.64,0.17,0.10])
W=sp.array([9.45,6.66,3.49,2.14])
errW=sp.array([0.09,1.06,0.29,0.16])

errEsqrt=0.5*E**(-0.5)*errE

plt.errorbar(W,E**0.5,fmt='o',capsize=3,yerr=errE,xerr=errW)

fit = sp.polyfit(W,E**0.5,1,w=errW)
p=sp.poly1d(fit)
T=sp.linspace(2,12,100)
plt.plot(T,p(T))
plt.xlabel('Width W',size=30)
plt.ylabel('square root of energy $E^{0.5}$',size=30)

plt.grid()


#peak at 122.89 is excluded since it has different broadening mechanism. The photons emitted by electron de-excitation can be in any direction, introducing another peak broadeing step.


#%%
#plt.plot(A,cCs/10)

AA=p(A)
plt.errorbar(AA,sCs1-bCs1,fmt='o',yerr=np.abs(sCs1-bCs1)**0.5,capsize=3)

A1=[]
sCsp1=[]
for i in range(0,len(A)):
    if 400<AA[i]<900:
        A1.append(AA[i])
        sCsp1.append(sCs1[i]-bCs1[i])
plt.errorbar(A1,sCsp1,fmt='o',yerr=np.abs(sp.array(sCsp1))**0.5,capsize=3)


'''
#group the data
gA1=[]
gsCsp1=[]
for i in range(0,len(A1)-8):
    gA1.append(sp.mean([A1[i],A1[i+1],A1[i+2],A1[i+3],A1[i+4],A1[i+5],A1[i+6],A1[i+7]]))
    gsCsp1.append(sp.mean([sCsp1[i],sCsp1[i+1],sCsp1[i+2],sCsp1[i+3],sCsp1[i+4],sCsp1[i+5],sCsp1[i+6],sCsp1[i+7]]))
plt.errorbar(gA1,gsCsp1,fmt='o',yerr=sp.array(gsCsp1)**0.5)


    '''    
    





fitsCs,covsCs = curve_fit(nor, A1, sCsp1,[ 5.95159322e+02,  3.73930801e+01,  9.32464930e+00, -3.72857332e-01,
  5.19971127e-04])
print(fitsCs)

for i in range(0,4):
    print(covsCs[i,i])
    
Ax=sp.linspace(0,1750,1000)

plt.plot(Ax,nor(Ax,*fitsCs))




fitsCs,covsCs = curve_fit(voi, A1, sCsp1,[ 4.19738061e+01, -4.14325173e+01,  3.31997291e+02,  5.94834892e+02,
  8.53485283e-01, -5.99601714e-04])
print(fitsCs)

for i in range(0,4):
    print(covsCs[i,i])
    
Ax=sp.linspace(0,1750,1000)

plt.plot(Ax,voi(Ax,*fitsCs))


fitsCs,covsCs = curve_fit(voi, A1, sCsp1,[ 4.19738061e+01, -4.14325173e+01,  3.31997291e+02,  5.94834892e+02,
  8.53485283e-01, -5.99601714e-04])
print(fitsCs)

for i in range(0,4):
    print(covsCs[i,i])
    
Ax=sp.linspace(0,1750,1000)

plt.plot(Ax,voi(Ax,*fitsCs))








#plt.plot(A,sCs2,'o')
#plt.plot(A,sCs3,'o')
plt.legend(['20','30','45'])
plt.grid()

print(0.49471538412305943**0.5+6.10)
print('the peak for 20 degree is at 610.01±6.80')





#%%

plt.errorbar(AA,sCs2-bCs2,fmt='o',yerr=(np.abs(sCs2-bCs2)**0.5),capsize=3)

A2=[]
sCsp2=[]
for i in range(0,len(AA)):
    if 400<AA[i]<750:
        A2.append(AA[i])
        sCsp2.append(sCs2[i]-bCs2[i])
plt.errorbar(A2,sCsp2,fmt='o',yerr=sp.sqrt(sp.array(np.abs(sCsp2))),capsize=3)



'''
gA2=[]
gsCsp2=[]
for i in range(0,len(A2)-8):
    gA2.append(sp.mean([A2[i],A2[i+1],A2[i+2],A2[i+3],A2[i+4],A2[i+5],A2[i+6],A2[i+7]]))
    gsCsp2.append(sp.mean([sCsp2[i],sCsp2[i+1],sCsp2[i+2],sCsp2[i+3],sCsp2[i+4],sCsp2[i+5],sCsp2[i+6],sCsp2[i+7]]))
plt.plot(gA2,gsCsp2,'o')
'''


        

        
fitsCs2,covsCs2 = curve_fit(nor, A2, sCsp2,[570,5,20,0,-0.05])
print(fitsCs2)

for i in range(0,4):
    print(covsCs2[i,i])
    
Ax=sp.linspace(0,1750,1000)

plt.plot(Ax,nor(Ax,*fitsCs2))





fitsCs2,covsCs2 = curve_fit(voi, A2, sCsp2,[10,10,20,550,0,0])
print(fitsCs2)

for i in range(0,4):
    print(covsCs2[i,i])
    
Ax=sp.linspace(0,1750,1000)

plt.plot(Ax,voi(Ax,*fitsCs2))


        
        
plt.plot(A2,sCsp2,'o')
#plt.plot(A,sCs2,'o')
#plt.plot(A,sCs3,'o')
#plt.legend(['20','30','45'])
plt.grid()




'''
plt.plot(sp.linspace(0,995,200),ampd)
plt.plot(gA2,gsCsp2,'o')

'''




        
print(0.5842594969439249**0.5+5.49)
print('the peak for 30 degree is at 549.13±6.25')


#%%
plt.errorbar(AA,sCs3-bCs3,fmt='o',yerr=np.abs(sCs3-bCs3)**0.5,capsize=3)






A3=[]
sCsp3=[]
for i in range(0,len(AA)):
    if 400<AA[i]<750:
        A3.append(AA[i])
        sCsp3.append(sCs3[i]-bCs3[i])
plt.errorbar(A3,sCsp3,fmt='o',yerr=np.abs(sCsp3)**0.5,capsize=3)

'''
        
gA3=[]
gsCsp3=[]
for i in range(0,len(A3)-8):
    gA3.append(sp.mean([A3[i],A3[i+1],A3[i+2],A3[i+3],A3[i+4],A3[i+5],A3[i+6],A3[i+7]]))
    gsCsp3.append(sp.mean([sCsp3[i],sCsp3[i+1],sCsp3[i+2],sCsp3[i+3],sCsp3[i+4],sCsp3[i+5],sCsp3[i+6],sCsp3[i+7]]))
plt.plot(gA3,gsCsp3,'o')
'''        
    
        
plt.plot(A3,sCsp3,'o')
        
fitsCs3,covsCs3 = curve_fit(nor, A3, sCsp3,[500,3,30,0,-0.05])
print(fitsCs3)

for i in range(0,4):
    print(covsCs3[i,i])
    
Ax=sp.linspace(0,1750,1000)

plt.plot(Ax,nor(Ax,*fitsCs3))
        
    
    
fitsCs2,covsCs2 = curve_fit(voi, A3, sCsp3,[50,50,300,470,0,0])
print(fitsCs2)

for i in range(0,4):
    print(covsCs2[i,i])
    
Ax=sp.linspace(0,1750,1000)

plt.plot(Ax,voi(Ax,*fitsCs2))


        

#plt.plot(A,sCs2,'o')
#plt.plot(A,sCs3,'o')
#plt.legend(['20','30','45'])
plt.grid()
print(0.5035778448149948**0.5+4.70)
print('the peak for 30 degree is at 470.21±5.41')



#%%


Eth=sp.array([662,609.41,549.43,470.41])*1.6e-16
errEth=sp.array([0.1,2.09,3.18,1.98])*1.6e-16
ang=sp.array([0,20,30,45])

plt.errorbar(sp.cos(ang*sp.pi/180),1/Eth,yerr=errEth/Eth**2,fmt='o',capsize=3)

fit= sp.polyfit(sp.cos(ang*sp.pi/180),1/Eth,1,w=errEth/Eth**2)
p=sp.poly1d(fit)
T=sp.linspace(0.67,1.05,1000)
plt.plot(T,p(T))
print(fit)

#plt.plot(sp.cos(ang*sp.pi/180),p(sp.cos(ang*sp.pi/180)))

dif=1/Eth-p(sp.cos(ang*sp.pi/180))
print(dif)

var=0
for i in range(len(dif)):
    var+=dif[i]**2
    
print('variance',var,'sd',var**0.5)

plt.grid()

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


n=55 #the number of points

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
plt.errorbar(sp.linspace(0,520,105)/0.512,ampd,'o')
plt.plot(sp.array(gA1),gsCsp1,'o')
#for a<Obf<=b:






fitcCs,covcCs = curve_fit(nor, sp.array(gA1),gsCsp1,[615,10,50,10,0])
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,1750,1000)
plt.plot(Ax,nor(Ax,*fitcCs))


xx1=sp.linspace(0,520,105)/0.512

    
fitcCs,covcCs = curve_fit(voi, sp.array(gA1),gsCsp1,[10,10,30,600,0,0])
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,1750,1000)
plt.plot(Ax,voi(Ax,*fitcCs))


xx1=sp.linspace(0,520,105)/0.512


fitcCs1,covcCs1 = curve_fit(nnor, xx1,ampd,[595,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)

for i in range(0,4):
    print(covcCs1[i,i])
    
plt.plot(Ax,nnor(Ax,*fitcCs1))
plt.title('Frequency-energy graph at 20 degree')



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
plt.plot(sp.linspace(0,520,105)/0.512,ampd,'o')
plt.plot(sp.linspace(0,520,105)/0.512,ampd,'o')
plt.plot(sp.array(gA2),gsCsp2,'o')
#for a<Obf<=b:


fitcCs,covcCs = curve_fit(nnor, sp.array(gA2),gsCsp2,[600,10,50,10])
print('experiment, U,S,A,M',fitcCs)

for i in range(0,4):
    print(covcCs[i,i])
    
Ax=sp.linspace(0,1750,1000)
plt.plot(Ax,nnor(Ax,*fitcCs))


fitcCs1,covcCs1 = curve_fit(nnor, sp.linspace(0,520,105)/0.512,ampd,[595,5,1000,100])
print('simulation, U,S,A,M',fitcCs1)

for i in range(0,4):
    print(covcCs1[i,i])
    
plt.plot(Ax,nnor(Ax,*fitcCs1))
plt.title('Frequency-energy graph at 30 degree')



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
plt.plot(sp.linspace(0,520,105)/0.512,ampd,'o')
plt.plot(sp.linspace(0,520,105)/0.512,ampd,'o')
plt.plot(sp.array(gA3),gsCsp3,'o')
#for a<Obf<=b:


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



#%%
#probability distribution over the angle
def sec(th,Eo):# define Klein-Nishina cross section
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
plt.ylabel('Cross section')
plt.grid()

plt.show()

plt.plot(E(x,662),y) 
plt.xlabel('Scattered photon energy')
plt.ylabel('Cross section')
plt.grid()

plt.show()

plt.plot(612-E(x,662),y)
plt.xlabel('Scattered electron energy')
plt.ylabel('Cross section')
plt.grid()


#%%
def sec(th,Eo):
    a=1/(1+Eo/511*(1-th))
    return 0.5*(2.8179e-15)**2*(a)**2*((a+1/a)-2*(1-th**2))/1.3350666520239514e-29 #normalised

def E(t,EE):
    return EE/(1+EE/511*(1-t))        

  
norm=spi.quad(sec,-1,1,args=662)
x=sp.arange(-1,1,0.01)


En=sp.arange(min(E(x,662)),max(E(x,662)),0.1)
theta=1-(662/En-1)*511/662
y=sec(theta,662)

#plt.plot(En,y)



#create a distribution from 0 to 1000

dispa=sp.arange(0,1000,0.1)
newdist=[0]*len(dispa)
for i in range(int(min(E(x,662))/0.1),int(max(E(x,662))/0.1)):
    newdist[i]+=y[i-int(min(E(x,662))/0.1)]
       

plt.plot(dispa,newdist) #scattered photon

#scattered electron distribution

#plt.plot(662-dispa,newdist)

disnew=[]
for i in range(len(newdist)):
    disnew.append(newdist[len(newdist)-i-1])
#plt.plot(dispa,disnew)

diste=[]
for i in range(len(disnew)):
    if disnew[i]!=0:
        diste.append(disnew[i])
for i in range(len(dispa)-len(diste)):
    diste.append(0)

print(len(diste),len(dispa))
plt.plot(dispa,sp.array(diste))

plt.plot(dispa,nor(dispa,11.5*3.53-17,(11.5*3.53-17)**1.5/343,1,0,0))
plt.plot(dispa,nor(dispa,28.5*3.53-17,(28.5*3.53-17)**1.5/343,1,0,0))
plt.plot(dispa,nor(dispa,662,662*0.075,1.5,0,0))

backscadel=sp.array(newdist)+sp.array(diste)*sp.exp(0.001*dispa-1)+nor(dispa,662,662**1.5/343,0.0,0,0)
plt.plot(dispa,backscadel*2000)



#normpeak1=nor(dispa,662,res*Ai/2/0.512,2,0,0)


cov=np.convolve(nor(dispa,662,662*0.075,1,0,0),sp.array(newdist)+sp.array(diste)*sp.exp(0.001*dispa-1)+nor(dispa,662,662*0.075,1,0,0))

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
plt.errorbar(A*3.53-17,cCs-cno,fmt='o',yerr=sp.sqrt(np.abs(cCs-cno)))


#%%

plt.plot(dispa,sp.array(newdist))
plt.plot(dispa,sp.array(diste))


for i in range(len(diste[300:])):
    if diste[300:][i]==0:
        print(i+300)
        break
        
cutting=dispa[777] 
print(cutting)

newdistex=[]

for i in range(len(newdist)):
    if dispa[i]<=477:
        newdistex.append(newdist[i])
    else:
        newdistex.append(0)
        
plt.plot(dispa,newdistex)

count=0

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
    if dispa[i]<=477:
        shift.append(0.5)
    elif 477<dispa[i]<=662:
        shift.append(0.25)
    else:
        shift.append(0)

plt.plot(dispa,disteex)
sim=sp.array(newdistex)+sp.array(disteex)+nor(dispa,11.5*3.53-17,(11.5*5)**1.5/343,1,0,0)+nor(dispa,27.01*3.53-17,(27.01*3.53-17)**1.5/343,0.5,0,0)+sp.array(shift)
plt.plot(dispa,sim)

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
ampd1=sp.array([0]*num)

for i in range(0,n):
    for j in range(0,num):
        if 5*j<=Obf[i]<5*j+5:
            ampd1[j]+=1
        
        
#print(ampd)
#print( sp.linspace(0,1000,201))   

plt.errorbar(sp.linspace(0,520,105)/0.512,ampd1,fmt='o',yerr=sp.sqrt(ampd),capsize=3)
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
#plt.plot(sp.linspace(0,520,105)/0.512,ampd,'o')
#plt.plot(sp.array(gA3),gsCsp3,'o')
#for a<Obf<=b:

'''
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

'''
plt.errorbar(A*3.53-17,cCs-cno,fmt='o',yerr=sp.sqrt(np.abs(cCs-cno)))
#plt.errorbar(sp.linspace(0,520,105)/0.512,ampd1,fmt='o',yerr=sp.sqrt(ampd),capsize=3)


