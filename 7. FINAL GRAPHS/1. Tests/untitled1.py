# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 11:15:49 2021

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

n=10*9000/1400    #45 #the number of points

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


n=3.5*900/1400 #35 #the number of points

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


n=20*9000/1400 #90 #the number of points

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
