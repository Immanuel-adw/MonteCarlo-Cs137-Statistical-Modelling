# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 07:13:01 2021

@author: imman
"""

#%%
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
    
#if the peak cannot be fitted with voigt profile.
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


'''
n = 8
params = {   
    'axes.labelsize': 15,
   'font.size': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'figure.figsize': [n,n]
   } 
plt.rcParams.update(params)
'''


params = {   
    'axes.labelsize': 20,
   'font.size': 20,
   'legend.fontsize': 20,
   'xtick.labelsize': 20,
   'ytick.labelsize': 20,
   'figure.figsize': [10,8]
   } 
plt.rcParams.update(params)


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




plt.plot(dispa,sp.array(diste), 'g-')  # Scattered Electron Cross-section
plt.plot(dispa,newdist, 'b-') # Scattered Photon Cross-section
plt.plot(dispa, sp.array(newdist) + sp.array(diste), 'r-')
plt.xlabel('Incident photon Energy, MeV')
plt.ylabel('Klein-Nishina Cross section')
plt.legend(['Electron Incident Energy vs Cross-Section', 'Photon Incident Energy vs Cross-Section', 'Energies vs Sum of Cross-Sections'])
plt.title('Klein-Nishina Incident Photon E v. Cross-section')
plt.savefig('Immanuel Klein-Nishina')
plt.grid()
plt.show()

'''
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
'''
#normpeak1=nor(dispa,662,res*Ai/2/0.512,2,0,0)


# Convolution Not used ! ! !
#cov=np.convolve(nor(dispa,662,662*0.075,1,0,0),sp.array(newdist)+sp.array(diste)*sp.exp(0.001*dispa-1)+nor(dispa,662,662*0.075,1,0,0))

#plt.plot(sp.linspace(0,1000,len(cov)),cov)

#cov1=np.convolve(nor(dispa,11.5,11.5**1.5/343,2,0,0),cov)
#plt.plot(sp.linspace(0,1000,len(cov1)),cov1)
#plt.plot(dispa,sp.array(newdist)*sp.array(newdiste)*10)



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

shift = []
n = 0.00298
o = 110
for i in range(len(dispa)):
    if 110<dispa[i]<=181:       # 130<dispa[i]<=181:
        shift.append(0.302-((n*(181-dispa[i]))**2))
    else:
        shift.append(0)
        
shift3=[]
n = 0.0008
for i in range(len(dispa)):
    if 500 < dispa[i]<=606:       # 130<dispa[i]<=181:
        shift3.append(0.148-((n*(499-dispa[i]))**2))
    else:
        shift3.append(0)


Pb = nor(dispa,77.66,3*77.67**0.5,0.435,0,0)
Pb0 = []
for i in range(len(dispa)):
    if dispa[i]<110:        #(o)
        Pb0.append(Pb[i])
    else:
        Pb0.append(0)
Pb0=np.array(Pb0)


klein_nishina = 1.2*(sp.array(newdistex)+sp.array(disteex))
#Pb0 = nor(dispa,77.66,3*77.67**0.5,0.435,0,0)
Pb1 = nor(dispa,75,3*75**0.5,0.435,0,0)
Pb2 = nor(dispa,73,75**0.5,0.435,0,0)
Pb3 = nor(dispa,85,75**0.5,0.435,0,0)
Photopeak1 = nor(dispa,662,1.26*(662**0.5),1.4,0,0)
Photopeak = []

for i in range(len(dispa)):
    if 607<=dispa[i]<=800:       # 130<dispa[i]<=181:
        Photopeak.append(Photopeak1[i])
    else:
        Photopeak.append(0)

# Why is ar multiplicative-height 0.5 ?
# What do newdist and diste represent
# Why are there 2 Steps in the Step-Function?
# Why do we cuttoff Newdist and diste at the points where we cut them off?
# Is Backscadel ever used again?
#
#sim = klein_nishina + nor(dispa,11.5*3.53-17,(11.5*5)**1.5/343,1,0,0)+nor(dispa,27.01*3.53-17,(27.01*3.53-17)**1.5/343,0.5,0,0)+sp.array(shift)
#sim3 = (sp.array(newdistex))+(sp.array(disteex))+nor(dispa,11.5*3.53-17,(11.5*5)**1.5/343,0.68,0,0) + nor(dispa,662,662*0.075,1.4,0,0)
sim2_unscaled = sp.array(shift) + sp.array(shift3) + klein_nishina + nor(dispa,11.5*3.53-17,(11.5*3.53-17)**0.5,0.68,0,0) + Photopeak + Pb0
sim2 = 1000*sim2_unscaled

#plt.plot(dispa,1000*np.array(newdistex))
#plt.plot(dispa,1000*np.array(disteex))
plt.plot(dispa,sim2)
plt.xlabel('Energy, KeV')
plt.ylabel('Frequency')
#plt.legend(['Photon E-Cross-Section','Electron E-Cross-Section','Cs-137 - Compton Analysis Compton Scatter Model']) # at 0 degree measurement
plt.legend(['Cs-137 - Gamma-Spectrum Model'])
plt.title('Cs-137 - Gamma-Ray Spectrum Model')
plt.savefig('Montecarlo Theory Alone Plot')
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
plt.savefig('Montecarlo Theory Plot')
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

r = 0.1 #0.1 # 90/14 #21/90 # 14/90 * 1.5

#for l in range(300,len(dispa)):
print('len dispa ', len(dispa))
#the number of points
for l in range(90,589):     # Check here for 10 factor
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
    #s = 2*gain
    s = 0.9*((Ei[l])**0.5)*gain
    #s=  Ei[l]**0.5 * gain
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

#Ep1 = sp.array([int((11.5*3.53)-17), 77.7, 662]) 
Ep1 = sp.array([int((11.5*3.53)-17), 74, 85, 662])    
#Ep  = 10*E(Ep1*0.1)      
Ai=Ep1*gain
           
ampd=sp.array([0]*num)
#the number of points
    #Cgood = (cCs-cno)
for l in range(4):
    if l == 0:
        r = np.sqrt((0.68/1.4)*90/14) #1 #90/14
        n=int(sim2[int(Ep1[l])*10]*r) #
        #print('sim Ep1 l is ', sim2[int(Ep1[l]*9)])
        #print('n is ', n)
        # random seeds
        randseed=[]
        for i in range(0,n):
            randseed.append(random.uniform(0, 1))
        # noise
        #s=res*Ai[l]/2
        #s=(1.26683723*(Ep1[l])**0.5+1.92435425)*gain  # What are these?
        s=0.55*((Ep1[l])**0.5)*gain
        s=((Ep1[l])**0.5)*gain
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
    elif l == 1:
        r = np.sqrt((0.43/1.4)*90/14) #0.5*0.43*90/14 #2
        n=int(sim2[int(Ep1[l]*10)]*r) #
        #print('sim Ep1 l is ', sim2[int(Ep1[l]*9)])
        print('n is ', n)
        # random seeds
        randseed=[]
        for i in range(0,n):
            randseed.append(random.uniform(0, 1))
        # noise
        #s=res*Ai[l]/2
        s=2*(1.26683723*(Ep1[l])**0.5+1.92435425)*gain  # What are these?
        s=((Ep1[l])**0.5)*gain
        #s=14 #13.9  #Calculated as diff+max+minpeak root75+root85+85-75                       #3*((Ep1[l])**0.5)*gain    
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
    elif l == 2:
        r = np.sqrt((0.3/1.4)*90/14) #0.5*(0.3)*90/14
        n=int(sim2[int(Ep1[l]*10)]*r) #
        #print('sim Ep1 l is ', sim2[int(Ep1[l]*9)])
        print('n is ', n)
        # random seeds
        randseed=[]
        for i in range(0,n):
            randseed.append(random.uniform(0, 1))
        # noise
        #s=res*Ai[l]/2
        #s=3*(1.26683723*(Ep1[l])**0.5+1.92435425)*gain  # What are these?
        s=((Ep1[l])**0.5)*gain
        #s=14 #13.9  #Calculated as diff+max+minpeak root75+root85+85-75                       #3*((Ep1[l])**0.5)*gain    
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
    elif l == 3:
        r = 90/14
        n=int(sim2[int(Ep1[l]*10)]*r) #
        #print('sim Ep1 l is ', sim2[int(Ep1[l]*9)])
        print('n is ', n)
        # random seeds
        randseed=[]
        for i in range(0,n):
            randseed.append(random.uniform(0, 1))
        # noise
        #s=res*Ai[l]/2
        #s=(1.26683723*(Ep1[l])**0.5+1.92435425)*gain  # What are these?
        s=((Ep1[l])**0.5)*gain
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
#ampd = ampd

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


nK = 270
nq = 150
plt.plot(A[0:nK]*3.53-17,(cCs-cno)[0:nK],'yo')
plt.errorbar(A[0:nK]*3.53-17,(cCs-cno)[0:nK],fmt='y ',yerr=sp.sqrt(np.abs((cCs-cno)[0:nK])))
plt.plot((sp.linspace(0,520,105)/0.512)[0:nq],((ampd)+ampd1)[0:nq],'bo')
#plt.xticks(np.arange(0,800,5))
plt.grid()
plt.xlabel('Energy, KeV')
plt.ylabel('Frequency')
plt.legend(['Measured Data','Simulated Model Cs 137']) # at 0 degree measurement
plt.title('Numerical Monte-Carlo Simulation')
plt.savefig('Immanuels Numerical Monte-Carlo Test')
plt.show()

