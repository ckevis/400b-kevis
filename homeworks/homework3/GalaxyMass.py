
# coding: utf-8

# In[2]:


#a program to determine the total mass of a galaxy
#we have to import the main functions to complete the tasks given
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from ReadFile import Read


# In[10]:


#we have to define a function to determine the total mass of the galaxy by first defining the component mass
#in this function we will input the mass information from a data set and return the total mass for each particle type
def ComponentMass(filename, partype):
    #values from the read file:
        #time - the time in seconds
        #total - the total number of particles
        #dat - the data pulled from the file
    time,total,dat = Read(filename)
    
    #now we can pull the information from the file using an index similar to ParticleProperties.py
    index = np.where(dat['type']==partype)
    
    #Now we need to get all of the masses for each particle type in units of 10*12 Msun
    #there are 3 types: Halo, Disk, and Bulge
    #the units are currently in terms of 10*10 Msun so we must divide it by 10*2 Msun
    
    # [number][10^10Msun] / [10^2] = [number][10^12Msun]
    compm = dat['m'][index]*10**-2*u.Msun
    
    #now we have to sum the mass of each component
    totmas = np.sum(compm)
    return totmas


# In[12]:


#now we are going to define the mass for M31, M33, and the Milky Way

#partype 1 = Halo particles
#partype 2 = Disk particles
#partype 3 = Bulge Particles

#for the data located in MW_000.txt we define the mass of each particle type as follows
MWHalo = np.around(ComponentMass('MW_000.txt',partype=1),3)
MWDisk = np.around(ComponentMass('MW_000.txt',partype=2),3)
MWBulge = np.around(ComponentMass('MW_000.txt',partype=3),3)

#for the data in M33_000.txt we define the mass of each particle type as follows
M33Halo = np.around(ComponentMass('M33_000.txt',partype=1),3)
M33Disk = np.around(ComponentMass('M33_000.txt',partype=2),3)
M33Bulge = np.around(ComponentMass('M33_000.txt',partype=3),3)

#for the data in M31_000.txt we define the mass of each particle type as follows
M31Halo = np.around(ComponentMass('M31_000.txt',partype=1),3)
M31Disk = np.around(ComponentMass('M31_000.txt',partype=2),3)
M31Bulge = np.around(ComponentMass('M31_000.txt',partype=3),3)


# In[13]:


#Now we need to calculate the total mass of the galaxies
MWmass = MWHalo + MWDisk + MWBulge
#the masses are still in units of 10^12 Msun

#we are going to print both the total component masses and the total galaxy mass to make sure the numbers are added properly
print (MWHalo, MWDisk, MWBulge)
print (MWmass)

#now we are going to do the same for M33 and M31
M33mass = M33Halo + M33Disk + M33Bulge
print (M33Halo, M33Disk, M33Bulge)
print (M33mass)

M31mass = M31Halo + M31Disk + M31Bulge
print (M31Halo, M31Disk, M31Bulge)
print (M31mass)

#this should return the component mass of each galaxy and the total mass


# In[15]:


#now we are going to compute the barion fraction fbar, which is the fraction of the galactic mass that is stellar
#fbar = TotalStellarMass/TotalMass
    #the halo mass is the dark mass, while the bulge and disk masses are the stellar masses
fMWbar = np.around((MWDisk+MWBulge)/MWmass,3)
fM33bar = np.around((M33Disk+M33Bulge)/M33mass,3)
fM31bar = np.around((M31Disk+M31Bulge)/M31mass,3)
#since fbar is a ratio, it does not have units

#now we have to print those values
print(fMWbar,fM33bar,fM31bar)


# In[17]:


#now we do the same calculation for the entire local group, which is made up of M33, M31, the Milky Way, and their satellite galaxies
#to find fbar for the local group, we have to sum the components of each galaxy
localHalo = MWHalo + M33Halo + M31Halo
localDisk = MWDisk + M33Disk + M31Disk
localBulge = MWBulge + M33Bulge + M31Bulge

#now we calculate the total mass of the local group in 10^12 Msun
localtot = localHalo + localDisk + localBulge

#now we calculate the fbar relation for the entire local group
flocbar = np.around((localDisk + localBulge)/localtot,3)
print (localHalo, localDisk, localBulge, localtot,flocbar)

 
