
# coding: utf-8

# In[50]:


#in this code we are going to define the properties of a particle from the ReadFile program
#importing the packs that we need to complete the following tasks
import numpy as np
import astropy.units as u
#here i am importing the function Read from part 1 of the homework
from ReadFile import Read


# In[90]:


#here we are going to define a function ParticleInfo
def ParticleInfo(filename,partype,parnumber):
    time, total, dat = Read(filename)
    #this will pull the data from the Read function
    
    #now we have to initialize the position and velocity function so they are 0 at the start and not
    #a random variable
    x=0
    y=0
    z=0
    vx=0
    vy=0
    vz=0
    
    #we have to define a way to pull the information from the Read file and import it into the code
    #to perform these tasks
    index = np.where(dat['type']==partype)
    
    #now we are defining the position and velocities using this index
    xnew = dat['x'][index]*u.kpc
    ynew = dat['y'][index]*u.kpc
    znew = dat['z'][index]*u.kpc
    vxnew = dat['vx'][index]*u.km/u.s
    vynew = dat['vy'][index]*u.km/u.s
    vznew = dat['vz'][index]*u.km/u.s
    
    #now we define the components of any given particle from the data set
    xcomp = xnew[parnumber]
    ycomp = ynew[parnumber]
    zcomp = znew[parnumber]
    vxcomp = vxnew[parnumber]
    vycomp = vynew[parnumber]
    vzcomp = vznew[parnumber]
    
    #using this information we can calculate the distance and velocity of the particle in 3D
    d = 0 #making sure the variable is initially set to 0 for proper calculations
    d = np.around((xcomp**2+ycomp**2+zcomp**2)**0.5,3) #the equation for distance is sqrt(x^2+y^2+z^2)
    
    #to convert from kpc to lyr
    print(np.around(d.to(u.lyr),3))
    
    #velocity is calculated the same
    v = 0
    v = np.around((vxcomp**2+vycomp**2+vzcomp**2)**0.5,3)
    
    #the mass of the particle is given, we just need to define it
    #the units are given in 10^10 solar masses
    m = dat['m'][index]*10**10*u.solMass
    mass = m[parnumber]
    
    return d,v,mass


# In[91]:


#now for the wanted outcome of the homework
#we are going to define the properties of the 100th disk particle
#we also have to convert the distance from kpc to ly
    
#here are the parameters
#they have different names than the functions listed above in the ParticleInfo function because they
#are a known variable
fname = "MW_000.txt"
ptype = 2 #the particle is a disk particle
partn = 99 #this will give us the 100th particle

#now we have to execute the code using the function given above
#this will add the variables we have listed above into the function above to produce a d,v, and mass
ParticleInfo(fname,ptype,partn)

