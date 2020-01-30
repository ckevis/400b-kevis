
# coding: utf-8

# # In Class Lab 1
# 
# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 

# In[3]:


import numpy as np


# In[4]:


#a function to compute the local standard of velocity
#4.75*mu*radius of the sun = VLSR + nu
#VLSR = 4.75*mu*Ro - vsun

def VLSR(Ro, mu=4.379, vsun=12.24):
    #Ro is the distance to the Sun from the Galactic Center (kpc)
    #mu is the proper motion of Sag A* (mas/yr)
    #vsun is the peculiar motion of the sun in the v direction (km/s)
    #VLSR, the local standard of rest (km/s)
    
    return 4.74*mu*Ro - vsun


# In[5]:


RoReid = 8.34   #distance to the Galactic center from Reid t al 2014 in kpc
RoGravity = 8.178 #distance to the galactic center from the GRAVITY collaboration Abuter 2019
RoSG = 7.9 #distance to the galactic center from the textbook by Sparke and Gallagher


# In[6]:


#compute VLSR using RoReid
VLSR(RoReid)


# In[8]:


#compute VLSR using RoGravity
VLSR(RoGravity)


# In[9]:


#compute VLSR using RoSG
VLSR(RoSG)


# ### b)
# 
# compute the orbital period of the sun using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr

# In[11]:


#orbital period of the Sun, using Ro from GRAVITY
#T = 2piRo/v [Gyr]
#v = vtan = VLSR + vsun
vtan = VLSR(RoGravity) + 12.24
T_grav = 2*np.pi*RoGravity/vtan
print(T_grav)


# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)

# In[12]:


#the number of rotations about the galactic center 
#age of the universe / orbital period
print(13.8/T_grav)


# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of $10^{10}$ M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4988e-6 kpc$^3$/Gyr$^2$/M$_\odot$
# 
# What about at 260 kpc (in units of 10$^{12}$ M$_\odot$) ? 

# In[19]:


#gravitational constant
G = 4.4988*10**-6 #[kpc^3Mo/Gyr^2]
#compute the mass enclosed within the solar radius assuming an isothermal sphere model
#density profile rho = VSLR*2/(4*pi*G*R^2)
#mass = integrate Rho dR
    #integrate rho 4*pi*r^2 dr
    #integrate VLSR*2/(4*pi*r^2)*(4*pi*r^2)dr
    #integrate VLSR*2/Gdr
def MassIso(r, VLSR=235):
    #VLSR the local standard of rest [km/s]
    #r is the distance from the galactic center [kpc]
    #returns mass enclosed [Mo]
    return VLSR**2/(G*r)


# In[20]:


#compute the mass enclosed within Ro
MIsoSolar = MassIso(RoGravity)
print(MIsoSolar)


# In[21]:


#compute the mass enclosed in 260 kpc
MIso = MassIso(260)
print(MIso)


# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of $10^{12}$ M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

# In[22]:


#potential for a Hernquist sphere
#phi = -GM/r^2
#using the Hernquist potential, the equation for the escape speed becomes vesc^2 = 2*|2GM/r+a|
#rearrange the escape speed equation vesc^2/2G(r+a) =M
#function that will determine the total halo mass needed to set a given escape speed at a given distance
#assuming a Hernquist profile
def Massesc(vesc,a,d):
    #input vesc (km/s) escape speed
    #input d (kpc) distance from the galactic center
    #imput a (kpc) Hernquist scale length
    #return Total Mass (Mo)
    return vesc**2/2/G*(a+d)


# In[23]:


##mass needed to keep Leo I bound, assuming a Hernquist Profile
MLeoI = Massesc(196,30,260)
print (MLeoI/1e12)


# In[24]:


comp = (MIso/1e12)/MLeoI
print (comp)

