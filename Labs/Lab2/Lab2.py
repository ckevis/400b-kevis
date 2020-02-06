
# coding: utf-8

# In[26]:


# ASTR 400 B 
# In Class Lab 2

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from scipy.integrate import quad # For integration
# Documentation and examples for quad : 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
# https://www.tutorialspoint.com/scipy/scipy_integrate.htm


# ## Part A:  Schechter Fxn
# 
# The galaxy luminosity function in the nearby universe is well described by a Schechter Function:
# 
# \begin{equation}
# \Phi(M)dM = ( 0.4 \, ln10 ) \, \phi_\ast \, 10^{0.4(M_\ast - M)(\alpha +1)} e^{-10^{0.4(M_\ast - M)}} dM
# \end{equation}
# 
# With the following parameters from Smith+2009 for Field Galaxies in SDSS at z$\sim$0.1 in the Kband:
# 
# 
#  $\phi_\ast$ =1.66 $  \times 10^{-2}$  $h^3$ Mpc$^{-3}$
# 
#  $\alpha$ =  -0.81 
# 
# 
#   M$_\ast$ =  M$_k^\ast$= -23.19  - 5*log($h$)
#   
#  $h$ = the Hubble constant in units of 100 km/s/Mpc . At z=0 this is 0.7. But we are going to ignore it here. Units will then be in "comoving" coordinates.
#   
#   This function is defined for you below:

# In[2]:


# Function for the Schechter Luminosity function (in terms of Magnitude)
def Schechter(M,phistar=0.0166,Mstar=-23.19,alpha=-0.81):
# Takes as input:  an array of Absolute Magnitudes M 
    # phistar, number density of galaxies (normalization)
           #0.0166*h^3 Mpc^-3  defaults from Smith+2009 in Kband
    # Mstar (Knee of the Schechter Fxn):
           #-23.19 - 5*log(h) defaults from Smith+2009 in Kband
    # alpha (Slope of the Fxn) 
           #-0.81 defaults from Smith+2009 in Kband
# Returns: number density of galaxies (comoving units) at that magnitude M.
    return 0.4*np.log(10)*phistar*10**(0.4*(Mstar - M)*(alpha +1.0))*np.exp(-10**(0.4*(Mstar - M)))


# # Q1 
# 
# Utilizing the defined function, plot the Schechter Function using the above parameter values over a magnitude range of -17 to -26. 
# Try to reproduce the black solid line in Smith+2009 MNRAS 397,868 [UKIDSS Survey] Figure below.
# 
# 
# ![Smith09.png](attachment:Smith09.png)

# # Q2 
# 
# Galaxies in the Virgo Cluster have different parameters, like $\alpha$=-1.35  (Ferrarese+2016 ApJ 824).
# 
# Overplot the Schechter Function with this new value of $\alpha$.  
# 
# Try a smaller value of $\alpha = -0.6$.
# 
# How does the function change?  What does this mean? 
# 

# In[3]:


# Create an array to store Kband Magnitudes from -26 to -17
MK = np.arange(-26,-17,0.1)


# In[11]:


# Plot the Schechter Function

fig = plt.figure(figsize=(10,10))  # sets the scale of the figure
ax = plt.subplot(111) 

# Plot the default values (y axis log)
# ADD HERE
ax.semilogy(MK,Schechter(MK),label='Smith09', color = 'blue', linestyle=':',linewidth=5)
ax.semilogy(MK,Schechter(MK,alpha=-0.6),label='alpha = -0.6', linestyle='--', color = 'orange',linewidth=5)
ax.semilogy(MK,Schechter(MK,alpha=-1.35),label='alpha = -1.35', color = 'green',linewidth=5)

# Q2 solutions: change alpha
# ADD HERE


# Add labels
plt.xlabel(r'M$_k$ + 5Log($h$)', fontsize=22)
plt.ylabel(r'$\Phi$ (Mpc$^{-3}h^3$/mag)', fontsize=22)

#set axis limits
plt.xlim(-17,-26)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

# Save to a file
#ax.set_rasterized(True)
#plt.savefig('Schechter.eps', rasterized=True, dpi=350)


# In[19]:


#schechner luminosity function in terms of luminosity
def SchechterL(L,nstar=8e-3,Lstar=1.4e10,alpha=-0.7):
    #L array of lumniosities (Lsun)
    #nstar the mean galaxy number density (Mpc*-3)
    #Lstar characteristic luminosity (Lsun)
    #alpha faint end slope
    #returns number density of galaxies (Moc*-3)
    
    A = (L/Lstar)**alpha #controls the faint end slope
    B = np.exp(-L/Lstar) #controls the high luminosity end
    C = nstar/Lstar #constants
    
    return A*B*C


# In[20]:


#understand lambda function
#short cut - defines and evaluates a function in one line

x = lambda a,b : a*b
print(x(5,6))


# In[27]:


#what fraction of the integrated luminosity density lies above Lstar?
Lupper = quad(lambda L: L*SchechterL(L),1.4e10,1e14)
print(Lupper)

Ltotal = quad(lambda L: L*SchechterL(L), 0.1,1e14)
print(Lupper[0]/Ltotal[0])


# ## Part B: IMF 
# 
# Create a function called {\it Salpeter} that defines the Salpeter IMF: 
# 
# \begin{equation}
# \xi(M) = \xi_0 (M/M_\odot)^{-\alpha}
# \end{equation}
# 
# $\alpha = 2.35$
# The function should take as input an array of stellar masses, M. 
# You will need to determine the normalization, $\xi_0$, by integrating this equation over mass from 0.1 to 120 M$_\odot$
# and setting the value to 1.  The function should then return $\xi(M)$, which will now represent the fractional number of stars. 
# 
# Integration: quad(lambda x:  fxn(x),xmin,xmax)
# 
# quad returns an array with 2 values. you want the first value. 
# Note I've used a "lambda" expression.   Python's lambda expressions allow a function to be created and passed around all in one line of code

# In[28]:


#create a function that evaluates the Salpeter Initial Mass Function
#from 0.1-120 Msun
def Salpeter (M,Mmin = 0.1,Mmax = 120):
    #M is the stellar mass (Msun)
    #Mmin is the minimum mass to be considered (Msun)
    #Mmax is the maximal mass to be considered (Msun)
    #returns the normalized fraction of stars between Mmin and Mmax
    
    alpha = 2.35 #power law for the Salpeter IMF
    #normalize the Salpeter IMF
    Norm = quad(lambda M: M**(-alpha), Mmin, Mmax)
    #normalizing to 1
    A = 1./Norm[0]
    
    return A*M**(-alpha)


# In[31]:


Test = quad (lambda M: Salpeter(M), 0.1, 120)
print(Test[0])


# ## Q1: 
# Integrate your normalized function to compute the fraction of stars with stellar masses greater than the sun and less 
# than 120 M$_\odot$.
# 
# Double Check: if you integrate your function from 0.1 to 120 you should return 1.0 
# 

# ## Q2:
# 
# How might you modify the above to return the fraction of MASS ? instead of numbers of stars.

# In[36]:


def SalpeterN(M,Mmin=0.1,Mmax = 120):
    #M is the mass (Msun)
    #Mmin is the minimum mass we consider (Msun)
    #Mmax is the maximum mass we consider (Msun)
    #returns the normalized Salpeter IMF *Mass
    
    alpha = 2.35
    #normalize the function
    Norm = quad(lambda M: M*M**(-alpha),Mmin,Mmax)
    A = 1./Norm[0]
    
    #returns Salpeter IMF*M normalized
    return A*M**M(-alpha)


# In[37]:


#determine the fraction of mass in stars that are more massive than the Sun
MassAboveSun = quad(lambda M: SalpeterN(M), 1, 120)
print(MassAboveSun[0])

