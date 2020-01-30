
# coding: utf-8

# The first program for homework 2, the read file function to collect data and sort.

# In[4]:


import numpy as np
import astropy.units as u


# This imports the packages that we need in order to complete the functions necessary to import data from a data file as well as assign each value units. Now we have to define the read function and how it reads the data.

# In[7]:


def Read(filename):
    
    #we have to define our variables for this section
    data=[]
    time = 0.
    totpart = 0.
    
    #now we open the file
    file = open(filename, 'r')
    #these lines import the data and assign units
    #the first column is time
    line1 = file.readline()
    label, tvalue = line1.split()
    time = float(tvalue)*10.0*u.Myr
    
    #the second column is the total number of particles
    line2 = file.readline()
    label, totpart = line2.split()
    #the number of particles does not have units
    
    file.close()
    
    #now we have to generate the data in order to print it
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    #now we have to return the data in order to have it print out 
    return time,totpart,data
    


# In[8]:


#these are the variables for printing to a file or code

#time
t=0
#total number of particles
p=0
#data
i=[]
#now we are pulling the values out of the code above
t,p,i=Read('MW_000.txt')

#printing the values
print(t)
print(p)
print(i[1])
#checking the particle type. This will print the type as defined by the data
print(i['type'][1])
