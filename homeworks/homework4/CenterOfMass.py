
# coding: utf-8

# In[ ]:


# Homework 4
# Center of Mass Position and Velocity
# Charlotte Kevis


# ### Keep in mind this is just a template; you don't need to follow every step and feel free to change anything.
# ### We also strongly encourage you to try to develop your own method to solve the homework.

# In[2]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[16]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        #the positions and velocities are defined by the index as x,y,z and vx,vy,vz
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)
    
    #define strings for each component
        xstr = []
        ystr = []
        zstr = []
        msum = np.sum(m)
    
            #the component for each vector is given as the sum of x*m / the sum of m
        #so we have to calculate the sums for all of these components
        #I'm not sure how to do this other than in a for loop
        for i in range (len(a)):
            #for the x component
            xstr.append(np.sum(a[i]*m[i]))
            #this will change the value of xstr based on the values imported from the data
            #now we repeat this for y and z
            ystr.append(np.sum(b[i]*m[i]))
            zstr.append(np.sum(c[i]*m[i]))
        xcom = np.sum(xstr)
        #this will sum all of the calculations made by the previous line
        ycom = np.sum(ystr)
        zcom = np.sum(zstr)
        #now we have the numerator value for the components from the equation listed before    
        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Acom = xcom/msum
        # ycomponent Center of mass
        Bcom = ycom/msum
        # zcomponent Center of mass
        Ccom = zcom/msum
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        RCOM = np.sqrt(XCOM**2+YCOM**2+ZCOM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        xNew = (self.x - XCOM)
        yNew = (self.y - YCOM)
        zNew = (self.z - ZCOM)
        RNEW = np.sqrt(xNew**2+yNew**2+zNew**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW <= RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2,m2)
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = np.sqrt(XCOM2**2+YCOM2**2+ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
       # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this.                                                                                              
            #print ("maxR", maxR)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            xNew = (x2-XCOM2)
            yNew = (y2-YCOM2)
            zNew = (z2-ZCOM2)
            RNEW = np.sqrt(xNew**2+yNew**2+zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units usint astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        return np.round(COMP*u.kpc,2)
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x*u.kpc - COMX
        yV = self.y*u.kpc - COMY
        zV = self.z*u.kpc - COMZ
        RV = np.sqrt(xV**2+yV**2+zV**2)
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV <= RVMAX)

        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew =  self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew,mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # we are rounding all of the values to 2 decimal places
        VXCOMr = np.around(VXCOM,2)
        VYCOMr = np.around(VYCOM,2)
        VZCOMr = np.around(VZCOM,2)
        # write your own code below
        COMV = [VXCOMr, VYCOMr, VZCOMr]

        # return the COM vector                                                                                        
        return COMV*(u.km/u.s)


# In[17]:


# Create a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MWCOM = CenterOfMass("MW_000.txt", 2)
#then we make the center of mass object in the same way for the other two galaxies
M31COM = CenterOfMass("M31_000.txt",2)
M33COM = CenterOfMass("M33_000.txt",2)


# In[18]:


# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM
MW_COMP = MWCOM.COM_P(0.1)
MW_COMV = MWCOM.COM_V(MW_COMP[0], MW_COMP[1], MW_COMP[2])
#now we do the same for M33 and M31
M33_COMP = M33COM.COM_P(0.1)
M33_COMV = M33COM.COM_V(M33_COMP[0], M33_COMP[1], M33_COMP[2])
M31_COMP = M31COM.COM_P(0.1)
M31_COMV = M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2])

#now we have to print all of our calculated values
print("MIlky Way")
print(MW_COMP)
print(MW_COMV)
print("M33")
print(M33_COMP)
print(M33_COMV)
print("M31")
print(M31_COMP)
print(M31_COMV)


# In[ ]:


# now write your own code to answer questions


# In[20]:


#determining the current separation and velocity between the MW and M31
#we can use the distance formula sqrt( (x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
#first we define the strings to calculate the answers
MWM31dif = []
MWM31vel = []

#now we define the final calues
MWM31distance = 0.
MWM31velocity = 0.

#now we compute using the distance formula
MWM31dif = np.subtract(MW_COMP,M31_COMP)**2
MWM31distance = np.sqrt(np.sum(MWM31dif))

#for the velocity we do the dame thing
MWM31vel = np.subtract(MW_COMV,M31_COMV)**2
MWM31velocity = np.sqrt(np.sum(MWM31vel))

#now we print thes numbers rounded to 2 decimal places
print(np.round(MWM31distance,2))
#we know the distance is ~770 kpc, which is what we get
print(np.round(MWM31velocity,2))
#we know the velocity is ~110 km/s. Our calculation is a little high but still accurate


# In[21]:


#now we are doing the same thing for the separation and velocity between M33 and M31
M3331dif = []
M3331vel = []

M3331distance = 0.
M3331velocity = 0.

M3331dif = np.subtract(M33_COMP,M31_COMP)**2
M3331distance = np.sqrt(np.sum(M3331dif))

print(np.round(M3331distance,2))
#i am not sure what the distance is

M3331vel = np.subtract(M33_COMV, M31_COMV)**2
M3331velocity = np.sqrt(np.sum(M3331vel))

print(np.round(M3331velocity,2))
#the answer from class is ~202 km/s. Here we are a little low but still very close


# Each component is constantly changing as the two galaxies approach each other. As they begin to merge, the components will shift and add into each other. Understanding how the galaxies components are moving will give us a good estimate/model of what the two galaxies are going to do in the future. This information will tell us when the two galaxies will begin to merge and at what distance. We can learn a lot about galaxy mergers from understanding how our two galaxies are behaving and how they may behave in the future.
