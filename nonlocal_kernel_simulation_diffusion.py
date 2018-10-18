
# coding: utf-8

# In[8]:


import numpy as np


# In[38]:


class simulator_1D():
    
    def __init__(self,E,rho,loadfunc,Ttotal,dt,Nnodes):
        self.E = E
        self.rho = rho
        self.loadfunc = loadfunc
        self.Ttotal = Ttotal
        self.dt = dt
        self.Nnodes = Nnodes
        
        
    def nonlocal_kernel_sim(self,kernel):
    
        # This is the function that use the given nonlocal kernel and the given material properties and the displacement
        # BC to do a simulation on 1-D elastic dynamic equation
        
        E = self.E
        rho = self.rho
        loadfunc = self.loadfunc
        Ttotal = self.Ttotal
        dt = self.dt
        Nnodes = self.Nnodes
    
        horizon = np.size(kernel)
        self.horizon = horizon
        kernel0 = 2*np.sum(kernel) # for symmetric kernel only
        Tsteps = int(Ttotal/dt)
        N = Nnodes+2*horizon # add some nodes at the boundaries for nonlocal BCs
        self.nonlocal_N = N
        u = np.zeros((N,Tsteps))
        t = np.linspace(0,Ttotal,Tsteps)
        self.t = t
    
        # some arrays for storing velocity and acceleration  
        ut1 = np.zeros(N)
        ut2 = np.zeros(N)


    
        for tt in range(Tsteps-1):
            for i in range(horizon,N-horizon):
                temp = -kernel0*u[i,tt] # temporary variable
            
                for k in range(horizon):
                    temp += kernel[k]*(u[i+k+1,tt]+u[i-k-1,tt])
            
                # use the explicit finite difference for time
                ut2[i] = E/rho*temp
                u[i,tt+1] = u[i,tt] + dt*(ut1[i]+ut2[i])/2
                ut1[i] = ut2[i]
        
            # BC at the fixed end
            for i in range(horizon):
                u[i,tt+1]= -u[2*horizon-1-i,tt+1]
        
            # BC with displacement load
            for i in range(N-horizon,N):
                P = loadfunc(t[tt])
                u[i,tt+1] = 2*P-u[2*(N-horizon)-1-i,tt+1]
    
        return u
    
    def nonlocal_kernel_middisplacement(self,kernel):
        
        # This function returns only the displacement at the midpoint of the rod
        
        u = self.nonlocal_kernel_sim(kernel)
        umid = (u[int(self.nonlocal_N/2),:]+u[int(self.nonlocal_N/2-1),:])/2
        
        return umid
        
    

