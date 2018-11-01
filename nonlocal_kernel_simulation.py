
# coding: utf-8



import numpy as np
import cffi
try:
    from explicit_time_api import ffi, lib
except ImportError:
    print("Error importing lib")
    pass




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
        loadfunc_array = np.zeros(Tsteps) # array for load at each time step
    
        # generate the load_array by the given load_function(fixed at one end, Neumann BC at the other end)   
        for tt in range(Tsteps):
            loadfunc_array[tt] = loadfunc(t[tt])
        
        # Use the CFFI to import the C code to do the simulation for all time steps by explicit time discretization 
        lib.explicit_time(ffi.cast("double *", ffi.from_buffer(u)), N, Tsteps, kernel0, ffi.cast("double *", ffi.from_buffer(kernel)), horizon, E ,rho, dt, ffi.cast("double *", ffi.from_buffer(loadfunc_array)))
         
    
        return u
    
    def nonlocal_kernel_middisplacement(self,kernel):
        
        # This function returns only the displacement at the midpoint of the rod
        
        u = self.nonlocal_kernel_sim(kernel)
        umid = (u[int(self.nonlocal_N/2),:]+u[int(self.nonlocal_N/2-1),:])/2
        
        return umid
        
    

