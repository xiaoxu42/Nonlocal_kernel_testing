import unittest
from nonlocal_kernel import Kernel_calculator as kc
from nonlocal_kernel_simulation import simulator_1D as sim1D
import numpy as np

class Testsimulator_1D(unittest.TestCase):
    #this is the test for nonlocal_kernel_simulation.py this test may take about more than half minute

    def test_nonlocalresult_with_FEMresult(self):
        
        #mechanics properties
        Eh = 200e9
        Es = 5e9
        rhoh = 8000
        rhos = 8000
        rhoave = (rhoh+rhos)/2
        # geometric properties
        L = 1
        l = 0.02
        A = 1e-4
        # simulation parameters
        Ttotal = 1e-3
        dt = 0.2e-9
        Nnodes = int(L/l)
        Tsteps = int(Ttotal/dt)

        nonlocal_kernel = kc(Eh,Es,rhoh,rhos)
        tolerance = 1e-4 
        fourth_order_kernel = nonlocal_kernel.kernel_generator(4,tolerance)

        displacement_load_array = np.zeros(Tsteps)


        #generate the BC condition array using the BC function defined later
        for tt in range(Tsteps-1):
            t = (tt+1)*dt
            displacement_load_array[tt+1] = displacement_load(t)

        nonlocal_kernel_result = sim1D(Eh,rhoave,displacement_load_array,Ttotal,dt,Nnodes)

        # calculating using fourth order nonlocal kernel
        u_4 = nonlocal_kernel_result.nonlocal_kernel_middisplacement(fourth_order_kernel)

        # This is the result calculated by using FEM with extremely fine mesh, we use it as the 'accurate' resulf for comparison
        ureference = np.loadtxt('FEM.dat')

        error_mid = 0
        for tt in range(Tsteps):
            mid_diff = abs(ureference[tt] - u_4[tt])
            if mid_diff > error_mid:
                error_mid = mid_diff
        
        assert error_mid<1.5


# this will be our displacement BC function
def displacement_load(t):
    T = 0.157e-3
    P0 = -50e3
    b0 = 1e46
    load = P0*b0*t**6*(t-T)**6*(1-np.heaviside(t-T,0))
    return load





if __name__ == '__main__':
    unittest.main()