# Nonlocal Kernel Calculator & 1-D Nonlocal Simulator

This is initially the python code that I use to generate some plots in the paper that I am currently working on. Basically, it can generate any order of discrete nonlocal kernel based on material's microstructure of a 1-D elastic composite rod. I generate the nonlocal kernels of different orders in a systematic way. However, when order gets higher, the numerical errors seems to be amplified exponentially which will make the matrix inside the lineq_solver in [nonlocal_kernel](nonlocal_kernel.py) to be nearly singular.  This limits that, right now, the highest order of the nonlocal kernel that the code can generate accurately is only 6. There are actually two critical numerical parameters (upperbound and number of samples) inside the fitting_kernel_formular in [nonlocal_kernel](nonlocal_kernel.py) that may greatly affect the accuracy of the resulf if not chosen properly. Therefore, I have some keywords arguments in these methods which can be used to possibly get better results for different problems.

After calculating the discrete nonlocal kernel, 1-D nonlocal elastic simulation is done using teh corresponding discrete kernel. And the simulation process has been accelerate by implementing a CFFI API solution to use C to do the for loops during simulation

Then we compare the numerical result calculated using different order of nonlocal kernel with the classical peridynamic kernel as well as the "accurate" result by using FEM with extremely fine mesh. It shows that the nonlocal kernel it generates is much more accurate than classical peridynamic kernel when simulating the composite rod with heterogeneous microstructure.

# GUI by tkinter

Run [Nonlocal_kernel_gui](Nonlocal_kernel_gui.py) to start the GUI.

Finally, a GUI is coded using tkinter modulus which allows users to input their elastic properties to generate their own discrete nonlocal kernel. Furthermore, the gui provides advance setting where users can modify the two critical numerical parameters to achieve better accuracy, and also the plots of the inverse fourier transform is ploted so that users can examine the quality of the nonlocal kernel. With nonlocal kernel generated, users can use the kernel to do 1-D nonlocal kernel by inputing some numerical parameters as well as boundary condition. The plot of displacement at midpoint can be generated and user can save the complete dispalcement result at each time step as a .dat file. The GUI also provide a option that users can directly do a nonlocal simulation by simply specifying the order of the nonlocal kernel you want to use.

One thing to notice is that if you are running [Nonlocal_kernel_gui](Nonlocal_kernel_gui.py) in a debugging environment, errors will occur when running the numpy.loadtxt. A detailed information about the bug can be referred to https://github.com/numpy/numpy/issues/11538

Notice that this is the code that I mainly use for some numerical testing and plots. So it is not robust right now. But it is promising that this code can be developed further if all the numerical problems have been solved. 
