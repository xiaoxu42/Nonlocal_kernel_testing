import matplotlib
matplotlib.use("TkAgg")
#import matplotlib.pyplot as plt

import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

from nonlocal_kernel import Kernel_calculator as kc
from nonlocal_kernel_simulation import simulator_1D as sim1D


import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

import os

LARGE_FONT = ("Verdana",12)
NORMAL_FONT = ("Verdana",10)

# funcction that convert array to string that will be used to display the discrete kernels
def DisplayArray(Array):
    Array_string = ""
    for x in Array:
        Array_string = Array_string + "     " + str(round(x,4)) + "," 
    
    return Array_string

def popupmsg(msg):
    popup = tk.Tk()

    popup.wm_title("Message")
    label = ttk.Label(popup, text =msg, font = NORMAL_FONT)
    label.pack(side = "top", fill ="x", pady =10)
    buttonn1 = ttk.Button(popup, text = "Got it", command = popup.destroy)
    buttonn1.pack()

    popup.mainloop()

class Mygui(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, *kwargs)

        tk.Tk.wm_title(self,"Testing")

        container = tk.Frame(self)
        container.pack(side="top", fill="both",expand= True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}
        self.discrete_kernel = None
        self.discrete_kernel_flag = 0 # flag that represents whether the discrete kernel has been assigned value or not
        self.displacement_load_data = None
        self.displacement_load_data_flag = 0 # flag that represents wheter the BC array has been assigned or not

        for F in (StartPage, PageOne, PageTwo, SimbyUserKernel):

            frame = F(container, self)
            frame.grid(row=0, column=0, sticky = "nsew")
            self.frames[F] = frame


        self.show_frame(StartPage)

        

        # This is the default material properties used for calculating the nonlocal kernel(default initialization)
        self.nonlocal_kernel = kc(200e9,5e9,8000,8000)
        
    def show_frame(self,index):

        frame = self.frames[index]
        frame.lift()
    
    def initialization(self,e1,e2,e3,e4):
        # This is used to re-initialize the basis 1-D elastic problem using nonlocal_kernel library from user's input
        E1 = float(e1.get())
        E2 = float(e2.get())
        rho1 = float(e3.get())
        rho2 = float(e4.get())

        if E1>E2:
            Eh, Es = E1, E2
            rhoh, rhos = rho1, rho2
        else:
            Eh, Es = E2, E1
            rhoh, rhos = rho2, rho1

        self.nonlocal_kernel = kc(Eh,Es,rhoh,rhos)
       

        




class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self,text = "Welcome, please specify the material properties of the composit", font=LARGE_FONT)
        label.grid(row=0,column = 0, columnspan = 3, pady=10,padx=10)

        label1 = tk.Label(self, text = "Elastic modulus of material 1 (Pa)")
        label1.grid(row=1,column=0)

        e1 = tk.Entry(self)
        e1.grid(row=1,column=1)
        e1.insert(10,"200e9")

        label2 = tk.Label(self,text = "Elastic modulus of material 2 (Pa)")
        label2.grid(row=2,column=0)

        e2 = tk.Entry(self)
        e2.grid(row=2,column=1)
        e2.insert(10,"5e9")

        label3 = tk.Label(self,text = "Density of material 1 (kg/m3)")
        label3.grid(row=3,column=0)

        e3 = tk.Entry(self)
        e3.grid(row=3,column=1)
        e3.insert(10,"8000")

        label4 = tk.Label(self,text = "Density of material 2 (kg/m3)")
        label4.grid(row=4,column=0)

        e4 = tk.Entry(self)
        e4.grid(row=4,column=1)
        e4.insert(10,"8000")

        

        button = ttk.Button(self, text = "OK", command=lambda: controller.initialization(e1,e2,e3,e4))
        button.grid(row=5,column = 2, sticky = "E", pady = 20)

        button1 = ttk.Button(self,text = "Kernel Generator", command=lambda: controller.show_frame(PageOne))
        button1.grid(row = 10,column=0,sticky = "W")

        button3 = ttk.Button(self,text = "Nonlocal 1-D Simulation using Kernel Generator", command=lambda: controller.show_frame(SimbyUserKernel))
        button3.grid(row = 11,column=0, sticky = "W")

        button2 = ttk.Button(self,text = "Quick Nonlocal 1-D Simulation", command=lambda: controller.show_frame(PageTwo))
        button2.grid(row = 12,column=0, sticky = "W")
    

class PageOne(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self,text = "Kernel Generator", font=LARGE_FONT)
        label.grid(row = 0,column = 1)

        label1  = tk.Label(self,text = "Calculating the kernel using")
        label1.grid(row=1,column=0)

        self.e1 = tk.Entry(self)
        self.e1.grid(row=1,column=1)
        self.e1.insert(10,"4")

        label2 = tk.Label(self,text = "order approximation")
        label2.grid(row=1,column=2,sticky ="W")

        label3 = tk.Label(self,text = "with error tolerance set to be")
        label3.grid(row=2,column=0)

        self.e2 = tk.Entry(self)
        self.e2.grid(row=2,column=1)
        self.e2.insert(10,"0.001")

        button3 = ttk.Button(self,text = "Generate discrete nonlocal kernel", command=lambda:self.GenerateKernel(controller))
        button3.grid(row=3,column=1)



        button1 = ttk.Button(self,text = "Quick Nolocal 1-D Simulation", command=lambda: controller.show_frame(PageTwo))
        button1.grid(row=20,column=0, sticky = "WS")

        button2 = ttk.Button(self,text = "Return", command=lambda: controller.show_frame(StartPage))
        button2.grid(row=30,column=0, sticky = "WS")

    def GenerateKernel(self,controller):
        self.order = int(self.e1.get())
        self.tolerance = float(self.e2.get())
        nonlocal_kernel = controller.nonlocal_kernel
        # Pops out a new window newwin
        newwin = tk.Toplevel(self)

        self.upperbound = None # This upperbound can be changed by user throught advanced setting
        
        kernel = nonlocal_kernel.kernel_generator(self.order, self.tolerance,upperbound = self.upperbound)
        self.discrete_kernel = kernel
        # Display the resulf of discrete kernel
        kernel_string = "Discrete Kernel:" + DisplayArray(kernel)

        # This is the label displaying the the resulf of dicrete kernel
        self.kernel_label = tk.Label(newwin, text = kernel_string)
        self.kernel_label.pack(side = tk.TOP,expand=True, pady = 15)

        # This is the button that will pass the calculated discrete kernel to controller so that it can be used in other pages for simluation later on
        buttonpasskernel = ttk.Button(newwin, text = "Use this kernel to do 1-D nonlocal simulation", command=lambda: self.PassDiscreteKernel(controller,newwin))
        buttonpasskernel.pack()


        #Also generate a plot of Fourier transform function which generates the discrete kernel,
        # this can be used to check if the result is accurate or not to some extent 

        
        (x,y) =  nonlocal_kernel.plot_test(self.order, upperbound = self.upperbound)
        self.f = matplotlib.figure.Figure(figsize=(5,4),dpi=200)
        self.a = self.f.add_subplot(111)
        self.a.plot(x,y)

        self.canvas1 = FigureCanvasTkAgg(self.f,newwin)
        self.canvas1.draw()
        self.canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(self.canvas1, newwin)
        toolbar.update()
        self.canvas1._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        button1 = ttk.Button(newwin, text="Advanced Setting", command=lambda: self.AdvancedSetting(controller))
        button1.pack(side=tk.RIGHT)

        button2 = ttk.Button(newwin, text="Quit", command=newwin.destroy)
        button2.pack(side =tk.BOTTOM)

    def AdvancedSetting(self,controller):

        
        newwin_ad = tk.Toplevel(self)

        label1 = tk.Label(newwin_ad, text="Upperbound = ")
        label1.grid(row=0,column=0,sticky = "W")

        e_ad = tk.Entry(newwin_ad)
        e_ad.grid(row=0,column=1,sticky = "W")
        e_ad.insert(10,"0.4")

        button1 = ttk.Button(newwin_ad, text = "Update the upperbound and Refresh the plot", command =lambda: self.SettingUpperboundandRefresh(e_ad,controller))
        button1.grid(row=1, column = 1, sticky = "E")


        button2 = ttk.Button(newwin_ad, text="Quit", command=newwin_ad.destroy)
        button2.grid(row = 10, sticky = "S")

    def PassDiscreteKernel(self, controller,newwin):
        controller.discrete_kernel = self.discrete_kernel
        controller.discrete_kernel_flag = 1
        newwin.destroy()
        controller.show_frame(SimbyUserKernel)

        

    def SettingUpperboundandRefresh(self,e,controller):
        # This function is used to update the plot on canvas1 using the parameter set by the user

        self.upperbound = float(e.get())
        nonlocal_kernel = controller.nonlocal_kernel

        kernel = nonlocal_kernel.kernel_generator(self.order, self.tolerance,upperbound = self.upperbound)
        self.kernel = kernel
        # Update the displayed discrete kernel result
        kernel_string = "Discrete Kernel:" + DisplayArray(kernel)

        self.kernel_label.configure(text = kernel_string)

        (x,y) =  nonlocal_kernel.plot_test(self.order, upperbound = self.upperbound)
        self.a.clear()
        self.a.plot(x,y)
        self.canvas1.draw()



        

class PageTwo(tk.Frame):
    # This is the page where you can do a quick nonlocal simlution by directly specify the number of oder of the kernel you want to use

    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent)

        label0 = tk.Label(self, text = "Order of the kernel:" )
        label0.grid(row=1,column = 0,sticky = "E")

        self.e0 = tk.Entry(self)
        self.e0.grid(row=1,column=1,sticky = "W")
        self.e0.insert(10,"4")

        label1 = tk.Label(self, text = "Total simulation time:" )
        label1.grid(row=2,column = 0,sticky = "E")

        self.e1 = tk.Entry(self)
        self.e1.grid(row=2,column=1,sticky = "W")
        self.e1.insert(10,"1e-3")
        
        label2 = tk.Label(self, text = "Time step:")
        label2.grid(row=3,column=0, sticky = "E")

        self.e2 = tk.Entry(self)
        self.e2.grid(row=3,column=1,sticky = "W")
        self.e2.insert(10,"1e-7")

        label3 = tk.Label(self, text = "Number of nodes" )
        label3.grid(row=4,column = 0,sticky = "E")

        self.e3 = tk.Entry(self)
        self.e3.grid(row=4,column=1,sticky = "W")
        self.e3.insert(10,"50")

        self.button_BC = ttk.Button(self, text = "Import displacement load data from file...", command=lambda: self.Opendisplacementfile(controller))
        self.button_BC.grid(row=5,column = 0, columnspan = 2)
        

        button_sim = ttk.Button(self, text = "Do nonlocal simulation!!", command=lambda: self.Nonlocal_sim(controller))
        button_sim.grid(row=6,column=2)
        

        button1 = ttk.Button(self,text = "Kernel Generator", command=lambda: controller.show_frame(PageOne))
        button1.grid(row=10,column=0,sticky = "W")

        button2 = ttk.Button(self,text = "Nonlocal 1-D Simulation using Kernel Generator", command=lambda: controller.show_frame(SimbyUserKernel))
        button2.grid(row=11,column=0,sticky = "W")

        button3 = ttk.Button(self,text = "Return", command=lambda: controller.show_frame(StartPage))
        button3.grid(row=12,column=0,sticky = "W")

    
    def Opendisplacementfile(self,controller):
        name  = filedialog.askopenfilename(initialdir = os.getcwd(), filetypes =(("Data File","*.dat"),("All files","*.*")), title = "Choose a file")

        if name != "":       # if user didn't click Cancel 
            # Use try just in case user types some unknown file or closes without choosing a file
            controller.displacement_load_data = np.loadtxt(name) 
            controller.displacement_load_data_flag = 1
            # This loadtxt command actually will cause problem in debug mode, but it works fine without debugging. Refer to https://github.com/Microsoft/ptvsd/issues/465 for details about the problem of numpy.loadtxt
            popupmsg("Load successfully!")
           
    def Nonlocal_sim(self, controller):

        order = int(self.e0.get())
        Ttotal = float(self.e1.get())
        dt = float(self.e2.get())
        Nnodes = int(self.e3.get())
        
        nonlocal_kernel = controller.nonlocal_kernel
        
        tolerance = 0.01    #default tolerance for kernel generator

        Eh = nonlocal_kernel.Eh

        rhoave = nonlocal_kernel.rhoh*2*nonlocal_kernel.alpha + nonlocal_kernel.rhos*nonlocal_kernel.beta

        if controller.displacement_load_data_flag == 0:
            popupmsg("Oops, you haven't defined your Boundary Condition yet, please use the \"import displacement load data\" to define your BC ")
        
        quick_kernel = nonlocal_kernel.kernel_generator(order,tolerance)    
        nonlocal_kernel_result = sim1D(Eh,rhoave,controller.displacement_load_data,Ttotal,dt,Nnodes)
        u_mid = nonlocal_kernel_result.nonlocal_kernel_middisplacement(quick_kernel)
        
        # Pops out a new window newwin
        newwin = tk.Toplevel(self)

        size1 = np.size(u_mid)
        t1 = np.linspace(0.0,Ttotal,num = size1)
        
        self.f = matplotlib.figure.Figure(figsize=(5,4),dpi=200)
        self.a = self.f.add_subplot(111)
        self.a.plot(t1,u_mid)
        self.a.set_xlabel('Time ($s$)')
        self.a.set_ylabel('Midpoint Displacement ($m$)')

        self.canvas1 = FigureCanvasTkAgg(self.f,newwin)
        self.canvas1.draw()
        self.canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(self.canvas1, newwin)
        toolbar.update()
        self.canvas1._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


class SimbyUserKernel(tk.Frame):

    # This is the page users can do their 1-D simulation using the kernel they calculate at Kernel Generator

    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self, parent)

        self.label = tk.Label(self,text = "Import the discrete kernel to see if you are using the corret one!", font=LARGE_FONT)
        self.label.grid(row=0, column = 0, columnspan = 4, pady=10,padx=10)

        button = ttk.Button(self, text = "Import the discrete kernel from kernel generator", command = lambda:self.DisplayDiscreteKernel(controller))
        button.grid(row=1,column = 1)

        label1 = tk.Label(self, text = "Total simulation time:" )
        label1.grid(row=2,column = 0,sticky = "E")

        self.e1 = tk.Entry(self)
        self.e1.grid(row=2,column=1,sticky = "W")
        self.e1.insert(10,"1e-3")
        
        label2 = tk.Label(self, text = "Time step:")
        label2.grid(row=3,column=0, sticky = "E")

        self.e2 = tk.Entry(self)
        self.e2.grid(row=3,column=1,sticky = "W")
        self.e2.insert(10,"1e-7")

        label3 = tk.Label(self, text = "Number of nodes" )
        label3.grid(row=4,column = 0,sticky = "E")

        self.e3 = tk.Entry(self)
        self.e3.grid(row=4,column=1,sticky = "W")
        self.e3.insert(10,"50")

        self.button_BC = ttk.Button(self, text = "Import displacement load data from file...", command=lambda: self.Opendisplacementfile(controller))
        self.button_BC.grid(row=5,column = 0, columnspan = 2)
        

        button_sim = ttk.Button(self, text = "Do nonlocal simulation!!", command=lambda: self.Nonlocal_sim(controller))
        button_sim.grid(row=6,column=2,sticky = "W")
        

        button1 = ttk.Button(self,text = "Kernel Generator", command=lambda: controller.show_frame(PageOne))
        button1.grid(row=10,column=0,sticky = "W")

        button2 = ttk.Button(self,text = "Return", command=lambda: controller.show_frame(StartPage))
        button2.grid(row=11,column=0,sticky = "W")

    def DisplayDiscreteKernel(self,controller):

        if controller.discrete_kernel_flag == 0:
            popupmsg("Oops, you haven't calculated the kernel yet, please got back to Kernel Generator to calculate your kernel")
        
        else:
            kernel_string = "Your discrete kernel is:" + DisplayArray(controller.discrete_kernel)
            self.label.configure(text = kernel_string)
    
    def Opendisplacementfile(self,controller):
        name  = filedialog.askopenfilename(initialdir = os.getcwd(), filetypes =(("Data File","*.dat"),("All files","*.*")), title = "Choose a file")

        if name != "":       # if user didn't click Cancel 
            # Use try just in case user types some unknown file or closes without choosing a file
            controller.displacement_load_data = np.loadtxt(name) 
            controller.displacement_load_data_flag = 1
            # This loadtxt command actually will cause problem in debug mode, but it works fine without debugging. Refer to https://github.com/Microsoft/ptvsd/issues/465 for details about the problem of numpy.loadtxt
            popupmsg("Load successfully!")
        

        
            
    def Nonlocal_sim(self, controller):

        Ttotal = float(self.e1.get())
        dt = float(self.e2.get())
        Nnodes = int(self.e3.get())
        
        nonlocal_kernel = controller.nonlocal_kernel
        

        Eh = nonlocal_kernel.Eh

        rhoave = nonlocal_kernel.rhoh*2*nonlocal_kernel.alpha + nonlocal_kernel.rhos*nonlocal_kernel.beta

        if controller.displacement_load_data_flag == 0:
            popupmsg("Oops, you haven't defined your Boundary Condition yet, please use the \"import displacement load data\" to define your BC ")
        
            
        if controller.discrete_kernel_flag == 0:
            popupmsg("Oops, you haven't generate a kernel you want to use yet, please go back to Kernel Generator to generate your kernel")

            
        nonlocal_kernel_result = sim1D(Eh,rhoave,controller.displacement_load_data,Ttotal,dt,Nnodes)
        u_mid = nonlocal_kernel_result.nonlocal_kernel_middisplacement(controller.discrete_kernel)
        
        # Pops out a new window newwin
        newwin = tk.Toplevel(self)

        size1 = np.size(u_mid)
        t1 = np.linspace(0.0,Ttotal,num = size1)
        
        self.f = matplotlib.figure.Figure(figsize=(5,4),dpi=200)
        self.a = self.f.add_subplot(111)
        self.a.plot(t1,u_mid)
        self.a.set_xlabel('Time ($s$)')
        self.a.set_ylabel('Midpoint Displacement ($m$)')

        self.canvas1 = FigureCanvasTkAgg(self.f,newwin)
        self.canvas1.draw()
        self.canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(self.canvas1, newwin)
        toolbar.update()
        self.canvas1._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)










myapp = Mygui()
#myapp.geometry("1000x500")
myapp.mainloop()