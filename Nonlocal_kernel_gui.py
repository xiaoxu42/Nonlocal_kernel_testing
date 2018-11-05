import matplotlib
matplotlib.use("TkAgg")

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

from nonlocal_kernel import Kernel_calculator as kc
from nonlocal_kernel_simulation import simulator_1D as sim1D


import tkinter as tk
from tkinter import ttk

LARGE_FONT = ("Verdana",12)

class Mygui(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, *kwargs)

        tk.Tk.wm_title(self,"Testing")

        container = tk.Frame(self)
        container.pack(side="top", fill="both",expand= True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (StartPage, PageOne, PageTwo):

            frame = F(container, self)
            frame.grid(row=0, column=0, sticky = "nsew")
            self.frames[F] = frame


        self.show_frame(StartPage)

    def show_frame(self,index):

        frame = self.frames[index]
        frame.lift()



class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self,text = "Welcome", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self,text = "Kernel Generator", command=lambda: controller.show_frame(PageOne))
        button1.pack()

        button2 = ttk.Button(self,text = "1-D Simulation", command=lambda: controller.show_frame(PageTwo))
        button2.pack()

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

        button3 = ttk.Button(self,text = "Generate discrete nonlocal kernel", command=lambda:self.GenerateKernel())
        button3.grid(row=3,column=1)



        button1 = ttk.Button(self,text = "1-D Simulation", command=lambda: controller.show_frame(PageTwo))
        button1.grid(row=20,column=1)

        button2 = ttk.Button(self,text = "Return", command=lambda: controller.show_frame(StartPage))
        button2.grid(row=30,column= 1)

    def GenerateKernel(self):
        order = int(self.e1.get())
        tolerance = float(self.e2.get())
        nonlocal_kernel = kc(200,5,8000,8000)
        # Pops out a new window newwin
        newwin = tk.Toplevel(self)
        
        self.upperbound = None # This upperbound can be changed by user throught advanced setting
        
        kernel = nonlocal_kernel.kernel_generator(order, tolerance,upperbound = self.upperbound)
        # Display the resulf of discrete kenerl
        kernel_string = "Discrete Kernel:"

        for x in kernel:
            kernel_string = kernel_string + "     " + str(round(x,4)) + "," 


        label = tk.Label(newwin, text = kernel_string)
        label.pack(side = "top")


        #Also generate a plot of Fourier transform function which generates the discrete kernel,
        # this can be used to check if the result is accurate or not to some extent 

        
        f =  nonlocal_kernel.plot_test(order, upperbound = self.upperbound)

        canvas = FigureCanvasTkAgg(f,newwin)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, newwin)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        button1 = ttk.Button(newwin, text="Advanced Setting", command=lambda: self.AdvancedSetting())
        button1.pack(side=tk.RIGHT)

        button2 = ttk.Button(newwin, text="Quit", command=newwin.destroy)
        button2.pack(side =tk.BOTTOM)

    def AdvancedSetting(self):
        newwin_ad = tk.Toplevel(self)

        label1 = tk.Label(newwin_ad, text="Upperbound = ")
        label1.grid(row=0,column=0,sticky = "W")

        e_ad = tk.Entry(newwin_ad)
        e_ad.grid(row=0,column=1,sticky = "W")
        e_ad.insert(10,"0.4")

        button1 = ttk.Button(newwin_ad, text = "Update the upperbound", command =lambda: self.SettingUpperbound(e_ad))
        button1.grid(row=1, column = 1, sticky = "E")


        button2 = ttk.Button(newwin_ad, text="Quit", command=newwin_ad.destroy)
        button2.grid(row = 10, sticky = "S")

    def SettingUpperbound(self,e):
        self.upperbound = float(e.get())        



        

class PageTwo(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self,text = "1-D Simulation", font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1 = ttk.Button(self,text = "Kernel Generator", command=lambda: controller.show_frame(PageOne))
        button1.pack()

        button2 = ttk.Button(self,text = "Return", command=lambda: controller.show_frame(StartPage))
        button2.pack()


myapp = Mygui()
myapp.geometry("1000x500")
myapp.mainloop()