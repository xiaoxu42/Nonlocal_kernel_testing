
# coding: utf-8

# In[3]:


from cffi import FFI
ffibuilder = FFI()

ffibuilder.set_source("explicit_time_api",
                     r"""
                     void explicit_time(double *u, int N, int Tsteps, const double kernel0, const double *kernel, const int horizon, const double E, const double rho, const double dt, const double *loadfunc);
                     """,sources=['./explicit_time.c']
                     )

ffibuilder.cdef("""
    void explicit_time(double *u, int N, int Tsteps, const double kernel0, const double *kernel, const int horizon, const double E, const double rho, const double dt, const double *loadfunc);
    """)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

