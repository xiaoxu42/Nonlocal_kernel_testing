#include <math.h>

void explicit_time(double *u, int N, int Tsteps, const double kernel0, const double *kernel, const int horizon, const double E, const double rho, const double dt, const double *loadfunc){
    double ut1[N], ut2[N], utt[N] ;
    double temp ;
    int tt, i, k;

    for (i = 0; i<N; ++i){
        ut1[i]=0;
    }

    for (tt = 0; tt< Tsteps-1; ++tt){
        for (i = horizon; i< N-horizon; ++i){
            temp = -kernel0*u[i*Tsteps + tt];
            for (k=0; k<horizon; ++k){
                temp = temp + kernel[k]*(u[(i+k+1)*Tsteps + tt] + u[(i-k-1)*Tsteps + tt]);
            }
            utt[i] = E/rho*temp;
            ut2[i] = ut1[i]+dt*utt[i];
            u[i*Tsteps + tt+1] = u[i*Tsteps + tt] + dt*(ut1[i]+ut2[i])/2;
            ut1[i] = ut2[i];
        
        }
    
        for (i=0; i<horizon; i++){
            u[i*Tsteps + tt+1] = - u[(2*horizon-1-i)*Tsteps + tt+1];
        }

        for (i=N-horizon; i<N; ++i){
            u[i*Tsteps + tt+1] = 2*loadfunc[tt+1] - u[(2*(N-horizon)-1-i)*Tsteps + tt+1];

        }
    }

    return;

}
