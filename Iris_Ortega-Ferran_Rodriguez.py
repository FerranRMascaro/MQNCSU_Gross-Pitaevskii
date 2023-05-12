import numpy as np 
import math as m
import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = 'Iris_Ortega-Ferran_Rodriguez.py', description = 'Computes')
    parser.add_argument('inpath', help = 'path')
    parser.add_argument('-ao', '--armonic_osci', action='store_true' , help='to reduce to the armonic oscillator (cequ=0)')

    arguments = parser.parse_args()
    inpath = arguments.inpath

    # Extracting data
    indata = pd.read_csv(inpath, sep = ' ', names=['a0','N_steps','step','N','alpha','time','iter'], header=None)
    alpha2_arr = np.array(indata['alpha'])**2
    cvar_arr = 2.0*np.sqrt(indata['alpha'])**3/m.sqrt(m.sqrt(m.pi))

    #  num_run = number of runs or lines in input file
    #  building the starting wave function R(r). Phi(r)=R(r)/r*Y00
    for num_run in range(len(indata['N'])):
        xr, frev, freo, fred, xmu, fren, den, u = np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000)

        [a0, N_steps, step, N, alpha, time, iter] = indata.iloc[num_run]
        alpha2 = alpha2_arr[num_run]
        cvar = cvar_arr[num_run]
        for i in range(int(N_steps)):
            xr[i] = step*i
            xr2 = xr[i]**2
            frev[i] = cvar*xr[i]*m.e**(-0.5*alpha2*xr2)
            freo[i] = frev[i]

    # starting the convergence process
    # ****************************************
    # to reduce to the armonic oscillator cequ=0
    # ************************************
        if not arguments.armonic_osci:
            cequ = a0*N
        else:
            cequ = 0
        as3n = N*a0**3
        # cequ = 0
        itw = 0
        for it in range(int(iter)):
            itw = itw + 1
            xnorm = 0
            ene0 = 0
            for i in range(2, N_steps):
                fred[i] = (freo[i-1]+freo[i+1]-2.0*freo[i])/(step*step)
            fred[N_steps]=(freo[N_steps-1]-2.0*freo[i])/(step**2)
            
        
    