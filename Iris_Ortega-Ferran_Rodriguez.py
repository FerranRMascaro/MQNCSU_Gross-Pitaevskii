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
        # N_steps=n1, N=aa
        [a0, N_steps, step, N, alpha, time, iter] = indata.iloc[num_run]
        alpha2 = alpha2_arr[num_run]
        cvar = cvar_arr[num_run]

        # Log where we will write the output
        log_path = '\\'.join(inpath.split('\\')[0:-1])
        log = open(log_path + str(N)  '.txt', 'w+')
        
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
            for i in range(1, N_steps - 1):
                fred[i] = (freo[i-1]+freo[i+1]-2.0*freo[i])/(step*step)
            fred[N_steps - 1]=(freo[N_steps-2]-2.0*freo[i])/(step**2)
            xmu[0]=0.0
            for i in range(N_steps):
                xr2=xr[i]**2
                if i != 0:
                    ene0 = 0.5*(ene0-freo[i]*fred[i] + xr2*freo[i]**2 + cequ*xr2*(freo[i]/xr[i])**4)
                    xmu[i] = 0.5*(xr2-fred[i]/freo[i]) + cequ*(freo[i]/xr[i])**2
                fren[i] = freo[i]-time*xmu[i]*freo[i]
                xnorm += fren[i]**2
            xnorm = m.sqrt(xnorm*step)
            ene0 = ene0*step
            if (itw % 200) == 0:
                print('ene0 =', ene0)      
            # I define the new wf.
            for i in range(N_steps):
                freo[i] = fren[i]/xnorm
            if it == iter:
              	log.write(xr[i], xmu[i], 'i=2', N_steps)
                
        log.close()
        
        # calculation ofthe radious, potential and kinetic energy, density and single particle potential
        for i in range(1, N_steps - 1):
            fred[i] = (freo[i-1] + freo[i + 1]-2*fro[i])/(step**2)
        fred[N_steps - 1] = (freo[N_steps-2]-2*freo[i])/(step**2)

        radious, xkin, potho, potself, chem, xaver, xnormden = 0, 0, 0, 0, 0, 0, 0
        
        for i in range(1, N_steps):
            xr2 = xr[i]**2
            radious = radious + xr2*freo[i]**2
            xkin = xkin + freo[i]*fred[i]
            poth0 = potho + xr2*freo[i]**2
            potself = potself + xr2*freo[i]**2
            potself = potself + xr2*(freo[i]/xr[i])**4
            chem = chem + xmu[i]*freo[i]**2
            u[i] = 0.5*xr2 + cequ*(freo[i]/xr[i])**2
            den[i] = (freo[i]/xr[i])**2
            xnormden = xnormden + den[i]*xr2
            xaver = xaver + (freo[i]**2)*as3n*den[i]
        radious2 = radious*step
        radious = m.sqrt(radious*step)
        xaver = xaver*step
        chem = chem*step
        xkin = xkin*step/2
        poth0 = poth0*step/2
        potself = potself*step*cequ/2
        pot = potself + poth0
        xnormden = xnormden*step
        
        