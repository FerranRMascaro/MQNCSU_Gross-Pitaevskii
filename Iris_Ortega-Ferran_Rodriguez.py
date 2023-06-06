import numpy as np 
import math as m
import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog = 'Iris_Ortega-Ferran_Rodriguez.py', description = 'Solves the Gross-Pitaevskii equation given some parameters for bosons in a spherical trap.')
    parser.add_argument('inpath', help = 'Path where the input file is stored')
    parser.add_argument('-ao', '--armonic_osci', action='store_true' , help='To reduce solution to the armonic oscillator (cequ=0)')

    arguments = parser.parse_args()
    inpath = arguments.inpath
    
    # Extracting data
    indata = pd.read_csv(inpath, sep = ' ', names=['a0','N_steps','step','N','time', 'alpha', 'iter'], header=None)
    alpha2_arr = np.array(indata['alpha'])**2
    cvar_arr = 2.0*np.sqrt(indata['alpha'])**3/m.sqrt(m.sqrt(m.pi))
    pd_res = pd.DataFrame(columns=['N', 'xnormden','ene', 'avg_chem_pot', 'kin_ene', 'total_pot', 'poth0', 'potint', 'radious', 'radious2'], index = [i for i in range(len(alpha2_arr))])

    #  num_run = number of runs or lines in input file
    #  building the starting wave function R(r). Phi(r)=R(r)/r*Y00

    log_path = '\\'.join(inpath.split('\\')[0:-1]) + '\\'
    if arguments.armonic_osci:
        log_path = log_path + 'ao_'
        cequ = 0

    for num_run in range(len(indata['N'])):
        xr, frev, freo, fred, xmu, fren, den, u = np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000), np.zeros(1000)
        # N_steps=n1, N=aa
        [a0, N_steps, step, N, time, alpha, iter] = indata.iloc[num_run]
        print(N, N_steps, step, a0, alpha, time, iter)
        alpha2 = alpha2_arr[num_run]
        cvar = cvar_arr[num_run]
        print('')
        print('---------------------------------------------------------')
        print('Computing for', int(N), 'number of particles.')
        N_steps, N = int(N_steps), int(N)

        # Log where we will write the output
        log_mu = open(log_path + str(N) + '_mu.txt', 'w+')
        log_den = open(log_path + str(N) + '_den.txt', 'w+')
        
        for i in range(N_steps):
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
        as3n = N*a0**3
        itw = 0
        for it in range(int(iter)):
            itw = itw + 1
            xnorm = 0
            ene0 = 0
            for i in range(1, N_steps - 1):
                fred[i] = (freo[i-1]+freo[i+1]-2.0*freo[i])/(step**2)
            fred[N_steps - 1] = (freo[N_steps-2]-2.0*freo[N_steps-1])/(step**2)
            xmu[0] = 0.0
            for i in range(N_steps):
                xr2=xr[i]**2
                if i != 0:
                    ene0 = ene0 + 0.5*(-freo[i]*fred[i] + xr2*freo[i]**2 + cequ*xr2*(freo[i]/xr[i])**4)
                    xmu[i] = 0.5*(xr2-fred[i]/freo[i]) + cequ*(freo[i]/xr[i])**2
                fren[i] = freo[i] - time*xmu[i]*freo[i]
                xnorm += fren[i]**2
            xnorm = m.sqrt(xnorm*step)
            ene0 = ene0*step
            if (itw % 200) == 0:
                print('ene0 =', ene0)      
            # I define the new wf.
            for i in range(N_steps):
                freo[i] = fren[i]/xnorm
        for i in range(1, N_steps):
            log_mu.write(str(xr[i]) + '\t' + str(xmu[i]) + '\n')
        log_mu.close()

        # calculation ofthe radious, potential and kinetic energy, density and single particle potential
        for i in range(1, N_steps - 1):
            fred[i] = (freo[i-1] + freo[i + 1]-2*freo[i])/(step**2)
        fred[N_steps - 1] = (freo[N_steps-2]-2*freo[i])/(step**2)

        radious, xkin, poth0, potself, chem, xaver, xnormden = 0, 0, 0, 0, 0, 0, 0
        
        for i in range(1, N_steps):
            xr2 = xr[i]**2
            radious = radious + xr2*freo[i]**2
            xkin = xkin + freo[i]*fred[i]
            poth0 = poth0 + xr2*freo[i]**2
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
        xkin = -xkin*step/2
        poth0 = poth0*step/2
        potself = potself*step*cequ/2
        pot = potself + poth0
        xnormden = xnormden*step
        for i in range(1, N_steps):
            log_den.write(str(xr[i]) + '\t' + str(den[i]) + '\n')
        log_den.close()
        
        pd_res.loc[num_run] = [N, xnormden, ene0, chem, xkin, pot, poth0, potself, radious, radious2]

    pd_res.to_csv(log_path + 'results.csv', header = True, index = False)