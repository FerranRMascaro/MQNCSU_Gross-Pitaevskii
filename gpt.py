import math

# Constants
n1 = 1000
xr = [0.0] * (n1 + 1)
frev = [0.0] * (n1 + 1)
freo = [0.0] * (n1 + 1)
fred = [0.0] * (n1 + 1)
xmu = [0.0] * (n1 + 1)
fren = [0.0] * (n1 + 1)
den = [0.0] * (n1 + 1)
u = [0.0] * (n1 + 1)

# Input values
a0, n1, step, aa, time, alpha, iter = map(float, input("Enter input values: ").split())

# Constants
pi = math.pi
piin = 1.0 / (4.0 * pi)
pi2in = math.sqrt(piin)
alpha2 = alpha * alpha
cvar = 2.0 * math.sqrt(alpha) ** 3 / math.sqrt(math.sqrt(pi))

# Building the starting wave function R(r). Phi(r) = R(r)/r * Y00
for i in range(1, n1 + 1):
    xr[i] = step * (i - 1)
    xr2 = xr[i] * xr[i]
    frev[i] = cvar * xr[i] * math.exp(-0.5 * alpha2 * xr2)
    freo[i] = frev[i]

# Starting the convergence process
# To reduce to the harmonic oscillator, set cequ=0
cequ = a0 * aa
as3n = aa * a0 * a0 * a0
cequ = 0.0  # Uncomment this line to recover the harmonic oscillator

itw = 0
for it in range(1, iter + 1):
    itw += 1
    xnorm = 0.0
    ene0 = 0.0

    # Calculation of the second derivative using finite difference method
    for i in range(2, n1):
        fred[i] = (freo[i-1] + freo[i+1] - 2.0 * freo[i]) / (step * step)

    fred[n1] = (freo[n1-1] - 2.0 * freo[n1]) / (step * step)

    for i in range(1, n1 + 1):
        xr2 = xr[i] * xr[i]
        if i == 1:
            xmu[i] = 0.0
        else:
            # Calculation of energies and xmu
            ene0 = ene0 - freo[i] * fred[i] * 0.5 + 0.5 * xr2 * freo[i] * freo[i] + 0.5 * cequ * xr2 * (freo[i] / xr[i]) ** 4
            xmu[i] = -0.5 * fred[i] / freo[i] + 0.5 * xr2 + cequ * (freo[i] / xr[i]) ** 2

        # Update the wave function
        fren[i] = freo[i] - time * xmu[i] * freo[i]
        xnorm += fren[i] * fren[i]

    xnorm = math.sqrt(xnorm * step)
    ene0 = ene0 * step

    if itw == 200:
        print('ene0')

# Update the wave function normalization
for i in range(1, n1 + 1):
    freo[i] = fren[i] / xnorm

if it == iter:
    with open('output.txt', 'w') as f:
        for i in range(1, n1 + 1):
            f.write(f"{xr[i]:15.5e} {xmu[i]:15.5e}\n")

# Calculation of radius, potential, kinetic energy, density, and single-particle potential
radious = 0.0
xkin = 0.0
poth0 = 0.0
potself = 0.0
chem = 0.0
xaver = 0.0
xnormden = 0.0

for i in range(2, n1 + 1):
    xr2 = xr[i] * xr[i]
    radious += xr2 * freo[i] * freo[i]
    xkin += freo[i] * fred[i]
    poth0 += xr2 * freo[i] * freo[i]
    potself += xr2 * (freo[i] / xr[i]) ** 4
    chem += xmu[i] * freo[i] * freo[i]
    u[i] = 0.5 * xr2 + cequ * (freo[i] / xr[i]) ** 2
    den[i] = (freo[i] / xr[i]) ** 2
    xnormden += den[i] * xr2
    xaver += freo[i] * freo[i] * as3n * den[i]

radious2 = radious * step
radious = math.sqrt(radious * step)
xaver = xaver * step
chem = chem * step
xkin = -xkin * step * 0.5
poth0 = 0.5 * poth0 * step
potself = potself * step * cequ * 0.5
pot = potself + poth0
xnormden = xnormden * step

print('xnormden =', xnormden)
print('ene0 =', ene0)
print('chemical =', chem)
print('kinetic energy =', xkin)
print('potential energy =', pot)
print('potho =', poth0)
print('potself =', potself)
print('radius =', radious)
print('radius2 =', radious2)

with open('density.txt', 'w') as f:
    for i in range(2, n1 + 1):
        f.write(f"{xr[i]:15.5e} {den[i]:15.5e}\n")