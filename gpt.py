import math

# Constants
num_points = 1000
x = [0.0] * (num_points + 1)
real = [0.0] * (num_points + 1)
imag = [0.0] * (num_points + 1)
d2r = [0.0] * (num_points + 1)
mu = [0.0] * (num_points + 1)
new_real = [0.0] * (num_points + 1)
density = [0.0] * (num_points + 1)
potential = [0.0] * (num_points + 1)

# Input values
a0, num_points, step_size, aa, time_step, alpha, iterations = map(float, input("Enter input values: ").split())

# Constants
pi = math.pi
pi_inv = 1.0 / (4.0 * pi)
sqrt_pi_inv = math.sqrt(pi_inv)
alpha_squared = alpha * alpha
c = 2.0 * math.sqrt(alpha) ** 3 / math.sqrt(math.sqrt(pi))

# Building the starting wave function R(r). Phi(r) = R(r)/r * Y00
for i in range(1, num_points + 1):
    x[i] = step_size * (i - 1)
    x_squared = x[i] * x[i]
    real[i] = c * x[i] * math.exp(-0.5 * alpha_squared * x_squared)
    imag[i] = real[i]

# Starting the convergence process
# To reduce to the harmonic oscillator, set cequ=0
cequ = a0 * aa
a_s3n = aa * a0 * a0 * a0
# Uncomment this line to recover the harmonic oscillator
# cequ = 0.0

itw = 0
for it in range(1, iterations + 1):
    itw += 1
    xnorm = 0.0
    ene0 = 0.0

    # Calculation of the second derivative using finite difference method
    for i in range(2, num_points):
        d2r[i] = (imag[i-1] + imag[i+1] - 2.0 * imag[i]) / (step_size * step_size)

    d2r[num_points] = (imag[num_points-1] - 2.0 * imag[num_points]) / (step_size * step_size)

    for i in range(1, num_points + 1):
        x_squared = x[i] * x[i]
        if i == 1:
            mu[i] = 0.0
        else:
            # Calculation of energies and xmu
            ene0 = ene0 - imag[i] * d2r[i] * 0.5 + 0.5 * x_squared * imag[i] * imag[i] + 0.5 * cequ * x_squared * (imag[i] / x[i]) ** 4
            mu[i] = -0.5 * d2r[i] / imag[i] + 0.5 * x_squared + cequ * (imag[i] / x[i]) ** 2

        # Update the wave function
        new_real[i] = imag[i] - time_step * mu[i] * imag[i]
        xnorm += new_real[i] * new_real[i]

    xnorm = math.sqrt(xnorm * step_size)
    ene0 = ene0 * step_size

    if itw == 200:
        print('ene0')

# Update the wave function normalization
for i in range(1, num_points + 1):
    imag[i] = new_real[i] / xnorm

if it == iterations:
    with open('output.txt', 'w') as f:
        for i in range(1, num_points + 1):
            f.write(f"{x[i]:15.5e} {mu[i]:15.5e}\n")

# Calculation of radius, potential, kinetic energy, density, and single-particle potential
radius = 0.0
kinetic_energy = 0.0
potential_harmonic = 0.0
potential_self = 0.0
chemical_potential = 0.0
average_x = 0.0
density_norm = 0.0

for i in range(2, num_points + 1):
    x_squared = x[i] * x[i]
    radius += x_squared * imag[i] * imag[i]
    kinetic_energy += imag[i] * d2r[i]
    potential_harmonic += x_squared * imag[i] * imag[i]
    potential_self += x_squared * (imag[i] / x[i]) ** 4
    chemical_potential += mu[i] * imag[i] * imag[i]
    density[i] = (imag[i] / x[i]) ** 2
    density_norm += density[i] * x_squared
    average_x += imag[i] * imag[i] * a_s3n * density[i]

radius_squared = radius * step_size
radius = math.sqrt(radius * step_size)
average_x = average_x * step_size
chemical_potential = chemical_potential * step_size
kinetic_energy = -kinetic_energy * step_size * 0.5
potential_harmonic = 0.5 * potential_harmonic * step_size
potential_self = potential_self * step_size * cequ * 0.5
potential_total = potential_self + potential_harmonic
density_norm = density_norm * step_size

print('Density normalization:', density_norm)
print('Ene0:', ene0)
print('Chemical potential:', chemical_potential)
print('Kinetic energy:', kinetic_energy)
print('Potential energy:', potential_total)
print('Potential harmonic:', potential_harmonic)
print('Potential self:', potential_self)
print('Radius:', radius)
print('Radius squared:', radius_squared)

with open('density.txt', 'w') as f:
    for i in range(2, num_points + 1):
        f.write(f"{x[i]:15.5e} {density[i]:15.5e}\n")