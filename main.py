import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pylab as plt
import math

h_bar = 1.05457180013e-34  # J*s
h = 6.62607004081e-34      # J*s
"CODATA Recommended Value"
amu = 1.66054e-27          # u to kg


# 1 hartree bohr−2 = 15.569141 mdyne Å−1
# Double the top number in the molpro output.
# Convert to mDyne/A by the above conversion.
# Convert to N/m by multiplying mDyne/A by 100.


def parse_input():
    f = open("input.txt", 'r')
    for lines in f.readlines():
        if 'Format' in lines:
            continue
        else:
            inputs = lines.split(', ')
            k = float(inputs[0]) * 2
            k = k * 15.569141 * 100
            mass_1 = float(inputs[1]) * amu
            mass_2 = float(inputs[2]) * amu
            bond_length = float(inputs[3]) * 1e-10

            return k, mass_1, mass_2, bond_length


def reduce_mass(m1, m2):
    """
    :param m1:
    Mass of atom 1 in kg.
    :param m2:
    Mass of atom 2 in kg.
    :return:
    Yields reduced mass in kg.
    """
    reduced = ((m1 * m2) / (m1 + m2))
    return reduced


# Force Constants to be in N/m.
# n should be 0.
def energy_levels(force_constant, reduced_mass, n):
    """
    :param force_constant:
    The force constant from the given data.
    :param reduced_mass:
    The reduced mass, calculated in the reduce_mass function.
    :param n:
    :return:
    """
    En = h_bar * (force_constant / reduced_mass) ** 1 / 2 * (n + 1 / 2)
    return En


def frequency(force_constant, reduced_mass):
    v = 1 / (2 * math.pi) * (force_constant / reduced_mass) ** 1 / 2
    return v


def disassociation_energy(v, En, n):
    """
    :param v:
    This is the frequency determined from frequency()
    :param En:
    These are the energy levels determined in energy_levels().
    :param n:
    n is 0, 1, 2, 3...
    :return:
    This will be the disassociation energy.
    """
    energy = (-(n + 1/2) * (h * v) ** 2) / (4 * (En - h * v * (n + 1/2)))
    return energy


def alpha(force_constant, d_energy):
    return force_constant / (2 * d_energy)


def morse_potential(De, alpha, x, u):
    function = De * (1 - np.exp(-alpha*(x - u)) ** 2)
    yield function


if __name__ == '__main__':
    k, m1, m2, eq_length = parse_input()
    rM = reduce_mass(m1, m2)
    levels = energy_levels(k, rM, 0)
    freq = frequency(k, rM)
    De = disassociation_energy(freq, levels, 0)
    a = alpha(k, De)

    x_values = []
    y_values = []
    for i in np.arange(eq_length-1e-10, eq_length+1e-10, 0.005e-10):
        x_values.append(i)

    print(x_values)

    for i in x_values:
        for y in morse_potential(De, a, i, eq_length):
            y_values.append(y)

    print(y_values)

    plt.plot(x_values, y_values)
    plt.show()
