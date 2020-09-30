import numpy as np
from scipy.optimize import optimize
import matplotlib.pylab as plt
import math

h_bar = 1.05457180013e-34
h = 6.62607004081e-34
"CODATA Recommended Value"

def parse_input():
    f = open("input.txt", 'r')
    for lines in f.readlines():
        if 'Format' in lines:
            continue
        else:
            print(lines)


def reduce_mass(m1, m2):
    reduced = ((m1*m2)/(m1+m2))
    return reduced


def energy_levels(force_constant, reduced_mass, n):
    """
    :param force_constant:
    The force constant from the given data.
    :param reduced_mass:
    The reduced mass, calculated in the reduce_mass function.
    :param n:
    :return:
    """
    En = h_bar * (force_constant/reduced_mass)**1/2 * (n + 1/2)
    return En


def frequency(force_constant, reduced_mass):
    v = 1/(2*math.pi) * (force_constant/reduced_mass)**1/2
    return v


def disassociation_energy(v, En, n):
    """
    :param v:
    This is the frequency determined from frequency()
    :param En:
    These are the energy levels determined in energy_levels().
    :param En:
    n is 0, 1, 2, 3...
    :return:
    This will be the disassociation energy.
    """
    energy = ((n+0.5)(h*v)**2)/(4*(h*v*(n+0.5)-En))
    return energy


def alpha(force_constant, d_energy):
    return force_constant/(2*d_energy)


def morse_potential(De, alpha, x, u):
    function = De*(1-np.exp(-alpha(x-u))**2)


parse_input()

