# -*- coding: utf-8 -*-
"""
Method of Characteristics Nozzle Design Code
@author: Peter Senior
email: pws1g13@soton.ac.uk
"""

import numpy
from numpy import sin, tan, arcsin, arctan, pi
import matplotlib.pyplot as plt

def mu(M):
    """Calculates and returns the value of mu for a given mach number"""
    return arcsin(1 / M)

def nu(M, gamma = 1.4):
    """Calculates and returns the value of nu (the Prantl-Meyer angle) based
    upon the Mach number at a point"""
    nux = (((gamma + 1) / (gamma - 1)) ** 0.5) * arctan(((gamma - 1) * (M ** 2
    - 1) / (gamma + 1)) ** 0.5) - arctan((M ** 2 - 1) ** 0.5) 
    return nux


def R_over_Re(M, Nu_exit, gamma=1.4):
    Rx = (1 - (((2 / (gamma + 1)) * (1 + (gamma - 1) * M ** 2 / 2)) ** 
    ((gamma + 1) / (2 * (gamma - 1))) * sin(phi(M))) / Eps) ** 0.5
    return Rx


def phi(M):
    return Nu_exit - nu(M) + mu(M)

N = 50 #solution resolution

#Constants
Mi = 1.
Me = 4.
#ht = 1.
gamma = 1.4
Re = 10 #Cowl Radius

M = numpy.linspace(1, Me, N)

Eps = (1 / Me) * ((2 / (gamma + 1)) * (1 + (gamma - 1) * Me ** 2 / 2)) ** ((
gamma + 1) / (2 * (gamma - 1)))

Nu_exit = nu(Me)
delta = pi / 2 - Nu_exit

#Re = ht * Eps * sin(delta) / (Eps - (Eps * (Eps - sin(delta))) ** 0.5)
ht = Re * Eps * sin(delta) / (Eps - (Eps * (Eps - sin(delta))) ** 0.5)

x = numpy.linspace(0, 1, N)

Min = (Me - 1) / N
Mx = 1 + x * Min

R = R_over_Re(M, Nu_exit)
X = Re * (1 - R) / tan(phi(M))

plt.plot(X / Re, R * 10 / Re)
plt.plot(0, 1, marker='o')
