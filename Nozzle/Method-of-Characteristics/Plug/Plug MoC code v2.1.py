# -*- coding: utf-8 -*-
"""
Method of Characteristics Nozzle Design Code
@author: Peter Senior
email: pws1g13@soton.ac.uk
"""

import numpy
from numpy import sin, tan, arcsin, arctan
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


def R_over_Re(M, Nu_exit, Eps, gamma=1.4):
    Rx = (1 - (((2 / (gamma + 1)) * (1 + (gamma - 1) * M ** 2 / 2)) ** 
    ((gamma + 1) / (2 * (gamma - 1))) * sin(phi(M, Nu_exit))) / Eps) ** 0.5
    return Rx


def phi(M, Nu_exit):
    return Nu_exit - nu(M) + mu(M)


def MoC_Plug_Nozzle_Geometry(Me, gamma=1.4, N=50):
    """Creates plug nozzle geometry using Me as the exit mach number, gamma
    defaulting to 1.4, and N as the solution resolution, defaulting to 50. The
    expansion point is located at 0, 1. The returned coordinates are
    nondimensional"""
    #Constants
    Me = float(Me)
    gamma = 1.4

    M = numpy.linspace(1, Me, N)

    Eps = (1 / Me) * ((2 / (gamma + 1)) * (1 + (gamma - 1) * Me ** 2 / 2)) ** (
    (gamma + 1) / (2 * (gamma - 1)))

    Nu_exit = nu(Me)

    R = R_over_Re(M, Nu_exit, Eps)
    X = (1 - R) / tan(phi(M, Nu_exit))
    return [list(X), list(R)]


def Plot_Geo(Me, gamma=1.4, N=50):
    """Plots the created geometry"""
    coords = MoC_Plug_Nozzle_Geometry(Me, gamma, N)
    plt.plot(coords[0], coords[1])
    plt.plot(0,1, marker='o')
    return
