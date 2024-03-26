#!/usr/bin/python

#############################################################
# jiri kvita 2018

from math import *

from Consts import *

from MaterialConsts import *
from Particles import *

if len(gMaterials) == 0:
    MakeMaterials()
if len(gParticles) == 0:
    MakeParticles()

#############################################################
# hm, X0 could be computed for electrons as a material property
# and modified for other particles by (M/me)^2 on the fly?

def dEdXBrems(E, particle, material, verbose = 0):
    rho = material.GetRho()
    A = material.GetA()
    Z = material.GetZ()
    if Z < 0:
        print("ERROR: seems like this is a compound material for which X0 has to be computed from weights and individual X0's, feature not supported yet!")
        return 0.
    X0 = GetX0(particle.GetM(), particle.GetZ(), Z, A)
    if verbose: print(' X0={:f} cm'.format(X0/rho))
    dedx = rho*E / X0
    if verbose: print(' K={:f} Lrad={:f} Lrad\'={:f} Brems. -dE/dX ={:f}'.format(Konst,Lrad(Z),LradPrime(Z),dedx))
    return dedx


