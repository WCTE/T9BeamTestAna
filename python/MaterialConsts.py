#!/usr/bin/python

from math import *

# jk 26.3.2018
# 3.4.2018: plan: add X0prime, X0 = X0prime / rho, similar for Lint_nucl
# add methods to give X0 as well as X0*rho?
# TODO: use for I etc:
# http://pdg.lbl.gov/2017/AtomicNuclearProperties/
# http://pdg.lbl.gov/2017/AtomicNuclearProperties/HTML/silicon_Si.html


### TODO:
# add X0, EC as a property
# need to remember solid/gas! E.g. bool isSolid?
# MakeMaterialFromCompound(chemical)
# and compute mixed averaged Z, A, X0...
# e.g. chemical = 'H2O'


from Consts import *


##########################################
# for radiation length:

def fce(Z):
    a = alpha*Z
    return a*a * (1./(1+a*a) + 0.20206 - 0.0369*a*a + 0.0083*pow(a,4) - 0.002*pow(a,6) )

def Lrad(Z):
    return log(184.15*pow(Z, -1./3.))

def LradPrime(Z):
    return log(1194*pow(Z, -2./3.))

# radiation length in g cm^-2
def GetX0(M, z, Z, A):
    Zfact =  ( Z*Z * ( Lrad(Z) - fce(Z) ) + Z*LradPrime(Z) )
    # ZfactSimple = (Z*Z * log(183./pow(Z,1./3.)) )
    inv = Konst * pow(alpha*hc/M,2) * pow(z,4) / A * Zfact
    #print 'Zfact/Z^2 = {:f}'.format(  ( Z*Z * ( Lrad(Z) - fce(Z) ) + Z*LradPrime(Z) ) / (Z*Z * log(183./pow(Z,1./3.)) ) )
    return 1./inv

##########################################
# more material helper functions:

# TODO: gas vs liquid?!
def GetEC(Z):
    return 610.4 / (Z + 1.24)

# simpler formula for electrons only:.
# radiation length in g cm^-2
def GetX0Ele(A, Z):
    return 716.4 * Z / ( Z*(Z+1)*math.log(287./math.sqrt(Z)) )


def GetHomega(rho, ZoverA):
    return sqrt(rho*ZoverA)*homega0*1.e-6 # MeV;)

##########################################
##########################################
##########################################


class MyMaterial:
    ### I and homega expected in eV, transforming here to MeV!!!
    def Set(self, name, Z, ZoverA, I, rho, homega, S0, S1, a, md, delta0):
        self._name = name
        self._Z = Z
        self._ZoverA = ZoverA
        self._I = I*1.e-6 ### I expected in eV, transforming here to MeV!!!
        self._rho = rho
        self._homega = homega*1.e-6 ### expected in eV, transforming here to MeV!!!
        if self._homega < 0.:
            self._homega = GetHomega(self._rho, self._ZoverA)
        self._S0 = S0
        self._S1 = S1
        self._a = a
        self._md = md
        self._delta0 = delta0

    def __init__(self, element):
        self.Set(element[0], element[1], element[2], element[3], element[4], element[5], element[6], element[7], element[8], element[9], element[10])

    def GetName(self)   : return self._name
    def GetZ(self)      : return self._Z
    def GetZoverA(self) : return self._ZoverA
    def GetA(self)      : return self._Z / self._ZoverA
    def GetI(self)      : return self._I
    def GetRho(self)    : return self._rho
    def GetHomega(self) : return self._homega

    # getters for density effect correction:
    def GetS0(self)     : return self._S0
    def GetS1(self)     : return self._S1
    def Geta(self)      : return self._a
    def GetMd(self)     : return self._md
    def GetDelta0(self) : return self._delta0

        

'''
from RADIATION-INTERACTION-IN-MATTER-AND-DETECTION.pdf
Claude Leroy
Universite de Montreal, Canada
Pier-Giorgio Rancoita
Instituto Nazionale di Fisica Nucleare, Milan, Italy
'''

# OLD notes:
# default parameters for copper: I=332 eV, rho=8.9 g/cm^3
# Bichsel:
# hbar*omega = 31.048 eV is the plasma frequency in Silicon


gElement = [
    # element, Z, Z/A, I eV!, rho g/cm3, plasma energy ev, S0, S1, a, md, delta0
    
['He', 2,  0.500, 41.8, 1.66e-4, 0.26, 2.202, 3.612, 0.134, 5.835, 0.00],
['Li', 3,  0.432, 40.0, 0.53, 13.84, 0.130, 1.640, 0.951, 2.500, 0.14],
['O',  8,  0.500, 95.0, 1.33e-3, 0.74, 1.754, 4.321, 0.118, 3.291, 0.00],
['Ne', 10, 0.496, 137.0, 8.36e-4, 0.59, 2.074, 4.642, 0.081, 3.577, 0.00],
['Al', 13, 0.482, 166.0, 2.70, 32.86, 0.171, 3.013, 0.080, 3.635, 0.12],
['Si', 14, 0.498, 173.0, 2.329, 31.06, 0.201, 2.872, 0.149, 3.255, 0.14],
['Ar', 18, 0.451, 188.0, 1.66e-4, 0.79, 1.764, 4.486, 0.197, 2.962, 0.00],
['Fe', 26, 0.466, 286.0, 7.87, 55.17, -0.001, 3.153, 0.147, 2.963, 0.12],
['Cu', 29, 0.456, 322.0, 8.96, 58.27, -0.025, 3.279, 0.143, 2.904, 0.08],
['Ge', 32, 0.441, 350.0, 5.32, 44.14, 0.338, 3.610, 0.072, 3.331, 0.14],
['Kr', 36, 0.430, 352.0, 3.48e-3, 1.11, 1.716, 5.075, 0.074, 3.405, 0.00],
['Ag', 47, 0.436, 470.0, 10.50, 61.64, 0.066, 3.107, 0.246, 2.690, 0.14],
['Xe', 54, 0.411, 482.0, 5.49e-3, 1.37, 1.563, 4.737, 0.233, 2.741, 0.0],
['Ta', 73, 0.403, 718.0, 16.65, 74.69, 0.212, 3.481, 0.178, 2.762, 0.14],
['W',  74, 0.403, 727.0, 19.30, 80.32, 0.217, 3.496, 0.155, 2.845, 0.14],
['Au', 79, 0.401, 790.0, 19.32, 80.22, 0.202, 3.698, 0.098, 3.110, 0.14],
['Pb', 82, 0.396, 823.0, 11.35, 61.07, 0.378, 3.807, 0.094, 3.161, 0.14],
['U',  92, 0.387, 890.0, 18.95, 77.99, 0.226, 3.372, 0.197, 2.817, 0.14],
]

gCompound = [
# material                          Zdummy, Z/A, I eV!, rho g/cm3, plasma energy ev, S0, S1, a, md, delta0dummy
#(dry) Air at sea level
['Air',                               -1,  0.499, 85.7, 1.21e-3, 0.71, 1.742, 4.276, 0.109, 3.399, -1],
###!!!['Water', ]
['Anthracene',                        -1,  0.527, 69.5, 1.28, 23.70, 0.115, 2.521, 0.147, 3.283, -1],
['Ethane',                            -1,  0.599, 45.4, 1.25e-3, 0.79, 1.511, 3.874, 0.096, 3.610, -1],
['Ethyl Alcohol',                     -1,  0.564, 62.9, 0.79, 19.23, 0.222, 2.705, 0.099, 3.483, -1],
['Freon-12',                          -1,  0.480, 143.0, 1.12, 21.12, 0.304, 3.266, 0.080, 3.463, -1],
['(lead) Glass',                      -1,  0.421, 526.4, 6.22, 46.63, 0.061, 3.815, 0.095, 3.074, -1],
['Kapton polyimide, film',            -1,  0.513, 79.6, 1.42, 24.59, 0.151, 2.563, 0.160, 3.192, -1],
['Lithium carbonate',                 -1,  0.487, 87.9, 2.11, 29.22, 0.055, 2.660, 0.099, 3.542, -1],
['Methane',                           -1,  0.623, 41.7, 6.67e-4, 0.59, 1.626, 3.972, 0.093, 3.626, -1],
['Methanol',                          -1,  0.562, 67.6, 0.79, 19.21, 0.253, 2.764, 0.090, 3.548, -1],
['Plastic scint. vinyltoluene',       -1,  0.541, 64.7, 1.03, 21.54, 0.146, 2.486, 0.161, 3.239, -1],
['Polyethylene',                      -1,  0.570, 57.4, 0.94, 21.10, 0.137, 2.518, 0.121, 3.429, -1],
['Propane',                           -1,  0.590, 47.1, 1.88e-3, 0.96, 1.433, 3.800, 0.099, 3.592, -1],
['Lucite',                            -1,  0.539, 74.0, 1.19, 23.09, 0.182, 2.668, 0.114, 3.384, -1],
['Silicon dioxide',                   -1,  0.499, 139.2, 2.32, 31.01, 0.139, 3.003, 0.084, 3.506, -1],
['Tissue, soft (ICRP)',               -1,  0.551, 72.3, 1.00, 21.39, 0.221, 2.780, 0.089, 3.511, -1],
['Tissue, soft (ICRP, four-comp.)',   -1,  0.550, 74.9, 1.00, 21.37, 0.238, 2.791, 0.096, 3.437, -1],
['Tissue-equiv. gas, (methane base)', -1,  0.550, 61.2, 1.06e-3, 0.70, 1.644, 4.140, 0.099, 3.471, -1],
['Tissue-equiv. gas, (propane base)', -1,  0.550, 59.5, 1.83e-3, 0.91, 1.514, 3.992, 0.098, 3.516, -1],
#    https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=226#:~:text=Mean%20Excitation%20Energy%20(eV)%20%3D,68.700000
['Polystyrene',                       -1,  0.53768, 68.7, 1.06, 68.7, 0, 0, 0, 0, -1]

]

# data are from [Sternheimer, Berger and Seltzer, (1984)]

gMaterials = {}
def MakeMaterials():
    print('Making elements...')
    for element in gElement:
        #print 'Adding ', element
        gMaterials[element[0]] =  MyMaterial(element) 
    for element in gCompound:
        gMaterials[element[0]] =  MyMaterial(element) 

def PrintMaterials():
    for mat in gMaterials:
        print(mat)
