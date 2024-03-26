#!/usr/bin/python

# jk 26.4.2018

from Losses import *
from Brems import *

print('*** Supported particles:')
PrintParticles()
print('*** Supported materials:')
PrintMaterials()


print('***************************************************')

pairs = [ #['Alpha', 'Si'],
          #['Deuteron', 'Si'],
          #['Proton', 'Si'],
          #['Kaon', 'Si'],
          #['Pion', 'Si'],
          #['Muon', 'Si'],
          #['Electron', 'Si'],
          #['Positron', 'Si'],
    
    ['Electron', 'Polystyrene'],
    ['Positron', 'Polystyrene'],
    ['Muon', 'Polystyrene'],
    ['Pion', 'Polystyrene'],
    ['Proton', 'Polystyrene'],
    ['Deuteron', 'Polystyrene'],

]

#KeepSameMomentum = False
KeepSameMomentum = True

for pair in pairs:
    pname = pair[0]
    mname = pair[1]
    particle = gParticles[pname]
    material = gMaterials[mname]

    beta = 0.
    gamma = 0.
    bg = 0.
    p = 0.
    M = particle.GetM()
    T = 0.

    useHigherCorrs = False
    SomeEnergy  = 540. # MeV
    
    if KeepSameMomentum:
        p = 1.*SomeEnergy # MeV
        E = sqrt(p*p + M*M)
        T = E - M
        # beta*gamma:
        bg = p/M
        beta = GetBetaFromBg(bg)
        gamma = GetGammaFromBg(bg)
    else:
        T = 1.*SomeEnergy
        gamma = 1. + T/M
        beta = sqrt(1. - 1./pow(gamma,2))
        gamma = GetGamma(beta)
        bg = beta*gamma
        p = bg*M

    pname = particle.GetName()
    print('*** {:} in {:} '.format(pname, material.GetName()))
    # print bg, beta, gamma
    print('    BEFORE:    T={:3.1f} MeV p={:3.1f} MeV, E={:3.1f}, beta={:1.4f}, gamma = {:3.3f}, beta*gamma={:3.3f}'.format(T, p, E, beta, gamma, beta*gamma))
    #print(beta, particle, material, useHigherCorrs)
    dedx, halflog = dEdX(beta, particle, material, useHigherCorrs)
    newE = E - dedx
    newp = sqrt(newE*newE - M*M)
    newT = newE - M
    fracEloss = dedx / E
    newbeta = newp / newE
    newgamma = newE/M
    print('    AFTER 1cm: T={:3.1f} MeV p={:3.1f} MeV, E={:3.1f}, beta={:1.4f}, gamma = {:3.3f}, beta*gamma={:3.3f}'.format(newT, newp, newE, newbeta, newgamma, newbeta*newgamma))
    print('    Ionization losses                                      : {:1.3f} MeV/cm'.format(dedx))
    #print('    New momentum after 1cm of the material                 : {:3.1f}'.format(newp))
    #print('    New energy after 1cm of the material                   : {:3.1f}'.format(newE))
    #print('    New kinetic energy after 1cm of the material           : {:3.1f}'.format(newT))

    print('    Energy fraction loss after 1cm of the material         : {:1.4f}'.format(dedx / E))
    print('    Kinetic energy fraction loss after 1cm of the material : {:1.4f}'.format(dedx / T))

    print('    New momentum fraction after 1cm of the material        : {:1.3f}'.format(newp/p))
    print('    New energy fraction after 1cm of the material          : {:1.3f}'.format(newE/E))
    print('    New kinetic energy fraction after 1cm of the material  : {:1.3f}'.format(newT/T))
    print('    New beta, deltaBeta, new gamma                         : {:1.3f}, {:1.3f}, {:1.3f}'.format(newbeta, beta-newbeta, newgamma))

    
    #E = gamma*M
    #print('    Radiation losses  : {:1.3f} MeV/cm'.format( dEdXBrems(E, particle, material), ))
