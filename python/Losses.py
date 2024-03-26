#!/usr/bin/python

# Jiri Kvita 2015-2018
#############################################################


from MaterialConsts import *
from Particles import *

if len(gMaterials) == 0:
    MakeMaterials()
if len(gParticles) == 0:
    MakeParticles()



#############################################################
# density effect:
def GetDelta(beta, gamma, material):
    delta = 0.
    
    I = material.GetI()
    homega = material.GetHomega()

    bg = beta*gamma
    bglog = log10(bg)
    bgln = log(bg)
    C = 2.*log(homega/I) - 1.
    
    if bglog  > material.GetS1():
        #print '    ... bg > S1'
        delta = 2.*bgln  + C
    elif bglog  > material.GetS0():
        #print '    ... S1 > bg > S0'
        delta = C + 2.*bgln + material.Geta()*pow( 1./log(10.)*log( pow(10,material.GetS1())/ bg),   material.GetMd()  )
    else:
        #print '    ... bg < S0'
        delta = material.GetDelta0() * pow(bg / pow(10, material.GetS0()), 2) 
    #print '      delta=', delta
    return delta

#############################################################
# ionization potential in MeV!
def GetI(Z):
    return pow(Z,0.9)*1e-6

#############################################################

def GetGamma(beta):
    gamma = 0
    if beta > 0 and beta < 1:
        gamma = 1./sqrt(1-beta*beta)
    return gamma

def GetBeta(gamma):
    beta = 0
    if gamma > 1.:
        beta = sqrt(1.-1./(gamma*gamma))
    return beta

def GetGammaFromBg(bg):
    gamma = 0
    if bg > 0:
        gamma = sqrt(1. + bg*bg)
    return gamma

def GetBetaFromBg(bg):
    beta = 0.
    if bg > 0:
        beta = bg /sqrt(1. + bg*bg)
    return beta

def GetTmax(beta, M):
    gamma = GetGamma(beta)
    tmax = 2.*me*pow(beta*gamma,2) / (1. + 2.*gamma*me/M + pow(me/M,2))
    return tmax

#############################################################

def dEdX(beta, particle, material, useHigherCorrs = True):
    #print('in dEdX')
    #print(beta, particle, material, useHigherCorrs)
    if particle.GetName() == 'Electron':
        #print('Calling Moller')
        return dEdXMoller(beta, particle, material)
    elif particle.GetName() == 'Positron':
        return dEdXBhabha(beta, particle, material)
    else:
        return dEdXHeavy(beta, particle, material, useHigherCorrs)

#############################################################

def dEdXHeavy(beta, particle, material, useHigherCorrs):
    if particle.GetName() == 'Electron' or particle.GetName() == 'Positron':
        print('WARNING, using Bethe-Bloch ionization losses, intended for heavy particles, for Electron or Positron, which is not applicable!')
    rho = material.GetRho()
    I = material.GetI()
    homega = material.GetHomega()

    #print 'rho, Z/A, I, homega:'
    #print rho, material.GetZoverA(), I, homega
    
    # protect against low beta, validity of Bethe-Bloch formula is only
    # down to beta of 0.1
    interpolateToZero = False
    betaorig = beta
    if beta < betaMin:
        interpolateToZero = True
        beta = betaMin
    gamma = GetGamma(beta)
    gamma2 = gamma*gamma
    beta2 = beta*beta
    Tmax = GetTmax(beta,particle.GetM())
    delta = GetDelta(beta, gamma, material)
    halflog = 0.5*log(2.*me*beta2*gamma2*Tmax/(I*I))
    higherCorrs = 0.
    if useHigherCorrs:
        higherCorrs =  - delta/2.
    dedx = rho*K*pow(particle.GetZ(),2)*material.GetZoverA()*1./beta2 * ( halflog - beta2 )

    if interpolateToZero:
        # linear:
        dedx = betaorig*dedx/betaMin
        # exponential:
        #alpha = 1.
        #dedx = dedx*(exp(alpha*betaorig) - 1.) / (exp(alpha*betaMin) -1)
    #print '   beta, gamma       : %f, %f' % (beta,gamma,)
    #print '   Argument logaritmu: %f' % (2*me*beta2*gamma2*Tmax/(I*I))
    #print '   Tmax/I            : %f' % (Tmax/I,)
    #print '   Tmax       [MeV]  : %f' % (Tmax,)
    #print '   Logaritmus        : %f' % (log(2*me*beta2*gamma2*Tmax/(I*I)),)
    #print '   Logaritmus        : %f' % (log(Tmax*Tmax/(I*I)),)
    #print '   delta             :', delta
    #print '   dedx              : ', dedx

    # returns energy loss in MeV/cm:
    return dedx, halflog


#############################################################
# losses for electrons!
def dEdXMoller(beta, particle, material):
    #print('In Moller')
    if particle.GetName() != 'Electron':
        print('WARNING, using Moller ionization losses intended for Electron for particle {:} instead!'.format(particle.GetName()))
    rho = material.GetRho()
    I = material.GetI()
    homega = material.GetHomega()
    #print('In Moller 2')
    # protect against low beta, validity of Bethe-Bloch formula is only
    # down to beta of 0.1
    interpolateToZero = False
    betaorig = beta
    if beta < betaMin:
        interpolateToZero = True
        beta = betaMin
    gamma = GetGamma(beta)
    gamma2 = gamma*gamma
    beta2 = beta*beta
    Tmax = GetTmax(beta,particle.GetM())
    delta = GetDelta(beta, gamma, material)
    halfLog =  log(me*beta2*gamma2*0.5*me*(gamma-1)/(I*I))
    dedx = rho*K/2.*material.GetZoverA()*1./beta2 * ( halfLog + 1 - beta2 - (2*gamma-1)/gamma2*log(2) + 1./8.*pow((gamma-1)/gamma, 2) - delta )
    #print('In Moller 3')
    if interpolateToZero:
        dedx = betaorig*dedx/betaMin

    # returns energy loss in MeV/cm:
    return dedx, halfLog



#############################################################
# losses for positrons!
def dEdXBhabha(beta, particle, material):
    if particle.GetName() != 'Positron':
        print('WARNING, using Bhabha ionization losses intended for Positron for particle {:} instead!'.format(particle.GetName()))
    rho = material.GetRho()
    I = material.GetI()
    homega = material.GetHomega()

    # protect against low beta, validity of Bethe-Bloch formula is only
    # down to beta of 0.1
    interpolateToZero = False
    betaorig = beta
    if beta < betaMin:
        interpolateToZero = True
        beta = betaMin
    gamma = GetGamma(beta)
    gamma2 = gamma*gamma
    beta2 = beta*beta
    Tmax = GetTmax(beta,particle.GetM())
    delta = GetDelta(beta, gamma, material)
    halfLog =  log(me*beta2*gamma2*0.5*me*(gamma-1)/(I*I))
    dedx = rho*K/2.*material.GetZoverA()*1./beta2 * (halfLog + 2.*log(2.) - beta2/12.*(23 + 14/(gamma+1) + 10/pow(gamma+1, 2) + 4/pow(gamma+1, 3) ) - delta )

    if interpolateToZero:
        dedx = betaorig*dedx/betaMin

    # returns energy loss in MeV/cm:
    return dedx, halfLog



#############################################################

def GetMostProbableLoss(x, beta, particle, material, useHigherCorrs = True):
    if particle.GetName() == 'Electron' or particle.GetName() == 'Positron':
        print('WARNING, using most probable ionization losses, intended for heavy particles, for Electron or Positron, which is not applicable!')
    rho = material.GetRho()
    I = material.GetI()
    homega = material.GetHomega()

    # protect against low beta, validity of Bethe-Bloch formula is only
    # down to beta of 0.1
    #interpolateToZero = False
    betaorig = beta
    if beta < betaMin:
        #interpolateToZero = True
        beta = betaMin
    gamma = GetGamma(beta)
    gamma2 = gamma*gamma
    beta2 = beta*beta
    xi = x/2.*rho*K*pow(particle.GetZ(),2)*material.GetZoverA()/beta2 # MeV
    delta = GetDelta(beta, gamma, material)
    higherCorrs = 0.
    if useHigherCorrs:
        higherCorrs = gkjterm + log(xi/I) - delta
    de =  xi * ( log(2*me*beta2*gamma2/I) - beta2 + higherCorrs ) # MeV
    #print 'x=', x*1e4, 'mum; xi=', xi, ' xi/I=', xi/I, ' delta=', delta, ' log(2*me*beta2*gamma2/I)=', log(2*me*beta2*gamma2/I)
    return de


#############################################################

#
# Obsolete:
#############################################################
# Landau; Claude Leroy, Pier-Georgio Rancoita
# default parameters for copper: I=332 eV, rho=8.9 g/cm^3
def GetMostProbableLossAlt(x, beta, particle, material, UseRelativisticDelta, useHigherCorrs = True):
    print('You should not be using this function!!!')
    rho = material.GetRho()
    A = material.GetA()
    Z = material.GetZ()
    I = material.GetI()
    homega = material.GetHomega()

    # protect against low beta, validity of Bethe-Bloch formula is only
    # down to beta of 0.1
    #interpolateToZero = False
    betaorig = beta
    if beta < betaMin:
        #interpolateToZero = True
        beta = betaMin
    gamma = GetGamma(beta)
    gamma2 = gamma*gamma
    beta2 = beta*beta
    xi = 0.1535*x*rho*pow(particle.GetZ(),2)*Z/A/beta2 # MeV
    delta = GetDelta(beta, gamma, material)
    dedxmean = dEdX(beta, particle, material, useHigherCorrs)
    Tmax = GetTmax(beta,particle.GetM())
    ded = 0.
    #if not UseRelativisticDelta:
    #    print 'OK, using NON-relativistic formula!'
    #print 'dedx   =', dedxmean
    #print 'de     =', dedxmean*x
    #print 'Landau =',  xi * ( beta2 + log(xi/Tmax) + 0.194 ) # MeV
    ded =  dedxmean*x + xi * ( beta2 + log(xi/Tmax) + 0.194 ) # MeV
    #else:
    #print 'OK, using relativistic formula!'
    ded = xi* ( log(2.*me*xi / (homega*homega)) + 0.194 )
    #print 'x=', x*1e4, 'mum; xi=', xi, ' xi/I=', xi/I, ' delta=', delta, ' log(2*me*beta2*gamma2/I)=', log(2*me*beta2*gamma2/I)
    return ded
    


#############################################################
