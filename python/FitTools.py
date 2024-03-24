#!/usr/bin/python

# jk 24.8.2021

import ROOT

from math import log, exp, pow

kBadMu = -0.1
kBadP0 = -1
kEpsilon = 1.e-5

gsigVals = []
gbgVals = []
gdataVals = []

# number of studied variables = free parameters
gN = 1

#################################
def getPoissLog(k, nu):
    plog = 0.
    # nu is the expected=true
    # k is the data
    # P(nu,k) = exp(-nu) nu^k / k!
    if nu <= 0.:
        return 0.
    plog = -nu + k*log(nu) # drop the constant data-dependent term (data is given!): -log(factorial(k))
    return plog

#################################
def getNegLogLhood(mu):
    lhoodLog = 0.
    for i in range(0, len(gsigVals)):
        lhoodLog = lhoodLog - getPoissLog(gdataVals[i], mu*gsigVals[i] + gbgVals[i])
    return lhoodLog


#################################
# for Minuit:
# https://root-forum.cern.ch/t/numerical-minimization-with-root-math-functor/13286/2

class MyFunction( ROOT.TPyMultiGenFunction ):
    def __init__( self ):
        print("CREATED")
        ROOT.TPyMultiGenFunction.__init__( self, self )

    def NDim( self ):
        print('PYTHON NDim called')
        return gN

    def DoEval(self, args):
        lhood = 0.
        # TO CHECK!!!
        mu = args[0]
        negloglhood = getNegLogLhood(mu)
        return negloglhood

#################################
def minimizeLhood(step, mumin, mumax):

    mu = 0.5 * (mumax + mumin)
    muerr = kBadMu
    # https://root-forum.cern.ch/t/numerical-minimization-with-root-math-functor/13286#p58684

    # not returning errors:
    # minimizer = ROOT.Math.Factory.CreateMinimizer("GSLMultiMin", "BFGS")
    minimizer = ROOT.Math.Factory.CreateMinimizer('Minuit2', 'Minuit2' )
    minimizer.SetMaxFunctionCalls(1000000)
    minimizer.SetMaxIterations(100000)
    minimizer.SetTolerance(0.001)
    minimizer.SetPrintLevel(1)
    
    myfun = MyFunction()
    # https://root.cern.ch/doc/master/classTMinuitMinimizer.html#a069cdb1a4cf6e2f373832a4cb094c6d2
    minimizer.SetLimitedVariable(0, "mu", mu, step, mumin, mumax)
    minimizer.SetFunction(myfun)

    # minimize the -log(L)!;-)
    minimizer.Minimize()
    status = minimizer.Status()

    # get the minimized parameter:
    if status == 0:
        print('minimizeLhood :: SUCCESSFUL MINIMIZATION! ;-)')
        mu =  minimizer.X()[0]
        muerr = minimizer.Errors()[0]
        #print('minimizeLhood :: Lhood fit result: mu={} +/- {}'.format(mu, muerr))
    else:
        print('minimizeLhood :: FAILED MINIMIZATION! :-(')
        mu = kBadMu
    return mu, muerr


#################################
# support both 1D and 2D versions on demand
def fitSignalStrength(hsig, hbg, hdata, is1D, step = 0.01, mumin = 0., mumax = 2.):

    del gsigVals[:]
    del gbgVals[:]
    del gdataVals[:]
    
    # using lists as storage ;-)
    # 1D:
    if is1D:
        for i in range(1, hdata.GetNbinsX()+1):
            gsigVals.append(hsig.GetBinContent(i))
            gbgVals.append(hbg.GetBinContent(i))
            gdataVals.append(hdata.GetBinContent(i))
    else:
        # 2D:
        for i in range(1, hdata.GetNbinsX()+1):
            for j in range(1, hdata.GetNbinsY()+1):
                gsigVals.append(hsig.GetBinContent(i,j))
                gbgVals.append(hbg.GetBinContent(i,j))
                gdataVals.append(hdata.GetBinContent(i,j))
            
    mu, muerr = minimizeLhood(step, mumin, mumax)
    return mu, muerr

#################################
def ComputeZeroCompatibility(mu, muerr):
    #hmu = dummmodel.Clone('muh')
    #hmu.Reset()
    #hmu.SetBinContent(1, mu)
    #hmu.SetBinError(1, muerr)
    #p0 = hmu.Chi2Test(dummmodel, 'WW') # by default returns the pval
    p0 = kBadP0
    if muerr > 0:
        p0 = ROOT.TMath.Prob(pow(mu/muerr, 2), 1)
    
    t0 = 999.
    logt0 = 999.
    if p0 > 0:
        t0 = -log(p0)
        if t0 > 0.:
            logt0 = log(t0)
    return p0, logt0
            
#################################
#################################
#################################
