#!/usr/bin/python

# jk 24.8.2021
# jk 24.3.2024

import ROOT

from math import log, exp, pow, sqrt
from ctypes import c_double
from array import array

from tofUtil import ms

kBadP0 = -1
kEpsilon = 1.e-5

#gPars = []
gInitPars = {}
gDataPoints = {}

# number of studied variables = free parameters
#gN = -1


# region names:
gTS0p = 'TS0p'
gTS0D = 'TS0D'
gTS0T = 'TS0T'
gTS1p = 'TS1p'
gTS1D = 'TS1D'
gTS1T = 'TS1T'

# not preferred:
gTSp = 'TSboth_p'
gTSD = 'TSboth_D'
gTST = 'TSboth_T'


# https://root.cern/manual/python/#alternative-for-tpymultigenfunction-and-tpymultigradfunction


##################################################################
def initGlobalPars():
    gInitPars['A'] = [1., 500.,]
    gInitPars['B'] = [1., 20.]
    gInitPars['Krel'] = [0.1, 10.]
    gInitPars['Conv'] = [1e-6, 1.e-2]
    return

##################################################################
class cPoint:
    def __init__(self,x,y,ey):
        self.x = x
        self.y = y
        self.ey = ey

##################################################################
# expect a dictionary of dE/dX graphs in fit regions
def prepareData(grs):
    print('Preparing fit data')
    for region in grs:
        gDataPoints[region] = []
    x = c_double(0.)
    y = c_double(0.)

    haveTS0 = False
    haveTS1 = False
    for region,gr in grs.items():
        print(f'  processing region {region}')
        if 'TS0' in region:
            haveTS0 = True
        if 'TS1' in region:
            haveTS1 = True
        for ip in range(0,gr.GetN()):
            gr.GetPoint(ip, x, y)
            ey = gr.GetErrorY(ip)
            gDataPoints[region].append(cPoint(x.value, y.value, 1.*ey))

            
    # basic one region fit, parameters A, B
    npars = 2
    if haveTS0 and haveTS1: # expect we fit split of TS0 and TS1 for one or more particle types
        npars = 4
    print('  initialized {} graph regions'.format(len(grs)))
    return npars

##################################################################
def getBaseFit(beta, A, B, debug = 0):
    fval = 0.
    if beta > 0. and beta < 1.:
        # Bethe-Bloch ansatz:
        fval = A / pow(beta,2) * (log(B*beta/sqrt(1-pow(beta,2))) - pow(beta,2))
    else:
        print('ERROR in getBaseFit: beta out of physically allowed range!')
    #if debug: print(f'fval={fval}')
    return fval

##################################################################
def getFitVal(region, x, npars, pars, debug = 0):
    #A = gPars[0]
    #B = gPars[1]
    A = pars[0]
    B = pars[1]

    beta = 1.*x
    Krel = 1.
    # dEdX congversion par
    C = 0.
    particle = ''

    # names: Trigger Scintillator TS 0 or 1
    # particles p, D, T:
    if npars > 3 and (region == gTS1p or region == gTS1D or region == gTS1T):
        # constant to equalize TS1 to TS0
        Krel = pars[2]
        # conversion from integrated charge in p.e. to MeV
        #print('converting')
        C = pars[3]
        if 'p' == region[-1]:
            particle = 'p'
        elif 'D' == region[-1]:
            particle = 'D'
        elif 'T' == region[-1]:
            particle = 'T'
        if debug: print(f'particle: "{particle}"')
        m = 0.
        try:
            m = ms[particle]
        except:
            print('Failed getting particle mass and correct the beta!')
        if debug: print(f'm={m:1.1f} A={A:1.1f} B={B:1.1f} C={C:1.1f}')
        if m > 0.:
            dE = C * getBaseFit(x, A, B, debug) / 2. # dividing by 2 to account for half material in TS0 compared to TS0+TS1!
            if dE > 0:
                gamma0 = 1./sqrt(1. - pow(beta,2))
                E0 = m*gamma0
                gamma1 = (E0 - dE) / m
                if debug: print(f'x={x} E0={E0:1.1f} dE={dE:1.1f} beta={beta:1.3f} gamma0={gamma0:1.3f} gamma1={gamma1:1.3f}')
                if gamma1 > 1.:
                    beta = sqrt( 1. - 1./pow(gamma1,2))
                else:
                    print('ERROR, negative gamma0!')
                if debug: print(f'  beta0={beta:1.3f}')
            else:
                if debug:
                    print('ERROR, negative energy correction!')
        
    fval = Krel*getBaseFit(beta, A, B, debug)
    return fval

##################################################################
def getChi2Term(region, dpoint, npars, pars, debug):
    term = 0.
    x = dpoint.x
    y = dpoint.y
    ey = dpoint.ey
    if ey <= 0.:
        return 0.
    fitVal = getFitVal(region, x, npars, pars, debug)
    term = pow( (y - fitVal)/ey , 2) 
    return term

##################################################################
#def getChi2(n, grad, x, par):
def getChi2(npars, pars, debug):
    chi2 = 0.
    for region in gDataPoints:
        for dpoint in gDataPoints[region]:
            chi2 = chi2 + getChi2Term(region, dpoint, npars, pars, debug)
    return chi2

##################################################################
# for Minuit:
# https://root-forum.cern.ch/t/numerical-minimization-with-root-math-functor/13286/2
#class MyFunction( ROOT.Math.IMultiGenFunction ):
class MyFunction( object ):
    
    # lass MyFunction( ROOT.Fit.Fitter.MinuitFCN_t ):
    # defint init with npars passed?
    def __init__(self, npars, debug):
        self.npars = npars
        self.debug = debug
    
    #def NDim( self ):
    #    print('PYTHON NDim called')
    #    return self.npars

    #def DoEval(self, x):
    def __call__(self, pars):
        # set the global parameters based on the actual args called
        #for i in range(0, self.NDim()):
        #    gPars[i] = x[i]
        chi2 = getChi2(self.npars, pars, self.debug)
        #if self.debug: print(f'chi2={chi2}')
        return chi2

    #def Clone( self ):
    #    x = MyMultiGenFCN()
    #    ROOT.SetOwnership(x, False)
    #    return x
##################################################################
def mySimpleFcn(npars, x):
    # set the global parameters based on the actual args called
    #for i in range(0, npars):
    #    gPars[i] = x[i]
    chi2 = getChi2(npars, x)
    return chi2


##################################################################
def minimizeChi2(npars, step = 0.01):

    # https://root-forum.cern.ch/t/numerical-minimization-with-root-math-functor/13286#p58684

    # not returning errors:
    # minimizer = ROOT.Math.Factory.CreateMinimizer("GSLMultiMin", "BFGS")
    # OLD:
    """
    minimizer = ROOT.Math.Factory.CreateMinimizer('Minuit2', 'Minuit2' )
    minimizer.SetMaxFunctionCalls(1000000)
    minimizer.SetMaxIterations(100000)
    minimizer.SetTolerance(0.001)
    minimizer.SetPrintLevel(1)
    """
    # OLD: myfun = MyFunction()
    # https://root.cern.ch/doc/master/classTMinuitMinimizer.html#a069cdb1a4cf6e2f373832a4cb094c6d2

    # NEW:
    fitter = ROOT.Fit.Fitter()
    
    pars = []
    # ipar in range(0,npars):
    ipar = -1
    for parName,parLimits in gInitPars.items():
        ipar = ipar + 1
        if ipar >= npars:
            break
        # set initial vals to the mean of limits:
        par = (parLimits[1] - parLimits[0]) / 2.
        print(f'Setting initial par {par}')
        pars.append(par)
        # OLD:
        # minimizer.SetLimitedVariable(0, parName, par, step, parLimits[0], parLimits[1])
        #fitter.SetLimitedVariable(0, parName, par, step, parLimits[0], parLimits[1])

    params = array('d', pars)

    # See:
    # $ROOTSYS/tutorials/fit/combinedFit.py
    #fitter.Config().MinimizerOptions().SetPrintLevel(0)

    fitter.Config().SetMinimizer("Minuit2", "Migrad")
    npars = len(params)
    debug = 1
    myfcn = MyFunction(npars, debug)

    # OLD:
    #minimizer.SetFunction(myfun)
    #minimizer.Minimize()
    #status = minimizer.Status()

    #fitter.FitFCN(myfcn)
    # NEW:
    print(f'Fit parameters: {npars}')
    # we can't pass the Python object globalChi2 directly to FitFCN.
    # It needs to be wrapped in a ROOT::Math::Functor.
    globalChi2Functor = ROOT.Math.Functor(myfcn, npars)
    #fitter.FitFCN(globalChi2Functor)
    dataSize = len(gDataPoints)

    print('Setting parameters limits')
    fitter.Config().SetParamsSettings(npars, params)
    ipar = -1
    for parName,parLimits in gInitPars.items():
        ipar = ipar + 1
        if ipar >= npars:
            break
        fitter.Config().ParSettings(ipar).SetLimits(parLimits[0], parLimits[1])

    fitter.Config().ParSettings(0).SetName("A")
    fitter.Config().ParSettings(1).SetName("B")
    if npars > 3:
        fitter.Config().ParSettings(2).SetName("Krel")
        fitter.Config().ParSettings(3).SetName("Conv")
    
    fitter.FitFCN(globalChi2Functor) #, 0, dataSize, True)

    # TO TRY:
    # fitter.Config().ParSettings(0).SetName("Parameter0").SetLimits(-10, 10)

    # back to minimizer...or not
    #minimizer = ROOT.Math.Minimizer(ROOT.Math.MinimizerOptions.DefaultMinimizerType())
    #minimizer.SetMaxFunctionCalls(10000)
    #minimizer.SetMaxIterations(1000)
    #minimizer.SetTolerance(0.001)
    #minimizer.SetFunction(myfcn)
    #minimizer.Minimize()
    
    #fitter.FitFCN(mySimpleFcn, len(params), params)
    #status = fitter.Result()
    # status = fitter.Result()
    #.Print(ROOT.std.cout, True)

 
    result = fitter.Result()
    result.Print(ROOT.std.cout)
    bestPars = result.Parameters()
    parErrs = result.Errors()
    # https://root.cern/doc/v610/classROOT_1_1Fit_1_1FitResult.html
    print('*** FIT RESULTS ***')
    chi2 = result.MinFcnValue() # .Chi2()
    ndf = 0 # #result.Ndf()
    for region,dpoints in gDataPoints.items():
        ndf = ndf + len(dpoints)
    ndf = ndf - npars
    pval = result.Prob()
    status = result.Status()
    print('Fit status {}, chi2/ndf={:1.3f}/{:}'.format(status, chi2, ndf))
    if ndf > 0:
        print('              chi2/ndf={:1.3f}'.format(chi2/ndf))
    for ipar in range(0, result.NPar()):
        pname = result.ParName(ipar)
        print('Parameter {:} {:1.3f} +/- {:1.3f}'.format(pname, bestPars[ipar], parErrs[ipar]))
    
    pars = []
    parerrs = []
    ## get the minimized parameter:
    #if status == 0:
    #    print('minimizeChi2 :: SUCCESSFUL MINIMIZATION! ;-)')
    #    ipar = -1
    #    Xs =  minimizer.X()
    #    Errs = minimizer.Errors()
    #    for parname in gInitPars:
    #        ipar = ipar+1
    #        par =  Xs[ipar]
    #        parerr = Errs[ipar]
    #        pars.append(par)
    #        parerrs.append(parerr)
    #        print('minimizeChi2 fit result: par {}: {} +/- {}'.format(parname, par, parerr))
    #else:
    #    print('minimizeLhood :: FAILED MINIMIZATION! :-(')

    return bestPars, parErrs


##################################################################
# support both 1D and 2D versions on demand
def doTheFit(grs, step = 0.001):

    gInitPars = {}
    #del gPars[:]

    initGlobalPars()
    npars = prepareData(grs)
    print('Will fit data points:')
    for region in gDataPoints:
        print(f'  Region {region}')
        for dpoint in gDataPoints[region]:
            print('    beta={} y={} ey={}'.format(dpoint.x, dpoint.y, dpoint.ey))
    pars, parerrs = minimizeChi2(npars, step)
    return pars, parerrs

##################################################################
def ComputeZeroCompatibility(val, sigma):
    p0 = kBadP0
    if sigma > 0:
        p0 = ROOT.TMath.Prob(pow(val/sigma, 2), 1)
    t0 = 999.
    logt0 = 999.
    if p0 > 0:
        t0 = -log(p0)
        if t0 > 0.:
            logt0 = log(t0)
    return p0, logt0
            
##################################################################
##################################################################
##################################################################
