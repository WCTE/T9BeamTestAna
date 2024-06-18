#!/usr/bin/python

# jk 24.8.2021
# jk 24.3.2024

import ROOT

from math import log, exp, pow, sqrt
from ctypes import c_double
from array import array

from graphTools import *
from tofUtil import ms
from Losses import *
#from Brems import *

# Define color codes
RED = "\033[31m"
LRED = "\033[91m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
CYAN = "\033[36m"
RESET = "\033[0m"
BOLD = "\033[1m"
LBLUE = "\033[94m"
   
gFitMeanLosses = True

kBadP0 = -1
kEpsilon = 1.e-5
gme = 0.511e6 # eV!
gScintiThicknes = 0.6 # cm!
gj = 0.200 # Landau most probable losses constant

#gPars = []
gInitPars = {}
gDataPoints = {}

# number of studied variables = free parameters
# doesn't work globally...
# gN = -1

# https://root.cern/manual/python/#alternative-for-tpymultigenfunction-and-tpymultigradfunction


##################################################################
def printMatrix(corr, npars):
  for i in range(0,npars):
    for j in range(0,npars):
        val = corr[i][j]
        if i == j:
            print(f'{BOLD}{GREEN}', end='')
        else:
            if val > 0:
                print(f'{LRED}', end='')
            else:
                print(f'{CYAN}', end='')
        print(f' {val:+1.5f}', end='')
        print(f'{RESET}', end='')
    print('')
  return

##################################################################
def printMatrixToFile(outfile, corr, npars, parNames):
  ll = ''
  for i in range(0,npars):
    ll = ll + 'l'
  outfile.write('\nCorrelation matrix:\n\n' + r'\begin{tabular}{l|' + ll + '} ' + '\n')
  outfile.write(' & ')
  for pn in parNames:
    outfile.write(pn)
    if pn != parNames[-1]:
      outfile.write(' & ')
  outfile.write(r' \\ \hline' + '\n')
    
  for i in range(0,npars):
    outfile.write(f' {parNames[i]} & ')
    for j in range(0,npars):
        val = corr[i][j]
        outfile.write(f'${val:+1.4f}$')
        if j < npars-1:
          outfile.write(' & ')
    outfile.write(r' \\ ' + '\n')
  outfile.write(r'\end{tabular}' + '\n')
  return


##################################################################
def initGlobalPars():
    gInitPars['A'] = [1., 50.,]
    gInitPars['I'] = [20., 10e3]
    gInitPars['C'] = [-50., 50.]
    #gInitPars['dXScale'] = [1., 2.]
    gInitPars['Krel'] = [1., 10.]
    gInitPars['Conv'] = [1e-6, 1. ] #1.e-1]
    return

##################################################################
class cPoint:
    def __init__(self,x,y,ey):
        self.x = x
        self.y = y
        self.ey = ey
        # fitted value to be computed after the fit
        # to be used for fit residuals later
        self.yfit = 0.
        # fitted and theory deltaE after passing PMTs of T0X, X = 0,1,2,3
        # in MeV
        self.dEfit = 0.
        self.dEtheory = 0.

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

            
    # basic one region fit, parameters A, I, C
    npars = 3
    if haveTS0 and haveTS1: # expect we fit split of TS0 and TS1 for one or more particle types
        npars = 5
    print('  initialized {} graph regions'.format(len(grs)))
    return npars

##################################################################
def getBaseMeanLossesFitVal(beta, A, I, C, debug = 0):
    fval = 0.
    if beta > 0. and beta < 1.:
        # Bethe-Bloch ansatz:
        fval = A / pow(beta,2) * ( log(2*gme/I*beta*beta/(1-pow(beta,2))) - pow(beta,2) ) + C
    else:
        print('ERROR in getBaseMeanLossesFitVal: beta out of physically allowed range!')
    #if debug: print(f'fval={fval}')
    return fval

##################################################################
def getBasePeakLossesFitVal(beta, A, I, C, dX, debug = 0):
    fval = 0.
    if beta > 0. and beta < 1.:
        # Bethe-Bloch ansatz:
        xi = A / 2. / pow(beta,2)
        fval = xi  * ( log(2*gme/I*beta*beta/(1-pow(beta,2))) + log(xi/I) + gj - pow(beta,2) ) + C
    else:
        print('ERROR in getBasePeakLossesFitVal: beta out of physically allowed range!')
    #if debug: print(f'fval={fval}')
    return fval

##################################################################
def getFitVal(region, x, npars, pars, debug = 0):
    A = pars[0]
    I = pars[1]
    C = pars[2]
    dXScale = 1. # cannot be a fit parametr as is completely correlated to Conv! # was: = pars[3]

    beta = 1.*x
    Krel = 1.
    # dEdX conversion par, for the moment zero, will be nonzero when needed for regions TS1X
    Conv = 0.
    particle = ''
    dE, theorydE = 0., 0.

    # names: Trigger Scintillator TS 0 or 1
    # particles p, D, T:
    if npars > 3 and 'TS1' in region:
        # constant to shift charges between TS1 to TS0
        Krel = pars[3]
        # conversion from integrated charge in p.e. to MeV using the fitted dE/dx
        # print('converting')
        Conv = pars[4]
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
        if debug: print(f'  m={m:1.1f} A={A:1.1f} I={I:1.1f} Krel={Krel:1.3f} C={C:1.4f} Conv={Conv:1.4f}')
        if m > 0.:
            dX = dXScale*gScintiThicknes
            # evaluate the fitted lossed in TS0
            dE = dX * Conv * getBaseMeanLossesFitVal(beta, A, I, C, debug) # ??? was NONSENSE???: / 2. # dividing by 2 to account for half material in TS0 compared to TS0+TS1!
            if dE > 0 and beta > 0. and beta < 1.:
                gamma0 = 1./sqrt(1. - pow(beta,2))
                E0 = m*gamma0
                p0 = sqrt(E0*E0 - m*m)
                gamma1 = (E0 - dE) / m
                fitFracEloss = dE / E0
                
                useHigherCorrs = False
                material = gMaterials['Polystyrene']
                #print(beta, particle, material, useHigherCorrs)
                fullName = 'e'
                if particle == 'p':
                    fullName = 'Proton'
                if particle == 'D':
                    fullName = 'Deuteron'
                if particle == 'T':
                    fullName = 'Tritium'

                # theory losses:
                theorydEdX, halflog = dEdX(beta, gParticles[fullName], material, useHigherCorrs)
                
                theorydE = theorydEdX*gScintiThicknes*dXScale
                #print(theorydE)
                newE = E0 - theorydE
                #newp = sqrt(newE*newE - m*m)
                #newT = newE - m
                theoryFracEloss = theorydE / E0
                #newbeta = newp / newE
                newgamma = newE / m
                newbeta = sqrt( 1. - 1./pow(newgamma,2))
                
                if debug:
                    print(f'  beta0={x:1.4} p0={p0:1.1f} E0={E0:1.1f} MeV; using dX={dX:1.2f}cm theorydE={theorydE:1.2f} MeV, actual dE={dE:1.2f} MeV, dE/theorydE ratio: {dE/theorydE:1.3f}; beta={beta:1.3f} gamma0={gamma0:1.3f} gamma1={gamma1:1.3f}')
                if gamma1 > 1.:
                    # HERE WE CHANGE THE BETA to adjust for energy losses after passing through T0!
                    # we use the new energy after fitted losses, which then define the new gamma factor at T1, therefore called gamma1
                    beta = sqrt( 1. - 1./pow(gamma1,2))
                else:
                    print('ERROR, negative gamma0!')
                if debug:
                    print(f'  beta1={beta:1.4f} theory dbeta = {x - newbeta:1.4f}, actual dbeta = {x-beta:1.4f}')
                    print(f'  fracEloss theory: {theoryFracEloss:1.5f}, actually fitted: {fitFracEloss:1.5f}')
            else:
                if debug:
                    print('ERROR, negative energy correction!')

    # EVALUATE THE ENERGY LOSSES, possibly with modified beta
    fval = 0.
    if gFitMeanLosses:
        fval = Krel*getBaseMeanLossesFitVal(beta, A, I, C, debug)
    else:
      # this is not ready yet and not so simple to translate from peak to mean losses!!
      # need a new function for this!
      #fval = Krel*getBasePeakLossesFitVal(beta, A, I, C, dX, debug)
      pass
    return fval, beta, dE, theorydE

##################################################################
def getChi2Term(region, dpoint, npars, pars, debug):
    term = 0.
    x = dpoint.x
    y = dpoint.y
    ey = dpoint.ey
    if ey <= 0.:
        return 0.
    fitVal, beta, dE, theorydE = getFitVal(region, x, npars, pars, debug)
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
def minimizeChi2(npars, nCalibCs, step = 0.01, debug = 0):

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
    parNames = []
    for parName,parLimits in gInitPars.items():
        parNames.append(parName)
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
    fitter.Config().ParSettings(1).SetName("I")
    fitter.Config().ParSettings(2).SetName("C")
    #fitter.Config().ParSettings(3).SetName("dxScale")
    if npars > 3:
        fitter.Config().ParSettings(3).SetName("Krel")
        fitter.Config().ParSettings(4).SetName("Conv")
    
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

    return fitter, result, npars, parNames

##################################################################
# June 2024

def AnalyzeFitResults(GrsBeta, fitter, result, nCalibCs, parNames, extraTag):

    npars = len(parNames)
    bestPars = result.Parameters()
    parErrs = result.Errors()  
    
    corr = ROOT.TMatrix(npars,npars)
    cov = ROOT.TMatrix(npars,npars)
    # https://root.cern/doc/v606/classROOT_1_1Fit_1_1FitResult.html
    result.GetCorrelationMatrix(corr)
    result.GetCovarianceMatrix(cov)
    result.Print(ROOT.std.cout)
    print('--- Correlation matrix:')
    printMatrix(corr, npars)
    print('--- Covariance matrix:')
    printMatrix(cov, npars)

    ### print fitted parameters and correlation matrix, also to a TeX file
    
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
    outtex = open('tex/fitpars_eRelCalib{}{}.tex'.format(nCalibCs, extraTag), 'w')
    outtex.write(r'\begin{tabular}{l|ll}' + '\n')
    outtex.write(r'Parameter & value & uncertainty \\ \hline' + '\n')
    for ipar in range(0, result.NPar()):
        pname = result.ParName(ipar)
        print('Parameter {:} {:1.3f} +/- {:1.3f}'.format(pname, bestPars[ipar], parErrs[ipar]))
        outtex.write(' {:} & {:1.3f} & {:1.3f} \\\\ \n'.format(pname, bestPars[ipar], parErrs[ipar]))
    if ndf > 0:
        outtex.write(r' \hline' + '\n')
        outtex.write(r' $\chi^2/\mathrm{ndf}$ & ' + '{:1.3f}/{:}'.format(chi2, ndf) + r' & \\' + '\n')
        outtex.write(r' $\chi^2/\mathrm{ndf}$ & ' + '{:1.3f}'.format(chi2/ndf) + r' & \\' + '\n')
    outtex.write(r'\end{tabular}' + '\n')
    printMatrixToFile(outtex, corr, npars, parNames)
    outtex.close()

    ### compute fit residuals:
    rmax = 50.
    hname = 'FitResiduals'
    htitle = ';fit residuals [N_{p.e.}];data points'
    nb = 100
    hres = ROOT.TH1D(hname, htitle, nb, -rmax, rmax)
    hres.SetLineWidth(2)
    hres.SetLineColor(ROOT.kRed)

    ### fill histograms of interest
    
    hname = 'dEratio'
    htitle = ';#DeltaE^{fit}/#DeltaE^{theory};data points'
    hdEFitOverTheory = ROOT.TH1D(hname, htitle, 50, 0., 2.)
    hdEFitOverTheory.SetLineWidth(2)
    hdEFitOverTheory.SetLineColor(ROOT.kBlue)
    
    canname = 'FitAnalysis{}'.format(nCalibCs)
    canres = ROOT.TCanvas(canname, canname, 0, 0, 1200, 600)
    canres.Divide(2,1)
    
    for region in gDataPoints:
      for dpoint in gDataPoints[region]:

        beta = dpoint.x
        debug = 1
        fitVal, beta, dE, theorydE = getFitVal(region, beta, len(bestPars), bestPars, 0)
        #print(fitVal, beta, dE, theorydE)
        #A, I, C, Krel, Conv = bestPars[0], bestPars[1], bestPars[2], bestPars[3], bestPars[4]
        dpoint.yfit = fitVal
        dpoint.dEfit = dE
        dpoint.dEtheory = theorydE
        if  theorydE  > 0.:
          dEfitOverTheory = dE/theorydE
          hdEFitOverTheory.Fill(dEfitOverTheory)
        residual = dpoint.y - dpoint.yfit
        hres.Fill(residual)

        
        # hard to compute total fit unc. b.c of the formula and all correlations...
        # dpoint.yfiterr = uff...;-)
        # totalErr = sqrt( pow(dpoint.yerr,2) + dpoint.yfiterr)
        # pull = residual / totalErr 
    canres.cd(1)
    hres.Draw()
    canres.cd(2)
    hdEFitOverTheory.Draw()
    canres.Update()

    ### Plot all data points and connect fit points in one canvas
    GrsFit = {}
    for region in gDataPoints:
      xs = []
      ys = []
      for dpoint in gDataPoints[region]:
        beta = dpoint.x
        fitVal, newbeta, dE, theorydE = getFitVal(region, beta, len(bestPars), bestPars, 0)
        xs.append(beta)
        ys.append(fitVal)
      GrsFit[region] = MakeGraphNoErrs(xs, ys, ROOT.kRed, 20, 0.)
      GrsFit[region].SetLineStyle(1)
      GrsFit[region].SetLineWidth(1)

    scint1 = 60.
    scint2 = 240.
    hn = 'hbtmp_data_fit'
    bmin = 0.4
    bmax = 0.8
    hb = ROOT.TH2D(hn, ';#beta;Rel. calib charge [N_{p.e.}]', 100, bmin, bmax, 100, scint1, scint2)
    hb.SetStats(0)
    
    canname = 'DataFitCmp{}'.format(nCalibCs)
    cancmp = ROOT.TCanvas(canname, canname, 200, 200, 1000, 900)
    leg = ROOT.TLegend(0.5, 0.5, 0.88, 0.88)
    leg.SetNColumns(2)
    cancmp.cd()
    hb.Draw()
    
    mstsf = [20, 21, 22, 23, 29, 33, 34, 45]
    mstso = [24, 25, 26, 32, 30, 27, 28, 44]
    ireg = -1

    chi2s = {}
    npts = {}
    for region in GrsFit:
      ireg = ireg+1
      gr = GrsBeta[region]
      grfit = GrsFit[region]
      col = ROOT.kBlack
      if region[-1] == 'D':
        col = ROOT.kBlue
      gr.SetMarkerColor(col)
      gr.SetMarkerSize(1.5)
      gr.SetLineColor(col)
      grfit.SetLineColor(col)
      mst = mstsf[ireg % 8]
      lst = 1
      if region[2] == '1':
        mst = mstso[ireg % 8]
        lst = 2
      gr.SetMarkerStyle(mst)
      grfit.SetLineStyle(lst)
      chi2s[region], npts[region] = getGraphsChi2(gr, grfit)
      chi2npts = 0.
      if npts[region] > 0:
        chi2npts =  chi2s[region] / npts[region] 
      
      leg.AddEntry(gr, region + ' data', 'P')
      leg.AddEntry(grfit, 'fit #chi^{2}/N_{pts}=' + f'{chi2npts:2.1f}', 'L')
      gr.Draw('P')
      grfit.Draw('L')

    leg.Draw()
    
    #pars = []
    #parerrs = []
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
    return canres, cancmp, leg, hres, hdEFitOverTheory, hb, GrsFit
##################################################################
# support both 1D and 2D versions on demand
def doTheFit(grs, nCalibCs, step = 0.001, debug = 0):

    gInitPars = {}
    #del gPars[:]

    initGlobalPars()
    npars = prepareData(grs)
    print('Will fit data points:')
    for region in gDataPoints:
        print(f'  Region {region}')
        for dpoint in gDataPoints[region]:
            print('    beta={} y={} ey={}'.format(dpoint.x, dpoint.y, dpoint.ey))
    fitter, result, npars, parNames = minimizeChi2(npars, nCalibCs, step, debug)
    return fitter, result, npars, parNames

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
