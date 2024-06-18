#!/snap/bin/pyroot

###  /usr/bin/python3

# jk
# 20/09/2022
# 14.7.2023
#  5.3.2024
#  6.3.2024
# 24.3.2024

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from labelTools import *
from tofUtil import *
from FitTools import *
#from graphTools import *


stuff = []
lines = []

####################################################################################

def PrintUsage(argv):
    print('Usage:')
    print('{} filename_plots.root [-b]'.format(argv[0]))
    print('Example:')
    print('{} output_300n_plots.root -b'.format(argv[0]))
    return

####################################################################################
def readInputFiles():
    dirname = 'histos/windowpe_analyzed/'
    filenames = [
        'peakAnalysed_timeCorr_windInt_000403_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000396_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000394_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000393_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000392_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000398_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000399_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000402_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000449_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000436_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000435_plots_f.root'
    ]
    rfiles = []
    for filename in filenames:
        rfile = ROOT.TFile(dirname + filename, 'read')
        if not rfile.IsZombie():
            rfiles.append(rfile)
    return rfiles

####################################################################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def singleFit(argv, rfiles, TStag, particle, calibOnly, calibCs, minEntries = 100., Opt1d = 'hist', Opt2d = 'box', drawElAnyway = False):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]
    cans = []
    reltag = ''
    if len(calibCs) > 0:
        print('calibCs:')
        print(calibCs)
        refmomentum = 1000
        reltag = '_eRelCalib{}'.format(len(calibCs[refmomentum]))

    ROOT.gStyle.SetOptFit(111)
    #ROOT.gStyle.SetStatW(0.8*ROOT.gStyle.GetStatW())
    #ROOT.gStyle.SetStatH(0.8*ROOT.gStyle.GetStatH())
    ROOT.gStyle.SetStatW(0.21)
    ROOT.gStyle.SetStatH(0.16)
    ROOT.gStyle.SetStatX(0.87)
    ROOT.gStyle.SetStatY(0.87)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    print('*** Settings:')

    bgmin = 0.4
    bgmax = 1.5
    bmin = 0.4
    bmax = 0.85
    
    #scint1 = 4.
    #scint2 = 22.
    scint1 = 0.
    scint2 = 500.
    t1 = 13.
    t2 = 35.
    if particle == 'D':
        t2 = 60.
        bgmin = 0.4
        bgmax = 0.7
        bmin = 0.4
        bmax = 0.6

    if particle == 'T':
        t2 = 85.
        bgmin = 0.2
        bgmax = 0.7
        bmin = 0.2
        bmax = 0.5
   
    # COMMON HACK!!!
    bgmin = 0.3
    bgmax = 1.5
    bmin = 0.3
    bmax = 0.85

    # ...except for electrons for calibration;)
    if particle == 'e':
        t1 = 10.
        t2 = 13.
        bgmin = 0.2
        bgmax = 2.7
        bmin = 0.9
        bmax = 1.1
        
    os.system('mkdir -p pdf png')

    # make these dictionaries keyed by momentum?
    hs = []
    momenta = []
    hscp = []
    projYs = []
    projYcps = []
    txts = []

    
    #basedir = 'TrigScint/'

    tstag = 'TS' + TStag
    if tstag == 'TS':
        tstag = 'bothTS'
    
    pbasedir = 'TrigScint_{}/'.format(particle)
    hname = 'hRef_TOF_TrigScint{}C_{}-like'.format(TStag, particle)
    
    can = None
    opt = ''
    ys = []
    eys = []
    betas = []
    betagammas = []
    ebetas = []
    leg = ROOT.TLegend(0.68, 0.5, 0.88, 0.88) # 0.12, 0.65, 0.40, 0.88)
    leg.SetBorderSize(0)
    stuff.append(leg)
    cols = [ROOT.kBlack, ROOT.kGreen+2, ROOT.kBlue, ROOT.kViolet, ROOT.kRed,
            ROOT.kOrange+1, ROOT.kYellow, ROOT.kYellow+2, ROOT.kCyan+1, ROOT.kTeal-7,
            ROOT.kAzure, ROOT.kGray+2, ROOT.kGreen+1, ROOT.kBlue, ROOT.kMagenta]
    ifile = -1
    pcans = {}
    popt = {}
    for rfile,col in zip(rfiles,cols):
        ifile = ifile + 1
        filename = rfile.GetName()
        momentum = None
        runindex = -1;
        srun = ''
        try:
            runindex = filename.index('run')
            srun = filename[runindex+6:runindex+9]
        except:
            runindex = filename.index('000')
            srun = filename[runindex+3:runindex+6]
        if momentum == None:
            momentum = getMomentum(srun)
        if momentum == None:
            momentum = getMergedMomentum(srun)
        print('---> Run {:}, p={:} MeV/c'.format(srun,momentum))
        momenta.append(momentum)
        Hs = []
        Txts = []
        ftag = filename.split('/')[-1].replace('output_','').replace('_plots_f.root','').replace('_plots.root','')

        h = rfile.Get(pbasedir + hname)
        try:
            print('   ...ok, got {} from file {}, I={}'.format(h.GetName(), rfile.GetName(), h.Integral()))
            tmp = h.GetName()
            #print('meanX: {1.2f}'.format(h.GetMean(1)))
            #print('meanY: {1.2f}'.format(h.GetMean(2)))
        except:
            print('ERROR getting histo {}{} from file {}!'.format(pbasedir,hname, rfile.GetName()))
            continue

        hs.append(h)
        if can == None and (drawElAnyway or not calibOnly):
            canname = 'WCTEJuly2023_BetheBloch_SingleFit_{}_{}{}'.format(particle, tstag, reltag)
            canname = canname.replace('_list_root','').replace('_ntuple','')
            can = ROOT.TCanvas(canname, canname, 0, 0, 1100, 800)
            cans.append(can)
            #can.Divide(8,4)
        
        h.SetStats(0)

        projY = h.ProjectionY()
        pname = f'projY_{tstag}_{ifile}_{momentum}'
        projY.SetName(pname)
        projYs.append(projY)
        if particle == 'D' or particle == 'T':
            #h.Rebin2D(2,2)
            projY.Rebin(4)
        if particle == 'p':
            #h.Rebin2D(2,2)
            projY.Rebin(2)

        beta = getBeta(ms[particle], momentum)
        betagamma = getBetaGamma(ms[particle], momentum)

        hcp = None
        if  drawElAnyway or not calibOnly:
            can.cd()
            alpha = 0.10
            h.SetFillStyle(1111)
            h.SetMarkerSize(2)
            h.SetMarkerStyle(20)
            h.SetMarkerColorAlpha(col, alpha)
            h.SetLineColorAlpha(col, alpha)
            h.SetFillColorAlpha(col,1)
            h.SetFillStyle(1111)
            h.SetMarkerSize(0.1)
            h.SetMarkerStyle(6)

            hcp = h.DrawCopy(Opt2d + opt)
            if opt == '':
                hcp.GetXaxis().SetRangeUser(t1, t2)
                hcp.GetYaxis().SetRangeUser(scint1, scint2)
            ROOT.gPad.Update()
            hscp.append(hcp)
            leg.AddEntry(h, 'p = {:4} MeV/c #beta={:1.2f}'.format(str(momentum), beta), 'F')

        projYcp = projY.Clone(projY.GetName() + f'_cp_{tstag}_{particle}_' + str(ifile))
        projYcps.append(projYcp)
        meanfull = projY.GetMean()
        sigmafull = projY.GetStdDev()

        x2 = projY.GetXaxis().GetXmax()
        x1 = max(projY.GetXaxis().GetXmin(), meanfull - sigmafull)
        projYcp.GetXaxis().SetRangeUser(x1, x2)


        if  drawElAnyway or not calibOnly:
            try:
                pcan = pcans[momentum]
            except:
                canname = 'WCTEJuly2023_LandauProfiles_{}_{}{}'.format(particle, momentum, reltag)
                #pcans[momentum] = ROOT.TCanvas(canname, canname, 300, 300, 1000, 800)
                popt[momentum] = ''
            #pcans[momentum].cd()
            #projY.Draw(popt[momentum] + ' hist')
            popt[momentum] = 'same'
        
        if projY.GetEntries() > minEntries:
            relcal = 1.
            # we expect we get a region like TS01p, so we remove the particle,
            # and add 'e' for which relative calibration was derived:
            if not calibOnly and len(calibCs) > 0:
                print('...ok, will try to use relative calibration constants!')
                #print('  ', calibCs)
                cregion = tstag + 'e'
                #print('   ', momentum, cregion)
                if momentum in calibCs and cregion in calibCs[momentum]:
                    try:
                        relcal = calibCs[momentum][cregion]
                        print('   ...using relcal={:1.4f}'.format(relcal))
                    except:
                        print(f'   ...Unable to use calibration constants for momentum {momentum}, TS {tstag} using {cregion}!')
                else:
                    print(f'   ...did NOT find calib constant for momentum {momentum}, TS {tstag} using {cregion}!')
                        
            mean = relcal*projYcp.GetMean()
            meanerr = relcal*projYcp.GetMeanError()
                        
            ys.append(mean)
            eys.append(meanerr)
            betas.append(beta)
            betagammas.append( betagamma)
            ebetas.append(0.)

        opt = 'same'
        #adjustStats(h)
        #ROOT.gPad.Update()
    if  drawElAnyway or not calibOnly:
        can.cd()
        leg.Draw()
        cnote, pnote = makePaperLabel(srun, momentum, 0.12, 0.92)
        #cnote.Draw()
        #pnote.Draw()
        pnote2 = makeMomentumLabel(srun, momentum, 0.12, 0.92)
        #pnote2.Draw()
        
    #parts = ['e', 'mu', 'pi', 'K', 'p', 'D', 'T']
    #lines = makeLines(h, 0., parts, momentum, True)
    #stuff.append(lines)
    #stuff.append([cnote, pnote, pnote2])

    ####################################################################################
    canp = None
    legp = None
    if drawElAnyway or not calibOnly:
        canname = 'WCTEJuly2023_LandauProfiles_{}_{}{}'.format(particle, tstag, reltag)
        canp = ROOT.TCanvas(canname, canname, 200, 200, 1000, 800)
        canp.cd()
        opt = ''
        legp = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
        legp.SetBorderSize(0)
        stuff.append(legp)
        for h,col,projY,momentum,beta in zip(hs,cols,projYs,momenta,betas):
            hname = h.GetName()
            #print(f'col: {col}')
            projY.SetLineColor(col)
            projY.SetFillColorAlpha(col, 0.3)
            #projY.SetFillStyle(1111)
            projY.SetLineWidth(2)
            projY.SetLineStyle(1)
            projY.Rebin(2)
            val = projY.Integral(0, projY.GetXaxis().GetNbins()+1)
            if val > 0.:
                projY.Scale(1./val)
            projY.SetMaximum(projY.GetMaximum()*1.2)
            if 'D' in hname:
                projY.SetMaximum(projY.GetMaximum()*1.2)
            projY.SetStats(0)
            projY.Draw(Opt1d + opt)
            if opt == '':
                projY.GetXaxis().SetRangeUser(scint1, scint2)
            legp.AddEntry(projY, 'p = {:4} MeV/c #beta={:1.2f}'.format(str(momentum), beta), 'F')
            opt = 'same'
        legp.Draw()
        canp.Update()
        cans.append(canp)

    gcan = None
    grb = None
    grbg = None
    if drawElAnyway or not calibOnly:
    ####################################################################################
        drawBetaOnly = True
        cw = 600
        dw = 50
        if not drawBetaOnly:
            cw = 2*cw + dw
        canname = 'BetaGraph_{}_{}{}'.format(particle, tstag, reltag)
        gcan = ROOT.TCanvas(canname, canname, 100, 100, cw, 600)
        if not drawBetaOnly:
            gcan.Divide(2,1)
        print('* betas, ebetas, ys, eys:')
        print(betas, ebetas, ys, eys)
        grb = MakeGraph(betas, ebetas, ys, eys)
        grbg = MakeGraph(betagammas, ebetas, ys, eys)

        hn = 'tmpbg' + tstag + particle 
        ht = hn + ';#beta#gamma;Mean trig. scint. charge [a.u.];'
        htmpbg = ROOT.TH2D(hn, ht, 100, bgmin, bgmax, 100, scint1, scint2)
        htmpbg.SetStats(0)
        htmpbg.GetXaxis().SetMoreLogLabels()

        if not drawBetaOnly:
            gcan.cd(1)
            htmpbg.Draw()
            ROOT.gStyle.SetOptTitle(0)
            #ROOT.gPad.SetLogx()
            ROOT.gPad.SetGridx(1)
            ROOT.gPad.SetGridy(1)
            grbg.Draw("P")
            #adjustStats(grbg)

        hn = 'tmpb' + tstag + particle
        ht = hn + ';#beta;Mean trig. scint. charge [a.u.];'
        htmpb = ROOT.TH2D(hn, ht, 100,bmin, bmax, 100, scint1, scint2)
        htmpb.SetStats(0)
        htmpb.GetXaxis().SetMoreLogLabels()

        if not drawBetaOnly:
            gcan.cd(2)
        htmpb.Draw()
        ROOT.gStyle.SetOptTitle(0)
        #ROOT.gPad.SetLogx()
        ROOT.gPad.SetGridx(1)
        ROOT.gPad.SetGridy(1)
        grb.Draw("P")
        #adjustStats(grb)
        #fun = ROOT.TF1('fun', '[0]/x^2 + [1]', 0.1, 1.)
        #fun.SetParameters(0.1, 1.)
        #fun = ROOT.TF1('fun', '[0]/x^2*(log([1]*x/sqrt(1-x*x)) - x^2) + [2]', 0.1, 1.)
        #fun.SetParameters(2., 10., 0.5)
        fun = ROOT.TF1('fun', '[0]/x^2*(log([3]*x*x/[1]/(1-x*x)) - x^2) + [2]', 0.1, 1.)
        fun.SetParameters(10, 700., 10.)
        fun.SetParName(0, 'A')
        fun.SetParName(1, 'I')
        fun.SetParName(2, 'C')
        fun.SetParName(3, 'me')
        fun.FixParameter(3, 2*0.511e6)
        #fun.SetParName(2, 'g')

        
        ##grb.Fit('fun', '', '0')
        ##fun.Draw('same')
        
        cans.append(gcan)
        stuff.append([momenta, grb, grbg, htmpb, htmpbg, projYs, projYcps, fun])
        gcan.Update()

    
    return momenta, grb, grbg, cans, pcans, projYs, projYcps, leg, legp

###################################
###################################
###################################

def main(argv):

    # execute only if run as a script"
    gBatch = False
    gTag=''

    pngdir = 'png_results/'
    pdfdir = 'pdf_results/'
    os.system(f'mkdir {pngdir}')
    os.system(f'mkdir {pdfdir}')

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    if len(argv) < 1:
        PrintUsage(argv)
        return

    #ROOT.gStyle.SetPalette(ROOT.kSolar)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    #ROOT.gStyle.SetPalette(1)

    allTStags = [ '00', '01', '02', '03', '10', '11', '12', '13' ]
    
    TStags = [
        # not supported anymore!
        #'0',
        #'1'
        #'' # both trigger scintillators
        # NOW supported:
        # individual TS PMTs:
        ##
        #'00',
        '01',
        ##
        #'02',
        '03',
        '10','11',
        '12',
        ##
        #'13',
    ]
    extraTag = ''
    if len(allTStags) != len(TStags):
        extraTag = '_rm'
        for atag in allTStags:
            if not atag in TStags:
                extraTag = extraTag + f'_{atag}'
    particles = [ 'e', # for calibration
                  'p', # protons
                  'D', # deuterons
                  #'T', # tritium,
                 ]
    rfiles = readInputFiles()
    stuff.append(rfiles)
    # graphs of E losses in p.e. as function of beta or beta*gamma
    GrsBeta = {}
    GrsBg = {}
    Cans = []
    projs = {}
    momenta = []
    calibCs = {}
    for particle in particles:
        calibOnly = False
        if particle == 'e':
            calibOnly = True
        for TStag in TStags:
            region = 'TS' + TStag + particle

            # based on chi2 / Npoints:
            #if region == 'TS00p' or region == 'TS03p' or region == 'TS12p':
            #    continue
            #if region == 'TS12p':
            #    continue
            
            
            print(f'Adding region {region}')
            momenta, grbeta, grbg, cans, pcans, projYs, projYcps, leg, legp = singleFit(sys.argv, rfiles, TStag, particle, calibOnly, calibCs)
            if not calibOnly:
                GrsBeta[region] = grbeta
            stuff.append([grbeta, grbg, cans, pcans, leg, legp])
            projs[region] = [projYs, projYcps]
            Cans.append(cans)

            if calibOnly:
                chargeVals = {}
                for region in projs:
                    for momentum,proj in zip(momenta,projs[region][0]):
                        val = proj.GetMean()
                        print('p={:1.0f} region {:} mean charge [Npe]: {:1.3f}'.format(momentum,region,val))
                        try:
                            chargeVals[momentum][region] = 1.*val
                        except:
                            chargeVals[momentum] = {}
                            chargeVals[momentum][region] = 1.*val
                #print(chargeVals)
                print('Working on relative TS calibrations constants...')
                for momentum in chargeVals:
                    aver = 0.
                    for region in chargeVals[momentum]:
                        val = chargeVals[momentum][region]
                        aver = aver + val
                    nn = len(chargeVals[momentum])
                    #print(f'nPMTs: {nn}')
                    if nn > 0:
                        aver = aver / nn
                    else:
                        aver = -1.
                    calibCs[momentum] = {}
                    for region in chargeVals[momentum]:
                        val = chargeVals[momentum][region]
                        if val > 0.:
                            calibCs[momentum][region] = aver/val
                        else:
                            calibCs[momentum][region] = -1
                #print(calibCs)
                print('"TOF" trigger scintillators relative calibration constants:')
                xs = {}
                ys = {}
                imomentum = -1
                for momentum in chargeVals:
                    imomentum = imomentum + 1
                    print(f'{momentum}: ', end='')
                    for region in chargeVals[momentum]:
                        c = calibCs[momentum][region]
                        if imomentum == 0:
                            xs[region] = []
                            ys[region] = []
                        xs[region].append(1.*momentum)
                        ys[region].append(1.*c)
                        print(' {:}: {:1.4f}'.format(region, c), end='')
                    print()
                cgrs = {}
                ireg = -1
                for region in xs:
                    ireg = ireg + 1
                    print(xs[region], ys[region])
                    col = ROOT.kBlue - ireg
                    if 'TS0' in region:
                        col = ROOT.kMagenta - ireg
                    cgrs[region] = MakeGraphNoErrs(xs[region], ys[region], col, 20 + (ireg % 4))

        if calibOnly:
            canname = 'WCTEJuly2023_CalibTS_eRelCalib{}'.format(len(cgrs))
            ccan = ROOT.TCanvas(canname, canname, 200, 200, 1100, 800)
            cans.append(ccan)
            ccan.cd()
            opt = 'PL'
            hname = 'calibGr'
            htitle = ';p [MeV/c];rel. calib const.'
            hh2 = ROOT.TH2D(hname, htitle, 100, 400., 1600., 100, 0., 1.7)
            hh2.SetStats(0)
            hh2.Draw()
            cleg = ROOT.TLegend(0.79, 0.23, 0.89, 0.89)
            cleg.SetBorderSize(0)
            stuff.append(cleg)
            for region,cgr in cgrs.items():
                print(f'plotting calib graph for {region}')
                cgr.Draw(opt)
                cleg.AddEntry(cgr, region, 'PL')
            cleg.Draw()
            ROOT.gPad.Update()
            # Then, using p and D, multi-fit individual 8 regions 00..13 usin the calibration, or apply calibration 'constants' in event Loop?
            # are these really constants? Don't thet include some delta electrons physics and momemntum dependence, too?
            # fit peaks, and not use means of projections?
            # End of relative calibration using electrons

    ##########################
    #  Now run the multifit! #
    ##########################
    step = 0.01
    debug = 1
    refmomentum = 1000
    nCalibCs = len(calibCs[refmomentum])
    fitter, result, npars, parNames = doTheFit(GrsBeta, nCalibCs, step, debug)

    # analyze the fitter parameters
    # plot individual subfits over data in each region!

    canres, cancmp, leg, hres, hdEFitOverTheory, hb, GrsFit = AnalyzeFitResults(GrsBeta, fitter, result, nCalibCs, parNames, extraTag)
    cans.append(canres)
    cans.append(cancmp)

    # End of multifit

    # And just print all canvases;)
    for cans in Cans:
        for can in cans:
            try:
                can.cd()
                if 'vs' in can.GetName():
                    pnote.Draw()            
                can.Update()
                can.Print(pngdir + can.GetName() + extraTag + '.png')
                can.Print(pdfdir + can.GetName() + extraTag + '.pdf')
            except:
                print('ERROR printing canvas!')
    
    if not gBatch:
        ROOT.gApplication.Run()

    
    print('DONE!')
    return

###################################
###################################
###################################





###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################

