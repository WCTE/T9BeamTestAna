#!/snap/bin/pyroot

#/usr/bin/python3

# jk
# 20/09/2022
# 14.7.2023
#  5.3.2024
# 15.5.2024

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from collections import OrderedDict

from labelTools import *

MeV = 1.

cans = []
stuff = []
lines = []

def SetStyle(gr, mst, msz, mc):
    gr.SetMarkerStyle(mst)
    gr.SetMarkerSize(msz)
    gr.SetMarkerColor(mc)
    gr.SetLineColor(mc)

####################################################################################

class cFitPeak:
    def __init__(self, p, x0, muerr, sigma, sigmaErr):
        self.p = p
        self.x0 = x0
        self.muerr = muerr
        self.sigma = sigma
        self.sigmaErr = sigmaErr
    
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
def PrintUsage(argv):
    print('Usage:')
    print('{} filename_plots.root [-b]'.format(argv[0]))
    print('Example:')
    print('{} output_300n_plots.root -b'.format(argv[0]))
    return

####################################################################################
####################################################################################
####################################################################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm

def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    pngdir = 'png_results/'
    pdfdir = 'pdf_results/'
    os.system(f'mkdir -p {pngdir}')
    os.system(f'mkdir -p {pdfdir}')

    opt2d = 'colz'

    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    #gBatch = True
    gBatch = False
    gTag=''
    print(argv[1:])
    try:
        # options that require an argument should be followed by a colon (:).
        opts, args = getopt.getopt(argv[2:], 'hbt:', ['help','batch','tag='])
        print('Got options:')
        print(opts)
        print(args)
    except getopt.GetoptError:
        print('Parsing...')
        print ('Command line argument error!')
        print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]]'.format(argv[0]))
        sys.exit(2)
    for opt,arg in opts:
        print('Processing command line option {} {}'.format(opt,arg))
        if opt == '-h':
            print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]'.format(argv[0]))
            sys.exit()
        elif opt in ("-b", "--batch"):
            gBatch = True
            print('OK, running in batch mode')
        elif opt in ("-t", "--tag"):
            gTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gTag,) )

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    if len(argv) < 1:
        PrintUsage(argv)
        return

    ROOT.gStyle.SetOptFit(111)
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    ROOT.gStyle.SetPadLeftMargin(0.15)
        
    #ROOT.gStyle.SetPalette(ROOT.kSolar)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    #ROOT.gStyle.SetPalette(1)

    rfiles = readInputFiles()
    stuff.append(rfiles)


    os.system('mkdir -p pdf png')

    # more histos NOT supported!
    hnames2d = [ 
        'hRef_pbC_TrigScintC',
    ]
    pbasedirs = [
        #'TrigScint_nonp/',
        'TrigScint_e/',
        #'TrigScint_nonp/',
    ]

    canname = 'WCTEJuly2023_Quick2D_all_PbG'
    cw = 1100
    ch = 800
    allcan = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
    cans.append(allcan)
    allsame = ''
    allleg = ROOT.TLegend(0.7, 0.6, 0.88, 0.88)

    canname = 'WCTEJuly2023_Quick2D_PbG_fits'
    can = ROOT.TCanvas(canname, canname, 200, 200, 1200, 800)
    can.Divide(4,3)
    cans.append(can)
    
    hs = []
    fitPeaks = {}
    ican = -1    
    for pbasedir in pbasedirs:
        
        suff = '-like'
        particle = ''
        if '_e/' in pbasedir:
            particle = 'e'
        if '_p/' in pbasedir:
            particle = 'p'
        if '_nonp/' in pbasedir:
            particle = 'nonp'
        if '_pi/' in pbasedir:
            particle = 'pi'
        if '_mu/' in pbasedir:
            particle = 'mu'
        if '_D/' in pbasedir:
            particle = 'D'
        if '_T/' in pbasedir:
            particle = 'T'
        suff = '_' + particle + suff

        allleg.SetHeader(f'Particle: {particle}')
        
        fitPeaks[particle] = {}
        for rfile in rfiles:
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
            print(srun,momentum)

            ftag = filename.split('/')[-1].replace('output_','').replace('_plots.root','')
      
            for hname in hnames2d:
                ican = ican + 1
                
                h = rfile.Get(pbasedir + hname + suff)
                try:
                    #print('ok, got ', h.GetName())
                    tmp = h.GetName()
                except:
                    print('ERROR getting histo {}{}!'.format(pbasedir,hname + suff))
                    continue

                #print('Pushing ', ich, hname)
                hs.append(h)


                """
                canname = 'WCTEJuly2023_Quick2D_{}_{}'.format(ftag, hname + suff)
                canname = canname.replace('_list_root','').replace('_ntuple','').replace('.root','')
                cw = 1100
                ch = 800
               
                can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
                cans.append(can)
                #can.Divide(8,4)
                #h.Rebin2D(2,2)
                opt = 'colz'
                is2d = True
                h.SetStats(0)
                h.Draw(opt)
                rho = h.GetCorrelationFactor()
                rtxt = ROOT.TLatex(0.76, 0.85, '#rho={:1.2f}'.format(rho))
                rtxt.SetNDC()
                rtxt.SetTextSize(0.04)
                rtxt.Draw()
                stuff.append(rtxt)
                """
                
                projX = h.ProjectionX(srun + hname + suff + '_projX')

                """
                canname = canname + '_projX'
                can = ROOT.TCanvas(canname, canname, ican*30, ican*30, 800, 600)
                cans.append(can)
                """
                
                can.cd(ican+1)
                projX.Draw('hist')
                allcan.cd()
                projXcp = projX.DrawCopy('hist plc' + allsame)
                projXcp.SetStats(0)
                projXcp.SetLineWidth(2)
                num = projXcp.Integral(0,projXcp.GetXaxis().GetNbins()+1)
                if num > 0:
                    projXcp.Scale(1./num)
                allleg.AddEntry(projXcp, 'Run {}, p={} MeV/c'.format(srun, momentum), 'L')
                allsame = 'same'
                projXcp.SetMaximum(0.085)
                stuff.append(projXcp)
                #ca.cd()
                can.cd(ican+1)
                ChargeCenter = 40. # some dummy val
                if particle == 'e':
                    # a*500 + b = 7
                    # a*1200 + b = 20
                    # ==>
                    # a*700 = 13 ==>
                    a = 13./700.
                    b = 20 - 1200*a
                    ChargeCenter = abs(momentum)*a + b
                    #if abs(momentum) < 600:
                    #    ChargeCenter = 5.
                    print(f'momentum: {momentum}, ChargeCenter={ChargeCenter}')
                    #projX.GetXaxis().SetRangeUser(ChargeCenter,projX.GetXaxis().GetXmax())
                if particle == 'p':
                    projX.GetXaxis().SetRangeUser(0., 200.)
                ibx = projX.GetMaximumBin()
                if particle == 'e':
                    projX.GetXaxis().SetRangeUser(0.,projX.GetXaxis().GetXmax())
                print(ibx)
                xmax = projX.GetBinCenter(ibx)
                #rms = projX.GetStdDev()
                rms = 4.
                x1 = ChargeCenter - rms
                x2 = ChargeCenter + rms
                print(f'INIT Charge limits: {x1} {ChargeCenter} {x2}')
                fitname = 'fit_{}_{}_{}_{}'.format(srun, momentum, hname, suff)
                fun = ROOT.TF1(fitname, '[0]*exp(-(x-[1])^2/(2*[2]^2))', projX.GetXaxis().GetXmin(), projX.GetXaxis().GetXmax())
                fun.SetNpx(1000)
                fun.SetParameters(projX.GetMaximum()/2, xmax, 2.)
                projX.Fit(fitname, '', '', x1, x2)
                x0 = fun.GetParameter(1)
                sigma = abs(fun.GetParameter(2))
                sf = 1.1
                if x0 > 50:
                    sf = 1.
                xx1, xx2 = x0 - sf*sigma, x0 + sf*sigma
                print(xx1,xx2)
                projX.Fit(fitname, '', '', xx1, xx2)
                fun.Draw('same')
                stuff.append(fun)

                x0 = fun.GetParameter(1)
                muerr = fun.GetParError(1)
                sigma = abs(fun.GetParameter(2))
                sigmaErr = abs(fun.GetParError(2))

                print(f'MOMENTUM {momentum} PART {suff} FITTED MAIN PEAK MEAN {x0} RMS {sigma}')
                fitPeaks[particle][momentum] = cFitPeak(momentum, x0, muerr, sigma, sigmaErr)
                
                adjustStats(projX)
                #ROOT.gPad.Update()
                cnote, pnote = makePaperLabel(srun, momentum, 0.12, 0.92)
                #cnote.Draw()
                #pnote.Draw()
                pnote2 = makeMomentumLabel(srun, momentum, 0.12, 0.92)
                pnote2.Draw()
                ROOT.gPad.Update()
                
                if 'TOF' in hname:
                    parts = ['e', 'mu', 'pi', 'K', 'p', 'D', 'T']
                    lines = makeLines(h, 0., parts, momentum, True)
                    stuff.append(lines)
                stuff.append([cnote, pnote, pnote2])

        
    ##################################
    #       plot all the canvas     #
    ##################################

    allcan.cd()
    allleg.Draw()

    canname = 'PbGResolution'
    cw = 1200
    ch = 400
    canreso = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
    canreso.Divide(3,1)
    cans.append(canreso)
    grsReso = {}
    grsResoRel = {}
    grsE = {}
    
    for particle in fitPeaks:
        grsReso[particle] = ROOT.TGraphErrors()
        grsResoRel[particle] = ROOT.TGraphErrors()
        grsE[particle] = ROOT.TGraphErrors()

        ip = -1
        momenta = []
        for momentum in fitPeaks[particle]:
            momenta.append(momentum)
            ip = ip+1
            fitPeak = fitPeaks[particle][momentum]
            print('  {' + '{}, {:.3f}, {:.3f}'.format(fitPeak.p, fitPeak.x0, fitPeak.sigma) + '}, ')

            # resolution of fitted charges
            grsReso[particle].SetPoint(ip, momentum, fitPeak.sigma)
            grsReso[particle].SetPointError(ip, 0., fitPeak.sigmaErr)
            SetStyle(grsReso[particle], 20, 1, ROOT.kBlue)

            # relative resolution, w.r.t. the nominal momentum
            #grsResoRel[particle].SetPoint(ip, momentum, fitPeak.sigma / momentum)
            #grsResoRel[particle].SetPointError(ip, 0., fitPeak.sigmaErr / momentum)
            grsResoRel[particle].SetPoint(ip, momentum, fitPeak.sigma / fitPeak.x0)
            grsResoRel[particle].SetPointError(ip, 0., fitPeak.sigmaErr / fitPeak.x0)
            SetStyle(grsResoRel[particle], 20, 1, ROOT.kGreen+2)

            # fitted charges -- linearity check as function as nominal momentum
            grsE[particle].SetPoint(ip, momentum, fitPeak.x0)
            grsE[particle].SetPointError(ip, 0., fitPeak.muerr)
            SetStyle(grsE[particle], 20, 1, ROOT.kRed)
            

    opt = 'P'
    helphs = {}
    dp = 60.*MeV
    p1, p2 = min(momenta)-dp, max(momenta)+dp
    calibFits = {}
    for particle in fitPeaks:

        helphs[particle] = []

        canreso.cd(1)
        helphs[particle].append(ROOT.TH2D('reso_h2_' + particle,';p [MeV/c];#sigma_{Charge} [N_{p.e.}]', 100, p1, p2, 100, 0, 25./10.))
        helphs[particle][-1].SetStats(0)
        helphs[particle][-1].Draw()
        grsReso[particle].Draw(opt)
        
        canreso.cd(2)
        helphs[particle].append(ROOT.TH2D('reso_h2_' + particle,';p [MeV/c];#sigma_{Charge} / Charge [-]', 100, p1, p2, 100, 0., 0.129))
        helphs[particle][-1].SetStats(0)
        helphs[particle][-1].Draw()
        grsResoRel[particle].Draw(opt)

        canreso.cd(3)        
        helphs[particle].append(ROOT.TH2D('reso_h2_' + particle,';p [MeV/c];Fitted mean charge [N_{p.e.}]', 100, p1, p2, 100, 0., 400/10.))
        helphs[particle][-1].SetStats(0)
        helphs[particle][-1].Draw()
        grsE[particle].Draw(opt)
        fun1 = ROOT.TF1('fun1_' + particle, '[0] + [1]*x', p1, p2)
        fun1.SetLineStyle(2)
        calibFits[particle] = fun1
        grsE[particle].Fit(fun1)
        opt = 'P'

            
            
    for can in cans:
        try:
            can.cd()
        except:
            print('ERROR printing canvas')
            continue
        if 'vs' in can.GetName():
            pnote.Draw()            
        can.Update()
        can.Print(pngdir + can.GetName() + '.png')
        ###can.Print(pdfdir + can.GetName() + '.pdf')
    
    if not gBatch:
        ROOT.gApplication.Run()
    return



###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################

