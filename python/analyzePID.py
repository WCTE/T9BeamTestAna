#!/snap/bin/pyroot

#/usr/bin/python3

# jk
# 20/09/2022
# 14.7.2023
# 5.3.2024

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from collections import OrderedDict

from labelTools import *

cans = []
stuff = []
lines = []

####################################################################################

class cFitPeak:
    def __init__(self, p, x0, sigma):
        self.p = p
        self.x0 = x0
        self.sigma = sigma
    


####################################################################################
def readInputFiles():
    dirname = 'histos/windowpe_analyzed/'
    filenames = [
        'peakAnalysed_timeCorr_windInt_000403_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000396_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000394_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000393_plots_f.root',
        'peakAnalysed_timeCorr_windInt_000392_plots_f.root',
        
        # 900:
        'peakAnalysed_timeCorr_windInt_000398_plots_f.root',
        #'peakAnalysed_timeCorr_windInt_000453_plots_f.root'
        # add 370??
        
        'peakAnalysed_timeCorr_windInt_000399_plots_f.root',

        # 700:
        'peakAnalysed_timeCorr_windInt_000402_plots_f.root',
        #'peakAnalysed_timeCorr_windInt_000438_plots_f.root',
        #'peakAnalysed_timeCorr_windInt_000457_plots_f.root',
        
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



def PrintUsage(argv):
    print('Usage:')
    print('{} filename_plots.root [-b]'.format(argv[0]))
    print('Example:')
    print('{} output_300n_plots.root -b'.format(argv[0]))
    return

##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    pngdir = 'png_results/'
    pdfdir = 'pdf_results/'
    os.system(f'mkdir {pngdir}')
    os.system(f'mkdir {pdfdir}')

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


    #ROOT.gStyle.SetPalette(ROOT.kSolar)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    #ROOT.gStyle.SetPalette(1)

    
    #filename = 'output_300n_plots.root'
    #filename = argv[1]
    #rfile = ROOT.TFile(filename, 'read')


    rfiles = readInputFiles()
    stuff.append(rfiles)


    os.system('mkdir -p pdf png')

    hnames2d = [
        'hRef_pbC_act23C',
        'hRef_pbC_act1C',
        'hRef_ACT1CACT23C'
    ]


    
    pbasedirs = [
        #'Charged/'
        'Charged_nonp/'
        #'TrigScint_p/',
        #'TrigScint_e/',
    ]


    canname = 'WCTEJuly2023_Quick2D_all_PbG'
    cw = 1100
    ch = 800
    allcan = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
    cans.append(allcan)
    allsame = ''
    allleg = ROOT.TLegend(0.7, 0.6, 0.88, 0.88)

    
    hs = []
    fitPeaks = {}
    for pbasedir in pbasedirs:

        pTag = pbasedir.replace('Charged_','').replace('Charged','').replace('/','')
        selTag = ''
        if len(pTag) > 0:
            selTag = '_' + pTag + '-like'
        print(f'pTag: {pTag} selTag: {selTag}')
        jcan = -1
                    
        for hname in hnames2d:
            jcan = jcan + 1
                

            canname = f'WCTEJuly2023_Quick2D_PID_{hname}_{selTag}'
            can = ROOT.TCanvas(canname, canname, jcan*100, jcan*100, 1200+200, 800)
            #can.Divide(4,3)
            can.Divide(2,1)
            cans.append(can)
        
       
            ican = -1
            for rfile in rfiles:
                filename = rfile.GetName()
                ican = ican + 1
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

                h = rfile.Get(pbasedir + hname + selTag)
                try:
                    #print('ok, got ', h.GetName())
                    tmp = h.GetName()
                except:
                    print('ERROR getting histo {}{}!'.format(pbasedir,hname))
                    continue

                #print('Pushing ', ich, hname)
                hs.append(h)


                """
                canname = 'WCTEJuly2023_Quick2D_{}_{}'.format(ftag, hname)
                canname = canname.replace('_list_root','').replace('_ntuple','').replace('.root','')
                cw = 1100
                ch = 800
               
                can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
                cans.append(can)
                #can.Divide(8,4)
                #h.Rebin2D(2,2)
                """

                
                can.cd(ican+1)
                opt = 'colz'
                is2d = True
                h.SetStats(0)
                if hname == 'hRef_ACT1CACT23C':
                    h.GetXaxis().SetRangeUser(0, 30.)
                    h.GetYaxis().SetRangeUser(0, 40.)
                else:
                    h.GetXaxis().SetRangeUser(0, 500.)
                    if hname == 'hRef_pbC_act1C':
                        h.GetYaxis().SetRangeUser(0, 30.)
                    else:
                        h.GetYaxis().SetRangeUser(0, 50.)
                    
                h.Draw(opt)
                ROOT.gPad.SetLogz(1)

                #ROOT.gPad.Update()
                cnote, pnote = makePaperLabel(srun, momentum, 0.12, 0.92)
                #cnote.Draw()
                #pnote.Draw()
                pnote2 = makeMomentumLabel(srun, momentum, 0.12, 0.92)
                pnote2.Draw()
                ROOT.gPad.Update()
                stuff.append([cnote, pnote, pnote2])
                
                rho = h.GetCorrelationFactor()
                rtxt = ROOT.TLatex(0.76, 0.85, '#rho={:1.2f}'.format(rho))
                rtxt.SetNDC()
                rtxt.SetTextSize(0.04)
                rtxt.Draw()
                stuff.append(rtxt)
                
                projX = h.ProjectionX(srun + hname + '_projX')
                projY = h.ProjectionX(srun + hname + '_projY')

                """
                canname = canname + '_projX'
                can = ROOT.TCanvas(canname, canname, ican*30, ican*30, 800, 600)
                cans.append(can)
                """

                """
                can.cd(ican+1)
                projX.Draw('hist')
                allcan.cd()
                projXcp = projX.DrawCopy('hist plc' + allsame)
                projXcp.SetStats(0)
                projXcp.SetLineWidth(2)
                projXcp.Scale(1./projXcp.Integral(0,projXcp.GetXaxis().GetNbins()+1))
                allleg.AddEntry(projXcp, 'Run {}, p={} MeV/c'.format(srun, momentum), 'L')
                allsame = 'same'
                projXcp.SetMaximum(0.085)
                stuff.append(projXcp)
                #ca.cd()
                can.cd(ican+1)
                ibx = projX.GetMaximumBin()
                #projX.GetXaxis().SetRangeUser(0.,projX.GetXaxis().GetXmax())
                print(ibx)
                xmax = projX.GetBinCenter(ibx)
                rms = projX.GetStdDev()
                x1 = xmax - rms/5.
                x2 = xmax + rms/5.
                fitname = 'fit_{}_{}_{}'.format(srun, momentum, hname)
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
                sigma = abs(fun.GetParameter(2))
                adjustStats(projX)
                #ROOT.gPad.Update()

               
                """

                if 'TOF' in hname:
                    parts = ['e', 'mu', 'pi', 'K', 'p', 'D', 'T']
                    lines = makeLines(h, 0., parts, momentum, True)
                    stuff.append(lines)
                

        
    ##################################
    #       plots all the canvas     #
    ##################################

    allcan.cd()
    allleg.Draw()
    
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

