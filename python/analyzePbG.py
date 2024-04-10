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
def readInputFiles():
    dirname = 'histos/windowpe_analyzed/'
    filenames = [
        'peakAnalysed_timeCorr_windInt_-16_45_000403_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000396_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000394_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000393_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000392_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000398_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000399_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000402_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000449_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000436_plots_f.root',
        'peakAnalysed_timeCorr_windInt_-16_45_000435_plots_f.root'
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

    if len(argv) < 2:
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

        Hs = []
        Txts = []
        ftag = filename.split('/')[-1].replace('output_','').replace('_plots.root','')

        os.system('mkdir -p pdf png')
        hs = []
        txts = []
        hnames2d = [ 
                     'hRef_pbC_TrigScintC',
                    ]
        pbasedirs = ['TrigScint_e/',
                     #'TrigScint_p/'
                    ]

        meanXmap = OrderedDict()

        for pbasedir in pbasedirs:

            suff = '-like'
            if '_e/' in pbasedir:
                suff = '_e' + suff 
            if '_p/' in pbasedir:
                suff = '_p' + suff 
            if '_pi/' in pbasedir:
                suff = '_pi' + suff 
            if '_mu/' in pbasedir:
                suff = '_mu' + suff
            if '_D/' in pbasedir:
                suff = '_D' + suff 
            if '_T/' in pbasedir:
                suff = '_T' + suff 

            for hname in hnames2d:
                h = rfile.Get(pbasedir + hname + suff)
                try:
                    #print('ok, got ', h.GetName())
                    tmp = h.GetName()
                except:
                    print('ERROR getting histo {}{}!'.format(pbasedir,hname + suff))
                    continue

                #print('Pushing ', ich, hname)
                hs.append(h)

                
                canname = 'WCTEJuly2023_Quick2D_{}_{}'.format(ftag, hname + suff)
                canname = canname.replace('_list_root','').replace('_ntuple','').replace('.root','')
                cw = 1100
                ch = 800
                """
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

                canname = canname + '_projX'
                can = ROOT.TCanvas(canname, canname, 300, 300, 800, 600)
                cans.append(can)
                can.cd()
                projX.Draw('hist')
                ibx = projX.GetMaximumBin()
                print(ibx)
                xmax = projX.GetBinCenter(ibx)
                rms = projX.GetStdDev()
                x1 = xmax - rms/5.
                x2 = xmax + rms/5.
                fitname = 'fit_{}_{}_{}_{}'.format(srun, momentum, hname, suff)
                fun = ROOT.TF1(fitname, '[0]*exp(-(x-[1])^2/(2*[2]^2))', projX.GetXaxis().GetXmin(), projX.GetXaxis().GetXmax())
                fun.SetNpx(1000)
                fun.SetParameters(projX.GetMaximum()/2, xmax, 2.)
                projX.Fit(fitname, '', '', x1, x2)
                x0 = fun.GetParameter(1)
                sigma = fun.GetParameter(2)
                sf = 1.5
                if x0 > 50:
                    sf = 2.
                xx1, xx2 = x0 - sf*sigma, x0 + sf*sigma
                print(xx1,xx2)
                projX.Fit(fitname, '', '', xx1, xx2)
                fun.Draw('same')
                stuff.append(fun)

                x0 = fun.GetParameter(1)
                sigma = fun.GetParameter(2)
                print(f'MOMENTUM {momentum} PART {suff} FITTED MAIN PEAK MEAN {x0} RMS {sigma}')

                #adjustStats(h)
                #ROOT.gPad.Update()
                cnote, pnote = makePaperLabel(srun, momentum, 0.12, 0.92)
                #cnote.Draw()
                #pnote.Draw()
                pnote2 = makeMomentumLabel(srun, momentum, 0.12, 0.92)
                pnote2.Draw()
                if 'TOF' in hname:
                    parts = ['e', 'mu', 'pi', 'K', 'p', 'D', 'T']
                    lines = makeLines(h, 0., parts, momentum, True)
                    stuff.append(lines)
                stuff.append([cnote, pnote, pnote2])

        
##################################
#       plots all the canvas     #
##################################

 
    
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

