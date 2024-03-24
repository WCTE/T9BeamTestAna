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
    filename = argv[1]
    rfile = ROOT.TFile(filename, 'read')

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
    hnames2d = [ #'hRef_TOFPbA
                 'hRef_TOFPbC',
                 
                 'hRef_TOF_TrigScintC',
                 'hRef_TOF_TrigScint0C',
                 'hRef_TOF_TrigScint1C',

                 'hRef_TOF_TrigScint0LC',
                 'hRef_TOF_TrigScint0RC',
                 'hRef_TOF_TrigScint1LC',
                 'hRef_TOF_TrigScint1RC',

   # L-L, R-R
                 'hRef_TrigScint0LC_TrigScint1LC',
                 'hRef_TrigScint0RC_TrigScint1RC',
                 'hRef_TrigScint0LC_TrigScint1LC_p-like',
                 'hRef_TrigScint0RC_TrigScint1RC_p-like',
                 
    #             ]
    #dummy = [

                 'hRef_pbC_TrigScintC',
                 #'hRef_pbC_TrigScintA',
                 'hRef_pbA_TrigScintC',
                 #'hRef_pbA_TrigScintA',
                 
                 'hRef_TrigScint0RC_TrigScint0LC',
                 'hRef_TrigScint1RC_TrigScint1LC',
                 'hRef_TrigScint0RC_TrigScint0LC_p-like',
                 'hRef_TrigScint1RC_TrigScint1LC_p-like',

                 'hRef_TrigScint0_Cweighted_xymap',
                 'hRef_TrigScint1_Cweighted_xymap',

                 'hRef_TrigScint0_Cweighted_xymap_p-like',
                 'hRef_TrigScint1_Cweighted_xymap_p-like',

                 'hRef_TrigScint0_Tweighted_xymap',
                 'hRef_TrigScint1_Tweighted_xymap',

                 'hRef_TrigScint0_Tweighted_xymap_p-like',
                 'hRef_TrigScint1_Tweighted_xymap_p-like',

                 # amplitudes:
                 #'hRef_act0LA_act0RA_nonZero',
                 #'hRef_act1LA_act1RA_nonZero',
                 #'hRef_act2LA_act2RA_nonZero',
                 #'hRef_act3LA_act3RA_nonZero',
                 #'hRef_act3LC_act3RC_nonZero',

                 # charges:
                 'hRef_act0LC_act0RC_nonZero',
                 'hRef_act1LC_act1RC_nonZero',
                 'hRef_act2LC_act2RC_nonZero',
                 'hRef_act3LC_act3RC_nonZero',

                 #'hRef_act0LC_minus_act0RC_nonZero',
                 #'hRef_act1LC_minus_act1RC_nonZero',
                 #'hRef_act2LC_minus_act2RC_nonZero',
                 #'hRef_act3LC_minus_act3RC_nonZero',

                 'hRef_act0Xmap_nonZero',
                 'hRef_act1Xmap_nonZero',
                 'hRef_act2Xmap_nonZero',
                 'hRef_act3Xmap_nonZero',
                 
                ]
    
    basedir = 'TrigScint/'
    pbasedir = 'TrigScint_p/'
    chbasedir = 'Charged/'

    meanXmap = OrderedDict()
    
    for hname in hnames2d:
        rdir = basedir
        if 'p-like' in hname:
            rdir = pbasedir
        elif 'nonZero' in hname or 'Xmap' in hname:
            rdir = chbasedir
        elif 'hRef_TOFPb' in hname:
            rdir = chbasedir
        h = rfile.Get(rdir + hname)
        try:
            #print('ok, got ', h.GetName())
            tmp = h.GetName()
        except:
            print('ERROR getting histo {}{}!'.format(rdir,hname))
            continue

        #print('Pushing ', ich, hname)
        hs.append(h)

        canname = 'WCTEJuly2023_Quick2D_{}_{}'.format(ftag, hname)
        canname = canname.replace('_list_root','').replace('_ntuple','').replace('.root','')
        cw = 1100
        ch = 800
        if 'xymap' in hname:
            cw = 1000
            ch = 900
        can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
        cans.append(can)
        #can.Divide(8,4)
        #h.Rebin2D(2,2)
        opt = 'colz'
        is2d = True
        h.SetStats(0)
        h.Draw(opt)
        if 'minus' in hname or 'Xmap' in hname:
            h.SetFillColor(ROOT.kMagenta)
            h.SetFillStyle(1111)
            opt = 'hist'
            is2d = False
            #h.SetStats(1)
            mean = h.GetMean()
            emean = h.GetMeanError()
            print(mean, emean)
            rtxt = ROOT.TLatex(0.72, 0.85, '#mu={:1.2f}#pm{:1.2f}'.format(mean, emean))
            rtxt.SetNDC()
            rtxt.SetTextSize(0.04)
            rtxt.Draw()
            stuff.append(rtxt)
            if 'Xmap' in hname:
                meanXmap[h.GetName()] = [mean, emean]
        if is2d:
            rho = h.GetCorrelationFactor()
            rtxt = ROOT.TLatex(0.76, 0.85, '#rho={:1.2f}'.format(rho))
            rtxt.SetNDC()
            rtxt.SetTextSize(0.04)
            rtxt.Draw()
            stuff.append(rtxt)

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
        x1,x2 = h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax()
        y1,y2 = h.GetYaxis().GetXmin(),h.GetYaxis().GetXmax()
        if ('LC' in hname and 'RC' in hname) or ( 'TrigScint0' in hname and 'TrigScint1' in hname) :
            diag = ROOT.TLine(x1, y1, x2, y2)
            diag.SetLineColor(ROOT.kMagenta)
            diag.SetLineStyle(1)
            diag.SetLineWidth(2)
            diag.Draw()
            stuff.append(diag)
        if 'xymap' in hname:
            hl = ROOT.TLine(x1, 0.5*(y1+y2), x2, 0.5*(y1+y2))
            hl.SetLineWidth(2)
            hl.SetLineColor(ROOT.kMagenta)
            hl.Draw()
            vl = ROOT.TLine(0.5*(x1+x2), y1, 0.5*(x1+x2), y2)
            vl.SetLineWidth(2)
            vl.SetLineColor(ROOT.kMagenta)
            vl.Draw()
            stuff.append([hl,vl])
            if 'Cweighted_xymap' in hname and not 'p-like' in hname:
                meanX = h.GetMean(1)
                emeanX = h.GetMeanError(1)
                meanY = h.GetMean(2)
                emeanY = h.GetMeanError(2)
                meanXmap[h.GetName()] = [meanX, emeanX]

        
##################################
#       plots all the canvas     #
##################################

    """
    srun = ''
    tokens = filename.split('_')
    momentum = None
    runindex = filename.index('run')
    srun = filename[runindex+6:runindex+9]
    if momentum == None:
        momentum = getMomentum(srun)
    if momentum == None:
        momentum = getMergedMomentum(srun)

    pnote = makeMomentumLabel(srun, momentum)
    stuff.append(pnote)
    """

    cw = 1100
    ch = 800
    canname = 'meanXmap'
    can = ROOT.TCanvas(canname, canname, 100, 100, cw, ch)
    cans.append(can)
    #grMeanXmap = ROOT.TGraphErrors()
    ip = -1
    nb = len(meanXmap)
    name = 'meanXmap'
    title = ';;#DeltaX [quasi mm]'
    hmean = ROOT.TH1D(name, title, nb, 0, nb)
    for key,vals in meanXmap.items():
        ip = ip+1
        mean = vals[0]
        emean = vals[1]
        #grMeanXmap.SetPoint(ip, ip, mean)
        #grMeanXmap.SetPointError(ip, ip, emean)
        hmean.SetBinContent(ip+1, mean)
        hmean.SetBinError(ip+1, emean)
        prettykey = key.replace('hRef_','').replace('_nonZero','').replace('_Cweighted_xymap','').replace('Xmap','').replace('TrigScint','TS').replace('act','ACT')
        hmean.SetStats(0)
        hmean.GetXaxis().SetBinLabel(ip+1, prettykey)

    can.cd()
    hmean.SetMarkerSize(2)
    hmean.SetMarkerStyle(20)
    hmean.SetMarkerColor(ROOT.kBlue)
    hmean.SetMinimum(-15.)
    hmean.SetMaximum(15.)
    hmean.Draw("P")
    onel = makeLine(0., 0., nb, 0.)
    onel.Draw()
    stuff.append(onel)
    ROOT.gPad.SetGridx(1)
    ROOT.gPad.SetGridy(1)
    
    for can in cans:
        can.cd()
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

