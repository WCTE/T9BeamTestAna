#!/usr/bin/python

#/snap/bin/pyroot

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
    hnames2d = [
        'hRef_pbC_act23C_nonp-like',
        'hRef_ACT1CACT23C_nonp-like',
                ]
    
    chbasedir = 'Charged_nonp/'

    for hname in hnames2d:
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
        can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
        cans.append(can)
        #can.Divide(8,4)
        #h.Rebin2D(2,2)
        opt = 'colz'
        is2d = True
        h.SetStats(0)
        h.Draw(opt)
        ROOT.gPad.SetLogz(1)

        if 'hRef_pbC_act23C_nonp-like' in h.GetName():
            h.GetXaxis().SetRangeUser(0., 20)
            h.GetYaxis().SetRangeUser(0., 60)            
        if 'hRef_ACT1CACT23C_nonp-like' in h.GetName():
            h.GetXaxis().SetRangeUser(0., 40)
            h.GetYaxis().SetRangeUser(0., 50)

            
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
    for can in cans:
        can.cd()
        if 'vs' in can.GetName():
            pnote.Draw()            
        can.Update()
        can.Print(pngdir + can.GetName() + '.png')
        can.Print(pdfdir + can.GetName() + '.pdf')
    
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

