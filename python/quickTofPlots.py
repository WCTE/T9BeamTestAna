#!/snap/bin/pyroot

#/usr/bin/python3

# jk
# 19.4.2024

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from labelTools import *
from tofUtil import *

cans = []
stuff = []
lines = []


def makeLine(x1, x2, y1, y2):
    line = ROOT.TLine(x1, y1, x2, y2)
    line.SetLineColor(ROOT.kGreen)
    line.SetLineWidth(2)
    line.Draw()
    return line


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
    gBatch = False
    if gBatch:
        ROOT.gROOT.SetBatch(1)

    if len(argv) < 2:
        PrintUsage(argv)
        return

    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetPalette(ROOT.kSolar)

    basedir = 'TOF/'
    
    filename = argv[1]
    rfile = ROOT.TFile(filename, 'read')

    srun = ''
    tokens = filename.split('_')

    momentum = None
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


    hs = []
    funs = []
    hnames = ['hTOFAllLow', 'hTOFElectronIDLow', 'hTOFMuonIDLow', 'hTOFPionIDLow']

    sameopt = ''
    opt = 'e1x0'
    cols = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]

    canname = f'WCTEJuly2023_QuickTOF_{srun}'
    canname = canname.replace('_list_root','').replace('_ntuple','')
    cw = 1000
    ch = 800
    can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
    cans.append(can)
    can.cd()
    ROOT.gPad.SetLogy(1)
    for col,hname in zip(cols,hnames):
        h = rfile.Get(basedir + hname)
        print('got {} I={:1.2f}'.format(h.GetName(), h.GetEntries()))
        h.SetStats(0)
        h.SetLineColor(col)
        h.SetMarkerColor(col)
        h.SetMarkerSize(1)
        h.SetMarkerStyle(20)
        h.SetLineWidth(2)
        h.Draw(opt + sameopt)
        if not 'All' in h.GetName():
            stddev = h.GetStdDev()
            sf = 2
            x1 = h.GetMean() - sf*stddev
            x2 = h.GetMean() + sf*stddev
            funname = 'gfit_{}'.format(h.GetName())
            fun = ROOT.TF1(funname, '[0]*exp(-(x-[1])^2/(2*[2]^2))', x1, x2)
            funs.append(fun)
            fun.SetParameters(h.GetMaximum(), h.GetMean(), stddev)
            h.Fit(funname, '0', '', x1, x2)
            fun.SetLineColor(col)
            fun.SetLineStyle(2)
            fun.Draw('same')
        sameopt = 'same'
        hs.append(h)



##################################
#       plots all the canvas     #
##################################

 
    pnote = makeMomentumLabel(srun, momentum, 0.12, 0.92)
    stuff.append(pnote)
    pnote.Draw()

    parts = ['e', 'mu', 'pi']
    eoff = 0.
    #te = getTof(ms['e'], momentum)
    #eoff = fit.GetParameter(1) - te

    lines = makeLines(hs[0], eoff, parts, momentum)
    stuff.append(lines)    

    fit_te = funs[0].GetParameter(1)
    fit_tmu = funs[1].GetParameter(1)
    fit_tpi = funs[2].GetParameter(1)

    fitlines = makeFitLines(hs[0], parts, [fit_te, fit_tmu, fit_tpi])
    stuff.append(fitlines)
    
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

