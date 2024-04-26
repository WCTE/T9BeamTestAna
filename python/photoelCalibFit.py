#!/snap/bin/pyroot

# /usr/bin/python3

# jk
# 20/09/2022
# 14.7.2023

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from labelTools import *

cans = []
stuff = []
lines = []

##########################################

def  getParEst(h, x1, x2):
    sigma = 0
    A = -1.
    xmax  = -1.
    sumsq = 0.
    mean = 0.
    sumw = 0
    i1 = h.GetXaxis().FindBin(x1)
    i2 = h.GetXaxis().FindBin(x2)
    for i in range(i1,i2):
        x = h.GetBinCenter(i)
        y = h.GetBinContent(i)
        sumw = sumw + y
        mean = mean + y*x
        sumsq = sumsq + y*x*x
        if y > A:
            A = 1.*y
            xmax = 1.*x
    if sumw > 0:
        mean = mean / sumw
        sumsq = sumsq / sumw
    sigma = sumsq - mean*mean
    if sigma > 0:
        sigma = sqrt(sigma)
    return A, xmax, sigma    

##########################################

def makeLine(x1, y1, x2, y2, c = ROOT.kGreen, lst = 1, lw = 2):
    line = ROOT.TLine(x1, y1, x2, y2)
    line.SetLineColor(c)
    line.SetLineStyle(lst)
    line.SetLineWidth(lw)
    line.Draw()
    return line

##########################################

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

    ChNames = ChNamesCharged
    
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


    ROOT.gStyle.SetPalette(ROOT.kSolar)

    basedir = 'General/'
    
    #filename = 'output_300n_plots.root'
    filename = argv[1]
    print(f"Opening file '{filename}'")
    rfile = ROOT.TFile(filename, 'read')
    hbasenames = {
        #'hRef_Time' : ROOT.kGreen,
        'hRef_Charge' : ROOT.kCyan,
        #'hRef_Voltage' : ROOT.kMagenta,
        #'hRef_nPeaksC' : ROOT.kYellow,
        #'hRef_nPeaksA' : ROOT.kGray,
    }
    
    nChannels = 19 # 32
    Hs = []
    Txts = []
  

    ftag = filename.split('/')[-1].replace('output_','').replace('_plots.root','')

    os.system('mkdir -p pdf png')
    
    for hbasename in hbasenames:

        hs = []
        txts = []
           
        for ich in range(0, nChannels):
            # hack just for old tofs
            #if not ( ich >= 8 and ich <= 15):
            #    continue
            hname = hbasename + str(ich)
            h = rfile.Get(basedir + hname)
            try:
                #print('ok, got ', h.GetName())
                tmp = h.GetName()
            except:
                print('ERROR getting histo {}!'.format(hname))
                continue

            #print('Pushing ', ich, hname)
            hs.append(h)
        Hs.append(hs)

        canname = 'WCTEJuly2023_Quick1D_{}_{}'.format(ftag, hbasename)
        canname = canname.replace('_list_root','').replace('_ntuple','')
        #can = ROOT.TCanvas(canname, canname, 0, 0, 1600, 800)
        #cans.append(can)
        #can.Divide(8,4)

        relCsDict = {}
        relCs = []
        for h in hs:
            try:
                #print('ok, got ', h.GetName())
                tmp = h.GetName()
            except:
                print('ERROR getting histo!')
                continue
            ich = hs.index(h)

            if ich % 8 == 0:
                idigi = ich / 8
                off = 60
                cw = 4*400 + 4*off
                ch = 2*400
                if idigi > 1:
                    cw = 2*400 + 2*off
                    ch = 400
                can = ROOT.TCanvas(canname + f'_digi{idigi}', canname + f'_digi{idigi}', 0, 0, cw, ch)
                if idigi < 2:
                    can.Divide(4,2)
                else:
                    can.Divide(3,1)
                cans.append(can)

            
            can.cd(ich % 8 + 1)
            h.SetStats(1)
            #if not 'Time' in h.GetName():
            #if 'nPeaks' in h.GetName():
            ROOT.gPad.SetLogy(1)
            if 'TOF' in h.GetName():
                h.Rebin(2)
            
            #h.GetYaxis().SetRangeUser(1.e-4, h.GetYaxis().GetXmax())
            h.SetFillColor(hbasenames[hbasename])
            h.SetFillStyle(1111)
            chname = ChNames[hs.index(h)]
            h.SetTitle(chname)
            
            h.Rebin(2)
            h.GetXaxis().SetTitle('Charge [uncalibrated p.e.]')
            h.Draw('hist')

            x1,x2 = 0.5, 2.5
            if ich == 10 or ich == 11:
                x1, x2 = 1.1, 3.5
            if ich == 1:
                x1, x2 = 1.5, 3.
            if ich == 16 or ich == 17:
                x1, x2 = 2., 6.
            if ich == 18:
                x1, x2 = 1.7, 5.

            # CLOSURE range HACK
            x1,x2 = 0.5, 1.5
            
            fname = 'fit_' + h.GetName()
            fun = ROOT.TF1(fname, "[0]*exp(-(x-[1])^2/(2*[2]^2))", x1, x2)
            A, mu, sigma = getParEst(h, x1, x2)
            if ich == 1:
                sigma = 1.
            fun.SetParameters(A, mu, sigma)
            stuff.append(fun)
            h.Fit(fname, "", "0", x1, x2)
            pars = []
            pars = []
            for i in range(0, fun.GetNpar()):
                pars.append(fun.GetParameter(i))
            sf = 1.5
            h.Fit(fname, "", "0", pars[1] - sf*pars[2], pars[1] + sf*pars[2])
            gain = fun.GetParameter(1)
            gainErr = fun.GetParError(1)
            y1, y2 = h.GetMinimum(), h.GetMaximum()
            lines = [ makeLine(gain, y1, gain, y2, ROOT.kRed, 1, 2),
                      makeLine(gain - gainErr, y1, gain - gainErr, y2, ROOT.kRed, 2, 2),
                      makeLine(gain + gainErr, y1, gain + gainErr, y2, ROOT.kRed, 2, 2),
                      makeLine(1., y1, 1., y2, ROOT.kBlue, 2, 2) ]
            stuff.append(lines)

            relCsDict[ChNames[hs.index(h)]] = gain
            relCs.append(gain)
            
            fun.Draw('same')

            ROOT.gPad.Update()
        oldCs = [  0.30446, 0.6183148, 0.2150617, 0.316799, 0.2275831, 0.171613, 0.2256, 0.295012,
                   0.04457263158, 0.049545, 0.04662082192, 0.03502604651, 0.03854842105, 0.06932761905, 0.06069396226, 0.06920426087,
                   0.03495, 0.04052, 1.]
        print(relCsDict)
        newCs = []
        for oldC,relC in zip(oldCs,relCs):
            newCs.append(oldC*relC)
        print('OLD CALIBRATION CONSTANTS:')
        print(oldCs)
        print('NEW CALIBRATION CONSTANTS:')
        print(newCs)
        for oC,nC in zip(oldCs,newCs):
            print(f'oldC: {oC:1.3f}  newC: {nC:1.3f}')
        
##################################
#       plots all the canvas     #
##################################

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

    pnote = makeMomentumLabel(srun, momentum)
    stuff.append(pnote)
    
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

