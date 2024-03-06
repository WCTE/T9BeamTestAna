#!/snap/bin/pyroot

#/usr/bin/python3

# jk
# 20/09/2022
# 14.7.2023
# 5.3.2024
# 6.3.2024

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from labelTools import *
from tofUtil import *

cans = []
stuff = []
lines = []

##########################################

def makeLine(x1, x2, y1, y2):
    line = ROOT.TLine(x1, y1, x2, y2)
    line.SetLineColor(ROOT.kGreen)
    line.SetLineWidth(2)
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

def MakeGraph(xs, exs, ys, eys, col = ROOT.kBlack, mst = 20, ms = 1.):
    gr = ROOT.TGraphErrors()
    for i in range(0,len(xs)):
        gr.SetPoint(i, xs[i], ys[i])
        gr.SetPointError(i, exs[i], eys[i])
    gr.SetMarkerColor(col)
    gr.SetMarkerStyle(mst)
    gr.SetMarkerSize(ms)
    return gr
    

##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    pngdir = 'png_results/'
    pdfdir = 'pdf_results/'
    os.system(f'mkdir {pngdir}')
    os.system(f'mkdir {pdfdir}')

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


    dirname = 'histos/'
    filenames = [ 'ntuple_000396_plots.root ',
                  'ntuple_000394_plots.root ',
                  'ntuple_000393_plots.root ',
                  'ntuple_000392_plots.root ',
                  'ntuple_000398_plots.root ',
                  'ntuple_000399_plots.root ',
                  'ntuple_000402_plots.root'
        

                 ]

    
    os.system('mkdir -p pdf png')
    hs = []
    momenta = []
    hscp = []
    projYs = []
    txts = []
    #basedir = 'TrigScint/'
    pbasedir = 'TrigScint_p/'
    hname = 'hRef_TOF_TrigScintC_p-like'
     #'hRef_TOFPbA',
     #'hRef_TOFPbC',
     
     
     #'hRef_TOF_TrigScint0C',
     #'hRef_TOF_TrigScint1C',
     
     #'hRef_TOF_TrigScint001C',
     #'hRef_TOF_TrigScint023C',
     #'hRef_TOF_TrigScint101C',
     #'hRef_TOF_TrigScint123C',
     
     #'hRef_pbC_TrigScintC',
     #'hRef_pbC_TrigScintA',
     #'hRef_pbA_TrigScintC',
     #'hRef_pbA_TrigScintA',



    can = None
    opt = ''
    rfiles = []
    ys = []
    eys = []
    betas = []
    betagammas = []
    ebetas = []
    leg = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0)
    stuff.append(leg)
    cols = [ROOT.kBlack, ROOT.kGray+1, ROOT.kGreen+2, ROOT.kBlue, ROOT.kViolet, ROOT.kRed, ROOT.kOrange+1, ROOT.kYellow]
    for filename,col in zip(filenames,cols):
        rfile = ROOT.TFile(dirname + filename, 'read')
        rfiles.append(rfile)
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
        momenta.append(momentum)
        Hs = []
        Txts = []
        ftag = filename.split('/')[-1].replace('output_','').replace('_plots.root','')

        h = rfile.Get(pbasedir + hname)
        try:
            print('ok, got {} from file {}'.format(h.GetName(), rfile.GetName()))
            tmp = h.GetName()
        except:
            print('ERROR getting histo {}{} from file {}!'.format(pbasedir,hname, rfile.GetName()))
            continue

        hs.append(h)
        if can == None:
            canname = 'WCTEJuly2023_protons2D_{}_{}'.format(ftag, hname)
            canname = canname.replace('_list_root','').replace('_ntuple','')
            can = ROOT.TCanvas(canname, canname, 0, 0, 1100, 800)
            cans.append(can)
            #can.Divide(8,4)
        h.SetStats(0)
        #h.Rebin2D(2,2)
        hcp = h.DrawCopy('scat' + opt)
        if opt == '':
            hcp.GetXaxis().SetRangeUser(13, 23)
            hcp.GetYaxis().SetRangeUser(4, 16)
        alpha = 0.1
        hcp.SetMarkerColorAlpha(col, alpha)
        hcp.SetLineColorAlpha(col, alpha)
        hcp.SetMarkerSize(0.1)
        hcp.SetMarkerStyle(6)
        ROOT.gPad.Update()
        h.SetFillColorAlpha(hcp.GetMarkerColor(),1)
        h.SetFillStyle(1111)
        h.SetMarkerSize(2)
        h.SetMarkerStyle(20)
        hscp.append(hcp)
        beta = getBeta(ms['p'], momentum)
        betagamma = getBetaGamma(ms['p'], momentum)
        leg.AddEntry(h, 'p = {:4} MeV/c #beta={:1.2f}'.format(str(momentum), beta), 'F')

        ys.append(h.GetMean(2))
        eys.append(h.GetMeanError(2))
        betas.append(beta)
        betagammas.append( betagamma)
        ebetas.append(0.)

        projY = h.ProjectionY()
        projYs.append(projY)
        
        opt = 'same'
        #adjustStats(h)
        #ROOT.gPad.Update()

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

    ##########################################
    canname = 'protonsLandau'
    canp = ROOT.TCanvas(canname, canname, 200, 200, 1000, 800)
    canp.cd()
    opt = ''
    legp = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
    legp.SetBorderSize(0)
    stuff.append(legp)
    for h,projY,momemntum in zip(hs,projYs,momenta):
        projY.SetLineColor(h.GetFillColor())
        projY.SetFillColorAlpha(h.GetFillColor(), 0.3)
        projY.SetLineWidth(2)
        projY.SetLineStyle(1)
        projY.Rebin(2)
        val = projY.Integral(0, projY.GetXaxis().GetNbins()+1)
        if val > 0.:
            projY.Scale(1./val)
        projY.SetMaximum(projY.GetMaximum()*1.2)
        projY.SetStats(0)
        projY.Draw('hist' + opt)
        if opt == '':
            projY.GetXaxis().SetRangeUser(4., 14.)
        legp.AddEntry(projY, 'p = {:4} MeV/c #beta={:1.2f}'.format(str(momentum), beta), 'F')
        opt = 'same'
    legp.Draw()
    canp.Update()
    cans.append(canp)

    
    ##########################################
    canname = 'protonsBetaGraph'
    gcan = ROOT.TCanvas(canname, canname, 100, 100, 1000, 400)
    gcan.Divide(2,1)
    print(betas, ebetas, ys, eys)
    grb = MakeGraph(betas, ebetas, ys, eys)
    grbg = MakeGraph(betagammas, ebetas, ys, eys)
    
    hn = 'tmpbg'
    ht = hn + ';#beta#gamma;Mean trig. scint. charge [a.u.];'
    htmpbg = ROOT.TH2D(hn, ht, 100, 0.5, 1.5, 100, 6., 12.)
    htmpbg.SetStats(0)
    htmpbg.GetXaxis().SetMoreLogLabels()

    gcan.cd(1)
    htmpbg.Draw()
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gPad.SetLogx()
    ROOT.gPad.SetGridx(1)
    ROOT.gPad.SetGridy(1)
    grbg.Draw("P")

    hn = 'tmpbg'
    ht = hn + ';#beta;Mean trig. scint. charge [a.u.];'
    htmpb = ROOT.TH2D(hn, ht, 100, 0.5, 0.85, 100, 6., 12.)
    htmpb.SetStats(0)
    htmpb.GetXaxis().SetMoreLogLabels()

    gcan.cd(2)
    htmpb.Draw()
    ROOT.gStyle.SetOptTitle(0)
    #ROOT.gPad.SetLogx()
    ROOT.gPad.SetGridx(1)
    ROOT.gPad.SetGridy(1)
    grb.Draw("P")
    #fun = ROOT.TF1('fun', '[0]/x^2 + [1]', 0.1, 1.)
    #fun.SetParameters(0.1, 1.)
    fun = ROOT.TF1('fun', '[0]/x^2*(log([1]*x/sqrt(1-x*x)) - x^2) + [2]', 0.1, 1.)
    fun.SetParameters(10., 1., 0.5)
    grb.Fit('fun')
    
    cans.append(gcan)
    stuff.append([grb, grbg])
    gcan.Update()


    
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

