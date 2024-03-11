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

from labelTools import *

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

    x1,x2 = -50, 50
    y1,y2 = -50, 50
    name = 'tmp'
    title = name + ';X [mm];Y[mm]'
    h = ROOT.TH2D(name, title, 100, x1, x2, 100, y1, y2)
    h.SetStats(0)
    ROOT.gStyle.SetOptTitle(0)

    canname = 'emptyXY'
    canname = canname.replace('_list_root','').replace('_ntuple','')
    cw = 1000
    ch = 1000

    can = ROOT.TCanvas(canname, canname, 0, 0, cw, ch)
    cans.append(can)
        
    h.Draw()
    
    x1,x2 = h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax()
    y1,y2 = h.GetYaxis().GetXmin(),h.GetYaxis().GetXmax()
        
    hl = ROOT.TLine(x1, 0.5*(y1+y2), x2, 0.5*(y1+y2))
    hl.SetLineWidth(1)
    hl.SetLineColor(ROOT.kMagenta)
    hl.Draw()
    vl = ROOT.TLine(0.5*(x1+x2), y1, 0.5*(x1+x2), y2)
    vl.SetLineWidth(1)
    vl.SetLineColor(ROOT.kMagenta)
    vl.Draw()
    stuff.append([hl,vl])

    diag = ROOT.TArrow(x1, y1, x2, y2, 0.02, '<|-|>')
    diag.SetLineColor(ROOT.kMagenta)
    diag.SetFillColor(ROOT.kMagenta)
    diag.SetLineStyle(1)
    diag.SetLineWidth(4)
    diag.Draw()
    adiag = ROOT.TArrow(x1, y2, x2, y1, 0.02, '<|-|>')
    adiag.SetLineColor(ROOT.kMagenta)
    adiag.SetFillColor(ROOT.kMagenta)
    adiag.SetLineStyle(1)
    adiag.SetLineWidth(4)
    adiag.Draw()
    stuff.append(adiag)

    
        
    ##################################
    #       plots all the canvas     #
    ##################################
    

    for can in cans:
        can.cd()
        if 'vs' in can.GetName():
            pnote.Draw()            
        can.Update()
        can.Print(can.GetName() + '.png')
        #can.Print(can.GetName() + '.pdf')
    
        
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

