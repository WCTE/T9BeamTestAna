#!/snap/bin/pyroot
# was: #!/usr/bin/python3
# Út 16. dubna 2024, 16:58:09 CEST

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []


####################################################################################
def MakeFitBox(fun, addtag, x1 = 0.15, y1 = 0.65):
    txt = ROOT.TLatex(x1, y1, addtag + '#alpha + #beta|p|, #alpha={:1.2f} #beta={:1.5f}'.format(fun.GetParameter(0), fun.GetParameter(1) ))
    txt.SetNDC()
    txt.SetTextColor(fun.GetLineColor())
    txt.Draw()
    return txt

####################################################################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
   

    canname = 'cutFitPars'
    can = ROOT.TCanvas(canname, canname)
    cans.append(can)

    stuff.append(can)

    # [momentum, b, c] where b, c are cut line intercepts of the y and x axes
    points = [
        [540., 16., 11.],
        [600., 18., 9.],
        [700., 13., 10.],
        [740., 17., 11.5],
        [800., 14., 9.],
        [860., 17.5, 12.5],
        [900., 15., 10.],
        [940., 12., 9.5],
        [1000., 15.5, 9.],
        [1060., 16.5, 9.],
        [1120., 16., 9.],
    ]

    gra = ROOT.TGraphErrors()
    grb = ROOT.TGraphErrors()
    grc = ROOT.TGraphErrors()

    for i in range(0, len(points)):
        p = points[i][0]
        b = points[i][1]
        c = points[i][2]
        a = -b/c
        gra.SetPoint(i, p, a)
        grb.SetPoint(i, p, b)
        grc.SetPoint(i, p, c)

        # manual guesstimate
        gra.SetPointError(i, 0., a*sqrt(2)*0.5)
        grb.SetPointError(i, 0., 0.5)
        grc.SetPointError(i, 0., 0.5)

        
    name = 'tmp'
    title = ';p [MeV/c];cut line fit pars b, c'
    ROOT.gStyle.SetOptTitle(0)
    p1 = 400.
    p2 = 1200.
    h2 = ROOT.TH2D(name, title, 100, p1, p2, 100, 0., 30.)
    h2.SetStats(0)
    h2.Draw()

    fb = ROOT.TF1('fb', '[0] + [1]*x', p1, p2)
    fb.SetParameters(15., 0.)
    fb.SetLineColor(ROOT.kBlue)
    fb.SetLineStyle(2)
    grb.Fit('fb')
    grb.Draw('PL')
    grb.SetMarkerColor(ROOT.kBlue)
    grb.SetLineColor(ROOT.kBlue)
    grb.SetMarkerStyle(20)
    grb.SetMarkerSize(1)

    fc = ROOT.TF1('fc', '[0] + [1]*x', p1, p2)
    fc.SetParameters(10., 0.)
    fc.SetLineColor(ROOT.kGreen+2)
    fc.SetLineStyle(2)
    grc.Fit('fc');
    grc.SetMarkerColor(ROOT.kGreen+2)
    grc.SetLineColor(ROOT.kGreen+2)
    grc.SetMarkerStyle(21)
    grc.SetMarkerSize(1)
    grc.Draw('PL')

    txtb = MakeFitBox(fb, 'b: ', 0.15, 0.82)
    txtc = MakeFitBox(fc, 'c: ', 0.15, 0.72)
    
    stuff.append([gra, grb, grc, txtb, txtc])

    can.Print(can.GetName() + '.png')
    can.Print(can.GetName() + '.pdf')
    
    ROOT.gPad.Update()
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

