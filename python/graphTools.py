#!/usr/bin/python
####################################################################################

import ROOT
from ctypes import c_double

####################################################################################

def makeLine(x1, x2, y1, y2):
    line = ROOT.TLine(x1, y1, x2, y2)
    line.SetLineColor(ROOT.kGreen)
    line.SetLineWidth(2)
    line.Draw()
    return line


####################################################################################
def MakeGraphNoErrs(xs, ys, col = ROOT.kBlack, mst = 20, ms = 1.):
    gr = ROOT.TGraph()
    for i in range(0,len(xs)):
        gr.SetPoint(i, xs[i], ys[i])
    gr.SetMarkerColor(col)
    gr.SetLineColor(col)
    gr.SetMarkerStyle(mst)
    gr.SetMarkerSize(ms)
    return gr

####################################################################################
def MakeGraph(xs, exs, ys, eys, col = ROOT.kBlack, mst = 20, ms = 1.):
    gr = ROOT.TGraphErrors()
    for i in range(0,len(xs)):
        gr.SetPoint(i, xs[i], ys[i])
        gr.SetPointError(i, exs[i], eys[i])
    gr.SetMarkerColor(col)
    gr.SetLineColor(col)
    gr.SetMarkerStyle(mst)
    gr.SetMarkerSize(ms)
    return gr

####################################################################################
def getGraphsChi2(gr, grfit):
    chi2 = 0.
    npts = 0
    x = c_double(0.)
    y = c_double(0.)
    xf = c_double(0.)
    yf = c_double(0.)
    for ip  in range(0, min(gr.GetN(), grfit.GetN())):
        gr.GetPoint(ip, x, y)
        grfit.GetPoint(ip, xf, yf)
        yerr = gr.GetErrorY(ip)
        if abs(x.value - xf.value) < 0.05 and yerr > 0:
            chi2 = chi2 + pow( (y.value-yf.value)/yerr, 2)
            npts = npts + 1
    return chi2, npts

        


####################################################################################
