#!/usr/bin/python
####################################################################################

import ROOT

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
