#!/usr/bin/python

import ROOT

from data_runs import *
from tofUtil import *

ChNamesCharged =  {0: 'ACT-00', 1: 'ACT-01', 2: 'ACT-10', 3: 'ACT-11', 4: 'ACT-20', 5: 'ACT-21', 6: 'ACT-30', 7: 'ACT-31',
                   8: 'TOF-00', 9: 'TOF-01', 10: 'TOF-10', 11: 'TOF-11', 12: 'TOF-20', 13: 'TOF-21', 14: 'TOF-30', 15: 'TOF-31',
                   16: 'HC-00', 17: 'HC-10', 18: 'LG-00',
                   19: 'X', 20: 'X', 21: 'X', 22: 'X', 23: 'X',
                   24: 'X', 25: 'X', 26: 'X', 27: 'X', 28: 'X', 29: 'X', 30: 'X', 31: 'X' }


ChNamesHodo =  {    0 : "ACT0L",    1: "ACT0R",
                    2: "ACT1L",    3: "ACT1R",
                    4: "ACT3L",    5: "ACT3R",
                    6: "TriggerScint",
                    7: 'TOF-00', 8: 'TOF-01', 9: 'TOF-10', 10: 'TOF-11', 11: 'TOF-20', 12: 'TOF-21', 13: 'TOF-30', 14: 'TOF-31',
                    15: 'LG-00',
                    16: "HD8",   17: "HD9",   18: "HD10",   19: "HD11",   20: "HD12", 21: "HD13",   22: "HD14",
                    23: "HD0",   24: "HD1",   25: "HD2",  26: "HD3",  27: "HD4",  28: "HD5",  29: "HD6",  30: "HD7"
            }

##########################################
def makeLine(x1, y1, x2, y2, col = ROOT.kBlack, lst = 1, lw = 1):
    line = ROOT.TLine(x1, y1, x2, y2)
    line.SetLineColor(col)
    line.SetLineStyle(lst)
    line.SetLineWidth(lw)
    line.Draw()
    return line

##########################################
def makeLines(h, eoff, parts, momentum, useYaxis = False):
    lines = []
    #y1 = h.GetYaxis().GetXmin()
    #y2 = h.GetYaxis().GetXmax()
    print("LINE MOMENTUM: ", momentum)
    y1 = 1.05*h.GetMaximum()
    y2 = h.GetMinimum()
    if useYaxis:
        y1 = h.GetYaxis().GetXmin()
        y2 = h.GetYaxis().GetXmax()
    te = getTof(ms['e'], momentum) + eoff
    for part in parts:
        dt = getTofDiff('e', part, momentum)
        #print(f'makeLines {part}: dt={dt} ns')
        #print('line coors: ', te + dt, y1, te + dt, y2)
        line = makeLine(te + dt, y1, te + dt, y2, pcols[part], 1, 1)
        line.Draw()
        lines.append(line)
    return lines

##########################################
def makeFitLines(h, parts, times):
    lines = []
    #y1 = h.GetYaxis().GetXmin()
    #y2 = h.GetYaxis().GetXmax()
    #print("LINE MOMENTUM: ", momentum)
    y1 = 1.05*h.GetMaximum()
    y2 = h.GetMinimum()
    for part,time in zip(parts,times):
        #print(f'makeFitLines {part}: time={time} ns')
        #print('line coors: ', time, y1, time, y2)
        line = ROOT.TLine(time, y1, time, y2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.SetLineColor(pcols[part])
        line.Draw()
        lines.append(line)
    return lines

##########################################
def makeOneLine(h, value, part):
    y1 = h.GetMaximum()
    y2 = h.GetMinimum()
    line = ROOT.TLine(value, y1, value, y2)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.SetLineColor(pcols[part])
    line.Draw()
    return line


##########################################
def makePaperLabel(srun, momentum, x = 0.12, y = 0.86, size = 0.04, addn = True, addSlit = True, dx = 0.45):
    #momentum = getMomentum(srun)
    tag = '+'
    if momentum < 0:
        tag = '-'
    ctxt = f'WCTE Collaboration'
    txt = f'p={tag}{abs(momentum)} MeV/c'
    if addn:
        n = -1
        try:
            n = runsRefractionIndexDict[int(srun)]
        except:
            print('ERROR getting the ACT n info for run {}'.format(srun))
        if n > 0:
            txt = txt + f' n={n}'
    if addSlit:
        slit = ''
        slitperc = -1
        try:
            slitperc = runsSlitDict[int(srun)]
        except:
            print('ERROR getting the slit information for run {}'.format(srun))
        if slitperc > 0:
            slit = f' slit {slitperc}%'
            txt = txt + slit
    pnote = ROOT.TLatex(x + dx, y, txt)
    pnote.SetTextSize(size)
    pnote.SetNDC()
    cnote = ROOT.TLatex(x, y, ctxt)
    cnote.SetTextSize(size)
    cnote.SetNDC()
    return cnote, pnote

##########################################
def makeMomentumLabel(srun, momentum, x = 0.12, y = 0.86, size = 0.04, addn = True, addSlit = True):
    #momentum = getMomentum(srun)
    tag = '+'
    if momentum < 0:
        tag = '-'
    txt = f'WCTE TB2023 run {srun}, p={tag}{abs(momentum)} MeV/c'
    if addn:
        n = -1
        try:
            n = runsRefractionIndexDict[int(srun)]
        except:
            print('ERROR getting the ACT n info for run {}'.format(srun))
        if n > 0:
            txt = txt + f' n={n}'
    if addSlit:
        slit = ''
        slitperc = -1
        try:
            slitperc = runsSlitDict[int(srun)]
        except:
            print('ERROR getting the slit information for run {}'.format(srun))
        if slitperc > 0:
            slit = f' slit {slitperc}%'
            txt = txt + slit
    pnote = ROOT.TLatex(x, y, txt)
    pnote.SetTextSize(size)
    pnote.SetNDC()
    return pnote

##########################################
def adjustStats(h):
    ROOT.gPad.Update()
    st = h.GetListOfFunctions().FindObject("stats")
    st = ROOT.gPad.GetPrimitive("stats")
    st.SetX1NDC(0.7)
    st.SetX2NDC(0.9)
    st.SetY1NDC(0.65)
    st.SetY2NDC(0.9)
    
