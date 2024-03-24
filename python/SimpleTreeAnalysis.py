#!/snap/bin/pyroot
# was: #!/usr/bin/python3
# jk Út 30. ledna 2024, 14:16:31 CET

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

##########################################
def InitTrees(rfile, trs):
    trees = {trn : rfile.Get(trn) for trn in trs}
    return trees

##########################################
def ReadEntry(trs, i):
    for trn, tr in trs.items():
        tr.GetEntry(i)

##########################################
def main(argv):

    fname = "output/ntuple_000409.root"
    rfile = ROOT.TFile(fname, 'read')

    treenames = [#"TOF00", "TOF01", "TOF02", "TOF03",
                 #"TOF10", "TOF11", "TOF12", "TOF13",
                 "ACT2L", "ACT2R", "ACT3L", "ACT3R",
                 "PbGlass"
                 ]
    refname = treenames[-1]

    print('Will analyze the following trees:')
    print(treenames)
    
    nbins = 50
    kPeakVoltage = "PeakVoltage"
    kSignalTime = "SignalTime"
    hbasenames = {  kPeakVoltage: [0., 2., ROOT.kYellow] } #, kSignalTime : [-150., 150., ROOT.kCyan] }
    histos = {}
    for treename in treenames:
        histos[treename] = {}
        for hbase,bins in hbasenames.items():
            hname = f'{treename}_{hbase}'
            htitle = hname + f';{hbase};events;'
            histos[treename][hbase] = ROOT.TH1D(hname, htitle, nbins, bins[0], bins[1])
            histos[treename][hbase].SetFillColor(bins[2])
            histos[treename][hbase].SetFillStyle(1111)


    histos2d = {}
    for i in range(0,len(treenames)):
        for j in range(0,len(treenames)):
            if j >= i: continue
            for hbase,bins in hbasenames.items():
                h2name = f'{treenames[j]}_vs_{treenames[i]}_{hbase}'
                htitle = h2name + f';{treenames[i]};{treenames[j]};'
                histos2d[h2name] = ROOT.TH2D(h2name, htitle, nbins, bins[0], bins[1], nbins, bins[0], bins[1])
                histos2d[h2name].SetStats(0)
            hname = f'{treename}_{hbase}'
            
    trees = InitTrees(rfile, treenames)
    Nevt = trees[refname].GetEntries()
    # HACK!!!
    #Nevt = int (Nevt / 100.)
    for ievt in range(0,Nevt):
        if ievt % 1000 == 0: print(f'processing {ievt}/{Nevt}')
        ReadEntry(trees, ievt)
        for treename in treenames:
            if eval('trees["{}"].nPeaks'.format(treename)) > 0: 
                for hbase,h in histos[treename].items():
                    # look at the zero's peak
                    h.Fill( eval('trees["{}"].{}[0]'.format(treename,hbase)) )
                    # alt:
                    #if hbase == kSignalTime:
                    #    h.Fill( trees[treename].SignalTime[0])
                    #elif hbase == kPeakVoltage:
                    #    h.Fill( trees[treename].PeakVoltage[0])
                        
        for i in range(0,len(treenames)):
            for j in range(0,len(treenames)):
                if j >= i: continue
                # look at the zero's peak
                if eval('trees["{}"].nPeaks'.format(treenames[i])) > 0 and eval('trees["{}"].nPeaks'.format(treenames[j])) > 0:
                    for hbase,bins in hbasenames.items():
                        h2name = f'{treenames[j]}_vs_{treenames[i]}_{hbase}'
                        histos2d[h2name].Fill(eval('trees["{}"].{}[0]'.format(treenames[i],hbase)), eval('trees["{}"].{}[0]'.format(treenames[j],hbase)))
                        # alt:
                        #if hbase == kSignalTime:
                        #    x = trees[treenames[i]].SignalTime[0]
                        #    y = trees[treenames[j]].SignalTime[0]
                        #    histos2d[h2name].Fill(x,y)
                        #elif hbase == kPeakVoltage:
                        #    x = trees[treenames[i]].PeakVoltage[0]
                        #    y = trees[treenames[j]].PeakVoltage[0]
                        #    histos2d[h2name].Fill(x,y)                    
                    
    for treename in treenames:
        canname = f'can_{treename}'
        print(f'Making canvas {canname}')
        can = ROOT.TCanvas(canname, canname)
        can.Divide(len(hbasenames),1)
        cans.append(can)
        i = 0
        for hbase,h in histos[treename].items():
            i = i + 1
            can.cd(i)
            h.Draw('hist')
            ROOT.gPad.Update()
        can.Update()

    canname = f'can2d'
    can2d = ROOT.TCanvas(canname, canname)
    can2d.Divide(4,4)
    cans.append(can2d)
    ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)
    i = 0
    for hname2d in histos2d:
        i = i + 1
        can2d.cd(i)
        histos2d[hname2d].Draw('colz')
        ROOT.gPad.SetLogz(1)
        ROOT.gPad.Update()
        
    for can in cans:
        try:
            can.Print(can.GetName() + '.png')
            can.Print(can.GetName() + '.pdf')
        except:
            print('Error')
        
    stuff.append(histos)
    #ROOT.gApplication.Run()
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

