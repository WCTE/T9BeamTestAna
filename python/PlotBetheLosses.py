#!/usr/bin/python

# jk 26.4.2018
# jk  3.4.2024

from Losses import *
#from Brems import *
from graphTools import *
from tofUtil import *

print('*** Supported particles:')
PrintParticles()
print('*** Supported materials:')
PrintMaterials()

print('***************************************************')
ROOT.gStyle.SetPadLeftMargin(0.2)

# material: dX in cm
myMat = {'Polystyrene' : 1.,
          'Al' : 1.e-2
         }

myParts = {
    #'Electron',
    'Positron' : pcols['e'],
    'Muon' : pcols['mu'],
    'Pion' : pcols['pi'],
    'Proton' : pcols['p'],
    'Deuteron' : pcols['D'],
}

# MeV/c:
p1 = 200.
p2 = 1200.
N = 200
dp = (p2-p1) / N
ps = [ p1 + i*dp for i in range(0,N+1) ]

GrsdE = {}
GrsdEdX = {}
GrsdBeta = {}
GrsdP= {}

for mname,dX in myMat.items():

    GrsdE[mname] = {}
    GrsdEdX[mname] = {}
    GrsdBeta[mname] = {}
    GrsdP[mname] = {}
    

    
    for pname,col in myParts.items():
    
        particle = gParticles[pname]
        material = gMaterials[mname]

        M = particle.GetM()

        Ps = []
        dEs = []
        dEdXs = []
        dPs = []
        dBetas = []
        
        useHigherCorrs = False
        # MIP:

        for p in ps:
            E = sqrt(p*p + M*M)
            T = E - M
            # beta*gamma:
            bg = p/M
            beta = GetBetaFromBg(bg)
            gamma = GetGammaFromBg(bg)

            #pname = particle.GetName()
            print('*** {:} in {:} using dX={:f}cm'.format(pname, material.GetName(), dX))
            # print bg, beta, gamma
            print('    BEFORE:    T={:3.1f} MeV p={:3.1f} MeV, E={:3.1f}, beta={:1.4f}, gamma = {:3.3f}, beta*gamma={:3.3f}'.format(T, p, E, beta, gamma, beta*gamma))
            #print(beta, particle, material, useHigherCorrs)
            dedx, halflog = dEdX(beta, particle, material, useHigherCorrs)
            dE = dedx*dX
            newE = E - dE

            if newE*newE > M*M:
                newp = sqrt(newE*newE - M*M)
                newT = newE - M
                fracEloss = dE / E
                newgamma = newE/M
                newbeta = newp / newE


                dP = p - newp
                dEs.append(dE)
                dEdXs.append(dedx)
                dPs.append(dP)
                dBetas.append(beta - newbeta)
                Ps.append(p)

                print('    AFTER dX={:}cm: T={:3.1f} MeV p={:3.1f} MeV, E={:3.1f}, beta={:1.4f}, gamma = {:3.3f}, beta*gamma={:3.3f}'.format(dX, newT, newp, newE, newbeta, newgamma, newbeta*newgamma))
                print('    Ionization losses                                      : {:1.3f} MeV/cm'.format(dedx))
                print('    Ionization losses                                      : {:1.3f} MeV'.format(dE))
                #print('    New momentum after 1cm of the material                 : {:3.1f}'.format(newp))
                #print('    New energy after 1cm of the material                   : {:3.1f}'.format(newE))
                #print('    New kinetic energy after 1cm of the material           : {:3.1f}'.format(newT))

                print('    Energy fraction loss after 1cm of the material         : {:1.4f}'.format(dedx / E))
                print('    Kinetic energy fraction loss after 1cm of the material : {:1.4f}'.format(dedx / T))

                print('    New momentum fraction after 1cm of the material        : {:1.3f}'.format(newp/p))
                print('    New energy fraction after 1cm of the material          : {:1.3f}'.format(newE/E))
                print('    New kinetic energy fraction after 1cm of the material  : {:1.3f}'.format(newT/T))
                print('    New beta, deltaBeta, new gamma                         : {:1.3f}, {:1.3f}, {:1.3f}'.format(newbeta, beta-newbeta, newgamma))


                #E = gamma*M
                #print('    Radiation losses  : {:1.3f} MeV/cm'.format( dEdXBrems(E, particle, material), ))
            else:
                print('new E is negative!')

        #print(ps, dEs)
        GrsdE[mname][pname]     = MakeGraphNoErrs(Ps, dEs, col)
        GrsdEdX[mname][pname]   = MakeGraphNoErrs(Ps, dEdXs, col)
        GrsdP[mname][pname]     = MakeGraphNoErrs(Ps, dPs, col)
        GrsdBeta[mname][pname]  = MakeGraphNoErrs(Ps, dBetas, col)

print(GrsdE)
cans = {}
stuff = []
legs = {}
imat = -1
N = 0
for mname,dX in myMat.items():
    imat = imat + 1
    canname = f'TheoryBethe_{mname}'
    cans[mname] = ROOT.TCanvas(canname, canname, 0, 0, 1200, 1000)
    cans[mname].Divide(2,2)
    opt = 'P'
    name = 'h2'
    title = ';p [MeV/c];'
    h2dE = ROOT.TH2D(name + 'dE', title + '#DeltaE [MeV]', 100, ps[0], ps[-1], 100, 0.1, 30.)
    h2dEdX = ROOT.TH2D(name + 'dEdX', title + 'dE/dX [MeV/cm]', 100, ps[0], ps[-1], 100, 0.1, 20.)
    h2dP = ROOT.TH2D(name + 'dP', title + '#Deltap [MeV/c]', 100, ps[0], ps[-1], 100, 0.1, 20.)
    h2dBeta = ROOT.TH2D(name + 'dBeta', title + '#Delta#beta', 100, ps[0], ps[-1], 100, 0.001, 0.05)
    hs = [h2dE, h2dEdX, h2dP, h2dBeta]
    N = len(hs)
    hcps = []
    for ih in range(0,len(hs)):
        cans[mname].cd(ih+1)
        hs[ih].SetStats(0)
        hcp = hs[ih].DrawCopy()
        hcps.append(hcp)
    ipart = -1
    for pname in myParts:
        ipart = ipart + 1
        if ipart == 0:
            legs[mname] = ROOT.TLegend(0.60, 0.55, 0.88, 0.88)
            legs[mname].SetHeader('{:} dX={:1.2f}cm'.format(mname,dX))
            legs[mname].SetBorderSize(0)
        legs[mname].AddEntry(GrsdE[mname][pname], pname, "L")
        print(pname, mname)
        grs = [GrsdE[mname][pname],
               GrsdEdX[mname][pname],
               GrsdP[mname][pname],
               GrsdBeta[mname][pname]
               ]
        for ih in range(0,len(grs)):
            cans[mname].cd(ih+1)
            grs[ih].SetLineWidth(2)
            grs[ih].Draw('L')
        
    for mname,can in cans.items():
        for ipad in range(0,N):
            can.cd(ipad+1)
            legs[mname].Draw()
            #ROOT.gPad.SetLogy(1)
            ROOT.gPad.SetGridy(1)
            ROOT.gPad.SetGridx(1)
            ROOT.gPad.Update()
            

    stuff.append([hs, hcps, legs])


# And just print all canvases;)
pngdir = 'png_results/'
pdfdir = 'pdf_results/'

for mname, can in cans.items():
    try:
        can.cd()
        can.Update()
        can.Print(pngdir + can.GetName() + '.png')
        can.Print(pdfdir + can.GetName() + '.pdf')
    except:
        print('ERROR printing canvas!')          

                
ROOT.gApplication.Run()
