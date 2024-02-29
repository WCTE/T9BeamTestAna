#!/snap/bin/pyroot

# was: #!/usr/bin/python3
# Pá 14. července 2023, 18:39:38 CEST
# devel: II/2024

from data_runs import *
from tofUtil import *
from labelTools import *

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

import os

##########################################

def PrintUsage(argv):
    print('Usage:')
    print('{} histos/output_list_root_run_XYZ_plots.root [-b]'.format(argv[0]))
    return

##########################################
def makeLines(h, eoff, parts, momentum):
    lines = []
    #y1 = h.GetYaxis().GetXmin()
    #y2 = h.GetYaxis().GetXmax()
    print("LINE MOMENTUM: ", momentum)
    y1 = 1.05*h.GetMaximum()
    y2 = h.GetMinimum()
    te = getTof(ms['e'], momentum) + eoff
    for part in parts:
        dt = getTofDiff('e', part, momentum)
        print(f'makeLines {part}: dt={dt} ns')
        print('line coors: ', te + dt, y1, te + dt, y2)
        line = ROOT.TLine(te + dt, y1, te + dt, y2)
        line.SetLineStyle(1)
        line.SetLineWidth(1)
        line.SetLineColor(pcols[part])
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
        print(f'makeFitLines {part}: time={time} ns')
        print('line coors: ', time, y1, time, y2)
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
# cts     ...  central peak time of assumed gauss
# w       ... width
# t1, t2: ... fit window

##########################################
def Fit(h, tag, momentum, ct, w, t1, t2, peaksf = 1.):
    fname = 'fit{}'.format(tag)
    hname = h.GetName()
    fit = ROOT.TF1(fname, '[0]*exp(-(x-[1])^2/(2*[2]^2))', t1, t2)
    fit.SetLineColor(ROOT.kBlack)
    fit.SetLineStyle(2)
    #ampl = h.GetMaximum() / peaksf
    ampl = h.GetBinContent(h.FindBin(ct)) / peaksf
    print('Amplitude initially ', ampl)
    fit.SetParameters(ampl, ct, w)
    fit.SetParLimits(2, 0., 0.8)
    for ip in range(0, fit.GetNpar()):
        print(fit.GetParameter(ip))
    #prefit = fit.DrawCopy('same')
    #stuff.append(prefit)
    h.Fit(fname, '', '', t1, t2)
    mean = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    #print(f'1) {hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
    sf = 2.
    h.Fit(fname, '', '', mean - sf*sigma, mean + sf*sigma)
    mean = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    print(f'{hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
    fit.SetLineColor(h.GetLineColor())
    fit.SetLineStyle(2)
    fit.Draw('same')

    te = getTof(ms['e'], momentum)
    eoff = fit.GetParameter(1) - te
    
    parts = ['p', 'd']
    lines = makeLines(h, eoff, parts, momentum)
    stuff.append(lines)    
    
    return fit


##########################################
# cts    ... central peak time of assuemd gauss
# w      ... widths
# t1, t2 ... full fit window
# sfs are 1/SFs for settign peak maginitues w.r.t. the main peak

def MultiFit(h, tag, momentum, cts, ws, t1, t2, peaksfs = [1., 1., 1.]):
    fname = 'fit{}'.format(tag)
    hname = h.GetName()
    fitform = ''
    ngpars = 3
    # possible future issue: same sigma of pi and mu?
    # fragile change to make...
    
    for ifit in range(0,len(cts)):
        #need to fix that the parameters for sigma mu and pi are the same
        #print(ifit, ngpars)
        # if ifit == len(cts)-1:
        #
        #     fitform = fitform + '[{}]*exp(-(x-[{}])^2/(2*[{}]^2))'.format(ifit*ngpars, ifit*ngpars + 1, ifit*(ngpars-2) + 2 )
        # else:
        fitform = fitform + '[{}]*exp(-(x-[{}])^2/(2*[{}]^2))'.format(ifit*ngpars, ifit*ngpars + 1, ifit*ngpars + 2 )

        if ifit < len(cts)-1:
            fitform = fitform + ' + '
    print('Fit formula: ', fitform)
    fit = ROOT.TF1(fname, fitform, t1, t2)
    # get electron peak applitude based on the expected electron arrival time!
    ampl = h.GetBinContent(h.FindBin(cts[0]))
    for ifit in range(0,len(cts)):
        fit.SetParameter(ifit*ngpars+1, ampl / peaksfs[ifit])
        fit.SetParameter(ifit*ngpars+1, cts[ifit])
        # constrain mu to mu \pm 3*sigma:
        fit.SetParLimits(ifit*ngpars+1, cts[ifit] - 3*ws[ifit],  cts[ifit] + 3*ws[ifit])
        # if ifit == len(cts)-1:
        fit.SetParameter(ifit*ngpars+2, ws[ifit])
        # this line effectively fixes the width parameters;)
        fit.SetParLimits(ifit*ngpars+2, 0.4, 0.1)
        # this line would only limit the widths, but fit is unstable:
        # fit.SetParLimits(ifit*ngpars+2, 0.1, 0.4)

    for ipar in range(0,fit.GetNpar()):
        print('par{} initially {:1.2f}'.format(ipar, fit.GetParameter(ipar)))
    print('*** Fitting!')
    h.Fit(fname, '', '' , t1, t2)
    fit.SetLineColor(ROOT.kBlack)
    fit.SetLineStyle(2)
    #fit.SetLineColor(h.GetLineColor())
    fit.SetLineStyle(2)
    fit.Draw('same')

    fits = [fit]
    for ifit in range(0,len(cts)):
        sfname = fname + str(ifit)
        fitform = '[{}]*exp(-(x-[{}])^2/(2*[{}]^2))'.format(0, 1, 2)
        single_fit = ROOT.TF1(sfname, fitform, t1, t2)
        for ip in range(0, ngpars):
            single_fit.SetParameter(ip, fit.GetParameter(ifit*ngpars + ip) )
        fits.append(single_fit)
    cols = [pcols['e'],
            pcols['mu'], pcols['pi'], ROOT.kBlack]
    for sfit,col in zip(fits[1:],cols):
        sfit.SetLineStyle(2)
        sfit.SetLineColor(col)
        sfit.Draw('same')


    te = getTof(ms['e'], momentum)
    print("TE", te, momentum, fit.GetParameter(1))

    # electrons tof offset, but for drawing, we want to draw
    # theoretical lines
    # and then fit lines
    # so we kep this ZERO!
    eoff = 0.

    fit_te = fit.GetParameter(1)
    fit_tmu = fit.GetParameter(4)
    fit_tpi = fit.GetParameter(7)
    
    parts = ['e', 'mu', 'pi']
    lines = makeLines(h, eoff, parts, momentum)
    fitlines = makeFitLines(h, parts, [fit_te, fit_tmu, fit_tpi])
    stuff.append(lines)
    stuff.append(fitlines)
    #makeOneLine()
    
    return fits


##########################################
##########################################
##########################################

# https://www.tutorialspoint.com/python/python_command_line_arguments.htm

def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    pngdir = 'png_results/'
    pdfdir = 'pdf_results/'
    os.system(f'mkdir -p {pngdir}')
    os.system(f'mkdir -p {pdfdir}')
    
    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    gBatch = False
    momentum = None
    n_spill = 1
    target = 'tun'
    #gBatch = False
    gTag=''
    print(argv[1:])
    # options that require an argument should be followed by a colon (:).
    #opts, args = getopt.getopt(argv[2:], 'hbtp:s:m:', ['help','batch','tag=', 'momentum', 'spill', 'target'])


    if gBatch:
        ROOT.gROOT.SetBatch(1)

    if len(argv) < 2:
        PrintUsage(argv)
        return


    fileName = argv[1]
    print(f'Opening input file "{fileName}"')
    inFile = ROOT.TFile(fileName, "READ")

    momentum = None
    runindex = fileName.index('run')
    srun = fileName[runindex+6:runindex+9]
    if momentum == None:
        momentum = getMomentum(srun)
    if momentum == None:
        momentum = getMergedMomentum(srun)

    print(f'Assuming run {srun} and momentum {momentum}')

    suff = ''
    if abs(momentum) < 500:
        suff = 'Low_act1cuts'
        print('Low momentum run, looking at zoomed version of tof histos!')
    hTOFOther = inFile.Get("hTOFOther" + suff)
    hTOFEl = inFile.Get("hTOFEl" + suff)


    signedmomentum = str(abs(momentum))
    if momentum > 0:
        signedmomentum = signedmomentum + 'Pos'
    else:
        signedmomentum = signedmomentum + 'Neg'
    canname = 'FitToF_run{}_{}'.format(srun, signedmomentum)
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    cans.append(can)

    #can.Divide(2,2)
    
    ROOT.gStyle.SetOptStat(0)
    
    hTOFOther.SetLineColor(ROOT.kBlack)
    hTOFOther.SetMarkerColor(hTOFOther.GetLineColor())
    hTOFOther.SetMarkerStyle(20)
    hTOFOther.SetMarkerSize(1)
    hTOFOther.SetLineWidth(1)
    # hack
    #hTOFOther.SetMaximum(1000)
    hTOFOther.GetYaxis().SetTitle('Events')
    hTOFOther.GetYaxis().SetMoreLogLabels()
    hTOFOther.GetXaxis().SetTitleOffset(1.15)
    hTOFOther.Draw('e1')

    off = 0.
    if 'uncalibrated' in inFile.GetName():
        off = 3.
    width = 0.2

    tofDiff_e_mu = getTofDiff('e','mu', momentum)
    tofDiff_e_pi = getTofDiff('e','pi', momentum)
    tofDiff_e_p = getTofDiff('e','p', momentum)
    tofDiff_e_d = getTofDiff('e','d', momentum)
    tofDiff_mu_pi = getTofDiff('mu','pi', momentum)

    print(f'ToF diffs for momentum {momentum}: mu-e: {tofDiff_e_mu:2.2f}, pi-e: {tofDiff_e_pi:2.2f}, p-e: {tofDiff_e_p:2.2f}, d-e: {tofDiff_e_d:2.2f}, mu-pi: {tofDiff_mu_pi:2.2f}')

    if abs(momentum) <= 300:
        ###################################################################################################
        print('Assuming low momentum run, will try to fit e/mu/pi')

        sigmat = 0.23
        tce = l/c*conv # see tofUtil

        # generic fit function arguments:
        fa = [tce, tce + tofDiff_e_mu, tce + tofDiff_e_pi,
              sigmat, sigmat, sigmat,
              1.0, 10.0, 30.,
              tce - 4*sigmat,  tce + tofDiff_e_pi + 3*sigmat]
        
        # -280:
        #if abs(momentum + 280) < 0.1:
        #    fa = [tce, 12.6, 13.1, sigmat, sigmat, sigmat, 1.0, 10.0, 30., 10.5, 14.5]

        if abs(abs(momentum) - 200) < 41:
            # smaller pion peak at low momenta:
            fa[7] = 50
            fa[8] = 200
            # and wider:
            #fa[5] = 0.24
            #fa[6] = 0.25
        
        # central times
        tcs = [fa[0], fa[1], fa[2]]
        # widths
        ws = [fa[3], fa[4], fa[5]]
        # scale factors for peak initial magnitudes
        sfs = [fa[6], fa[7], fa[8]]
        # fit range
        begin_fit = fa[9]
        end_fit = fa[10]

        # MULTIFIT!
        fits = MultiFit(hTOFOther, '_mupi', abs(momentum), tcs, ws, begin_fit, end_fit, sfs)
        stuff.append(fits)

        print("Integral fits[1]", fits[1].Integral(10, 20)/ hTOFOther.GetBinWidth(1))
        i1 = 0 
        i2 = hTOFOther.GetXaxis().GetXmax() 
        n_e  = fits[1].Integral(i1, i2) / hTOFOther.GetBinWidth(1)
        n_mu = fits[2].Integral(i1, i2) / hTOFOther.GetBinWidth(1)
        n_pi = fits[3].Integral(i1, i2) / hTOFOther.GetBinWidth(1)

        print("Integral fits[2]", fits[2].Integral(0, 100)/ hTOFOther.GetBinWidth(1))
        print("Integral fits[3]", fits[3].Integral(0, 100)/ hTOFOther.GetBinWidth(1))

        print("mu/all = ", n_mu/(n_e+n_mu+n_pi))
        print("pi/all = ", n_pi/(n_e+n_mu+n_pi))

        # fitted tof peak positions:
        t_e = fits[1].GetParameter(1)
        t_mu = fits[2].GetParameter(1)
        t_pi = fits[3].GetParameter(1)

        # relative unc. in peak integrals:
        err_e = fits[0].GetParError(0)/fits[0].GetParameter(0);
        err_mu = fits[0].GetParError(3)/fits[0].GetParameter(3);
        err_pi = fits[0].GetParError(6)/fits[0].GetParameter(6);

        # fitted widths, but can be fixed in fit!
        sig_e = fits[1].GetParameter(2)
        sig_mu = fits[2].GetParameter(2)
        sig_pi = fits[3].GetParameter(2)

        # uncertainties in fitted tof peak positions:
        sigmate = fits[0].GetParError(1)
        sigmatmu = fits[0].GetParError(4)
        sigmatpi = fits[0].GetParError(7)

        print('width mu: %.2fns, '% sig_mu, 'width pi %.2fns' % sig_pi, ', width mu/width pi %.2f'%(fits[0].GetParameter(5)/fits[0].GetParameter(8)))

        txtSize = 0.035
        xx = 0.55
        
        tex = ROOT.TLatex(xx, 0.85, 'e: t=' + '{:1.2f} #pm {:1.2f}'.format(t_e, sigmate) + 'ns, N=' + '{:1.0f}'.format(n_e) + ' #pm ' + '{:1.0f}'.format(n_e*err_e))
        # + ', n/spill= ' + '{:1.0f}'.format(n_e/n_spill) + '#pm' + '{:1.0f}'.format(n_e/n_spill*err_e))

        tex3 = ROOT.TLatex(xx, 0.80, '#mu: t=' + '{:1.2f} #pm {:1.2f}'.format(t_mu, sigmatmu) + 'ns, N=' + '{:1.0f}'.format(n_mu) + ' #pm ' + '{:1.0f}'.format(n_mu*err_mu))
        # + ', n/spill=' + '{:1.0f}'.format(n_mu/n_spill)  + '#pm' + '{:1.0f}'.format(n_mu/n_spill*err_mu))

        tex2 = ROOT.TLatex(xx, 0.75, '#pi: t=' + '{:1.2f} #pm {:1.2f}'.format(t_pi, sigmatpi) + 'ns, N=' + '{:1.0f}'.format(n_pi) + ' #pm ' + '{:1.0f}'.format(n_pi*err_pi))
        # + ', n/spill=' + '{:1.0f}'.format(n_pi/n_spill)  + '#pm' + '{:1.0f}'.format(n_pi/n_spill*err_pi))

        #use the separation to t_e instead of the absolute value
        #mom_pred_mu = TofToMomentum(t_mu - t_e + getTof(ms['e'], momentum), ms['mu'])
        #mom_pred_pi = TofToMomentum(t_pi - t_e + getTof(ms['e'], momentum), ms['pi'])
        
        mom_pred_mu, mom_pred_mu_err = TofDiffToMomentum(t_mu - t_e, ms['mu'], sigmate, sigmatmu)
        mom_pred_pi, mom_pred_pi_err = TofDiffToMomentum(t_pi - t_e, ms['pi'], sigmate, sigmatpi)
        
        xx = 0.60
        delta = 0.01
        yy = 0.68
        Delta = 0.1
        yspace = 0.05
        
        tex4 = ROOT.TLatex(xx, yy, '#hat{p} from t_{#mu}: ' + '{:1.1f}'.format(mom_pred_mu) + ' #pm {:1.1f} MeV/c'.format(mom_pred_mu_err))
        tex5 = ROOT.TLatex(xx, yy - 1*yspace, '#hat{p} from t_{#pi}: ' + '{:1.1f}'.format(mom_pred_pi) + ' #pm {:1.1f} MeV/c'.format(mom_pred_pi_err))
        tex4b = ROOT.TLatex(xx, yy - 2*yspace - 1*delta, '#hat{p}-p from t_{#mu}: ' + '{:+1.1f}'.format(mom_pred_mu-abs(momentum)) + ' #pm {:1.1f} MeV/c'.format(mom_pred_mu_err))
        tex5b = ROOT.TLatex(xx, yy - 3*yspace - 1*delta, '#hat{p}-p from t_{#pi}: ' + '{:+1.1f}'.format(mom_pred_pi-abs(momentum)) + ' #pm {:1.1f} MeV/c'.format(mom_pred_pi_err))
        tex4c = ROOT.TLatex(xx, yy - 4*yspace - 2*delta, '#hat{p}/p from t_{#mu}: ' + '{:1.3f}'.format(mom_pred_mu/abs(momentum)) + ' #pm {:1.3f}'.format(mom_pred_mu_err/abs(momentum)))
        tex5c = ROOT.TLatex(xx, yy - 5*yspace - 2*delta, '#hat{p}/p from t_{#pi}: ' + '{:1.3f}'.format(mom_pred_pi/abs(momentum)) + ' #pm {:1.3f}'.format(mom_pred_pi_err/abs(momentum)))

        texs = [tex, tex2, tex3, tex4, tex5, tex4b, tex5b, tex4c, tex5c]
        for tx in texs:
            tx.SetNDC()
            tx.SetTextSize(txtSize)
            # https://root.cern.ch/doc/master/classTAttText.html
            if tx == tex or tx == tex2 or tx == tex3:
                tx.SetTextFont(62)
            else:
                tx.SetTextFont(42)
            tx.Draw()
        stuff.append(texs)

        outdir='fitres'

        os.system(f'mkdir -p {outdir}')
        outfile = open(f'{outdir}/fitres_{momentum}.txt', 'a')
        outfile.write('Momentum {:1.0f}MeV/c, Target: {}, Run:{}, Ne: {:1.0f}:{:1.0f}, Nmu:{:1.0f}:{:1.0f}, Npi:{:1.0f}:{:1.0f}, n_spill:{:1.0f}, sig_e:{:.2f}, sig_mu:{:.2f}, sig_pi:{:.2f}'.format(momentum,target, srun,n_e, err_e*n_e, n_mu, err_mu*n_mu, n_pi, err_pi*n_pi, n_spill, sig_e, sig_mu, sig_pi) + '\n')
        outfile.write('Fit parameters:{}'.format(fa)+ '\n'+ '\n')
        # print also the mmentum bias etc
        outfile.write('Fitted momentum using mu, with an error:{:3.2f}, sigmap:{:3.2f}'.format(mom_pred_mu, mom_pred_mu_err) + '\n')
        outfile.write('Fitted momentum using pi, with an error:{:3.2f}, sigmap:{:3.2f}'.format(mom_pred_pi, mom_pred_pi_err) + '\n')
        outfile.close()

    elif abs(momentum) <= 500:
        ###################################################################################################
        print('Assuming medium momentum run, will try to fit e + mu/pi')

        #fitting array to be saved for each configuration
        fa = [11.64, 12.6, 0.4, 0.3, 1., 1., 10.5, 15.] #380MeV #340

        tcs = [fa[0], fa[1]]
        ws = [fa[2], fa[3]]
        sfs = [fa[4], fa[5]]
        begin_fit = fa[6]
        end_fit = fa[7]

        fits = MultiFit(hTOFOther, '_mupi', momentum, tcs, ws, begin_fit, end_fit, sfs)
        stuff.append(fits)
        t_mu = max(fits[1].GetParameter(1), fits[2].GetParameter(1))
        t_e = min(fits[1].GetParameter(1), fits[2].GetParameter(1))

        if abs(t_mu - fits[1].GetParameter(1)) < 10e-4 :
            n_mu = fits[1].Integral(0, 100)/ hTOFOther.GetBinWidth(1)
            n_e = fits[2].Integral(0, 100)/ hTOFOther.GetBinWidth(2)
            err_e = fits[0].GetParError(3)/fits[0].GetParameter(3);
            err_mu = fits[0].GetParError(0)/fits[0].GetParameter(0);
            print("error ", err_mu, err_e)
        else:
            n_mu = fits[2].Integral(0, 20)/ hTOFOther.GetBinWidth(1)
            n_e = fits[1].Integral(0, 20)/ hTOFOther.GetBinWidth(2)
            err_mu = fits[0].GetParError(3)/fits[0].GetParameter(3);
            err_e = fits[0].GetParError(0)/fits[0].GetParameter(0);
            print("error :", err_mu, err_e)

        #mom_pred = TofToMomentum(t_mu+t_e-getTof(ms['e'], momentum),ms['mu'])
        mom_pred = TofDiffToMomentum(t_mu-t_e,ms['mu'])

        # #integrate the gaussians
        print("Integral fits[1]", fits[1].Integral(10, 20)/ hTOFOther.GetBinWidth(1))
        print("Integral fits[2]", fits[2].Integral(10, 20)/ hTOFOther.GetBinWidth(1))

        print("Momentum from muon ", mom_pred, "Mometum = ", momentum)
        #tex = ROOT.TLatex(0.4, 0.8, 't_{e} = '+ '{:1.2f}'.format(t_e) + ', t_{mu}=' + '{:1.2f}'.format(t_mu))
        #tex2 = ROOT.TLatex(0.4, 0.7, 'Mom. from tof mu = ' + '{:1.1f}'.format(mom_pred) + 'MeV/c')

        #tex3 = ROOT.TLatex(0.4, 0.6, 'Nmu+Npi = ' + '{:1.0f}'.format(n_mu) + ', Ne = ' + '{:1.0f}'.format(n_e))


        txtSize = 0.035

        tex = ROOT.TLatex(0.35, 0.8, 'e: t= ' + '{:1.2f}'.format(t_e) + 'ns, N= ' + '{:1.0f}'.format(n_e) + '#pm' + '{:1.0f}'.format(n_e*err_e) + ', n/spill= ' + '{:1.1f}'.format(n_e/n_spill) + '#pm' + '{:1.1f}'.format(n_e/n_spill*err_e))
        tex2 = ROOT.TLatex(0.35, 0.7, '#mu+#pi: t= ' + '{:1.2f}'.format(t_mu) + 'ns, N= ' + '{:1.0f}'.format(n_mu) + '#pm' + '{:1.0f}'.format(n_mu*err_mu) + ', n/spill= ' + '{:1.1f}'.format(n_mu/n_spill)  + '#pm' + '{:1.2f}'.format(n_mu/n_spill*err_mu))

        tex3 = ROOT.TLatex(0.35, 0.6, 'Mom_{pred} from t_{#mu} = ' + '{:1.1f}'.format(mom_pred) + ' MeV/c')
        tex2.SetTextSize(txtSize)
        tex3.SetTextSize(txtSize)
        
        tex.SetNDC()
        tex2.SetNDC()
        tex3.SetNDC()
        tex.Draw()
        tex2.Draw()
        tex3.Draw()
        stuff.append(tex)
        stuff.append(tex2)
        stuff.append(tex3)

        #
        # width = 0.24
        # fite = Fit(hTOFOther, '_el', momentum, 11.6, width, 10., 14., 1.)
        # fitmu = Fit(hTOFOther, '_p', momentum,  11.63 - tofDiff_e_mu, width*1.1, 10 - tofDiff_e_mu,  14 - tofDiff_e_mu, .5)




    else:
        ###################################################################################################
        width = 0.24
        print('Assuming high momentum run, will try to fit e/p/d')
        print(tofDiff_e_d)
        fite = Fit(hTOFOther, '_el', momentum, 11.6, width, 10., 14., 1.)
        
        fitp = Fit(hTOFOther, '_p', momentum,  11.63 - tofDiff_e_p, width, 10 - tofDiff_e_p,  14 - tofDiff_e_p, 1.)
        fitd = Fit(hTOFOther, '_d', momentum,  25.5, width+2, 24,  27, 0.8)
        stuff.append([fitp, fitd])
        
        tdif_e_p = fitp.GetParameter(1) - fite.GetParameter(1)

        mom_pred = TofToMomentum(fitp.GetParameter(1),ms['p'])

        a = (fitd.GetParameter(1)*c)/(l*conv)
        print(a)
        m3 = mom_pred*sqrt((a)**2-1)

        if momentum > 0:
            print("Mometum from proton ", TofToMomentum(fitp.GetParameter(1),ms['p']))
            tex = ROOT.TLatex(0.3, 0.8, 't_{e} = '+ '{:1.2f}'.format(fite.GetParameter(1)) + ', t_{p}=' + '{:1.2f}'.format(fitp.GetParameter(1)) + ', t_{3}=' + '{:1.2f}'.format(fitd.GetParameter(1)) + ', m_{3}=' + '{:1.0f}'.format(m3) + 'MeV')
        else:
            tex = ROOT.TLatex(0.4, 0.8, 't_{e} =' + '{:1.2f}'.format(fite.GetParameter(1)))
        tex2 = ROOT.TLatex(0.4, 0.7, 'p from proton tof  = ' + '{:1.1f}'.format(mom_pred) + 'MeV/c')

        tex.SetNDC()
        tex2.SetNDC()
        tex.Draw()
        tex2.Draw()
        stuff.append(tex)
        stuff.append(tex2)

        #tdif_e_d = fitd.GetParameter(1) - fite.GetParameter(1)
        #tex = ROOT.TLatex(0.7, 0.7, 't_{d}-t_{e}=' + '{:1.2f}'.format(tdif_e_d))
        #tex.SetNDC()
        #tex.Draw()
        #stuff.append(tex)
    
    #hTOFEl.SetLineColor(ROOT.kRed)
    #hTOFEl.SetLineWidth(2)
    #hTOFEl.Draw('hist same')

    ROOT.gPad.Update()

    cnote, pnote = makePaperLabel(srun, momentum, 0.12, 0.92)
    ROOT.gPad.Update()

    stuff.append(pnote)
    for can in cans:
        can.cd()
        cnote.Draw()
        pnote.Draw()
        can.Print(pngdir + can.GetName() + '_liny.png')
        can.Print(pdfdir + can.GetName() + '_liny.pdf')
        ROOT.gPad.SetLogy(1)
        can.Print(pngdir + can.GetName() + '_logy.png')
        can.Print(pdfdir + can.GetName() + '_logy.pdf')

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

