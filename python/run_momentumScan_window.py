#!/usr/bin/python
# jk 16.11.2023, II/2024
# 26.4.2024

import data_runs

import os

#isDryRun = False
isDryRun = True

runs=[ '514_515_516',
       '517_518_520',
       '339',
       '477', # +380 MeV/c
       '352', # -280
       '507_508_509',
       '510_511_512',
       '342',
       '341', '353',
       '504_505_506',
       '498_499_500_501',
       '521_522_523',
       '524_525_526'
      ]

selection = 'f'
#noAct1Cuts = 'true'
noAct1Cuts = 'false'


os.system('mkdir -p pdf/')

dirname='windowpe_analyzed/'
prefix='peakAnalysed_timeCorr_windInt_'
useWindowIntCharge = 'true'


cmd='mkdir -p histos/windowpe_analyzed/'
os.system(cmd)

for srun in runs:
    print('Processing {}'.format(srun))
    onerun = srun.split('_')[0]
    p = data_runs.getMomentum(onerun)
    thisprefix = prefix  + ''
    #if len() > 1:
    #    thisprefix = thisprefix + '000'
    cmd = 'root -l -b -q "macros/runMakeAllDataPlots.C(\\"{}{}000{}.root\\", {}, false, {}, \\"{}\\", {})"'.format(dirname, thisprefix, srun, p, noAct1Cuts, selection, useWindowIntCharge)

    print(cmd)
    if not isDryRun:
           os.system(cmd)
    
    outfile = 'peakAnalysed_timeCorr_windInt_000{}_plots.root'.format(srun)
    
    if len(selection) > 0:
        outfile = 'peakAnalysed_timeCorr_windInt_000{}_plots_{}.root'.format(srun, selection)
        
    #cmd = './python/quickPlots1d.py histos/' + outfile
    cmd = './python/slowPlots1d.py histos/windowpe_analyzed/' + outfile
    print(cmd)
    #if not isDryRun:
    #   os.system(cmd)
    
    cmd = './python/fitToF.py histos/windowpe_analyzed/' + outfile
    print(cmd)
    if not isDryRun:
        os.system(cmd)
    
