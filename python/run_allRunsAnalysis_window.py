#!/usr/bin/python
# jk 16.11.2023, II/2024

from data_runs import *

import os


#isDryRun = False
isDryRun = True


selection = 'f'
#noAct1Cuts = 'true'
noAct1Cuts = 'false'


os.system('mkdir -p pdf/')

#dirname='output/'
#prefix='ntuple_'

dirname='windowpe_analyzed/'
prefix='peakAnalysed_timeCorr_windInt_'
useWindowIntCharge = 'true'


cmd='mkdir -p histos/windowpe_analyzed/'
os.system(cmd)

rpath='/scratch/windowpe_analyzed/'

for xfname in os.popen('cd {} ; ls peakAnalysed_timeCorr_windInt_000???.root'.format(rpath)).readlines():
    filename = xfname[:-1]
    srun = ''
    momentum = None
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
    #print('---> Run {:}, p={:} MeV/c'.format(srun,momentum))

    cmd = 'root -l -b -q "macros/runMakeAllDataPlots.C(\\"{}{}000{}.root\\", {}, false, {}, \\"{}\\", {})"'.format(dirname, prefix, srun, momentum, noAct1Cuts, selection, useWindowIntCharge)

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
    #if not isDryRun:
    #    os.system(cmd)
    
