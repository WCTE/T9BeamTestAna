#!/usr/bin/python
# jk 16.11.2023

import data_runs

import os

isDryRun = False
isDryRun = True

runs=[ #'514_515_516', '517_518_520',
      #'339',
      #  '477', # +380 MeV/c
    '352'
        ]

selection = 'f'
#noAct1Cuts = 'true'
noAct1Cuts = 'false'


os.system('mkdir -p pdf/')

#dirname='output/'
#prefix='ntuple_'

dirname='windowpe_analyzed/'
prefix='windowPE_-16ns_45ns_run'
useWindowIntCharge = 'true'


cmd='mkdir -p histos/windowpe_analyzed/'
os.system(cmd)

for srun in runs:
    print('Processing {}'.format(srun))
    onerun = srun.split('_')[0]
    p = data_runs.getMomentum(onerun)
    cmd = 'root -l -b -q "macros/runMakeAllDataPlots.C(\\"{}{}000{}.root\\", {}, false, {}, \\"{}\\", {})"'.format(dirname, prefix, srun, p, noAct1Cuts, selection, useWindowIntCharge)
    print(cmd)
    if not isDryRun:
           os.system(cmd)
    outfile = 'windowPE_-16ns_45ns_run000{}_plots.root'.format(srun)
    if len(selection) > 0:
        outfile = 'windowPE_-16ns_45ns_run000{}_plots_{}.root'.format(srun, selection)
        
    #cmd = './python/quickPlots1d.py histos/' + outfile
    cmd = './python/slowPlots1d.py histos/windowpe_analyzed/' + outfile
    print(cmd)
    if not isDryRun:
        os.system(cmd)
    
    cmd = './python/fitToF.py histos/windowpe_analyzed/' + outfile
    print(cmd)
    if not isDryRun:
        os.system(cmd)
    
