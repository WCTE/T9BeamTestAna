#!/usr/bin/python
# jk 16.11.2023

import data_runs

import os

runs = [ # not available?? '381', '382',
         '383', '388' ]

selection = 'f'

os.system('mkdir -p pdf/')


for srun in runs:
    print(srun)
    p = data_runs.getMomentum(srun)
    cmd = 'root -l -b -q "macros/runMakeAllDataPlots.C(\\"data/ntuple_000{}.root\\", {}, false, true, \\"{}\\")"'.format(srun, p, selection)
    print(cmd)
    #os.system(cmd)
    outfile = 'ntuple_000{}_plots.root'.format(srun)
    if len(selection) > 0:
        outfile = 'ntuple_000{}_plots_{}.root'.format(srun, selection)
        
    #cmd = './python/quickPlots1d.py histos/' + outfile
    cmd = './python/slowPlots1d.py histos/' + outfile
    print(cmd)
    #os.system(cmd)
    
    #./python/fitToF.py histos/merged_ntuple_940p_n1p01_plots.root
