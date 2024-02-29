#!/usr/bin/python
# jk 16.11.2023

import data_runs

import os

runs = [ '407', '408', '409']

selection = 'f'

os.system('mkdir -p pdf/')

#dirname='output/'
#prefix='ntuple_'

dirname='windowpe_analyzed/'
prefix='windowPE_-16ns_45ns_run'
useWindowIntCharge = 'true'


cmd='mkdir -p histos/windowpe_analyzed/'
os.system(cmd)

for srun in runs:
    p = data_runs.getMomentum(srun)
    cmd = 'root -l -b -q "macros/runMakeAllDataPlots.C(\\"{}{}000{}.root\\", {}, false, true, \\"{}\\", {})"'.format(dirname, prefix, srun, p, selection, useWindowIntCharge)
    print(cmd)
    #os.system(cmd)
    outfile = 'windowPE_-16ns_45ns_run000{}_plots.root'.format(srun)
    if len(selection) > 0:
        outfile = 'windowPE_-16ns_45ns_run000{}_plots_{}.root'.format(srun, selection)
        
    #cmd = './python/quickPlots1d.py histos/windowpe_analyzed/' + outfile
    cmd = './python/slowPlots1d.py histos/windowpe_analyzed/' + outfile
    print(cmd)
    #os.system(cmd)
    
    #./python/fitToF.py histos/merged_ntuple_940p_n1p01_plots.root
