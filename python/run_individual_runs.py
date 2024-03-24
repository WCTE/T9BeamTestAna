#!/usr/bin/python
# jk 25.8.2023

from data_runs import *
from data_badruns import *

import os, sys

os.system('mkdir -p output')

dryrun = True # default
#dryrun = False # Careful!!

if not dryrun:
    print("WARNING, this will runMakeAllDataPlots for all runs in your output/, do you really wish to continue? Y/n")
    a = input()
    if not (a == 'Y' or a == 'y'):
        exit(1)
    
for xlistname in  os.popen('cd output/ ; ls ntuple_*.root'):
    
    print('######################################################')

    rfilename = xlistname[:-1]
    base = rfilename.replace('.root','')
    #print('base: ', base)
    srun = ''
    tokens = base.split('_')
    srun = tokens[-1].replace('000','')
    #print(srun)
    
    run = srun
    momentum = getMomentum(int(srun))
    #print(srun, momentum)

    if momentum == None:
        print('ERROR getting the momentum, continuing!')
        continue


    
    isHodo = isHodoscopeRun(int(run))
    CisHodo = 'false'
    if isHodo:
        isHodo = 'true'


    if int(run) <= 253:
        # commissioning
        continue
    
    # HACK!
    if isHodo:
        continue

    
    noAct1Cuts = 'true'
    cmd = 'root -l -b -q "macros/runMakeAllDataPlots.C(\\"output/{}\\", {}, {}, {})"'.format(rfilename,momentum,CisHodo, noAct1Cuts)
    print(cmd)
    if not dryrun:
        os.system(cmd)

    hfilename = rfilename.replace('.root','_plots.root')
    
    cmd='./python/quickPlots1d.py histos/{}'.format(hfilename)
    print(cmd)
    if not dryrun:
        cmd = cmd + ' -b'
        os.system(cmd)
        
    cmd='./python/slowPlots1d.py histos/{}'.format(hfilename)
    print(cmd)
    if not dryrun:
        cmd = cmd + ' -b'
        os.system(cmd)
        
    cmd='./python/fitToF.py histos/{}'.format(hfilename)
    print(cmd)
    #if not dryrun:
    #    cmd = cmd + ' -b'
    #    os.system(cmd)
        





