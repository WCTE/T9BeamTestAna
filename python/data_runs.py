#!/usr/bin/python

# jiri kvita 2023

# from collections import OrderedDict

from data_runs_dicts import *


####################################################################

def isHodoscopeRun(run):
    return run >= 580

####################################################################
# https://sba.web.cern.ch/sba/targets/TargetNorth.html

def getTarget(run):
    if run >= 537:
        return 'Target3' # 200mm Al
    else:
        return 'Target1' # 200mm Be + 3mm W

def isBeryllium(run):
    if run >= 537:
        return False
    else:
        return True

####################################################################

def getMergedMomentum(srun):
    # we assume srun is actually an expression like 240p or 1000n:
    momentum = None
    print(f'getMergedMomentum: Trying to extract momentum from {srun}')
    if srun[-1:] == 'p' or srun[-1:] == 'n':
        try:
            momentum = int(srun[:-1])
            if srun[-1:] == 'n':
                momentum = momentum*(-1)
            return momentum
        except:
            print('still an issue to deduce momenyum from ', srun)
    return momentum
####################################################################

def getMomentum(srun):
    momentum = None
    run = -999
    #print(f'getMomentum: Trying to extract momentum from {srun}')
    #if srun[-1] == 'p' or srun[-1] == 'n':
    #    print('Will try to extract the momentum from assuming this is a merged positive or negative file.')
    #    srun = srun[:-1]
    try:
        run = int(srun)
    except:
        print(f'ERROR converting run {srun} to int!')
    try:
        momentum = runsDict[run]
    except:
        print(f'ERROR getting momentum information for run {run} from python/data_runs.py !')
    return momentum

####################################################################
def getListOfRuns(momentum):
    runs = []
    try:
        runs = momentaDict[momentum]
    except:
        print(f'ERROR getting the list of runs for momentum {momentum}')
    return runs

####################################################################

def getAllRuns():
    runs = []
    for run in runsDict:
        runs.append(run)
    return runs

####################################################################
