#!/usr/bin/python

import os, sys

from data_runs import *


def main(argv):
    if len(argv) < 2:
        print('Usage: {} runnumber'.format(argv[0]))
    smomentum = argv[1]
    momentum = int(smomentum)
    runs =  getListOfRuns(momentum)
    print('momentum: ', momentum, ' Runs: ', runs)
    for irun in runs:
        srun = str(irun)
        n = runsRefractionIndexDict[irun]
        target = getTargetFullName(irun)
        slit = getSlitPerc(srun)
        print('  run: {}, momentum: {} MeV/c, n: {}, target: {}'.format(srun, momentum, n, target))



###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################


