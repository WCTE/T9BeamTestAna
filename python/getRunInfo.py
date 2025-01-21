#!/usr/bin/python

import os, sys

from data_runs import *


def main(argv):
    if len(argv) < 2:
        print('Usage: {} runnumber'.format(argv[0]))
    srun = argv[1]
    momentum = getMomentum(srun)
    irun = int(srun)
    n = runsRefractionIndexDict[irun]
    target = getTargetFullName(irun)
    slit = getSlitPerc(srun)
    print('Run: {}, momentum: {} MeV/c, n: {}, target: {} slit: {}%'.format(srun, momentum, n, target, slit))



###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################


