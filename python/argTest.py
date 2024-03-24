#!/usr/bin/python

import os, sys

argv = sys.argv

print(argv)

print('Running script {}, enjoy!'.format(argv[0]))

noDraw = False

if len(argv) > 1:
    print('Got arguments:')
    for arg in argv[1:]:
        print(arg)
        if arg == 'NoDraw':
            noDraw = True

print('Configured as:')
print('noDraw: {}'.format(noDraw))


