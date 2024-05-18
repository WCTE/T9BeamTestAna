#!/usr/bin/python

# JK 20.2.2023
# 26.4.2024

import os

multiruns = [
    # n=1.15
    [514, 515, 516],
    [517, 518, 520],
    [521, 522, 523],
    [524, 525, 526],
    [507, 508, 509],
    [510, 511, 512],
    [504, 505, 506],
    [498, 499, 500, 501],
    # n=1.11
    [344],
    [343],
    [342],
    [341],
    [340],
    [347, 348, 349, 350],
    [351],
    [352],
    [353],
    [354],
    # n=1.06
    [476],
    [477],
    [479],
    [480],
    [481],
    [485],
    [486],
    [487, 488, 489],
    [490],
    [491],
    # n=1.047:
    [430],
    [431],
    [432],
    [433],
    [434],
    [427],
    [426],
    [425],
    [424],
    [423],
    # n=1.03
    [414],
    [413],
    [412],
    [411],
    [410],
    [417],
    [418],
    [419],
    [420],
    [421],
    # n=1.02
    [436],
    [437, 449],
    [438],
    [439],
    [440],
    [447],
    [446],
    [445],
    [444],
    [443]
]

dirname='histos/windowpe_analyzed/'
base = 'peakAnalysed_timeCorr_windInt_000'
suff = '_plots_f.root'

#os.system('mkdir -p merged')

for mrun in multiruns:
    cmd = ''
    runs = ''
    if len(mrun) > 1:
        for run in mrun:
            cmd = cmd + dirname + base + f'{run}{suff} '
            runs = runs + f'{run}'
            if run != mrun[-1]:
                runs = runs + '_'

        merged = dirname + base[:-3] + f'{runs}{suff}'
        cmd = f'hadd {merged} ' + cmd
        print(cmd)
        os.system(cmd)
    #else:
    #    fname = base + f'{run}.root'
    #    cmd = f'cd  ; ln -s ../{fname} . ; cd -'
    #    print(cmd)
    #    #os.system(cmd)
