#!/usr/bin/python


fname = 'plots2D.tex'
outf = open(fname, 'w')

runs = [403, 396, 394, 393, 392, 398, 399, 402, 449, 436, 435]

for run in runs:
    outf.write(r'\clearpage' + '\n')
    outf.write(r'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' + '\n')
    outf.write(r'\begin{figure}[p]' + '\n')
    outf.write(r'\centering' + '\n')
    outf.write(r'\begin{tabular}{cc}' + '\n')
    outf.write(r'    \includegraphics[width=0.49\textwidth]{' + f'../pdf_results/WCTEJuly2023_Quick2D_peakAnalysed_timeCorr_windInt_000{run}_plots_f_hRef_TOF_TrigScint0LC.pdf' + r'} &' + '\n')
    outf.write(r'    \includegraphics[width=0.49\textwidth]{' + f'../pdf_results/WCTEJuly2023_Quick2D_peakAnalysed_timeCorr_windInt_000{run}_plots_f_hRef_TOF_TrigScint0RC.pdf' + r'} \\' + '\n')
    outf.write(r'    \includegraphics[width=0.49\textwidth]{' + f'../pdf_results/WCTEJuly2023_Quick2D_peakAnalysed_timeCorr_windInt_000{run}_plots_f_hRef_TOF_TrigScint1LC.pdf' + r'} &' + '\n')
    outf.write(r'    \includegraphics[width=0.49\textwidth]{' + f'../pdf_results/WCTEJuly2023_Quick2D_peakAnalysed_timeCorr_windInt_000{run}_plots_f_hRef_TOF_TrigScint1RC.pdf' + r'} \\' + '\n')
    
    outf.write(r'\end{tabular}' + '\n')
    outf.write(r'\caption{The charge vs. time-of-flight for the Left and Right PMTs in TS0 (top) and TS1 (bottom) for run ' + f'{run}.' + '}' + '\n')
    outf.write(r'\label{fog:2d:' + f'{run}' + '}' + '\n')
    outf.write(r'\end{figure}' + '\n\n')
