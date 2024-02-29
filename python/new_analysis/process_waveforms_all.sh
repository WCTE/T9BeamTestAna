#!/bin/bash
shopt -s extglob

# DEFAULT testbeam daq
#datadir=/media/wcte/T7/WCTE/data/2023/
spath=/home/${user}/np/T9BeamTestAna/python/

user=$USER
datadir=/media/${user}/T7/WCTE/data/2023/

# Jiri 20.11.2023
spath=/home/qitek/work/gitlab/Matej/T9BeamTestAna/

dir=root_files/
for i in `cd $datadir/$dir ; ls root_run_*.root` ; do
  outfile="${datadir}/ntuple_files/${i/root_run/ntuple}"
#  if ! [ -f $outfile ] ; then
    script="python3 ${spath}/python/new_analysis/process_waveform_analysis.py"
    config="${spath}/config/config_hodoscope.json"
    run_number=${i//[^0-9]/}
    if [ $run_number -lt 579 ] ; then
      config="${config/hodoscope/noHodoscope}"
    fi
    echo "${script} ${datadir}/${dir}/${i} ${config} ${outfile}"
    $script "${datadir}/${dir}/${i}" "${config}" "${outfile}"
    
#  else
#    echo "NOT processing as output ${outfile} exist!"
#  fi
done
