#! /bin/sh

for file in ToF_beam dEdxStudy ; do
  pdflatex ${file}.tex

done
